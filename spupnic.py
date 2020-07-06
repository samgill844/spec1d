########################################
# IMPORTS
########################################
import matplotlib.pyplot as plt
import numpy as np, os, sys
import argparse
from astropy.io import fits 
from astropy.visualization import simple_norm
from scipy.signal import find_peaks
from scipy.optimize import curve_fit, differential_evolution
from scipy.interpolate import interp2d
from tqdm import tqdm 
from astropy.modeling import models
import emcee, corner
from scipy.interpolate import interp1d
from multiprocessing import Pool
from astropy.table import Table
import matplotlib.patches as patches
from astropy.table import Table, Column, Row
import glob 
from scipy.ndimage import median_filter, maximum_filter
from lmfit import minimize, Parameters, Minimizer, printfuncs, conf_interval
import lmfit 
from lmfit.models import GaussianModel, PolynomialModel, VoigtModel,LorentzianModel
plt.style.use('dark_background')

########################################
# Argumant parser
########################################
description = '''Spupnic reduction code'''


parser = argparse.ArgumentParser('spupnic', description=description)

parser.add_argument('--create_log', action="store_true", default=False, help="Create the log file")

parser.add_argument('-a', 
                    '--extract_CuNe',
                     help='Extract spectrum from a fits file', type=str, default='no')

parser.add_argument('-b', 
                    '--extract_science',
                     help='Extract science spectra from a log file', type=str, default='no')


parser.add_argument('-c', 
                    '--path_to_ref_spectra',
                     help='The path to the reference spectra.', type=str, default='/home/sam/Software/SpUpNIC/CuNe_ref_spectra')

parser.add_argument('-d', 
                    '--threads',
                     help='The number of threads to use.', type=int , default=12)


parser.add_argument('--extract_all_CuNe', action="store_true", default=False, help="Extract all CuNe lamp spectra")
parser.add_argument('--extract_all_science', action="store_true", default=False, help="Extract all science spectra")
parser.add_argument('--process_folder', action="store_true", default=False, help="Process entire night of data.")


########################################
# Useful functions
########################################

def box(x, p):
    height, center, width = p
    return height*(center-width/2 < x)*(x < center+width/2)

def gaus(x,p):
    a,x0,sigma = p
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def sum_of_gaussians(x, thetas):
    flux_ = np.zeros(x.shape[0])
    for i in range(len(thetas)//3) : flux_ = flux_ + gaus(x, thetas[3*i:3*i+3])
    return flux_

def gaus_lnlike(theta, x, flux, bounds):
    for i in range(len(theta)):
        if (theta[i] < bounds[i][0]) or (theta[i] > bounds[i][1]) : return -np.inf
    
    model = sum_of_gaussians(x, theta)
    return -0.5*np.sum( (flux - model)**2 )

class wavespreader:
    def __init__(self, datax,datay, refx,refy):
        self.datax, self.datay, self.refx, self.refy = datax,datay, refx,refy
        
    def __call__(self, x, pars, return_wave=False): return self.residuals_lmfit(pars, return_model = True, x=x, return_wave=return_wave)

    def residuals_lmfit(self, pars, return_model = False, x = False, return_wave=False) :
        if return_model : 
            if return_wave : return np.interp(x , self.refx, self.refx - float(pars['c0']) + float(pars['c2'])*(self.refx-float(pars['c1'])) + float(pars['c3'])*(self.refx-float(pars['c1']))**2)
            else : return np.interp(x , self.refx - float(pars['c0']) + float(pars['c2'])*(self.refx-float(pars['c1'])) + float(pars['c3'])*(self.refx-float(pars['c1']))**2, self.refy) 
        else : return np.abs(self.datay - np.interp(self.datax , self.refx - float(pars['c0']) + float(pars['c2'])*(self.refx-float(pars['c1'])) + float(pars['c3'])*(self.refx-float(pars['c1']))**2, self.refy))

def align_CuNe_spectra(ref_spectra, traced_flux):
    # First, sort out the two data sets
    refx = np.arange(ref_spectra.shape[0])
    refy = ref_spectra

    datax = np.arange(traced_flux.shape[0])
    datay = traced_flux

    spreader = wavespreader(datax,datay, refx,refy)

    dw = np.argmax(datay) - np.argmax(refy)

    theta0 = np.array([30 - dw,1500-dw,0.,0.])

    params = Parameters()
    params.add('c0', value=theta0[0],     min=theta0[0] - 50, max=theta0[0] + 50)
    params.add('c1', value=theta0[1],     min=theta0[1] - 10, max=theta0[1] + 10)
    params.add('c2', value=theta0[2],     min=theta0[2] - 1e-1, max=theta0[2] + 1e-1)
    params.add('c3', value=theta0[3],     min=theta0[3] - 1e-5, max=theta0[3] + 1e-5)


    mini = Minimizer(spreader.residuals_lmfit, params, nan_policy='propagate')
    lmfit_minimize_result = mini.minimize(method='differential_evolution')

    print('LMFIT report')
    printfuncs.report_fit(lmfit_minimize_result, min_correl=0.5)

    return spreader, lmfit_minimize_result


def extract_all_CuNe_worker(i):
    os.system('spupnic --extract_CuNe {:} --ref_spectra ../arc.dat'.format(t['FILE'][i], args.ref_spectra))

def extract_all_science_worker(i):
    os.system('spupnic --extract_science {:}'.format(t['FILE'][i]))

if __name__=="__main__":
    # First, parse the args 
    args = parser.parse_args()

    # Create the logfile
    if args.create_log:
        files = glob.glob('a*.fits')
        files.sort() 

        time = [] # HJD-OBS
        target = [] # OBJECT
        #ra = [] # TARG-RA
        #dec = [] # TARG-DEC
        arc_lamp = [] # ARC-LAMP
        exptype = [] # EXPTYPE
        for i in range(len(files)):
            header = fits.open(files[i])[0].header 
            time.append(float(header['HJD-OBS']))
            target.append(header['OBJECT'])
            arc_lamp.append(header['ARC-LAMP'])
            exptype.append(header['EXPTYPE'])
    
        t = Table()
        t.add_column(Column(files, name='FILE'))
        t.add_column(Column(np.array(time) + 2400000, name='HJD'))
        t.add_column(Column(target, name='TARGET'))
        t.add_column(Column(arc_lamp, name='ARC-LAMP'))
        t.add_column(Column(exptype, name='EXPTYPE'))

        t.write('log.fits', overwrite=True)




    if args.extract_CuNe!='no':

        # First, extract the image
        h = fits.open(args.extract_CuNe)
        grating = h[0].header['GRATING']
        gratingangle = int(str(h[0].header['GR-ANGLE']).split('.')[0])

        norm = simple_norm(h[0].data, 'sqrt', min_percent=5, max_percent=95)

        # Now plot the image
        f1 = plt.figure()
        ax1 = plt.gca()
        plt.imshow(h[0].data, aspect='auto', norm=norm, origin='lower', cmap='Greys')
        plt.axhline(40, c='r', alpha = 0.3, ls = '--')
        plt.axhline(100, c='r', alpha = 0.3, ls = '--')
        ax1.set(xlabel='X', ylabel='Y', title='{:} [grating {:} Angle {:}]'.format(args.extract_CuNe, grating, gratingangle))


        # Get first spectra corrosponding to the bounds in the above plot
        # This is only a rough spectra, and not properly traced. 
        initial_spectra = np.sum(h[0].data[40:100,:], axis=0)[::-1]
        initial_spectra = initial_spectra / initial_spectra.max() 

        # To trace the spectral line, we will find peaks in the initial_spectra
        # and fit a box to it. 
        peaks, _ = find_peaks(initial_spectra, height = np.std(initial_spectra)) 
        peaks = np.sort(peaks)
        if len(peaks) > 50 : 
            print('Too many spectral lines... problem?')
            f1.savefig('{:}_image.png'.format(args.extract_CuNe.split('.')[0]), bbox='tight')
            exit()

        f2, (ax2, ax3, ax4) = plt.subplots(nrows=3, ncols=1, figsize=(15,10), sharex=True)
        ax2.set(ylabel = 'Normalised Counts', xlabel='Pixel', title = args.extract_CuNe)
        ax3.set(ylabel = 'Normalised Counts', xlabel='Pixel')
        ax4.set(ylabel = 'Normalised Counts', xlabel='Pixel')

        ax2.plot(np.arange(initial_spectra.shape[0]), initial_spectra, c='g', alpha = 0.4, ls='--', label='Initial spectra')
        ax2.plot(np.arange(initial_spectra.shape[0])[peaks], initial_spectra[peaks], 'r+', label='Initial spectra peaks')


        # Now we get the vertical profile at each peak
        # This will allow us to fit a box to the slit 
        # for each line, and thus trace it efficiently. 
        upper_line, lower_line = np.zeros(len(peaks)), np.zeros(len(peaks))
        for j in range(len(peaks)):
            vertical_spectra = h[0].data[:, h[0].data.shape[1] - peaks[j]] # extract vertical profile
            vertical_spectra = vertical_spectra - np.percentile(vertical_spectra, 20) # slide the baseline down a bit 
            x = np.arange(vertical_spectra.shape[0]) # create an "x" axis to fit the box to
            theta = [np.percentile(vertical_spectra, 80), 80, 80] # create a starting guess with height, center and width
            res = differential_evolution(lambda p: np.sum((box(x, p) - vertical_spectra)**2), [[i-20, i+20] for i in theta])
            upper_line[j] = res.x[1] + res.x[2]/2 # use differential evolution to fit the box
            lower_line[j] = res.x[1] - res.x[2]/2

            # Add patch to image tracing the spectral line
            rect = patches.Rectangle((-15 + initial_spectra.shape[0]-np.arange(initial_spectra.shape[0])[peaks][j],upper_line[j]),30,lower_line[j] - upper_line[j],linewidth=1,edgecolor='g',facecolor='none', alpha = 0.8)
            ax1.add_patch(rect)


        # Now we will fit a line to the top and bottom of the boxes to trace the 
        # spectral lines
        z_upper = np.polyfit(initial_spectra.shape[0]-np.arange(initial_spectra.shape[0])[peaks],  upper_line, 1)
        z_lower = np.polyfit(initial_spectra.shape[0]-np.arange(initial_spectra.shape[0])[peaks],  lower_line, 1)
        p_upper = np.poly1d(z_upper)
        p_lower = np.poly1d(z_lower)
        avg_m = (z_upper[0] + z_lower[0]) / 2.
        dc = abs(z_upper[1] - z_lower[1])
        dist = dc / np.sqrt(avg_m**2 + 1)
        xs = np.sqrt(dc**2 - dist**2)

        p_x = initial_spectra.shape[0]-np.arange(initial_spectra.shape[0])
        ax1.plot(p_x, p_upper(p_x), 'g',linewidth=1, ls='--')
        ax1.plot(p_x, p_lower(p_x), 'g',linewidth=1, ls='--')

        # Thats it for the image, so let's save and close
        ax1.set_xlim(0,h[0].data.shape[1])
        ax1.set_ylim(0,h[0].data.shape[0])

        f1.tight_layout()
        f1.savefig('{:}_image.png'.format(args.extract_CuNe.split('.')[0]), bbox='tight')
        plt.close(f1)

        # Now we need to extract a spectrum from inbetween the green lines
        # which mark the spectral lines. 
        image_interpolater = interp2d(np.arange(h[0].data.shape[1]), np.arange(h[0].data.shape[0]), h[0].data)
        npoints = 100
        traced_flux = np.zeros(initial_spectra.shape[0])
        for i in range(initial_spectra.shape[0]):
            y = np.linspace(p_lower(i), p_upper(i), npoints)
            x = np.linspace(i, i + xs, npoints )
            traced_flux[i] = np.sum(image_interpolater(x,y))*dist / npoints
        traced_flux = traced_flux[::-1] / traced_flux.max() # normalise and reverse flux
        traced_flux = traced_flux - np.percentile(traced_flux, 20)
        ax2.plot(np.arange(initial_spectra.shape[0]), traced_flux, c='y', alpha = 0.4, ls='--', label = 'Traced spectra')


        # Now we need to load in the reference spectra 
        # that corrosponds to the mask        
        ref_spectra = np.loadtxt( '{:}/{:}_{:}_ref_spectra.dat'.format(args.path_to_ref_spectra, grating, gratingangle))
        ref_spectra = ref_spectra[::-1] / ref_spectra.max()
        ax3.plot(np.arange(ref_spectra.shape[0]), ref_spectra, c='r', alpha = 0.5, ls='-', label='reference spectra')
        ax3.plot(np.arange(traced_flux.shape[0]), traced_flux , c='y', alpha = 0.5, ls='-', label='Trace spectra')

        # Now we need to get the stretch and slide to match out spectrum to 
        # the reference spectrum.
        # We want a transorm of the pixel axis to the referece pixel axis
        # so we can easily match spectral lines
        np.save('ref_traced', np.array([ref_spectra, traced_flux]))
        spreader, lmfit_minimize_result  = align_CuNe_spectra(ref_spectra, traced_flux)
        ax3.plot(np.arange(ref_spectra.shape[0]), spreader(np.arange(ref_spectra.shape[0]), lmfit_minimize_result.params), 'b--', label='Distorted model', alpha = 0.5)
        calibration = np.array([spreader( np.load('{:}/{:}_{:}_ref_spectra.npy'.format(args.path_to_ref_spectra, grating, gratingangle))[:,1], lmfit_minimize_result.params, return_wave=True), np.load('{:}/{:}_{:}_ref_spectra.npy'.format(args.path_to_ref_spectra, grating, gratingangle))[:,0]]).T # the pixel positions of the wavelength calibration

        # Now the problim with calibration is that some peaks are there and some are not. 
        # Lets do a peak find of the traced flux to see what peaks match up with the calibration
        traced_flux_peaks, _ = find_peaks(traced_flux, height = 0.6*np.std(traced_flux)) 
        ax4.plot(np.arange(traced_flux.shape[0]), traced_flux, c='y', alpha = 0.8, ls='--', label='Trace spectra')
        ax4.scatter(np.arange(traced_flux.shape[0])[traced_flux_peaks], traced_flux[traced_flux_peaks], c='y', marker = "2", label='Traced spectra peaks')
        for cal in range(len(calibration)) : ax4.axvline(calibration[cal,0], ls = '-', c='b', alpha = 0.5, label = 'Calibration peaks' if cal ==0 else None)

        # Find the ones we have peaks for 
        # We need to match calibration[cal,0] with np.arange(traced_flux.shape[0])[traced_flux_peaks]
        mask = np.zeros(calibration.shape[0], dtype = np.bool)
        for i in range(len(calibration)):
            if np.sum(np.abs(calibration[i,0] -   np.arange(traced_flux.shape[0])[traced_flux_peaks]) < 5) ==1 : mask[i] = True
        calibration = calibration[mask]
        ax4.scatter(calibration[:,0], np.interp(calibration[:,0], np.arange(traced_flux.shape[0]), traced_flux), s=80, facecolors='none', edgecolors='r')

        # we know which peaks to calibrate off, we still need to fit a Gaussian to them so we can get the exact peak position instead of just the highest pixel

        pars = []
        Model = None
        i = 0
        gaussmodel =  GaussianModel(prefix='f{:}_'.format(i))
        pars = gaussmodel.guess(np.arange(traced_flux.shape[0]), traced_flux)
        pars['f{:}_center'.format(i)].set(vary=True, value = calibration[i,0], min = calibration[i,0]-5 , max = calibration[i,0]+5)
        pars['f{:}_amplitude'.format(i)].set(vary=True,value = 0.5,)
        pars['f{:}_sigma'.format(i)].set(vary=True, value = 2, min=0.5, max = 10)

        for i in range(1,calibration.shape[0]):
            gaussmodel_ = GaussianModel(prefix='f{:}_'.format(i))
            pars_ = gaussmodel_.guess(np.arange(traced_flux.shape[0]), traced_flux)
            pars_['f{:}_center'.format(i)].set(vary=True, value = calibration[i,0], min = calibration[i,0]-5 , max = calibration[i,0]+5)
            pars_['f{:}_amplitude'.format(i)].set(vary=True, value = 0.5)
            pars_['f{:}_sigma'.format(i)].set(vary=True, value = 2, min=0.5, max = 10)
            gaussmodel = gaussmodel +  gaussmodel_
            pars = pars + pars_
        out = gaussmodel.fit(traced_flux, pars, x=np.arange(traced_flux.shape[0]))
        print(out.fit_report())

        ax4.plot(np.arange(traced_flux.shape[0]), out.best_fit, c='r', alpha = 0.8, ls='--', label='Fitted Gaussians')
        for i in range(calibration.shape[0]) : calibration[i,0] = float(out.params['f{:}_center'.format(i)].value)
        for cal in range(len(calibration)) : ax4.axvline(calibration[cal,0], ls = '--', c='b', alpha = 0.5, label = 'Re-calibrated peaks' if cal ==0 else None)

        ax4.legend()
        
        print('Calibration:')
        bbox = dict(boxstyle="round", fc="0.8")
        arrowprops = dict(
            arrowstyle = "->",
            connectionstyle = "angle,angleA=0,angleB=90,rad=10")

        for cal in calibration:
            print('\t{:} -> {:} nm'.format(*cal))
            s = '{:.1f} nm'.format(cal[1])
            xy = (cal[0], np.interp(cal[0], np.arange(initial_spectra.shape[0]), traced_flux))
            xytext = (cal[0] , np.interp(cal[0], np.arange(initial_spectra.shape[0]), traced_flux) + 0.1)
            ax2.annotate(s, xy=xy, xytext=xytext, bbox=bbox, arrowprops=arrowprops)
        ax2.legend()
        ax3.legend()

        np.save('{:}_calibration'.format(args.extract_CuNe.split('.')[0]), calibration)

        ax2.set_xlim(0, len(traced_flux))
        ax2.set(title='{:} [grating {:} Angle {:}]'.format(args.extract_CuNe, grating, gratingangle))

        f2.tight_layout()
        f2.savefig('{:}_CuNe.png'.format(args.extract_CuNe.split('.')[0]), bbox='tight')

        plt.close()


    if args.extract_all_CuNe:
        # Load the log file
        t = Table.read('log.fits')

        # now filter all those with CuNe lamp spectra
        t = t[(t['EXPTYPE']=='ARC') & (t['ARC-LAMP']=='CuNe')]

        # Now multiprocess
        with Pool(args.threads) as pool:
            pool.map(extract_all_CuNe_worker, range(len(t)))



    if args.extract_science != 'no':
        h = fits.open(args.extract_science)
        norm = simple_norm(h[0].data, 'sqrt', min_percent=5, max_percent=95)

        f1 = plt.figure()
        plt.imshow(h[0].data, aspect='auto', norm=norm, origin='lower', cmap = 'Greys')
        plt.axhline(40, c='r', alpha = 0.3, ls = '--')
        plt.axhline(100, c='r', alpha = 0.3, ls = '--')
        plt.xlabel('X')
        plt.ylabel('Y')        

        # Now find ThAr within specified time limit
        t = Table.read('log.fits') 
        time_of_obs = t['HJD'][np.where(t['FILE']==args.extract_science)[0][0]]
        arc_table = t[(t['EXPTYPE']=='ARC') & (t['ARC-LAMP']=='CuNe')]
        arc_table = arc_table[np.abs(arc_table['HJD'] - time_of_obs) < 1/24] 
        if len(arc_table) == 0 : 
            #raise ValueError('No arc calibrations within 1 hr')
            print('No arcs within 1 hour of obs :(')
            exit()

        # Ideally, we want to take the average of the one before and the one after.
        # The backup is to use the closest in time
        FAILED_WAVECAL = False
        arc_table_before = arc_table[arc_table['HJD'] < time_of_obs]
        arc_table_after= arc_table[arc_table['HJD'] > time_of_obs]
        print('Number of calibrations before {:}'.format(len(arc_table_before)))
        print('Number of calibrations adter  {:}'.format(len(arc_table_after)))

        if (len(arc_table_before)!=0) and (len(arc_table_after) != 0):
            print('Using the average wavelength calibraration for the reference file before and after. ')
            before_filename = arc_table_before['FILE'][np.argmax(arc_table_before['HJD'])].split('.')[0] + '_calibration.npy'
            after_filename = arc_table_after['FILE'][np.argmin(arc_table_after['HJD'])].split('.')[0] + '_calibration.npy'
            print('\tnearest calibration before : {:}'.format(before_filename))
            print('\tnearest calibration after  : {:}'.format(after_filename))

            xx_before, yy_before = np.load(before_filename).T
            xx_after, yy_after = np.load(after_filename).T
            if xx_before.shape[0]!=xx_after.shape[0]:
                print('Shapes do not match in the calibration peaks suggesting that different peaks are used. Defaulting to using the nearest.')
                FAILED_WAVECAL=True
            elif False in (yy_before==yy_after):
                print('Different calibration peaks were used suggesting in the nearest calibrations. Defaulting to using the nearest.')
                FAILED_WAVECAL=True
            else:
                # This is when it is sucessful
                xcalibration = (xx_before + xx_after) / 2.
                ycalibration = np.copy(yy_before)
        else:
            FAILED_WAVECAL=True

        # Use the nearest instead
        if FAILED_WAVECAL:
            nearest_filename = arc_table['FILE'][np.argmin(arc_table['HJD'] - time_of_obs)].split('.')[0] + '_calibration.npy'
            print('Using nearest calibration file : {:}'.format(nearest_filename))
            xcalibration, ycalibration = np.load(nearest_filename).T
        

        # Now trace the thin line        
        upper_line, lower_line = [], []
        for j in range(50,1950)[::25]:
            vertical_spectra = h[0].data[:, j]
            vertical_spectra = vertical_spectra - np.percentile(vertical_spectra, 20)
            #plt.plot(vertical_spectra, 'k')
            x = np.arange(vertical_spectra.shape[0])
            theta = [np.max(vertical_spectra), np.argmax(vertical_spectra), 3]
            bounds = [[theta[0] - 100, theta[0] + 100],   [theta[1] - 5, theta[1] + 5],    [theta[2] - 2, theta[2] + 2]]
            #plt.plot(x, box(x, theta), 'r')
            res = differential_evolution(lambda p: np.sum((box(x, p) - vertical_spectra)**2), bounds)
            #plt.plot(x, box(x, res.x), 'r')
            upper_line.append(res.x[1] + res.x[2]/2 + 2)
            lower_line.append(res.x[1] - res.x[2]/2 - 2)
            #plt.show()
        
        lower_line, upper_line = np.array(lower_line), np.array(upper_line) 
        lower_line[(lower_line < (np.median(lower_line) - 5)) | (lower_line > (np.median(lower_line) + 5))] = np.median(lower_line)
        upper_line[(upper_line < (np.median(upper_line) - 5)) | (upper_line > (np.median(upper_line) + 5))] = np.median(upper_line)
        
        plt.scatter(np.arange(50,1950)[::25],  upper_line, c = 'r', s = 2)
        plt.scatter(np.arange(50,1950)[::25], lower_line, c='r', s = 2)
        
        # Now fit upper line 
        z_upper = np.polyfit(np.arange(50,1950)[::25],  upper_line, 1)
        z_lower = np.polyfit(np.arange(50,1950)[::25],  lower_line, 1)
        #z_upper, z_lower = np.loadtxt('a{:}{:}_bounds.csv'.format(arc['RUN-NO'][0], arc['FRAME'][closest_arc_idx] ), delimiter=',')
        m_avg = (z_upper[0] + z_lower[0])/2.
        z_upper[0] = m_avg
        z_lower[0] = m_avg
        p_upper = np.poly1d(z_upper)
        p_lower = np.poly1d(z_lower)

        avg_m = (z_upper[0] + z_lower[0]) / 2.
        dc = abs(z_upper[1] - z_lower[1])
        dist = dc / np.sqrt(avg_m**2 + 1)
        xs = np.sqrt(dc**2 - dist**2)

        p_x = h[0].data.shape[1]-np.arange(h[0].data.shape[1])
        plt.plot(p_x, p_upper(p_x), 'g')
        plt.plot(p_x, p_lower(p_x), 'g')

        plt.plot(p_x, p_upper(p_x)-20, 'r')
        plt.plot(p_x, p_lower(p_x)-20, 'r')

        plt.savefig('{:}_{:}_image.png'.format(args.extract_science.split('.')[0], h[0].header['OBJECT']), bbox='tight')
        plt.close()


        # now create the iamge interpolater
        image_interpolater = interp2d(np.arange(h[0].data.shape[1]), np.arange(h[0].data.shape[0]), h[0].data)
        npoints = 100
        flux = np.zeros(h[0].data.shape[1])
        for i in range(h[0].data.shape[1]):
            y = np.linspace(p_lower(i), p_upper(i), npoints)
            y_back = np.linspace(p_lower(i)-20, p_upper(i)-20, npoints)
            x = np.linspace(i, i + xs, npoints )
            flux[i] = np.sum(image_interpolater(x,y)) -  np.sum(image_interpolater(x,y_back))*dist / npoints
        
        figg, axs = plt.subplots(nrows=2, ncols=1, figsize=(15,10))
        x = np.arange(h[0].data.shape[1])[::-1]

        # Now we need to re-interpolate onto a wavelength axis
        mask = (x > xcalibration[0]) & (x < xcalibration[-1])
        wavelength = np.interp(x, xcalibration, ycalibration) # convert pixel to wavelength
        wavelength = wavelength[mask]
        flux = flux[mask]
        flux_err = np.ones(flux.shape[0])*10
        median_max = maximum_filter(median_filter(flux,150),150)

        axs[0].plot(wavelength, flux, 'k')
        axs[0].set(xlabel='Wavelength [nm]', ylabel='Counts')
        axs[1].plot(wavelength, flux/median_max, 'k')
        axs[1].set(xlabel='Wavelength [nm]', ylabel='Counts')
        plt.savefig('{:}_{:}_spectra.png'.format(args.extract_science.split('.')[0], h[0].header['OBJECT']), bbox='tight')
        plt.close()

        tmp = np.array([wavelength.tolist(),flux.tolist(),flux_err.tolist()]).T[::-1]
        np.savetxt('{:}_{:}_spectra.dat'.format(args.extract_science.split('.')[0], h[0].header['OBJECT']), tmp)
        tmp = np.array([wavelength.tolist(),(flux/median_max).tolist(),flux_err.tolist()]).T[::-1]
        np.savetxt('{:}_{:}_spectra_normlised.dat'.format(args.extract_science.split('.')[0], h[0].header['OBJECT']), tmp)


    if args.extract_all_science:
        # Load the log file
        t = Table.read('log.fits')

        # now filter all those with Science frames
        t = t[t['EXPTYPE']=='SCIENCE']

        # Now multiprocess
        with Pool(args.threads) as pool:
            pool.map(extract_all_science_worker, range(len(t)))


    if args.process_folder:
        # Get the log file
        os.system('spupnic --create_log')

        # Load the log file
        t = Table.read('log.fits')
        t = t[(t['EXPTYPE']=='ARC') & (t['ARC-LAMP']=='CuNe')]
        with Pool(args.threads) as pool:
            pool.map(extract_all_CuNe_worker, range(len(t)))

        t = Table.read('log.fits')
        t = t[t['EXPTYPE']=='SCIENCE']
        with Pool(args.threads) as pool:
            pool.map(extract_all_science_worker, range(len(t)))

