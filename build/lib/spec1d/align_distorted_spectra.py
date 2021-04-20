import numpy as np
from lmfit import minimize, Parameters, Minimizer, printfuncs, conf_interval
import lmfit 
from lmfit.models import GaussianModel, PolynomialModel, VoigtModel,LorentzianModel
import matplotlib.pyplot as plt  

class wavespreader:
    def __init__(self, datax,datay, refx,refy):
        self.datax, self.datay, self.refx, self.refy = datax,datay, refx,refy
        
    def __call__(self, x, pars, return_wave=False): return self.residuals_lmfit(pars, return_model = True, x=x, return_wave=return_wave)

    def residuals_lmfit(self, pars, return_model = False, x = False, return_wave=False) :
        if return_model : 
            if return_wave : return np.interp(x , self.refx, self.refx - float(pars['c0']) + float(pars['c2'])*(self.refx-float(pars['c1'])) + float(pars['c3'])*(self.refx-float(pars['c1']))**2)
            else : return np.interp(x , self.refx - float(pars['c0']) + float(pars['c2'])*(self.refx-float(pars['c1'])) + float(pars['c3'])*(self.refx-float(pars['c1']))**2, self.refy) 
        else : return np.abs(self.datay - np.interp(self.datax , self.refx - float(pars['c0']) + float(pars['c2'])*(self.refx-float(pars['c1'])) + float(pars['c3'])*(self.refx-float(pars['c1']))**2, self.refy))

def align_distorted_spectra(ref_spectra, traced_flux, verbose=False):
    # First, sort out the two data sets
    refx = np.arange(ref_spectra.shape[0])
    refy = ref_spectra

    datax = np.arange(traced_flux.shape[0])
    datay = traced_flux

    #plt.close()
    #plt.plot(datax,datay, 'r')
    #plt.plot(refx,refy, 'b')
    #plt.show()

    spreader = wavespreader(datax,datay, refx,refy)

    dw = np.argmax(datay) - np.argmax(refy)

    theta0 = np.array([ dw,3000,0.,0.])

    params = Parameters()
    params.add('c0', value=theta0[0],     min=theta0[0] - 50, max=theta0[0] + 50)
    params.add('c1', value=theta0[1],     min=theta0[1] - 10, max=theta0[1] + 10)
    params.add('c2', value=theta0[2],     min=theta0[2] - 1e-1, max=theta0[2] + 1e-1)
    params.add('c3', value=theta0[3],     min=theta0[3] - 1e-5, max=theta0[3] + 1e-5)


    mini = Minimizer(spreader.residuals_lmfit, params, nan_policy='propagate')
    lmfit_minimize_result = mini.minimize(method='differential_evolution')

    if verbose:
        print('Align distorted spectr fit report')
        printfuncs.report_fit(lmfit_minimize_result, min_correl=0.5)

    return spreader, lmfit_minimize_result


class spready:

    def __init__(self, x,y,refx,refy):
        self.x = x 
        self.y = y 
        self.refx = refx 
        self.refy = refy 

    def transform(self, c1, c2):
        x = self.x + c1 
        x = x + (x - 1500)*c2
        return x

    def lmfit_residuals(self, pars):
        # First, transform
        x_transformed = self.transform(pars['c1'].value, pars['c2'].value)

        # Now compare on common axis
        x_to_compare = self.refx[(self.refx > max(np.min(self.refx), np.min(x_transformed))) & (self.refx < min(np.max(self.refx), np.max(x_transformed)))]
        residuals = np.interp(x_to_compare, x_transformed, self.y) - np.interp(x_to_compare, self.refx, self.refy)
        distance = np.max(x_to_compare) - np.min(x_to_compare)
        return residuals # / distance

    def fit_spready(self,):
        params = Parameters()
        c1 = self.refx[np.argmax(self.refy)] - self.x[np.argmax(self.y)]
        params.add('c1', value=c1,     min=c1 - 50, max=c1 + 50)
        params.add('c2', value=0.5e-2,     min=1e-3, max=1e-2)

        mini = Minimizer(self.lmfit_residuals, params, nan_policy='propagate')
        self.lmfit_minimize_result = mini.minimize(method='differential_evolution')

        printfuncs.report_fit(self.lmfit_minimize_result, min_correl=0.5)

    def plot_solution(self,):
        f, ax = plt.subplots(nrows=1, ncols=1, constrained_layout=True, figsize=(15,5))
        ax.plot(self.refx, self.refy, c='k', label='Reference spectra')
        ax.plot(self.x, self.y, c='b', label='Initial spectra', alpha = 0.3)
        ax.plot(self.transform(self.lmfit_minimize_result.params['c1'], self.lmfit_minimize_result.params['c2']), self.y, c='b', ls='--', label='Fitted diostortion')
        ax.set(xlabel='Pixel', ylabel='Counts')
        plt.legend()
        plt.show()

    def convert_reference_points_to_spectra(self, refx_points):
        x_transformed = self.transform(self.lmfit_minimize_result.params['c1'], self.lmfit_minimize_result.params['c2'])
        return np.interp(refx_points, x_transformed, self.x)


'''
spreader = spready(x,y,refx,refy)
spreader.fit_spready()
print(spreader.convert_reference_points_to_spectra([493, 494]))
spreader.plot_solution()
'''
