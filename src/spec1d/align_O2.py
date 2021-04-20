from lmfit import minimize, Parameters, Minimizer, printfuncs, conf_interval
import lmfit 
from lmfit.models import GaussianModel, PolynomialModel, VoigtModel,LorentzianModel, ExpressionModel , StepModel
import numpy as np , os, sys
import matplotlib.pyplot as plt 
from scipy.ndimage import minimum_filter, median_filter

################################################################################
#--- iSpec directory -------------------------------------------------------------
ispec_dir = '/Users/sam/iSpec/'
sys.path.insert(0, os.path.abspath(ispec_dir))
import ispec


def normalize_whole_spectrum(star_spectrum):
    """
    Use the whole spectrum, strategy 'median+max'
    """
    #--- Continuum fit -------------------------------------------------------------
    model = "Splines" # "Polynomy"
    degree = 2
    nknots = None # Automatic: 1 spline every 5 nm
    from_resolution = 5000

    # Strategy: Filter first median values and secondly MAXIMUMs in order to find the continuum
    order='median+max'
    median_wave_range=0.05
    max_wave_range=1.0

    star_continuum_model = ispec.fit_continuum(star_spectrum, from_resolution=from_resolution, \
                                nknots=nknots, degree=degree, \
                                median_wave_range=median_wave_range, \
                                max_wave_range=max_wave_range, \
                                model=model, order=order, \
                                automatic_strong_line_detection=True, \
                                strong_line_probability=0.5, \
                                use_errors_for_fitting=True)

    #--- Continuum normalization ---------------------------------------------------
    normalized_star_spectrum = ispec.normalize_spectrum(star_spectrum, star_continuum_model, consider_continuum_errors=False)
    # Use a fixed value because the spectrum is already normalized
    return normalized_star_spectrum, star_continuum_model



def lmfit_func(pars, normed_spectra, O2_spectra, return_model=False):
    model = ispec.correct_velocity(np.copy(normed_spectra), pars['rv'])
    model['flux'] = pars['zp']*model['flux']
    o2_model = 1 - (1 - np.copy(O2_spectra['flux']))*pars['strength']

    if return_model : return model, o2_model
    else : return np.interp(O2_spectra['waveobs'], model['waveobs'], model['flux']) - o2_model


def align_with_O2(spectra, O2_spectra, R=5000):
    # Convert to ispec structures
    spectra = ispec.create_spectrum_structure(spectra.T[0], spectra.T[1], spectra.T[2])
    O2_spectra = ispec.create_spectrum_structure(O2_spectra.T[0], O2_spectra.T[1], O2_spectra.T[2])

    # Now normalise
    normed_spectra, spectra_continuum_model = normalize_whole_spectrum(spectra)
    #wfilter = ispec.create_wavelength_filter(normed_spectra, wave_base=np.min(O2_spectra['waveobs']) - 3, wave_top=np.max(O2_spectra['waveobs']) + 3)
    wfilter = ispec.create_wavelength_filter(normed_spectra, wave_base=684, wave_top=690)
    normed_spectra = normed_spectra[wfilter]

    params = Parameters()
    params.add('rv', value=0, min=-100, max=100)
    params.add('strength', value=1, min=0.7, max = 1.3)
    params.add('zp', value=1., min=0.9, max = 1.1, vary=False )
    mini = Minimizer(lmfit_func, params, fcn_args=(np.copy(normed_spectra), np.copy(O2_spectra), False))
    lmfit_minimize_result = mini.minimize()
    printfuncs.report_fit(lmfit_minimize_result, min_correl=0.5)
    
    corrected_spectra, o2_model = lmfit_func(lmfit_minimize_result.params, normed_spectra, O2_spectra, return_model=True)

    f = plt.figure()
    plt.plot(normed_spectra['waveobs'], normed_spectra['flux'], 'k', alpha = 0.3)
    plt.plot(corrected_spectra['waveobs'], corrected_spectra['flux'], 'k', label='Corrected solution')
    plt.plot(O2_spectra['waveobs'], O2_spectra['flux'], 'r', alpha = 0.3)
    plt.plot(O2_spectra['waveobs'], o2_model, 'r', label='O2 model')
    plt.gca().set(title = 'RV = {:.3f} +- {:.3f}'.format(lmfit_minimize_result.params['rv'].value , lmfit_minimize_result.params['rv'].stderr), xlabel='Wavelength [nm]', ylabel='Normalized flux')
    plt.legend()
    plt.xlim(685,690)
    corrected = ispec.correct_velocity(np.copy(spectra), lmfit_minimize_result.params['rv'].value)
    return f, np.array([corrected['waveobs'], corrected['flux'], corrected['err']]).T , lmfit_minimize_result.params['rv'].value , lmfit_minimize_result.params['rv'].stderr