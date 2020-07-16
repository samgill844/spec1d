import numpy
from lmfit import minimize, Parameters, Minimizer, printfuncs, conf_interval
import lmfit 
from lmfit.models import GaussianModel, PolynomialModel, VoigtModel,LorentzianModel

class wavespreader:
    def __init__(self, datax,datay, refx,refy):
        self.datax, self.datay, self.refx, self.refy = datax,datay, refx,refy
        
    def __call__(self, x, pars, return_wave=False): return self.residuals_lmfit(pars, return_model = True, x=x, return_wave=return_wave)

    def residuals_lmfit(self, pars, return_model = False, x = False, return_wave=False) :
        if return_model : 
            if return_wave : return np.interp(x , self.refx, self.refx - float(pars['c0']) + float(pars['c2'])*(self.refx-float(pars['c1'])) + float(pars['c3'])*(self.refx-float(pars['c1']))**2)
            else : return np.interp(x , self.refx - float(pars['c0']) + float(pars['c2'])*(self.refx-float(pars['c1'])) + float(pars['c3'])*(self.refx-float(pars['c1']))**2, self.refy) 
        else : return np.abs(self.datay - np.interp(self.datax , self.refx - float(pars['c0']) + float(pars['c2'])*(self.refx-float(pars['c1'])) + float(pars['c3'])*(self.refx-float(pars['c1']))**2, self.refy))

def align_distorted_spectra(ref_spectra, traced_flux):
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