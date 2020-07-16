from lmfit import minimize, Parameters, Minimizer, printfuncs, conf_interval
import lmfit 
from lmfit.models import GaussianModel, PolynomialModel, VoigtModel,LorentzianModel, ExpressionModel , StepModel
import numpy as np 
import matplotlib.pyplot as plt 

#def box_lmfit(pars, x, p):
#    height, center, width = p
#    return height*(center-width/2 < x)*(x < center+width/2)

box_model = ExpressionModel('height*(center-width/2 < x)*(x < center+width/2) + off')


def fit_box_to_step_function(y,initial_guess, x = None, ):
    # Create the x axis if needed
    if x is None : x = np.arange(y.shape[0])

    
    # Create the parameters
    params = Parameters()
    params.add('height', value=initial_guess[0], min=0)
    params.add('center', value=initial_guess[1], min=np.min(x), max = np.max(x))
    params.add('width', value=initial_guess[2], min=1e-3, max = 0.5*(np.max(x) - np.min(x)) )
    params.add('off', value=0, min=np.min(y), max = 0.4*np.max(y) )


    # now fit
    out = box_model.fit(y, params, x=x)

    left = float(out.params['center'].value) - 0.5*float(out.params['width'].value)
    right = float(out.params['center'].value) + 0.5*float(out.params['width'].value)


    
    plt.figure()
    plt.plot(x,y)
    plt.plot(x, out.init_fit, 'k--', label='initial fit')
    plt.plot(x, out.best_fit, 'r-', label='best fit')
    plt.legend()
    plt.show(block=True)
    return  float(out.params['height'].value), float(out.params['center'].value), float(out.params['width'].value), float(out.params['off'].value), left, right

    