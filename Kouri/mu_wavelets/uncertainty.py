
import scipy.integrate
from math import *

hermite = scipy.special.hermite





class mu_wavelet:
    
    def __init__(self, n ):
        self.sqrt_n = sqrt(float(n))
        self.n = n
        self.H_n = hermite(n)
    
    def __call__(self, x ):
        exp(log(self.H_n(x))-self.sqrt_n*x*x)
    

















"""
def Var_x_integrand( x, n ):

def Var_k_integrand( k, x ):
    



scipy.integrate.dblquad(cylinder_integrand, 0., _2pi, lambda x:0., lambda x:H, args=(x,) )

"""




