
import scipy.integrate
from math import *
from numpy import *

hermite = scipy.special.hermite

_2pi = 6.2831853071795865
j2pi = _2pi *1.0j
sqrt_2pi = 1.772453850905516


def int_reals( f, a=-Inf, b=Inf, args=() ):
    return scipy.integrate.quadpack.quad( f, a,b. args )

def int_reals2( f, a=-Inf,b=Inf, d=-Inf,c=Inf, args=() )
    return scipy.integrate.dblquad( f, a, b, lambda x:c, lambda x:d, args )

def FT( f ):
    return lambda k: int_reals( lambda x: f(x)*exp(x*k*1.0j ) )[0] /sqrt_2pi



class mu_wavelet:
    
    def __init__(self, n ):
        self.sqrt_n = sqrt(float(n))
        self.n = n
        self.H_n = hermite(n)
        self.N, err = int_reals( lambda x: exp(log(self.H_n(x)**2)-2.0*self.sqrt_n*x*x) )
        self.N = sqrt(self.N)
    
    def __call__(self, x ):
        return self.H_n(x)*exp(-self.sqrt_n*x*x)/self.N






if __name__ == "__main__":
    
    from numpy import *
    import pylab as p
    
    X0 = 3
    X = arange(-X0,X0,0.03)
    
    for n in range(1,10):
        
        u_n = mu_wavelet(n)
        
        int_reals( lambda x: (x*u_n(x))**2 )
        int_reals2( lambda x,k: (k*u_n(x))**2 )






