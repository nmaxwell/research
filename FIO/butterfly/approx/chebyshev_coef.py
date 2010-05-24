

from numpy import *
import pylab as p
from math import *
norm = linalg.norm

def choose(n,k):
    return float(factorial(n))/(factorial(n-k)*factorial(k))


class Polynomial:
    a=[]
    
    def __init__(self, a ):
        self.a = [x for x in a]
    
    def __call__(self, x):
        if len(a) == 0:
            return 0.0
        if len(a) == 1:
            return a[0]
        
        s = x*a[len(a)-1]
        for k in range(len(self.a)-1):
            a *= 
        return s
    

def T(n ):
    a = [0.0 for k in range(n+1)]
    
    for k in range( int(floor(float(n)/2.0))+1 ):
        for j in range(k+1):
            a[n-2*j] += choose(n,2*k)*choose(k,j)*(-1)**j
    
    return Polynomial(a)


h=0.01
errors = []

for N in range(50):
    
    err=0.0
    for n in range(N):
        q = T(n)
        e = norm(array([q(x)-cos(acos(x)*n) for x in arange(-1.0, 1.0, h)]))
        try:x
            err += e/norm( array([cos(acos(x)*n) for x in arange(-1.0, 1.0, h)]) )
        except:
            err += e

    print N,err
    try:
        errors.append(log10(err))
    except:
        errors.append(err)

p.plot(errors)
p.show()

