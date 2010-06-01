
from numpy import *
import pylab as p
import random
from math import *

def L(x, k, X ):
    s=1.0
    for j in range(len(X)):
        if j != k:
            s *= float(x-X[j])/(X[k]-X[j])
    return s



N = 12

f = lambda x: float(x <= 0.0)
#(cos(x*2.0*pi))

X_reg = [ -1.0 + 2.0*float(k)/(N-1)  for k in range(N) ]
f_reg = [f(x) for x in X_reg ]

theta = [ float(k)*pi/N for k in range(N+1) ]
X_cheb = [ cos(t) for t in theta ]
f_cheb = [ f(x) for x in X_cheb ]

X = [ -1.0 + 2.0*float(k)/(500-1)  for k in range(500) ]
p_reg = lambda x: sum([ L(x, k, X_reg )*f_reg[k] for k in range(len(X_reg)) ])
Y_reg = [ p_reg(x) for x in X ]

p_cheb = lambda x: sum([ L(x, k, X_cheb )*f_cheb[k] for k in range(len(X_cheb)) ])
Y_cheb = [ p_cheb(x) for x in X ]

p.plot(X, Y_reg, 'r' )
p.plot(X_reg, f_reg, 'ro')

p.plot(X, Y_cheb , 'b')
p.plot(X_cheb, f_cheb, 'bo')

p.plot(X,[f(x) for x in X], 'g')
p.show()


