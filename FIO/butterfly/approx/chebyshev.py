

from numpy import *
import pylab as p
from math import *

def T(x, n ):
    return cos(acos(x)*n)

h = 0.01
X = arange(-1.0+h, 1.0, h)

n = 8
q = 8


Z = [ 1.0*cos(float(k)*pi/(q-1)) for k in range(q) ]


#p.plot( Z, zeros(len(Z)), 'ro' )


p.plot(X, [T(x,n) for x in X])
p.show()












