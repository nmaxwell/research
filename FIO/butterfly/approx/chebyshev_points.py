
from math import *
from numpy import *

import pylab as p



W = arange(-1,1,0.002)



def L(x, k, X ):
    s=1.0
    for j in range(len(X)):
        if j != k:
            s *= float(x-X[j])/(X[k]-X[j])
    return s

Q = lambda x: sum([ L(x,k,X)*(-1)**k for k in range(len(X)) ])


N = 4

theta = [ float(k)*pi/N for k in range(N+1) ]
X = [ cos(t) for t in theta ]

#p.plot( W, [  Q(x) for x in W] )
p.plot( W, [ Q(x) - cos(acos(x)*N) for x in W ])
p.plot( X, zeros(len(X)), 'ro' )

p.show()








