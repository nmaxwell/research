

from math import *
from numpy import *
import numpy
import pylab as p
import random
import itertools as itt

def L(x, k, X ):
    s=1.0
    for j in range(len(X)):
        if j != k:
            s *= float(x-X[j])/(X[k]-X[j])
    return s


q=20

[theta1, theta2] = mgrid[0.0: pi+pi/(q-1): pi/(q-1), 0.0: pi+pi/(q-1): pi/(q-1) ]
z1 = numpy.vectorize(lambda t1,t2: 0.5*cos(t1))
z2 = numpy.vectorize(lambda t1,t2: 0.5*cos(t2))
Z1 = z1(theta1, theta2)
Z2 = z2(theta1, theta2)

theta = arange( 0.0, pi+pi/(q-1), pi/(q-1) )
Z = [ 0.5*cos(t) for t in theta]

f = lambda x,y: cos(x*4.0*pi)*sin(y*4.0*pi)

f_ = numpy.vectorize(f)
F = f_(Z1,Z2)

L2 = lambda z1,z2,i1,i2: L(z1, i1, Z )*L(z2, i2, Z )  
g = lambda z1,z2: sum([ L2(z1,z2,i1,i2)*F[i1][i2] for (i1,i2) in itt.product(range(q),range(q)) ])


X = [ (random.uniform(-0.5,0.5),random.uniform(-0.5,0.5)) for k in range(100) ]

Y1 = [ f(z1,z2) for (z1,z2) in X ]
Y2 = [ g(z1,z2) for (z1,z2) in X ]


err = []
for k in range(len(Y1)):
    try:
        err.append( log10(abs(Y1[k] -Y2[k])/abs(Y1[k])+1E-17) )
    except:
        pass

import pylab as p
p.hist(err, bins=30)
p.show()







