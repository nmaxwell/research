

from numpy import *
import pylab as p
import random

def L(x, k, X ):
    s=1.0
    for j in range(len(X)):
        if j != k:
            s *= float(x-X[j])/(X[k]-X[j])
    return s



N = 7


q = lambda x: (x-0.5)**2

X = array([ float(x)/(N-1) for x in range(N) ])
F = [ random.uniform(-1,1) for x in X ]
F = [ q(x)  for x in X]

def G(x, F):
    s = 0.0
    for k,f in enumerate(F):
        s += f*L(x, k, X)
    return s

Z = array([ float(x)/1000 for x in range(1000) ])

Y = [ G(z, F) for z in Z ]
p.plot(X,zeros(len(X)), 'go' )
p.plot(X,F, 'bo' )

p.plot(Z,Y)
p.plot(Z,[q(z) for z in Z], 'g' )

p.show()





"""
X = array([ float(x)/(N-1) for x in range(N) ])

Z = array([ float(x)/1000 for x in range(1000) ])

for k in range(len(X)):
    Y = [ L(z,k, X) for z in Z ]
    
    p.plot(X,zeros(len(X)), 'go' )
    p.plot( (X[k],) , (L(X[k],k, X), ), 'ro' )
    p.plot(Z,Y)

p.show()
"""
















