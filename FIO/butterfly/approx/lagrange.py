

from numpy import *
import pylab as p



X = [ 0.0, 0.2, 0.5, 0.7, 1.3, 1.8, 2.0 ]



def L(x, k, X ):
    s=1.0
    for j in range(len(X)):
        if j != k:
            s *= float(x-X[j])/(X[k]-X[j])
    return s















