
from math import *
from numpy import *
import numpy
import pylab as p
import random
import itertools as itt




N = 8
phi = lambda x,k: dot(x,k)

X = [ (float(i1)/N, float(i2)/N) for (i1,i2) in itt.product(range(N),range(N)) ]
K = [ (float(i1), float(i2)) for (i1,i2) in itt.product(range(-N/2,N/2),range(-N/2,N/2)) ]
P = [ ( sqrt(k1*k1+k2*k2)*2.0/(sqrt(2.0)*N),  ) for (k1,k2) in K ]

















