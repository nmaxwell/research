

from math import *
from pysig import *
from pysig_filtering import *
import numpy

_2pi = 6.2831853071795865

def rms(X):
    sum =0.
    for k in range(len(X)):
        sum += abs(X[k])**2
    return sqrt(sum/len(X))


def frequency_response(filter, n, N ):
    
    s = pysig([ cos((_2pi*k*n)/N) for k in range(N)])
    
    u = filter(s)
    
    return rms(u)/rms(s)


def frequency_spread(filter, N):
    
    freq = [ float(k)/N for k in range(0,N/2+1) ]
    h_hat = [ frequency_response(filter, k,N ) for k in range(0,N/2+1) ]
    
    h_hat /= numpy.sum(h_hat) 
    sum = 0.
    for k in range(len(freq)):
        sum += freq[k]*h_hat[k]
    
    return sum, (freq,h_hat)

def time_spread(filter, N ):
    
    h = [ i for i in filter.impulse_response ]
    
    time = [ float(k)/N for k in range(len(h)) ]
    
    h /= numpy.sum(h)
    sum = 0.
    for k in range(len(time)):
        sum += time[k]*h[k]
    
    return sum, (time,h)
    



def fermi_filter( N, cutoff, slope, eps = 0.001 ):
    
    N0 = int(math.ceil(N*(slope*math.log(1.0/eps-1)/4.0+cutoff)))
    
    h = [ 1./(1+math.exp(4.*(float(k)/N-cutoff)/slope))  for k in range(N0+1)]
    
    h /= numpy.sum(h)
    
    return pyfilter( h )

def box_filter( N, cutoff ):

    N0 = cutoff*N+5
    
    h = [ float((float(k)/N) <= cutoff)  for k in range(N0+1)]
    
    h /= numpy.sum(h)
    
    return pyfilter( h )





if __name__ == "__main__":
    
    import random
    import pylab as p
    from numpy import *
    
    N = 256
    
    #h = [ 1.0, 0.8, 0.5,0.2,0.05  ]
    #h= ones(0.25*N)
    
    
    F = fermi_filter( N, .1, .05 )
    #F = box_filter( N, .1 )
    
    sf, fgraph = frequency_spread(F, N)
    print sf
    st, tgraph = time_spread(F, N)
    print st
    
    print sf*st
    
    p.plot(fgraph[0],fgraph[1], 'r.' )
    p.plot(tgraph[0],tgraph[1], 'b.' )
    p.show()
    
    
    
    """
    X = range(2,N) 
    #[ 2**k for k in range(1,10) ]
    Y = [ frequency_response(F,k) for k in X ]
    X = [1./x for x in X]
    p.plot(X,Y, 'bo')
    p.show()
    """ 
    
    
    
    
    
    
    
    
    
    
    """
    m = 6
    F = pyfilter( ones(m)/m )
    
    N = 128
    for n in range(N/2+1):
        
        
        r,s = frequency_response(F, n, N )
        
        p.plot(s.data, 'bo' )
        p.savefig( "/workspace/output/scratch/" + "%03d" % n )
        p.close()
    """

