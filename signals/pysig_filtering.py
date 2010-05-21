
from pysig import *






class pyfilter:
    
    def __init__(self, impulse_response=[], center=0 ):
        """
        center is index of center tap, center=0 corresponds to causal filter
        """
        
        self.impulse_response = numpy.array(impulse_response)
        self.center = 0
    
    def __call__(self, input ):
            
        result = copy.copy(input)
        result.data = numpy.zeros(len(input))
        h = self.impulse_response
        j0 = self.center
        
        for k in range(len(input)):
            
            sum = 0.
            for j in range(len(h)):
                sum += input[k-j]*h[j-j0]
            
            result[k] = sum
        
        return result







if __name__ == "__main__":
    import random
    import pylab as p
    import numpy as n
    
    a = pysig(20, 0, 1)
    a.data = numpy.array([ random.random() for x in a.data ])
    
    F = pyfilter( n.ones(20)/20 )
    
    b = F(a)
    
    print a.data
    print b.data
    
    #print '\n\n', a.__dict__

    p.plot(a.data, 'bo')
    p.plot(b.data, 'ro')
    p.show()

