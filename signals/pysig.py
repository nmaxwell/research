
import sys
import math
import numpy
import scipy
import pylab
import copy


floor = math.floor
ceil = math.ceil

inf = float("infinity")
_2pi = 6.2831853071795865



"""
class zero_padding():
    
    #print "non default __getitem__"
    def __call__(self, )
                if 0<= k or k < len(data):
                    return self.data[k]
                return 0.
            
            self.__getitem__ = f
            print self.__getitem__




        if extension_mode in ["per", "periodic", "cyclic" ]:
            def __getitem__(self, k):
                return self.data[k%len(data)]


"""




class pysig:
    "a basic signal (processing) class"
    
    def __init__(self):
        self.data = numpy.array([])
        self.sampling_period = 0
        self.initial_time = 0
        self.extension_mode = None
        self.define_indexing()
    
    def __init__(self, data=[], initial_time=0, sampling_period=1):
        self.data = numpy.array(data)
        self.sampling_period = sampling_period
        self.initial_time = initial_time
        self.define_indexing()
    
    def __getitem__(self, k):
        #print "default __getitem__"
        return self.data[k%len(self.data)]
    
    def __setitem__(self, k, x):
        self.data[k%len(self.data)] = x
    
    def define_indexing(self, extension_mode = None):
        
        return None
    """
        self.extension_mode = extension_mode
        
        if extension_mode == None or extension_mode in ["zero", "zero padding", "zero_padding" ]:
        
            
            
        if extension_mode in ["per", "periodic", "cyclic" ]:
            def __getitem__(self, k):
                return self.data[k%len(data)]
    
    """
    
    __contains__ = lambda self,x: x in data
    __len__ = lambda self: len(self.data)
    
#    def read_lvm(self, file_string):
#        [self.data, self.sampling_period,  self.initial_time ] = read_lvm.read_lvm(file_string)
    
    def sub(self, start=0, end=inf ):
        if ( isinstance(start,list) ):
            S = pysig()
            S.data = [self.data[i] for i in start ]
            S.initial_time = start[0]*self.initial_time+self.initial_time
            S.sampling_period = self.sampling_period
            return S
        
        if (end == inf):
            end = len(self.data)
        if ( isinstance(start,float) ):
            start = int((start-self.initial_time)/(self.sampling_period))
        if ( isinstance(end,float) ):
            end = int((end-self.initial_time)/(self.sampling_period))
        S = pysig()
        S.data = self.data[start:end]
        S.initial_time = self.time(start)
        S.sampling_period = self.sampling_period
        return S
    
    def time(self,i):
        #conver index to time
        if ( isinstance(i,int) ):
            return self.sampling_period*i+self.initial_time
        if ( isinstance(i,list) ):
            return numpy.array([self.sampling_period*j+self.initial_time for j in i])
    
    def final_time(self ):
        return self.sampling_period*len(data)+self.initial_time
    
    def total_time(self ):
        return self.sampling_period*len(self.data)
    
    def time_list(self):
        return self.time(range(0,len(self.data)))
    
    def time_array(self):
        return self.time(range(0,len(self.data)))
    
    def index(self, t):
        #convert time to nearest index
        if ( isinstance(t,list) ):
            return [round((t_i-self.initial_time)/self.sampling_period) for t_i in t]
        else:
            return round((t-self.initial_time)/self.sampling_period)
    
    def index_down(self, t):
        #convert time to index, rounded down
        if ( isinstance(t,list) ):
            return [int(floor((t_i-self.initial_time)/self.sampling_period)) for t_i in t ]
        else:
            return int(floor((t-self.initial_time)/self.sampling_period))
    
    def index_up(self, t):
        #convert time to index, rounded up
        if ( isinstance(t,list) ):
            return [int(ceil((t_i-self.initial_time)/self.sampling_period)) for t_i in t ]
        else:
            return int(ceil((t-self.initial_time)/self.sampling_period))
    
    def interpolate_linear(self, t):
        # linear interpolation
        if ( isinstance(t,list) ):
            return [ self.interpolate_linear(t_i) for t_i in t ]
        else:
            i1 = self.index_down(t)
            i2 = self.index_up(t)
            if i1 < 0 or i1 >= len(self.data) or i2 < 0 or i2 >= len(self.data):
                return 0
            if i1 == i2:
            	return self.data[i1]
            t1 = self.time(i1)
            t2 = self.time(i2)
            y1 = self.data[i1]
            y2 = self.data[i2]
            return (y2-y1)*(t-t1)/(t2-t1)+y1
    
    
    
    
    
