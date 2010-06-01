import time as time_module
from math import *
from numpy import *
import numpy
from waveprop import *
import pickle
import scipy.io
import pickle

norm = numpy.linalg.norm

wall_time_0 = time_module.clock()
wall_time = lambda : time_module.clock()-wall_time_0

class polygon2D:
    
    def __init__(self, vertices ):
        
        self.x_points = [ float(v[0]) for v in vertices ]
        self.y_points = [ float(v[1]) for v in vertices ]
        self.n_points = len(self.x_points)
    
    def interior_test(self, x,y):
        
        i=j=0
        c = bool(False)
        j=self.n_points-1
        while i < self.n_points:
            if ( self.y_points[i] > y ) != ( self.y_points[j] > y ) :
                if x < (self.x_points[j]-self.x_points[i]) * (y-self.y_points[i]) / ( self.y_points[j]-self.y_points[i] ) + self.x_points[i]:
                    c = not c
            j = i
            i += 1
        
        return c
    
    def __contains__(self, x):
        return self.interior_test(x[0],x[1])

class union:
    def __init__(self, set1=None, set2=None, set3=None, set4=None  ):
        self.sets=[]
        for E in [set1,set2,set3,set4]:
            if E != None:
                self.sets.append(E)
    
    def __contains__(self, x):
        for E in self.sets:
            if x in E:
                return True
        return False

class complement:
    def __init__(self, complement_set  ):
        self.complement_set = complement_set
    
    def __contains__(self, x):
        return not ( x in self.complement_set )

class simpleFunction:
    
    def __init__(self, sets_values ):
        self.sets = [ x[0] for x in sets_values ]
        self.values = [ x[1] for x in sets_values ]
        if len(self.sets) != len(self.values):
                print "error in simple_function: len(self.sets) != len(self.values):"
    
    def __call__(self, x ):
        sum = (self.values[0])*0.0
        for k,E in enumerate(self.sets):
            if x in E:
                sum += self.values[k]
        return sum




