
from math import *
from numpy import *
import numpy

norm = numpy.linalg.norm

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



#-----------------------


from hdaf import grid2d,hdafLaplacian2D,write_png
from ode import *


grid = grid2d(n1=256, n2=256, a1=-10, a2=-10, b1=10, b2=10 )


n_vertices = 20
vertices=[]
for k in range(n_vertices):
    t = 2.0*pi*float(k)/n_vertices
    vertices.append( array((cos(t),sin(t)))*(cos((t-atan2(-1,-1))*3)+2)*3.0 )

D = polygon2D(vertices)
chi_D = simpleFunction([(D,1.0)])
chi_D = grid.evaluate( lambda x,y: chi_D((x,y)) )
chi_Dc = grid.ones()-chi_D


f = lambda x,y: 5.0*exp(-norm(array((x,y))-array((-5,-5)))**2/9)
F = chi_Dc* grid.evaluate( f )

write_png( chi_D, "chi_D.png", major_scale=std(chi_D)*1, center=mean(chi_D) )
write_png( F, "F.png", major_scale=std(F)*1, center=mean(F) )


#----------------------

class terminator:
    
    def __init__(self, ode_instance ):
        self.ode_instance = ode_instance
        self.last = None
    
    def terminate(self, t, v, dt):
        
        v_norm = norm(v)
        
        if v_norm > 1E10:
            print 'Blew up!'
            return dt, False, True
        
        if self.last is not None:
            e = norm(self.last-v)/(dt*v_norm)
            if e < 0.1:
                return dt, False, True
            print t, '\t', v_norm, '\t', e
        
        self.last = v.copy()
        
        return dt, True, False


Del2 = hdafLaplacian2D()
Del2.set(grid=grid, m1=8,m2=8,gamma1=0.8,gamma2=0.8 )

def rhs(t,v):
    r = grid.zeros()
    Del2(v,r)
    r *= chi_D
    return r


heat = ode()
heat.solver = RK4solver()
heat.step_size = 0.002
heat.terminator = terminator(heat)
heat.rhs = rhs

t,u = heat.run((0.0, F))

#u *= chi_D

#write_png( u, "u.png", major_scale=std(u)*1, center=mean(u) )

#a = u - F
#write_png( a, "out.png", major_scale=std(a)*1, center=mean(a) )





#-----------
"""

from dirichlet_laplace import *

x = (-1.0, -1.0)

x = grid.discretize(x)
i,j = grid.index(x)
n_samples=500
(w,sig_w) = solution( x=x, f=lambda x: f(x[0],x[1]), n_samples=n_samples, dt=0.02, drift=lambda t:array((0.0, 0.0)), region=lambda x: chi_D(x)>0.5 )
print u[i][j], w, sig_w


def draw(u, fname, s=1.0):
    write_png( u, fname, major_scale=std(u)*s, center=mean(u) )

"""



