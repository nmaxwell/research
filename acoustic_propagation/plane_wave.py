
import sys
sys.path.append( "/workspace/research/acoustic_propagation" )
from math import *  
from solve import solve1 as solve1
from solve import solve2 as solve2
import pickle
import scipy.io
import waveprop as wp
import numpy
norm = numpy.linalg.norm

def render( data, fname, major_scale=None, center=None ):
    if major_scale is None:
        major_scale = std(data)*10.0
    if center is None:
        center=mean(data)
    wp.write_png( data, fname, center=center, major_scale=major_scale, red_params=(0.33,0,0,1), green_params=(0.5,0,1,0), blue_params=(0.66,1,0,0,) )


grid = wp.grid2d( a1 = 0.0, b1 = 1, a2=0, b2=1, n1=256, n2=256  )

c=2.0
velocity = grid.evaluate( lambda x,y: c**2 )

k1 = 2.0*2.0*pi
k2 = -3.0*2.0*pi
u0 = lambda x1,x2: cos(x1*k1+x2*k2)

time_step=0.002
results1 = solve2(grid, initial_wave=u0, time_step=time_step, final_time=1.0, damping=None, velocity=velocity, expansion_order=2 )
#results2 = solve2(grid, initial_wave=u0, time_step=time_step, final_time=1.0, damping=None, velocity=velocity )


for step in range(len(results1['files'])):
    data = scipy.io.loadmat( results1['files'][step] )
    u1 = data['m'].copy('C')
    #data = scipy.io.loadmat( results2['files'][step] )
    #u2 = data['m'].copy('C')
    
    k=sqrt(k1*k1+k2*k2)
    t = float(step)*time_step
    u2 = grid.evaluate( lambda x1,x2: cos(x1*k1+x2*k2)*cos(c*t*k) )
    
    render( u1-u2, "/workspace/output/acoustic_propagation/" + "%05d.png"%step , major_scale=1E-1, center=0.0 )
    print step, log10( norm( u1-u2 )/norm(u2) + 1E-17)
    
