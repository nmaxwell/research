
import sys
sys.path.append( "/workspace/research/acoustic_propagation" )
from math import *  
from solve import solve1 as solve
import waveprop as wp
from setup import *


def render( data, fname, major_scale=None, center=None ):
    ar = data.copy()
    if major_scale is None:
        major_scale = std(ar)*10.0
    if center is not None:
        center=mean(ar)
    wp.write_png_720p( ar, fname, center=center, major_scale=major_scale, red_params=(0.33,0,0,1), green_params=(0.5,0,1,0), blue_params=(0.66,1,0,0,) )


grid = wp.grid2d( a1 = -500, b1 = 3500, a2=-300, b2=1200, n1=1024, n2=384)

"""
region_1 = polygon2D([ (0,375), (1500,375), (1500,438), (0,438) ])
region_2 = polygon2D([ (1500,438), (3000,373), (3000,438), (1500,495) ])
region = union( region_1, region_2 )
velocity = simpleFunction([ (region,4500.0**2), (complement(region),2000.0**2) ])
"""

c=2000.0
velocity = grid.evaluate( lambda x,y: c**2 )
damping = grid.zeros()
#wp.cosh_damping( grid, x1=grid.a1, x2=grid.b1, y1=grid.a2, y2=grid.b2, xs=50.0, ys=50.0, amp=1.0E4 )


k1 = 0.01*2.0*pi
k2 = 0.01*2.0*pi

u0 = lambda x1,x2: cos(x1*k1+x2*k2) 
#*(1.0+(1.15*(x1-1500.0)/2000.0)**48)**-1 *(1.0+(1.4*(x2-450)/750.0)**48)**-1

time_step=0.001

u0 = grid.evaluate( u0 )
u0 = u0 * exp( -damping * time_step )

results = solve(grid, initial_wave=u0, time_step=time_step, final_time=  1.0, damping=damping, velocity=velocity )


for step in range(len(results['files'])):
    data = scipy.io.loadmat( results['files'][step] )
    u1 = data['m'].copy('C')
    #data = scipy.io.loadmat( results2['files'][step] )
    #u2 = data['m'].copy('C')
    
    k=sqrt(k1*k1+k2*k2)
    t = float(step)*time_step
    u2 = grid.evaluate( lambda x1,x2: cos(x1*k1+x2*k2)*cos(c*t*k) )
    
    render( u1-u2, "/workspace/output/acoustic_propagation/" + "%05d.png"%step , major_scale=1E-1, center=0.0 )
    print step, norm( u1-u2 )/norm(u2)



"""
for step, fname in enumerate(results['files']): 
    print step, "/", len(results['files'])
    data = scipy.io.loadmat( fname )
    data = data['m'].copy('C')
    render( data, "/workspace/output/acoustic_propagation/" + "%05d.png"%step , major_scale=1.0, center=0.0 )


output = open('results', 'wb')
pickle.dump(results, output)
output.close()


def load( fname='results' ):
    input = open(fname, 'rb')
    results = pickle.load(input)
    input.close()
    return results

"""









"""
freq = 47.110022310144

FT( session, freq )
"""




"""
render( session['velocity'], "velocity.png" )
render_all(session)
"""













