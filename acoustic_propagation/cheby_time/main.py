
from setup import *


region1 = polygon2D([ (-2.0,-2.0), (-4.0,8.0), (-3.5,7.0), (0.0,0.0)  ])
region2 = circle( array((0.0,0.0)), 2.0 )
back = complement(union(region1, region2))

velocity = simpleFunction([ (region1,20.0**2), (region2,15.0**2), (back,10.0**2) ])

hdaf_order=24
hdaf_gamma=0.8
grid = grid2d( a1=-10.0, b1=10.0, a2=-10.0, b2=10.0, n1=256, n2=256 )

P = propagator1()

velocity = grid.evaluate( lambda x,y: velocity((x,y)) )

write_png( velocity, "vel.png", center=mean(velocity), major_scale=std(velocity), red_params=(0.33,0,0,1), green_params=(0.5,0,1,0), blue_params=(0.66,1,0,0,) )


print 'initializing propagator'
P.set( grid, velocity=velocity, damping=None,  expansion_order=2, hdaf_order=hdaf_order, hdaf_gamma=hdaf_gamma )

x0 = array((7.,7.))
u0 = lambda x,y: exp( -norm( array((x,y)) - x0 )**2/(2.0*(0.3)**2 )) 

u = grid.evaluate( u0 )
v = grid.zeros()

fname = temp_directory + str(time_module.time())
scipy.io.savemat( fname, {'m':u} )
u_fnames = [fname ]

time_step = 0.002
stop_time = 3.0

T0 = wall_time()
"""
step=0
time = 0.0
while time<=stop_time:
    
    print time
    
    u_norm = norm(u)*grid.dx1*grid.dx2
    
    if u_norm > 1E4:
        print "Blew up!"
        quit()
    
    err = P( time_step, u,v, u,v, )
        
    if err is not None:
        print err
        break
    
    step += 1
    time = time_step*step
    
    fname = temp_directory + str(time_module.time())
    scipy.io.savemat( fname, {'m':u} )
    u_fnames.append(fname)

T_full = wall_time()-T0

print "full P time", T_full
"""


P = propagator2()
P.set( grid, velocity=velocity, damping=None,  expansion_order=2, hdaf_order=hdaf_order, hdaf_gamma=hdaf_gamma )

u1 = grid.evaluate( u0 )
u2 = grid.zeros()
err = P( time_step, u1, u2 )
if err is not None:
    print err
    quit()
u3 = grid.zeros()

fname = temp_directory + str(time_module.time())
scipy.io.savemat( fname, {'m':u1} )
u_fnames = [fname ]

fname = temp_directory + str(time_module.time())
scipy.io.savemat( fname, {'m':u2} )
u_fnames.append( fname )

T0 = wall_time()

step=1
time = time_step*step

while time<=stop_time:
    
    print time
    
    u_norm = norm(u1)*grid.dx1*grid.dx2
    
    if u_norm > 1E4:
        print "Blew up!"
        quit()
    
    err = P( time_step, u2, u3 )
    
    if err is not None:
        print err
        break
    
    u3 = u3*2.0-u1
    u1 = u2.copy()
    u2 = u3.copy()
    
    step += 1
    time = time_step*step
    
    fname = temp_directory + str(time_module.time())
    scipy.io.savemat( fname, {'m':u2} )
    u_fnames.append(fname)

T_cheby = wall_time()-T0


print "cheby P time", T_cheby


center = 0.0
major_scale = 0.2
directory="/workspace/output/acoustic_propagation/"

for step,fname in enumerate(u_fnames):
    print step
    data = scipy.io.loadmat( fname )
    u = data['m'].copy('C')
    
    write_png( u, directory+"%05d.png"%step, center=center, major_scale=major_scale, red_params=(0.33,0,0,1), green_params=(0.5,0,1,0), blue_params=(0.66,1,0,0,) )












