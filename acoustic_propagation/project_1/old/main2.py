
from setup import *

dir = "/tmp/research/"

region_1 = polygon2D([ (0,375), (1500,375), (1500,438), (0,438) ])
region_2 = polygon2D([ (1500,438), (3000,373), (3000,438), (1500,495) ])
region = union( region_1, region_2 )
velocity = simpleFunction([ (region,4500.0**2), (complement(region),2000.0**2) ])


def damping(x,y, grid):
    scale = 50.0
    return 10000*( cosh(abs(x-grid.a1)/scale)**-2 + cosh(abs(x-grid.b1)/scale)**-2 + cosh(abs(y-grid.a2)/scale)**-2 + cosh(abs(y-grid.b2)/scale)**-2  )

def ricker_wavelet(t, sigma, gamma, tau ):
    return -sqrt(2.0/pi)*sigma*gamma*(sigma-2.0*sigma*gamma*(sigma*t-tau)**2)*exp(-gamma*(sigma*t-tau)**2)


session = {}
session['grid'] = grid2d( a1 = -500, b1 = 3500, a2=-300, b2=1200, n1=1024, n2=384)
session['damp_index1'] = session['grid'].index1( 1500 )
session['damp_index2'] = session['grid'].index2( 800 )
session['driving_function'] = lambda t: ricker_wavelet(t=t, sigma=1.5*20, gamma=8, tau=1 )
#lambda t: float( t<=1.2*(session['time_step']) )
session['velocity_function'] = lambda x,y: velocity((x,y))
session['damping_function'] = damping
session['velocity'] = array([])
session['damping'] = array([])
session['driving'] = []
session['expansion_order'] = 2
session['m1'] = 8
session['m2'] = 8
session['gamma1'] = 0.8
session['gamma2'] = 0.8
session['time_step'] = 0.001
session['final_time'] = 0.7
session['u_fnames'] = []
session['time'] = []

P = propagator1()

def init(session, P ):
    print 'initializing, evaluating velocity'
    session['velocity'] = session['grid'].evaluate( session['velocity_function'] )
    print 'initializing, evaluating damping'
    session['damping'] = session['grid'].evaluate( lambda x,y: (session['damping_function'])(x,y,session['grid']) )
    print 'initializing propagator'
    P.set( session['grid'], velocity=session['velocity'], damping=session['damping'],  expansion_order=session['expansion_order'],   hdaf_m1=session['m1'], hdaf_m2=session['m2'], hdaf_gamma1=session['gamma1'], hdaf_gamma2=session['gamma2'] )


def run( session ):
    
    grid = session['grid']
    
    u = grid.zeros()
    v = grid.zeros()
    driving_1 = grid.zeros()
    driving_2 = grid.zeros()
    i0 = session['damp_index1']
    j0 = session['damp_index2']
    fname = dir  + str(time_module.time())
    scipy.io.savemat( fname, {'m':u} )
    #results[str(step)] = { 'time': time, 'step':step, 'fname':fname }
    u_fnames = [fname]
    driving = []
    time_step = session['time_step']
    final_time = session['final_time']
    step=0
    time = [step]
    
    while time[step] <= final_time:
        
        print "step:", step, "\ttime:", time[step], "\twall time", wall_time(), "\tmagnitude:", norm(u)
        
        driving_1[i0,j0] = (session['driving_function'])(time[step])
        driving_2[i0,j0] = (session['driving_function'])(time[step]+time_step)
        driving.append(driving_1[i0,j0])
        
        P( time_step, u,v, u,v, driving_1, driving_2 )
        
        fname = dir  + str(time_module.time())
        scipy.io.savemat( fname, {'m':u} )
        u_fnames.append(fname)
        
        #results[str(step)] = { 'data':'wave', 'time': time, 'step':step, 'fname':fname }
        #results['final_step'] = step
        
        step += 1
        time.append(time_step*step)
    
    session['u_fnames'] = u_fnames
    session['driving'] = driving
    session['time'] = time


def FT( session, frequency ):
    
    grid = session['grid']
    time_step = session['time_step']
    FT_sum_real = grid.zeros()  
    FT_sum_imag = grid.zeros()
    
    for step,time in enumerate(session['time']):
        print time
        data = scipy.io.loadmat( session['u_fnames'][step] )
        u = data['m']
        
        FT_sum_real += u*cos(-2.0*pi*time*frequency)*time_step
        FT_sum_imag += u*sin(-2.0*pi*time*frequency)*time_step
    
    return (FT_sum_real, FT_sum_imag)


    

def save( session, fname ):
    del(session['velocity_function'])
    del(session['damping_function'])
    del(session['driving_function'])
    output = open(fname, 'wb')
    pickle.dump(session, output)
    output.close()


def load( fname ):
    input = open(fname, 'rb')
    session = pickle.load(input)
    input.close()
    return session

def draw( data, fname ):
    write_png( data, fname, center=mean(data), major_scale=std(data)*3, red_params=(0.33,0,0,1), green_params=(0.5,0,1,0), blue_params=(0.66,1,0,0,) )


def render_all( session, dir="/workspace/output/scratch/" ):
    
    for step,time in enumerate(session['time']):
        print time
        data = scipy.io.loadmat( session['u_fnames'][step] )
        u = data['m']
        
        draw( u, dir+"%04d.png"%step )



if __name__ == "__main__":
    
    print "running main"
    
    init(session, P)
    run(session)
    save(session, 'session')




















"""



(ft_real, ft_imag) = FT( session, 10)
phase = numpy.vectorize( lambda x,y: atan2(x,y) )
A = sqrt( ft_real*ft_real + ft_imag*ft_imag)
phi = phase(ft_real, ft_imag)

dir = "/workspace/output/acoustic_propagation/"
draw(A, dir+ "A.png" )
draw(ft_real, dir+ "ft_real.png" )
draw(ft_imag, dir+ "ft_imag.png" )
write_png( phi, dir+ "phase.png", center=0.0, major_scale=1.5, red_params=(0.33,0,0,1), green_params=(0.5,0,1,0), blue_params=(0.66,1,0,0,) )









(FTr, FTi) = ft = FT(results, 47.1)
write_png( FTr, dir+ "FT_real.png", center=mean(FTr), major_scale=std(FTr)*3, red_params=(0.33,0,0,1), green_params=(0.5,0,1,0), blue_params=(0.66,1,0,0,) )
write_png( FTi, dir+ "FT_imag.png", center=mean(FTi), major_scale=std(FTi)*3, red_params=(0.33,0,0,1), green_params=(0.5,0,1,0), blue_params=(0.66,1,0,0,) )



def FT( results ):
    
    
    

#write_png( velocity_grid, dir+ "vlocity.png", center=mean(velocity_grid), major_scale=std(velocity_grid), red_params=(0.33,0,0,1), green_params=(0.5,0,1,0), blue_params=(0.66,1,0,0,) )
#write_png( damping_grid, dir+ "damping.png", major_scale=1.0, red_params=(0.33,0,0,1), green_params=(0.5,0,1,0), blue_params=(0.66,1,0,0,) )




    
#write_png( u, dir+ "%04d.png"%count, center=mean(u), major_scale=std(u)*3, red_params=(0.33,0,0,1), green_params=(0.5,0,1,0), blue_params=(0.66,1,0,0,) )





FT_sum = grid.zeros()
kernel = lambda t,f: exp(-2.0j*pi*t*f)
frequency = 47.110022310144111  
K = grid.ones()*kernel(time,frequency)
FT_sum += u*K*time_step

FT_sum_real = real(FT_sum)
FT_sum_imag = imag(FT_sum)


"""
"""
T = 0.7*array([ float(i)/1000.0 for i in range(1000) ])
F = [ driving_function_time(t) for t in T ]

import pylab as p
p.plot( T, F)
p.show()
"""
