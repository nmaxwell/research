
from setup import *

temp_directory = "/workspace/tmp/research/"


region_1 = polygon2D([ (0,375), (1500,375), (1500,438), (0,438) ])
region_2 = polygon2D([ (1500,438), (3000,373), (3000,438), (1500,495) ])
region = union( region_1, region_2 )
velocity = simpleFunction([ (region,4500.0**2), (complement(region),2000.0**2) ])


def ricker_wavelet(t, sigma, gamma, tau ):
    return -sqrt(2.0/pi)*sigma*gamma*(sigma-2.0*sigma*gamma*(sigma*t-tau)**2)*exp(-gamma*(sigma*t-tau)**2)


grid = grid2d( a1 = -500, b1 = 3500, a2=-300, b2=1200, n1=1024/2, n2=384/2)

session = {}
session['grid'] = grid
session['drive_index1'] = session['grid'].index1( 1500 )
session['drive_index2'] = session['grid'].index2( 800 )
session['driving_function'] = lambda t: ricker_wavelet(t=t, sigma=1.5*20, gamma=8, tau=1 )
session['velocity_function'] = lambda x,y: velocity((x,y))
session['velocity'] = array([])
session['damping'] = cosh_damping( grid, x1=grid.a1, x2=grid.b1, y1=grid.a2, y2=grid.b2, xs=50.0, ys=50.0, amp=1.0E4 )
session['driving'] = []
session['expansion_order'] = 8
session['hdaf_order'] = 8
session['hdaf_gamma'] = 0.8
session['time_step'] = 0.002
session['final_time'] = 0.7
session['u_fnames'] = []
session['time'] = []
session['skip_save']=int(0.01/session['time_step'])

P = propagator1()

def init(session, P ):
    print 'initializing, evaluating velocity'
    session['velocity'] = session['grid'].evaluate( session['velocity_function'] )
    print 'initializing propagator'
    P.set( session['grid'], velocity=session['velocity'], damping=session['damping'],  expansion_order=session['expansion_order'],   hdaf_order=session['hdaf_order'], hdaf_gamma=session['hdaf_gamma'] )


def run( session ):
    
    grid = session['grid']
    
    u = grid.zeros()
    v = grid.zeros()
    driving_1 = grid.zeros()
    driving_2 = grid.zeros()
    i0 = session['drive_index1']
    j0 = session['drive_index2']
    fname = temp_directory  + str(time_module.time())
    scipy.io.savemat( fname, {'m':u} )
    #results[str(step)] = { 'time': time, 'step':step, 'fname':fname }
    u_fnames = [fname]
    driving = []
    time_step = session['time_step']
    final_time = session['final_time']
    
    step=0
    save_times = []
    time = 0.0
    
    start_wall_time = wall_time()
    
    while time <= final_time:
        
        u_norm = norm(u)*grid.dx1*grid.dx2
        
        if u_norm > 1E6:
            print "Blew up!"
            break
        
        try:
            rate = ((wall_time()-start_wall_time)/time)
            eta = rate*(session['final_time']-time)
            print "step:", step, "\ttime:", time, "\twall time", wall_time(), "\teta:", eta ,"\tmagnitude:", u_norm, "\trate:", rate
        except:
            print "step:", step, "\ttime:", time, "\twall time", wall_time()
        
        driving_1[i0,j0] = (session['driving_function'])(time)
        driving_2[i0,j0] = (session['driving_function'])(time+time_step)
        driving.append(driving_1[i0,j0])
        
        err = P( time_step, u,v, u,v, driving_1, driving_2 )
        
        if err is not None:
            print err
            break
        
        if step%session['skip_save'] == 0:
            fname = temp_directory  + str(time_module.time())
            scipy.io.savemat( fname, {'m':u} )
            u_fnames.append(fname)
            save_times.append(time)
            
            session['u_fnames'] = u_fnames
            session['time'] = save_times
        
        step += 1
        time = time_step*step
    

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
    try:
        del(session['velocity_function'])
    except:
        pass
    try:
        del(session['driving_function'])
    except:
        pass
    output = open(fname, 'wb')
    pickle.dump(session, output)
    output.close()


def load( fname ):
    input = open(fname, 'rb')
    session = pickle.load(input)
    input.close()
    return session

def render( data, fname ):
    ar = data.copy()
    write_png( ar, fname, center=mean(ar), major_scale=std(ar)*3, red_params=(0.33,0,0,1), green_params=(0.5,0,1,0), blue_params=(0.66,1,0,0,) )

def render_all( session, directory="/workspace/output/acoustic_propagation/",  skip=1, new_grid=None ):
    
    resampler=None
    
    if new_grid is not None:
        grid_in=session['grid']
        grid_out = new_grid
        resampler = linearResampler2d( in_grid=grid_in, out_grid=new_grid )
    
    skipped=0
    
    for step,time in enumerate(session['time']):
        if skipped%skip == 0:
            print time
            data = scipy.io.loadmat( session['u_fnames'][step] )
            u = None
            if resampler is not None:
                u = resampler(data['m'].copy('C'))
            else:
                u = data['m'].copy('C')
            
            render( u, directory+"%05d.png"%step )
        skipped += 1

if __name__ == "__main__":
    
    #grid=grid2d( a1=grid_in.a1, a2=grid_in.a2,  b1=grid_in.b1, b2=grid_in.b2, n1=grid_in.n1/2, n2=grid_in.n2/2 )
    
    
    print "running main"
    
    init(session, P)
    render( session['velocity'], "velocity.png" )
    run(session)
    save(session, 'session')
    
    render_all(session)







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
