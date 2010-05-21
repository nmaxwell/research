
from setup import *

temp_dir = "/workspace/tmp/research/"


grid = grid2d( a1 = -2, b1 = 2, a2=-1, dx1=.03125, dx2=0.03125, n2=512 )
grid.debug()

region = polygon2D([ (-1.0,0.0), (-1.0,1E3), (1.0,1E3),  (1.0,0.0) ])

velocity = simpleFunction([ (region,100.0**2), (complement(region),200.0**2) ])

model = grid.evaluate( lambda x,y: velocity((x,y)) )


session = {}
session['grid'] = grid
session['damping'] = cosh_damping( grid, x1=grid.a1, x2=grid.b1, y1=grid.a2, y2=grid.b2, xs=0.1, ys=0.1, amp=1.0E2 )
session['initial_condition'] = lambda x,y: 0.0
#exp(-norm( array((x,y)) - array((0.0, 0.0)) )**2/(2.0*0.03**2))
session['velocity'] = model
session['driving'] = []
session['expansion_order'] = 2
session['hdaf_order'] = 8
session['hdaf_gamma'] = 0.80
session['time_step'] = 0.0001
session['final_time'] = 1.0
session['u_fnames'] = []
session['time'] = None
session['skip_save']=0.002/session['time_step']

P = propagator1()


write_png(session['damping'], "damp.png", major_scale=std(session['damping'])*3, center=mean(session['damping']), ordering='rm' )
write_png(model, "model.png", major_scale=std(model)*3, center=mean(model), ordering='rm' )

driving_function = lambda t: 100.0*sin(6.283*t*120.0)
i0 = grid.index1(0.0)
j0 = grid.index2(0.0)

"""
import pylab as p
Y = grid.X2[256]
Z = model[256]
damp = session['damping'][256]

p.plot(Y,damp)
p.plot(Y,Z)
p.show()
quit()
"""

"""
import pylab as p
X = transpose(grid.X1)[128]
Z = transpose(model)[128]
damp = transpose(session['damping'])[128]

p.plot(X,damp)
p.plot(X,Z)
p.show()
quit()
"""



def init(session, P ):
    print 'initializing propagator'
    P.set( session['grid'], velocity=session['velocity'], damping=session['damping'],  expansion_order=session['expansion_order'],   hdaf_order=session['hdaf_order'], hdaf_gamma=session['hdaf_gamma'] )
    



def run( session ):
    
    grid = session['grid']
    
    u = grid.evaluate( session['initial_condition'] )
    v = grid.zeros()
    
    fname = temp_dir  + str(time_module.time())
    scipy.io.savemat( fname, {'m':u} )
    u_fnames = [fname]
    
    time_step = session['time_step']
    final_time = session['final_time']
    step=0
    save_times = []
    time = 0.0
    
    driving_1 = grid.zeros()
    driving_2 = grid.zeros()
    
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
        
        
        driving_1[i0][j0] = driving_function(time)
        driving_2[i0][j0] = driving_function(time+time_step)    
        
        err = P( time_step, u,v, u,v, driving_1, driving_2 )
        
        if err is not None:
            print err
            break
        
        if step%session['skip_save'] == 0:
            fname = temp_dir  + str(time_module.time())
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
    del(session['initial_condition'])
    output = open(fname, 'wb')
    pickle.dump(session, output)
    output.close()


def load( fname ):
    input = open(fname, 'rb')
    session = pickle.load(input)
    input.close()
    return session

def render( data, fname, ordering='rm' ):
    ar = data.copy()
    write_png( ar, fname, center=mean(ar), major_scale=std(ar)*3, red_params=(0.33,0,0,1), green_params=(0.5,0,1,0), blue_params=(0.66,1,0,0,), ordering=ordering )

def render_all( session, directory="/workspace/output/acoustic_propagation/",  skip=0 ):
    
    skipped=0
    
    for step,time in enumerate(session['time']):
        if skipped%skip == 0:
            print time
            data = scipy.io.loadmat( session['u_fnames'][step] )
            u = data['m']
            
            render( u, directory+"%05d.png"%step )
        skipped += 1

if __name__ == "__main__":
    
    print "running main"
    
    init(session, P)
    run(session)
    save(session, 'session')
    
    render_all(session, skip=1)






