
import waveprop as wp
from numpy import *
import time as time_module
from math import *
import scipy.io
import numpy
norm = numpy.linalg.norm

def solve1( grid, out_directory="/workspace/tmp/research/", initial_wave=None, initial_time=0.0, time_step=0.01, final_time=1.0, velocity=lambda x,y:1.0, damping=None,  expansion_order=2, hdaf_order=12, hdaf_gamma=0.8 ):
    
    results={}
    
    P = wp.propagator2()
    P.set( grid, velocity=velocity, damping=damping,  expansion_order=expansion_order, hdaf_order=hdaf_order, hdaf_gamma=hdaf_gamma )
    
    try:
        u1 = grid.evaluate( initial_wave )
    except:
        try:
            u1 = array(initial_wave)
        except:
            u1=grid.zeros()
    
    u2 = grid.zeros()
    err = P( time_step, u1, u2 )
    if err is not None:
        print "solve1, error: ", err
        return results
    
    results["files"]=[]    
    fname = out_directory + str(time_module.time())
    scipy.io.savemat( fname, {'m':u1} )
    results["files"].append( fname ) 
    fname = out_directory + str(time_module.time())
    scipy.io.savemat( fname, {'m':u2} )
    results["files"].append( fname )
    
    step = 1
    time = initial_time+step*time_step
    temp = grid.zeros()
    
    while time<=final_time:
        
        err = P( time_step, u2, temp )
        if err is not None:
            print "solve1, error: ", err
            return results
        
        temp = temp*2.0-u1
        u1 = u2
        u2 = temp
        
        step += 1
        time = initial_time+step*time_step
        
        u_norm = norm(u2)*grid.dx1*grid.dx2
        if u_norm >= 1E6:
            print "solve1, blew up!"
            return results
        
        print "step: ", "%05d"%step, "time: ", "%6.3f"%time, "norm: ", "%03.2e"%u_norm
        
        fname = out_directory + str(time_module.time())
        scipy.io.savemat( fname, {'m':u2} )
        results["files"].append( fname ) 
        
    print "success!"
    return results


    

def solve2( grid, out_directory="/workspace/tmp/research/", initial_wave=None, initial_time=0.0, time_step=0.01, final_time=1.0, velocity=lambda x,y:1.0, damping=None, driving=None, expansion_order=2, hdaf_order=12, hdaf_gamma=0.8 ):
    
    results={}
    
    P = wp.propagator1()
    P.set( grid, velocity=velocity, damping=damping,  expansion_order=expansion_order, hdaf_order=hdaf_order, hdaf_gamma=hdaf_gamma )
    
    u=None
    try:
        u = grid.evaluate( initial_wave )
    except:
        try:
            u = array(initial_wave)
        except:
            u=grid.zeros()
    v=grid.zeros()
    
    results["files"]=[]    
    fname = out_directory + str(time_module.time())
    scipy.io.savemat( fname, {'m':u} )
    results["files"].append( fname )
    
    if driving is not None:
        pass
    
    step = 0
    time = initial_time+step*time_step
    temp = grid.zeros()
    
    while time<=final_time:
        
        err = P( time_step, u,v, u,v )
        if err is not None:
            print "solve2, error: ", err
            return results
        
        step += 1
        time = initial_time+step*time_step
        
        uv_norm = (norm(u)+norm(v))*grid.dx1*grid.dx2
        if uv_norm >= 1E6:
            print "solve2, blew up!"
            return results
        
        print "step: ", "%05d"%step, "time: ", "%6.3f"%time, "norm: ", "%03.2e"%uv_norm
        
        fname = out_directory + str(time_module.time())
        scipy.io.savemat( fname, {'m':u} )
        results["files"].append( fname ) 
        
    print "success!"
    return results
    












