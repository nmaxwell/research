
import os
from pylab import *
sys.path.append("/workspace/mathlib/tools/python/")
from read_file import *

dir = "/workspace/output/SDE/test/"

def rms(a):
    return sqrt(sum(a*a)/len(a))


def run_sde(n_runs = 100, n_steps=100, X0=0.0, stop_time=1.0, drift=lambda t: 0.0 ):
    
    time = (arange(n_steps)/float(n_steps))*stop_time
    
    drift = array([ drift(t) for t in time ])
    
    fl = open(dir + 'drift', 'w+' )
    for k in range(len(drift)):
        fl.write(str(drift[k]) + '\n' )
    fl.close()
    
    os.system( "./a.out " + str(n_runs) + " " + str(n_steps) + " " + str(X0) + " " + str(stop_time) + " " + "drift" )
    
    X = array(read_file( dir + "X" ))
    M = array(read_file( dir + "M" ))
    PX = array(read_file( dir + "PX" ))
    QX = array(read_file( dir + "QX" ))
    VX = array(read_file( dir + "VX" ))
    
    return { 'time': time, 'X': X, 'M': M, 'drift': drift, 'PX': PX, 'QX': QX, 'VX': VX }







"""

n_runs = len(X)
n_steps = len(X[0])

time = array(time)
w_slices= array(X)
t_slices = w_slices.transpose()

x_vars = array([ var(t_slices[k]) for k in range(n_steps) ])


err = (x_vars-time)/time
err[0] = 0

print "variance of sample paths, err = (x_var-time)/time"
print "var(err)= ", var(err),"\nmean(err)= ",mean(err), "\nrms(err)= ",rms(err)






"""


