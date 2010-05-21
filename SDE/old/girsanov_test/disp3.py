
from pylab import *
sys.path.append("/workspace/mathlib/tools/python/")
from read_file import *

dir = "/workspace/output/SDE/test/"



time = read_file( dir + "time" )
Xm = read_file( dir + "X_mean" )

plot( time, Xm )
show()

