
from pylab import *
sys.path.append("/workspace/mathlib/tools/python/")
from read_file import *

dir = "/workspace/output/SDE/test/"



time = read_file( dir + "time" )

QX = read_file( dir + "QX" )

plot( time, QX)
show()

