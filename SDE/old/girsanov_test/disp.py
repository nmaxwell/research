
from pylab import *
sys.path.append("/workspace/mathlib/tools/python/")
from read_file import *

dir = "/workspace/output/SDE/test/"



X = read_file( dir + "X" )
M = read_file( dir + "M" )
time = read_file( dir + "time" )
Xm = read_file( dir + "X_mean" )
Y = [ t*t/2 for t in time ]

XM = read_file( dir + "XM" )
QX = read_file( dir + "QX" )

"""
for k in range(4):
    plot( time, X[k])

plot( time, Xm)
plot( time, Y)
show()
"""

for k in range(1):
    plot( time, X[k])
plot( time, QX)
plot( time, Xm)
show()

