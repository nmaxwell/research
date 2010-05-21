
from pylab import *
sys.path.append("/workspace/mathlib/tools/python/")
from read_file import *

dir = "/workspace/output/SDE/project1/hitting_distribution/"


dist = read_file( dir + "distribution" )[0]
dist_nodrift = read_file( dir + "distribution_nodrift" )[0]
dist_nodrift_nochange = read_file( dir + "distribution_drift_nochange" )[0]


plot(dist, 'ro' )
plot(dist_nodrift, 'bo' )
plot(dist_nodrift_nochange, 'go' )

show()

