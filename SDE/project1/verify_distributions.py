
from dirichlet_laplace import *
from scipy import stats


class interval:
    "interval (a,b] of the real numbers."
    
    def __init__(self, a,b ):
        self.a = a
        self.b = b
    
    def __contains__(self, x):
        return self.a < x and x <= self.b



def check_distribution_nball_quadrant():
    
    dim = 6
    
    measure_sets = [ (k,) for k in range(2**dim) ]
    region = lambda x: l2norm(x)<=1.0
    f = lambda x : sum([(2.0**k)*a for k,a in enumerate([xi>=0 for xi in x])])
    
    dt = 0.03
    n_samples = 20000
    x = zeros(dim)
    drift = lambda t: array((0.,0.,.2,0.,-0.9,0.))
    
    distribution, distribution_nocom = hitting_value_distribuion( n_samples, measure_sets, x, f, dt, drift, region )
    
    drift = lambda t: zeros(dim)
    distribution_nodrift, blah = hitting_value_distribuion( n_samples, measure_sets, x, f, dt, drift, region )
    
    import matplotlib.pyplot as plt
    left = range(2**dim+1)
    
    print len(left), len(distribution_nodrift)
    
    distribution_nodrift.append(distribution_nodrift[-1])
    distribution.append(distribution[-1])
    distribution_nocom.append(distribution_nocom[-1])
    plt.plot(left, distribution_nodrift, 'g', drawstyle='steps-post'  )
    plt.plot(left, distribution, 'b', drawstyle='steps-post' )
    plt.plot(left, distribution_nocom, 'r',  drawstyle='steps-post' )
    plt.legend( ('a=0','w.r.t. Q', 'w.r.t. P') )
    plt.xlim( (0,2**dim) )
    plt.ylabel("probability")
    plt.xlabel("hitting position, quadrant")
    plt.title("hitting distributions")
    
    plt.savefig('figure6.eps')





def check_distribution_R2_unitball():
    
    f = lambda x: math.atan2(x[0],x[1])
    region=lambda x: l2norm(x)<=1.0
    
    n_samples = 1000
    dt = 0.002
    n_measure_sets = 32
    
    measure_sets = [ interval( 2.*pi*k/n_measure_sets-pi, 2.*pi*(k+1)/n_measure_sets-pi ) for k in range(n_measure_sets) ]
    
    x = (0.2, 0.0 )
    
    drift=lambda t: array(( 0.0,0.0 ))
    
    distribution_nodrift, blah = hitting_value_distribuion( n_samples, measure_sets, x, f, dt, drift, region )
    
    drift=lambda t: array(( -3.0*sin(t*pi*3),cos(t*pi*3) ))*(t+1)
    
    distribution, distribution_nocom = hitting_value_distribuion( n_samples, measure_sets, x, f, dt, drift, region )
    
    print "sums:"
    print sum(distribution), sum(distribution_nocom), sum(distribution_nodrift)
    
    import matplotlib.pyplot as plt
    left = [  (2.*k/n_measure_sets-1)*pi  for k in range(n_measure_sets+1) ]
    
    distribution_nodrift.append(distribution_nodrift[-1])
    distribution.append(distribution[-1])
    distribution_nocom.append(distribution_nocom[-1])
    plt.plot(left, distribution_nodrift, 'g', drawstyle='steps-post'  )
    plt.plot(left, distribution, 'b', drawstyle='steps-post' )
    plt.plot(left, distribution_nocom, 'r',  drawstyle='steps-post' )
    plt.legend( ('a=0','w.r.t. Q', 'w.r.t. P') )
    
    import fractions
    
    def pi_labeler(x):
        if float(x) in [-1., 0., 1.]:
            x = int(x)
            if x == -1:
                return '-' + u"\u03C0"
            if x == 0:
                return '0'
            if x == 1:
                return u"\u03C0"
        else:
            return str(x) + u"\u03C0"
    
    n_ticks = 8
    ticks = [ (2.0*k/n_ticks-1)*pi for k in range(n_ticks+1) ]
    labels = [ pi_labeler(x) for x in [ fractions.Fraction(2*k-n_ticks, n_ticks) for k in range(n_ticks+1) ] ]
    plt.xticks( ticks, labels )
    plt.xlim( (-3.3,3.3) )
    plt.ylabel("probability")
    plt.xlabel("hitting position")
    plt.title("hitting distributions")
    plt.savefig('figure3.eps')
    
    
    
    
    
    
    
    
    
check_distribution_nball_quadrant()

#check_distribution_R2_unitball()




