

//#define ML_NOT_USE_FFTW_MALLOC
 
#include <mathlib/math/random/ml_random.h>
#include <mathlib/math/SDE/BM.h>
#include <mathlib/math/grids/grid1D.h>
#include <mathlib/math/grids/grid2D.h>
//#include <mathlib/math/grids/extra/plots.cpp>
 
#include <mathlib/link.cpp>
#include <mathlib/non-lib_things.h>
 
//#include <boost/math/special_functions/erf.hpp>
 
 
 
 
/*
 
 
 
 
import math
import pylab as p
import sys
sys.path.append("/workspace/mathlib/tools/python/")
from read_file import *
 
for n in range(1,10):
p.plot( read_file( "out_" + str(n) ) )
 
p.plot ( read_file( "U_mean" ) )
p.show()
 
Y = read_file( "U_mean" )
X = [ float(k) / len(Y) for k in range(len(Y))]
Y = [ math.exp(x) for x in X ]
p.plot(X,Y,X, read_file( "U_mean" ) )
p.show()
 
 
 
 
 
*/
 
 
 
int main()
{
    /*
* Ito and Stratonovich integrals of W dW
*
*
*/
    
    std_setup();
    
    int n_runs = 1E1;
    int n_steps = 1E3;
    double **W=0;
    
    double stop_time = 1.0;
    double dt = stop_time/n_steps;
    
    gen_BM( dt, n_steps, W, n_runs, BM_mode_2 );
    
    
    for (int run=0; run<n_runs; run++)
    {
        double ito_sum = 0.0;
        
        for (int step=0; step<n_steps-1; step++)
            ito_sum += W[run][step]*(W[run][step+1]-W[run][step]);
        
        cout << "ito: " << ito_sum << "\t" << 0.5*( W[run][n_steps-1]*W[run][n_steps-1] - stop_time ) << "\t" << fabs( ito_sum - 0.5*(W[run][n_steps-1]*W[run][n_steps-1] - stop_time) ) <<endl;
        
        
        double strat_sum = 0.0;
        
        for (int step=0; step<n_steps-1; step++)
            strat_sum += (W[run][step] + W[run][step+1] )*(W[run][step+1]-W[run][step])/2;
        
        cout << "strat: " << strat_sum << "\t" << 0.5*( W[run][n_steps-1]*W[run][n_steps-1] ) << "\t" << fabs( strat_sum - 0.5*(W[run][n_steps-1]*W[run][n_steps-1] ) ) <<endl;
        
    }
    
    ml_free( W, n_runs );
    
    std_exit();
}
 
 
 
 
