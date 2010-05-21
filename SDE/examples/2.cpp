

//#define ML_NOT_USE_FFTW_MALLOC

#include <mathlib/math/SDE/BM.h>

#include <mathlib/link.cpp>
#include <mathlib/non-lib_things.h>




/*



import math
import pylab as p
import sys
sys.path.append("/workspace/mathlib/tools/python/")
from read_file import *

n = 10

for k in range(1,n):
    print p.std( read_file( "out_" + str(k) ) ), math.sqrt( float(k)/n ), p.mean( read_file( "out_" + str(k) ) )



*/



int main()
{
    std_setup();
    
    int n_runs = 1E4;
    int n_steps = 1E3;
    double **W=0;
    
    double stop_time = 1.0;
    double dt = stop_time/n_steps;
    
    double t1,t2;
    
    t1 = get_real_time();
    
    gen_BM( dt, n_steps, W, n_runs, BM_mode_3 );
    
    t2 = get_real_time();
    
    cout << "time: " << t2-t1 << endl;
        
    if (0)
    {
        double * M = ml_alloc<double > (n_runs );
        
        for ( int step=0; step<n_steps; step++ )
        {
            for (int run=0; run<n_runs; run++)
                M[run] = W[run][step];
            
            sprintf(fname, "/workspace/output/temp/out_%d", step);
            output( M, n_runs, fname );
        }
    }
    
    ml_free( W, n_runs );
    
    std_exit();
}


