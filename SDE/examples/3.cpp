
#define ML_NOT_USE_FFTW_MALLOC

#include <mathlib/math/random/ml_random.h>
#include <mathlib/math/SDE/BM.h>
#include <mathlib/math/grids/grid1D.h>
#include <mathlib/math/grids/grid2D.h>
//#include <mathlib/math/grids/extra/plots.cpp>

#include <mathlib/link.cpp>
#include <mathlib/non-lib_things.h>

//#include <boost/math/special_functions/erf.hpp>




template< class T >
void output( T * data, int n, const char * fname, const char * delim1="\t", const char * delim2="\n" )
{
    ofstream out;
    out.open(fname, fstream::out);
    assert(out.good() && out.is_open());
    
    for ( int k=0; k<n; k++ )
        out << data[k] << delim2;
    
    out.close();
}

template< class T1, class T2 >
void output( T1 * data1, T2 * data2, int n, const char * fname, const char * delim1="\t", const char * delim2="\n" )
{
    ofstream out;
    out.open(fname, fstream::out);
    assert(out.good() && out.is_open());
    
    for ( int k=0; k<n; k++ )
        out << data1[k] << delim1 << data2[k] << delim2;
    
    out.close();
}



/*



import math
import pylab as p
import sys
sys.path.append("/workspace/mathlib/tools/python/")
from read_file import *

for n in range(1,10):
    p.plot( read_file( "out_" + str(n) ) )

p.plot ( read_file( "U_mean" )  )
p.show()

Y = read_file( "U_mean" )
X = [ float(k) / len(Y) for k in range(len(Y))]
Y = [ math.exp(x) for x in X ]
p.plot(X,Y,X, read_file( "U_mean" )  )
p.show()





*/



int main()
{
    /*
     * functions along a brownian path
     *  exp( t + 0.5 Wt )
     * 
     */
    
    std_setup();
    
    int n_runs = 1E4;
    int n_steps = 1E3;
    double **W=0;
    
    double stop_time = 1.0;
    double dt = stop_time/n_steps;
    
    gen_BM( dt, n_steps, W, n_runs, BM_mode_2 );
    
    
    double * U_mean = ml_alloc<double > (n_steps);
    
    for ( int step=0; step<n_steps; step++ )
    {
        U_mean[step] = 0;
        
        for (int run=0; run<n_runs; run++)
            U_mean[step] += exp( step*dt + 0.5*W[run][step] )/n_runs;
    }
    
    output( U_mean, n_steps, "/workspace/output/temp/U_mean" );    
    
    double * M = ml_alloc<double > (n_steps );
    
    for (int run=0; run<min(n_runs,50); run++)
    {
        for ( int step=0; step<n_steps; step++ )
            M[step] = exp( step*dt + 0.5*W[run][step] );
        
        sprintf(fname, "/workspace/output/temp/out_%d", run);
        output( M, n_steps, fname );
    }
    
    ml_free( U_mean );
    ml_free( W, n_runs );
    
    std_exit();
}



