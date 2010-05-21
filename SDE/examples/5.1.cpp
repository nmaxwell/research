
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

time = read_file( "time" )

n= 1

for k in range(1,n):
    p.plot( time, read_file( "X_" + str(k) ) )

p.plot ( time, read_file( "X_mean" )  )

p.plot ( time, read_file( "X_true" )  )
p.show()





*/


double lambda = 2.0;
double mu = 1.0;


double f( double x, double t )
{
    return lambda*x;
}

double g( double x, double t )
{
    return mu*x;
    //return sqrt(fabs(x)+1.0)+1.0;
}




int main()
{
    /*
     * Euler-Maruyama Method, with the 'R' factor
     * 
     * 
     */
    
    std_setup();
    
    double X0 = 1.0;
    int R = 100;
    int n_steps = 500;
    int n_runs = 1;
    double **W=0, *time=0, *X_true=0;
    
    double stop_time = 2.0;
    double dt = stop_time/n_steps;
    
    gen_BM( dt/R, n_steps*R, W, n_runs, BM_mode_1 );
    time = reg_discr( 0, stop_time, n_steps );
    
    double **X = ml_alloc<double> ( n_runs, n_steps );
    
    for ( int run=0; run<n_runs; run++ )
    {
        X[run][0] = X0;
        
        for ( int j=1; j<n_steps; j++ )
        {
            X[run][j] = X[run][j-1] + f(X[run][j-1], dt*j)*dt + g(X[run][j-1], dt*j)*(W[run][j*R]-W[run][(j-1)*R]);
        }
    }
    
    double * X_mean = ml_alloc<double> (n_steps);
    
    for ( int step=0; step<n_steps; step++ )
        X_mean[step] = 0;
    
    for ( int run=0; run<n_runs; run++ )
    for ( int step=0; step<n_steps; step++ )
        X_mean[step] += X[run][step];
    
    for ( int step=0; step<n_steps; step++ )
        X_mean[step] /= n_runs;
    
    X_true = ml_alloc<double> (n_steps);
    for ( int k=0; k<n_steps; k++ )
        X_true[k] = X[0][0]*exp( (lambda-0.5*mu*mu)*time[k] + mu*W[0][k*R] );
    
    if (1)
    {
        output( X_mean, n_steps, "/workspace/output/temp/X_mean" );
        
        for (int run=0; run<min(n_runs,20); run++)
        {
            sprintf(fname, "/workspace/output/temp/X_%d", run);
            output( X[run], n_steps, fname );
        }
        
        output( time, n_steps, "/workspace/output/temp/time" );
        
        output( X_true, n_steps, "/workspace/output/temp/X_true" );
    }
    
    ml_free( X, n_runs );
    ml_free( X_mean );
    ml_free( X_true );
    ml_free( W, n_runs );
    ml_free( time );
    
    std_exit();
}









