#include <mathlib/math/std_math.h>

#include <mathlib/math/SDE/BM.h>

#include <mathlib/link.cpp>
#include <mathlib/non-lib_things.h>




double drift( double t )
{
	return (0.7*t+cos(t*_2pi))*1.0;
}

double volatility( double X, double t )
{
	return 1.0;
	
}



int main()
{
	std_setup();
    
    
    // setup:
    
    int n_steps = 500;
    int n_runs = 10000;
    
    double X0 = 0.0;
    double stop_time = 2.0;
    
    double dt = stop_time/n_steps;
    
    double * time = ml_alloc<double > ( n_steps );
    for (int k=0; k<n_steps; k++ )
        time[k] = dt*k;
    
    
    
    double **W=0, *X_true=0;
    double **X = ml_alloc<double> ( n_runs, n_steps );
    double * X_mean = ml_alloc<double> (n_steps);
    double **M = ml_alloc<double> ( n_runs, n_steps );
    
    
    
    // generate the BMs:
    
    gen_BM( dt, n_steps, W, n_runs, BM_mode_2 );
    
    
    // integrate the SDE
    
    for ( int run=0; run<n_runs; run++ )
    {
        X[run][0] = X0;
			
        for ( int j=1; j<n_steps; j++ )
        {
			X[run][j] = X[run][j-1] + ( drift(dt*j) + drift(dt*(j-1)) )*dt/2 + volatility(X[run][j-1], dt*(j-1))*(W[run][j]-W[run][j-1]);
        }
    }
    
    // compute the Girsanov change of measure at every t
    
    for ( int run=0; run<n_runs; run++ )
    for ( int k=0; k<n_steps-1; k++ )
    {
        double sum1 = 0;
        for (int j=0; j<=k; j++)
            sum1 += drift(dt*j)*drift(dt*j)*dt;
        
        double sum2 = 0;
        for (int j=0; j<=k; j++)
            sum2 += drift(dt*j)*(W[run][j+1]-W[run][j]);
        
        M[run][k] = exp(-sum2 -0.5*sum1);
    }
    
    // compute expectation of X, output results so far.
    
    for ( int step=0; step<n_steps; step++ )
        X_mean[step] = 0;
    
    for ( int run=0; run<n_runs; run++ )
    for ( int step=0; step<n_steps; step++ )
        X_mean[step] += X[run][step]/n_runs;
    
    output( X, max(n_runs,min(100,n_runs)), n_steps, "/workspace/output/SDE/test/X" );
    output( X_mean, n_steps, "/workspace/output/SDE/test/X_mean" );
    output( time, n_steps, "/workspace/output/SDE/test/time" );
    output( M, max(n_runs,min(100,n_runs)), n_steps, "/workspace/output/SDE/test/M" );
    
    
    // compute expectation of X, with respect to new measure, output results.
    
    for ( int run=0; run<n_runs; run++ )
    for ( int step=0; step<n_steps; step++ )
        X[run][step] *= M[run][step];
    
    for ( int step=0; step<n_steps; step++ )
        X_mean[step] = 0;
    
    for ( int run=0; run<n_runs; run++ )
    for ( int step=0; step<n_steps; step++ )
        X_mean[step] += X[run][step]/n_runs;
    
    
    output( X, max(n_runs,min(100,n_runs)), n_steps, "/workspace/output/SDE/test/XM" );
    output( X_mean, n_steps, "/workspace/output/SDE/test/QX" );
    
    
    
    
    
    
    ml_free( X, n_runs );
    ml_free( X_mean );
    ml_free( W, n_runs );
    ml_free( time );
    ml_free( X_mean );
    
    std_exit();
}




