
#include <mathlib/math/std_math.h>


#include <mathlib/math/SDE/BM.h>





#include <mathlib/link.cpp>
#include <mathlib/non-lib_things.h>




double drift( double X, double t )
{
	return 0.0;
	
	
}

double volatility( double X, double t )
{
	return 1.0;
	
}



int main()
{
	std_setup();
    
    
    int R = 100;
    int n_steps = 500;
    int n_runs = 100;
    double **W=0, *time=0, *X_true=0;
    
    double stop_time = 2.0;
    double dt = stop_time/n_steps;
    
    gen_BM( dt/R, n_steps*R, W, n_runs, BM_mode_1 );
    
    time = ml_alloc<double > ( n_steps );
    for (int k=0; k<n_steps; k++ )
        time[k] = dt*k;
    
    double **X = ml_alloc<double> ( n_runs, n_steps );
    
    
    double X0 = 0.0;
    
    for ( int run=0; run<n_runs; run++ )
    {
        X[run][0] = X0;
        
        for ( int j=1; j<n_steps; j++ )
        {
            X[run][j] = X[run][j-1] + drift(X[run][j-1], dt*j)*dt + volatility(X[run][j-1], dt*j)*(W[run][j*R]-W[run][(j-1)*R]);
        }
    }
    
    
    
    
    double * X_mean = ml_alloc<double> (n_steps);
    
    for ( int step=0; step<n_steps; step++ )
        X_mean[step] = 0;
    
    for ( int run=0; run<n_runs; run++ )
    for ( int step=0; step<n_steps; step++ )
        X_mean[step] += X[run][step]/n_runs;
    
    output( X_mean, n_steps, "/workspace/output/SDE/test/X_mean" );
    
    output( X, n_runs, n_steps, "/workspace/output/SDE/test/X" );
    
    output( time, n_steps, "/workspace/output/SDE/test/time" );
    
    ml_free( X, n_runs );
    ml_free( X_mean );
    ml_free( W, n_runs );
    ml_free( time );
    
    std_exit();
}




