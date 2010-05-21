
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

for n in range(1,10):
    p.plot( time, read_file( "X_" + str(n) ) )

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
     * Euler-Maruyama Method, convergence test
     * 
     * 
     */
    
    std_setup();
    
    mu = 1.0;
    lambda = 2.0;
    double X0 = 1.0;   
    double stop_time = 1.0;
    int n_steps = 512;
    int n_runs = 1000;
    int n_R = 7;
    
    double Xerr[n_R][n_runs];
    double X_true[n_R][n_runs];
    
    for (int p=0; p<n_R; p++)
    {
        cout << p << "\t" << n_steps*pow(2,p)*n_runs << endl;
        int R = pow(2,p);
        double **W=0, Xt, dt = stop_time/n_steps;
        
        gen_BM( dt/R, n_steps*R, W, n_runs, BM_mode_3 );
        
        Xt = 0;
        
        for ( int run=0; run<n_runs; run++ )
        {
            Xt = X0;
            
            for ( int j=1; j<n_steps; j++ )
                Xt += f(Xt, dt*j)*dt + g(Xt, dt*j)*(W[run][j*R]-W[run][(j-1)*R]);
            
            X_true[p][run] = X0*exp( (lambda-0.5*mu*mu)*stop_time + mu*W[run][R*n_steps-1] );
            
            //Xerr[p][run] = 100.0*fabs((( Xt - X_true[p][run] )/X_true[p][run]));
            Xerr[p][run] = fabs( Xt - X_true[p][run] );
        }
        
        ml_free( W, n_runs );
    }
    
    cout << endl;
    
    
    for (int p=0; p<n_R; p++)
    {
        double mean_err = 0;
        for (int k=0; k<n_runs; k++)
            mean_err += Xerr[p][k]/n_runs;
        
        cout << log10((stop_time/n_steps)/(pow(2,p))) << "\t" << log10(mean_err) << endl;
    }
    
    
    
    
    std_exit();
}









