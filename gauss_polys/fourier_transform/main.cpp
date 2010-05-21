
#include <mathlib/math/std_math.h>
#include <mathlib/math/polynomials/ml_poly.h>


#include <mathlib/link.cpp>

#define include_integration

#include <mathlib/non-lib_things.h>


struct transform_args
{
    transform_args( double a, int n, double xi ): a(a), n(n), xi(xi) {}
    
    double a;
    int n;
    double xi;
};



double real_part(double x, void * arg_ptr )
{
    double a = ((transform_args *) arg_ptr)->a;
    double xi = ((transform_args *) arg_ptr)->xi;
    int n = ((transform_args *) arg_ptr)->n;
    
    if ( fabs(x) < 1.0 )
        return cos(x*xi)*exp(-a*x*x)*pow(x,n);
    
    if (n%2)
        return cos(x*xi)*exp( log(fabs(x))*n -a*x*x )*sign(x);
    else
        return cos(x*xi)*exp( log(fabs(x))*n -a*x*x );
}

double imag_part(double x, void * arg_ptr )
{
    double a = ((transform_args *) arg_ptr)->a;
    double xi = ((transform_args *) arg_ptr)->xi;
    int n = ((transform_args *) arg_ptr)->n;
    
    if ( fabs(x) < 1.0 )
        return sin(x*xi)*exp(-a*x*x)*pow(x,n);
    
    if (n%2)
        return sin(x*xi)*exp( log(fabs(x))*n -a*x*x )*sign(x);
    else
        return sin(x*xi)*exp( log(fabs(x))*n -a*x*x );
}



/*




from math import *
import pylab as p
import sys
sys.path.append("/workspace/mathlib/tools/python/")
from read_file import *


F = read_file( "F" )
FTF = read_file( "FTF" )

p.plot( [ line[0] for line in F ], [ line[1] for line in F ] )
p.plot( [ line[0] for line in FTF ], [ line[1] for line in FTF ] )
p.plot( [ line[0] for line in FTF ], [ line[2] for line in FTF ] )

p.show()




*/


int main()
{
    std_setup();
    
    int n=6;
    double a=1.0;
    transform_args args( a,n,0 );
    
    {
        double xi_f = 10.0;
        int N = 500;
        double h = xi_f/N;
        
        double *rp = ml_alloc<double> (2*N+1);
        double *ip = ml_alloc<double> (2*N+1);
        double *xi = ml_alloc<double> (2*N+1);
        
        for ( int j=-N; j<=N; j++ )
        {
            args.xi = h*j;
            
            xi[j+N] = args.xi;
            rp[j+N] = integrate_R ( real_part, (void *)(&args) );
            ip[j+N] = integrate_R ( imag_part, (void *)(&args) );
        }
        
        output( xi, rp, ip, 2*N+1, "/workspace/output/temp/FTF" );
    }
    
    //---------------------
    
    {
        double xf = 10.0;
        int N = 500;
        double h = xf/N;
        
        double *F = ml_alloc<double> (2*N+1);
        double *X = ml_alloc<double> (2*N+1);
        
        for ( int j=-N; j<=N; j++ )
        {
            double x = h*j;
            
            X[j+N] = x;
            F[j+N] = exp(-a*x*x)*pow(x,n);
        }
        
        output( X, F, 2*N+1, "/workspace/output/temp/F" );
    }
    
    //-----------------------
    
    {
        double xi_f = 10.0;
        int N = 500;
        double h = xi_f/N;
        
        double *FTR = ml_alloc<double> (2*N+1);
        double *FTI = ml_alloc<double> (2*N+1);
        double *Xi = ml_alloc<double> (2*N+1);
        
        ml_poly<double > P0(n); P0 = 0.0;
        ml_poly<double > P1(n); P1 = 0.0;
        ml_poly<double > P2(n); P2 = 0.0;
        ml_poly<double > P3(n); P3 = 0.0;
        
        for ( int k=0; k<=n; k++ )
        {
            if (k%4 == 0) P0[k] = binom(n,k)*dfact2n_n[n-k]*sqrtpi*pow(a,-n)*pow(2.0,k-2*n);
            if (k%4 == 1) P1[k] = binom(n,k)*dfact2n_n[n-k]*sqrtpi*pow(a,-n)*pow(2.0,k-2*n);
            if (k%4 == 2) P2[k] = binom(n,k)*dfact2n_n[n-k]*sqrtpi*pow(a,-n)*pow(2.0,k-2*n);
            if (k%4 == 3) P3[k] = binom(n,k)*dfact2n_n[n-k]*sqrtpi*pow(a,-n)*pow(2.0,k-2*n);
        }
        
        ml_poly<double > PR(n); P0 = 0.0;
        ml_poly<double > PI(n); P1 = 0.0;
        
        PR = P0; PR -= P2;
        PI = P1; PI-= P3;
        
        
        for ( int j=-N; j<=N; j++ )
        {
            double xi = h*j;
            
            Xi[j+N] = xi;
            
            FTR[j+N] = exp(-xi*xi/(a*4))*PR(xi);
            FTI[j+N] = exp(-xi*xi/(a*4))*PI(xi);
        }
        
        output( Xi, FTR, FTI, 2*N+1, "/workspace/output/temp/test" );
    }
    
    
    
    
    std_exit();
}




