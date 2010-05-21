
#include <mathlib/math/std_math.h>
#include <mathlib/tools/graphing/plot2D.h>

#include <mathlib/math/random/ml_random.h>
#include <mathlib/math/grids/grid2D.h>
#include <mathlib/math/grids/extra/plots.cpp>
#include <mathlib/math/PDE/laplacian_hdaf.h>

#include <mathlib/link.cpp>
#include <mathlib/non-lib_things.h>

#include <gsl/gsl_integration.h>

#include "Polygon2D.cpp"




double integrand(double t, void * args)
{
    double * dargs = (double *)args;
  /*  double x = dargs[0]-cos(_2pi*t)*dargs[2];
    double y = dargs[1]-sin(_2pi*t)*dargs[2];
    if ( fabs(x) >= 1E-12 and fabs(y) >= 1E-12 )
        return 20.0*log(x*x+y*y)/(-_2pi*2);
    else return 0.0;*/
    
  /*  double x_eval = dargs[0];
    double y_eval = dargs[1];
    
    double x1 = -2.0;
    double y1 = 2.0;
    double x2 = 7.0;
    double y2 = 7.0;
    
    x_eval -= x1;
    y_eval -= y1;
    x2 -= x1;
    y2 -= y1;
    
    double L = sqrt(x2*x2+y2*y2);
    if (fabs(L) >= 1E-12)
    {
        x2 /= L;
        y2 /= L;
    }
    else return 0;
    
    double x = (x_eval*y2-y_eval*x2);
    double y = (x_eval*x2+y_eval*y2);
    
    y -= L*t;
    
    if ( fabs(x) >= 1E-12 and fabs(y) >= 1E-12 )
        return 1.0*log(x*x+y*y)/(-_2pi*2);
    else return 0.0;*/
    
    
    
    double x = dargs[0];
    double y = dargs[1];
    
    y -= 3.0*t;
    
    if ( fabs(x) >= 1E-12 and fabs(y) >= 1E-12 )
        return 1.0*log(x*x+y*y)/(-_2pi*2);
    else return 0.0;
    
    
}




double V_1(double x, double y)
{
    double * dargs = ml_alloc<double> (10);
    void * arg_ptr = (void *)dargs;
    
    dargs[0] = x;
    dargs[1] = y;
    
    double epsabs=1E-6;
    double epsrel=1E-10;
    
    double result,real_abserr,imag_result,imag_abserr;
    static gsl_integration_workspace * workspace = gsl_integration_workspace_alloc (1000000);
    
    gsl_function F;
    F.function = integrand;
    F.params = arg_ptr;
    
    gsl_integration_qag (
        &F,
        0,
        1,
        epsabs,
        epsrel,
        200000,
        6,
        workspace,
        &result,
        &real_abserr );
    
    return result;
}


double antideriv(double x, double z)
{
    return z*log(z*z+x*x)/2 + atan(z/fabs(x))-z;
}

double V_2(double x, double y)
{
/*    double x_eval = x;
    double y_eval = y;
    
    double x1 = -2.0;
    double y1 = 2.0;
    double x2 = 7.0;
    double y2 = 7.0;
    
    x_eval -= x1;
    y_eval -= y1;
    x2 -= x1;
    y2 -= y1;
    
    double L = sqrt(x2*x2+y2*y2);
    if (fabs(L) >= 1E-12)
    {
        x2 /= L;
        y2 /= L;
    }
    else return 0;
    
    x = (x_eval*y2-y_eval*x2);
    y = (x_eval*x2+y_eval*y2);
    
    return (antideriv(x,y)-antideriv(x,y-L))/(-_2pi);*/
    
   /* double s = sqrt(x*x+y*y);
    if ( fabs(s) >= 1E-14 )
        return log(s)/(-_2pi);
    else return 0.0;*/
    
    return (antideriv(x,y)-antideriv(x,y-3.0))/(-_2pi);
}







ml_color color_map_green( double z )
{
	float x = atan((double)z)/pi+0.5;
	return ml_color(0,x,0);
}

ml_color color_map_blue( double z )
{
	float x = atan((double)z)/pi+0.5;
	return ml_color(0,0,x);
}

ml_color color_map_red( double z )
{
	float x = atan((double)z)/pi+0.5;
	return ml_color(x,0,0);
}

ml_color color_map_error( double z )
{
    //cout << z << endl;
	float x = float(  (log10(fabs(z))+12)/12  );
	return ml_color(x,x,x);
}


int main()
{
    std_setup();
    
    ml_random rng;
    
    
    
    {
		// output analytic solution
		
        grid2D<double, double > F1( 500, -10.0, 10.0, 500, -10.0, 10.0 );
        grid2D<double, double > F2( F1 );
        grid2D<double, double > G1( F1 );
        grid2D<double, double > G2( F1 );
        
        F1 = V_1;
        F2 = V_2;
        
        laplacian_2d_hdaf Del2;
        Del2.init( F1.n1, F1.n2, F1.b1-F1.a1, F1.b2-F1.a2, 24, 24, .5, .5 );
        
        Del2.execute( F1.array, G1.array );
        Del2.execute( F2.array, G2.array );
        
        
        plotGrid2D_1( F1,  "/workspace/output/scratch/F1.png", color_map_green );
        plotGrid2D_1( F2,  "/workspace/output/scratch/F2.png", color_map_green );
        plotGrid2D_1( G1,  "/workspace/output/scratch/G1.png", color_map_error );
        plotGrid2D_1( G2,  "/workspace/output/scratch/G2.png", color_map_error );
	}
    
    
    
    
    
    
    std_exit();
}











