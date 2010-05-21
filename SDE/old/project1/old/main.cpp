
#include <mathlib/math/std_math.h>
#include <mathlib/tools/graphing/plot2D.h>

#include <mathlib/math/random/ml_random.h>
#include <mathlib/math/grids/grid2D.h>
#include <mathlib/math/grids/extra/plots.cpp>
#include <mathlib/math/PDE/laplacian_hdaf.h>

#include <mathlib/link.cpp>
#include <mathlib/non-lib_things.h>



#include "Polygon2D.cpp"


double antideriv(double x, double z)
{
    return z*log(z*z+x*x)/2 + atan(z/x)-z;
}




double V_linecharge(double x_eval, double y_eval, double x1, double y1, double x2, double y2, double rho )
{
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
    
    //return (antideriv(x,y)-antideriv(x,y-L))/(-_2pi);
    
    //double L = 1.7;
    //return (antideriv(x,y)-antideriv(x,y-L))/(-_4pi*L);
    //if (!is_number( log((y-L)*(y-L)+x*x)) or !is_number(log(y*y+x*x)) )
    //    return 0.0;
    
    return (y*log(y*y+x*x) + atan(y/x)*x*2.0-2.0*y-(y*log((y-L)*(y-L)+x*x) + atan((y-L)/x)*x*2.0-2.0*(y-L)))/(-_4pi*L);
}


float boundary_def(float & x, float & y, float & f, double t)
{
	// t \in [0,1]
	
	x = cos(t*_2pi)*(cos(t*_2pi*5)*2.5+6);
	y = sin(t*_2pi)*(cos(t*_2pi*5)*2.5+6);
	
	f = cos(_2pi*t*15)*2.0 + sin(_2pi*t*13)*3.0+sin(_2pi*t*60)*7.0 ;
    f *= 2.0;
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
    


int main()
{
    std_setup();
    
    ml_random rng;
    
    
    int n_bd=500;
    Polygon2D boundary;
    float * f = 0;
    
    
    if (1)
    {
		// setup boundary
		
		//boundary.set_inscribed_circle(3, 7);
		
		f = ml_alloc<float> (n_bd);
		boundary.init(n_bd+1);
		
		for (int k=0; k<n_bd; k++)
		{
			double t1 = double(k)/n_bd;
			double t2 = double(k+1)/n_bd;
			
			boundary_def( boundary.x_point(k),boundary.y_point(k), f[k], t1 );
			boundary_def( boundary.x_point(k+1),boundary.y_point(k+1), f[k+1], t2 );
		}
	}
    else
    {
        n_bd = 4;
        
        f = ml_alloc<float> (n_bd);
        boundary.init(n_bd+1);
        
        f[0] = 0.0;
        f[1] = 7.0;
        f[2] = 0.0;
        f[3] = 0.0;
        
        boundary.x_point(0) = -8.0;
        boundary.y_point(0) = +8.0;
        
        boundary.x_point(1) = +8.0;
        boundary.y_point(1) = +8.0;

        boundary.x_point(2) = +8.0;
        boundary.y_point(2) = -8.0;
        
        boundary.x_point(3) = -8.0;
        boundary.y_point(3) = -8.0;
        
        boundary.x_point(4) = -8.0;
        boundary.y_point(4) = +8.0;
    }
    
    
    
    
    {
		// output analytic solution
		
		plot2D plt( -10.0, 10.0, -10.0, 10.0, 512, 512 );
        grid2D<double, double > F1( plt.NX, -10.0, 10.0, plt.NY, -10.0, 10.0 );
        grid2D<double, double > F2( F1 );
		plt = ml_white;
		F1 = 0.0;
        
		for (int k=0; k<n_bd; k++)
			plt.ptLineThick( boundary.x_point(k),boundary.y_point(k), boundary.x_point(k+1), boundary.y_point(k+1), 5, color_map_red(f[k]) );
            //plt.ptLine( boundary.x_point(k),boundary.y_point(k), boundary.x_point(k+1), boundary.y_point(k+1), color );
        
		for (int i=0; i<plt.NX; i++)
        {
            cout << i << endl;
            for (int j=0; j<plt.NY; j++)
            {
                float x = plt.to_x(i);
                float y = plt.to_y(j);
                
                if (boundary.interior_test(x,y))
                {
                    double V = 0.0;
                    for (int k=0; k<n_bd; k++)
                        V += V_linecharge( x,y, boundary.x_point(k),boundary.y_point(k), boundary.x_point(k+1), boundary.y_point(k+1), f[k] );
                    F1(i,j) = exp(V);
                }
            }
        }
        
        double mn = min(F1);
        double mx = max(F1);
        
        for (int i=0; i<plt.NX; i++)
        for (int j=0; j<plt.NY; j++)
        {
            float x = plt.to_x(i);
            float y = plt.to_y(j);
            if (boundary.interior_test(x,y))
            {
                float s = (F1(i,j)-mn)/(mx-mn);
                plt.set_px(i,j, ml_color(s,s,s) );
            }
        }
        
		plt.png("/workspace/output/scratch/analytic.png");
        
	}
    
    
    
    
    
    
    std_exit();
}


