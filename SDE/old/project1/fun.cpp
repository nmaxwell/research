
#include <mathlib/math/std_math.h>
#include <mathlib/tools/graphing/plot2D.h>

#include <mathlib/math/random/ml_random.h>

#include <mathlib/link.cpp>
#include <mathlib/non-lib_things.h>



#include "Polygon2D.cpp"



float boundary_def(float & x, float & y, double t)
{
	x = cos(t*_2pi)*cos(t*_2pi*1.5)*7;
	y = sin(t*_2pi)*cos(t*_2pi*1.5)*7;
}

	




int main()
{
    std_setup();
    
    ml_random rng;
    
    
    
    Polygon2D boundary;
    
    //P.set_inscribed_circle(3, 7);
    
    boundary.init(100);
    
    {
		int n=boundary.n_points-1;
		
		for (int k=0; k<n; k++)
		{
			double t1 = double(k)/n;
			double t2 = double(k+1)/n;
			
			boundary_def( boundary.x_point(k),boundary.y_point(k), t1 );
			boundary_def( boundary.x_point(k+1),boundary.y_point(k+1), t2 );
		}
	}
    
    
    
    
    
    plot2D plt( -10.0, 10.0, -10.0, 10.0, 1000, 1000 );
    plt = ml_white;
    
    plot( boundary, plt, ml_black );
    
    for (int k=0; k<100000; k++)
    {
        float x = rng.gen_double()*20-10.0;
        float y = rng.gen_double()*20-10.0;
        
        if (boundary.interior_test(x,y))
            plt.ptDot(x,y,1, ml_blue);
        else
            plt.ptDot(x,y,1, ml_red);
    }
    
    
    plt.png("/workspace/output/scratch/out.png");
    
    std_exit();
}





