
#include <mathlib/math/std_math.h>
#include <mathlib/tools/graphing/plot2D.h>

#include <mathlib/math/random/ml_random.h>
#include <mathlib/math/grids/grid2D.h>
//#include <mathlib/math/grids/grid2D_Del2_FD.h>
#include <mathlib/math/grids/extra/plots.cpp>
#include <mathlib/math/PDE/laplacian_hdaf.h>

#include <mathlib/link.cpp>
#include <mathlib/non-lib_things.h>


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




float distance_squared( float const & X11, float const & X12, float const & X21, float const & X22, float const & Y1, float const & Y2 )
{
	float y1 = Y1 - X11;
	float y2 = Y2 - X12;
	float x1 = X21 - X11;
	float x2 = X22 - X12;
	
	float dot = x1*y1+x2*y2;
	
	if ( dot < 0.f )
	    return y1*y1+y2*y2;
	
	float lx2 = x1*x1+x2*x2;
	
	if ( dot > lx2 )
		return (Y1-X21)*(Y1-X21)+(Y2-X22)*(Y2-X22);
    
	return y1*y1+y2*y2 - dot*dot/lx2;
}







float X11;
float X12;
float X21;
float X22;





ml_color cmap( double z )
{
    float x = sqrt(z)/.5;
    return ml_color(0,x,0);
}



ml_color f(double x, double y)
{
    return cmap( distance( X11, X12, X21, X22  ,x,y  ) );
}


int main()
{
    std_setup();
    
    ml_random rng;
    
    int N = 1024;
    
    plot2D plt( -10,10,-10,10, N, N);
    plt = ml_white;
    
    
    for (int k=0; k<1000; k++)
    {
        cout << k << endl;
        
        X11 = rng.gen_double()*20-10;
        X12 = rng.gen_double()*20-10;
        X21 = rng.gen_double()*20-10;
        X22 = rng.gen_double()*20-10;
        
        
        plt = f;
        
        plt.ptLine(  X11, X12, X21, X22, ml_red );
        
        sprintf(fname, "/workspace/output/SDE/project1/test/%04d.png", k );
        plt.png(fname);
    }
    
    std_exit();
}


