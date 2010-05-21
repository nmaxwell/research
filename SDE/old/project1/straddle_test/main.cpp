
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



bool straddle_test( float const & x1, float const & x2, float const & x3, float const & x4, float const & y1, float const & y2, float const & y3, float const & y4 )
{
    if ( max(x1,x2) < min(x3,x4) )
        return false;
    if ( max(x3,x4) < min(x1,x2) )
        return false;
    if ( max(y1,y2) < min(y3,y4) )
        return false;
    if ( max(y3,y4) < min(y1,y2) )
        return false;
    
    float z1 = (x3-x1)*(y2-y1)-(y3-y1)*(x2-x1);
    float z2 = (x4-x1)*(y2-y1)-(y4-y1)*(x2-x1);
    
    //if ( z1==0.f && z2==0.f ) return true;
    
    if ( (z1<0) xor (z2<0) ) return true;
    
    return false;
}

int intersection( double & x, double & y, double X11, double X12, double X21, double X22, double Y11, double Y12, double Y21, double Y22  )
{
    double y11 = Y11-X11;
    double y12 = Y12-X12;
    double y21 = Y21-X11;
    double y22 = Y22-X12;
    
    double x1 = X21-X11;
    double x2 = X22-X12;
    
    double y1 = y21-y11;
    double y2 = y22-y12;
    
    double lx = sqrt(x1*x1+x2*x2);
    x1 /= lx;
    x2 /= lx;
    
    double ly = sqrt(y1*y1+y2*y2);
    y1 /= ly;
    y2 /= ly;
    
    double det = -x1*y2 + y1*x2;
    
    if ( det == 0.f )
        return -1;
    
    double sx = (y1*y12-y2*y11)/det;
    double sy = (x1*y12-x2*y11)/det;
    
    x = sx*x1+X11;
    y = sx*x2+X12;
    
    if ( min(X11,X21) <= x and x <= max(X11,X21) and min(X12,X22) <= y and y <= max(X12,X22) and
         min(Y11,Y21) <= x and x <= max(Y11,Y21) and min(Y12,Y22) <= y and y <= max(Y12,Y22) )
        return 0;
    else
        return 1;
}





int main()
{
    std_setup();
    
    
    int N = 300;
    
    grid2D<> u( N, 0, 1, N, 0, 1 );
    
    
    
    
    plot2D plt( 0,1,0,1, N, N);
    plt = ml_white;
    ml_random rng;
    
    int count = 0;
    
    
    for (int k=0; k<100000; k++)
    {
        cout << k << endl;
        
        plt = ml_white;
        
        double x1 = rng.gen_double();
        double x2 = rng.gen_double();
        double x3 = rng.gen_double();
        double x4 = rng.gen_double();
        double y1 = rng.gen_double();
        double y2 = rng.gen_double();
        double y3 = rng.gen_double();
        double y4 = rng.gen_double();
        
        double x,y;
        
        int res = intersection( x,y, x1, y1, x2, y2, x3, y3, x4, y4 );
        bool straddle = straddle_test( x1,x2,x3,x4,y1,y2,y3,y4 );
        
        plt.ptDot( x,y, 5, ml_green );
        
        if ( ( straddle and res == 1 ) or ( !straddle and res == 0 ) )
        {
            plt.ptLine( x1,y1,x2,y2, ml_black );
            plt.ptLine( x3,y3,x4,y4, ml_black );
            
          // sprintf(fname, "/workspace/output/SDE/project1/test/%04d.png", k );
          // plt.png(fname);
            
            count ++;
        }
        
        
    }
    
    cout << count << endl;
    
    
    std_exit();
}


