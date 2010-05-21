
#include <mathlib/math/grids/grid2D.h>		
#include <mathlib/math/grids/grid3D.h>
#include <mathlib/math/grids/extra/plots.cpp>
#include <mathlib/math/grids/grid_file.cpp>

#include <mathlib/link.cpp>

#include <mathlib/non-lib_things.h>


/*
void plot(const char * fname, grid2D< > & G1, grid2D< > & G2, ml_color (*cmap)(double, double ) )
{	
	plot2D plot(G1.a1,G1.b1,G1.a2,G1.b2, G1.n1,G1.n2);
	
	for (int i = 0; i<plot.NX; i++)
	for (int j = 0; j<plot.NY; j++)
	{
		plot.set_px(i,j, cmap(G1(i,j),G2(i,j)) );
	}
	
	plot.png(fname);
}*/

ml_color cmap(double x)
{
	x *= 10;
    
	float s = atan(x)/pi+0.5;
    
	return ml_color(s,s,s);
	
}


class euVec
{
public:
	double _1;
	double _2;
	double _3;
	
	euVec( double a, double b, double c ):_1(a),_2(b),_3(c) {}
	
	euVec operator* (double s)
	{
		return euVec( _1*s, _2*s, _3*s );
	}
	
	euVec operator+ (euVec v)
	{
		return euVec( _1+ v._1, _2+ v._2, _3+ v._3 );
	}
	
};

void slice( grid3D<> &F, grid2D<> &s, euVec offset, euVec y_hat, euVec x_hat )
{
	for (int i=0; i<s.n1; i++)
	for (int j=0; j<s.n2; j++)
	{
		euVec x = x_hat*s.x1(i)+y_hat*s.x2(j)+offset;
		
		s(i,j) = F( x._1, x._2, x._3 );
	}
}

double f(double x, double y, double z)
{
	double r2 = x*x+y*y+z*z;
	
	return exp(-r2);
}



int main()
{
	std_setup();
	
	
	/*
	grid3D<> u( 512, -10,10, 512, -10,10, 512, -10,10 );
	
	grid2D<> s( 512, -10,10, 512, -10,10  );
	
	u = f;
	
	euVec offset(0,0,0);
	euVec x_hat(1,0,0);
	euVec y_hat(0,1,0);
	
	slice( u, s, offset, y_hat, x_hat );
	
	plotGrid2D_1(s,"/workspace/output/elastic_propagate_3d/out.png",cmap);
	*/
	
	
	grid3D<> u1;
    grid3D<> u2;
    grid3D<> u3;
    
    grid2D<> s1( 500, -10,10, 500, -10,10  );
    grid2D<> s2(s1);
    grid2D<> s3(s1);
	
	euVec offset(0,0,0);
	euVec x_hat(1,0,0);
	euVec y_hat(0,1,0);
	
    
	for( int n = 0; n<9000; n+=1)
	{
		sprintf(fname,"/workspace/output/elastic_propagate_3d/out_dat/u1_%05d.dat" ,n);
        cout << fname << endl;
		if ( readFile(u1,fname) ) break;
		
        sprintf(fname,"/workspace/output/elastic_propagate_3d/out_dat/u2_%05d.dat" ,n);
        cout << fname << endl;
		if ( readFile(u2,fname) ) break;
		
		sprintf(fname,"/workspace/output/elastic_propagate_3d/out_dat/u3_%05d.dat" ,n);
        cout << fname << endl;
		if ( readFile(u3,fname) ) break;
         
        slice( u1, s1, offset, y_hat, x_hat );
        slice( u2, s2, offset, y_hat, x_hat );
        slice( u3, s3, offset, y_hat, x_hat );        
        
		sprintf(fname,"/workspace/output/elastic_propagate_3d/out_png/%05d_1.png" ,n);	    
		plotGrid2D_1(s1,fname,cmap);
		sprintf(fname,"/workspace/output/elastic_propagate_3d/out_png/%05d_2.png" ,n);	    
		plotGrid2D_1(s2,fname,cmap);
		sprintf(fname,"/workspace/output/elastic_propagate_3d/out_png/%05d_3.png" ,n);	    
		plotGrid2D_1(s3,fname,cmap);
		
	}
	
}


