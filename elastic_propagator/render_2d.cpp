		
#include <mathlib/math/grids/grid2D.h>
#include <mathlib/math/grids/extra/plots.cpp>
#include <mathlib/math/grids/grid_file.cpp>

#include <mathlib/link.cpp>

#include <mathlib/non-lib_things.h>



void plot(const char * fname, grid2D< > & G1, grid2D< > & G2, ml_color (*cmap)(double, double ) )
{	
	plot2D plot(G1.a1,G1.b1,G1.a2,G1.b2, G1.n1,G1.n2);
	
	for (int i = 0; i<plot.NX; i++)
	for (int j = 0; j<plot.NY; j++)
	{
		plot.set_px(i,j, cmap(G1(i,j),G2(i,j)) );
	}
	
	plot.png(fname);
}

ml_color cmap(double x1, double x2)
{
	x1 *= 10;
    x2 *= 10;
    
	float s1 = atan(x1)/pi+0.5;
    float s2 = atan(x2)/pi+0.5;
    
	return ml_color(s1,0,s2);
	
}

int main()
{
	std_setup();
	
	grid2D<> G1;
    grid2D<> G2;
	
	/*sprintf(fname,"/workspace/output/anisotropic_propagate_2d/out_dat/C2.dat" );
	readFile(G,fname);
	sprintf(fname,"/workspace/output/anisotropic_propagate_2d/out_png/C2.png" );
	plotGrid2D_1(G,fname,cmap);*/
	
    
	for( int n = 0; n<9000; n+=1)
	{
		sprintf(fname,"/workspace/output/elastic_propagate_2d/out_dat/u1_%05d.dat" ,n);
        cout << fname << endl;
		if ( readFile(G1,fname) ) break;
		
        sprintf(fname,"/workspace/output/elastic_propagate_2d/out_dat/u2_%05d.dat" ,n);
        cout << fname << endl;
		if ( readFile(G2,fname) ) break;
        
		sprintf(fname,"/workspace/output/elastic_propagate_2d/out_png/%05d.png" ,n);	    
		plot( fname, G1, G2, cmap );
	}	
}


