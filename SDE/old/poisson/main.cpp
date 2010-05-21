
#include <mathlib/math/grids/grid2D.h>
#include <mathlib/math/grids/grid2D_Del2_FD.h>

#include <mathlib/math/grids/extra/plots.cpp>

#include <mathlib/link.cpp>


#define DEL2( G,i,j ) (G(i+1,j)+G(i-1,j)+G(i,j+1)+G(i,j-1)-4.0*G(i,j))/(G.dx1()*G.dx1())
#define DEl2_center_coefficient (-2.0)

color3 cmap(double x)
{
	float s = atan(x)/pi+0.5;
	return color3(s,s,s);
	
}


template<class T >
void relaxation_method_factors( grid2D<T,T > & factors, int order, int BD_order  )
{
    // factors(i,j) = h_1^-2 delta(i,j)_1,0 + h_2^-2 delta(i,j)_2,0
    // 
    
    static FD_Del2_LR LR;
    int M = order;
    
    factors = T(0.0);
    
    for (int i=0; i<factors.n1; i++)
    for (int j=0; j<factors.n2; j++)
    {
        {
            int L = i;
        	int R = factors.n1-i-1;
        	
        	if (L>=M && R>=M) { L = M; R = M; }
        	if (L>BD_order) { L = BD_order; }
        	if (R>BD_order) { R = BD_order; }
        	
            factors(i,j) = LR(L,R,0)/(factors.dx1()*factors.dx1());
        }
        
        {
            int L = j;
        	int R = factors.n2-j-1;
        	
        	if (L>=M && R>=M) { L = M; R = M; }
        	if (L>BD_order) { L = BD_order; }
        	if (R>BD_order) { R = BD_order; }
        	
            factors(i,j) += LR(L,R,0)/(factors.dx2()*factors.dx2());
        }
    }
}




template<class T >
void solve_poisson( grid2D<T,T > & phi, grid2D<T,T > & f, grid2D<T,T > & chi_D, int order, T eps_max= 1E-6 )
{
    // solves Posson's equation in 2D. Phi is the 'initial condution',
    // and chi_D is the characeristic funciton of the Dirichlet boundary.
    // Where chi_D != 0.0, the boundary condition is taken as phi at those points.\
    //
    // R = f - Del2 phi 
    
    int BD_order = 8;
    
    grid2D<T,T > factors,R;
    factors.copy(phi);
    R.copy(phi);
    
    relaxation_method_factors(factors, order, BD_order );
    
    Del2_FD(
	    phi,
	    R,
	    order,
	    GRID_NONPER_BD_1,
	    BD_order );
        	
	LC(R, R,-1.0, f,+1.0 ); 
	
    int count = 1;
    
    while ( maxNorm(R) >= eps_max && count < 5000 )
    {
        cout << count << "\t" << maxNorm(R) <<  endl;
    	count++;
    	    	
    	for (int i=0; i<phi.n1; i++)
    	for (int j=0; j<phi.n2; j++)
    	{
    	    if (chi_D(i,j) == 0.0)
    	        phi(i,j) += R(i,j)/factors(i,j);
    	}
    	
        Del2_FD(
    	    phi,
    	    R,
    	    order,
    	    GRID_NONPER_BD_1,
    	    BD_order );
            	
    	LC(R, R,-1.0, f,+1.0 );
    }
}


template<class T >
void solve_poisson2( grid2D<T,T > & phi, grid2D<T,T > & f, grid2D<T,T > & chi_D, int order, T eps_max= 1E-6 )
{    
    grid2D<T,T > R;
    R.copy(phi);
        
    for (int i=0; i<phi.n1; i++)
    for (int j=0; j<phi.n2; j++)
        R(i,j) = f(i,j)-DEL2(phi,i,j); 
    
    int count = 1;
    
    while ( L2norm(R) >= eps_max && count < 5000 )
    {
        cout << count << "\t" << L2norm(R) <<  endl;
    	count++;
    	    	
    	for (int i=0; i<phi.n1; i++)
    	for (int j=0; j<phi.n2; j++)
    	{
    	    if (chi_D(i,j) == 0.0)
              phi(i,j) += R(i,j) /
              (DEl2_center_coefficient/(phi.dx1()*phi.dx1())
               + DEl2_center_coefficient/(phi.dx2()*phi.dx2()) );
    	}
    	
        for (int i=0; i<phi.n1; i++)
        for (int j=0; j<phi.n2; j++)
            R(i,j) = f(i,j)-DEL2(phi,i,j); 
    }
}







int main()
{
    std_setup();
    
    int n1=200 ,n2=300;
    grid2D<double,double > grid( n1,-5,5, n2,-5,5 );
    grid = 0.0;
    
    grid2D<double,double > f(grid);
    f = 0.0;
    
    grid2D<double,double > phi(grid);
    phi = 0.0;
    
    grid2D<double,double > chi_D(grid);
    chi_D = 0.0;
    
    chi_D(n1/2,n2/2) = 1.0;
    phi(n1/2,n2/2) = 10.0;
    chi_D(n1/2+1,n2/2) = 1.0;
    phi(n1/2+1,n2/2) = 10.0;
    chi_D(n1/2,n2/2+1) = 1.0;
    phi(n1/2,n2/2+1) = 10.0;
    chi_D(n1/2+1,n2/2+1) = 1.0;
    phi(n1/2+1,n2/2+1) = 10.0;    
    
    
    for (int i=0; i<n1; i++)
    {
        chi_D(i,0) = 1.0;
        phi(i,0) = 0.0;
        chi_D(i,1) = 1.0;
        phi(i,1) = 0.0;
        chi_D(i,2) = 1.0;
        phi(i,3) = 0.0;
        
        chi_D(i,n2-1) = 1.0;
        phi(i,n2-1) = 0.0;
        chi_D(i,n2-2) = 1.0;
        phi(i,n2-2) = 0.0;
        chi_D(i,n2-3) = 1.0;
        phi(i,n2-3) = 0.0;         
    }
    
    for (int i=0; i<n2; i++)
    {
        chi_D(0,i) = 1.0;
        phi(0,i) = 0.0;
        chi_D(1,i) = 1.0;
        phi(1,i) = 0.0;
        chi_D(2,i) = 1.0;
        phi(2,i) = 0.0;
        
        chi_D(n1-1,i) = 1.0;
        phi(n1-1,i) = 0.0;
        chi_D(n1-2,i) = 1.0;
        phi(n1-2,i) = 0.0;
        chi_D(n1-3,i) = 1.0;
        phi(n1-3,i) = 0.0;
    }
    
    solve_poisson2( phi, f, chi_D, 2,  1E-6 );
    
    for (int i=0; i<n1; i++)
    {
        cout << phi(i,n2/2) << endl;
    }
    
    grid2D<double,double > D2phi(phi);
    D2phi = 0.0;
    
    Del2_FD(
	    phi,
	    D2phi,
	    10,
	    GRID_NONPER_BD_1,
	    10 );
    
	sprintf(fname,"%s/out/phi.png",pwd);
    plotGrid2D_1(phi,fname,cmap);
    
	sprintf(fname,"%s/out/Del2_phi.png",pwd);
    plotGrid2D_1(D2phi,fname,cmap);
}









