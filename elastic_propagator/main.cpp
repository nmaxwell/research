
/*

mencoder 'mf://*.png' -mf  fps=25:type=png -ovc copy -oac copy -o output.avi

*/

// #define FFTW_PLAN_MODE FFTW_PATIENT

#define N_FFT_THREADS  3

//#define ML_TRACK_MEMORY
//#define ML_FFT_3D_R2C_TELL

#include <mathlib/math/std_math.h>

#include <mathlib/math/grids/grid3D.h>
#include <mathlib/math/grids/grid_file.cpp>
#include <mathlib/math/grids/extra/plots.cpp>

#include <mathlib/math/PDE/grad_hdaf.h>

#include <mathlib/link.cpp>
#include <mathlib/non-lib_things.h>



#include "3d.h"


void output(double t, grid_vec & u, grid_vec & v)
{
	static int n = 0;
    
    sprintf(fname,"/workspace/output/elastic_propagate_3d/out_dat/u1_%05d.dat", n);
	assert(!writeFile(u._1,fname));
    
    sprintf(fname,"/workspace/output/elastic_propagate_3d/out_dat/u2_%05d.dat", n);
	assert(!writeFile(u._2,fname));
    
    sprintf(fname,"/workspace/output/elastic_propagate_3d/out_dat/u3_%05d.dat", n);
	assert(!writeFile(u._3,fname));
	
	n++;
}




class C1_func : public functor3<double,double >
{
public:
	double operator() (double const & x, double const & y, double const & z) const 
	{
        return 1.0;
	}
};

class C2_func : public functor3<double,double >
{
public:
	double operator() (double const & x, double const & y, double const & z) const 
	{
        return 0.0;
	}
};

class C3_func : public functor3<double,double >
{
public:
	double operator() (double const & x, double const & y, double const & z) const 
	{
        return 0.0;
	}
};

class C4_func : public functor3<double,double >
{
public:
	double operator() (double const & x, double const & y, double const & z) const 
	{
        return 0.0;
	}
};

class C5_func : public functor3<double,double >
{
public:
	double operator() (double const & x, double const & y, double const & z) const 
	{
        return 0.0;
	}
};

class C6_func : public functor3<double,double >
{
public:
	double operator() (double const & x, double const & y, double const & z) const 
	{
        return 0.0;
	}
};

class C7_func : public functor3<double,double >
{
public:
	double operator() (double const & x, double const & y, double const & z) const 
	{
        return 1.0;
	}
};

class C8_func : public functor3<double,double >
{
public:
	double operator() (double const & x, double const & y, double const & z) const 
	{
        return 0.0;
	}
};

class C9_func : public functor3<double,double >
{
public:
	double operator() (double const & x, double const & y, double const & z) const 
	{
        return 0.0;
	}
};

class C10_func : public functor3<double,double >
{
public:
	double operator() (double const & x, double const & y, double const & z) const 
	{
        return 0.0;
	}
};
class C11_func : public functor3<double,double >
{
public:
	double operator() (double const & x, double const & y, double const & z) const 
	{
        return 0.0;
	}
};

class C12_func : public functor3<double,double >
{
public:
	double operator() (double const & x, double const & y, double const & z) const 
	{
        return 1.0;
	}
};

class C13_func : public functor3<double,double >
{
public:
	double operator() (double const & x, double const & y, double const & z) const 
	{
        return 0.0;
	}
};

class C14_func : public functor3<double,double >
{
public:
	double operator() (double const & x, double const & y, double const & z) const 
	{
        return 0.0;
	}
};

class C15_func : public functor3<double,double >
{
public:
	double operator() (double const & x, double const & y, double const & z) const 
	{
        return 0.0;
	}
};

class C16_func : public functor3<double,double >
{
public:
	double operator() (double const & x, double const & y, double const & z) const 
	{
        return 0.0;
	}
};

class C17_func : public functor3<double,double >
{
public:
	double operator() (double const & x, double const & y, double const & z) const 
	{
        return 0.0;
	}
};

class C18_func : public functor3<double,double >
{
public:
	double operator() (double const & x, double const & y, double const & z) const 
	{
        return 0.0;
	}
};

class C19_func : public functor3<double,double >
{
public:
	double operator() (double const & x, double const & y, double const & z) const 
	{
        return 0.0;
	}
};

class C20_func : public functor3<double,double >
{
public:
	double operator() (double const & x, double const & y, double const & z) const 
	{
        return 0.0;
	}
};

class C21_func : public functor3<double,double >
{
public:
	double operator() (double const & x, double const & y, double const & z) const 
	{
        return 0.0;
	}
};




class u1_func : public functor3<double,double >
{	
public:
	double operator() (double const & x, double const & y, double const & z) const 
	{
		return 0;
	}
};

class u2_func : public functor3<double,double >
{	
public:
	double operator() (double const & x, double const & y, double const & z) const 
	{
		return 0;
	}
};

class u3_func : public functor3<double,double >
{	
public:
	double operator() (double const & x, double const & y, double const & z) const 
	{
		return 0;
	}
};



class driving_func_1 : public functor3<double,double >
{
public:
	double t;
	driving_func_1(double t):t(t) {}
	
public:
	double operator() ( double const & x, double const & y, double const & z) const 
	{
        double x1 = x - 0;
        double x2 = y - 0;
        double x3 = z - .5;
        
		double s = sqrt(x1*x1+x2*x2+x3*x3);
        double a = 0.5;
        
		if (s <= 0.5*a)
			return (cos(ml_2pi*s/a)+1.0)*cos(_2pi*t+0*_2pi);
		else return 0.0;
	}
};

class driving_func_2 : public functor3<double,double >
{
public:
	double t;
	driving_func_2(double t):t(t) {}
	
public:
	double operator() ( double const & x, double const & y, double const & z) const 
	{
        double x1 = x - .5;
        double x2 = y - 0;
        double x3 = z - 0;
        
		double s = sqrt(x1*x1+x2*x2+x3*x3);
        double a = 0.5;
        
		if (s <= 0.5*a)
			return (cos(ml_2pi*s/a)+1.0)*cos(_2pi*t+.3*_2pi);
		else return 0.0;
	}
};

class driving_func_3 : public functor3<double,double >
{
public:
	double t;
	driving_func_3(double t):t(t) {}
	
public:
	double operator() ( double const & x, double const & y, double const & z) const 
	{
        double x1 = x - 0;
        double x2 = y - .5;
        double x3 = z - 0;
        
		double s = sqrt(x1*x1+x2*x2+x3*x3);
        double a = 0.5;
        
		if (s <= 0.5*a)
			return (cos(ml_2pi*s/a)+1.0)*cos(_2pi*t+.6*_2pi);
		else return 0.0;
	}
};



class damping_func : public functor3<double,double >
{
	
public:
	double operator() ( double const & x, double const & y, double const & z ) const 
	{
        double ax = 0.4;
        double ay = 0.4;
        double az = 0.4;
        
        return  75.0*( 1.0/(cosh((x-10)/ax)) +  1.0/(cosh((x+10)/ax)) + 1.0/(cosh((y-10)/ay)) +  1.0/(cosh((y+10)/ax)) + 1.0/(cosh((z-10)/az)) +  1.0/(cosh((z+10)/az))  );
	}
};


int main()
{
    std_setup();
    
    int n = 128;
    double x0 = 10; 
    int expansion_order = 5;
    int hdaf_order = 12;
    double hdaf_gamma = 0.6;
    
    double tf = 30.0;
    double dt = 0.01;
    
    grid G( n,-x0,x0, n,-x0,x0, n,-x0,x0 );
    grid_vec u,v;
    
    G = 0;
    
    u._1 = G;
    u._2 = G;
    u._3 = G;
    v._1 = G;
    v._2 = G;
    v._3 = G;
    
    grid C[21];
    for (int k=0; k<21; k++ )
        C[k] = G;
    
    C[0] = C1_func();
    C[1] = C2_func();
    C[2] = C3_func();
    C[3] = C4_func();
    C[4] = C5_func();
    C[5] = C6_func();
    C[6] = C7_func();
    C[7] = C8_func();
    C[8] = C9_func();
    C[9] = C10_func();
    C[10] = C11_func();
    C[11] = C12_func();
    C[12] = C13_func();
    C[13] = C14_func();
    C[14] = C15_func();
    C[15] = C16_func();
    C[16] = C17_func();
    C[17] = C18_func();
    C[18] = C19_func();
    C[19] = C20_func();
    C[20] = C21_func();
    
    u._1 = u1_func();
    u._2 = u2_func();
    u._3 = u3_func();
    v = 0.0;
    
    anisotropic_propagator P;
    P.init( G, expansion_order, C, hdaf_order, hdaf_order, hdaf_order, hdaf_gamma, hdaf_gamma, hdaf_gamma );
    
    driver F;
    F.init( G, expansion_order, C, hdaf_order, hdaf_order, hdaf_order, hdaf_gamma, hdaf_gamma, hdaf_gamma );
    
    grid_vec f1,f2;
    f1._1 = G;
    f1._2 = G;
    f1._3 = G;
    f2._1 = G;
    f2._2 = G;
    f2._3 = G;
    
    grid D;
    D = G;
    D = damping_func();
    
    for (int i=0; i<G.n1; i++ )
    for (int j=0; j<G.n2; j++ )
    for (int k=0; k<G.n3; k++ )
		D(i,j,k) = exp(-dt*D(i,j,k));    
    
    double t = 0.0;
    
    double T1,T2,T0 = get_real_time();
    
    while (t <= tf)
    {
        double mag = L2norm(u._1) + L2norm(u._2) + L2norm(u._3) + L2norm(v._1) + L2norm(v._2) + L2norm(v._3);     
        
        if ( mag > 1E10 ) 
        {
            cout << "Exploded; maginitude is " << mag << ", quitting..." << endl;
            exit(0);
        }
        
        cout << "Step: " << t << endl;
        
        cout << "Current magnitude: " << mag << endl;
        
        T1 = get_real_time();
        
        output(t,u,v);
        
        f1._1 = driving_func_1(t);
        f1._2 = driving_func_2(t);
        f1._3 = driving_func_3(t);
        
        f2._1 = driving_func_1(t+dt);
        f2._2 = driving_func_2(t+dt);
        f2._3 = driving_func_3(t+dt);
        
        T2 = get_real_time();
        
        cout << "Total compute time: " << T2-T0 << "\tother: " << T2-T1 << endl;
                
        T1 = get_real_time();
        
        P( dt, u,v, u,v );
        
        T2 = get_real_time();
        
        cout << "Total compute time: " << T2-T0 << "\tpropagate: " << T2-T1 << endl;
        
        T1 = get_real_time();
        
        F( dt, f1, f2, u,v );
        
        u *= D;
        v *= D;
        
        T2 = get_real_time();
        
        cout << "Total compute time: " << T2-T0 << "\tforce & damp: " << T2-T1 << endl;
        
        t += dt;
        
        cout << endl;
    }
    
    output(t,u,v);
}
