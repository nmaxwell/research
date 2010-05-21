
#include <mathlib/math/random/ml_random.h>
#include <mathlib/tools/graphing/plot2D.h>

#include <mathlib/link.cpp>

#include <boost/math/special_functions/erf.hpp>

const double sqrt12 = sqrt(12.0);



double N_RV(double omega,double mu = 0.0,double sigma = 1)
{        
    return sigma*sqrt2*boost::math::erf_inv(2.0*omega-1.0)+0.0;
}


double phi( double x )
{
    if ( x >=0 && x < 0.5 )
        return x*sqrt12;
    if ( x >=0.5 && x < 1.0 )
        return (1.0-x)*sqrt12;
    return 0;
    
}

double phi( double x, int k, int j )
{
    return phi( x*pow(2.0,j) -k )*pow(sqrt2,j);
}

class BM
{
public:
    double * A;
    int j_max;
    
    BM(int j_max):A(0),j_max(j_max) {
    
        A = new double [2*(2<<j_max)]; 
        ml_random rng;
        for (int k = 0; k<= (2*(2<<j_max-1)) ; k++)
            A[k] = N_RV(rng.gen_double(),0,1 );
        
    };
    
    ~BM() { if (A) delete [] A; A = 0; }
    
    double operator() (double x) const {
    
        double sum = 0;
        int i=0;
        
        for (int j=0; j<=j_max; j++)
        for (int k=0; k<(1<<j); k++)
        {
            //cout << A[i++] << endl;
            sum += A[i++]*phi(x,k,j);
        }
        
        return sum*pow(2,-0.5*j_max); 
    }
};


int main()
{
    std_setup();
    
    for (int j_max=0; j_max<=20; j_max++)
    {   
        BM bm(j_max);
                
        plot2D plot(-0.1, 1.1,-5,5, 900,600 );
        plot = c3_creme;
        plot.axes_1(0,-7, 1, c3_creme,c3_black,plot2D_stdFormat,c3_black, 0.1,10,1,10  );
        
        plot.plot ( bm ,c3_black);
        //bm(0);
        
        sprintf(fname,"%s/out/%03d.png",pwd,j_max);
        plot.png(fname);
        
    }
}








