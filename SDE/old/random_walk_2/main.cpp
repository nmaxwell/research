
#include <mathlib/math/random/ml_random.h>
#include <mathlib/tools/graphing/plot2D.h>

#include <mathlib/link.cpp>

#include <boost/math/special_functions/erf.hpp>

class gaussian : functor<double, double >
{
public:
    double mu,sig;
    
    gaussian(double mu, double sig ):mu(mu),sig(sig) {};
    
    double operator() (double const & x) const
    {
        return exp(-(x-mu)*(x-mu)/(2.0*sig*sig))/(sqrt2pi*sig);
    }
};

int main()
{
    std_setup();
    
    if (0)
    {
        plot2D plot(-4.0,4.0,-4.0,4.0, 600,600 );
        plot = c3_creme;
        plot.axes_1(0,0, 1, c3_creme,c3_black,plot2D_stdFormat,c3_black, 10,10,1,10  );
        
        
       // plot.plot( plot_schauder ,c3_blue);
        
        sprintf(fname,"%s/out/out.png",pwd);
        plot.png(fname);    
    }
    
    if(1)
    {
        int n_steps = 1000;
        float stop_time = 1.0;
        float time_step = stop_time/n_steps;
        float x_step = 0.1;
        
        for (int trial =0; trial<100; trial++)
        {
            double W[n_steps];
            double t[n_steps];
            
            for (int step = 0; step<n_steps; step++)
            {
                t[step] = time_step*step;
                W[step] = 0;
            }
            
            ml_random rng;
            
            int n_series = (n_steps/32+1);
            bool series[n_series*32];
            
            int series32[n_series];
            for (int k=0; k<n_series; k++)
                series32[k] = rng.gen_int();
            
            int32_to_bool_decomp( n_series, series32, series );
            
            W[0] = 0.0;
            
            for (int step = 1; step<n_steps; step++)
                if ( series[step] )
                    W[step] = W[step-1] + x_step;
                else
                    W[step] = W[step-1] - x_step;
            
                 
            
              
            plot2D plot(-0.1,stop_time*1.1,-10,10, 900,600 );
            plot = c3_creme;
            plot.axes_1(0,-7, 1, c3_creme,c3_black,plot2D_stdFormat,c3_black, 0.1,10,1,10  );
            
            plot.plot<double, double >(&(t[0]),&(W[0]),n_steps ,c3_black);
            for (int step = 0; step<n_steps; step++)
                W[step] = W[step]*W[step]-t[step];
            plot.plot<double, double >(&(t[0]),&(W[0]),n_steps ,c3_red);
                
            sprintf(fname,"%s/out/%03d.png",pwd,trial);
            plot.png(fname);
        }
    }
    
    
}




