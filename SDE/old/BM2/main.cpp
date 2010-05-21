
#include <mathlib/math/random/ml_random.h>
#include <mathlib/tools/graphing/plot2D.h>

#include <mathlib/link.cpp>

#include <boost/math/special_functions/erf.hpp>



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
    
    
    double stop_time = 1.0;
    int j = 3;
    
    for (j=0; j<=7; j++)
    {
        int k_max = stop_time*(1 << j);
        
        int n_steps = 400;
        
        double W[n_steps];
        double t[n_steps];
        
     /*   ml_random rng;
        
        for (int step = 0; step<n_steps; step++)
        {
            t[step] = (stop_time*step)/n_steps;
            
            double sum = 0;
            for (int k=0; k<=k_max; k++)
                sum += N_RV(rng.gen_double())*phi(t[step],k,j);
            
            W[step] = sum;
        }
        */
        
                
        plot2D plot(-0.1,stop_time*1.1,-10,10, 900,600 );
        plot = c3_creme;
        plot.axes_1(0,-7, 1, c3_creme,c3_black,plot2D_stdFormat,c3_black, 0.1,10,1,10  );
        
        //plot.plot<double, double >(&(t[0]),&(W[0]),n_steps ,c3_black);
        
        sprintf(fname,"%s/out/.png",pwd,j);
        plot.png(fname);
        
    }
    
    

}








