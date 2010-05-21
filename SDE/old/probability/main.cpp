
#include <mathlib/math/random/ml_random.h>
#include <mathlib/tools/graphing/plot2D.h>

#include <mathlib/link.cpp>

#include <boost/math/special_functions/erf.hpp>


double mu = 0.0;
double sigma = 1.0;

double X(double omega)
{        
    return sigma*sqrt2*boost::math::erf_inv(2.0*omega-1.0)+0.0;
    //return omega;
}

double N(double x)
{
    return exp(-(x-mu)*(x-mu)/(sigma*sigma*2))/(sqrt2pi*sigma);
}


int main()
{
    std_setup();
    
    int n_trials = 1E6;
    int n = 40;
    
    for ( n_trials = 10; n_trials <=1E7; n_trials *= 2)
    {        
        int64 bins[n+1];
        double breaks[n];
        
        // Divide R into n+2 intervals,
        // I_0 = (-infty,break[0] ]
        // I_1 = (break[0], break[1] ]
        // I_i = (break[i-1], break[i] ]
        // I_(n-1) = (break[n-2], break[n-1] ]
        // I_n = (break[n-1], +infty )
        
        for (int i=0; i<n; i++)
            breaks[i] = (2.0*double(i)/double(n)-1.0)*4.0;
        
        for (int k=0; k<n+1; k++)
            bins[k] = 0;
        
        ml_random rng;
        
        for (int j=0; j<n_trials; j++)
        {
            double omega = rng.gen_double();            
            double X_w = X(omega);
            
           // cout << omega << "\t" << X_w << endl;
            
            int k = 0; // put in kth bin, corresponding ot kth interval
            
            if ( X_w <= breaks[0] ) bins[0]++;
            
            for ( k=1; k<n; k++)
                if ( breaks[k-1] < X_w && X_w <= breaks[k] )
                    bins[k]++;
            
            if ( X_w > breaks[n-1] ) bins[n]++;
        }
        
        
        
        cout << 0 << "\t(inf," << breaks[0] << "]\t\t" << (double)(bins[0])/n_trials << endl;
        
        for (int k=1; k<n; k++)
            cout << k << "\t(" << breaks[k-1] << "," << breaks[k] << "]\t\t" << (double)(bins[k])/n_trials << endl;
        
        cout << n << "\t(" << breaks[n-1] << ",inf)\t\t" << (double)(bins[n])/n_trials << endl;
        
        double sum = 0;
        for (int k=0; k<=n; k++)
            sum += (double)(bins[k])/n_trials;
        cout << sum << endl;
        
        plot2D plot(-10,10,-0.1,1.1, 900,600 );
        plot = c3_creme;
        plot.axes_1(0,0, 1, c3_creme,c3_black,plot2D_stdFormat,c3_black, 1,10,0.1,10  );
        
        //plot.ptLine( -1000,(double)(bins[0])/n_trials, breaks[0],(double)(bins[0])/n_trials, c3_red);
        
        for (int k=1; k<n; k++)
        {
            double y = ((double)(bins[k])/n_trials)/(breaks[k]-breaks[k-1]);
            plot.ptLine( breaks[k-1],y, breaks[k],y, c3_red);
        }
        //plot.ptLine( breaks[n],(double)(bins[n])/n_trials, 1000,(double)(bins[n])/n_trials, c3_red);
        
        plot.plot(N,c3_blue);
        
        sprintf(fname,"%s/out/%09d.png",pwd,n_trials);
        plot.png(fname);
        
    }
    
    
    

}








