


#include <mathlib/math/std_math.h>
#include <mathlib/math/transforms/fft.h>

#include <mathlib/link.cpp>
#include <mathlib/non-lib_things.h>


// only a mockup...

void direct_convolution( double * result, double * kernel, double * data, int n_data, int n_kernel )
{
    for (int k=0; k<n_data; k++)
    {
        result[k]=0;
        for (int j=0; j<n_kernel; j++)
            result[k] += data[ (k+j)%n_data ]*kernel[j];
    }
}



int main()
{
    std_setup();
    
    int n_trials = 20;
    double T1,T2,T3;
    double rate1=0,rate2;
    
    int n_kernel = 4;
    int n_data = 1024*64;
    
    int p=13;
    int q=13;
    for (p=1; p<14; p++)
    for (q=1; q<14; q++)
    {
        n_data = pow(2,p);
        n_kernel = pow(2,q);
        
        cout << "W:\t\t"  << n_kernel << endl;
        cout << "N:\t\t"  << n_data << endl;
        cout << "logN:\t\t"  << log(n_kernel) << endl;
        
        double * data= ml_alloc<double> (n_data);
        double * kernel= ml_alloc<double> (n_data);
        double * result= ml_alloc<double> (n_data);
        
        complex<double> * fft_kernel = ml_alloc<complex<double> > (n_data);
        complex<double> * fft_data = ml_alloc<complex<double> > (n_data);
        
        fft(kernel, fft_kernel, n_data);
        ifft(fft_kernel, result, n_data);
        
        
        
        T1 = get_real_time();
        for (int trial=0; trial<n_trials; trial++)
        {
            fft(data, fft_data, n_data);
            
            for(int k=0; k<n_data/2+1; k++)
                fft_data[k] *= fft_kernel[k];
            
            ifft(fft_kernel, result, n_data);
        }
        T2 = get_real_time();
        rate1= (T2-T1)/n_trials;
        cout << "fft: \t\t" << rate1 << "\t" << 1E9*rate1/n_data << endl;
        
        T1 = get_real_time();
        for (int trial=0; trial<n_trials; trial++)
        {
            direct_convolution( result, kernel, data, n_data, n_kernel );
        }
        T2 = get_real_time();
        rate2= (T2-T1)/n_trials;
        cout << "direct: \t" << rate2 << "\t" <<  1E9*rate2/(n_data*n_kernel) << endl;
        
        cout << "ratio: \t" << (rate1)/(rate2/n_kernel) << endl;
        
        
        cout << endl;
    }
    
    std_exit();
}




