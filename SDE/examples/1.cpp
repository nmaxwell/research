
#define ML_NOT_USE_FFTW_MALLOC

#include <mathlib/math/random/ml_random.h>
#include <mathlib/math/SDE/BM.h>
#include <mathlib/math/grids/grid1D.h>
#include <mathlib/math/grids/grid2D.h>
//#include <mathlib/math/grids/extra/plots.cpp>

#include <mathlib/link.cpp>
#include <mathlib/non-lib_things.h>

//#include <boost/math/special_functions/erf.hpp>




template< class T >
void output( T * data, int n, const char * fname, const char * delim1="\t", const char * delim2="\n" )
{
    ofstream out;
    out.open(fname, fstream::out);
    assert(out.good() && out.is_open());
    
    for ( int k=0; k<n; k++ )
        out << data[k] << delim2;
    
    out.close();
}

template< class T1, class T2 >
void output( T1 * data1, T2 * data2, int n, const char * fname, const char * delim1="\t", const char * delim2="\n" )
{
    ofstream out;
    out.open(fname, fstream::out);
    assert(out.good() && out.is_open());
    
    for ( int k=0; k<n; k++ )
        out << data1[k] << delim1 << data2[k] << delim2;
    
    out.close();
}


/*



import math
import pylab as p
import sys
sys.path.append("/workspace/mathlib/tools/python/")
from read_file import *


print len(read_file( "out_0" ))

for n in range(1,10):
    print p.plot( read_file( "out_" + str(n) ) )

p.show()




*/



//     gen_BM( double h, int n_steps, float **& W, int n_runs, int mode )
//     gen_BM( double h, int n_steps, float *& W, int mode )


int main()
{
    std_setup();
    
    int n_runs = 10;   //1E2*32;
    int n_steps = 1000;
    double **W=0;
    
    double stop_time = 1.0;
    double dt = stop_time/n_steps;
    
    double t1,t2;
    
    t1 = get_real_time();
    
    gen_BM( dt, n_steps, W, n_runs, BM_mode_3 );
    
    t2 = get_real_time();
    
    cout << "time: " << t2-t1 << endl;
    
    for (int run=0; run<n_runs; run++)
    {
        sprintf(fname, "/workspace/output/temp/out_%d", run);
        output( W[run], n_steps, fname );
    }
    
    ml_free( W, n_runs );
    
    std_exit();
}



