
#include <mathlib/math/std_math.h>
#include <mathlib/tools/graphing/plot2D.h>

#include <mathlib/math/random/ml_random.h>
#include <mathlib/math/grids/grid2D.h>
//#include <mathlib/math/grids/grid2D_Del2_FD.h>
#include <mathlib/math/grids/extra/plots.cpp>
#include <mathlib/math/PDE/laplacian_hdaf.h>

#include <mathlib/link.cpp>
#include <mathlib/non-lib_things.h>

#include "Polygon2D.cpp"








bool straddle_test( float const & x1, float const & x2, float const & x3, float const & x4, float const & y1, float const & y2, float const & y3, float const & y4 )
{
    if ( max(x1,x2) < min(x3,x4) )
        return false;
    if ( max(x3,x4) < min(x1,x2) )
        return false;
    if ( max(y1,y2) < min(y3,y4) )
        return false;
    if ( max(y3,y4) < min(y1,y2) )
        return false;
    
    float z1 = (x3-x1)*(y2-y1)-(y3-y1)*(x2-x1);
    float z2 = (x4-x1)*(y2-y1)-(y4-y1)*(x2-x1);
    
    //if ( z1==0.f && z2==0.f ) return true;
    
    if ( (z1<0) xor (z2<0) ) return true;
    
    return false;
}

float distance_squared( float const & X11, float const & X12, float const & X21, float const & X22, float const & Y1, float const & Y2 )
{
	float y1 = Y1 - X11;
	float y2 = Y2 - X12;
	float x1 = X21 - X11;
	float x2 = X22 - X12;
	
	float dot = x1*y1+x2*y2;
	
	if ( dot < 0.f )
	    return y1*y1+y2*y2;
	
	float lx2 = x1*x1+x2*x2;
	
	if ( dot > lx2 )
		return (Y1-X21)*(Y1-X21)+(Y2-X22)*(Y2-X22);
    
	return y1*y1+y2*y2 - dot*dot/lx2;
}


void drift( float & x, float & y, float t )
{
    //x = cos(_2pi*t)*0.05;
    //y = sin(_2pi*t)*0.05;
    
    x = .04;
    y = .04;
    
   // x = y = 0;
}





int hitting_distribution_nodrift( double *distribution,  float x, float y, Polygon2D &boundary, double dt, int n_runs, int &failed_runs )
{
    static ml_random rng;
    
    int n_bd = boundary.n_points-1;
    float sqrt_dt = sqrt(dt);
    int max_steps = 1E7;
    int max_bad_runs = 100;
    float dW1,dW2,W1,W2;
    
    for (int k=0; k<n_bd; k++)
        distribution[k] = 0;
    
    int bad_runs = 0;
    int run = 0;
    double sum = 0.0;
    double time = 0;
    
    while ( run < n_runs )
    {
        int hitting_index = -1;
        int step = 0;
        
        W1 = x;
        W2 = y;
        
        while( true )
        {
            rng.std_normal_rv(dW1,dW2);
            
            if ( boundary.interior_test(W1+dW1*sqrt_dt, W2+dW2*sqrt_dt) )
            {
            	W1 += dW1*sqrt_dt;
            	W2 += dW2*sqrt_dt;
            }
            else
            {
            	float nearest = 1.E10;
                hitting_index = 0;
                
                for (int k=0; k<n_bd; k++ )
                {
                	float d2 = distance_squared( boundary.x_point(k), boundary.y_point(k), boundary.x_point(k+1), boundary.y_point(k+1), W1, W2 );
                	
                	if ( d2 < nearest )
                	{
                		hitting_index = k;
                		nearest = d2;
                	}
                }
                
                break;
            }
            
            step++;
            time += dt;
            
            if ( step == max_steps )
                break;
        }
        
        if ( hitting_index >=0 and hitting_index<n_bd )
        {
            run++;
            distribution[hitting_index] += 1.0/n_runs;
        }
        else
        {
        	bad_runs++;
        	
        	if ( bad_runs == max_bad_runs )
        	{
        		failed_runs += bad_runs;
        		return 1;
        	}
        }
    }
    
    failed_runs += bad_runs;
    return 0;
}


int hitting_distribution( double *distribution,  float x, float y, Polygon2D &boundary, double dt, int n_runs, int &failed_runs )
{
    static ml_random rng;
    
    int n_bd = boundary.n_points-1;
    float sqrt_dt = sqrt(dt);
    int max_steps = 1E7;
    int max_bad_runs = 100;
    float dW1,dW2,W1,W2, V1,V2;
    
    for (int k=0; k<n_bd; k++)
        distribution[k] = 0;
    
    int bad_runs = 0;
    int run = 0;
    double sum = 0.0;
    double time = 0;
    double M=0;
    
    while ( run < n_runs )
    {
        int hitting_index = -1;
        int step = 0;
        
        W1 = x;
        W2 = y;
        M=0;
        
        while( true )
        {
            rng.std_normal_rv(dW1,dW2);
            
            drift( V1, V2, time );
            
            if ( boundary.interior_test(W1+V1*dt+dW1*sqrt_dt, W2+V2*dt+dW2*sqrt_dt) )
            {
            	W1 += V1*dt+ dW1*sqrt_dt;
            	W2 += V2*dt+ dW2*sqrt_dt;
                
                M += -(V1*dW1+V2*dW2) + -.5*(V1*V1+V2*V2)*dt;
            }
            else
            {
            	float nearest = 1.E10;
                hitting_index = 0;
                
                for (int k=0; k<n_bd; k++ )
                {
                	float d2 = distance_squared( boundary.x_point(k), boundary.y_point(k), boundary.x_point(k+1), boundary.y_point(k+1), W1, W2 );
                	
                	if ( d2 < nearest )
                	{
                		hitting_index = k;
                		nearest = d2;
                	}
                }
                
                break;
            }
            
            step++;
            time += dt;
            
            if ( step == max_steps )
                break;
        }
        
        if ( hitting_index >=0 and hitting_index<n_bd )
        {
            run++;
            distribution[hitting_index] += exp(M)/n_runs;
        }
        else
        {
        	bad_runs++;
        	
        	if ( bad_runs == max_bad_runs )
        	{
        		failed_runs += bad_runs;
        		return 1;
        	}
        }
    }
    
    failed_runs += bad_runs;
    return 0;
}



int hitting_distribution_drift_nochange( double *distribution,  float x, float y, Polygon2D &boundary, double dt, int n_runs, int &failed_runs )
{
    static ml_random rng;
    
    int n_bd = boundary.n_points-1;
    float sqrt_dt = sqrt(dt);
    int max_steps = 1E7;
    int max_bad_runs = 100;
    float dW1,dW2,W1,W2, V1,V2;
    
    for (int k=0; k<n_bd; k++)
        distribution[k] = 0;
    
    int bad_runs = 0;
    int run = 0;
    double sum = 0.0;
    double time = 0;
    
    while ( run < n_runs )
    {
        int hitting_index = -1;
        int step = 0;
        
        W1 = x;
        W2 = y;
        
        while( true )
        {
            rng.std_normal_rv(dW1,dW2);
            
            drift( V1, V2, time );
            
            if ( boundary.interior_test(W1+V1*dt+dW1*sqrt_dt, W2+V2*dt+dW2*sqrt_dt) )
            {
            	W1 += V1*dt+ dW1*sqrt_dt;
            	W2 += V2*dt+ dW2*sqrt_dt;
            }
            else
            {
            	float nearest = 1.E10;
                hitting_index = 0;
                
                for (int k=0; k<n_bd; k++ )
                {
                	float d2 = distance_squared( boundary.x_point(k), boundary.y_point(k), boundary.x_point(k+1), boundary.y_point(k+1), W1, W2 );
                	
                	if ( d2 < nearest )
                	{
                		hitting_index = k;
                		nearest = d2;
                	}
                }
                
                break;
            }
            
            step++;
            time += dt;
            
            if ( step == max_steps )
                break;
        }
        
        if ( hitting_index >=0 and hitting_index<n_bd )
        {
            run++;
            distribution[hitting_index] += 1.0/n_runs;
        }
        else
        {
        	bad_runs++;
        	
        	if ( bad_runs == max_bad_runs )
        	{
        		failed_runs += bad_runs;
        		return 1;
        	}
        }
    }
    
    failed_runs += bad_runs;
    return 0;
}



void boundary_def(float & x, float & y, double t)
{
    // t \in [0,1]
    
    //x = cos(t*_2pi)*(cos(t*_2pi*5)*2.5+6);
    //y = sin(t*_2pi)*(cos(t*_2pi*5)*2.5+6);
    
    x = cos(t*_2pi)*5.0;
    y = sin(t*_2pi)*5.0;
}



int main()
{
    std_setup();
    
    
    int n_bd=300;
    Polygon2D boundary;
    
    {
        // setup boundary
        boundary.init(n_bd+1);
        
        for (int k=0; k<n_bd; k++)
        {
            double t1 = double(k)/n_bd;
            double t2 = double(k+1)/n_bd;
            
            boundary_def( boundary.x_point(k),boundary.y_point(k), t1 );
            boundary_def( boundary.x_point(k+1),boundary.y_point(k+1), t2 );
        }
    }
    
    double *distribution = ml_alloc<double> (n_bd);
    double *distribution_nodrift = ml_alloc<double> (n_bd);
    double *distribution_drift_nochange = ml_alloc<double> (n_bd);
    int n_runs=10000;
    double dt = 0.1;
    int failed_runs=0;
    
    float x=0;
    float y=0;
    
    hitting_distribution( distribution, x, y, boundary, dt, n_runs, failed_runs );
    hitting_distribution_nodrift( distribution_nodrift, x, y, boundary, dt, n_runs, failed_runs );
    hitting_distribution_drift_nochange( distribution_drift_nochange, x, y, boundary, dt, n_runs, failed_runs );
    
    //for (int k=0; k<n_bd; k++)
    //    cout << k << "\t" << (double)distribution[k]  << "\t" << (double)distribution_nodrift[k] << endl;
    
    cout << "failed_runs: " << failed_runs << endl;
    
    output( distribution, n_bd, "/workspace/output/SDE/project1/hitting_distribution/distribution" );
    output( distribution_nodrift, n_bd, "/workspace/output/SDE/project1/hitting_distribution/distribution_nodrift" );
    output( distribution_drift_nochange, n_bd, "/workspace/output/SDE/project1/hitting_distribution/distribution_drift_nochange" );
    
    
    int N = 512;
    plot2D plt( -10.0, 10.0, -10.0, 10.0, N, N);
    plt = ml_white;
    
    for ( int i=0; i<N; i++ )
    for ( int j=0; j<N; j++ )
    {
        if ( boundary.interior_test(plt.to_x(i), plt.to_y(j)) )
            plt(i,j) = ml_blue;
    }
    
    for (int k=0; k<n_bd; k++)
        plt.ptLine( boundary.x_point(k),boundary.y_point(k), boundary.x_point(k+1),boundary.y_point(k+1), ml_red );
    
    plt.png( "/workspace/output/SDE/project1/hitting_distribution/interiod_test.png" );
    
    
    std_exit();
}


