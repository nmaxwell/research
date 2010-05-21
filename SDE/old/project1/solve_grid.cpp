
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


ml_color color_map_green( double z )
{
    float x = atan((double)z)/pi+0.5;
    return ml_color(0,x,0);
}

ml_color color_map_blue( double z )
{
    float x = atan((double)z)/pi+0.5;
    return ml_color(0,0,x);
}

ml_color color_map_red( double z )
{
    float x = atan((double)z)/pi+0.5;
    return ml_color(x,0,0);
}




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
    //x = 0.2;
    //y = 0.2;
    
    x = 0;
    y = 0;
}


int solve_laplace( double &u, float x, float y, Polygon2D &boundary, double *f, double dt, int n_runs, int &failed_runs )
{
    static ml_random rng;
    
    int n_bd = boundary.n_points-1;
    float sqrt_dt = sqrt(dt);
    int max_steps = 10000;
    int max_bad_runs = 100;
    float dW1,dW2,W1,W2, V1,V2;
    
    u = 0;
    
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
            sum += f[hitting_index]*exp(M);
        }
        else
        {
        	bad_runs++;
        	
        	if ( bad_runs == max_bad_runs )
        	{
        		failed_runs += bad_runs;
        			u = 0.0;
        		
        		return 1;
        	}
        }
    }
    
    failed_runs += bad_runs;
    u = sum/n_runs;
    return 0;
}






struct thread_args
{
    int i;
    grid2D<> u;
    Polygon2D boundary;
    double * f;
    double dt;
    int n_runs;
    
    int eval_count;
    int failed_runs;
};

void *compute_line(void *vargs)
{
    thread_args &args = *((thread_args*)vargs);
    
    grid2D<> &u = args.u;
    int i = args.i;
    
    int N = u.n2;
    
    for (int j=0; j<N; j++)
    {
        double x = u.x1(i);
        double y = u.x2(j);
        
        if ( args.boundary.interior_test(x,y) )
        {
            int err = solve_laplace( u(i,j), x,y, args.boundary, args.f, args.dt, args.n_runs, args.failed_runs );
            args.eval_count++;
        }
    }
    
    return NULL;
}




void boundary_def(float & x, float & y, double & f, double t)
{
    // t \in [0,1]
    
    x = cos(t*_2pi)*(cos(t*_2pi*5)*2.5+6);
    y = sin(t*_2pi)*(cos(t*_2pi*5)*2.5+6);
    
    f = cos(_2pi*t*15)*2.0 + sin(_2pi*t*13)*3.0+0*sin(_2pi*t*60)*7.0 ;
    f *= 2.0;
}



int main()
{
    std_setup();
    
    
    
    
    int n_bd=100;
    Polygon2D boundary;
    double * f = 0;
    
    {
        // setup boundary
        
        f = ml_alloc<double> (n_bd);
        boundary.init(n_bd+1);
        
        for (int k=0; k<n_bd; k++)
        {
            double t1 = double(k)/n_bd;
            double t2 = double(k+1)/n_bd;
            
            boundary_def( boundary.x_point(k),boundary.y_point(k), f[k], t1 );
            boundary_def( boundary.x_point(k+1),boundary.y_point(k+1), f[k+1], t2 );
        }
    }
    
    
    
    
    int n_runs=5;
    
    //for ( n_runs=1; n_runs<=50; n_runs++ )
    {
        
        double dt = 0.01;
        int N = 512;
        
        grid2D<> u( N, -10, 10, N, -10, 10 );
        u = 0.0;
        
        int eval_count = 0;
        int failed_runs = 0;
        
        int n_threads = 8;
        
        thread_args args[n_threads];
        pthread_t pth[n_threads];
        bool threads_created[n_threads];
        
        for (int k=0; k<n_threads; k++)
        {
            args[k].u = u;
            args[k].boundary = boundary;
            args[k].f = ml_alloc<double> (n_bd);
            for (int j=0; j<n_bd; j++)
                args[k].f[j] = f[j];
            args[k].dt = dt;
            args[k].n_runs = n_runs;
            
            args[k].eval_count =0;
            args[k].failed_runs =0;
            
            threads_created[k] = false;
        }
        
        
        int thread=0;
        for (int i=0; i<N; i++)
        {
            thread = thread%n_threads;
            if (threads_created[thread])
            {
                int err = pthread_join( pth[thread], NULL );
                if (err) { cout << "pthread_join error: " << err << "\n"; exit(0); }
                threads_created[thread] = false;
            }
            args[thread].i = i;
            cout << i << "\t" << thread << endl;
            
            int err = pthread_create( &(pth[thread]), NULL, compute_line, (void*) &(args[thread]) );
            if (err ) { cout << "pthread_create error: " << err << "\n"; exit(0); }
            threads_created[thread] = true;
            thread ++;
        }
        
        for (thread=0; thread<n_threads; thread++)
        {
            if (threads_created[thread])
            {
                int err = pthread_join( pth[thread], NULL );
                if (err) { cout << "pthread_join error: " << err << "\n"; exit(0); }
            }
            
            u += args[thread].u;
            
            eval_count += args[thread].eval_count;
            failed_runs += args[thread].failed_runs;
        }
        
        
        int n_smooth = 0;
        grid2D<> u_smooth( u );
        
        for ( int i=0; i<u.n1; i++ )
        for ( int j=0; j<u.n2; j++ )
        {
            if ( boundary.interior_test(u.x1(i),u.x2(j)))
            {
                double sum = 0;
                int n = 0;
                
                for ( int k=-n_smooth; k<=n_smooth; k++ )
                for ( int l=-n_smooth; l<=n_smooth; l++ )
                {
                    if (boundary.interior_test( u.x1(i+k) , u.x2(j+l) ))
                    {
                        sum += u(i+k,j+l);
                        n++;
                    }
                }
                
                if (n>0)
                    u_smooth(i,j) = sum/n;
                else
                    u_smooth(i,j) = u(i,j);
            }
            else
                u_smooth(i,j) = u(i,j);
        }
        
        u = u_smooth;
        
        
        
        
        
        grid2D<> Du(u);
        
        laplacian_2d_hdaf Del2;
        Del2.init( N,N, u.b1-u.a1, u.b2-u.a2, 12, 12, 0.5, 0.5 );
        
        Del2.execute( u.array, Du.array );
        
        {
            double sum = 0.0;
            int count = 0;
            
            for (int i=0; i<N; i++)
            for (int j=0; j<N; j++)
            {
                double x = u.x1(i);
                double y = u.x2(j);
                
                if ( -1<=x and x<=1 and -1<=y and y<= 1 )
                {
                    sum += (Du(i,j));
                    count ++;
                }
            }
            
            cout << "result: " << sum/count << endl;
        }
        
        
        
        
        
        // --------------------
        
        //plotGrid2D_1( Du, "/workspace/output/SDE/project1/Du.png", color_map_red );
        
        
        plot2D plt( -10.0, 10.0, -10.0, 10.0, N, N);
        plt = ml_white;
        
        for (int k=0; k<n_bd; k++)
            plt.ptLineThick( boundary.x_point(k),boundary.y_point(k), boundary.x_point(k+1), boundary.y_point(k+1), 5, color_map_red(f[k]) );
        
        for (int i=0; i<N; i++)
        for (int j=0; j<N; j++)
        {
            if ( u(i,j) != 0.0 )
                plt.set_px( i,j, color_map_blue( u(i,j) ) );
        }
        
        char str[20];
        sprintf(str, "n_runs:%03d", n_runs );
        
        plt.text_1(10, 10, str, ml_white, ml_black, 20,0 ); 
        
        sprintf(fname, "/workspace/output/SDE/project1/n_runs:%03d", n_runs );
        
        plt.png( fname );
        
    }
    
    
    /*
    //cout << "failure rate: " << 100.0*(double)fails/(eval_count*n_runs) << "%\n";
        
    cout << "eval_count: " << eval_count << endl;
    cout << "failed_runs: " << failed_runs << endl;	
    cout << "fail rate: " << (100.0*failed_runs)/(n_runs*eval_count) << "%\n";
        */  

    
    
    std_exit();
}


