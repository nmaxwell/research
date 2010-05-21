





typedef grid3D<double, double, double > grid;

struct grid_vec
{
    grid _1;
    grid _2;
    grid _3;
    
    void copy( grid_vec & rhs )
    {
        _1.copy( rhs._1 );
        _2.copy( rhs._2 );
        _3.copy( rhs._3 );
    }
    
    void copy_data( grid_vec & rhs )
    {
        _1.copyData( rhs._1 );
        _2.copyData( rhs._2 );
        _3.copyData( rhs._3 );
    }
    
    void operator= ( grid_vec & rhs )
    {
        this->copy_data( rhs );
    }
    
    void operator *= ( double rhs )
    {
        _1 *= rhs;
        _2 *= rhs;
        _3 *= rhs;
    }
    
    void operator = ( double rhs )
    {
        _1 = rhs;
        _2 = rhs;
        _3 = rhs;
    }
    
    void operator *= ( grid & rhs )
    {
        _1 *= rhs;
        _2 *= rhs;
        _3 *= rhs;
    }
    
};

struct grid_mat
{
    grid _11;
    grid _12;
    grid _13;
    grid _21;
    grid _22;
    grid _23;
    grid _31;
    grid _32;
    grid _33;
    
};

struct grid_symmat
{
    grid _11;
    grid _22;
    grid _33;
    grid _23;
    grid _13;
    grid _12;
};

class strain_tensor
{
public:
    grad_3d_hdaf grad;
    
public:
    
    void init( grid & g, int m1, int m2, int m3, double gam1, double gam2, double gam3 )
    {
        grad.init( g.n1, g.n2, g.n3, g.b1-g.a1, g.b2-g.a2, g.b3-g.a3, m1, m2, m3, gam1, gam2, gam3 );
    }
    
    void execute( grid_vec & u, grid_symmat & del_u )
    {
        del_u._11.copy(u._1);
        del_u._22.copy(u._1);
        del_u._33.copy(u._1);
        del_u._23.copy(u._1);
        del_u._13.copy(u._1);
        del_u._12.copy(u._1);
        
        grad.linear_strain_tensor( u._1.array, u._2.array, u._3.array, del_u._11.array, del_u._12.array, del_u._13.array, del_u._22.array, del_u._23.array, del_u._33.array );
    }
};


class stress_tensor
{
public:
    strain_tensor del;
    grid_symmat eps;
    grid stiffness_tensor[21];
    
public:
    
    void init( grid & g, grid * C, int m1, int m2, int m3, double gam1, double gam2, double gam3 )
    {
        del.init( g, m1, m2, m3, gam1, gam2, gam3 );
        
        for (int k=0; k<21; k++ )
            stiffness_tensor[k] = C[k];
    }
    
    void execute( grid_vec & u, grid_symmat & sig )
    {
        // sig = C : eps
        
        int N = u._1.n1*u._1.n2;
        
        sig._11.copy(u._1);
        sig._12.copy(u._1);
        sig._13.copy(u._1);
        sig._22.copy(u._1);
        sig._23.copy(u._1);
        sig._33.copy(u._1);
        
        del.execute( u, eps );
        
        for ( int n=0; n<N; n++ )
            sig._11[n] =
                stiffness_tensor[0][n]*eps._11[n] +
                stiffness_tensor[1][n]*eps._22[n] +
                stiffness_tensor[2][n]*eps._33[n] +
                stiffness_tensor[3][n]*eps._23[n]*2 +
                stiffness_tensor[4][n]*eps._13[n]*2 +
                stiffness_tensor[5][n]*eps._12[n]*2;
        
        for ( int n=0; n<N; n++ )
            sig._22[n] =
                stiffness_tensor[1][n]*eps._11[n] +
                stiffness_tensor[6][n]*eps._22[n] +
                stiffness_tensor[7][n]*eps._33[n] +
                stiffness_tensor[8][n]*eps._23[n]*2 +
                stiffness_tensor[9][n]*eps._13[n]*2 +
                stiffness_tensor[10][n]*eps._12[n]*2;
        
        for ( int n=0; n<N; n++ )
            sig._33[n] =
                stiffness_tensor[2][n]*eps._11[n] +
                stiffness_tensor[7][n]*eps._22[n] +
                stiffness_tensor[11][n]*eps._33[n] +
                stiffness_tensor[12][n]*eps._23[n]*2 +
                stiffness_tensor[13][n]*eps._13[n]*2 +
                stiffness_tensor[14][n]*eps._12[n]*2;
        
        for ( int n=0; n<N; n++ )
            sig._23[n] =
                stiffness_tensor[3][n]*eps._11[n] +
                stiffness_tensor[8][n]*eps._22[n] +
                stiffness_tensor[12][n]*eps._33[n] +
                stiffness_tensor[15][n]*eps._23[n]*2 +
                stiffness_tensor[16][n]*eps._13[n]*2 +
                stiffness_tensor[17][n]*eps._12[n]*2;
        
        for ( int n=0; n<N; n++ )
            sig._13[n] =
                stiffness_tensor[4][n]*eps._11[n] +
                stiffness_tensor[9][n]*eps._22[n] +
                stiffness_tensor[13][n]*eps._33[n] +
                stiffness_tensor[16][n]*eps._23[n]*2 +
                stiffness_tensor[18][n]*eps._13[n]*2 +
                stiffness_tensor[19][n]*eps._12[n]*2;
        
        for ( int n=0; n<N; n++ )
            sig._12[n] =
                stiffness_tensor[5][n]*eps._11[n] +
                stiffness_tensor[10][n]*eps._22[n] +
                stiffness_tensor[14][n]*eps._33[n] +
                stiffness_tensor[17][n]*eps._23[n]*2 +
                stiffness_tensor[19][n]*eps._13[n]*2 +
                stiffness_tensor[20][n]*eps._12[n]*2;
        
    }
};


class anis_diff_op
{
public:
    stress_tensor del_sig;
    grad_3d_hdaf del;
    grid_symmat sig;
    
public:
    
    void init( grid & g, grid * C, int m1, int m2, int m3, double gam1, double gam2, double gam3 )
    {
        del_sig.init( g, C, m1, m2, m3, gam1, gam2, gam3  );
        del.init( g.n1, g.n2, g.n3, g.b1-g.a1, g.b2-g.a2, g.b3-g.a3, m1, m2, m3, gam1, gam2, gam3 );
        
        sig._11.copy(g);
        sig._12.copy(g);
        sig._13.copy(g);
        sig._22.copy(g);
        sig._23.copy(g);
        sig._33.copy(g);
    }
    
    void execute( grid_vec & u, grid_vec & del_u )
    {
        del_sig.execute( u, sig );
        del.grad_dot_sym( sig._11.array, sig._12.array, sig._13.array, sig._22.array, sig._23.array, sig._33.array, del_u._1.array, del_u._2.array, del_u._3.array );
    }
};


class anisotropic_propagator
{
public:
    grid_vec * U;
    grid_vec * V;
    int exp_order;
    anis_diff_op Del;
    
    anisotropic_propagator( ):exp_order(0),U(0),V(0),Del() {}
    
    ~anisotropic_propagator() {
        if (U) delete [] U;
        if (V) delete [] V;
        U=0; V=0; exp_order=0; }
    
public:
    
    void init( grid3D<> & g, int new_order, grid * C, int m1, int m2, int m3, double gam1, double gam2, double gam3 )
    {
        exp_order = new_order;
        
        if (U) delete [] U;
        if (V) delete [] V;
        U = new grid_vec [exp_order+2];
        V = new grid_vec [exp_order+1];
        
        for (int k=0; k<exp_order+2; k++ )
        {
            U[k]._1 = g;
            U[k]._2 = g;
            U[k]._3 = g;
        }
        
        for (int k=0; k<exp_order+1; k++ )
        {
            V[k]._1 = g;
            V[k]._2 = g;
            V[k]._3 = g;
        }
        
        Del.init( g, C, m1, m2, m3, gam1, gam2, gam3  );
    }
    
    void operator() ( double t, grid_vec & u0, grid_vec & v0, grid_vec & u1, grid_vec & v1 )
    {
        propagate( t, u0, v0, u1, v1 );
    }
    
    void propagate ( double t, grid_vec & u0, grid_vec & v0, grid_vec & u1, grid_vec & v1 )
    {
        u1.copy(u0);
        v1.copy(v0);
        
        U[0] = u0;
        V[0] = v0;
        
        for (int n = 1; n<=exp_order+1; n++)
        {
            Del.execute( U[n-1], U[n] );
            
            U[n] *= t*t/(4*n*n-2*n);
        }
        
        for (int n = 1; n<=exp_order; n++)
        {
            Del.execute( V[n-1], V[n] );
            
            V[n] *= t*t/(4*n*n-2*n);
        }
        
        
        u1 = 0.0;
        v1 = 0.0;
        
        for (int n = 0; n<=exp_order; n++)
        GRID2D_LINLOOP( u0._1 )
        {
            u1._1[i] += U[n]._1[i];
            u1._1[i] += V[n]._1[i]* t/(2.0*n+1.0);
            
            v1._1[i] += U[n+1]._1[i]*((4.0*n+6.0)*n+2.0)/(t*(2.0*n+1.0));
            v1._1[i] += V[n]._1[i];
            
            u1._2[i] += U[n]._2[i];
            u1._2[i] += V[n]._2[i]* t/(2.0*n+1.0);
            
            v1._2[i] += U[n+1]._2[i]*((4.0*n+6.0)*n+2.0)/(t*(2.0*n+1.0));
            v1._2[i] += V[n]._2[i];
        }
    }
};




class driver
{
public:
    grid_vec * F;
    int exp_order;
    anis_diff_op Del;
    
    driver( ):exp_order(0),F(0),Del() {}
    
    ~driver() {
        if (F) delete [] F;
        F=0; exp_order=0; }
    
public:
    
    void init( grid3D<> & g, int new_order, grid * C, int m1, int m2, int m3, double gam1, double gam2, double gam3 )
    {
        exp_order = new_order;
        
        if (F) delete [] F;
        F = new grid_vec [exp_order+2];
        
        for (int k=0; k<exp_order+2; k++ )
        {
            F[k]._1 = g;
            F[k]._2 = g;
            F[k]._3 = g;
        }
        
        Del.init( g, C, m1, m2, m3, gam1, gam2, gam3  );
    }
    
    void operator() ( double dt, grid_vec & f1, grid_vec & f2, grid_vec & u1, grid_vec & v1 )
    {
        drive( dt, f1, f2, u1, v1 );
    }
    
    void drive( double dt, grid_vec & f1, grid_vec & f2, grid_vec & u1, grid_vec & v1 )
    {
	    F[0] = f1;
        
        for (int n = 1; n<=exp_order+1; n++)
        {
            Del.execute( F[n-1], F[n] );
            
            F[n] *= dt*dt/(4*n*n-2*n);
        }
        
        for (int n = 0; n<=exp_order; n++)
        GRID3D_LINLOOP( u1._1 )
        {
            u1._1[i] += (dt/2)*F[n]._1[i];
            u1._2[i] += (dt/2)*F[n]._2[i];
            
            v1._1[i] += F[n+1]._1[i]*(double)((4*n+6)*n+2)/(2*(2*n+1));            
            v1._2[i] += F[n+1]._2[i]*(double)((4*n+6)*n+2)/(2*(2*n+1));
        }
        
        GRID3D_LINLOOP( u1._1 )
        {
	        u1._1[i] += (dt/2)*f2._1[i];
            u1._2[i] += (dt/2)*f2._2[i];
        }
        
    }
};




