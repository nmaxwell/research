#include <mathlib/math/std_math.h>




class Polygon2D
{
public:
    float * x_points;
    float * y_points;
    int n_points;
    
    Polygon2D():x_points(0),y_points(0),n_points(0) {}
    ~Polygon2D();
    
    void deallocate();
    void init(int n_points);
    void debug(bool display_data=false);
    
public:
    float & x_point(const int & k);
    float & y_point(const int & k);
    
public:
    void set_inscribed_circle(int n, float radius, float x_center=0, float y_center=0, float theta_0=0 );

public:
    int interior_test(const float & x, const float & y);
};



int Polygon2D::interior_test(const float & x, const float & y)
{
    // from http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
    
    int i, j, c = 0;
    for (i = 0, j = n_points-1; i < n_points; j = i++)
        if ( ((y_points[i]>y) != (y_points[j]>y)) && (x < (x_points[j]-x_points[i]) * (y-y_points[i]) / (y_points[j]-y_points[i]) + x_points[i]) )
			c = !c;
    return c;
}


#ifdef PLOT2D_H

void plot(Polygon2D & pg, plot2D & plt, ml_color const & color )
{
    for (int k=0; k<pg.n_points-1; k++)
    {
        plt.ptLine( pg.x_point(k),pg.y_point(k), pg.x_point(k+1), pg.y_point(k+1), color );
    }
}

#endif



void Polygon2D::debug(bool display_data)
{
    std::cout << "Polygon2D instance: \n";
    std::cout << "\taddress:\t" << this << endl;
    std::cout << "\tn_points:\t" << n_points << endl;
    std::cout << "\tx_points:\t" << x_points << endl;
    std::cout << "\ty_points:\t" << y_points << endl;
    
    if (display_data and x_points and y_points)
    {
        std::cout << "data:\n";
        for (int k=0; k<n_points; k++)
            std::cout << "\t" << k << "\t" << this->x_point(k) << "\t" << this->y_point(k) << endl;
    }
    
    std::cout << "end of " << this << " debug.\n";
}

void Polygon2D::set_inscribed_circle(int n, float radius, float x_center, float y_center, float theta_0 )
{
    // n sided regular polygon
    // n_poibnts = n+1
    
    int n_points = n+1;
    this->init(n_points);
    
    for (int k=0; k<n_points-1; k++)
    {
        double theta1 = (ml_2pi*k)/(n_points-1)+theta_0;
        double theta2 = (ml_2pi*(k+1))/(n_points-1)+theta_0;
        
        this->x_point(k) = cos(theta1)*radius+x_center;
        this->x_point(k+1) = cos(theta2)*radius+x_center;
        this->y_point(k) = sin(theta1)*radius+y_center;
        this->y_point(k+1) = sin(theta2)*radius+y_center;
    }
}


float & Polygon2D::x_point(const int & k)
{
    if (k>=n_points or x_points==NULL)
        std::cerr << "error in Polygon2D::x_point\n";
    
    return x_points[k];
}

float & Polygon2D::y_point(const int & k)
{
    if (k>=n_points or y_points==NULL)
        std::cerr << "error in Polygon2D::y_point\n";
    
    return y_points[k];
}



void Polygon2D::init(int n )
{
    this->deallocate();
    
    x_points = ml_alloc<float> (n);
    y_points = ml_alloc<float> (n);
    n_points = n;
}






void Polygon2D::deallocate()
{
    if (x_points) ml_free(x_points);
    if (y_points) ml_free(y_points);
    x_points = 0;
    y_points = 0;
}

Polygon2D::~Polygon2D()
{
    if (x_points) ml_free(x_points);
    if (y_points) ml_free(y_points);
    x_points = 0;
    y_points = 0;
}
















