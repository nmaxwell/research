#include "bfio.hpp"
#include "serialize.hpp"

using std::istringstream;
using std::ifstream;
using std::ofstream;
using std::set;
using std::queue;
using std::cerr;

extern "C"
{
  void vdsqrt_(int* n, double*, double*);
  void vdsincos_(int* n, double*, double*, double*);
}

//---------------------------------------
int serialize(const Entry& e, ostream& os, const vector<int>& mask)
{
  iC( serialize(e._grid, os, mask) );
  iC( serialize(e._mats, os, mask) );
  iC( serialize(e._dir, os, mask) );
  return 0;
}

int deserialize(Entry& e, istream& is, const vector<int>& mask)
{
  iC( deserialize(e._grid, is, mask) );
  iC( deserialize(e._mats, is, mask) );
  iC( deserialize(e._dir, is, mask) );
  return 0;
}

//---------------------------------------
int BFIO::setup(map<string,string>& opts)
{
  map<string,string>::iterator mi;
  vector<int> all(1,1);
  //mi = opts.find("-" + prefix() + "N");
  //if(mi!=opts.end()) {    istringstream ss((*mi).second);    ss>>_N;  }
  mi = opts.find("-" + prefix() + "EPS");
  if(mi!=opts.end()) {    istringstream ss((*mi).second);    ss>>_EPS;  }
  mi = opts.find("-" + prefix() + "fi");
  if(mi!=opts.end()) {    istringstream ss((*mi).second);    ss>>_fi;  }
  mi = opts.find("-" + prefix() + "datfile");
  char datfile[100];
  if(mi!=opts.end()) {    istringstream ss((*mi).second);    ss>>datfile;  }
  ifstream fin(datfile);
  iC( deserialize(_e2dmap, fin, all) );
  cerr<<_EPS<<" "<<_e2dmap.size()<<endl;
  return 0;
}







int BFIO::kernel(int N, const vector<Point2>& xi, const vector<Point2>& ki, CpxNumMat& result)
{
	vector<Point2> xs(xi);
	vector<Point2> ks(ki);
	
	if (_fi != 2)
		iA(0);
	
	int m = xs.size();
	int n = ks.size();
	for(int i=0; i<m; i++)
		xs[i] = xs[i]/double(N); //LEXING: IMPORTANT nick: why?
	for(int j=0; j<n; j++)
		ks[j] = ks[j]/double(N); //LEXING: IMPORTANT nick: why?
	
	DblNumMat phase(m,n);
	
	for(int j=0; j<n; j++)
	{
		double w = 1.0 + ks[j](0)*3.0;
		double s = ks[j](1);
			
		for(int i=0; i<m; i++)
		{
			double x1 = xs[i](0);
			double x2 = xs[i](1);
			
			phase(i,j) = -2.0*N*( w*sqrt( (s-x1)*(s-x1) + x2*x2 + 1 ));
		}
	}
	
	DblNumMat imag_part(m,n), real_part(m,n);
	result.resize(m,n);
	
	for(int j=0; j<n; j++)
	for(int i=0; i<m; i++)
	{
		sincos( phase(i,j), &(imag_part(i,j)), &(real_part(i,j)));
		result(i,j) = cpx( real_part(i,j), imag_part(i,j) );
	}
	
	return 0;
}





//---------------------------------------
int BFIO::check(const CpxNumMat& f, const CpxNumMat& u, int NC, double& relerr)
{
  int N = f.m();
  vector<Point2> src;
  for(int j=0; j<N; j++)
    for(int i=0; i<N; i++)
      src.push_back( Point2(i, j) );      //src.push_back( Point2(i-N/2, j-N/2) );
  vector<cpx> app(NC);
  vector<cpx> dir(NC);

  for(int g=0; g<NC; g++) {
    int x1 = int( floor(drand48()*N) );
    int x2 = int( floor(drand48()*N) );
    //
    app[g] = u(x1,x2);
    //
    vector<Point2> trg;    trg.push_back( Point2(x1,x2) );
    CpxNumMat res(1,N*N); iC( kernel(N, trg, src, res) );
    CpxNumMat resaux(N,N,false,res.data());
    cpx ttl(0,0);
    for(int j=0; j<N; j++)
      for(int i=0; i<N; i++)
	ttl = ttl + resaux(i,j) * f(i,j);
    dir[g] = ttl;
  }
  vector<cpx> err(NC);
  for(int g=0; g<NC; g++)
    err[g] = app[g] - dir[g];
  double dn = 0;
  double en = 0;
  for(int g=0; g<NC; g++) {
    dn += abs(dir[g])*abs(dir[g]);
    en += abs(err[g])*abs(err[g]);
  }
  dn = sqrt(dn);
  en = sqrt(en);
  relerr = en/dn;
  //
  return 0;
}
