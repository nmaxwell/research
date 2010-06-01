#ifndef _BFIO_HPP_
#define _BFIO_HPP_

#include "comobject.hpp"
#include "numtns.hpp"
#include "offtns.hpp"
#include "vec2t.hpp"
#include "vecmatop.hpp"

using std::vector;
using std::pair;
using std::map;

class Entry
{
public:
  DblNumVec _grid;
  NumVec<CpxNumMat> _mats;
  CpxNumMat _dir;
public:
  Entry() {;}
  ~Entry() {;}
  DblNumVec& grid() { return _grid; }
  NumVec<CpxNumMat>& mats() { return _mats; }
  CpxNumMat& dir() { return _dir; }
};

int serialize(const Entry&, ostream&, const vector<int>&);
int deserialize(Entry&, istream&, const vector<int>&);

//-----------------------------
class BFIO: public ComObject
{
public:
  int _EPS;
  int _fi;
  map<int, Entry> _e2dmap;
public:
  BFIO(const string& p): ComObject(p) {;}
  ~BFIO() {;}
  int& EPS() { return _EPS; }
  int& fi() { return _fi; }
  map<int, Entry>& e2dmap() { return _e2dmap; }
  //
  int setup(map<string,string>& opts);
  int eval(const CpxNumMat& f, CpxNumMat& u);
  int kernel(int N, const vector<Point2>& trg, const vector<Point2>& src, CpxNumMat& res);
  int check(const CpxNumMat& f, const CpxNumMat& u, int NC, double& relerr);
};

#endif

