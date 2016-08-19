#ifndef BLOCK1

#include "Grid.hpp"
#include <functional>

class Geo1 {
public:
  double s0, s1;
  int_8 nx;
  std::function<Vec3 (double)> f;
  Geo1() {
    s0=0.0; s1=1.0;
    f = [&] (double t) -> Vec3 {return Vec3(0.0) + t*Vec3(1.0, 1.0);};
  }
  virtual ~Geo1() {};
};

class geoLine: public Geo1 {
public:
  geoLine():Geo1() {};
  geoLine(Vec3 x0, Vec3 x1) {
    f = [=] (double t) -> Vec3 {return x0 + t*(x1-x0);}; 
  };
};

class geoSine: public Geo1 {
public:
  geoSine():Geo1() {};
  geoSine(Vec3 x0, Vec3 x1, double a, double freq) {
    Vec3 norm(x0[1]-x1[1], x1[0]-x0[0]);
    norm = norm/norm.abs();
    Vec3 del = x1-x0; double pi = 4*atan(1.0); 
    f = [=] (double t) -> Vec3 {return x0 + t*del + a*sin(freq*t*pi)*norm;}; 
  };
};

class geoCircle: public Geo1 {
public:
  geoCircle():Geo1() {};
  geoCircle(Vec3 x0, double r0) {
    double pi = 4*atan(1.0); 
    s0 = 0; s1 = 2*pi;
    f = [=] (double t) -> Vec3 {return x0 + r0*Vec3(cos(t), sin(t));}; 
  };
};


class Block1: public Grid {
public:
  Block1():Grid() {};
  Block1(Vec3 n1, Vec3 n2, int_4 nx);
  Block1(Geo1* a, int_4 nx); 
  Block1(initializer_list<double> n1, initializer_list<double> n2, int_4 nx):Block1((Vec3)n1, Vec3(n2), nx){}
  
  Block1(string geo, initializer_list<initializer_list<double> > pt, int_4 nx); 
  Block1(double s0, double s1, int_8 nx
	 , std::function<Vec3 (double)> f); 
  
  void resolve(double del); 

  void add(Geo1* a, int_4 nx); 
  void add(Block1& o);
  void add( double s0, double s1, int_8 nx
	    , std::function<Vec3 (double)> f);
  
  void subtract(Block1& o);

  void attach( double s0, double s1, int_8 nx
	       , std::function<Vec3 (double)> f);
  void attach(Block1& o); 
  
};




#endif
