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

class Geo1Line: public Geo1 {
public:
  Geo1Line():Geo1() {};
  Geo1Line(Vec3 x0, Vec3 x1) {
    f = [=] (double t) -> Vec3 {return x0 + t*(x1-x0);}; 
  };
};

class Geo1Sine: public Geo1 {
public:
  Geo1Sine():Geo1() {};
  Geo1Sine(Vec3 x0, Vec3 x1, double a, double freq) {
    Vec3 norm(x0[1]-x1[1], x1[0]-x0[0]);
    norm = norm/norm.abs();
    Vec3 del = x1-x0; double pi = 4*atan(1.0); 
    f = [=] (double t) -> Vec3 {return x0 + t*del + a*sin(freq*t*pi)*norm;}; 
  };
};

class Geo1Circle: public Geo1 {
public:
  Geo1Circle():Geo1() {};
  Geo1Circle(Vec3 x0, double r0, double a0=0.0, double a1=360.0) {
    double pi = 4*atan(1.0);
    if (a0 == a1) a1 = a0 + 360;
    if (a1 > a0+360) a1 = a0+360;
    if (a1 < a0-360) a1 = a0-360; 
    s0 = a0/180.0*pi; s1 = a1/180.0*pi;
    f = [=] (double t) -> Vec3 {return x0 + r0*Vec3(cos(t), sin(t));}; 
  };
};

// CONTENTS in src/grid/Block1.cpp
class Block1: public Grid {
public:
  Block1():Grid() {};
  Block1(Vec3 n1, Vec3 n2, int_4 nx);
  Block1(Geo1* a, int_4 nx); 
  Block1(initializer_list<double> n1, initializer_list<double> n2, int_4 nx):Block1(Vec3(n1), Vec3(n2), nx){}
  
  Block1(initializer_list<initializer_list<double> > pt, double del); 
  Block1(double s0, double s1, int_8 nx
	 , std::function<Vec3 (double)> f); 
  
  void resolve(double del); 

  void add(Geo1* a, int_4 nx); 
  void add(Block1& o);
  void add( double s0, double s1, int_8 nx
	    , std::function<Vec3 (double)> f);

  // mergers
  void addVol(Block1& o);
  void subVol(Block1& o);

  // modify self
  void trans(Vec3 dir, double dist);
  void rot(Vec3 x0, Vec3 norm, double angle); 

  // Create 1D higher
  Block2 rotate(Vec3 x0, Vec3 norm, int_8 nx);
  Block2 extrude(Vec3 norm, double dist, int_8 nx);  
  
};




#endif
