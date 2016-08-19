

#include <grid/Block1.hpp>

//  Block1():Grid() {}; //in include
Block1::Block1(Vec3 n1, Vec3 n2, int_4 nx):Block1::Block1() {   
  Vec3 del = (n2 - n1)/nx;
  for (auto i = 0; i < nx+1; ++i) {
    addVertex(n1[0] + i * del);
  }
  for (auto i = 0; i < nx; ++i) {
    addCell({i, i+1}); 
  }
  setCurrentLevels(); 
  //makeFace(); 
  cout << "Block1: Cells: " << listCell.size() << endl; 
  return;
}

Block1::Block1(Geo1* a, int_4 nx) {
  add(a, nx); 
}

void Block1::add(Geo1* a, int_4 nx) {
  double del = (a->s1 - a->s0)/(double)nx;
  double eps = (a->f(a->s0 + del) - a->f(a->s0)).abs();
  int_8 nVertex=listVertex.size(), i0 = 0; 
  if (nVertex > 0 && (a->f(a->s0) - **listVertex.rbegin()).abs() < eps*0.1) {i0 = 1; --nVertex;}
  cout << "First point: " << a->f(a->s0)<<endl; 
  bool loop = false; 
  for (auto i = i0; i < nx+1; ++i) {
    // basic closure catch;
    auto pt = a->f(a->s0 + (double)(i) * del);
    if (i > 0 && (pt - *listVertex[nVertex]).abs() < 0.9*eps) {
      cout << " loop true " << endl; 
      loop = true;
      break; 
    }
    addVertex(pt); 
  }
  if (loop) {
    addCell({(*listVertex.rbegin())->id, nVertex});    
  }
  for (auto i = nVertex; i < listVertex.size()-1; ++i) {
    addCell({i, i+1});
  }
}

//   Block1(initializer_list<double> n1, initializer_list<double> n2, int_4 nx):Block1((Vec3)n1, Vec3(n2), nx){} // in include

Block1::Block1(string geo, initializer_list<initializer_list<double> > pt, int_4 nx):Block1::Block1() {
  if (geo.compare("poly") == 0) {
    double del = 0; bool loop = false; 
    for (auto it = pt.begin(); it != pt.end()-1; ++it) {
      del += ((Vec3)(*(it+1)) - (Vec3)(*it)).abs();
      addVertex(*it);
    }
    if (del == 0) return;
    del = del/(double)nx; 
    if ((Vec3)*(pt.end()-1) == (Vec3)*(pt.begin())) {loop = true; }
    else {addVertex(*(pt.end()-1));}
    for (auto i = 0; i < listVertex.size()-1; ++i) {
      addCell({i, i+1});
    }
    if (loop) {
      addCell({(int_8)(listVertex.size()-1), 0});
    }      
    setCurrentLevels();
    resolve(del);
    //makeFace(); 
  } else if (geo.compare("arc") == 0) {
    auto pi = 4.0*atan(1.0);
    if (pt.size() < 2) {return;}
    auto x = (Vec3)(*pt.begin()); 
    double r = ((Vec3)(*(pt.begin()+1))).abs(); 
    Vec3 d = (pt.size() >= 3) ? (Vec3)(*(pt.begin()+2)) : Vec3(0, 360);
    if (d[1]-d[0] > 360 || d[1] == d[0]) d[1] = d[0] + 360;
    if (d[0]-d[1] > 360) d[1] = d[0] - 360; 
    bool loop = (abs(d[1] - d[0]) == 360) ? true : false;
    d = pi/180.0*d;
    double del = (d[1] - d[0])/nx;
    for (auto i = 0; i < nx; ++i) {
      addVertex(x + Vec3({r*cos(d[0] + (double)i*del), r*sin(d[0] + (double)i*del)}));
    }
    if (!loop) addVertex(x + Vec3(r*cos(d[0] + (double)(nx)*del),
				  r*sin(d[0] + (double)(nx)*del)));
    for (auto i = 0; i < nx-1; ++i) {
      addCell({i, i+1});	
    }
    if (loop) addCell({(int_8)(listVertex.size()-1), 0});
    else addCell({(int_8)listVertex.size()-2, (int_8)listVertex.size()-1}); 
    setCurrentLevels();
  }    
}

Block1::Block1(double s0, double s1, int_8 nx
	       , std::function<Vec3 (double)> f) {
  add(s0, s1, nx, f); 
}

void Block1::resolve(double del) {
  bool isadapt=true; 
  while (isadapt) {
    isadapt = false; int nold = listCell.size(); 
    for (auto c: listCell) {
      if (c->level[0] == levelHighBound[0]) continue; 
      if (c->vol().abs() > del) {
	c->adapt[0] = 1;
	isadapt = true;
      }
    }
    if (isadapt) adapt();
    if (nold == listCell.size()) break; 
  }
  setCurrentLevels(); 
}

void Block1::add(Block1& o) {
  // SIMPLE add;
  auto noldv = listVertex.size(); 
  for (auto v : o.listVertex) {
    addVertex(*v);
  }
  for (auto c : o.listCell) {
    addCell({c->node[0] + (int_8)noldv, c->node[1] + (int_8)noldv});
  }
}
// Parametric add; 
void Block1::add(double s0, double s1, int_8 nx
		 , std::function<Vec3 (double)> f) {
  double del = (s1 - s0)/(double)nx;
  int ncell = listCell.size(); 
  for (auto i = 0; i < nx+1; ++i) {
    addVertex(f(s0 + (double)i * del));
  }
  if (listCell.size() > ncell) addCell({ncell-1, ncell}); 
  for (auto i = 0; i < nx; ++i) {
    addCell({ncell+i, ncell+i+1});
  }
}

void Block1::attach(double s0, double s1, int_8 nx
		    , std::function<Vec3 (double)> f) {
  int_8 nVertex = listVertex.size(); 
  Block1 tmp(s0, s1, nx, f);
  add(tmp);
  if (tmp.listVertex.size() > 0) addCell({nVertex-1, nVertex});
  tmp.~Block1(); 
}
