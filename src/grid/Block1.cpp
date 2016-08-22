

#include <grid/Block1.hpp>

//  Block1():Grid() {}; //in include
Block1::Block1(Vec3 n1, Vec3 n2, int_4 nx):Block1::Block1() {
  Geo1* tmp = new Geo1Line(n1, n2); 
  add(tmp, nx);
  delete(tmp); 
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
  bool loop = false; 
  for (auto i = i0; i < nx+1; ++i) {
    // basic closure catch;
    auto pt = a->f(a->s0 + (double)(i) * del);
    if (i > 0 && (pt - *listVertex[nVertex]).abs() < 0.9*eps) {
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

Block1::Block1(initializer_list<initializer_list<double> > pt, double del):Block1::Block1() {
  auto itend = pt.end()-1; bool loop = false; 
  if (Vec3(*pt.begin()) == Vec3(*(pt.end()-1))) {loop = true; --itend;} 
  for (auto it = pt.begin(); it != itend; ++it) {
    Geo1* tmp = new Geo1Line(Vec3(*it), Vec3(*(it+1))); 
    add(tmp, 1);
    delete(tmp);
  }
  if (loop) addCell({(int_8)listVertex.size()-1, 0}); 
  resolve(del);
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


