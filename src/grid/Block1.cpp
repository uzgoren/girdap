

#include <grid/Block1.hpp>

Block1::Block1(Vec3 n1, Vec3 n2, int_4 nx):Block1::Block1() {
  add(make_shared<Geo1>( Geo1Line(n1, n2) ), nx);
  return;
}

Block1::Block1(shared_ptr<Geo1> a, int_4 nx) {
  add(a, nx); 
}

void Block1::add(shared_ptr<Geo1> a, int_4 nx) {
  add(a->s0, a->s1, nx, a->f); 
}

Block1::Block1(initializer_list<initializer_list<double> > pt):Block1::Block1() {
  auto itend = pt.end()-1; bool loop = false; 
  if (Vec3(*pt.begin()) == Vec3(*(pt.end()-1))) {loop = true; --itend;} 
  for (auto it = pt.begin(); it != itend; ++it) {
    add(make_shared<Geo1>( Geo1Line(Vec3(*it), Vec3(*(it+1))) ), 1);
  }
  if (loop) addCell({(int_8)listVertex.size()-1, 0}); 
  //  resolve(del);
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
  double eps = (f(s0 + del) - f(s0)).abs();
  int_8 nVertex=listVertex.size(), i0 = 0; 
  if (nVertex > 0 && (f(s0) - **listVertex.rbegin()).abs() < eps*0.1) {
    i0 = 1; --nVertex;
  }
  bool loop = false; int_8 i1 = nx+1; 
  if ((f(s0) - f(s1)).abs() < eps*0.5) {loop = true; --i1;}
  for (auto i = i0; i < i1; ++i) {
    auto pt = f(s0 + (double)(i) * del);
    if (i > 0 && (pt - **listVertex.rbegin()).abs() < eps*0.1) continue; 
    addVertex(pt); 
  }
  if (loop) {
    addCell({(*listVertex.rbegin())->id, nVertex});    
  }
  for (auto i = nVertex; i < listVertex.size()-1; ++i) {
    addCell({i, i+1});
  }
}


