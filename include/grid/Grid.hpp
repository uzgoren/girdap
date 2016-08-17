/*
 ***************************************************************************
 *   Copyright (C) 2015, Eray Uzgoren                                      *
 *                                                                         *
 *   This program is free software: you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation, either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>. *
 ***************************************************************************
*/
#ifndef GRID
#define GRID

#define SL 0
#define FOU 1

#include <memory>
#include <typeinfo>
#include <cmath>
#include <time.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>

#include "../base/Interp.hpp"
#include "../field/Var.hpp"
#include "Vertex.hpp"
#include "Cell.hpp"

class LinSys;
class triLinSys;  

//class Vertex; 
//class Var; 

class Grid {
private:
  unsigned int filecnt; 

public:
  vector<shared_ptr<Vertex > > listVertex; 
  vector<shared_ptr<Cell > > listCell; 
  vector<shared_ptr<Cell > > listFace; 
  vector<shared_ptr<Var > > listVar;
  vector<int_8 > otherVertex; 

  vector<bool> selectCell; 

  //vector<shared_ptr<Boundary > > listBNDR; 
  shared_ptr<Var > thisVar; 
  double dt, dt0, cfl; 

  int_8 nCellSize; 
  int nFace; 

  Interp int3D; // interpolation scheme

  int_2 levelLowBound[3], levelHighBound[3]; 
  int_2 levelMin[3], levelMax[3]; 

  double minD, maxD, meanD; 
  // vector<shared_ptr<BNDR> > 

  Grid() { nCellSize = 0; for (auto i = 0; i<3; ++i) {levelLowBound[i] = 0; levelHighBound[i] = 4; cfl = 0.5; }; filecnt=0;}

  Grid(initializer_list<initializer_list<double > > pts):Grid() {
    addVertex(pts);
  }
  
  Grid(Vec3 pts): Grid() {
    addVertex(pts); 
  }

  Grid(initializer_list<initializer_list<double> > pts, 
       initializer_list<initializer_list<int_8> > cell): Grid(pts) {
    addCell(cell);
    // correct connectivity ! missing (hanging nodes) for the next step; 
    // makeFace(); 
  }

  // ~Grid() {
  //   // for (auto i = 0; i < otherVertex.size(); ++i) otherVertex[i].reset(); 
  //   // for (auto i = 0; i < listFace.size(); ++i) listFace[i].reset(); 
  //   // for (auto i = 0; i < listCell.size(); ++i) listCell[i].reset(); 
  //   // for (auto i = 0; i < listVertex.size(); ++i) listVertex[i].reset(); 
    
  //   if (listVertex.size() > 0 ) cout << listVertex[0].use_count() << endl; 

  //   otherVertex.clear(); cout << "OtherVertex size = "<< otherVertex.size() << endl; 
  //   listFace.clear(); cout << "listFace size = "<< listFace.size() << endl; 
  //   listCell.clear(); cout << "listCell size = "<< listCell.size() << endl; 
  //   listVertex.clear(); cout << "listVertex size = "<< listVertex.size() << endl; 
  // }


  // Add vertex to the list; 
  void addVertex(Vec3 v) {
    listVertex.emplace_back( shared_ptr<Vertex>(new Vertex(v)) ); 
    (*listVertex.rbegin())->id = listVertex.size()-1;
    (*listVertex.rbegin())->grid = this; 
    otherVertex.emplace_back( -1 ); 
  }     

  void addVertex(initializer_list<double > pt) {
     addVertex(Vec3(pt)); 
  }

  void addVertex(vector<Vec3 > pts) {
    for (auto pt : pts) { addVertex(pt); }
  }

  void addVertex(initializer_list<initializer_list<double > > pts) {
     for (auto i = pts.begin(); i != pts.end(); ++i) { addVertex(*i); }    
  }

  // // //Add cell to the list
  void addCell(initializer_list<int_8 > l) {    
    if (l.size() == 2) {
      listCell.emplace_back( shared_ptr<Cell > (new Line(l)) ); 
 //     //listCell.emplace_back(shared_ptr<Cell > (new Line(tmp)) );
    } else if (l.size() == 3) {
  //     //listCell.emplace_back(shared_ptr<Cell > (new Tri(tmp)) ); 
    } else if (l.size() == 4) {
      listCell.emplace_back( shared_ptr<Cell > (new Quad(l)) ); 
    } else if (l.size() == 8) {
       //listCell.emplace_back( shared_ptr<Cell > (new Hexa(tmp)) ); 
    } else {
      cout << "Cell type not understood!"<<endl; 
      exit(1); 
    }
    nCellSize += l.size() + 1;
    auto c = *(listCell.rbegin()); 

    c->id = listCell.size()-1;
    c->grid = this; 
    c->assignCelltoNode();
    for (auto v: listVar) {
      if (v->loc == 0) {
	v->data.push_back(0); 
	v->data.uncompress();
      }
    }
    c->masterx.resize(levelHighBound[0]+1, false); 
    c->mastery.resize(levelHighBound[1]+1, false); 
    c->masterz.resize(levelHighBound[2]+1, false); 
  }

  void addCell(initializer_list<initializer_list<int_8 > > link) {
    for (auto it = link.begin(); it != link.end(); ++it) {
      addCell(*it);
    }
  }

  void setDt(double cdt) {
    auto U = getVar("u"); auto V = getVar("v"); auto W = getVar("w");     
    auto maxU = abs(min(U->data.min(), min(V->data.min(), W->data.min()))); 
    maxU = max(maxU, max(U->data.max(), max(V->data.max(), W->data.max())));
    double dx = 1e10;  
    for (auto c : listCell) {
      double cdx, cdy, cdz; 
      cdx = c->dx(); cdy = c->dy(); cdz = c->dz(); 
      if (cdx > 0) dx = min(dx, cdx); 
      if (cdy > 0) dx = min(dx, cdy); 
      if (cdz > 0) dx = min(dx, cdz); 
    }
    cout << " SetDt: " << dt << " Vel: " << maxU << endl; 
    dt0 = dt;
    if (maxU < 1e-6 || dx < 1e-6)
      dt = cdt;
    else 
      dt = min(cdt, cfl*dx/maxU);
  }

  void setIntCoef() {
    int_8 sum = 0; 
    for (auto v : listVertex) {
      //if (v->coefUpdate) { 
	v->setInterpCoef(); 
	sum++;
	//}
    }
  }

  void makeFace() {
    setIntCoef(); 

    listFace.clear(); 
    for (auto c : listCell)  
      c->face.clear(); 
    
    for (auto v : listVertex) {
      if (v->cell.size() == 2) {
	listFace.emplace_back(new Cell({v->id})); 
	(*listFace.rbegin())->next = (v->cell[1] < 0) ? -1 : v->cell[1]; 
	(*listFace.rbegin())->prev = (v->cell[0] < 0) ? -2 : v->cell[0];
	(*listFace.rbegin())->grid = this;
	(*listFace.rbegin())->orient = 0; 
	(*listFace.rbegin())->id = listFace.size()-1; 
	if (v->cell[1] >= 0) listCell[v->cell[1]]->face.push_back((*listFace.rbegin())->id) ;
	if (v->cell[0] >= 0) listCell[v->cell[0]]->face.push_back((*listFace.rbegin())->id); 
	return; 
      } else if (v->cell.size() == 4) { 
	auto n = v->ngbr(1); 
	if (n && *n) {
	  listFace.emplace_back(shared_ptr<Cell>(new Line({v->id, (*n)->id}))); 
	  (*listFace.rbegin())->next = (v->cell[2] < 0) ? -2 : v->cell[2];
	  (*listFace.rbegin())->prev = (v->cell[1] < 0) ? -1 : v->cell[1];
	  if (v->cell[2] < 0) v->cell[2] = -2; 
	  if (v->cell[1] < 0) v->cell[1] = -1;
	  (*listFace.rbegin())->grid = this;
	  (*listFace.rbegin())->orient = 0; 
	  (*listFace.rbegin())->id = listFace.size()-1; 
	  if (v->cell[1] >= 0) listCell[v->cell[1]]->face.push_back((*listFace.rbegin())->id); 
	  if (v->cell[2] >= 0) listCell[v->cell[2]]->face.push_back((*listFace.rbegin())->id); 
	} 
	if (v->ngbr(-1)) { 
	  if (v->cell[0] < 0) v->cell[0] = -1; 
	  if (v->cell[3] < 0) v->cell[3] = -2; 
	}
	n = v->ngbr(2); 
	if (n && *n) { 
	  listFace.emplace_back(shared_ptr<Cell>(new Line({(*n)->id, v->id}))); 
	  (*listFace.rbegin())->next = (v->cell[2] < 0) ? -4 : v->cell[2];
	  (*listFace.rbegin())->prev = (v->cell[3] < 0) ? -3 : v->cell[3];
	  if (v->cell[2] < 0) v->cell[2] = -4; 
	  if (v->cell[3] < 0) v->cell[3] = -3; 
	  (*listFace.rbegin())->grid = this;
	  (*listFace.rbegin())->orient = 1; 
	  (*listFace.rbegin())->id = listFace.size()-1; 
	  if (v->cell[2] >= 0) listCell[v->cell[2]]->face.push_back((*listFace.rbegin())->id); 
	  if (v->cell[3] >= 0) listCell[v->cell[3]]->face.push_back((*listFace.rbegin())->id); 
	} 
	if (v->ngbr(-2)) { 
	  if (v->cell[0] < 0) v->cell[0] = -3; 
	  if (v->cell[1] < 0) v->cell[1] = -4; 
	}
      }	
    }
    nFace = listFace.size();

    for (auto v : listVar) {
      if (v->loc == 2) v->data.resize(nFace); 
    }
  }

  void addVar(std::string n, int t = 0); 
  void addVec(std::string n, int t = 0); 
  void addVar(initializer_list<std::string> nl, int t = 0);
  void addVec(initializer_list<std::string> nl, int t = 0);
  shared_ptr<Var> getVar(std::string n);   
  void lockBC(shared_ptr<Var> v);
  void unlockBC(); 


  VecX<double> getCoord(int dir) {
    VecX<double> val(listCell.size()); val.uncompress(); 
    for (auto i = 0; i < listCell.size(); ++i) {
      double sum=0; double icnt = 0; 
      for (auto j = 0; j < listCell[i]->node.size(); ++j) {
	if (listCell[i]->node[j] < 0) continue; 
	auto k = listCell[i]->getVertex(j); 
	sum = sum + (**k)[dir]; icnt++; 
      }
      val[i] = sum/icnt; 
      //cout << sum  << " " << icnt << " " << val[i] << endl; 
    }
    return val; 
  }

  VecX<Vec3> getVel() { 
    VecX<Vec3> val(listCell.size()); val.uncompress(); 
    auto u = *(listVar.begin()); auto v = *(listVar.begin()+1); auto w = *(listVar.begin()+2);
    for (auto i = 0; i< listCell.size(); ++i) {
      val[i].set({u->get(i), v->get(i), w->get(i)});
    }
    return val; 
  }
  
  VecX<Vec3> getVecVertex(std::string a) {
    VecX<Vec3> val(listVertex.size()); 
    shared_ptr<Var > u, v, w; 
    if (a == "u") {
      u = getVar("u"); v = getVar("v"); w = getVar("w");
    } else {
      u = getVar(a+"x"); v = getVar(a+"y"); w = getVar(a+"z"); 
    }
    if (!u || !v || !w) {
      cout << "Vector not found! " << a << endl; 
      exit(-1); 
    }

    for (auto i = 0; i<listVertex.size(); ++i) {
      val[i][0] = (u->loc==1) ? u->get(i): listVertex[i]->evalPhi(u); //.eval(u); 
      val[i][1] = (v->loc==1) ? v->get(i): listVertex[i]->evalPhi(v); //.eval(v); 
      val[i][2] = (w->loc==1) ? w->get(i): listVertex[i]->evalPhi(w);//.eval(w); 
    }
    return val; 
  }

  VecX<double> getPhiVertex(shared_ptr<Var> phi) {
    VecX<double> val(listVertex.size()); 
    //    for (auto i = 0; i<listVertex.size(); ++i) {
      // double sum=0; double icnt = 0;  
      // for (auto j = 0; j < listVertex[i]->cell.size(); ++j) {
      // 	if (listVertex[i]->cell[j] < 0) continue; 
      // 	auto k = listVertex[i]->cell[j]; 
      // 	sum = sum + phi->data[k]; 
      // 	++icnt; 
      // }
      //val[i] = listVertex[i]->phi(phi).eval(phi); //sum/icnt; 
    for (auto v : listVertex) {
      auto w = v->getIntWeight(); 
      for (auto i = 0; i < v->cell.size() ; ++i) 
	val[v->id] += w[i]*(phi->data[v->cell[i]]); 
    }
    return val; 
  }
  
#include "Adapt.hpp"

  void checkGeo() { 
    auto l0 = 0.0, l1 = 0.0; 
    for (auto i = cbegin(); i !=cend(); ++i)
      l0 += (*i)->vol().abs(); 
    for (auto i = vbegin(); i != vend(); ++i) { 
      auto j = (*i)->ngbr(1); 
      if (j != nullptr) {l1 += ((**i)^(**j))*Vec3(1,1,1); }
    }
    l1 *= 0.5; 
    cout << " L0 = " << l0 << " L1 = "<< l1 << " ";  
  }

  void checkGrid() {
    for (auto i = cbegin(); i != cend(); ++i) {
      if (((*i)->id < 0) || ((*i)->id) >= listCell.size()) {
	cout << "Cell "<< (*i)->id << " is not valid " << endl; 
	cout << "Cell list size = " << listCell.size() << endl; 
	exit(1);
      }
      for (auto j = 0; j < (*i)->node.size(); ++j) {
	if (((*i)->node[j]) < 0 || ((*i)->node[j]) > listVertex.size()) {
	  cout << "Cell "<< (*i)->id << " has problems in its vertex list @" << j ; 
	  cout << " is " << (*i)->node[j] << endl; 
	  cout << "Vertex list size = " << listVertex.size() << endl; 
	  exit(1);
	}
      }
    }
    if (otherVertex.size() != listVertex.size() && otherVertex.size() > 0) {
      cout << "Size of otherVertex does not match with listVertex size" << endl; 
      cout << "Vertex size = " << listVertex.size() << endl; 
      cout << "OtherVertex size = " << otherVertex.size() << endl;
      exit(1); 
    }
    for (auto i = vbegin(); i != vend(); ++i) {
      if (((*i)->id) < 0 || ((*i)->id) > listVertex.size()) {
	cout << "Vertex "<< (*i)->id << " is not valid " << endl; 
	  cout << "Vertex list size = " << listVertex.size() << endl; 
	exit(1);
      }
      for (auto j = 0; j < (*i)->cell.size(); ++j) {
	if (int_8((*i)->cell[j]) >= int_8(listCell.size())) {
	  cout << "Vertex "<<(*i)->id << " has problems in its cell list @" << j ; 
	  cout << " is " <<(*i)->cell[j] << endl;
	  cout << "Cell list size = " << listCell.size() << endl; 
	  exit(1);
	}
      }
    }    
  }
	
  void debug() {
    for (auto v = vbegin(); v!=vend(); ++v) {
      cout << (*v)->id << " ("; 
      for (auto i = 0; i < (*v)->cell.size(); ++i) {
	if ((*v)->cell[i] < 0) {cout << "x "; continue;}
  	cout << (*v)->cell[i] << " "; 
      }
      cout << ") "<<endl; 
    }
  }

  // INTERPOLATION at a point
  double interp(shared_ptr<Var> &phi, Vec3 x, int_8 i0=0) {
    //    cout << " Interp: " << x << endl; 
    if (i0 < 0 || i0 >= listVertex.size()) i0=0; 
    // cout << " IN search "<< endl; 
    i0 = searchVertexbyCoords(x, i0); // minimum distance; 
    //cout << " OUT search " << endl; 
    if (i0 < 0) { cout << " coefs not found! after search point! " << endl; exit(1);  return 0;} 
    return listVertex[i0]->evalPhi(phi, &x);    
  }

  Vec3 interpVec(shared_ptr<Var> &phi, Vec3 x, int_8 i0=0) { 
    //   cout << " in " << endl; 
    if (!phi->isVec) {cout << "Nonvector interpolation" << endl; return Vec3(0); }
    if (i0 < 0 || i0 >= listVertex.size()) i0=0; 
    i0 = searchVertexbyCoords(x, i0); 
    if (i0 < 0) return 0; //BC
    shared_ptr<Var > u, v, w; 
    if (phi->name == "u") {
      u = getVar("u"); v = getVar("v"); w = getVar("w");
    } else {
      auto a = phi->name; 
      a.pop_back(); 
      u = getVar(a+"x"); v = getVar(a+"y"); w = getVar(a+"z"); 
    } 
    auto val0 = listVertex[i0]->evalPhi(u, &x); 
    auto val1 = listVertex[i0]->evalPhi(v, &x); 
    auto val2 = listVertex[i0]->evalPhi(w, &x); 
    
    //    cout << " out "<< endl; 
    return Vec3(val0, val1, val2); 

  }
  /**

  VecX<double> expDdt();
  VecX<double> expDiv();
  VecX<double> expLaplace();
  VecX<double> expGrad();
  
  triLinSys impDdt(); 
  triLinSys impDiv();
  triLinSys impLaplace();
  triLinSys impGrad(); 
  **/
  void advanceDiv(shared_ptr<Var> &phi, VecX<Vec3> &vel, 
		  int &rko, double aveflux[3], int &tri, bool isflux[3]); 

  void laplace(LinSys &axb, shared_ptr<Cell > f, double const &c);
  LinSys laplace(double c, initializer_list<double> n={}); 
  LinSys laplace(VecX<double> c, initializer_list<double> n={}); 
  LinSys source(double c, double a, initializer_list<double> n={}); 
  LinSys source(double c, VecX<double> a, initializer_list<double> n={}); 

  void div(LinSys &axb, shared_ptr<Cell > f, Vec3 const &c); 
  LinSys div(VecX<Vec3> vel, double c, initializer_list<double> n={}); 
  LinSys divRK2E(VecX<Vec3> vel, double c); 
  LinSys divRK4E(VecX<Vec3> vel, double c); 
  LinSys div(VecX<Vec3> vel, VecX<double> c, initializer_list<double> n={}); 
  void getBCPhi(shared_ptr<Cell > f, shared_ptr<Var> phi, double &a, double &b);
  void getPhi_val(shared_ptr<Cell > f, shared_ptr<Var> phi, double &out); 
  VecX<double> valDiv(VecX<Vec3> vel);
  VecX<Vec3> valGrad(shared_ptr<Var > phi); 
  void correctVel(double coef);

  LinSys ddt(double c); 
  LinSys ddt(VecX<double> c); 
  void timeScheme(LinSys &axb, initializer_list<double> &n, VecX<double> &prev); 
  
  void triLaplace(triLinSys &axb, shared_ptr<Cell > f, double const &c);
  triLinSys laplace2(double c, initializer_list<double> n);
  triLinSys laplace2(VecX<double>& c, initializer_list<double> n);
  triLinSys source2(double c, double a, initializer_list<double> n);
  triLinSys source2(double c, VecX<double>& a, initializer_list<double> n);

  void triDiv(triLinSys &axb, shared_ptr<Cell > f, Vec3 const &c);
  triLinSys div2(VecX<Vec3>& vel, double c, initializer_list<double> n);
  triLinSys div2(VecX<Vec3>& vel, VecX<double> c, initializer_list<double> n);

#include "IO_grid.hpp"

#include "Connect_grid.hpp"

#include "Search_grid.hpp"  

  vector<shared_ptr<Vertex > >::iterator vbegin() { return listVertex.begin(); }
  vector<shared_ptr<Vertex > >::iterator vend() {return listVertex.end(); }
  
  vector<shared_ptr<Cell > >::iterator  cbegin() {return listCell.begin(); }
  vector<shared_ptr<Cell > >::iterator  cend() {return listCell.end(); }

};


class Block2: public Grid {
public:
  Block2():Grid() {};
  Block2(initializer_list<double> n1, initializer_list<double > n2, int_4 nx, int_4 ny):Grid() {  
    Vec3 node1 = n1; 
    Vec3 node2 = n2; 
    Vec3 del = (node2 - node1);  
    meanD = min(del[0]/double(nx), del[1]/double(ny)); 
    addVertex({
	 {node1[0], node1[1], node1[2]}
	,{node2[0], node1[1], node1[2]} 
	,{node2[0], node2[1], node1[2]} 
	,{node1[0], node2[1], node1[2]}
      }); 
    addCell({0,1,2,3});
    (*listCell.rbegin())->convertToSimpleBlock({nx,ny}); 
    setCurrentLevels(); 
    makeFace(); 
    //setQuadBoundary(); 
    cout << "Block2: Cells: " << listCell.size(); 
    cout << " Faces: " << nFace << endl; 
    addVec("u"); //, "v", "w"});
  }

  


}; 


class Block3: public Grid {
public:
  Block3():Grid() {};
  Block3(initializer_list<double> n1, initializer_list<double > n2, int_4 nx, int_4 ny, int_4 nz):Grid() {    
    Vec3 node1 = n1; 
    Vec3 node2 = n2; 
    addVertex({ 
	 {node1[0], node1[1], node1[2]}
	,{node2[0], node1[1], node1[2]}
	,{node2[0], node2[1], node1[2]}
	,{node1[0], node2[1], node1[2]}
	,{node1[0], node1[1], node2[2]}
	,{node2[0], node1[1], node2[2]}
	,{node2[0], node2[1], node2[2]}
	,{node1[0], node2[1], node2[2]}
      });
    //addCell({0,1,2,3,4,5,6,7}); 
    //if (listCell.capacity() < (listCell.size() + nx*ny*nz)) {listCell.reserve(listCell.size() + nx*ny*nz);}
    //(*listCell.rbegin())->split({nx,ny, nz}); 
    //    listCell.rbegin()->splitHexa(nx,ny,nz); 
  }

}; 


// template <typename T>
// inline bool equals(const std::weak_ptr<T>& t, const std::weak_ptr<T>& u)
// {
//     return !t.owner_before(u) && !u.owner_before(t);
// }

// template <typename T>
// inline bool equals(const std::weak_ptr<T>& t, const std::shared_ptr<T>& u)
// {
//     return !t.owner_before(u) && !u.owner_before(t);
// }


#endif
