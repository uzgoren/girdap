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
#ifndef CELL
#define CELL

#include "Vertex.hpp"

class Var; 
//class Grid;
//class Vertex; 

class Cell {
protected:
  short int type;   
public:
  int_8 id; 
  Grid* grid; 
  bool isAlive;   
  vector<int_8> node;
  vector<int_8> face; 
  int_2 orient; // x-0; y-1; z-2; no particular dir-3;  
  int_8 next;
  int_8 prev;
  vector<int_2> level; 
  vector<int_2> adapt; 
  vector<bool> masterx, mastery, masterz; 
  bool cmx, cmy, cmz; 
    
  Cell(initializer_list<int_8> l, initializer_list<int_2> a={0,0,0}) { 
    node.assign(l.begin(), l.end()); 
    level.assign(a.begin(), a.end());
    adapt.assign(3, 0);     
    next = -1; prev = -1; 
    isAlive = true; 
  }

  short int getType() { return type;}
  void setType(short int a) { type = a;} 

  shared_ptr<Vertex> *getVertex(int_2 i); 
  void setVertex(int_2 i, int_8 v) { node[i] = v; }
    
  void reset() {node.clear();}
  void reset(initializer_list<int_8> l) {node.assign(l.begin(), l.end()); assignCelltoNode();} 

  void virtual assignCelltoNode() {}
  bool virtual refine() {return false;}
  void virtual refine(int dir) {return;}
  bool virtual coarsen(int dir) {return false;}  
  void virtual convertToSimpleBlock(initializer_list<int> n, bool debug=false) {};
  int_8 virtual hangingVertexOnFace(int_2 i, int_2 j=0) { return -1;}

  //  Vec3 virtual grad(shared_ptr<Var > phi) {return Vec3(); }
  void virtual getBCPhi(shared_ptr<Var> phi, double &a, double &b) {};  
  //Scheme virtual phi(int_2 bias = 0) {};

  double phiVal_iface(VecX<double> &phi, int_2 const &bias); 
  double phiVal_bcface(VecX<double> &phi, vector<shared_ptr<Boundary> > const &bc, int_2 const &bias); 
  double phiVal(VecX<double> &phi, vector<shared_ptr<Boundary> > const &bc, int_2 bias=0); 
  double phiVal(VecX<double> &phi, int_2 bias=0); 
  double phiVal(shared_ptr<Var> &var, int_2 bias=0); 

  Scheme<double> phi(vector<shared_ptr<Boundary> > const &bc, int_2 bias=0); 
  Scheme<Vec3> grad(vector<shared_ptr<Boundary> > const &bc); 
  Scheme<double> normGrad(vector<shared_ptr<Boundary> > const &bc); 

  Scheme<double> phi(shared_ptr<Var> var, int_2 bias=0); 
  Scheme<Vec3> grad(shared_ptr<Var> var);
  Scheme<double> normGrad(shared_ptr<Var> var);

  // Scheme virtual getGrad() {};
  //  void virtual makeSimpleBlock(initializer_list<int_8> n) { }

  // void virtual split(initializer_list<int> n, bool debug=false) { }
  // void virtual splitAndExpand(initializer_list<int> n, bool debug=false) { }
  // void virtual extrude(int dir, vector<shared_ptr<Vertex > > jnd) { } 
  Vec3 virtual edge(int_2 i=0) { return Vec3(0);  } 
  // int_2 virtual ngbrLevel(int_2 dir=0) { return 0; }
  // vector<shared_ptr<Vertex> > virtual edgeVertexList(int_2 i) {
  //   vector<shared_ptr<Vertex> > l; return l;
  // }
  // shared_ptr<Cell > virtual ngbrCell(int_2 i, int_2 j) { 
  //   shared_ptr<Cell > c; return c;
  // }
  vector<int_8> virtual ngbrCellList() { return vector<int_8>(); }
  void virtual checkIslandLevels(int dir) {};
  void virtual checkNgbrLevel(int dir) {};

  Vec3 virtual vol() { return Vec3(1,0,0); }
  double virtual dx() { return 0; }
  double virtual dy() { return 0; }
  double virtual dz() { return 0; }
  bool virtual isPointIn(Vec3 a) { return false; }

  Vec3 getCoord() { 
    Vec3 sum(0,0,0); double icnt = 0; 
    for (auto i = 0; i<node.size(); ++i) {
      if (auto v = getVertex(i)) {
	sum += (**v); icnt++; 
      }
    }
    return sum/icnt; 
  }

  friend ostream &operator<<(ostream &out, shared_ptr<Cell > const &a){
    if (a) {
      out << (a)->node.size() << " "; 
      for (auto e = (a)->node.begin(); e != (a)->node.end(); ++e) {
	if (*e >= 0) out << *e << " "; 
      }
      out << endl; 
      return out; 
    } else {
      return out << "NULL" << endl;
    }
  };
  virtual ~Cell() {}; 
};

#include "Line.hpp"

class Tri: public Cell {
public:
  Tri(initializer_list<int_8> l):Cell(l) {
    setType(5); 
  }
};

#include "Quad.hpp"

class Hexa: public Cell {
protected: 
  const int_2 map[8] = {6, 7, 4, 5, 2, 3, 0, 1};
public:
  //  short int type;   
  Hexa(initializer_list<int_8> l):Cell(l) {
    setType(12); 
  }
  

};




#endif
