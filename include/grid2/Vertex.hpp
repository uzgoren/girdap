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
#ifndef VERTEX
#define VERTEX

#include <grid2/Scheme.hpp>

// Classes to be define later, but needed
class Grid; 
class Cell;

// Vertex
class Vertex : public Vec3 {
protected:
    
public:
  vector<int_8 > cell; 
  // Cell contains structured array
  // many constains node shared by many cells; unstructured
  vector<int_8 > many;
  // --- NOT IMPLEMENTED -- //
  // When both cell and many has elements it is a mixed node;
  // Mixed nodes should be dealt specifically
  // --- NOT IMPLEMENTED -- //
  
  // Coefficients of a bilinear equation - used for mapping from x to xhat
  //    i.e. for quad: 
  //    x = xcoef(0) + xcoef(1)*xhat0 + xcoef(2)*xhat1 + xcoef(3)*xhat0*xhat1
  bool coefUpdate; 
  VecX<double> xcoef, ycoef, zcoef; 

  // Weights are computed using trasformed coordinate xhat; 
  Vec3 xhat;  
  //    i.e. for quad:
  //    p_v = (1-xhat1)*(1-xhat0)*p_0 + xhat0*(1-xhat1)*p_1+..
  //          xhat0*xhat1*p_2 + (1-xhat0)*xhat1*p_3; 

  // For gradients we need to fit p into a different bilinear equation; 
  //     but we only need weigths, which are computed from xhat again;  
  // OR we can simply calculate using gauss formula! Not sure yet; 

  int_2 bndr; 

  int_8 id; 
  Grid* grid; 
  
  Vertex():Vec3() { coefUpdate = true;}
  Vertex(double const & a, double const & b, double const & c) :  Vec3 (a, b, c) {coefUpdate = true;}
  Vertex(initializer_list<double> a) : Vec3( a ) {coefUpdate = true; }
  Vertex(Vec3 a) : Vec3(a) {coefUpdate = true;}

  // 00 getCoord 
  Vec3 getCoord() { return Vec3(*this); }

  // 01 RESET node's cell list
  void reset() {cell.clear(); coefUpdate = true;}
  // 01a RESET only one cell; 
  void reset(int_2 s) {cell.assign(s, -1);coefUpdate = true;}
  // 01b RESET and replace connectivity; 
  void reset(initializer_list<int_8> l) {cell.assign(l.begin(), l.end());coefUpdate = true;}

  // 02 Change number of cells (structure of a node)
  void cellResize(int_2 size); 
  // 03a assigns cell at a prescribed loc; Quad 0-1-2-3 (CCW) Hexa 0-1..7-8 (CCW-CCW)
  //     given that it doesn't have a meaningful value (<0)
  void setCell(int_2 ind, int_8 value);
  // 03b replaces whole structure based on list given; 
  void setCell(initializer_list<int_8> l) {
    cellResize(l.size()); 
    for (auto i = 0; i < l.size(); ++i) 
      {if (cell[i] < 0) cell[i] = *(l.begin()+i); }
  }
  
  // 03c replace ind with another cell -- Not a duplicate of setcell 
  //                                     (it overrides existing cell)
  void replaceCell(int_2 ind, int_8 value);

  // 04 gets a pointer to the cell object at ind
  shared_ptr<Cell > *getCell(int_2 ind, bool debug=0); 

  // 05 gets a pointer to an immediate neighboring vertex (IMPORTANT); 
  //     d stands for direction; and 1 represents x1, 2 for x2 and 3 for x3; 
  //     negative values of 1, 2, 3 reverses the direction.
  //     NOTE: ngbr should not be used during grid generation/adaptation
  //           it is only made available after faces are created!
  shared_ptr<Vertex > *ngbr(int_2 d);
  shared_ptr<Cell> getFace(int_2 order);

  // 06 interpolation scheme of a Field variable 
  Scheme<double> phi(vector<shared_ptr<Boundary> > const &bc, int_2 bias=0);  
  Scheme<double> phi(shared_ptr<Var> var, int_2 bias=0);  
  
  double evalPhi(shared_ptr<Var> &var, Vec3 *x=NULL);
  double evalPhi(VecX<double> &phi, vector<shared_ptr<Boundary> > const &bc, Vec3* x=NULL); 
 
  bool setInterpCoef(); 
  //vector<double> getIntWeight();
  vector<double> getIntWeight(Vec3* x = NULL); 

  //  double getIntVal(shared_ptr<Var> var); 
  // double getIntVal(shared_ptr<Var> var, Vec3 x); 

  

};

#endif
