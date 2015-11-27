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
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <time.h>
#include <iomanip>

#include <base/Interp.hpp>
#include <field/Var>
#include <grid2/Grid>

int main() {
  std::string flname; 
  int filecnt = 0; 
  ofstream myfile; 

  Interp interp; 

  VecX<double> x = {-1, -1, 13, -4, -1, -1, 13, -4};
  VecX<double> y = {-1, -1, 11,  8, -1, -1, 11,  8};
  VecX<double> z = { 0,  0,  0,  0,  0,  0,  0.5,  1.3};

  auto xcoef = interp.getCoef(x); 
  auto ycoef = interp.getCoef(y); 
  auto zcoef = interp.getCoef(z); 
  
  Vec3 pthat(0.7, 0.6, 1);
  Vec3 pt; 
  pt[0] = interp.linFunc(pthat, xcoef); 
  pt[1] = interp.linFunc(pthat, ycoef); 
  pt[2] = interp.linFunc(pthat, zcoef); 
  
  auto xhat = interp.findXhat(pt, xcoef, ycoef, zcoef); 

  Grid* grid = new Block2({0, 0, 0}, {1, 1, 0}, 50, 50); 

  // for (auto i = 0; i < grid->listVertex.size(); ++i) {
  //   auto w = grid->listVertex[i]->getIntWeight(); 
  //   double sum = 0; 
  //   cout << i << " : "; 
  //   for (auto v : w) { cout << v << ", "; sum += v; } 
  //   cout << "= " <<sum << endl; 
  // }

  grid->addVar({"T"}); 
  
  auto T = grid->getVar("T");
  T->setBC("west", "val", 0.8);
  T->setBC("south", "val", 1); 
  T->setBC("east", "val", -1); 
  T->setBC("north", "val", -0.5); 
  
  double pi = 4.0*atan(1); 
  for (auto i = 0; i < grid->listCell.size(); ++i) {
    auto c = grid->listCell[i]; 
    auto x = c->getCoord(); 
    T->set(i, sin(x[0]*4*pi)*cos(x[1]*2*pi)); // + cos(x[0]*4*pi)*sin(x[1]*2*pi)); 
  }

  cout << pt << endl; 
  cout << xhat << endl; 

     myfile.open("grid.vtk"); 
    myfile << grid << endl;
    myfile.close(); 

  return 0; 
}
