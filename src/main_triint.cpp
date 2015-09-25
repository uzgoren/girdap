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
//#include <field/Var>
//#include <grid2/Grid>

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
  
  Vec3 pthat(1, 1, 1);
  Vec3 pt; 
  pt[0] = interp.linFunc(pthat, xcoef); 
  pt[1] = interp.linFunc(pthat, ycoef); 
  pt[2] = interp.linFunc(pthat, zcoef); 
  
  auto xhat = interp.findXhat(pt, xcoef, ycoef, zcoef); 

  cout << pt << endl; 
  cout << xhat << endl; 

  return 0; 
}
