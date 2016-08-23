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

#include <linSolve/LinSys>
#include <grid2/Grid>
#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>
#include <time.h>
#include <field/Var>

int main() {
  double pi = 4.0*atan(1); 

  Block2* grid = new Block2({0, 0, 0}, {1, 1, 0}, 10, 10);
  
  grid->addVar({"T"}); 
  auto T = grid->getVar("T");
  T->set(0); 

  // Indicator for a circle; 
  for (auto j = 0; j < 5; ++j) {
    for (auto i = 0; i < grid->listCell.size(); ++i) {
      auto x = grid->listCell[i]->getCoord(); 
      auto r = (grid->listCell[i]->getCoord() - Vec3(0.5, 0.5)).abs(); 
      T->set(i, 1.0/(1.0 + exp(-2.0*80*(0.15-r)))); 
    }
    grid->solBasedAdapt2(grid->getError(T));
    grid->adapt(); 
  }

  // Create a surface grid at T = 0.5
  auto surf = contour(T, 0.5); 

  // Solves for Laplacian(x) = 0; and Laplacian(y) = 0; 
  // interface location and adapts; 
  grid->subtract(surf); 

  T->setBC("east", "val", 100);
  T->setBC("west", "val", 160); 
  T->setBC("north", "grad", 250*20, -250); 
 
  for (auto itt = 0; itt< 3; ++itt) {
    grid->lockBC(T); 
    T->solve(grid->laplace(2.0) + grid->source(0,1000)); 
    grid->unlockBC(); 

    grid->writeVTK("reg"); 
    
    grid->solBasedAdapt2(grid->getError(T));
    grid->adapt(); 
  }
  grid->writeVTK("reg"); 

  delete(grid); 

  return 0; 
}
