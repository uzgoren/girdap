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
#include <math.h>
#include <fstream>
#include <sstream>
#include <time.h>
#include <iomanip>

#include <field/Var>
#include <grid2/Grid>


int main() {
  auto k = 2.0; auto qdot = 5e3; auto h = 50; auto Tinf = 20;
  Block2* grid = new Block2({0, 0, 0}, {1, 1, 0}, 10, 10); 
  grid->levelHighBound = {2, 2}; 
  grid->addVar("T"); 

  Block2* uni = new Block2({0, 0, 0}, {1, 1, 0}, 40, 40); 
  uni->levelHighBound[0] = {0, 0}; 
  uni->addVar("T");   

  auto Tref = uni->getVar("T");
  Tref->solver = "BiCGSTAB";
  Tref->itmax = 1000; 
  
  Tref->set(100); 
  
  Tref->setBC("south", "grad", 0); 
  Tref->setBC("north", "grad", h/k*Tinf, -h/k);
  Tref->setBC("east", "val", 200); 
  Tref->setBC("west", "val", 100); 
  
  uni->lockBC(Tref); 					 // for boundary conditions
  Tref->solve( uni->laplace(k)  			 // thermal conductivity as gamma
	       + uni->source(0, qdot) ); 
  uni->unlockBC();
  uni->writeVTK("heat_ref"); 

  auto T = grid->getVar("T");
  auto u = grid->getVar("u"); 
  auto v = grid->getVar("v"); 
  T->solver = "BiCGSTAB";
  T->itmax = 1000; 

  T->set(100); 
  
  T->setBC("south", "grad", 0); 
  T->setBC("north", "grad", h/k*Tinf, -h/k);
  T->setBC("east", "val", 200); 
  T->setBC("west", "val", 100); 
  
  for (auto i = 0; i< 7; ++i) {
    grid->solBasedAdapt2(grid->getError2(T), 2e-3, 2e-1);
    grid->adapt();  
    
    grid->lockBC(T); 					 // for boundary conditions
    T->solve( grid->laplace(k)  			 // thermal conductivity as gamma
	      + grid->source(0, qdot) ); 
    grid->unlockBC();

    auto sum = 0; 
    for (auto v : grid->listVertex) {
      auto T0 = uni->interp(Tref, *v);
      auto T1 = v->evalPhi(T);
      sum += pow(T1-T0, 2); 
    }
    cout << " Error: " << sqrt(sum)/grid->listVertex.size() << endl; 
    auto err = grid->getError2(T); 
    u->set(err.comp(0)); 
    v->set(err.comp(1)); 
    
    grid->writeVTK("heat"); 
  }
    
  delete(grid); 
  
  return 0; 
}
