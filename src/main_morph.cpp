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
  auto t = clock(); 
  //Block2* grid = new Block2({0, 0, 0}, {0.5, 0.3, 0}, 50, 30);
  Block2* grid = new Block2({0, 0, 0}, {1, 1, 0}, 10, 10);

  //grid->adaptCriteria(); 
  grid->addVar({"T"}); 
  auto T = grid->getVar("T");

  T->set(0); 
  double pi = 4.0*atan(1); 
  
  for (auto j = 0; j < 5; ++j) {
    for (auto i = 0; i < grid->listCell.size(); ++i) {
      auto x = grid->listCell[i]->getCoord(); // - Vec3(0.5, 0.5); 
      auto r = (grid->listCell[i]->getCoord() - Vec3(0.5, 0.5)).abs(); 

      T->set(i, 1.0/(1.0 + exp(-2.0*80*(0.15-r)))); 
    }
    //    auto gt = grid->valGrad(T); 
    //grid->solBasedAdapt(gt); 
    grid->solBasedAdapt2(grid->getError(T));
    grid->adapt(); 
  }


  grid->addVar({"dx", "dy"}); 
  auto dx = grid->getVar("dx"); 
  auto dy = grid->getVar("dy"); 

  for (auto c: grid->listCell) {    
    if ((c->getCoord() - Vec3(0.5, 0.5)).abs() < 0.14) {
      for (auto i = 0; i < c->node.size() ; ++i) 
	grid->listVertex[c->node[i]]->cell[(i+2)%4] = -1; 
 
      c->isAlive = false; 
    }
  }

  grid->cleanGrid(); 
  grid->checkGrid(); 
  grid->setCurrentLevels();
  grid->makeFace(); 

  for (auto f: grid->listFace) {
    auto fr = f->getCoord(); 
    auto x = fr.x(); 
    auto y = fr.y(); 
    auto theta = atan2(y-0.5, x-0.5); 
    auto r = 0.15;

    if (f->prev < 0) {
	
      dx->listBC.push_back(shared_ptr<Boundary>(new Boundary()));
      if (abs((fr-Vec3(0.5,0.5)).abs()-r) < 0.1) {
	(*(dx->listBC.rbegin()))->setBC(0, (r*cos(theta)+0.5));
      } else {
	(*(dx->listBC.rbegin()))->setBC(0, x); //(r*cos(theta)+0.5));
      }
      f->prev = -(dx->listBC.size()); 
	  
      dy->listBC.push_back(shared_ptr<Boundary>(new Boundary()));
      if (abs((fr-Vec3(0.5,0.5)).abs()-r) < 0.1) {
	(*(dy->listBC.rbegin()))->setBC(0, (r*sin(theta)+0.5));
      } else {
	(*(dy->listBC.rbegin()))->setBC(0, y); //(r*sin(theta)+0.5));
      }
      f->prev = -(dy->listBC.size());
	  
    } else if (f->next < 0) { 
      dx->listBC.push_back(shared_ptr<Boundary>(new Boundary()));
      if (abs((fr-Vec3(0.5,0.5)).abs()-r) < 0.1) {
	(*(dx->listBC.rbegin()))->setBC(0, (r*cos(theta)+0.5));
      } else {      
	(*(dx->listBC.rbegin()))->setBC(0, x); //(r*cos(theta)+0.5));
      }
      f->next = -(dx->listBC.size()); 
      
      dy->listBC.push_back(shared_ptr<Boundary>(new Boundary()));
      if (abs((fr-Vec3(0.5,0.5)).abs()-r) < 0.1) {
	(*(dy->listBC.rbegin()))->setBC(0, (r*sin(theta)+0.5));
      } else {
	(*(dy->listBC.rbegin()))->setBC(0, y); //(r*sin(theta)+0.5));
      }
      f->next = -(dy->listBC.size());
    }

  }

  dx->set(0.0); 
  // dx->setBC("north", "val", 0); 
  // dx->setBC("south", "val", 0); 
  dy->set(0.0); 
  // dy->setBC("east", "val", 0); 
  // dy->setBC("west", "val", 0); 
  dx->solver = "CG"; 
  dy->solver = "CG"; 
  dx->itmax = 1000; 
  dx->tol = 1e-6; 
  dy->itmax = 1000; 
  dy->tol = 1e-6; 


  grid->lockBC(dx); 
  dx->solve(
	    grid->laplace(1) 
	    ); 
  grid->unlockBC(); 

  grid->lockBC(dy); 
  dy->solve(
  	    grid->laplace(1) 
  	    ); 
  grid->unlockBC(); 

  // update vertex 
  auto icount = 0; 
  for (auto v: grid->listVertex) {
    double sumx = 0, sumy = 0; double icnt = 0; double jcnt = 0; 
    for (auto c : v->cell) {
      if (c >=0 ) {
	sumx += dx->data[c]; ++icnt; 
	sumy += dy->data[c]; ++jcnt;
      }
    }
    v->data[0] = sumx/icnt; 
    v->data[1] = sumy/jcnt; 
    auto thet = atan2(v->data[1]-0.5, v->data[0]-0.5); 
    for (auto ic = 0; ic < 4; ++ic) {
      if (v->cell[ic] < 0) {
    	if (abs((v->getCoord()-Vec3(0.5, 0.5)).abs()-0.15)< 0.2) {
    	  v->data[0] = 0.5 + 0.15*cos(thet); 
    	  v->data[1] = 0.5 + 0.15*sin(thet); 
    	  break; 
    	}
      }
    }
  }
  
  // update hanging nodes; 

  for (auto v: grid->listVertex) { 
    for (auto iv = 0; iv < 4; ++iv) { 
      if (v->cell[iv] >= 0 && v->cell[iv] == v->cell[(iv+1)%4]) {
	auto vb = grid->listCell[v->cell[iv]]->getVertex((iv+3)%4); 
	auto vn = grid->listCell[v->cell[iv]]->getVertex((iv+2)%4); 
	if (vn && vb) {
	  v->set(0.5*(**vn +  **vb)); 
	}
      }
    }
  }

  grid->listVar.clear(); 

  grid->makeFace(); 


  grid->addVar({"T"}); 
  T = grid->getVar("T"); 
  T->solver = "CG"; 
  T->itmax =1000; 
  T->tol = 1e-6; 

  T->setBC("east", "val", 100);
  T->setBC("west", "val", 160); 
  T->setBC("north", "grad", 250*20, -250); 
 
  for (auto itt = 0; itt< 3; ++itt) {

    grid->lockBC(T); 
    T->solve(grid->laplace(2.0) + grid->source(0,1000)); 
    grid->unlockBC(); 
    
  // grid->lockBC(dx); 
  // dx->solve(
  // 	    grid->laplace(1) 
  // 	    //	     - grid->source(0, 10000)
  // 	    ); 
  // grid->unlockBC(); 

    grid->writeVTK("reg"); 
    
    grid->solBasedAdapt2(grid->getError(T));
    grid->adapt(); 

  }
  grid->writeVTK("reg"); 

  // ofstream myfile; 
  // myfile.open("reg.vtk"); 
  // myfile << grid << endl;
  // myfile.close();    
 

  delete(grid); 

  return 0; 
}
