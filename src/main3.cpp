
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
  Block2* grid = new Block2({0, 0, 0}, {1, 1, 0}, 50, 50);

  grid->adaptCriteria(); 

  grid->addVar({"dx", "dy"}); 
  auto dx = grid->getVar("dx"); 
  auto dy = grid->getVar("dy"); 

  for (auto v: grid->listVertex) {
    for (auto i =0; i< 2; ++i) {
      if (v->face[i]) {
	auto x = v->face[i]->getCoord().x(); 
	auto y = v->face[i]->getCoord().y(); 
	auto theta = atan2(y-0.5, x-0.5); 
	auto r = sqrt(2.0); 
	if (v->face[i]->prev < 0) {
	  dx->listBC.push_back(shared_ptr<Boundary>(new Boundary()));
	  (*(dx->listBC.rbegin()))->setBC(0, (r*cos(theta)+0.5));
	  v->face[i]->prev = -(dx->listBC.size()); 
	  
	  dy->listBC.push_back(shared_ptr<Boundary>(new Boundary()));
	  (*(dy->listBC.rbegin()))->setBC(0, (r*sin(theta)+0.5));
	  v->face[i]->prev = -(dy->listBC.size());
	  
	} else if (v->face[i]->next < 0) { 
	  dx->listBC.push_back(shared_ptr<Boundary>(new Boundary()));
	  (*(dx->listBC.rbegin()))->setBC(0, (r*cos(theta)+0.5));
	  v->face[i]->next = -(dx->listBC.size()); 
	  
	  dy->listBC.push_back(shared_ptr<Boundary>(new Boundary()));
	  (*(dy->listBC.rbegin()))->setBC(0, (r*sin(theta)+0.5));
	  v->face[i]->next = -(dy->listBC.size());
	}
      }
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
	    //	     - grid->source(0, 10000)
	    ); 
  grid->unlockBC(); 

  grid->lockBC(dy); 
  dy->solve(
  	    grid->laplace(1) 
  	    //	     - grid->source(0, 10000)
  	    ); 
  grid->unlockBC(); 

  for (auto v: grid->listVertex) {
    double sumx = 0, sumy = 0; double icnt = 0; double jcnt = 0; 
    for (auto c : v->cell) {
      if (c >=0 ) {
	sumx += dx->data[c]; ++icnt; 
	sumy += dy->data[c]; ++jcnt;
      }
    }
    v->data[0] = sumx/icnt; 
    v->data[1] = sumy/icnt; 
  }
  //grid->listVar.clear(); 

  grid->makeFace(); 
  grid->addVar({"T"}); 
  auto T = grid->getVar("T"); 
  T->solver = "CG"; 
  T->itmax =1000; 
  T->tol = 1e-6; 

  T->setBC("east", "val", 100);
  T->setBC("west", "val", 160); 
  T->setBC("north", "grad", 250*20, -250); 
 
  grid->lockBC(T); 
  T->solve(grid->laplace(2.0) + grid->source(0,1000)); 
  grid->unlockBC(); 
  
  // grid->lockBC(dx); 
  // dx->solve(
  // 	    grid->laplace(1) 
  // 	    //	     - grid->source(0, 10000)
  // 	    ); 
  // grid->unlockBC(); 

  ofstream myfile; 
  myfile.open("reg.vtk"); 
  myfile << grid << endl;
  myfile.close();    
 

  delete(grid); 

  return 0; 
}
