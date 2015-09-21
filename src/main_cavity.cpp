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

  //CAVITY FLOW 
  auto dt = 0.01; 
  auto t = clock(); 

  // GRID 
  Block2* grid = new Block2({0, 0, 0}, {1.0, 1.0, 0}, 10, 10);
  // grid->adaptCriteria(); 

  // CONST variables; 
  double rho = 1; double mu = 0.01; 

  // FIELD variables;
  grid->addVar({"p", "vor"}); 
    
  // initial and bc values; 
  auto u = grid->getVar("u"); 
  u->set(0.0);
  u->setBC("west", "val", 0); u->setBC("east", "val", 0);
  u->setBC("south", "val", 0); u->setBC("north", "val", 1); 
  auto v = grid->getVar("v"); 
  v->set(0.0);
  v->setBC("west", "val", 0); v->setBC("east", "val", 0);
  v->setBC("south", "val", 0); v->setBC("north", "val", 0); 

  auto p = grid->getVar("p");
  p->set(0.0); 

  auto vor = grid->getVar("vor"); 

  // solver behavior
  u->solver = "BiCGSTAB"; u->itmax = 20; u->tol = 1e-4; 
  v->solver = "BiCGSTAB"; v->itmax = 20; v->tol = 1e-4; 
  p->solver = "BiCGSTAB"; p->itmax = 2000; p->tol = 1e-6; 

  auto gp = grid->valGrad(p); 
  
 
  // Time control 
  grid->setDt(dt); 
  double time= 0; double endTime = 30; //dt*50; 
  int filecnt = 0; int it = 0, writeInt = 1; 
  ofstream myfile; 
  while (time < endTime) {
    cout << setiosflags(ios::scientific) << setprecision(2); 
    cout << "------------- Processing TIME = "<< time << " ------------------"<<endl; 

    auto vel = grid->getVel();
    auto diverge = grid->valDiv(vel);
    cout << "+ div=("<<diverge.min()<<":" << diverge.max()<<") "; 
    cout << " u=(" << vel.comp(0).min()<<":"<<vel.comp(0).max()<<")"; 
    cout << " v=(" << vel.comp(1).min()<<":"<<vel.comp(1).max()<<")"<< endl; 
    //    cout << grid->valDiv(vel) << endl; 

    vor->set(grid->valGrad(v).comp(0) - grid->valGrad(u).comp(1)); 
    auto gu = grid->valGrad(u); 
    auto gv = grid->valGrad(v);

    grid->lockBC(u); 
    u->solve(
             grid->ddt(rho) 
   	     + grid->div(vel, rho, {0.5}) 
  	     - grid->laplace(mu, {0.5}) 
	     - grid->source(0, -gp.comp(0)) // gradX is not defined; 
             );
    grid->unlockBC();

    grid->lockBC(v); 
    v->solve(
	     grid->ddt(rho) 
   	     + grid->div(vel, rho, {0.5})
   	     - grid->laplace(mu, {0.5})
   	     - grid->source(0, -gp.comp(1))
	     );
    grid->unlockBC();     

    //    grid->solBasedAdapt(vor->data);
    //    grid->solBasedAdapt(gp);
    if ((it == 1) || (it % writeInt) == 0) { 
      //grid->solBasedAdapt(gu, 1.0);
      //      grid->solBasedAdapt(grid->valGrad(vor), 0.4);
      grid->solBasedAdapt(vor->data, 0.9); 
      // for (auto i = 0; i < grid->listCell.size(); ++i)
      // 	for (auto j = 0; j< 3; ++j) 
      // 	  grid->listCell[i]->checkNgbrLevel(j); 
      grid->refine2(); 
    }

    // if (periodic) grid->adapt();
    auto velstar = grid->getVel();    

    auto vsdiv = grid->valDiv(velstar); 
    cout << "+ div=("<<vsdiv.min()<<":" << vsdiv.max()<<") "<<endl; 
    //cout << grid->valDiv(velstar) << endl; 
    //cin.ignore().get(); 

    grid->lockBC(p); 
    p->solve(grid->laplace(dt/rho) - grid->source(0, grid->valDiv(velstar)));
    double pref = p->data.min(); 
    for (auto i = 0; i < grid->listCell.size(); ++i) {      
       p->data[i] = p->data[i] - pref; 
    }
    grid->unlockBC(); 

    gp = grid->valGrad(p); 
    //cout << gp <<endl; 
    u->set(velstar.comp(0)-dt/rho*gp.comp(0));  // 
    v->set(velstar.comp(1)-dt/rho*gp.comp(1));  //
    //grid->correctVel(velstar, dt/rho); 

    grid->setDt(dt); 
    time += grid->dt; 

    //grid->solBasedAdapt(vor->data);
    //    grid->solBasedAdapt(gv);
    //grid->refine2(); 
    if (( it++ % writeInt) == 0) {
      //grid->writeFace("face"+std::to_string(filecnt)+".vtk"); 
      std::string flname="cav"+std::to_string(filecnt++)+".vtk"; 
      myfile.open(flname); 
      myfile << grid << endl;
      myfile.close(); 
            
    } 

    cout << "---------------------------------------------------"<<endl; 
    cin.ignore().get(); 
  }

 
  delete(grid); 

  return 0; 
}
