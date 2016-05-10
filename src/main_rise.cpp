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
  auto dt = 0.0005; 
  auto t = clock(); 

  // GRID 
  //  Block2* grid = new Block2({0, 0, 0}, {1.0, 2.0, 0}, 10, 20);
  Block2* grid = new Block2({0, 0, 0}, {1.0, 1.0, 0}, 10, 10);
  // grid->adaptCriteria(); 

  // CONST variables; sv->specific volume
  double sv1 = 0.1; double mu1 = 0.067;   double g = 0; //-1;
  double sv2 = 1;   double mu2 = 0.00067; double sigma = 1;

  // FIELD variables;
  grid->addVar({"p", "vor", "I"}); 
    
  // initial and bc values; say this is rho*u then divide it by rho 
  auto u = grid->getVar("u"); 
  u->set(0.0);
  u->setBC("west", "val", 0); u->setBC("east", "val", 0);
  u->setBC("south", "val", 0); u->setBC("north", "val", 0); 
  auto v = grid->getVar("v"); 
  v->set(0.0);
  v->setBC("west", "val", 0); v->setBC("east", "val", 0);
  v->setBC("south", "val", 0); v->setBC("north", "val", 0); 


  auto p = grid->getVar("p");
  p->setBC("south", "val", 0);  
  p->setBC("north", "val", 0);  
  p->setBC("west", "val", 0);  
  p->setBC("east", "val", 0);  
  p->set(0.0); 

  grid->levelHighBound[0] = 3;
  grid->levelHighBound[1] = 3; 
  grid->levelHighBound[2] = 0; 

  auto I = grid->getVar("I"); 
  I->set(0.0); 
  double pi = 4.0*atan(1);  
  for (auto j = 0; j < 4; ++j) {
    for (auto i = 0; i < grid->listCell.size(); ++i) {
      auto x = grid->listCell[i]->getCoord(); // - Vec3(0.5, 0.5); 
      double r = 0.125 - (x - Vec3(0.5, 0.5)).abs();       
      I->set(i, 1.0/(1.0 + exp(-2.0*80*(r)))); 
      r = 0.125 - (x - Vec3(0.5, 0.625)).abs(); 
      // I->set(i, max(0.0, min(1.0, I->get(i) + 1.0/(1.0 + exp(-2.0*80*(r))))));
    }    
    grid->solBasedAdapt2(grid->getError(I)); 
    grid->adapt(); 
  }
  auto gI = grid->valGrad(I); 
  auto gp = grid->valGrad(p);
  // grid->addVec("st"); 
  // auto stx = grid->getVar("stx"); 
  // auto sty = grid->getVar("sty"); 
  grid->writeVTK("rise"); 

  // exit(1); 
  auto vor = grid->getVar("vor"); 

  // solver behavior
  u->solver = "Gauss"; u->itmax = 200; u->tol = 1e-4; 
  v->solver = "Gauss"; v->itmax = 200; v->tol = 1e-4; 
  p->itmax = 2000; p->tol = 1e-5; 

  I->solver = "Gauss"; I->itmax = 100; I->tol = 1e-6; 

  VecX<double> rho; 
  VecX<double> mu; 
   
  dt=0.0001; 
  // Time control 
  grid->setDt(dt);  
  double time= 0; double endTime = 10.0; 
  int filecnt = 0; int it = 0, writeInt = 1; auto adaptInt = 10; 
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

    // Interface move
    grid->lockBC(I); 
    I->solve (grid->ddt(1.0) 
     	      +grid->divRK2E(vel, 1.0)
	       ); 
    grid->unlockBC(); 

    gI = grid->valGrad(I);
    // for (auto i = 0; i < gI.comp(0).size(); ++i) {
    //   double mag = 1.0/gI.abs(); 
    //   gI[i][0] *= mag; 
    //   gI[i][1] *= mag; 
    //   gI[i][2] *= mag; 
    // }
    
    // auto gK = grid->valDiv(gI);

    // double mingk = 100; double maxgk = -100; 
    // for (auto i = 0; i < gK.size(); ++i) {
    //   if (gK[i] > maxgk) maxgk = gK[i]; 
    //   if (gK[i] < mingk) mingk = gK[i]; 
    // }

    // cout << " curv=(" << mingk<<":"<<maxgk<<")"<< endl; 
    // cout << " nx=(" << gI.comp(0).min()<<":"<<gI.comp(0).max()<<")"<< endl; 
    // cout << " ny=(" << gI.comp(1).min()<<":"<<gI.comp(1).max()<<")"<< endl; 
 
    // for (auto i = 0; i < grid->listCell.size(); ++i) {
    //   gI[i][0] *= gK[i]; 
    //   gI[i][1] *= gK[i]; //sigma*(sv1 + I->get(i) * (sv2 - sv1)); 
    // }

    vor->set(grid->valGrad(v).comp(0) - grid->valGrad(u).comp(1)); 
    auto gu = grid->valGrad(u); 
    auto gv = grid->valGrad(v);

    // rho*d(u)/dt + rho*u*d(u)/dx = -dp/dx + mu*d2u/dx2 + rho*g + sigma*nx; 
    // u : rhou & v : rhov

    for (auto i = 0; i < grid->listCell.size(); ++i) {
      rho[i] = 1.0/(sv1 + I->get(i) * (sv2 - sv1)); 
      vel[i][0] *= rho[i]; 
      vel[i][1] *= rho[i]; 
    }  
    mu = mu1 + I->data * (mu2 - mu1); 

    grid->lockBC(u); 
    u->solve(
             grid->ddt(rho) 
   	     + grid->div(vel, 1, {0.5}) 
  	     - grid->laplace(mu, {0.5}) 
	     - grid->source(0, sigma*gI.comp(0)) //*(sv1 + I->data *(sv2 - sv1))) // gradX is not defined; 
             );
    grid->unlockBC();

    grid->lockBC(v); 
    v->solve(
	     grid->ddt(rho) 
   	     + grid->div(vel, 1, {0.5})
   	     - grid->laplace(mu, {0.5})
   	     - grid->source(0, sigma*gI.comp(1) - rho*g) //*(sv1 + I->data *(sv2 - sv1)))
	     );
    grid->unlockBC();     

    //    grid->solBasedAdapt(vor->data);
    //    grid->solBasedAdapt(gp);
    if ((it == 1) || (it % adaptInt) == 0) { 
      //grid->solBasedAdapt(gu, 1.0);
      //      grid->solBasedAdapt(grid->valGrad(vor), 0.4);
      grid->solBasedAdapt2(grid->getError(I)); 
      grid->solBasedAdapt(vor->data, 0.9); 
      // for (auto i = 0; i < grid->listCell.size(); ++i)
      // 	for (auto j = 0; j< 3; ++j) 
      // 	  grid->listCell[i]->checkNgbrLevel(j); 
      grid->adapt(); 
    }

    // if (periodic) grid->adapt();
    auto velstar = grid->getVel();    

    auto vsdiv = grid->valDiv(velstar); 
    cout << "+ div=("<<vsdiv.min()<<":" << vsdiv.max()<<") "<<endl; 
    //cout << grid->valDiv(velstar) << endl; 
    //cin.ignore().get(); 

    // d(rhou)/dt = -dp/dx
    // div(rhou(n+1))/dt-div(u(n))/dt = -d2p/dx2
    // dt/den*d2p/dx2 = 
    grid->lockBC(p); 
    p->solve(grid->laplace(dt*(sv1 + I->data *(sv2 - sv1)))
			   - grid->source(0, grid->valDiv(velstar))
			   );
    grid->unlockBC(); 

    gp = grid->valGrad(p); 
    for (auto i = 0; i < grid->listCell.size(); ++i) {
      gp[i][0] *= dt*(sv1 + I->get(i) *(sv2 - sv1) ); 
      gp[i][1] *= dt*(sv1 + I->get(i) *(sv2 - sv1) ); 
    }
    //cout << gp <<endl; 
    u->set(velstar.comp(0)-gp.comp(0));  // 
    v->set(velstar.comp(1)-gp.comp(1));  //
    //grid->correctVel(velstar, dt/rho1); 

    grid->setDt(dt); 
    time += dt; 

    //grid->solBasedAdapt(vor->data);
    //    grid->solBasedAdapt(gv);
    //grid->refine2(); 
    if (( it++ % writeInt) == 0) {
      //grid->writeFace("face"+std::to_string(filecnt)+".vtk"); 
      std::string flname="rise"+std::to_string(filecnt++)+".vtk"; 
      myfile.open(flname); 
      myfile << grid << endl;
      myfile.close(); 
            
    } 

    cout << "---------------------------------------------------"<<endl; 
    //cin.ignore().get(); 
  }

 
  delete(grid); 

  return 0; 
}
