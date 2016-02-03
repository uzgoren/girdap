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
  auto dt = 0.5; auto writeTime = 0.02; 
  auto t = clock(); int iter = 0; 
  Block2* grid = new Block2({0, 0, 0}, {1, 1, 0}, 10, 10);
  //Block2* grid = new Block2({0, 0, 0}, {0.02, 0.02, 0}, 10, 1);
  double time= 0; double endTime = 8; //dt*50; 

  // Field variables; 
  double rho=10000, cp=1000, k=1; 

  grid->addVar({"T"}); 
  grid->addVar({"msx", "msy"}); 
    
  auto T = grid->getVar("T");
  auto u = grid->getVar("u");
  auto v = grid->getVar("v");
  auto msx = grid->getVar("msx"); 
  auto msy = grid->getVar("msy");   

  T->solver = "Gauss"; 
  
  T->set(0); 
    double pi = 4.0*atan(1); 
  
  for (auto j = 0; j < 5; ++j) {
    for (auto i = 0; i < grid->listCell.size(); ++i) {
      auto x = grid->listCell[i]->getCoord(); // - Vec3(0.5, 0.5); 
      u->set(i, -2*sin(pi*x[1])*cos(pi*x[1])*sin(pi*x[0])*sin(pi*x[0]));
      v->set(i, 2*sin(pi*x[0])*cos(pi*x[0])*sin(pi*x[1])*sin(pi*x[1])); 

      auto r = (grid->listCell[i]->getCoord() - Vec3(0.5, 0.75)).abs(); 

      T->set(i, 1.0/(1.0 + exp(-2.0*80*(0.15-r)))); 
    }
    //    auto gt = grid->valGrad(T); 
    //grid->solBasedAdapt(gt); 
    grid->solBasedAdapt2(grid->getError(T));
    grid->adapt(); 
  }

  double mass0=0; double mass=0; 
  for (auto i=0; i < grid->listCell.size(); ++i) {
    mass0 += grid->listCell[i]->vol().abs()*T->get(i); 
  }

  for (auto i = 0; i < grid->listCell.size(); ++i) {
    auto xlevel = grid->listCell[i]->level[0]; 
    auto ylevel = grid->listCell[i]->level[0];
    msx->set(i, grid->listCell[i]->masterx[xlevel]); 
    msy->set(i, grid->listCell[i]->mastery[ylevel]); 
  }

  int filecnt = 0; int it = 0, writeInt = 1; 
  ofstream myfile;   
  std::string flname="heat"+std::to_string(filecnt++)+".vtk"; 
  myfile.open(flname); 
  myfile << grid << endl;
  myfile.close();   

  T->setBC("west", "val", 0);
  T->setBC("south", "val", 0); 
  T->setBC("east", "val", 0); 
  T->setBC("north", "val", 0); 
  // T->setBC("south", "grad", -20, 100);   
  // T->setBC("north", "grad", -20, 100);  
  T->itmax = 1000; 
  T->tol = 1e-6;  

  // u->set(1); 
  // v->set(1); 


  auto writeCnt = writeTime; 
  while (time < endTime) {
    iter++; 
    for (auto i = 0; i < grid->listCell.size(); ++i) {
      auto x = grid->listCell[i]->getCoord(); // - Vec3(0.5, 0.5); 
      u->set(i, -2*sin(pi*x[1])*cos(pi*x[1])*sin(pi*x[0])*sin(pi*x[0])*cos(pi*time/endTime));
      v->set(i, 2*sin(pi*x[0])*cos(pi*x[0])*sin(pi*x[1])*sin(pi*x[1])*cos(pi*time/endTime)); 
    }
    grid->setDt(2.0*dt); 

 
    cout << setiosflags(ios::fixed) << setprecision(6); 
    cout << "------------- Processing TIME = "<< time << " ------------------"<<endl; 
    cout << "  Mass: " << mass << " Initial: " << mass0 << endl; 

    auto vel = grid->getVel();
    
    grid->lockBC(T); 
    T->solve (grid->ddt(1.0) 
     	      +grid->divRK2E(vel, 1.0)
	      //- grid->laplace(k/cp/rho) 
	      //- grid->source(0, 100000/cp/rho)
	      //- grid->laplace(1.0) 
	      //- grid->source(-25, 25*20)
	       ); 
    grid->unlockBC(); 

    // remove over-shoot and undershoot
    //    for (auto i = 0; i < grid->listCell.size(); ++i) { 
    //}

    mass=0; double part=0; 
    for (auto i=0; i < grid->listCell.size(); ++i) {
      double vol = grid->listCell[i]->vol().abs();
      //      T->set(i, max(min(T->get(i)+1e-5, 1.0), 0.0)); 
      T->set(i, min(T->get(i)+1e-4, 1.0));
      T->set(i, max(T->get(i)-1e-4, 0.0)); 
      double Tval = T->get(i); 
      mass += vol*Tval;
      if (Tval > 0) part += vol; 
    }

    // double corr = (mass - mass0)/part; mass = 0; 
    // for (auto i=0; i < grid->listCell.size(); ++i) { 
    //   double vol = grid->listCell[i]->vol().abs();
    //   if (T->get(i) > 0) 
    // 	T->set(i, max(min(T->get(i)+corr*vol, 1.0), 0.0)); 
    //   mass += vol*T->get(i); 
    // }
	     
    
    //    auto gt = grid->valGrad(T); 
    if (iter%1 == 0) {
      grid->solBasedAdapt2(grid->getError(T)); 
      grid->adapt(); 
    }

    for (auto i = 0; i < grid->listCell.size(); ++i) {
      auto xlevel = grid->listCell[i]->level[0]; 
      auto ylevel = grid->listCell[i]->level[0];
      msx->set(i, grid->listCell[i]->masterx[xlevel]); 
      msy->set(i, grid->listCell[i]->mastery[ylevel]); 
    }

    time += grid->dt; 
    writeCnt -= grid->dt; 
    grid->setDt(dt); 

    if (writeCnt <= 0 || time >= endTime) {
      std::string flname="heat"+std::to_string(filecnt++)+".vtk"; 
      myfile.open(flname); 
      myfile << grid << endl;
      myfile.close();   
      writeCnt = writeTime; 
    } 

    cout << "---------------------------------------------------"<<endl; 
  }

  delete(grid); 


  // LinSys ls(6); 
  // VecX<double> x(6); x.uncompress();
  // ls.A = {
  //   { 2, -1,  0,  0,  0,  0}, 
  //   {-1,  2, -1,  0,  0,  0}, 
  //   { 0, -1,  2, -1,  0,  0}, 
  //   { 0,  0, -1,  2, -1,  0}, 
  //   { 0,  0,  0, -1,  2, -1},
  //   { 0,  0,  0,  0, -1,  2}    
  // };
  // // ls.A[5] = 
  // // //vel.A[5][4] = -0.5; 

  // ls.A.info(); 
  // ls.b.info(); 

  // x = {1, -2, 3, -4, 9, 3}; 
  // cout << ls.A << endl; 
  // cout << ls.b << endl; 
  // cout << x<<endl; 
  // ls.b = (ls.A)*(x);
  // cout << ls.b << endl;
  // x = {0, 0, 0, 0, 0, 0}; 

  // ls.x = &x;
  // ls.setLimits(1e-8); 

  // // //vel.Gauss_Seidel();
  // ls.BiCGSTAB();
  // cout <<"x: "<< x<< endl;



  return 0; 
}
