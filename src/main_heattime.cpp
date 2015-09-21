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
  auto dt = 1; 
  auto t = clock(); 
  //Block2* grid = new Block2({0, 0, 0}, {0.5, 0.3, 0}, 50, 30);
  Block2* grid = new Block2({0, 0, 0}, {1.0, 1.0, 0}, 5, 5);
  grid->setDt(dt); 
  double time= 0; double endTime = 2; //dt*50; 

  // Field variables; 
  double rho=10000, cp=1000, k=10; 

  grid->addVar({"T"}); 
    
  auto T = grid->getVar("T");
  auto u = grid->getVar("u");
  auto v = grid->getVar("v");
  
  T->set(0.0); 

  T->setBC("west", "val", 100);
  T->setBC("east", "grad", 100); 
  T->setBC("south", "grad", -250*20, 250);   
  T->itmax = 1000; 
  T->tol = 1e-6;

  // u->set(0.01); 
  // v->set(0.03); 

  int filecnt = 0; int it = 0, writeInt = 1; 
  ofstream myfile;   
  // while (time < endTime) {
    cout << setiosflags(ios::fixed) << setprecision(6); 
    cout << "------------- Processing TIME = "<< time << " ------------------"<<endl; 

    auto vel = grid->getVel();
    
    grid->lockBC(T); 
    T->solve(
	     //grid->ddt(1.0) 
	     //	     + grid->div(vel, 1.0)
	     - grid->laplace(1.0) 
	     // - grid->source(0, 10000/cp/rho)
	     ); 
    grid->unlockBC();     

    // auto gt = grid->valGrad(T); 
    // grid->solBasedAdapt(gt, 1); 
    // grid->refine2(); 

    time += grid->dt; 
    //    grid->setDt(dt); 

    //if (( it % writeInt) == 0) {
      std::string flname="heat"+std::to_string(filecnt++)+".vtk"; 
      myfile.open(flname); 
      myfile << grid << endl;
      myfile.close();   
      //} 

    cout << "---------------------------------------------------"<<endl; 
    //  }


  // myfile.open("reg.vtk"); 
  // myfile << grid << endl;
  // myfile.close();    
 
  //grid->writeFace(); 
 
//  grid->writeFace(); 

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
