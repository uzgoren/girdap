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
  auto dt = 0.5; 
  auto t = clock(); 
  //Block2* grid = new Block2({0, 0, 0}, {0.5, 0.3, 0}, 50, 30);
  Block2* grid = new Block2({0, 0, 0}, {0.02, 0.02, 0}, 5, 5);
  grid->setDt(dt); 
  double time= 0; double endTime = 10; //dt*50; 

  // Field variables; 
  double rho=10000, cp=1000, k=10; 

  // auto y = grid->getCoord(1); 
  // for (auto i = 0; i < grid->listCell.size(); ++i) {
  //   if (y[i] < 0.05) {
  //     grid->listCell[i]->adapt[1] = 2; 
  //   } else if (y[i] < 0.1) {
  //     grid->listCell[i]->adapt[1] = 1; 
  //   }
  // }

  // grid->adaptCriteria(); 

  // cout<< "pt: 1430 dir 0 and 1" << endl; // west
  // grid->debugFace(1430, 0);   grid->debugFace(1430, 1); 
  // //cout<< "pt: 1431 dir 0 and 1" << endl; // internal
  // //grid->debugFace(1431, 1); 
  // cout<< "pt: 25 dir 0 and 1" << endl; // south
  // grid->debugFace(25, 0);   grid->debugFace(25, 1);
  // cout<< "pt: 1480 dir 0 and 1" << endl; // east
  // grid->debugFace(1480, 0);   grid->debugFace(1480, 1);  
  // // cout<< "pt: 3302 dir 0 and 1" << endl; // problem
  // // grid->debugFace(3302, 0);   grid->debugFace(3302, 1); 
  // // cout<< "pt: 3301 dir 0 and 1" << endl; // problem
  // // grid->debugFace(3301, 0);   grid->debugFace(3301, 1); 
  // //cout<< "pt: 2760 dir 0 and 1" << endl; // problem
  // //grid->debugFace(2760, 0);   grid->debugFace(2760, 1); 

  grid->addVar({"T", "p"}); 
    
  auto u = grid->getVar("u"); 
  auto v = grid->getVar("v"); 
  auto w = grid->getVar("w"); 
  auto p = grid->getVar("p");
  auto T = grid->getVar("T");

  T->set(250); 

  //  T->setBC("west", "val", 250); 
  T->setBC("west", "val", 0);
  //  T->setBC("north", "grad", 250*20, -250); 
  T->setBC("south", "grad", -250*20, 250); 
  // T->setBC("north", "grad", 0); 
  // T->setBC("south", "grad", 0); 
  T->itmax = 10; 

  u->set(0.01); 
  v->set(0.03); 
  

  // if (auto u = grid->getVar("u")) {
  //   u->set({"x", "y"}, "sqrt(x^2 + y^2)*sin(atan2(y-0.5, x-0.5))"); 
  // }

  // if (auto v = grid->getVar("v")) {
  //   v->set({"x", "y"}, "-sqrt(x^2 + y^2)*cos(atan2(y-0.5, x-0.5))"); 
  // }


  // for (auto c: grid->listCell) {c->adapt[0] = 1; c->adapt[1] = 1;}
  // grid->refine();
  while (time < endTime) {
    cout << setiosflags(ios::fixed) << setprecision(2); 
    cout << "------------- Processing TIME = "<< time << " ------------------"<<endl; 

    auto vel = grid->getVel();
    //    cout << grid->valDiv(vel) << endl; 


  //   if (ibmGrid) {
  //     vel2 = grid->passTo(ibmGrid, vel);
  //     ibmGrid->move(vel2, 'RK4'); 
  //     surfTen = ibmGrid->getNorm(); 
  //     force = ibmGrid->passTo(grid, surfTen);
  //   } else {
  //     force = 0;
  //   }
  // auto axb = grid->laplace(1.0); 
  // cout << axb.A << endl; 
  // cout << axb.b << endl;
  //  u->solve(
  //           grid->ddt(rho) 
  // 	     + grid->div(rho*vel, {0, 1.5, -0.5}) 
  //	     grid->laplace(1.0) //rho/mu), {0.5, 0.5}) 
  // 	     + grid->source(0, gradX(p) - rho*aVec + force.x()) 
  //           );
    
  //   v->solve(grid->ddt(rho) 
  // 	     + (grid->div(rho*vel), {0, 1.5, -0.5}) 
  // 	     - grid->laplace(rho/mu, {0.5, 0.5}) 
  // 	     + matSource(0, gradY(p) - rho*aVec + force.y()) );
    
  //   if (periodic) grid->adapt();
    
  //   p->solve(grid->laplace(1./rho) - grid->div(rho));
    
  //   u->update(p);
  //   v->update(p); 
    
    grid->lockBC(T); 
    // T->solver="CG"; 
    T->solve(
	     grid->ddt(1.0) 
	     + grid->div(vel, 1.0)
	     - grid->laplace(k/cp/rho) 
	     - grid->source(0, 10000/cp/rho)
	     ); 
    grid->unlockBC();     
    time += grid->dt; 
    grid->setDt(dt); 
    // cin.ignore().get(); 
    cout << "---------------------------------------------------"<<endl; 
  }
  //  grid->adaptCriteria(); 

  ofstream myfile; 
  myfile.open("reg.vtk"); 
  myfile << grid << endl;
  myfile.close();    
 
  //grid->writeFace(); 
 
  grid->writeFace(); 

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
