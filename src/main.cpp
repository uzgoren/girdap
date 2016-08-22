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
#include <fstream>

#include <girdap>


int main() {
  Block2 grid({0, 0, 0}, {1, 1, 0}, 20, 20);
  Block1 lin({0.2, 0.3, 0}, {0.8, 0.6, 0}, 20); //checks
  lin.writeVTK("lin");
  
  Block1 poly({{0.1, 0.1}, {0.2, 0.3}, {0.6, 0.2}
      , {0.3, 0.2}, {0.4, 0.15}, {0.3, 0.09}, {0.1,0.1}}, 0.02); // check resolve
  poly.writeVTK("poly");
  
  Block1 circ(new Geo1Circle(Vec3(0.5, 0.75), 0.15, 360, 0), 4); // OK
  circ.writeVTK("circ");
  
  poly.add(circ); // OK
  poly.writeVTK("poly");
  
  // //--- Circular arc ----
  Block1 arc(new Geo1Circle(Vec3(0.5, 0.75), 0.15, 30, 120), 40); // OK
  arc.writeVTK("arc");
  
  // //--- Rotated polygon "equal sided"
  Block1 rotsqr(new Geo1Circle(Vec3(0.5,0.75), 0.15, 30, 400), 4); // OK
  rotsqr.resolve(0.15/sqrt(2)*0.1);
  rotsqr.writeVTK("rotsqr"); 
  
  Block1 rotpenta(new Geo1Circle(Vec3(0.5,0.75), 0.15, 30, 400), 5); //OK
  rotpenta.resolve(0.01);
  rotpenta.writeVTK("rotpenta");

  lin.add(0.0, 1.0, 10
	       , [](double t) -> Vec3 {return Vec3(0.0) + t*Vec3(0.2, -0.1)/10;} );
  lin.writeVTK("lin");

  Block1 sine;
  sine.add(0.0, 1.0, 50
		, [](double t) -> Vec3 {return Vec3(t, 0.1*sin(5*t*3.14), 0);} );
  sine.writeVTK("sine"); 

  Block1 arc2;
  arc2.add(0, 360, 20
	      , [](double t) -> Vec3 {return Vec3(0.5,0.75) + 0.15*Vec3(cos(t*3.14/180), sin(t*3.14/180)); } ); 
  arc2.writeVTK("ccc");

  Block1 line2(new Geo1Sine(Vec3(0.2, 0.4), Vec3(0.5, 0.5), 0.1, 5), 50);
  line2.add(new Geo1Line(Vec3(0.5, 0.5), Vec3(0.8, 0.6)), 20); 
  line2.writeVTK("sdsd"); 

  Geo1Circle* g0 = new Geo1Circle(Vec3(0.5, 0.5), 0.2);  
  Block1 arc3(g0, 30);
  delete g0; 
  arc3.writeVTK("newCircle"); 
  
  // grid->addVertex({ {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0} }); 

  // grid->addCell( {0, 1, 2, 3} ) ; 

  // for (auto i =0; i<3; ++i) {
  //   grid->listCell[0]->adapt = {1, 1}; 
  //   grid->adapt(); 
  //   grid->writeVTK("myFirstGrid_"); 
  // }
  
  //delete(grid); 

  // Block2* volgrid = new Block2({0,0,0}, {1,1,0}, 10, 10); 
  // Grid* surf = new Grid(); 

  // surf->addVertex({{0.5,0.4}, {0.6,0.5}, {0.5,0.6}, {0.4,0.5}}); 
  // surf->addCell({{0,1}, {1,2}, {2,3}, {3,0}}); 

  // volgrid->writeVTK("vol"); 
  // surf->writeVTK("surf"); 

  // delete(volgrid); 
  // delete(surf); 
  //double pi = 4*atan(1.0); 

  // Block2* volgrid = new Block2({0,0,0}, {1,1,0}, 5, 5); 
  
  // // add a new variable
  // volgrid->addVar("f"); 
  // auto f = volgrid->getVar("f"); // variable handle
  
  // f->setBC("east", "grad", 0);   // This is the default
  // f->setBC("north", "val", 1);   //   

  // for (auto i=0; i < 4; ++i) { 
  //   for (auto c : volgrid->listCell) {
  //     auto x = c->getCoord(); // cell-centers
  //     f->set(c->id, sin(3*pi*x[0])*cos(2*pi*x[1]));
  //   }
  //   volgrid->solBasedAdapt2(volgrid->getError(f)); 
  //   volgrid->adapt(); 
  //   volgrid->writeVTK("field_"); 
  // }
  
  // delete(volgrid); 


  // Block2* volgrid = new Block2({0,0,0}, {1,1,0}, 50, 50); 

  // // Velocity field
  // auto uv = volgrid->getVar("u"); auto vv = volgrid->getVar("v"); 
  // uv->set(1.0); // set velocity
  // vv->set(-0.5); // set velocity
  // // New variable at cell center
  // volgrid->addVar("f"); auto f = volgrid->getVar("f"); 

  // Grid* surf = new Grid(); 

  // surf->addVertex({{0.55,0.32}, {0.58,0.5}, {0.45,0.68}, {0.42,0.46}}); 
  // surf->addCell({{0,1}, {1,2}, {2,3}, {3,0}}); 
  // // Refine cell; 
  // for (auto i=0; i<4; ++i) {
  //   for (auto c: surf->listCell) if (c->vol().abs() > 0.02) c->adapt[0] = 1;
  //   surf->adapt(); 
  // }
  // volgrid->updateOtherVertex(surf);
  // // mark location of this surface
  // volgrid->indicator(surf, f);

  // // Assign velocity variables to surface at vertex  
  // surf->addVec("u",1);

  // // Get velocity on the surface
  // auto us = surf->getVar("u"); auto vs = surf->getVar("v");   
  // volgrid->passVar(surf, uv, us); 
  // volgrid->passVar(surf, vv, vs);   

  // volgrid->writeVTK("vol"); 
  // surf->writeVTK("surf"); 

  // delete(volgrid); 
  // delete(surf); 

  // // Problem parameters
  // auto k = 2.0; auto qdot = 5e3; auto h = 50; auto Tinf = 20;
  // // Grid
  // Block2* grid = new Block2({0, 0, 0}, {1, 1, 0}, 10, 10); 
  // grid->levelHighBound[0] = 2; 
  // grid->levelHighBound[1] = 2; 
  // grid->addVar("T"); 
  // // Variables
  // auto T = grid->getVar("T");
  // // Linear solver
  // T->solver = "BiCGSTAB";
  // T->itmax = 1000; 
  // T->set(100); 
  // // Boundary conditions
  // T->setBC("south", "grad", 0); 
  // T->setBC("north", "grad", h/k*Tinf, -h/k);
  // T->setBC("east", "val", 200); 
  // T->setBC("west", "val", 100); 
  
  // for (auto i = 0; i< 4; ++i) {
  //   grid->solBasedAdapt2(grid->getError2(T), 2e-3, 2e-1);
  //   grid->adapt();  
    
  //   // Equation 
  //   grid->lockBC(T); 	
  //   T->solve( grid->laplace(k) 
  // 	      + grid->source(0, qdot) ); 
  //   grid->unlockBC();
    
  //   grid->writeVTK("heat"); 
  // }
    
  // delete(grid); 
  
  // Block2* grid = new Block2({0, 0, 0}, {1, 1, 0}, 10, 10);
  // double time= 0; double endTime = 1; //dt*50; 
  // // Problem constants
  // auto k = 2.0; auto qdot = 5e4; auto h = 20; auto Tinf = 20;
  // auto rho=1000, cp=4000;
  // // Field variables, velocity already defined
  // grid->addVar("T"); 

  // auto u = grid->getVar("u"); 
  // auto v = grid->getVar("v"); 
  // auto T = grid->getVar("T"); 

  // T->set(100); 
  // T->setBC("west", "val", 20);
  // T->setBC("south", "val", 20);
  // T->itmax = 50; 

  // u->set(1); 
  // v->set(0.2); 

  // grid->cfl = 0.5; 
  // int it = 0; 
  // while (time < endTime) {
  //   grid->setDt(0.5); // CFL condition     
  //   auto vel = grid->getVel(); // freeze velocity! 

  //   // Advection-diffusion equation ---
  //   grid->lockBC(T); 
  //   T->solve(
  // 	     grid->ddt(1.0) 
  // 	     + grid->div(vel, 1.0)
  // 	     - grid->laplace(k/rho/cp) 
  // 	     - grid->source(0, qdot/rho/cp)
  // 	     ); 
  //   grid->unlockBC();     
  //   grid->writeVTK("heattime_"); 
    
  //   if (it++ % 5 == 0) {
  //     grid->solBasedAdapt(grid->valGrad(T));
  //     grid->adapt(); 
  //   }

  //   time += grid->dt; 
  // }


  // // GRID 
  // Block2* grid = new Block2({0, 0, 0}, {10, 1, 0}, 20, 20);

  // // CONST variables; 
  // double rho = 1; double mu = 0.01; 

  // // FIELD variables; Vorticity to appear in the output
  // grid->addVar({"p", "vor"}); 
    
  // // initial and bc values; 
  // auto u = grid->getVar("u"); 
  // u->set(0.0);
  // u->setBC("west", "val", 1.0); u->setBC("east", "val", 0);
  // u->setBC("south", "val", 0); u->setBC("north", "val", 0); 
 
  // auto v = grid->getVar("v"); 
  // v->set(0.0);
  // v->setBC("west", "val", 0); v->setBC("east", "val", 0);
  // v->setBC("south", "val", 0); v->setBC("north", "val", 0); 

  // auto p = grid->getVar("p");
  // p->set(0.0);

  // auto vor = grid->getVar("vor");

  // // Solver behavior
  // u->solver = "Gauss"; u->itmax = 20; u->tol = 1e-4; 
  // v->solver = "Gauss"; v->itmax = 20; v->tol = 1e-4; 
  // p->solver = "BiCGSTAB"; p->itmax = 10000; p->tol = 1e-6; 

  // auto gp = grid->valGrad(p);   
 
  // // Time control 
  // grid->setDt(1.0); 
  // double time= 0; double endTime = 10; 
  // int it = 0, writeInt = 4; 

  // while (time < endTime) {
  //   auto vel = grid->getVel();
  //   // Compute vorticity; 
  //   vor->set(grid->valGrad(v).comp(0) - grid->valGrad(u).comp(1)); 
    
  //   // Solve for u*
  //   grid->lockBC(u); 
  //   u->solve(grid->ddt(rho) + grid->div(vel, rho) - grid->laplace(mu) );
  //   grid->unlockBC();

  //   // Solve for v*
  //   grid->lockBC(v); 
  //   v->solve(grid->ddt(rho) + grid->div(vel, rho) - grid->laplace(mu) ); 
  //   grid->unlockBC();     

  //   // Adapt grid using vorticity (0.9 sigma up/down for refine/coarsen); 
  //   if (it == 1 || (it % writeInt == 0)) { 
  //     grid->solBasedAdapt(vor->data, 0.9); 
  //     grid->adapt(); 
  //   }

  //   // Get Vel*
  //   auto velstar = grid->getVel();    

  //   // Solve for pressure Poisson equation
  //   grid->lockBC(p); 
  //   p->solve(grid->laplace(1.0/rho) - grid->source(0, grid->valDiv(velstar)/dt));
  //   grid->unlockBC(); 

  //   // Correct velocities; 
  //   u->set(velstar.comp(0)-dt/rho*gp.comp(0));  // 
  //   v->set(velstar.comp(1)-dt/rho*gp.comp(1));  //

  //   // Set new time step 
  //   grid->setDt(dt); 
  //   time += grid->dt; 

  //   //Write output at intervals; 
  //   if ( ( it++ % writeInt) == 0) {
  //     grid->writeVTK(); //
  //   } 

  // }

 
  // delete(grid);
  


  return 0; 
};
