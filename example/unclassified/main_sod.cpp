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

//*************************************************************
// Codelet: Performance analysis of Linear system construction  
//*************************************************************

#include <linSolve/MatY.hpp>
#include <linSolve/LinSys.hpp>
#include <grid2/Grid>
#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>
#include <time.h>
#include <field/Var>
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */



double heavy(double x, double e) {
  if (x/e <= -1)
    return 0; 
  else if (x/e <= 0)
    return 0.5 + x/e + 0.5*pow(x/e,2); 
  else if (x/e <= 1)
    return 0.5 + x/e - 0.5*pow(x/e,2); 
  else 
   return 1; 
}

double heavy(Vec3 x, Vec3 xp, double e) {
  return heavy((x - xp).abs(), e); 
}


int main(int argc,char *argv[]) {
  int_8 N, l; 
  float cfl = 0.5; 
  clock_t t; 
  
  if (argc == 1) { N = 20; l = 0;}
  else if (argc == 2) {N = atoi(argv[1]); l = 0;}
  else if (argc == 3) {N = atoi(argv[1]); l = atoi(argv[2]);} 
  else {N = atoi(argv[1]); l = atoi(argv[2]); cfl = atof(argv[3]); }  
  
  Block2* grid = new Block2({0, 0, 0}, {1, 1, 0}, N, 1); 
  grid->cfl = cfl; 
  grid->levelHighBound[0] = l; 
  grid->levelHighBound[1] = 0;
  grid->levelHighBound[2] = 0; 
  
  grid->addVar({"rho", "E", "p", "pu"}); 
  
  auto rho = grid->getVar("rho"); 
  auto E = grid->getVar("E"); 
  auto u = grid->getVar("u"); 
  auto p = grid->getVar("p"); 
  auto pu = grid->getVar("p"); 
 
  rho->setBC("west", "val", 1); 
  rho->setBC("east", "val", 0.125);
  p->setBC("west", "val", 1); 
  p->setBC("east", "val", 0.1); 
  E->setBC("west", "val", 2.5); 
  E->setBC("east", "val", 0.25);
  u->setBC("west", "val", 0); 
  u->setBC("east", "val", 0); 
  pu->setBC("west", "val", 0); 
  pu->setBC("east", "val", 0); 
  
  pu->set(0); 
  u->set(0); 

  for (auto j = 0; j < l; ++j) {
    for (auto i = 0; i < grid->listCell.size(); ++i) {
      auto x = grid->listCell[i]->getCoord();
      if (x[0] <= 0.5) { 
	rho->set(i, 1.0); 
	p->set(i, 1.0); 
	E->set(i, 2.5); 
      } else {
	rho->set(i, 0.125); 
	p->set(i, 0.1); 
	E->set(i, 0.25);  
      }
    } 
    grid->solBasedAdapt(grid->valGrad(rho)); 
    grid->adapt(); 
  }

  grid->writeVTK("H"); 

  double endTime = 0.1; //0.8/u->data.max(); 

  auto vel = grid->getVel(); 
  double time = 0; 
  auto writeTime = 0.0; 
  auto writeCnt = writeTime; 
  while (time <= tend) {
    grid->setDt(0.001);       
    
    time += grid->dt;    
    cout << "Time: " << time << " " << grid->dt << " --> " << "xp: " << xp << "xm: "<< xm << endl;
    // Density
    grid->advectDiv(vel, 1); 

    grid->lockBC(u); 
    u->solve( grid->ddt(1) + grid->divRK2E(vel*rho, 1) 
	      + grid->source(0, grid->valGrad(p).comp(0) )  );
    grid->unlockBC(); 

    pu->set(p->data*u->data);

    grid->lockBC(E); 
    E->solve( grid->ddt(1.0) + grid->divRK2E(vel, 1.0) 
	      + grid->source(0, grid->valGrad(pu).comp(0) ) );
    grid->unlockBC(); 

    grid->solBasedAdapt(grid->valGrad(rho)); 
    grid->solBasedAdapt(grid->valGrad(E)); 
    grid->solBasedAdapt(grid->valGrad(u)); 
    grid->adapt(); 

    for (auto c: grid->listCell) {
      p->set(c->id, 0.4*(E->get(c->id) 
			 - 0.5*rho->get(c->id)*pow(u->get(c->id), 2.0)));
    }

    if (writeCnt <= 0 || time >= endTime) {
      grid->writeVTK("sod_");  
      writeCnt = writeTime; 
    }

  }

  // NEW
  // for (auto f:grid->listFace) {
  //   int n = f->next; 
  //   int p = f->prev; 
    
  //   auto xf = f->getCoord(); 
  //   double uf, phif; 
  //   auto nid = grid->searchVertexbyCoords(xf, f->node[0]); 
  //   uf = grid->listVertex[nid]->evalPhi(u, &xf);
  //   Vec3 x0 = xf - Vec3(0.5*grid->dt*uf,0,0); 
  //   nid = grid->searchVertexbyCoords(x0, nid); 
  //   phif = grid->listVertex[nid]->evalPhi(phi, &x0); 
  //   cout << phif << " "; 
  // }



  // // Speed test: Matrix construction and solution old and new techniques  
  // LinSys old(N); VecX<double> x0(N); VecX<double> err0(N); 
  // old.x = &x0; 
  // old.error = &err0; 

  // t = clock(); 
  // for (auto i =0; i < N; ++i) {
  //   old.A[i][i] = 4; 
  //   if (i > 0) old.A[i][i-1] = -1; 
  //   if (i < N-1) old.A[i][i+1] = -1; 
  //   if (i-N >= 0) old.A[i][i-N] = -1; 
  //   if (i+N < N) old.A[i][i+N] = -1; 
  //   old.b[i] = 3*i+1; 
  // }
  // old.setVECX(old.A*old.b); 
  // t = clock() - t; 
  // float setup0 = float(t)/CLOCKS_PER_SEC; 

  // t = clock(); 
  // old.BiCGSTAB(); 
  // t = clock() - t; 
  // float solve0 = float(t)/CLOCKS_PER_SEC; 

  // //cout << *old.x << endl; 
  // cout << old.error->abs() << endl; 

  // triLinSys sp(N); VecX<double> x1(N); VecX<double> err1(N); 
  // sp.x = &x1; 
  // sp.error = &err1; 

  // t = clock(); 
  // for (auto i =0; i < N; ++i) {
  //   sp.A += {(double)i, (double)i, 4}; //old.A[i][i] = 2; 
  //   if (i > 0) sp.A += {(double)i, (double)i-1, -1}; //old.A[i][i-1] = -1; 
  //   if (i < N-1) sp.A += {(double)i, (double)i+1, -1}; //old.A[i][i+1] = -1; 
  //   if (i-N >= 0) sp.A += {(double)i, (double)i-N, -1}; //old.A[i][i-1] = -1; 
  //   if (i+N < N) sp.A += {(double)i, (double)i+N, -1}; //old.A[i][i+1] = -1; 

  //   sp.b[i] = 3*i+1; 
  // }
  // sp.setMat(sp.A); 
  // sp.b = sp.S*sp.b; 
  // t = clock()-t; 
  // float setup1 = float(t)/CLOCKS_PER_SEC; 

  // t = clock(); 
  // sp.BiCGSTAB(); 
  // t = clock()-t; 
  // float solve1 = float(t)/CLOCKS_PER_SEC; 
  // cout << " ----------- " << endl ; 
  // //cout << *sp.x << endl; 
  // cout << sp.error->abs() << endl; 

  // cout << "OLD: setup-> " << setup0 << " sec, solve-> " << solve0 << "sec."<<endl; 
  // cout << "NEW: setup-> " << setup1 << " sec, solve-> " << solve1 << "sec."<<endl; 

  // // cout << "Abs" << endl; 
  // Triplets vt, bt; 
  
  // triLinSys Ab; 
  // vt += {{0,0,4}, {0,1,-1}, {1,2,-1}, {2,1,-1}}; 
  // vt += {{2,2,4}, {3,2,-1}, {3,3,4}, {3,4,-1}}; 
  // bt += {{4,4,4}, {4,3,-1}, {1,1,4}, {2,3,-1}};
  // bt += {{1,0,-1}, {2,1,-2}, {1,4,-6}, {2,2,1}};
  // bt(2,2) = 7; 

  // Ab.setMat(-(2*vt-bt)); 
  
  // for (auto i = 0; i < Ab.S.aij.size(); ++i) { 
  //    cout << i << " " << Ab.S.aij[i] << " " << Ab.S.val[i] << endl; 
  // }
  
  // // VecX<double> b(5); 
  // // cout << "----"<< endl<< Ab.A.rank << " " << b.size() << endl; 
  // // for (auto i = 0; i < Ab.A.rank; ++i) b[i] = i+1; 
  
  // // Ab.b = Ab.A*b;
  // // old.setVECX(Ab.b); 
  
  // // VecX<double> x(b.size());
  // // Ab.error = &b; 
  // // Ab.x = &x;
  // // Ab.BiCGSTAB(); 
  
  // // cout << *Ab.x << endl; 
  // // cout << Ab.error->abs() << endl; 

  // // cout << old.A<<endl; 
  // // cout << old.b<<endl; 
  // // old.error = &b; 
  // // old.x = &x;

  // // old.BiCGSTAB(); 
  
  // // cout << *old.x << endl; 
  // // cout << old.error->abs() << endl; 
  
  

  return 0; 
};




