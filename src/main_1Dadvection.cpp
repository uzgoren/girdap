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
  
  double L = 20; 

  Block2* grid = new Block2({0, 0, 0}, {L, 0.1*L, 0}, N, 1); 
  grid->cfl = cfl; 
  grid->levelHighBound[0] = l; 
  grid->levelHighBound[1] = 0;
  grid->levelHighBound[2] = 0; 
  
  grid->addVar({"phi", "phi2", "ref"}); 
  
  auto u = grid->getVar("u"); 
  auto v = grid->getVar("v"); 
 
  
  auto phi = grid->getVar("phi"); 
  auto phi2 = grid->getVar("phi2"); 
  auto ref = grid->getVar("ref"); 

  phi->set(0); 
  u->set(0.1); 
  v->set(0.0); 

  VecX<double> u0(grid->listCell.size(), 0.0); 
  auto pi= 4*atan(1.0); 
  for (auto i = 0; i < grid->listCell.size(); ++i) {
    auto x = grid->listCell[i]->getCoord(); 
    u0[i] = 1 + sin(pi*x[0]/L);
  } 
  u->set(u0); 

  phi->itmax = 10; 
  phi->tol = 1e-8; 
  phi->solver = "Gauss";

  double xp = 1.8;
  double xm = 6.0; 

  double e = 1.0/N/pow(2,grid->levelMax[0])*(2); 
  for (auto c:grid->listCell) {
    auto h= heavy(c->getCoord().x()-xp, e)-heavy(c->getCoord().x()-xm, e); 
    phi->set(c->id, h);
    phi2->set(c->id, h); 
    ref->set(c->id, h); 
  }
  
  

  for (auto j = 0; j < l; ++j) {
    grid->solBasedAdapt2(grid->getError(ref), 1e-5, 1e-3); 
    //grid->valAdapt(phi, 0.1, 0.9); 
    grid->adapt(); 
    double e = 1.0/N/pow(2,grid->levelMax[0])*(2); 
    for (auto c:grid->listCell) { 
      auto h= heavy(c->getCoord().x()-xp, e)-heavy(c->getCoord().x()-xm, e); 
      phi->set(c->id, h);
      phi2->set(c->id, h);
      ref->set(c->id, h);  
    }    
  }

    grid->writeVTK("H"); 

    double endTime = L; //0.8/u->data.max(); 
    double tend = endTime; 
    //    endTime = 1; 

    auto vel = grid->getVel(); 
    double time = 0; 
    auto writeTime = 0.1; 
    auto writeCnt = writeTime; 
    while (time <= tend) {
      for (auto i = 0; i < grid->listCell.size(); ++i) {
	auto x = grid->listCell[i]->getCoord(); 
	u->set(i, (1.0 + sin(pi*x.x()/L))*cos(2*pi*time/tend));
      }            
      cout << u->data.max() << " " << u->data.min() <<endl; 
      grid->setDt(2.0);       

      time += grid->dt; 

      auto up = grid->interp(u, xp); 
      auto um = grid->interp(u, xm); 

      xp += up*grid->dt; 
      xm += um*grid->dt; 
     
      for (auto c:grid->listCell) { 
	auto h= heavy(c->getCoord().x()-xp, e)-heavy(c->getCoord().x()-xm, e); 
	ref->set(c->id, h);  
      }    

      for (auto i = 0; i < grid->listCell.size(); ++i) {
	auto x = grid->listCell[i]->getCoord(); 
	phi2->set(i, grid->interp(phi2, x - Vec3(u->get(i), v->get(i), 0)*grid->dt)); 
      }

      writeCnt -= grid->dt; 
      
    auto vel = grid->getVel();
   
    cout << "Time: " << time << " " << grid->dt << " --> " << "xp: " << xp << "xm: "<< xm << endl;
    grid->advanceDiv(phi, vel); //, true); 
    // FOU
    // grid->lockBC(phi); 
    // phi->solve(grid->ddt(1.0)+grid->divRK2E(vel,1.0)); 
    // grid->unlockBC(); 

    // for (auto f:grid->listFace) {
    //   int n = f->next; 
    //   int p = f->prev;  

    //   Vec3 norm = f->vol();
    //   auto xf = f->getCoord(); 
    //   double phif; 
    //   auto nid = grid->searchVertexbyCoords(xf, f->node[0]); 
    //   auto uf1 = Vec3(grid->listVertex[nid]->evalPhi(u, &xf),0,0);
    //   // Vec3 x1 = xf - 0.125*grid->dt*uf1;
    //   // nid = grid->searchVertexbyCoords(x1, nid); 
    //   // auto uf2 = Vec3(grid->listVertex[nid]->evalPhi(u, &x1),0,0);
    //   // Vec3 x2 = xf - 0.250*grid->dt*uf2;
    //   // nid = grid->searchVertexbyCoords(x2, nid); 
    //   // auto uf3 = Vec3(grid->listVertex[nid]->evalPhi(u, &x2),0,0);
    //   // Vec3 x3 = xf - 0.5*grid->dt*(uf1 - 2*uf2 + 2*uf3); 
    //   // nid = grid->searchVertexbyCoords(x3, nid); 
    //   // auto uf4 = Vec3(grid->listVertex[nid]->evalPhi(u, &x3),0,0);
    //   // Vec3 x0 = xf - grid->dt/12*(uf1 + 4*uf3 - uf4); 
    //   Vec3 x0 = xf - grid->dt*0.5*uf1; 
    //   auto uf = uf1; 

    //   double flux = uf*norm;
    //   //      if (p>=0 && n>=0) {   
    // 	// FOU
    // 	// if (flux > 0) 
    // 	//   flux *= phi->get(p); 
    // 	// else
    // 	//   flux *= phi->get(n);    

    //   //      nid = grid->searchVertexbyCoords(x0, nid); 
    //   phif = max(0.0, min(1.0, grid->interp(phi, x0, nid))); //listVertex[nid]->evalPhi(phi, &x0)));	

    //   auto invvoln = (n < 0) ? 1.0 : 1.0/grid->listCell[n]->vol().abs(); 
    //   auto invvolp = (p < 0) ? 1.0 : 1.0/grid->listCell[p]->vol().abs();      
      
    //   if (n >= 0) phi->data[n] += flux*phif*grid->dt*invvoln; 
    //   if (p >= 0) phi->data[p] -= flux*phif*grid->dt*invvolp; 
    //   //      } 
    // }

    grid->solBasedAdapt2(grid->getError(phi)); //, 1e-5, 1e-3);// grid->valGrad(phi)); //
    //grid->valAdapt(phi, 0.1, 0.9); 
    grid->adapt(); 
    if (writeCnt <= 0 || time >= endTime) {
      grid->writeVTK("H");  
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




