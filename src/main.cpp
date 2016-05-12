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
#include <linSolve/MatY.hpp>
#include <linSolve/LinSys.hpp>
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */


int main(int argc,char *argv[]) {
  int_8 N; 
  clock_t t; 

  if (argc == 1) N = 5; 
  else N = atoi(argv[1]); 

  cout << N << endl; 

  // Speed test: Matrix construction and solution old and new techniques  
  LinSys old(N); VecX<double> x0(N); VecX<double> err0(N); 
  old.x = &x0; 
  old.error = &err0; 

  t = clock(); 
  for (auto i =0; i < N; ++i) {
    old.A[i][i] = 4; 
    if (i > 0) old.A[i][i-1] = -1; 
    if (i < N-1) old.A[i][i+1] = -1; 
    if (i-N >= 0) old.A[i][i-N] = -1; 
    if (i+N < N) old.A[i][i+N] = -1; 
    old.b[i] = 3*i+1; 
  }
  old.setVECX(old.A*old.b); 
  t = clock() - t; 
  float setup0 = float(t)/CLOCKS_PER_SEC; 

  t = clock(); 
  // old.BiCGSTAB(); 
  t = clock() - t; 
  float solve0 = float(t)/CLOCKS_PER_SEC; 

  //cout << *old.x << endl; 
  cout << old.error->abs() << endl; 

  triLinSys sp(N); VecX<double> x1(N); VecX<double> err1(N); 
  sp.x = &x1; 
  sp.error = &err1; 

  t = clock(); 
  for (auto i =0; i < N; ++i) {
    sp.vt += {(double)i, (double)i, 4}; //old.A[i][i] = 2; 
    if (i > 0) sp.vt += {(double)i, (double)i-1, -1}; //old.A[i][i-1] = -1; 
    if (i < N-1) sp.vt += {(double)i, (double)i+1, -1}; //old.A[i][i+1] = -1; 
    if (i-N >= 0) sp.vt += {(double)i, (double)i-N, -1}; //old.A[i][i-1] = -1; 
    if (i+N < N) sp.vt += {(double)i, (double)i+N, -1}; //old.A[i][i+1] = -1; 

    sp.b[i] = 3*i+1; 
  }
  sp.setMat(sp.vt); 
  sp.b = sp.A*sp.b; 
  t = clock()-t; 
  float setup1 = float(t)/CLOCKS_PER_SEC; 

  t = clock(); 
  //sp.BiCGSTAB(); 
  t = clock()-t; 
  float solve1 = float(t)/CLOCKS_PER_SEC; 
  cout << " ----------- " << endl ; 
  //cout << *sp.x << endl; 
  cout << sp.error->abs() << endl; 

  cout << "OLD: setup-> " << setup0 << " sec, solve-> " << solve0 << "sec."<<endl; 
  cout << "NEW: setup-> " << setup1 << " sec, solve-> " << solve1 << "sec."<<endl; 

  // cout << "Abs" << endl; 
  // vector<Triplet> vt, bt; 
  
  // triLinSys Ab; 
  // vt += {{0,0,4}, {0,1,-1}, {1,2,-1}, {2,1,-1}}; 
  // vt += {{2,2,4}, {3,2,-1}, {3,3,4}, {3,4,-1}}; 
  // bt += {{4,4,4}, {4,3,-1}, {1,1,4}, {2,3,-1}};
  // bt += {{1,0,-1}, {2,1,-2}, {1,4,-6}, {2,2,1}};

  //  Ab.setMat(vt+bt); 
  
  // for (auto i = 0; i < Ab.A.aij.size(); ++i) { 
  //   cout << i << " " << Ab.A.aij[i] << " " << Ab.A.val[i] << endl; 
  // }
  
  // VecX<double> b(5); 
  // cout << "----"<< endl<< Ab.A.rank << " " << b.size() << endl; 
  // for (auto i = 0; i < Ab.A.rank; ++i) b[i] = i+1; 
  
  // Ab.b = Ab.A*b;
  // old.setVECX(Ab.b); 
  
  // VecX<double> x(b.size());
  // Ab.error = &b; 
  // Ab.x = &x;
  // Ab.BiCGSTAB(); 
  
  // cout << *Ab.x << endl; 
  // cout << Ab.error->abs() << endl; 

  // cout << old.A<<endl; 
  // cout << old.b<<endl; 
  // old.error = &b; 
  // old.x = &x;

  // old.BiCGSTAB(); 
  
  // cout << *old.x << endl; 
  // cout << old.error->abs() << endl; 
  
  

  return 0; 
};
