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
#include <girdap>
#include <linSolve/MatY.hpp>
#include <time.h>

//------------------------------------------------------------------//
// ---------------- LAPLACE -- CELL CENTER -------------------------//
//------------------------------------------------------------------//
// (1) The following function is the core others use this function  //
//------------------------------------------------------------------//
void Grid::triLaplace(triLinSys &axb, shared_ptr<Cell > f, double const &c) {
  if (!thisVar) {cout << "Laplace0: Variable is not locked!!"<< endl; return; }
  //  cout << "LAPLACE !!!! "<< endl; 
  //  cout << f->vol() << endl; 
  Scheme<double> sch = f->normGrad(thisVar); 
  auto area = f->vol().abs(); 
  int n = f->next; 
  int p = f->prev;
  for (auto i = 0; i < sch.size(); ++i) { 
    //    cout << "[" << n << ", "<< p << "] "<< c << " : " << sch.val[i] << " : " << area; 
    auto flux = c*sch.val[i]*area;
    //cout << " : " << flux << endl; 
    //cin.ignore().get(); 
    if (n >= 0) axb.A(n, sch.ind[i]) -= flux/listCell[n]->vol().abs(); 
    if (p >= 0) axb.A(p, sch.ind[i]) += flux/listCell[p]->vol().abs(); 
  }
  if (n >= 0) axb.b[n] += c*sch.c*area/listCell[n]->vol().abs(); 
  if (p >= 0) axb.b[p] -= c*sch.c*area/listCell[p]->vol().abs(); 
  return;  
  // int n = f->next; 
  // int p = f->prev;
  // auto c0 = (p>=0) ? *(listCell.begin() + p) : f;
  // auto c1 = (n>=0) ? *(listCell.begin() + n) : f;
  // Vec3 dx = c1->getCoord() - c0->getCoord(); 
  // Vec3 faceArea = f->vol(); 
  // Vec3 norm = faceArea/faceArea.abs();
  // if (p>=0 && n>=0) {
  //   double flux = c*faceArea.abs()/(dx*norm);
  //   axb.A[p][p] -= flux; 
  //   axb.A[n][n] -= flux; 
  //   axb.A[p][n] += flux; 
  //   axb.A[n][p] += flux;
  //   return; 
  // } else {
  //   if (!thisVar) {cout<<"laplace0: Boundary is not locked!"<<endl; return;}
  //   auto bndr = (p >= 0) ? -n-1 : -p-1; 
  //   auto row = (p >= 0) ? p : n; 
  //   if (bndr < 0 || bndr >= thisVar->listBC.size() || !(thisVar->listBC[bndr])) {
  //     cout << "wrong with bndr : Value "<< bndr <<endl; 
  //     return; 
  //   }
  //   auto bc = thisVar->listBC[bndr];
  //   if (bc->type == 0) {
  //     double flux = c*faceArea.abs()/(dx*norm);
  //     axb.A[row][row] -= flux*(1.0 - bc->a_val);
  //     axb.b[row] -= flux*(bc->b_val); 
  //   } else if (bc->type == 1) { 
  //     double flux = c*faceArea.abs();
  //     if (n >= 0) flux = -flux; 
  //     axb.A[row][row] += flux*(bc->a_val);
  //     axb.b[row] -= flux*(bc->b_val);
  //   } else {
  //     cout << "Type is not recognized!"<< endl; 
  //   }
  // }	
  return; 
}


//------------------------------------------------------------------//
// (2) Laplace with a constant Gamma throughout                     //
//------------------------------------------------------------------//
triLinSys Grid::laplace2(double c, initializer_list<double> n) {
  triLinSys axb;
  auto t = clock(); auto icnt = 0; 
  for (auto f : listFace) {
    triLaplace(axb, f, c); 
  }
  //  timeScheme(axb, n, thisVar->triLaplace);  
  t = clock()  -t; 
  //  cout << "Laplace is prepared in " << t/ (double) CLOCKS_PER_SEC << " secs"<< endl; 
  return axb;   
}; 

//------------------------------------------------------------------//
// (3) Laplace with a field variable Gamma (multiphase)             //
//------------------------------------------------------------------//
triLinSys Grid::laplace2(VecX<double>& c, initializer_list<double> n) {
  triLinSys axb;
  auto t = clock(); auto icnt = 0; 
  for (auto f : listFace) {
    auto n = f->next; 
    auto p = f->prev; 
    double cbar; 
    if (n >= 0 && p >= 0) cbar = 0.5*(c[n] + c[p]); 
    else if (n >= 0) cbar = c[n]; 
    else if (p >= 0) cbar = c[p]; 
    else { 
      cout << "Error: no cell is attached to a face ("<<f->id<<","; 
      cout<<"). next: "<<n<<" prev: "<<p<<endl; 
      exit(1); 
    }
    triLaplace(axb, f, cbar); 
  } 
  //timeScheme(axb, n, thisVar->triLaplace);  
  t = clock()  -t; 
  //cout << "Laplace is prepared in " << t/ (double) CLOCKS_PER_SEC << " secs"<< endl; 
  return axb; 
} 


//------------------------------------------------------------------//
// ---------------- SOURCE  -- CELL CENTER -------------------------//
//------------------------------------------------------------------//
// (1) The following function is the core others use this function  //
//------------------------------------------------------------------//
triLinSys Grid::source2(double c, double a, initializer_list<double> n) {
  triLinSys axb;
  auto t = clock(); 
  for (auto i = 0; i < listCell.size(); ++i) {
    // double vol = listCell[i]->vol().abs(); 
    if (c != 0) {
      axb.A(i,i) += c;//*vol; 
    }
    axb.b[i] -= a;//*vol;
  }
  t = clock() - t; 
  //cout << "Source is prepared in " << t/ (double) CLOCKS_PER_SEC << " secs"<< endl; 

  return axb; 
}

//------------------------------------------------------------------//
// (2) The following function is the core others use this function  //
//------------------------------------------------------------------//
triLinSys Grid::source2(double c, VecX<double>& a, initializer_list<double> n) {
  triLinSys axb;
  //auto t = clock(); 
  for (auto i = 0; i < listCell.size(); ++i) {
    //double vol = listCell[i]->vol().abs(); 
    axb.A(i,i) += c;//*vol; 
    axb.b[i] -= a[i];//*vol;
  }
  //t = clock() - t; 
  //cout << "Source is prepared in " << t/ (double) CLOCKS_PER_SEC << " secs"<< endl; 

  return axb; 
}


//------------------------------------------------------------------//
// ---------------- DIVERGENCE -- CELL CENTER - FOU-----------------//
//------------------------------------------------------------------//
// (1) The following function is the core others use this function  //
//------------------------------------------------------------------//
void Grid::triDiv(triLinSys &axb, shared_ptr<Cell > f, Vec3 const &c) {
  int n = f->next; 
  int p = f->prev;

  Vec3 norm = f->vol(); 
  double flux = c*norm;
  auto xf = f->getCoord(); 
  auto xphi = xf - c*dt; 
  norm /= norm.abs(); 
  if (p>=0 && n>=0) {
    auto phif = f->phi(thisVar).eval(thisVar); 
    auto row = (flux > 0) ? p : n; 
    auto invvoln = 1.0/listCell[n]->vol().abs(); 
    auto invvolp = 1.0/listCell[p]->vol().abs(); 
    // auto down = (row == p) ? n : p; 

    // auto t = xf - xphi; //t = t/t.abs(); 
    // auto d0 = thisVar->grad[row]*t; //d0 /= t.abs(); 
    // auto d1 = thisVar->grad[down]*t;
    // //    d1 /= t.abs(); 
    // //limiter 
    // if (d0*d1 < 0) {
    //   if (abs(d0) < abs(d1)) d0 = 0; 
    //   else d1 = 0; 
    // }
    auto dx0 =(xphi - listCell[row]->getCoord()); 
    double corr = 0.5*(phif + dx0*thisVar->grad[row]); // + (d0-d1)/12.0; 

    //    auto phix = thisVar->data[row] + corr; 

   // if (thisVar->data[n] > 0 && thisVar->data[n] < 1) {
   //   cout << "n: " << n << " p: "<< p << " Orient: " << f->orient; 
   //   cout << " vel: " << c << " xo: " << xphi; 
   //   cout << " xu: "<< listCell[row]->getCoord(); 
   //   cout << " del: "<< (xphi - listCell[row]->getCoord()); 
   //   cout << " xf : " << f->getCoord(); 
   //   cout << " xx : " << - 0.5*c*dt << endl; 
   //   cin.ignore().get(); 
   // }

    // if (phix < phiL) {
    //   corr = phiL - thisVar->data[row]; 
    // } else if (phix > phiU) {
    //   corr = phiU - thisVar->data[row]; 
    // }
    // if (n == row) {
    //   axb.A[p][p] += 0.25*flux; 
    //   axb.A[n][p] -= 0.25*flux; 
    // } else {
    //   axb.A[p][n] += 0.25*flux; 
    //   axb.A[n][n] -= 0.25*flux; 
    // }
    // if (p == row) axb.b[n] += 0.01*flux*thisVar->data[n]; 

    axb.A(p, row) += 0.5*flux*invvolp;
    axb.A(n, row) -= 0.5*flux*invvoln; 

    axb.b[p] -= flux*corr*invvolp; 
    axb.b[n] += flux*corr*invvoln; 

    return; 
  } else {
    if (!thisVar) {cout<<"Boundary is not locked!"<<endl; return;}
    double a, b; 
    auto row = (p >= 0) ? p : n; 
    auto sign = (p >= 0) ? 1.0 : -1.0; 
    if (flux*sign > 0) { //FOU
      axb.A(row, row) += flux*sign/listCell[row]->vol().abs();  
    } else {
      auto bndr = (p >= 0) ? -n-1 : -p-1; 
      auto invvol = (p >= 0) ? 1.0/listCell[p]->vol().abs()
	: 1.0/listCell[n]->vol().abs();  

      if (thisVar->listBC[bndr]->type == 0) 
	axb.b[row] -= flux*sign*thisVar->listBC[bndr]->b_val*invvol; 
    //   Scheme<double> s = f->phi(thisVar); 
    //   //      f->getBCPhi(thisVar, a, b); //
    //   //getBCPhi(f, thisVar, a, b);
    //   for (auto i=0; i < s.ind.size(); ++i) {
    // 	axb.A[row][s.ind[i]] += sign*s.val[i]; 
    //   }
    //   axb.b[row] -= sign*s.c; 
    //   axb.A[row][row] += sign*flux*a; 
    //   axb.b[row] -= sign*flux*b; 
    }
  }	
  return; 
}

//----------------------------------------------------------------------//
// (2) The following uses the core function to evaluate the whole array //
//----------------------------------------------------------------------//
triLinSys Grid::div2(VecX<Vec3>& vel, double c, initializer_list<double> n) {  
  triLinSys axb; 

  auto U = getVar("u"); 
  auto V = getVar("v"); 
  auto W = getVar("w"); 
  //cout << axb.A << endl; 
  auto t = clock(); auto icnt = 0; 
  //double flux = 1; 
  for (auto f : listFace) {
    auto n = f->next; 
    auto p = f->prev; 
    Vec3 vbar; 
    double ubc, vbc, wbc; 
    ubc = f->phi(U->listBC).eval(vel.comp(0)); //
    //getPhi_val(f, vel.comp(0), ubc); 
    vbc = f->phi(V->listBC).eval(vel.comp(1)); //
    //getPhi_val(f, V, vbc); 
    wbc = f->phi(W->listBC).eval(vel.comp(2)); //
    //getPhi_val(f, W, wbc); 
    vbar = {ubc, vbc, wbc}; 
    triDiv(axb, f, c*vbar); 
  }
  //timeScheme(axb, n, thisVar->div);  
  t = clock()  -t; 
  //cout << "Div is prepared in " << t/ (double) CLOCKS_PER_SEC << " secs"<< endl; 

  return axb;   
}; 

//----------------------------------------------------------------------//
// (3) The following uses the core function to evaluate the whole array //
//----------------------------------------------------------------------//
triLinSys Grid::div2(VecX<Vec3>& vel, VecX<double> c, initializer_list<double> n) {
  triLinSys axb;
  auto U = getVar("u"); 
  auto V = getVar("v"); 
  auto W = getVar("w"); 

  //cout << axb.A << endl; 
  auto t = clock(); auto icnt = 0; 
  //double flux = 1; 
  for (auto f : listFace) {
    auto n = f->next; 
    auto p = f->prev; 
    double cbar; Vec3 vbar; 
    if (n >= 0 && p >= 0) cbar = 0.5*(c[n] + c[p]); 
    else if (n >= 0) cbar = c[n]; 
    else if (p >= 0) cbar = c[p]; 
    double ubc, vbc, wbc; 
    ubc = f->phi(U->listBC).eval(vel.comp(0)); //
    //getPhi_val(f, vel.comp(0), ubc); 
    vbc = f->phi(V->listBC).eval(vel.comp(1)); //
    //getPhi_val(f, V, vbc); 
    wbc = f->phi(W->listBC).eval(vel.comp(2)); //
    // getPhi_val(f, U, ubc); 
    // getPhi_val(f, V, vbc); 
    // getPhi_val(f, W, wbc); 
    vbar = {ubc, vbc, wbc}; 
    triDiv(axb, f, cbar*vbar); 
  }
  //  timeScheme(axb, n, thisVar->div);  
  t = clock()  -t; 
  //cout << "Div is prepared in " << t/ (double) CLOCKS_PER_SEC << " secs"<< endl; 

  return axb;   
}; 

// triLinSys Grid::divRK2E(VecX<Vec3> &vel, double c) {
//   //back up data; 
//   auto phi = thisVar->data; 
//   auto grad = thisVar->grad; 
//   auto k1 = triDiv(vel, c, {0});
//   thisVar->data = phi + 0.5*k1.b*dt;
//   thisVar->grad = valGrad(thisVar); 
//   auto k2 = triDiv(vel, c, {0}); 
//   thisVar->data = phi; 
//   thisVar->grad = grad; 
//   return k2; 
// }

// triLinSys Grid::divRK4E(VecX<Vec3> &vel, double c) {
//   //back up data; 
//   auto phi = thisVar->data; 
//   auto grad = thisVar->grad; 
//   auto k1 = triDiv(vel, c, {0});
//   thisVar->data = phi + 0.5*k1.b*dt;
//   thisVar->grad = valGrad(thisVar); 
//   auto k2 = triDiv(vel, c, {0});
//   thisVar->data = phi + 0.5*k2.b*dt;
//   thisVar->grad = valGrad(thisVar); 
//   auto k3 = triDiv(vel, c, {0});
//   thisVar->data = phi + k3.b*dt;
//   thisVar->grad = valGrad(thisVar); 
//   auto k4 = triDiv(vel, c, {0});
//   k4.b = 1.0/6.0*(k1.b + 2.0*k2.b + 2.0*k3.b + k4.b); 
//   thisVar->data = phi; 
//   thisVar->grad = grad; 
//   return k4; 
// }



// void Grid::getBCPhi(shared_ptr<Cell > f, shared_ptr<Var> phi, double &a, double &b) {
//   int n = f->next; 
//   int p = f->prev;
//   auto bndr = (p >= 0) ? -n-1 : -p-1; 
//   auto row = (p >= 0) ? p  : n; 
//   if (bndr < 0 || bndr >= phi->listBC.size() || !(phi->listBC[bndr])) {
//     cout << "Something wrong with bndr : Value "<< bndr <<endl; 
//     return; 
//   }
//   if (phi->listBC[bndr]->type == 0) {
//     a = phi->listBC[bndr]->a_val; 
//     b = phi->listBC[bndr]->b_val; 
//   } else if (phi->listBC[bndr]->type == 1) { 
//     auto c0 = (p>=0) ? *(listCell.begin() + p) : *(listCell.begin()+n);	
//     auto norm = f->vol(); norm = norm/norm.abs(); 
//     double dx = norm*(f->getCoord() - c0->getCoord()); 
//     a = (1.0 + dx*phi->listBC[bndr]->a_val);
//     b = dx*phi->listBC[bndr]->b_val; 
//   }
//   return;
// }


// void Grid::getPhi_val(shared_ptr<Cell > f, shared_ptr<Var> phi, double &out) {
//   int n = f->next; 
//   int p = f->prev;
//   if (n >= 0 && p >= 0) {
//     out = 0.5*(phi->data[p] + phi->data[n]); 
//   } else {
//     double a=0, b=0; 
//     f->getBCPhi(phi,a,b); //
//     //getBCPhi(f, phi, a, b); 
//     out = (p >= 0) ? a*phi->data[p] + b : a*phi->data[n] + b; 
//   } 
//   return; 
// }

// VecX<double> Grid::valDiv(VecX<Vec3> vel) { 
//   VecX<double> a(vel.rank); 

//   auto U = getVar("u"); 
//   auto V = getVar("v"); 
//   auto W = getVar("w"); 

//   for (auto f : listFace) {
//     auto n = f->next; 
//     auto p = f->prev; 
//     Vec3 vbar; 
//     double ubc, vbc, wbc; 
//     ubc = f->phi(U->listBC).eval(vel.comp(0)); //
//     vbc = f->phi(V->listBC).eval(vel.comp(1)); //
//     wbc = f->phi(W->listBC).eval(vel.comp(2)); //
//     // getPhi_val(f, U, ubc); 
//     // getPhi_val(f, V, vbc); 
//     // getPhi_val(f, W, wbc); 
//     vbar = {ubc, vbc, wbc}; 
//     double flux = vbar*f->vol(); 
//     if (p>=0) a[p] += flux; 
//     if (n>=0) a[n] -= flux; 
//   }
//   for (auto c : listCell) {
//     a[c->id] = a[c->id]/c->vol().abs(); 
//   }
//   return a;   
// }

// VecX<Vec3> Grid::valGrad(shared_ptr<Var > phi) { 
//   VecX<Vec3> a(listCell.size(), {0,0,0}); 
//   for (auto f : listFace) {
//     double phif = f->phi(phi->listBC).eval(phi); 
//     //getPhi_val(f, phi, phif); 
//     Vec3 flux = (f->vol())*phif; 
//     auto n = f->next; 
//     auto p = f->prev;
//     if (n >= 0)  a[n] -= flux; 
//     if (p >= 0)  a[p] += flux; 
//   }
//   for (auto c : listCell) {
//     a[c->id] = a[c->id]/c->vol().abs(); 
//   }
//   return a;   
// }


// void Grid::correctVel(double coef) {
//   auto U = getVar("u");   auto V = getVar("v");   auto W = getVar("w"); auto P = getVar("p"); 
//   VecX<Vec3> sum(listCell.size()); 
//   VecX<Vec3> factor(listCell.size());
  
//   for (auto f : listFace) {
//     auto n = f->next, p = f->prev; 
//     auto norm = f->vol(); 
//     auto area = norm.abs(); norm /= area; 
//     Vec3 gpf, usf, uf; 
//     double Uf = f->phi(U).eval(U); //(p >=0) ? a*vel[p].x() + b : a*vel[n].x() + b; 
//     double Vf = f->phi(V).eval(V); //(p >=0) ? a*vel[p].y() + b : a*vel[n].y() + b; 
//     double Wf = f->phi(W).eval(W); //(p >=0) ? a*vel[p].z() + b : a*vel[n].z() + b; 
//     usf = Vec3(Uf, Vf, Wf); 
//     gpf = f->grad(P).eval(P); 
//     uf = (usf*norm)*norm - coef*(gpf*norm);
//     if (n >= 0) { sum[n] += area*uf; factor[n] += area*norm;}
//     if (p >= 0) { sum[p] += area*uf; factor[p] += area*norm;}
//   }
//   for (auto i = 0; i < listCell.size(); ++i) {
//     U->data[i] = sum[i].x()/factor[i].x(); 
//     V->data[i] = sum[i].y()/factor[i].y(); 
//     if (factor[i].z() != 0) W->data[i] = sum[i].z()/factor[i].z(); 
//   }
//   return; 
// }





// VecX<double> Grid::valDiv(VecX<Vec3> c) {
//   VecX<double> a; 
//   for (auto v : listVertex) {
//     for (auto i = 0; i< 2; ++i) {
//       if (v->face[i]) {
// 	auto n = v->face[i]->next; 
// 	auto p = v->face[i]->prev; 
// 	Vec3 cbar; 
// 	if (n >= 0 && p >= 0) cbar = 0.5*(c[n] + c[p]); 
// 	else {

	  

// 	  if (n >= 0) {
// 	  // condition for velocity 	  
// 	  cbar = c[n]; 
// 	} else if (p >= 0) {
// 	  // condition for velocity
// 	  cbar = c[p]; 
// 	}
// 	else { 
// 	  cout << "Error: no cell is attached to a face ("<<v->id<<","<<i; 
// 	  cout<<"). next: "<<n<<" prev: "<<p<<endl; 
// 	  exit(1); 
// 	}
// 	div(axb, v->face[i], cbar); 
//       }
//     }
//   }

// }



// void Grid::limiter(LinSys &axb, shared_ptr<Vertex > v, shared_ptr<Cell > f, Vec3 const &c) {
//   Vec3 flux = f->vol()*c; 
//   double r = 0; 
//   // find ff 
//   if (f->next>=0 && f->prev>=0) {
//     if (flux > 0) {
//       auto D = *(listCell.begin() + f->next);
//       auto U = *(listCell.begin() + f->prev);
//       auto UU = U; 
//       if (f->node[0]->id == f->node[0]) {//north south
//         UU = (*(listCell.begin() + f->prev))->node[0]->cell[3];
// 	if (UU >= 0) { r = (thisVar->data[U]-thisVar->data[UU])/
//       } else { // east west
// 	UU = (*(listCell.begin() + f->prev))->node[0]->cell[1]; 
//       }
//       if (UU->id == U->id) {
// 	cout << "UU not found" << endl; 
//       }
// }

// //------------------------------------------------------------------//
// // ---------------- DDT -- CELL CENTER -----------------------------//
// //------------------------------------------------------------------//
// // (1) The following function is the core others use this function  //
// //------------------------------------------------------------------//
// LinSys Grid::ddt(double c) {
//   LinSys axb(listCell.size());
//   if (thisVar) thisVar->dt = dt; 
//   auto t = clock(); 
//   for (auto i = 0; i < listCell.size(); ++i) {
//     //    double vol = listCell[i]->vol().abs(); 
//     if (c != 0) {
//       axb.A(i, i) += c/dt; 
//     }
//     axb.b[i] += c/dt*(thisVar->data[i]);
//   }
//   t = clock() - t; 
//   //cout << "DDT (Euler) is prepared in " << t/ (double) CLOCKS_PER_SEC << " secs"<< endl; 
//   return axb; 
// }

// LinSys Grid::ddt(VecX<double> c) {
//   LinSys axb(listCell.size());
//   if (thisVar) thisVar->dt = dt; 
//   auto t = clock(); 
//   for (auto i = 0; i < listCell.size(); ++i) {
//     //    double vol = listCell[i]->vol().abs(); 
//     if (c[i] != 0) {
//       axb.A(i, i) += c[i]/dt; 
//     }
//     axb.b[i] += c[i]/dt*(thisVar->data[i]);
//   }
//   t = clock() - t; 
//   //cout << "DDT (Euler) is prepared in " << t/ (double) CLOCKS_PER_SEC << " secs"<< endl; 
//   return axb; 
// }


// void Grid::timeScheme(LinSys &axb, initializer_list<double> &n, VecX<double> &prev) {  
//   if (n.size()==0) {
//   } else if (n.size() == 1) { 
//     double n0 = *(n.begin()); 
//     axb.b = (axb.b - (1.0-n0)*(axb.A*(thisVar->data))); 
//     axb.A = n0*axb.A;  
//   } else if (n.size() == 2) { 
//     double n0 = *(n.begin()); 
//     //cout << prev.size() <<  " " << dt0 << endl; 
//     if (prev.size() == 0 || dt0 == 0) {
//       prev = axb.b - axb.A*(thisVar->data); 
//       //cout << " ++--> "<< prev << endl; 
//       axb.b = (axb.b - (1.0-n0)*(axb.A*(thisVar->data))); 
//       axb.A = n0*axb.A; 
//     } else {
//       if (dt != dt0 && n0 != 0.0) { cout << "variable time step is only allowed with Adam-Bashfort scheme"<< endl; exit(1); }
//       double n1 = *(n.begin()+1);
//       double n2 = (1.0 - n1 - n0);
//       //cout << "--> " << n0 << " "<< n1 << " "<< n2 << endl; 
//       auto prevNew = axb.b - axb.A*(thisVar->data); 
//       //cout << " ::==> " << prevNew << endl; 
//       //cout << " ::--> " << prev << endl; 
//       if (n1 > 1.0) {
// 	axb.b = n0*(axb.b) + (1.0+(n1-1.0)*dt/dt0)*prevNew + n2*dt/dt0*(prev); 
//       } else {
// 	axb.b = n0*(axb.b) + (n1*dt/dt0)*prevNew + n2*dt/dt0*(prev); 
//       }
//       //      cout << " ::~~> " << axb.b << endl; 
//       axb.A = n0*(axb.A); 
//       prev = prevNew; 
//     }
//   }
//   return; 
// }

