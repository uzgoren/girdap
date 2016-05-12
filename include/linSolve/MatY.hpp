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
#ifndef MATY
#define MATY

//-------------------------------------------
// Matrix storage scheme
//    Defines: basic matrix operations (not all)
//             Vec*Matrix addition
//-------------------------------------------

#include <VecX>

class Triplet {
public:
  int_8 row; 
  int_8 col; 
  double val; 

  Triplet(int_8 i, int_8 j, double a) { row = i; col = j; val = a;}

  // Operators
  // sum; 
};

vector<Triplet>& operator+=(vector<Triplet> &b, initializer_list<double> a) {
  b.emplace_back(Triplet((int_8)*(a.begin()), (int_8)*(a.begin()+1), 
			 *(a.begin()+2))); 
  return b; 
}

vector<Triplet>& operator+=(vector<Triplet> &b, 
 			   initializer_list<initializer_list<double> > a) {
  for (auto it = a.begin(); it < a.end(); ++it) 
    b += *it; 
  return b; 
}

vector<Triplet>& operator+(vector<Triplet> &a, vector<Triplet> &b) {
  a.reserve(a.size() + b.size()); 
  std::move(b.begin(), b.end(), std::back_inserter(a));
  b.clear();
  return a;
}

vector<Triplet>& operator*(double a, vector<Triplet> &b) {
  for (auto c=b.begin(); c < b.end(); ++c)  
    c->val = a*c->val; 
  return b;
}

vector<Triplet>& operator*(vector<Triplet> &b, double a) {
  b = a*b; 
  return b; 
}

vector<Triplet>& operator/(double a, vector<Triplet> &b) {
  b = (1.0/a)*b; 
  return b;
}

vector<Triplet>& operator/(vector<Triplet> &b, double a) {
  b = (1.0/a)*b; 
  return b;
}

vector<Triplet>& operator-(vector<Triplet> &a) {
  a = (-1.0)*a; 
  return a; //(-1.0)*a; 
}

vector<Triplet> operator-(vector<Triplet> &a, vector<Triplet> &b) {
   return a + (-b); 
}

class Sparse {
public:
  vector<int_8> aij; 
  vector<double> val; 
  int_8 rank; 

  VecX<double> operator*(VecX<double> &right) {    
    VecX<double> tmp(right.size()); 
    for (auto i = 0; i < tmp.size(); ++i) {      
      tmp[i] = val[i]*right[i];
      //      cout << i << " " << aij[i] << "-" << aij[i+1]-1 << endl; 
      for (auto j = aij[i]; j < aij[i+1]; ++j) { 
	auto col = aij[j]; 
	tmp[i] += val[j]*right[col];
      }
    }
    return tmp; 
  }

}; 


class triLinSys {
private:
  int_8 N; 
  double tol;
  int itmax; 

public:
  vector<Triplet> vt; 
  Sparse A; 
  VecX<double> b;
  VecX<double> *x; 
  VecX<double> *error; 

  double err; 
  int it; 
  
  triLinSys() { setLimits(); A.rank = 0;}; 
  triLinSys(int_8 const &N) { A.rank = N; b.resize(N); setLimits(); };
  

  void setLimits(double t=1e-6, unsigned int imax = 1000) {
    if (t>0) tol = t; 
    if (imax>0) itmax = imax; 
  }


  void setMat(vector<Triplet> &vtr) {
    if (vtr.empty()) {
      if (A.rank > 0) {
	A.aij.assign(A.rank+1, A.rank+1);
	A.val.assign(A.rank+1, 0); 
      }
      return; 
    }
    // first sort a (1st col then row); 
    sort(vtr.begin(), vtr.end(), 
	 [](const Triplet & a, const Triplet & b) -> bool
	 { 
	   return (a.row < b.row) || ((a.row == b.row) && (a.col < b.col)); 
	 });
    A.rank = max(A.rank, vtr.rbegin()->row + 1); 
    
    //cout << " ------------------ " <<endl;     
    // for (auto c : vtr) { 
    //   cout << c.row << " " << c.col << " " << c.val <<endl; 
    // }
    // cout << " ------------------ " <<endl; 

    A.aij.reserve(30*A.rank); 
    A.val.reserve(30*A.rank); 

    A.aij.assign(A.rank+1, A.rank+1);
    A.val.assign(A.rank+1, 0); 

    int_8 oldrow = vtr.begin()->row;

    for (auto c : vtr) {
      // c.row == c.col
      if (c.row == c.col) {
	A.val[c.row] += c.val;
	continue; 
      }
      if (c.row != oldrow) {
	A.aij[c.row+1] = A.aij[c.row]; 
	oldrow = c.row; 
      }

      // check if col exist at the last location
      if (*(A.aij.rbegin()) == c.col) {
	*(A.val.rbegin()) += c.val; 
	continue; 
      }

      A.aij.emplace_back(c.col); 
      A.val.emplace_back(c.val); 
      A.aij[c.row+1]++; 

    }
    vtr.clear(); 
    //cout << " Size: " << A.aij.size() << " " << A.val.size() <<endl; 


  } 


  void BiCGSTAB() {
    int iter = 0; 
    double rho, alpha, beta, omega, rho_new; 
    VecX<double> r(A.rank), r0(A.rank), p(A.rank);
    VecX<double> t(A.rank), v(A.rank), s(A.rank);

    err = 0; 
    //check if t is initialized as zero; 
    //cout<< t << endl; 

    *error = (b - A*(*x));
    //    cout << error->abs() << endl; 
    if (error->abs() < tol) { err = tol; it = 0; return; }
    r0 = *error; 
    rho = 1; alpha = 1; omega = 1; 
    for (; iter < itmax; ++iter) {
      rho_new = r0*(*error); 
      beta = (rho_new/rho)*(alpha/omega);
      rho = rho_new; 
      p = (*error) + beta*(p - omega*v); 
      v = A*p; 
      //cout << " r0v: " << r0*v << " p: " << p << endl; 
      alpha = rho/(r0*v);
      s = *error - alpha*v; 
      t = A*s; 
      //cout << " t: " << t << " s: " << s << endl; 
      omega = (t*s)/(t*t); 
      (*x) = (*x) + alpha*p + omega*s; 
      *error = s - omega*t; 
      err = sqrt(*error*(*error)); 
      if (err < tol) break; 
      //if (iter%50 == 0) 
	//cout << "BiCGstab: " << iter << ": "<< log10(err) << endl; 
    }  
    it = iter; 
    //err = sqrt(rsnew); 
    //cout << "BiCGstab: " << iter << ": " << log10(err) << " DONE!" << endl;
   }

};




#endif
