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

class Triplets {
public:
  vector<Triplet> data; 
 
  Triplets() {}; 
  
  Triplets& operator+=(const initializer_list<double> a) {
    data.emplace_back(Triplet((int_8)*(a.begin()), (int_8)*(a.begin()+1), 
			 *(a.begin()+2))); 
    return *this;
  } 

  Triplets& operator+=(const initializer_list<initializer_list<double> > a) {
    for (auto it = a.begin(); it < a.end(); ++it) 
      *this += *it; 
    return *this; 
  }

  Triplets& operator+=(Triplets &b) {
    this->data.reserve(this->data.size() + b.data.size()); 
    std::move(b.data.begin(), b.data.end(), std::back_inserter(this->data));
    b.data.clear();
    return *this;
  }

  Triplets& operator+(vector<Triplet> &b) {
    this->data.reserve(this->data.size() + b.size()); 
    std::move(b.begin(), b.end(), std::back_inserter(this->data));
    b.clear();
    return *this;
  }
  
  Triplets& operator+(Triplets &b) {
    *this = *this + b.data; 
    return *this; 
  }
  
  Triplets& operator*(double a) {
    for (auto c=data.begin(); c < data.end(); ++c)  
      c->val = a*c->val; 
    return *this;
  }
  
  friend Triplets& operator*(double a, Triplets& b) {
    return b*a; 
  }
  
//  friend Triplets& operator*(double a, vector<

  Triplets& operator/(double a) {
    return *this*(1/a); 
  }

  friend Triplets& operator/(double a, Triplets &b) {
   return b*(1.0/a); 
 }

  Triplets& operator-() {
    return *this*(-1); 
  }

  Triplets& operator-(Triplets &b) {
    return *this + (-b); 
  }

  double& operator()(int_8 i, int_8 j) {
    data.emplace_back(Triplet(i, j, 0));
    return (data.rbegin()->val); 
  }
};


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
  Triplets A; 
  Sparse S; 
  VecX<double> b;
  VecX<double> *x; 
  VecX<double> *error; 

  double err; 
  int it; 
  
  triLinSys() {}; // setLimits(); S.rank = 0;}; 
  triLinSys(int_8 const &N) { S.rank = N; b.resize(N); setLimits(); };
  

  void setLimits(double t=1e-6, unsigned int imax = 1000) {
    if (t>0) tol = t; 
    if (imax>0) itmax = imax; 
  }


  void setMat(Triplets &vtr) {
    if (vtr.data.empty()) {
      if (S.rank > 0) {
	S.aij.assign(S.rank+1, S.rank+1);
	S.val.assign(S.rank+1, 0); 
      }
      return; 
    }
    // first sort a (1st col then row); 
    sort(vtr.data.begin(), vtr.data.end(), 
	 [](const Triplet & a, const Triplet & b) -> bool
	 { 
	   return (a.row < b.row) || ((a.row == b.row) && (a.col < b.col)); 
	 });
    S.rank = max(S.rank, vtr.data.rbegin()->row + 1); 
    
    //cout << " ------------------ " <<endl;     
    // for (auto c : vtr) { 
    //   cout << c.row << " " << c.col << " " << c.val <<endl; 
    // }
    // cout << " ------------------ " <<endl; 

    S.aij.reserve(30*S.rank); 
    S.val.reserve(30*S.rank); 

    S.aij.assign(S.rank+1, S.rank+1);
    S.val.assign(S.rank+1, 0); 

    int_8 oldrow = vtr.data.begin()->row;

    for (auto c : vtr.data) {
      // c.row == c.col
      if (c.row == c.col) {
	S.val[c.row] += c.val;
	continue; 
      }
      if (c.row != oldrow) {
	S.aij[c.row+1] = S.aij[c.row]; 
	oldrow = c.row; 
      }

      // check if col exist at the last location
      if (*(S.aij.rbegin()) == c.col) {
	*(S.val.rbegin()) += c.val; 
	continue; 
      }

      S.aij.emplace_back(c.col); 
      S.val.emplace_back(c.val); 
      S.aij[c.row+1]++; 

    }
    vtr.data.clear(); 
    //cout << " Size: " << A.aij.size() << " " << A.val.size() <<endl; 


  } 


  void BiCGSTAB() {
    int iter = 0; 
    double rho, alpha, beta, omega, rho_new; 
    VecX<double> r(S.rank), r0(S.rank), p(S.rank);
    VecX<double> t(S.rank), v(S.rank), s(S.rank);

    err = 0; 
    //check if t is initialized as zero; 
    //cout<< t << endl; 

    *error = (b - S*(*x));
    //    cout << error->abs() << endl; 
    if (error->abs() < tol) { err = tol; it = 0; return; }
    r0 = *error; 
    rho = 1; alpha = 1; omega = 1; 
    for (; iter < itmax; ++iter) {
      rho_new = r0*(*error); 
      beta = (rho_new/rho)*(alpha/omega);
      rho = rho_new; 
      p = (*error) + beta*(p - omega*v); 
      v = S*p; 
      //cout << " r0v: " << r0*v << " p: " << p << endl; 
      alpha = rho/(r0*v);
      s = *error - alpha*v; 
      t = S*s; 
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

  triLinSys& operator+(triLinSys &T) {
    A += T.A; 
    return *this; 
  }

  triLinSys& operator-(triLinSys &T) {
    T.A = -T.A; 
    A += T.A; 
    return *this; 
  }



};




#endif
