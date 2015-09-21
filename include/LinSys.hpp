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
#ifndef LINSYS
#define LINSYS

#include <MatX>

class LinSys {
private:
  typedef unsigned long int uint_8; 
  uint_8 N; 
  double tol;
  unsigned int itmax; 

public:
  MatX<double> A; 
  VecX<double> b; 
  VecX<double> *x; 

  //-------------------------
  // Constructor
  //-------------------------
  LinSys() { setLimits();  }

  LinSys(uint_8 const &a) {
    setLimits();
    N = a; 
    resize(a); 
  }
  

  //-------------------------
  // Methods
  //-------------------------
  void resize(uint_8 const &a) {
    N = a; 
    A.resize(a); 
    b.resize(a); 
  }  
  void clear() {
    A.empty(); 
    b.empty(); 
  }
  void setLimits(double t=1e-6, unsigned int imax = 1000) {
    tol = t; itmax = imax; 
  }
  
 void setMATX(MatX<double>const &MAT){
  A=MAT;
 }

 void setVECX(VecX<double>const &VEC){
  b=VEC;
 }
 

  //------------------------
  // Operators
  //------------------------
  VecX<double> solve() {
    return b;
  }

  void CG() {
    unsigned int iter; 
    double alpha, rsold, rsnew; 
    VecX<double> r(N), p(N), q(N);
    
    cout << "askdkjasda" <<endl; 

    r = (b - A*(*x));
    cout << "----"<<endl<<A<<"-----"<<endl; 
    cout << x<<endl;
    cout << b<<endl; 
    cout << r<<endl; 
    p = r; 
    rsold = r*r;
    for (iter=0; iter < 10; ++iter) {
	q = A*p;
	alpha = rsold/(q*p);
	*x += alpha*p;
	r -= alpha*q;
	rsnew = r*r; 
	if (sqrt(rsnew) < tol) break;
	p = r + rsnew/rsold*p; 
	rsold = rsnew; 
      }  
    cout << "CG: " << iter << ", " << log10(rsnew) << endl;
   }
    


  
    void Gauss_Seidel() {
      unsigned int iter; 
      double sigma, err; 
    
      iter = 0; 
      while (iter<itmax){
	uint_8 r = 0; 
	err = 0;
	for (MatX<double>::rowiter rit=A.data.begin(); rit!=A.data.end(); ++rit, ++r) {
	  sigma = 0; 
	  for (VecX<double>::it_cmap cit=rit->px.begin(); cit!=rit->px.end(); ++cit) {
	    if (cit->first == r) continue;
	    sigma += cit->second * (*x)[cit->first]; 
	  }
	  sigma = (b[r] - sigma)/rit->px[r];
	  err += ((sigma-(*x)[r])); err *= err; 
	  (*x)[r] = sigma; 
	}
	err = sqrt(err/A.rows); 
	if (err < tol) break;
	//	cout << "GS: " << iter << "->" << log10(err) <<endl;
	iter+=1;
      }
      cout << "GS: " << iter << ", " << log10(err) << endl;
    }

 
};


#endif
