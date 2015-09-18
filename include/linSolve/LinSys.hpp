#ifndef LINSYS
#define LINSYS

#include "MatX.hpp"

class LinSys {
private:
  int_8 N; 
  double tol;
  int itmax; 

public:
  MatX<double> A; 
  VecX<double> b; 
  VecX<double> *x; 
  VecX<double> *error; 

  double err; 
  int it; 

  //-------------------------
  // Constructor
  //-------------------------
  LinSys() { setLimits();  }

  LinSys(int_8 const &a) {
    setLimits();
    N = a; 
    resize(a); 
  }
  

  //-------------------------
  // Methods
  //-------------------------
  void resize(int_8 const &a) {
    N = a; 
    A.resize(a); 
    b.resize(a); 
  }  
  void clear() {
    A.empty(); 
    b.empty(); 
  }
  void setLimits(double t=1e-6, unsigned int imax = 1000) {
    if (t>0) tol = t; 
    if (imax>0) itmax = imax; 
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
    int iter = 0; 
    double alpha, rsold, rsnew; 
    VecX<double> p(A.data.size()), q(A.data.size());
 
    *error = (b - A*(*x));
    p = *error; 
    rsold = *error*(*error);
    for (; iter < itmax; ++iter) {
	q = A*p;
	alpha = rsold/(q*p);
	*x += alpha*p;
	*error -= alpha*q;
	rsnew = (*error)*(*error); 
	if (sqrt(rsnew) < tol) break;
	p = *error + rsnew/rsold*p; 
	rsold = rsnew; 
	//if (iter%50 == 0) 
	  // cout << "CG: " << iter << ": "<< log10(rsnew) << endl; 
      }  
    it = iter; 
    err = sqrt(rsnew); 
    //cout << "CG: " << iter << ": " << log10(rsnew) << " DONE!" << endl;
   }
    
  
  void Gauss_Seidel() {
    int iter;
    double sigma; 
    
    iter = 0; 
    while (iter<itmax){
      int_8 r = 0; 
      err = 0;
      for (MatX<double>::rowiter rit=A.data.begin(); rit!=A.data.end(); ++rit, ++r) {
	sigma = 0; 
	for (VecX<double>::it_cmap cit=rit->px.begin(); cit!=rit->px.end(); ++cit) {
	  if (cit->first == r) continue;
	  sigma += cit->second * (*x)[cit->first]; 
	}
	sigma = (b[r] - sigma)/rit->px[r];
	err += ((sigma-(*x)[r])); //err *= err; 
	(*x)[r] = sigma; 
      }
      err = sqrt(err/A.rows); 
      if (err < tol) break;
      //	cout << "GS: " << iter << "->" << log10(err) <<endl;
      iter+=1;
    }
    it = iter; 
    //      cout << "GS: " << iter << ", " << log10(err) << endl;
  }

  void BiCGSTAB() {
    int iter = 0; 
    double rho, alpha, beta, omega, rho_new; 
    VecX<double> r(A.data.size()), r0(A.data.size()), p(A.data.size());
    VecX<double> t(A.data.size()), v(A.data.size()), s(A.data.size());
    
    err = 0; 
    //check if t is initialized as zero; 
    //cout<< t << endl; 

    *error = (b - A*(*x));
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




  LinSys operator+(LinSys const &a) {
    LinSys c; 
    c.A = A + a.A;
    c.b = b + a.b; 
    return c;
  }

  LinSys operator-(const LinSys& a) {
    LinSys c; 
    c.A = A - a.A;
    c.b = b - a.b; 
    return c;
  }

  LinSys operator-() {
    A = -A; 
    b = -b; 
    return *this; 
  }
 
};


#endif
