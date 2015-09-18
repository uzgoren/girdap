//#include "SparseMtrx.h"
//#include <Case>
#include <Fvm>


int main() {  
  // Read case files
  // 1. read domain (geometry, interface)
  // 2. read bcs
  // 3. create grid (based on bcs ?? )
  // 4. read init file
  // 5. read solver related conditions
  //Case case(); 
  //Grid grid(case); 
  //Eqn eqn(grid);   

  // Start Time iteration
  // Solve for Ax = b

  // End simulation

  Fvm vel(9); 
  Fvm p; 
  p.resize(9); 

  VecX<double> a(9); 
  MatX<double> x(9), y(5), z(9);
  a = {3, 2, 0, 3, 1, 2, 3, 4, 5};
  a += 1;
  cout << a << endl;
  VecX<double> b(9); 
  b = a*2; 
  cout << b <<endl; 
  b /= a; 
  cout << b << endl;
  a.info(); 
  z.insertDiag(a);
  cout << z << endl;
  a = -1; 
  z.insertDiag(a, -1); z.insertDiag(a, 1); 

  y.resize(9); 

  a = {1, 2, 3, 4, 5, 6, 7, 8, 9};

  cout << z << endl; 
  cout << a<< endl;
  cout << z*a <<endl;
  cout << z.getDiag() <<endl;


  return 0; 
};
