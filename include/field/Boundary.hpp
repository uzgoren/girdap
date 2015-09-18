#ifndef BOUNDARY
#define BOUNDARY


class Boundary {
public:
  int type; // 0: Dirichlet and 1: Neumann
  double a_val, b_val; // Linearized a_val*phi_p + b_val; 
 
  Boundary():type(1), b_val(0), a_val(0) {}
  Boundary(int t) :type(t), b_val(0), a_val(0) {}
  Boundary(int t, double b) : type(t), b_val(b), a_val(0) {}
  Boundary(int t, double b, double a) : type(t), b_val(b), a_val(a){}

  void setBC() {
    type = 1; a_val = 0; b_val = 0; 
  };

  void setBC(int t) {
    type = t; b_val = 0; a_val = 0; 
  }; 
 
  void setBC(int t, double b, double a=0) {
    type = t; b_val = b; a_val = a; 
  };
     
};


#endif
