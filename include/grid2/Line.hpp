#ifndef LINE
#define LINE

class Line: public Cell {
protected:
  const int_2 map[2] = {1, 0};  
public:
  Line(initializer_list<int_8> l):Cell(l) {
    setType(3); 
  }
  Vec3 vol() { 
    Vec3 a = edge(); 
    return Vec3(-a.data[1], a.data[0], 0); 
  }
  Vec3 edge(int_2 i=0) { return (**getVertex(1)) - (**getVertex(0));} 
  //  Vec3 grad(shared_ptr<Var > phi);
  //void getBCPhi(shared_ptr<Var> phi, double &a, double &b);  
};


#endif
