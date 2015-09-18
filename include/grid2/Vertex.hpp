#ifndef VERTEX
#define VERTEX

#include <VecX>
#include <grid2/Scheme.hpp>

class Grid; 
class Cell;

class Vertex : public Vec3 {
protected:
    
public:
  vector<int_8 > cell; 
  // Cell contains structured array
  // many constains node shared by many cells; unstructured
  vector<int_8 > many;
  // When both cell and many has elements it is a mixed node;
  // Mixed nodes should be dealt specifically
  // Keep at most 3 faces (attached to a vertex) for structured
  vector<shared_ptr<Cell> > face; 
  //  vector<shared_ptr<Vector> > ngbr; 
  int_2 bndr; 

  int_8 id; 
  Grid* grid; 
  
  Vertex():Vec3() {face.resize(3);}
  Vertex(double const & a, double const & b, double const & c) :  Vec3 (a, b, c) {face.resize(3);}
  Vertex(initializer_list<double> a) : Vec3( a ) { face.resize(3); }
  Vertex(Vec3 a) : Vec3(a) {face.resize(3); }

  void reset() {cell.clear();}
  void reset(int_2 s) {cell.assign(s, -1);}
  void reset(initializer_list<int_8> l) {cell.assign(l.begin(), l.end());}

  void cellResize(int_2 size); 
  void setCell(int_2 ind, int_8 value);
  void setCell(initializer_list<int_8> l) {for (auto i = 0; i < l.size(); ++i) {if (cell[i] < 0) cell[i] = *(l.begin()+i); }}
  void replaceCell(int_2 ind, int_8 value);
  shared_ptr<Cell > *getCell(int_2 ind, bool debug=0); 
  shared_ptr<Vertex > *ngbr(int_2 d);
  // void addFace();
  
  Scheme<double> phi(vector<shared_ptr<Boundary> > const &bc, int_2 bias=0);  
  Scheme<double> phi(shared_ptr<Var> var, int_2 bias=0);  
    
  
  // shared_ptr<Cell > isHanging(int_2 dir) {
  //    if (cell.size() == 2) return NULL;     
  //   int_2 ind[4] = {dir, static_cast<int_2>((dir+1)%4)
  // 		    , static_cast<int_2>(dir+4)
  // 		    , static_cast<int_2>(((dir+1)%4)+4)};
  //   if (cell.size() == 4) {
  //     if (cell[ind[0]]) { // one combination
  // 	if (cell[ind[0]].get() == cell[ind[1]].get()) 
  // 	  return cell[ind[0]]; 
  //     } 
  //   } else if (cell.size() == 8) {
  //     if (cell[ind[0]]) { // many combinations 
  // 	if (cell[ind[0]].get() == cell[ind[1]].get() 
  // 	    || cell[ind[0]].get() == cell[ind[2]].get())
  // 	  return cell[ind[0]];
  //     } 
  //     if (cell[ind[3]]) {
  // 	if (cell[ind[3]].get() == cell[ind[1]].get()
  // 	    || cell[ind[3]].get() == cell[ind[2]].get())
  // 	  return cell[ind[3]]; 
  //     }
  //   }
  //   return NULL; 
  // }  
  //  void updateFaceLink(int i, int j, int k, int_8 c0, int_8 c1);  

};

#endif
