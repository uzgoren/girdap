#ifndef FACE
#define FACE

class Cell; 

class Face {
public:
  typedef unsigned long int uint_8;  
  vector<Vertex*  > node;
  Cell* parent; 
  
  Face(Cell* c, unsigned short int k) {
    parent = c; 
  }

};

#endif
