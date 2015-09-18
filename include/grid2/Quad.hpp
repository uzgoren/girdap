#ifndef QUAD
#define QUAD

class Quad: public Cell {
protected: 
  const int_2 map[4] = {2, 3, 0, 1};  
public:
  //  short int type;   
  Quad(initializer_list<int_8> l):Cell(l) {
    setType(9); 
  }
  
  void assignCelltoNode(); 
  void convertToSimpleBlock(initializer_list<int> n, bool debug=false);
  bool refine();
  void refine(int dir);
  bool coarsen(int dir); 
  // // initNodeList assigns this cell to all vertices (hanging nodes hang)
  // void initNodeList(shared_ptr<Cell > c) {
  // }
  // // assignHangingNode updates the 
  int_8 hangingVertexOnFace(int_2 i, int_2 j=0) { 
    shared_ptr<Vertex > *v0, *v1;
    int i0=(i+1)%4; int i1=(i0+1)%4;
    v0 = getVertex(i); v1 = getVertex(i0); 
    if (v0 && v1) {
      if ((*v0)->cell[i0] >= 0 && ((*v0)->cell[i0] != (*v1)->cell[i])) {
        return (*((*v0)->getCell(i0)))->node[i1]; 
      }
    }  
    return -1;
  }
  Vec3 edge(int_2 i=0) {
    int_2 j, k; 
    if      (i == 0) {j = 1; k = 0; }
    else if (i == 1) {j = 2; k = 1; }
    else if (i == 2) {j = 2; k = 3; }
    else if (i == 3) {j = 3; k = 0; }
    auto v1 = getVertex(j); 
    auto v0 = getVertex(k); 
    if ((*v1) && (*v0)) return **v1 - **v0;
    cout << "Error: Edge do not exist for Quad: "<< i <<endl; 
    exit(1); 
  }
  Vec3 vol() { 
    return 0.5*((**getVertex(2)) - (**getVertex(0)))
      ^((**getVertex(3)) - (**getVertex(1)));
  }
  double dx() {
    return 0.5*(edge(0).abs() + edge(2).abs()); 
  }
  double dy() {
    return 0.5*(edge(1).abs() + edge(3).abs()); 
  }
  double dz() { return 0; }
  // Scheme phi(int_2 bias=0); 
  vector<int_8> ngbrCellList(); 
  void checkIslandLevels(int dir);
  void checkNgbrLevel(int dir);
  // // shared_ptr<Cell >  makeFace(int_2 i=0) { 
  // //   vector<shared_ptr<Vertex > > tmp;
  // //   if (i >= 4) {cout << "Quad: makeFace(i): i is out of bounds "<<i<<endl; exit(1);}
  // //   int_2 j = i+1; if (j == 4) j = 0;
  // //   tmp.emplace_back(node.at(i)); 
  // //   tmp.emplace_back(node.at(j)); 
  // //   shared_ptr<Cell > c = shared_ptr<Cell >(new Line(tmp));
  // //   return c;
  // //}
};




#endif
