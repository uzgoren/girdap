
#include <linSolve/LinSys>
#include <grid2/Grid>
#include <iostream>
#include <fstream>

int main() {
  // Grid grid;
  // cout << "1" << &grid << endl; 
  // grid.addVertex({ {0, 0, 0}, {1, 0, 0} 
  //     , {1, 1, 0}, {0, 1, 0}, {0, 2, 0}, {1, 2, 0}}); 
  // cout << "2" << &grid << endl; 

  // grid.addCell({ { 3, 2, 5, 4} }); 
  // cout << "3" << &grid << endl; 
  // grid.addVertex({
  //     {2, 0, 0}, {2, 1, 0}, {2.1, 2.2, 0}, {0.5, 0, 0}
  //     , {0, 0.5, 0}, {1, 0.5, 0}
  //     , {0.5, 1, 0}, {0.5, 0.5, 0}
  //     , {1.5,1,0},{2,1.5,0},{1.5,2,0}
  //     ,{1,1.5,0}, {1.4, 1.6, 0} });
  // cout << "4" << &grid << endl; 
  // grid.addCell({ {2,14,18,17}, {14,7,15,18}, {1,6, 7, 2}
  //     , {18,15,8,16}, {17,18, 16,5}
  //     , {0, 9, 13, 10} , {9, 1, 11, 13}, {10, 13, 12, 3}, {13, 11, 2, 12} });
  // cout << "5" << &grid << endl; 
  //Grid* dummy = new Grid();
  //Grid* dummy2 = new Grid(); 
  // dummy = new shared_ptr(&grid); 
  //dummy2 = dummy; 
  // dummy2.reset(); 
  // cout << &dummy << " " << dummy->cbegin()->grid << endl;
  //cout << &grid << " " << grid.cbegin()->grid << endl; 
  //Grid dummy3 = grid; 
  //cout << &dummy3 << " " << dummy3.cbegin()->grid << endl;
  
  cout << "1" << endl; 
  Grid* block = new Block3({0,0,0}, {1,1,1}, 10, 9, 5); 

  if (block->listCell.capacity() < (block->listCell.size() + 30)) {block->listCell.reserve(max(block->listCell.capacity()*2,block->listCell.size() + 30));}
  (block->cbegin())->splitHexa(5,3,2); 
  //cout << "3" << endl; 
  //(block->cbegin())->splitHexa(2,2,2); 
  // cout << "4" << endl;
  
  // Grid* d2 = new Grid(); 
  // d2->addVertex({ {0,0,0}, {1,0,0}, {1,1,0}, {0, 1, 0}});
  // d2->addCell({0,1,2,3}); 
  // (d2->cbegin())->splitQuad(10,10); 
  
  ofstream myfile; 
  myfile.open("grid3d.vtk"); 
  myfile << *block << endl;
  myfile.close(); 

  // myfile.open("block.vtk"); 
  // myfile << d2 <<endl; 
  // myfile.close(); 
 
  // (grid.cbegin())->splitQuad(6,10); 
  // (grid.cbegin()+3)->splitQuad(2, 4); 
  // (grid.cbegin()+1)->splitQuad(2, 10); 
  // (grid.cbegin()+2)->splitQuad(2, 10); 
  // (grid.cbegin()+4)->splitQuad(2, 10); 
  // (grid.cbegin()+5)->splitQuad(2, 10); 
  
  // (grid.cbegin()+6)->splitQuad(3, 1); 
  // (grid.cbegin()+7)->splitQuad(3, 1); 
  // (grid.cbegin()+8)->splitQuad(3, 1); 
  // (grid.cbegin()+9)->splitQuad(6, 1);
  

  //(*(grid->listFace.begin()))->slice(0, 4);//, 4); 
  short int map1[4][2]={{1,2},{0,3},{2,3},{0,1}};

  // cout << "Vertex List" << endl; 
  // for (auto p = grid.vbegin(); p != grid.vend(); ++p) {
  //   cout << *p << " ["; 
  //   for (const auto& l : p->node) {
  // 	cout <<  l << " ";
  //   }
  //   cout << "] " <<endl;
  // }
  // cout << endl;
  
  cout << "Cell list" << endl; 
  vector<Cell> faces;
    for (auto f = block->cbegin(); f != block->cend(); ++f) {
    cout << *f << " "; 
    cout << "faces "; 
    faces = f->ngbrAllFaces(); 
    cout << "[" << faces.size() << "] ("; 
    for (auto c = faces.begin(); c != faces.end(); ++c) {
       cout << *c << " ("<<c->Area()<<") - "; 
    }
    cout << ") "<< endl;
    faces.clear(); 
  }
  cout << endl; 

  myfile.open("grid2d.vtk"); 
  myfile << grid <<endl; 
  myfile.close(); 


  // auto c = grid->cend()-1; 
  // cout << "-----------"<< endl << *(c) << endl;   
  // faces = c->ngbrFaces(2); 
  // cout<< *(c->ngbr(0,3)) << endl << *(c->ngbr(3,0)) << endl; 
  // for (auto f = faces.begin(); f != faces.end(); ++f) {
  //   cout << *f << endl; 
  // }
  // cout << "===========" <<endl; 


  // MatX<double> ma(5000); //ma.uncompress(); 
  // ma[4999][250] = 1; 
    
  // ma.info(); 

  LinSys vel(6); 
  VecX<double> x(6); x.uncompress();
  vel.A = {
    { 2, -1,  0,  0,  0,  0}, 
    {-1,  2, -1,  0,  0,  0}, 
    { 0, -1,  2, -1,  0,  0}, 
    { 0,  0, -1,  2, -1,  0}, 
    { 0,  0,  0, -1,  2, -1}   
  };
  vel.A[5] = { 0,  0,  0,  0, -1,  2} ;
  //vel.A[5][4] = -0.5; 

  //vel.A.info(); 

  x = {1, -2, 3, -4, 9, 3}; 
  vel.b = vel.A*x;
  cout << vel.b << endl;
  x = {0, 0, 0, 0, 0, 0}; 

  vel.x = &x;
  vel.setLimits(1e-8); 

  //vel.Gauss_Seidel();
  vel.CG();
  cout <<"x: "<< x<< endl;


  Vec3 a({1, 1, 0}); 
  Vec3 b({-1, 1, 0}); 
  cout << endl; 
  cout << 0.5*abs(a^b) << endl; 

  return 0; 
};
