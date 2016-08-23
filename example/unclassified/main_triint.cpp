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
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <time.h>
#include <iomanip>

#include <grid2/Grid>

int main(int argc, char* argv[]) {
  double x=0.6, y=0.0, z=0.0; 
  Grid* grid = new Block2({0, 0, 0}, {1, 1, 0}, 2, 2);
  
  if (argc >= 2) x = atof(argv[1]); 
  if (argc >= 3) y = atof(argv[2]); 

  auto pt = Vec3(x, y, z); 

  grid->listCell[0]->adapt[0] = 1; 
  grid->listCell[0]->adapt[1] = 1; 
  grid->listCell[1]->adapt[0] = 1; 
  grid->listCell[2]->adapt[1] = 1; 
  grid->adapt(); 
  grid->levelHighBound[0]=6; 
  grid->levelHighBound[1]=6; 

  for (auto c : grid->listCell) {
    if (c->isPointIn(pt)) 
      cout << "cell " << c->id << " contains pt= " << pt << endl; 
  }
  // exit(1); 

  auto i0 = grid->searchVertexbyCoords(pt); 
  cout << "Vertex " << i0 << " : " << *grid->listVertex[i0] << endl; 
  for (auto ci : grid->listVertex[i0]->cell) {
    if (ci < 0) continue; 
    cout << "Cells: " << ci << endl; 
  }

  //exit(1); 


  grid->addVar({"T"}); 
  
  auto T = grid->getVar("T");
  T->setBC("west", "grad", 1);
  T->setBC("south", "grad", 2); 
  T->setBC("east", "grad", 1); 
  T->setBC("north", "grad", 2); 
  
  double pi = 4.0*atan(1); 
  for (auto i = 0; i < grid->listCell.size(); ++i) {
    auto c = grid->listCell[i]; 
    auto x = c->getCoord(); 
    T->set(i, x[0] + 2*x[1]); 
    //T->set(i, sin(x[0]*2*pi)*sin(x[1]*2*pi)); // + cos(x[0]*4*pi)*sin(x[1]*2*pi)); 
  }
  
  for (auto k = 0; k < 2; ++k) {
    for (auto i = 0; i < grid->listCell.size(); i++) {
      grid->listCell[i]->adapt[0] = 1; 
      grid->listCell[i]->adapt[1] = 1; 
    }
    grid->adapt(); 
  }  


  Grid* other = new Block2({0.05, 0.05, 0}, {0.95, 0.95, 0}, 51, 78); 
  other->addVar({"T"}); 
  auto Tn = other->getVar("T"); 

  for (auto c:other->listCell) {    
    auto x = c->getCoord();
    Tn->set(c->id, grid->interp(T, x)); 
  }  

  grid->writeVTK("int2_"); 
  other->writeVTK("oth_"); 

  //exit(1); 

   for (auto k = 0; k < 2; ++k) {
    for (auto i = 0; i < grid->listCell.size(); i++) {
      grid->listCell[i]->adapt[0] = -1; 
      //      grid->listCell[i]->adapt[1] = -1; 
    }
    grid->adapt(); 
  }   

  // for (auto i = 0; i < grid->listCell.size(); ++i) {
  //   auto c = grid->listCell[i]; 
  //   auto x = c->getCoord(); 
  //   //T->set(i, x[0] + 2*x[1]); 
  //   T->set(i, sin(x[0]*2*pi)*sin(x[1]*2*pi)); // + cos(x[0]*4*pi)*sin(x[1]*2*pi)); 
  // }

  for (auto c:other->listCell) {
    auto x = c->getCoord();
    auto i0 = grid->searchVertexbyCoords(x); //grid->searchInterpPoint(x, grid->searchVertexbyCoords(x)); 
    if (i0 < 0) {
      // cout << "search of << " << x << " failed " << endl;
      // cout << "result of init search " << grid->searchVertexbyCoords(x) << endl; 
    } else {
      Tn->set(c->id, grid->listVertex[i0]->evalPhi(T, &x)); 
    }
  }  

  grid->writeVTK("int2_"); 
  other->writeVTK("oth_");

  ofstream out;
  out.open("interp.vtk"); 
 

  int_8 nCell = grid->listVertex.size();
  int t; // = listCell[0]->node.size()+1;
  out << "# vtk DataFile Version 2.0" << endl; 
  out << "Unstructure Grid" << endl; 
  out << "ASCII"<< endl; 
  out << endl << "DATASET UNSTRUCTURED_GRID"<< endl; 
  out << "POINTS " << 4*nCell << " float" << endl;
  for (auto v : grid->listVertex) {
    out << grid->int3D.linFunc(Vec3(0,0,0), v->xcoef) << " "; 
    out << grid->int3D.linFunc(Vec3(0,0,0), v->ycoef) << " "; 
    out << grid->int3D.linFunc(Vec3(0,0,0), v->zcoef) << endl; 
    out << grid->int3D.linFunc(Vec3(1,0,0), v->xcoef) << " "; 
    out << grid->int3D.linFunc(Vec3(1,0,0), v->ycoef) << " "; 
    out << grid->int3D.linFunc(Vec3(1,0,0), v->zcoef) << endl; 
    out << grid->int3D.linFunc(Vec3(1,1,0), v->xcoef) << " "; 
    out << grid->int3D.linFunc(Vec3(1,1,0), v->ycoef) << " "; 
    out << grid->int3D.linFunc(Vec3(1,1,0), v->zcoef) << endl; 
    out << grid->int3D.linFunc(Vec3(0,1,0), v->xcoef) << " "; 
    out << grid->int3D.linFunc(Vec3(0,1,0), v->ycoef) << " "; 
    out << grid->int3D.linFunc(Vec3(0,1,0), v->zcoef) << endl;     
  }
  auto icnt = 0; 
  out << endl << "CELLS "<< nCell << " " << nCell*5 << endl; 
  for (auto v : grid->listVertex) { 
    out << "4 " << icnt << " " << icnt+1 << " " << icnt+2 << " " << icnt+3 << endl; 
    icnt += 4; 
  }
  out << endl << "CELL_TYPES " << nCell <<endl; 
  for (auto v : grid->listVertex) {
    out << "9 "<< endl; 
  }
  out.close(); 
  return 0; 
}
