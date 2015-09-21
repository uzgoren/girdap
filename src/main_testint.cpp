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
#include <math.h>
#include <fstream>
#include <sstream>
#include <time.h>
#include <iomanip>

#include <field/Var>
#include <grid2/Grid>


int main() {
  std::string flname; 
  int filecnt = 0; 
  ofstream myfile; 

  // GRID 
  Block2* grid = new Block2({0, 0, 0}, {1.0, 1.0, 0}, 20, 20);

  // FIELD variables;
  grid->addVar({"p", "gpx", "gpy", "cgpx", "cgpy"}); 
  grid->addVar("fgpx", 2); 
  grid->addVar("fgpy", 2); 

  auto p = grid->getVar("p"); 
  auto gpx = grid->getVar("gpx"); 
  auto gpy = grid->getVar("gpy");
  auto cgpx = grid->getVar("cgpx"); 
  auto cgpy = grid->getVar("cgpy"); 
  auto u = grid->getVar("u"); 
  auto v = grid->getVar("v"); 
  auto w = grid->getVar("w"); 
  auto fgpx = grid->getVar("fgpx"); 
  auto fgpy = grid->getVar("fgpy"); 

  // set the field variable
  // write face, node based and cell based to a file; 
  auto t = clock(); 
  for (auto i = 0; i < grid->listCell.size(); ++i) { 
   auto c = grid->listCell[i]; 
    auto x = c->getCoord(); 
    double pi = 3.1456; 
    p->set(i, sin(x[0]*4*pi)*cos(x[1]*2*pi)); // + cos(x[0]*4*pi)*sin(x[1]*2*pi)); 
    gpx->set(i, 4*pi*cos(x[0]*4*pi)*cos(x[1]*2*pi)); 
    gpy->set(i, 2*pi*sin(x[0]*4*pi)*sin(x[1]*2*pi)); 
  }
  t = clock() - t; 
  cout << " Exact time: "<< t/(double) CLOCKS_PER_SEC << " secs"<< endl; 

  t = clock(); 
  for (auto j = 0; j< 2; ++j) {
    grid->solBasedAdapt(gpx->data, 1.0); 
    grid->solBasedAdapt(gpy->data, 1.5); 
    grid->refine2();
    for (auto i = 0; i < grid->listCell.size(); ++i) { 
      auto c = grid->listCell[i]; 
      auto x = c->getCoord(); 
      double pi = 3.1456; 
      p->set(i, sin(x[0]*4*pi)*cos(x[1]*2*pi)); // + cos(x[0]*4*pi)*sin(x[1]*2*pi)); 
      gpx->set(i, 4*pi*cos(x[0]*4*pi)*cos(x[1]*2*pi)); 
      gpy->set(i, -2*pi*sin(x[0]*4*pi)*sin(x[1]*2*pi)); 
    } 
  }
  t = clock() - t; 
  cout << " Adaptation time: "<< t/(double) CLOCKS_PER_SEC << " secs"<< endl; 


  auto ix = 186; 
  cout << "exact:  " << gpx->get(ix) << " " << gpy->get(ix) <<  " 0 " << endl; 
  auto res = grid->listCell[ix]->grad(p).eval(p);
  cout << "method: " << res << endl;

  t = clock(); 
  for (auto i = 0; i < grid->listCell.size(); ++i) {
    auto c = grid->listCell[i]; 
    auto gp =  c->grad(p).eval(p); 
    cgpx->set(i, gp[0]); 
    cgpy->set(i, gp[1]); 
    //    w->set(c->id, gp[2]); 
    if (i == ix) cout << "cgps: " << gp << endl;
  }
  t = clock() - t; 
  cout << " Cell grad time: "<< t/(double) CLOCKS_PER_SEC << " secs"<< endl; 


  t = clock(); 
  grid->lockBC(p); 
  for (auto i = 0; i < grid->listFace.size(); ++i) {
    auto f = grid->listFace[i]; 
    auto gp =  f->grad(p).eval(p); 
    fgpx->set(i, gp[0]); fgpy->set(i, gp[1]); 
  }
  t = clock() - t; 
  cout << " Face grad time: "<< t/(double) CLOCKS_PER_SEC << " secs"<< endl; 


  grid->writeFace("face.vtk"); 


  cout << " -------------------------- " << endl; 
  t=clock(); 
  // for (auto i = 0; i < grid->listFace.size(); ++i) {
  //   auto f = grid->listFace[i]; 
  //   auto pf = f->phi(p).eval(p); 
  //   auto flux = pf*f->vol();
  //   if (f->next >= 0) {
  //     u->data[f->next] -= flux[0]; 
  //     v->data[f->next] -= flux[1]; 
  //     w->data[f->next] -= flux[2];
  //     if (f->next == ix) {
  // 	cout << " Face: " << f->id << " orient: " << f->orient << " next: " << f->next; 
  // 	cout << " prev: " << f->prev << " flux " << -flux << endl; 
  //     }
  //   }
  //   if (f->prev >= 0) {
  //     u->data[f->prev] += flux[0]; 
  //     v->data[f->prev] += flux[1]; 
  //     w->data[f->prev] += flux[2];
  //     if (f->prev == ix) {
  // 	cout << " Face: " << f->id << " orient: " << f->orient << " next: " << f->next; 
  // 	cout << " prev: " << f->prev << " flux " << flux << endl; 
  //     }
  //   }
  // }

  for (auto i = 0; i <grid->listCell.size(); ++i) {
    auto c = grid->listCell[i]; 
    double phix=0, phiy=0, sx=0, sy=0; 
    for (auto j :c->face) {
      auto f = grid->listFace[j]; 
      auto area = f->vol().abs(); 
      if (f->orient == 1) {	
	phix += (f->grad(p).eval(p).x())*area; 
	sx += area; 
      } else if(f->orient == 0) { 
	phiy += (f->grad(p).eval(p).y())*area; 
	sy += area; 
      }
    }
    u->data[i] = phix/sx; 
    v->data[i] = phiy/sy; 
  }
  t = clock() - t; 
  cout << " thru average time: "<< t/(double) CLOCKS_PER_SEC << " secs"<< endl; 


  cout << "cool: " << u->data[ix] << " " << v->data[ix] << " " << w->data[ix] << endl;
  grid->unlockBC(); 


  // calculate the gradient of p; 
  // write face, node based and cell based to a file; 

  

  // adapt grid randomly and redo above; 

  flname = "cav"+std::to_string(filecnt++)+".vtk"; 
  myfile.open(flname); 
  myfile << grid << endl;
  myfile.close(); 
  
 
  delete(grid); 

  return 0; 
}
