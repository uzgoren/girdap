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
#ifndef GRID
#define GRID

#include <memory>
#include <typeinfo>
#include <cmath>
#include <time.h>
#include <iostream>
#include <sstream>
#include <fstream>

#include <field/Var>
#include <grid2/Vertex>
#include <grid2/Cell>

class LinSys;  

//class Vertex; 
//class Var; 

class Grid {
public:
  vector<shared_ptr<Vertex > > listVertex; 
  vector<shared_ptr<Cell > > listCell; 
  vector<shared_ptr<Cell > > listFace; 
  vector<shared_ptr<Var > > listVar;

  // Matrix for bilinear coordinate transformation
  MatX<double> AI; 

  //vector<shared_ptr<Boundary > > listBNDR; 
  shared_ptr<Var > thisVar; 
  double dt, dt0, cfl; 

  int_8 nCellSize; 
  int nFace; 

  int_2 levelLowBound[3], levelHighBound[3]; 
  int_2 levelMin[3], levelMax[3]; 

  double minD, maxD, meanD; 
  // vector<shared_ptr<BNDR> > 

  Grid() { nCellSize = 0; for (auto i = 0; i<3; ++i) {levelLowBound[i] = 0; levelHighBound[i] = 4; cfl = 0.5; }; setAI(); }

  Grid(initializer_list<initializer_list<double > > pts):Grid() {
    addVertex(pts);
  }
  
  Grid(Vec3 pts): Grid() {
    addVertex(pts); 
  }

  Grid(initializer_list<initializer_list<double> > pts, 
       initializer_list<initializer_list<int_8> > cell): Grid(pts) {
    addCell(cell);
    // correct connectivity ! missing (hanging nodes) for the next step; 
    // makeFace(); 
  }

  // 00 Set transformation matrix; 
  void setAI() { 
    AI = { { 1, 0, 0, 0, 0, 0, 0, 0}, 
	   {-1, 1, 0, 0, 0, 0, 0, 0},
           {-1, 0, 0, 1, 0, 0, 0, 0},
           {-1, 0, 0, 0, 1, 0, 0, 0},
           { 1,-1, 1,-1, 0, 0, 0, 0},
           { 1,-1, 0, 0,-1, 1, 0, 0},
           { 1, 0, 0,-1,-1, 0, 0, 1},
           {-1, 1,-1, 1, 1,-1, 1,-1} };
  }
  
  // Add vertex to the list; 
  void addVertex(Vec3 v) {
    listVertex.emplace_back( shared_ptr<Vertex>(new Vertex(v)) ); 
    (*listVertex.rbegin())->id = listVertex.size()-1;
    (*listVertex.rbegin())->grid = this; 
  }     

  void addVertex(initializer_list<double > pt) {
     addVertex(Vec3(pt)); 
  }

  void addVertex(vector<Vec3 > pts) {
    for (auto pt : pts) { addVertex(pt); }
  }

  void addVertex(initializer_list<initializer_list<double > > pts) {
     for (auto i = pts.begin(); i != pts.end(); ++i) { addVertex(*i); }    
  }

  // void forCells(vector<shared_ptr<Cell > >::iterator cit, function fn()) {
  //   vector<bool> visited(listCell.size(), false); 
  //   vector<int_8> queue(1, (*cit)->id); 
  //   auto qit = queue.begin(); 
  //   while (qit != queue.end()) {
  //     c = listCell[*qit]; 
  //     auto ngbrList = c->ngbrCellList(); 
      
  //   }
  // }

  // // //Add cell to the list
  void addCell(initializer_list<int_8 > l) {    
    if (l.size() == 2) {
  //     //listCell.emplace_back(shared_ptr<Cell > (new Line(tmp)) );
    } else if (l.size() == 3) {
  //     //listCell.emplace_back(shared_ptr<Cell > (new Tri(tmp)) ); 
    } else if (l.size() == 4) {
      listCell.emplace_back( shared_ptr<Cell > (new Quad(l)) ); 
    } else if (l.size() == 8) {
       //listCell.emplace_back( shared_ptr<Cell > (new Hexa(tmp)) ); 
    } else {
      cout << "Cell type not understood!"<<endl; 
      exit(1); 
    }
    nCellSize += l.size() + 1;
    auto c = *(listCell.rbegin()); 
    c->id = listCell.size()-1;
    c->grid = this; 
    c->assignCelltoNode();
    for (auto v: listVar) {v->data.push_back(0); v->data.uncompress(); }
    c->masterx.resize(levelHighBound[0]+1, false); 
    c->mastery.resize(levelHighBound[1]+1, false); 
    c->masterz.resize(levelHighBound[2]+1, false); 
  }

  void addCell(initializer_list<initializer_list<int_8 > > link) {
    for (auto it = link.begin(); it != link.end(); ++it) {
      addCell(*it);
    }
  }

  void setDt(double c) {
    auto U = getVar("u"); auto V = getVar("v"); auto W = getVar("w");     
    auto maxU = abs(min(U->data.min(), min(V->data.min(), W->data.min()))); 
    cout << maxU << " ";
    maxU = max(maxU, max(U->data.max(), max(V->data.max(), W->data.max())));
    cout << maxU << endl; 
    double dx = 1e10;  
    for (auto c : listCell) {
      double cdx, cdy, cdz; 
      cdx = c->dx(); cdy = c->dy(); cdz = c->dz(); 
      if (cdx > 0) dx = min(dx, cdx); 
      if (cdy > 0) dx = min(dx, cdy); 
      if (cdz > 0) dx = min(dx, cdz); 
    }
    dt0 = dt;
    if (maxU < 1e-6 || dx < 1e-6)
      dt = c;
    else 
      dt = min(c, cfl*dx/maxU);
    cout << dt << " " << dx << " " << maxU << endl; 
  }

  void makeFace() {
    listFace.clear(); 
    for (auto c : listCell)  
      c->face.clear(); 
    
    for (auto v : listVertex) {
      if (v->cell.size() == 2) {
	listFace.emplace_back(new Cell({v->id})); 
	(*listFace.rbegin())->next = (v->cell[1] < 0) ? -1 : v->cell[1]; 
	(*listFace.rbegin())->prev = (v->cell[0] < 0) ? -2 : v->cell[0];
	(*listFace.rbegin())->grid = this;
	(*listFace.rbegin())->orient = 0; 
	(*listFace.rbegin())->id = listFace.size()-1; 
	if (v->cell[1] >= 0) listCell[v->cell[1]]->face.push_back((*listFace.rbegin())->id) ;
	if (v->cell[0] >= 0) listCell[v->cell[0]]->face.push_back((*listFace.rbegin())->id); 
	return; 
      } else if (v->cell.size() == 4) { 
	auto n = v->ngbr(1); 
	if (n && *n) {
	  listFace.emplace_back(new Line({v->id, (*n)->id})); 
	  (*listFace.rbegin())->next = (v->cell[2] < 0) ? -2 : v->cell[2];
	  (*listFace.rbegin())->prev = (v->cell[1] < 0) ? -1 : v->cell[1];
	  if (v->cell[2] < 0) v->cell[2] = -2; 
	  if (v->cell[1] < 0) v->cell[1] = -1;
	  (*listFace.rbegin())->grid = this;
	  (*listFace.rbegin())->orient = 0; 
	  (*listFace.rbegin())->id = listFace.size()-1; 
	  if (v->cell[1] >= 0) listCell[v->cell[1]]->face.push_back((*listFace.rbegin())->id); 
	  if (v->cell[2] >= 0) listCell[v->cell[2]]->face.push_back((*listFace.rbegin())->id); 
	} 
	if (v->ngbr(-1)) { 
	  if (v->cell[0] < 0) v->cell[0] = -1; 
	  if (v->cell[3] < 0) v->cell[3] = -2; 
	}
	n = v->ngbr(2); 
	if (n && *n) { 
	  listFace.emplace_back(new Line({(*n)->id, v->id})); 
	  (*listFace.rbegin())->next = (v->cell[2] < 0) ? -4 : v->cell[2];
	  (*listFace.rbegin())->prev = (v->cell[3] < 0) ? -3 : v->cell[3];
	  if (v->cell[2] < 0) v->cell[2] = -4; 
	  if (v->cell[3] < 0) v->cell[3] = -3; 
	  (*listFace.rbegin())->grid = this;
	  (*listFace.rbegin())->orient = 1; 
	  (*listFace.rbegin())->id = listFace.size()-1; 
	  if (v->cell[2] >= 0) listCell[v->cell[2]]->face.push_back((*listFace.rbegin())->id); 
	  if (v->cell[3] >= 0) listCell[v->cell[3]]->face.push_back((*listFace.rbegin())->id); 
	} 
	if (v->ngbr(-2)) { 
	  if (v->cell[0] < 0) v->cell[0] = -3; 
	  if (v->cell[1] < 0) v->cell[1] = -4; 
	}
      }	
    }
    nFace = listFace.size();

    for (auto v : listVar) {
      if (v->loc == 2) v->data.resize(nFace); 
    }
  }

  void addVar(std::string n, int t = 0); 
  void addVar(initializer_list<std::string> nl);
  shared_ptr<Var> getVar(std::string n);   
  void lockBC(shared_ptr<Var> v);
  void unlockBC(); 


  VecX<double> getCoord(int dir) {
    VecX<double> val(listCell.size()); val.uncompress(); 
    for (auto i = 0; i < listCell.size(); ++i) {
      double sum=0; double icnt = 0; 
      for (auto j = 0; j < listCell[i]->node.size(); ++j) {
	if (listCell[i]->node[j] < 0) continue; 
	auto k = listCell[i]->getVertex(j); 
	sum = sum + (**k)[dir]; icnt++; 
      }
      val[i] = sum/icnt; 
      //cout << sum  << " " << icnt << " " << val[i] << endl; 
    }
    return val; 
  }

  VecX<Vec3> getVel() { 
    VecX<Vec3> val(listCell.size()); val.uncompress(); 
    auto u = *(listVar.begin()); auto v = *(listVar.begin()+1); auto w = *(listVar.begin()+2);
    for (auto i = 0; i< listCell.size(); ++i) {
      val[i].set({u->get(i), v->get(i), w->get(i)});
    }
    return val; 
  }
  
  VecX<Vec3> getVelVertex() {
    VecX<Vec3> val(listVertex.size()); 
    auto u = getVar("u"); 
    auto v = getVar("v"); 
    auto w = getVar("w"); 

    for (auto i = 0; i<listVertex.size(); ++i) {
      // Vec3 sum(0,0,0); double icnt = 0;  
      // for (auto j = 0; j < listVertex[i]->cell.size(); ++j) {
      // 	if (listVertex[i]->cell[j] < 0) continue;
      // 	auto k = listVertex[i]->cell[j]; 
      // 	sum = sum + Vec3(u[k], v[k], w[k]); 
      // 	++icnt; 
      // }
      val[i][0] = listVertex[i]->phi(u).eval(u); 
      val[i][1] = listVertex[i]->phi(v).eval(v); 
      val[i][2] = listVertex[i]->phi(w).eval(w); 

    }
    return val; 
  }

  VecX<double> getPhiVertex(shared_ptr<Var> phi) {
    VecX<double> val(listVertex.size()); 
    for (auto i = 0; i<listVertex.size(); ++i) {
      // double sum=0; double icnt = 0;  
      // for (auto j = 0; j < listVertex[i]->cell.size(); ++j) {
      // 	if (listVertex[i]->cell[j] < 0) continue; 
      // 	auto k = listVertex[i]->cell[j]; 
      // 	sum = sum + phi->data[k]; 
      // 	++icnt; 
      // }
      val[i] = listVertex[i]->phi(phi).eval(phi); //sum/icnt; 
    }
    return val; 
  }


  void solBasedAdapt(VecX<double> phi, double a=1.0, double b = 0.5) {
    double icnt=0, sumx=0, varx=0, sumy=0, vary=0, sumz=0, varz=0; 
    //double a = 0.5; 
    if (phi.size() != listCell.size()) {
      cout << "solBasedAdapt: CANNOT process a variable with a different size than cell size!"<<endl;
      exit(1); 
    }
    for (auto i = 0; i < listCell.size(); ++i) { 
      auto qx = abs(phi[i])*pow(listCell[i]->dx(),1.5);  
      auto qy = abs(phi[i])*pow(listCell[i]->dy(),1.5);  
      auto qz = abs(phi[i])*pow(listCell[i]->dz(),1.5);  
      ++icnt; 
      sumx += qx; varx += qx*qx; 
      sumy += qy; vary += qy*qy; 
      sumz += qz; varz += qz*qz; 
    }
    varx = sqrt((varx - sumx*sumx/icnt)/(icnt-1)); 
    vary = sqrt((vary - sumy*sumy/icnt)/(icnt-1)); 
    varz = sqrt((varz - sumz*sumz/icnt)/(icnt-1)); 
    sumx /=icnt; 
    sumy /=icnt; 
    sumz /=icnt; 
    
    for (auto i = 0; i < listCell.size(); ++i) {
      auto c=listCell.at(i);
      auto qx = abs(phi[i])*pow(c->dx(),1.5);  
      auto qy = abs(phi[i])*pow(c->dy(),1.5);  
      auto qz = abs(phi[i])*pow(c->dz(),1.5);  

      if ((qx - sumx) > a*varx)
	c->adapt[0] = min(1, levelHighBound[0] - c->level[0]);
      if ((qy - sumy) > a*vary)
	c->adapt[1] = min(1, levelHighBound[1] - c->level[1]);
      if ((qz - sumz) > a*varz)
	c->adapt[2] = min(1, levelHighBound[2] - c->level[2]);
    } 
  }

  void solBasedAdapt2(VecX<Vec3> phi) {
    // refine if err > 0.03; coarsen if err < 1e-6; 
    for (auto i = 0; i < listCell.size(); ++i) {
      auto c = listCell[i]; 
      auto qx = abs(phi[i].x())*listCell[i]->dx();  
      auto qy = abs(phi[i].y())*listCell[i]->dy();  
      auto qz = abs(phi[i].z())*listCell[i]->dz();  
      
      if (qx > 8e-5) 
	c->adapt[0] = min(1, levelHighBound[0] - c->level[0]);
      else if (qx < 5e-6) 
	c->adapt[0] = max(-1, levelLowBound[0] - c->level[0]);
      if (qy > 8e-5) //      if (0.5*qy > sumy)
	c->adapt[1] = min(1, levelHighBound[1] - c->level[1]);
      else if (qy < 5e-6)
	c->adapt[1] = max(-1, levelLowBound[1] - c->level[1]);
      if (qz > 8e-5) //if (0.5*qz > sumz)
	c->adapt[2] = min(1, levelHighBound[2] - c->level[2]);
      else if (qz < 5e-6)
	c->adapt[2] = max(-1, levelLowBound[2] - c->level[2]);      
    }       
  }

  void solBasedAdapt(VecX<Vec3> phi, double a=1.0, double b=0.1) {
    double icnt=0, sumx=0, varx=0, sumy=0, vary=0, sumz=0, varz=0;
    // double a = 0.5; 
    if (phi.size() != listCell.size()) {
      cout << "solBasedAdapt: CANNOT process a variable with a different size than cell size!"<<endl;
      exit(1); 
    }
    for (auto i = 0; i < listCell.size(); ++i) {      
      auto qx = abs(phi[i].x())*pow(listCell[i]->dx(),1.5);  
      auto qy = abs(phi[i].y())*pow(listCell[i]->dy(),1.5);  
      auto qz = abs(phi[i].z())*pow(listCell[i]->dz(),1.5);  
      ++icnt; 
      sumx += qx; varx += qx*qx; 
      sumy += qy; vary += qy*qy; 
      sumz += qz; varz += qz*qz; 
    }
    varx = sqrt((varx - sumx*sumx/icnt)/(icnt-1)); 
    vary = sqrt((vary - sumy*sumy/icnt)/(icnt-1)); 
    varz = sqrt((varz - sumz*sumz/icnt)/(icnt-1)); 
    sumx /=icnt; 
    sumy /=icnt; 
    sumz /=icnt; 

    for (auto i = 0; i < listCell.size(); ++i) {
      auto c=listCell.at(i); 
      auto qx = abs(phi[i].x())*pow(c->dx(),1.5);  
      auto qy = abs(phi[i].y())*pow(c->dy(),1.5);  
      auto qz = abs(phi[i].z())*pow(c->dz(),1.5);  
      if (qx - sumx > a*varx)
	c->adapt[0] = min(1, levelHighBound[0] - c->level[0]);
      else if (qx - sumx < b*varx)
	c->adapt[0] = max(-1, levelLowBound[0] - c->level[0]);
      if (qy - sumy > a*vary)
	c->adapt[1] = min(1, levelHighBound[1] - c->level[1]);
      else if (qy - sumy < b*vary)
	c->adapt[1] = max(-1, levelLowBound[1] - c->level[1]);
      if (qz - sumz > a*varz)
	c->adapt[2] = min(1, levelHighBound[2] - c->level[2]);
      else if (qz - sumz < b*varz)
	c->adapt[2] = max(-1, levelLowBound[2] - c->level[2]);
    }
  }

  void setCurrentLevels() { 
    if (listCell.empty()) {
      cout << "setCurrentLevels: There is no cell to set levels!"<< endl; 
      exit(1); 
    }
    for (auto j = 0; j < 3; ++j) { 
      levelMin[j] = listCell.at(0)->level[j]; 
      levelMax[j] = listCell.at(0)->level[j];
    } 
    for (auto i = 0; i < listCell.size(); ++i) {
      for (auto j = 0; j < 3; ++j) {
	levelMin[j] = min(levelMin[j], listCell[i]->level[j]); 
	levelMax[j] = max(levelMax[j], listCell[i]->level[j]); 
      }
    }
  }


  void updateLevels() {
    // in x-direction; 
    for (auto j = 0; j < 3; ++j) {
      for (auto lx = levelMax[j]; lx > levelMin[j]-1; --lx) {
	for (auto i = 0; i < listCell.size(); ++i) {
	  auto level0 = listCell[i]->level[j]; 
	  auto adapt0 = listCell[i]->adapt[j]; 
	  if (level0 + adapt0 != lx + 1) continue; // true maximum 
	  auto ngbrList = listCell[i]->ngbrCellList(); 
	  for (auto n : ngbrList) {
	    auto level1 = listCell[n]->level[j]; 
	    auto adapt1 = listCell[n]->adapt[j];
	    if (level1 + adapt1 > level0 + adapt0 - 2) continue; 
	    listCell[n]->adapt[j] = adapt0+level0-1;
	  }
	}
      }
    }
  }
  

  void adapt() { 

    auto t = clock(); 
    // [0] FIND cells that needs to be adapted for continuous levels
    for (auto j = 0; j < 3; ++j) {
      for (auto l = levelMax[j]; l > levelMin[j]-1; --l) {
	for (auto i = 0; i < listCell.size(); ++i) {
	  auto c = listCell[i]; 
	  if (c->level[j] != l) continue; 
	  c->checkNgbrLevel(j); 
	}
      }
    }
    //
    for (auto j = 0; j < 3; ++j) {
      for (auto l = levelMax[j]; l > levelMin[j]-1; --l) {
	for (auto i = 0; i < listCell.size(); ++i) {
	  auto c = listCell[i]; 
	  if (c->level[j] != l) continue; 
	  c->checkIslandLevels(j); 
	}
      }
    }
    // [2] Refine cells;
    for (auto j = 0; j < 3; ++j) {
      auto lmin = levelMin[j]; auto lmax = levelMax[j]+1;
      auto l = lmin; 
      while (l < lmax) {
	auto icnt = 0; 
	auto nCell = listCell.size();
    	for (auto i = 0; i < nCell; ++i) {
	  auto c = listCell[i]; 
    	  if (c->level[j] != l) continue;
    	  if (c->adapt[j] <= 0) continue;
	  c->refine(j); 
	}
	++l;
      }
    }

    bool isCellRemoved = false; 
    // [3] Coarsen cells (index directions - 3 to be defined)
    for (auto pass =0; pass < 3; ++pass) { 
      for (auto j = 0; j < 2; ++j) {	
	auto lmin = levelMin[j]; auto lmax = levelMax[j];
	for (auto l = lmax; l > lmin-1; --l) {
	  for (auto i = 0; i < listCell.size(); ++i) {
	    auto c = listCell[i]; 
	    if (!c->isAlive) continue; 
	    if (c->adapt[j] >= 0) continue;
	    if (c->level[j] != l) continue; 
	    if (c->coarsen(j) && !isCellRemoved) { 
	      isCellRemoved = true;
	    }
	  }
	}
      }
    }

    // [4] correct cell lists;     
    if (isCellRemoved) { 
      vector<int_8> o2n_cell(listCell.size(), -1); int_8 icell = 0; 
      for (auto it = cbegin(); it != cend(); ) {
	if ((*it)->isAlive) {
	  o2n_cell[(*it)->id] = icell++; 
	  ++it; 
	} else { 
	  if ((*it).use_count() > 1) {
	    cout << "Cannot clear cell list: USED a lot "<<  endl; 
	    exit(1); 
	  }
	  it = listCell.erase(it); 
	}	
      }
      vector<int_8> o2n_node(listVertex.size(), -1); int_8 inode = 0; 
      for (auto it = vbegin(); it != vend(); ) {
	//      if (false) {
	if (((*it)->cell[0] == (*it)->cell[1] && (*it)->cell[2] == (*it)->cell[3]) ||
	    ((*it)->cell[0] == (*it)->cell[3] && (*it)->cell[1] == (*it)->cell[2])) {
	  if ((*it).use_count() > 1) {
	    cout << " vertex count exceeds one "<< endl; 
	    exit(1); 
	  }
	  it = listVertex.erase(it); 
	} else {      
	  o2n_node[(*it)->id] = inode++;
	  ++it; 
	}
      }

      for (auto j = 0; j < listCell.size(); ++j) {
	auto c = listCell[j]; 
	for (auto var : listVar) {
	  var->data[j] = var->data[c->id]; 
	}
	c->id = j; 
	for (auto i = 0; i < c->node.size(); ++i) {
	  if (c->node[i] >= 0) c->node[i] = o2n_node[c->node[i]]; 
	}
      }

      for (auto j = 0; j < listVertex.size(); ++j) {
	auto v = listVertex[j]; 
	v->id = j; 
	for (auto i = 0; i < v->cell.size(); ++i) {
	  if (v->cell[i] >= 0) v->cell[i] = o2n_cell[v->cell[i]]; 
	} 
      }
      for (auto v : listVar) { 
	if (v->loc == 0) v->data.resize(listCell.size()); 
      }
      o2n_cell.clear(); 
      o2n_node.clear(); 
    }
 
    setCurrentLevels();
    makeFace(); 

    for (auto j = 0; j < 3; ++j) {
      for (auto i = 0; i < listCell.size(); ++i) {
	auto c= listCell[i]; 
	c->adapt[j] = 0;
      }
    }

    t = clock() - t;   
 
    cout << " Cells: " << listCell.size(); 
    cout << " Faces: " << nFace << " in " << ((float)t)/CLOCKS_PER_SEC << " secs."<< endl; 
  }



  void refine2() {
    auto t = clock(); 
    for (auto i = 0; i < listCell.size(); ++i) { 
      for (auto j = 0; j < 3; ++j) {
    	listCell[i]->checkNgbrLevel(j); 
      }
    }
    for (auto i = 0; i < listCell.size(); ++i) { 
      for (auto j = 0; j < 3; ++j) {
	listCell[i]->checkIslandLevels(j); 
      }
    }
    for (auto j = 0; j < 3; ++j) {
      auto lmin = levelMin[j]; auto lmax = levelMax[j]+1;
      auto l = lmin; 
      while (l < lmax) {
	auto icnt = 0; 
	auto nCell = listCell.size();
	for (auto c:listCell) {
	  c->checkNgbrLevel(j); 
	}

    	for (auto i = 0; i < nCell; ++i) {
	  auto c = listCell[i]; 
    	  if (c->level[j] != l) continue;
    	  if (c->adapt[j] == 0) continue;
	  c->refine(j); 
	  icnt++; 
	}
	if (icnt == 0) ++l; 
	setCurrentLevels();
	makeFace(); 
      }
    }
    t = clock() - t;   
 
    cout << " Cells: " << listCell.size(); 
    cout << " Faces: " << nFace << " in " << ((float)t)/CLOCKS_PER_SEC << " secs."<< endl; 

    return; 
  }
    
  bool coarsen(int_2 dir) {
    bool isAnyCoarsened = false; 
    auto lmin = levelMin[dir]-1; auto lmax = levelMax[dir];
    for (auto l = lmax; l > lmin; --l) {
      for (auto i = 0; i < listCell.size(); ++i) { 
	auto c = listCell[i]; 
	if (c->level[dir] != l) continue; 
	if (dir == 0 && !c->masterx[l]) continue; 
	if (dir == 1 && !c->mastery[l]) continue; 
	if (c->coarsen(dir) && !isAnyCoarsened) {
	  isAnyCoarsened=true;
	} 
      }
    }
    return true; 
  }

  void refine() { 
    auto t = clock(); 
    bool isFinished = false; 
    auto maxLevel = 0; auto l = 0;  auto minLevel=0; 
    vector<int > levels; 
    cout << "Refine: "; 
    while (!isFinished) {
      isFinished = true;
      for (auto v : listVertex) {
    	auto c1 = v->getCell(1); 
    	auto c2 = v->getCell(2); 
    	auto c3 = v->getCell(3); 
    	if (c1 && c2 && c1!=c2) {
	  l = (*c1)->adapt[0] + (*c1)->level[0]; if (levels.size() < l+1) {levels.resize(l+1); levels[l]=1;} else ++levels[l]; 
	  l = (*c2)->adapt[0] + (*c2)->level[0]; if (levels.size() < l+1) {levels.resize(l+1); levels[l]=1;} else ++levels[l]; 
	  maxLevel = max(maxLevel, (*c1)->adapt[0] + (*c1)->level[0]);
	  maxLevel = max(maxLevel, (*c2)->adapt[0] + (*c2)->level[0]);
    	  if (((*c1)->adapt[0])+((*c1)->level[0])-((*c2)->adapt[0])-((*c2)->level[0]) > 1 ) {
    	    (*c2)->adapt[0] = ((*c1)->level[0])+((*c1)->adapt[0])-((*c2)->level[0])-1; 
    	    isFinished = false;
	  } else if (((*c2)->adapt[0])+((*c2)->level[0])-((*c1)->adapt[0])-((*c1)->level[0]) > 1 ) {
    	    (*c1)->adapt[0] = ((*c2)->level[0])+((*c2)->adapt[0])-((*c1)->level[0])-1;
    	    isFinished = false; 
    	  }
    	}
    	if (c2 && c3 && c2!=c3) {
	  l = (*c3)->adapt[1] + (*c3)->level[1]; if (levels.size() < l+1) {levels.resize(l+1); levels[l]=1;} else ++levels[l]; 
	  l = (*c2)->adapt[1] + (*c2)->level[1]; if (levels.size() < l+1) {levels.resize(l+1); levels[l]=1;} else ++levels[l]; 
	  maxLevel = max(maxLevel, (*c3)->adapt[1] + (*c3)->level[1]);
	  maxLevel = max(maxLevel, (*c2)->adapt[1] + (*c2)->level[1]);
    	  if (((*c3)->adapt[1])+((*c3)->level[1])-((*c2)->adapt[1])-((*c2)->level[1]) > 1 ) {
    	    (*c2)->adapt[1] = ((*c3)->level[1])+((*c3)->adapt[1])-((*c2)->level[1])-1;
    	    isFinished = false; 
    	  } else if (((*c2)->adapt[1])+((*c2)->level[1])-((*c3)->adapt[1])-((*c3)->level[1]) > 1 ) {
    	    (*c3)->adapt[1] = ((*c2)->level[1])+((*c2)->adapt[1])-((*c3)->level[1])-1;
    	    isFinished = false; 
    	  }	
    	}
     }
    }
    for (auto i=0; i<levels.size(); ++i) if (levels[i] > 0) { minLevel = i; break; }
    cout << " Levels: ("<< minLevel << "-" << maxLevel << ") : "; 

    isFinished = false; 
    auto nCell = 0; 
    while (!isFinished && listCell.size() != nCell) {
       vector<bool > visited(listCell.size(), false);
      nCell = listCell.size(); 
      
      isFinished = true; 
      for (auto l = minLevel; l < maxLevel+1; ++l) {
	int nlevel=0;      
	for (auto i = 0; i< nCell; ++i) {
	  if (((*(cbegin()+i))->level[0] > l) || ((*(cbegin()+i))->level[1] > l)|| visited[i]) continue; 
	  visited[i] = true; 
	  if ((*(cbegin()+i))->refine()) isFinished = false;
	}
      }    
    }
    setCurrentLevels(); 
    makeFace(); 
    t = clock() - t;   
  
    cout << " Cells: " << listCell.size(); 
    cout << " Faces: " << nFace << " in " << ((float)t)/CLOCKS_PER_SEC << " secs."<< endl; 

  }

  VecX<Vec3> getError(shared_ptr<Var> const &a) {
    VecX<Vec3> err(a->data.size()); 
    VecX<double> sum(a->data.size()); 
    auto grad = valGrad(a); 
    for (auto j = 0; j < listFace.size(); ++j) {
      auto f = listFace[j]; 
      auto xf = f->getCoord(); auto area = f->vol(); 
      auto pf1 = f->phi(a).eval(a); 
      auto n = f->next; auto p = f->prev; 
      if (n >= 0) {
	auto xc = listCell[n]->getCoord(); 
	auto pf0 = a->data[n] + grad[n]*(xf - xc); 
	err[n] += abs(pf1 - pf0)*area;
	sum[n] += area.abs(); 
      } 
      if (p >= 0) {
	auto xc = listCell[p]->getCoord(); 
	auto pf0 = a->data[p] + grad[p]*(xf - xc); 
	err[p] += abs(pf1 - pf0)*area;
	sum[p] += area.abs(); 
      } 	
    }   
    for (auto j = 0; j < sum.size(); ++j) {     
      err[j] /= sum[j]; 
    }

    return err; 
  }
	// if (c->adapt[dir] == -1) {
	//   // check east first; 
	//   if (c->node[1] < 0 || c->node[2] < 0) continue; 
	  
	//   // check hanging node at first; 
	//   auto e2 = (*c->getVertex(0)); 
	//   auto e3 = (*c->getVertex(3)); 
	//   if (e2->cell[0]>=0 && (e2->cell[0] == e2->cell[1])) continue; 
	//   if (e3->cell[2]>=0 && (e3->cell[2] == e3->cell[3])) continue; 
	  
	//   auto e0 = (*c->getVertex(1))->cell[2]; 
	//   auto e1 = (*c->getVertex(2))->cell[1]; 
	  
	//   if (e0 >= 0 && e1 == e0) { // tackle only with e0
	//     // CHECK LEVELS;
	//     if (grid->listCell[e0]->adapt[0] >= 0) continue; 
	//     if (grid->listCell[e0]->level[1] != c->level[1]) continue;
	//     e2 = *grid->listCell[e0]->getVertex(1); 
	//     e3 = *grid->listCell[e1]->getVertex(2); 
	//     if (e2->cell[0]>=0 && (e2->cell[0] == e2->cell[1])) continue; 
	//     if (e3->cell[2]>=0 && (e3->cell[2] == e3->cell[3])) continue; 	    
	    
	//     cout << " --> POSSIBLE MATCH : cell "<< i << " -> "<< e0 << " " <<e1 << endl; 

	    
	//     for (auto v : grid->listVar) {
	//       auto vol0 = c->vol().abs(); 
	//       auto vol1 = grid->listCell[e0]->vol().abs(); 
	//       v->set(i, (v->data[i]*vol0 + v->data[e0]*vol1)/(vol0+vol1)); 
	//     }
	    
	//     //check for hanging node; 
	//     auto hv = grid->listCell[e0]->hangingVertexOnFace(1); 	   
	//     if (hv >= 0) {
	//       grid->listVertex[hv]->cell[0] = i; 
	//       grid->listVertex[hv]->cell[3] = i;
	//     } 

	//     for (auto j = 0; j < grid->listCell[e0]->node.size(); ++j) {
	//       auto v = *grid->listCell[e0]->getVertex(j); 
	//       v->cell[(j+2)%4] = i; 
	//     }	    
	    
	//     grid->listCell[e0]->isAlive = false; //node.assign(4, -1); 
	//     grid->listCell[e1]->isAlive = false; //node.assign(4, -1); 
	//     grid->listCell[e0]->node.assign(4,-1); 
	//     grid->listCell[e1]->node.assign(4,-1); 
	    
	//     (*c->getVertex(1))->cell[2] = i; 
	//     (*c->getVertex(2))->cell[1] = i; 
	    
	//     c->node[1] = e2->id; //grid->listCell[e0]->node[1]; 
	//     c->node[2] = e3->id; //grid->listCell[e1]->node[2]; 
	//     c->adapt[0] = 0; 
	//     c->level[0] -= 1; 
	    
	    // cout << c << endl; 	   
	    // for (auto n : c->node) {
	    //   auto v = grid->listVertex[n]; 
	    //   cout << v->cell[0] << " " << v->cell[1] << " " << v->cell[2] << " " << v->cell[3] << endl; 
	    // }
  // 	  }
  // 	}
  //     }
  //   }    
  // }
  

  // void updateHangingVertex() {
  //   for (auto cw: listCell) { 
  //     if (auto c = cw.lock()) c->assignHangingNode();
  //   }
  // }

  // void removeUnusedVertexFromList() {
  //   int_8 icnt = 0;
  //   for (auto v=listVertex.begin(); v!=listVertex.end(); ) {
  //     if ((*v)->isActive()) {
  // 	(*v)->id = icnt; 
  // 	++icnt; ++v; 
  //     }	else {
  // 	listVertex.erase(v); 
  //     }
  //   }
  // }

  // void removeCell(vector<weak_ptr<Cell > >::iterator cit) {
  //   //    (*cit)->reset(); 
  //   listCell.erase(cit); 
  // }
  void adaptCriteria() { // For now it is a circle; 
    Vec3 x0 = Vec3(0.5, 0.5, 0); 
    double tol = 1./10.; 
    for (auto i = 0; i < 2 ; ++i) {
      for (auto v : listVertex) {
	if (abs((Vec3(*v)-x0).abs() - 0.2) <= tol) {
	  auto c = (v)->getCell(0); 
	  if (c && *c) {
	    ++((*c)->adapt[0]); 
	    ++((*c)->adapt[1]); 
	  }
	}
      }
      tol = tol/2; 
      refine();       
    }
  }

  void debug() {
    for (auto v = vbegin(); v!=vend(); ++v) {
      cout << (*v)->id << " ("; 
      for (auto i = 0; i < (*v)->cell.size(); ++i) {
	if ((*v)->cell[i] < 0) {cout << "x "; continue;}
  	cout << (*v)->cell[i] << " "; 
      }
      cout << ") "<<endl; 
    }
  }

  void laplace(LinSys &axb, shared_ptr<Cell > f, double const &c);
  LinSys laplace(double c, initializer_list<double> n={}); 
  LinSys laplace(VecX<double> c, initializer_list<double> n={}); 
  LinSys source(double c, double a, initializer_list<double> n={}); 
  LinSys source(double c, VecX<double> a, initializer_list<double> n={}); 

  void div(LinSys &axb, shared_ptr<Cell > f, Vec3 const &c); 
  LinSys div(VecX<Vec3> vel, double c, initializer_list<double> n={}); 
  LinSys divRK2E(VecX<Vec3> vel, double c); 
  LinSys divRK4E(VecX<Vec3> vel, double c); 
  LinSys div(VecX<Vec3> vel, VecX<double> c, initializer_list<double> n={}); 
  void getBCPhi(shared_ptr<Cell > f, shared_ptr<Var> phi, double &a, double &b);
  void getPhi_val(shared_ptr<Cell > f, shared_ptr<Var> phi, double &out); 
  VecX<double> valDiv(VecX<Vec3> vel);
  VecX<Vec3> valGrad(shared_ptr<Var > phi); 
  void correctVel(double coef);

  LinSys ddt(double c); 
  void timeScheme(LinSys &axb, initializer_list<double> &n, VecX<double> &prev); 
  

  //   cout << "Hanging nodes on faces of a cell" << endl; 
  //   for (auto c : listCell) { 
  //     if (auto cs = c.lock()) {
  // 	  cout << cs->id << "-> "; 
  // 	  cout << "S: "; 
  // 	  if (cs->hangingNodeOnFace(0)) 
  // 	    cout << cs->hangingNodeOnFace(0)->id<<" ";
  // 	  cout << "E: "; 
  // 	  if (cs->hangingNodeOnFace(1)) 
  // 	    cout << cs->hangingNodeOnFace(1)->id<<" ";
  // 	  cout << "N: "; 
  // 	  if (cs->hangingNodeOnFace(2)) 
  // 	    cout << cs->hangingNodeOnFace(2)->id<<" ";
  // 	  cout << "W: "; 
  // 	  if (cs->hangingNodeOnFace(3)) 
  // 	    cout << cs->hangingNodeOnFace(3)->id<<" ";
  // 	  cout << endl; 
  //     }
  //   }    

  friend ostream &operator<<(ostream &out, Grid* const &a) {
    if (a == NULL) {out << " "; return out; }
    int_8 nCell = a->listCell.size(); //a->nCellSize = 0; 
    // for (auto i =0; i < a->listCell.size(); ++i) {
    //   if (a->listCell[i]->isAlive) {
    // 	nCell++; 
    // 	a->nCellSize += a->listCell[i]->node.size() + 1; 
    //   }
    // }
    out << "# vtk DataFile Version 2.0" << endl; 
    out << "Unstructure Grid" << endl; 
    out << "ASCII"<< endl; 
    out << endl << "DATASET UNSTRUCTURED_GRID"<< endl; 
    out << "POINTS " << a->listVertex.size() << " float" << endl; 
    for (auto p = a->vbegin(); p != a->vend(); ++p) {
      out << **p <<  endl ; 
    }
    out << endl << "CELLS "<< nCell << " " << nCell*5 << endl; 
    for (auto c = a->cbegin(); c != a->cend(); ++c) {
      if ((*c)->isAlive) out << *c ; 
   }
    out << endl << "CELL_TYPES " << nCell <<endl; 
    for (auto c = a->cbegin(); c != a->cend(); ++c) {
      if ((*c)->isAlive) out << (*c)->getType() << endl; 
    }
    out << endl << "POINT_DATA " << a->listVertex.size() << endl;
    for (auto i = 3; i < a->listVar.size(); ++i) {
      auto var = a->listVar[i]; 
      // cout << var->name << " " << var->loc << " "<< var->data.size() << " " << a->listCell.size() << endl; 
      // if (var->data.size() != a->listCell.size()) continue; 
      out << "SCALARS " << var->name << " float 1"<<endl; 
      out << "LOOKUP_TABLE default"<<endl; auto icnt = 0; 
      for (auto v: a->listVertex) {
	auto val = v->phi(var).eval(var); 
	out << ((abs(val) < 1e-10) ? 0 : val) << endl; 
      }
      //auto d = a->getPhiVertex(v);
      //for (auto i = 0; i<d.size(); ++i) out << d[i] << endl; 
      //      for (auto s = v->data.begin(); s != v->data.end(); ++s) {
      //	out << *s << " "; ++icnt; }
      out << endl; 
    }
    auto v = a->getVelVertex();
    out << "VECTORS vel float" << endl; 
    for (auto i = 0; i<v.size(); ++i) out << float(v[i][0]) << " " << float(v[i][1]) << " " << float(v[i][2]) << endl;

    // extra -remove
    // out << endl << "CELL_DATA " << nCell << endl;
    // for (auto j = 0; j < 3; ++j) {
    //   out << "SCALARS adapt"<<j<<" int 1"<<endl; 
    //   out << "LOOKUP_TABLE default"<<endl; auto icnt = 0; 
    //   for (auto i = 0; i < a->listCell.size(); ++i)
    // 	if (a->listCell[i]->isAlive) out << a->listCell[i]->adapt[j] << endl;
    // }
    // for (auto j = 0; j < 3; ++j) {
    //   out << "SCALARS level"<<j<<" int 1"<<endl; 
    //   out << "LOOKUP_TABLE default"<<endl; auto icnt = 0; 
    //   for (auto i = 0; i < a->listCell.size(); ++i)
    // 	if(a->listCell[i]->isAlive) out << a->listCell[i]->level[j] << endl;
    // } 
    // out << "SCALARS masterx int 1"<<endl; 
    // out << "LOOKUP_TABLE default"<<endl; 
    // for (auto i = 0; i < a->listCell.size(); ++i)
    //   if(a->listCell[i]->isAlive) out << a->listCell[i]->masterx[a->listCell[i]->level[0]] << endl;
    // out << "SCALARS mastery int 1"<<endl; 
    // out << "LOOKUP_TABLE default"<<endl; 
    // for (auto i = 0; i < a->listCell.size(); ++i)
    //   if(a->listCell[i]->isAlive) out << a->listCell[i]->mastery[a->listCell[i]->level[1]] << endl;


    out<<endl; 
    return out;
  }
 
  
  void writeFace(std::string a = "face.vtk") { 
    ofstream out; 
    out.open(a);
    out << "# vtk DataFile Version 2.0" << endl; 
    out << "Unstructure Grid" << endl; 
    out << "ASCII"<< endl; 
    out << endl << "DATASET UNSTRUCTURED_GRID"<< endl; 
    out << "POINTS " << listVertex.size() << " float" << endl; 
    for (auto p = vbegin(); p != vend(); ++p) {
      out << **p <<  endl ; 
    }
    out << endl << "CELLS "<< nFace << " " << 3*(nFace) << endl; 
    for (auto f: listFace) 
      if (f) out << f; 
  
    out << endl << "CELL_TYPES " << nFace <<endl; 
    for (auto f: listFace) 
      if (f) out << f->getType() << endl; 

    out << endl << "CELL_DATA " << nFace << endl; 
    for (auto i = 3; i < listVar.size(); ++i) {
      auto var = listVar[i]; 
      if (var->data.data.size() != nFace) continue; 
      out << "SCALARS " << var->name << " float 1"<<endl; 
      out << "LOOKUP_TABLE default"<<endl; 
      out << var->data << endl; 
      out << endl; 
    }
    
    // out << "SCALARS next" << " int 1"<<endl; 
   
    // for (auto f: listFace) 
    //   if (f) out << f->next << endl; 

    // out << "SCALARS prev" << " int 1"<<endl; 
    // out << "LOOKUP_TABLE default"<<endl; 
    // for (auto f: listFace) 
    //   if (f) out << f->prev << endl; 
  }

  vector<shared_ptr<Vertex > >::iterator vbegin() { return listVertex.begin(); }
  vector<shared_ptr<Vertex > >::iterator vend() {return listVertex.end(); }
  
  vector<shared_ptr<Cell > >::iterator  cbegin() {return listCell.begin(); }
  vector<shared_ptr<Cell > >::iterator  cend() {return listCell.end(); }

};


class Block2: public Grid {
public:
  Block2():Grid() {};
  Block2(initializer_list<double> n1, initializer_list<double > n2, int_4 nx, int_4 ny):Grid() {  
    Vec3 node1 = n1; 
    Vec3 node2 = n2; 
    Vec3 del = (node2 - node1);  
    meanD = min(del[0]/double(nx), del[1]/double(ny)); 
    addVertex({
	 {node1[0], node1[1], node1[2]}
	,{node2[0], node1[1], node1[2]} 
	,{node2[0], node2[1], node1[2]} 
	,{node1[0], node2[1], node1[2]}
      }); 
    addCell({0,1,2,3});
    (*listCell.rbegin())->convertToSimpleBlock({nx,ny}); 
    setCurrentLevels(); 
    makeFace(); 
    //setQuadBoundary(); 
    cout << "Block2: Cells: " << listCell.size(); 
    cout << " Faces: " << nFace << endl; 
    addVar({"u", "v", "w"});
  }
}; 


class Block3: public Grid {
public:
  Block3():Grid() {};
  Block3(initializer_list<double> n1, initializer_list<double > n2, int_4 nx, int_4 ny, int_4 nz):Grid() {    
    Vec3 node1 = n1; 
    Vec3 node2 = n2; 
    addVertex({ 
	 {node1[0], node1[1], node1[2]}
	,{node2[0], node1[1], node1[2]}
	,{node2[0], node2[1], node1[2]}
	,{node1[0], node2[1], node1[2]}
	,{node1[0], node1[1], node2[2]}
	,{node2[0], node1[1], node2[2]}
	,{node2[0], node2[1], node2[2]}
	,{node1[0], node2[1], node2[2]}
      });
    //addCell({0,1,2,3,4,5,6,7}); 
    //if (listCell.capacity() < (listCell.size() + nx*ny*nz)) {listCell.reserve(listCell.size() + nx*ny*nz);}
    //(*listCell.rbegin())->split({nx,ny, nz}); 
    //    listCell.rbegin()->splitHexa(nx,ny,nz); 
  }

}; 


// template <typename T>
// inline bool equals(const std::weak_ptr<T>& t, const std::weak_ptr<T>& u)
// {
//     return !t.owner_before(u) && !u.owner_before(t);
// }

// template <typename T>
// inline bool equals(const std::weak_ptr<T>& t, const std::shared_ptr<T>& u)
// {
//     return !t.owner_before(u) && !u.owner_before(t);
// }

#endif
