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
#include "Grid.hpp" 
#include <iostream>

void Quad::assignCelltoNode() {
  // check size first; if this is the first time assign only current nodes; 
  for (int_2 i = 0; i<4; ++i) {
    if (auto v = getVertex(i)) {  
      if ((*v)->cell.size() != 4) (*v)->reset(4); 
      (*v)->replaceCell(int(i+2)%4, id); 
    }
  }
}

void Quad::convertToSimpleBlock(initializer_list<int> n, bool debug) {
  int n1=1; int n2=1;


  // First find the existing vertex if non then proceed; 
  if (n.size() >= 2) {n1 = *(n.begin()); n2 = *(n.begin()+1);}
  else if (n.size() == 1) {n1 = *(n.begin());}
  else {  return; }
  if (debug) cout << "Converting cell: "<<id << " to block ("<< n1<< ", "<<n2<<")"<<endl; 

  // auto xcold = getCoord(); 
  // vector<double> oldv(grid->listVar.size()); 
  // vector<Vec3> oldg(grid->listVar.size()); 
  // for (auto i = 0; i < grid->listVar.size(); ++i) {
  //   auto var = grid->listVar[i]; 
  //   oldv[i] = var->get(id); 
  //   oldg[i] = grad(var).eval(var); 
  // }  

  vector<vector<int_8 > > ind; 
  // vector<vector<vector<double > > > ovar; 
  // vector<int_8 > intv; 
  ind.resize(n1+1); 
  for (int_8 i=0; i < n1+1; ++i) ind[i].assign(n2+1,-1);
 
  // variables; 
  //ovar.resize(n1+1); 
  // for (int_8 i=0; i < n1+1; ++i) {
  //   ovar[i].resize(n2+1); 
  //   for (int_8 j = 0; j < n2+1; ++j) {
  //     ovar[i][j].assign(grid->listVar.size(), 0); 
  //   }
  // }
  
  ind[0][0] = node[0];   
  ind[n1][0] = node[1];
  ind[n1][n2] = node[2];
  ind[0][n2] = node[3];

  if (n1 >= 2) {
    ind[n1/2][0] = hangingVertexOnFace(0); 
    ind[n1/2][n2] = hangingVertexOnFace(2); 
  }
  if (n2 >= 2) {
    ind[n1][n2/2] = hangingVertexOnFace(1); 
    ind[0][n2/2] = hangingVertexOnFace(3); 
  }

  Vec3 dx1, dx2, dy1, dx, x0;
  dx1 = edge(0)/double(n1); 
  dx2 = edge(2)/double(n1); 
  dy1 = edge(3)/double(n2);     
  x0 = Vec3(**getVertex(0)); 

  // for (auto j = 0; j < n2+1 ; ++j) {
  //   for (auto i = 0; i < n1+1 ; ++i) {
  //     if (ind[i][j] < 0) continue;
  //     if (grid->listVertex[ind[i][j]]->xcoef.data.size() > 0) {
  // 	intv.emplace_back(ind[i][j]); 
  //     }
  //   }
  // }
  // cout << id <<  " : " ; 
  // for (auto j : intv) {
  //   cout << j << " " ; 
  // }
  // cout << endl; 

  // Fill in non-existing Vertices
  for (auto j = 0; j < n2+1 ; ++j) {
    for (auto i = 0; i < n1+1 ; ++i) {
      if (debug) cout << "Processing i = "<< i << " and j = " << j << " using index "<< ind[i][j] << endl; 
      if (ind[i][j]>=0) {
	if (debug)
	  cout<< "Existing vertex: " << ind[i][j] << " ("<< *(*(grid->listVertex.begin()+ind[i][j])) <<")"<< endl; 
	//(*(grid->listVertex.begin()+ind[i][j]))->reset(4); 
	continue;
      }
      dx = dx1 + (dx2 - dx1)*double(j)/double(n2);	
      grid->addVertex(x0 + double(i)*dx + double(j)*dy1);      
      ind[i][j] = grid->listVertex.size()-1;
      (*(grid->listVertex.rbegin()))->reset(4);
      if (debug) 	  
	cout<< "New vertex: " << ind[i][j] << " ("<< *(*(grid->listVertex.begin()+ind[i][j])) <<")"<< endl; 
       // boundary
      //      (*(grid->listVertex.rbegin()))->setCell(0, (j-1)*n1+i-1+nCell); 
    }
  }

  // for (auto j = 0; j < n2+1 ; ++j) {
  //   for (auto i = 0; i < n1+1 ; ++i) {
  //     int_8 i0 = -1; 
  //     for (auto k : intv) {
  // 	if (grid->listVertex[k]->isIn(*(grid->listVertex[ind[i][j]]))) {
  // 	  i0 = k; 
  // 	  break; 
  // 	}
  //     }
  //     if (i0 < 0) {
  // 	for (auto k = 0; k < grid->listVar.size(); ++k) { 
  // 	  ovar[i][j][k] = grid->listVar[k]->get(id); 
  // 	}
  //     } else {
  // 	for (auto k = 0; k < grid->listVar.size(); ++k) { 
  // 	  auto x = *(grid->listVertex[ind[i][j]]); 
  // 	  ovar[i][j][k] = grid->listVertex[i0]->evalPhi(grid->listVar[k], &x); 
  // 	}
  //     }
  //   }
  // }
  //cin.ignore().get(); 

  // 0 th cell; xhat(0-3) =  0.5*(0.5*(1+1/n1), 0.5*(1+1/n2))
  //                   + 0.5*(0.5*(1+1/n1), 0.5*(1-1/n2))
  // 1 th node; xhat = (0.5*(1+1/n1
  

  auto nCell = grid->listCell.size()-1; 
  for (auto j = 0; j < n2+1; ++j) {
    for (auto i = 0; i < n1+1; ++i) { 
      auto v = (*(grid->listVertex.begin()+ind[i][j]));
      int_8 i0 = (j-1)*n1 + i-1; i0 = (i0 == 0) ? id : i0+nCell; 
      int_8 i1 = (j-1)*n1 + i;   i1 = (i1 == 0) ? id : i1+nCell; 
      int_8 i2 = j*n1 + i;       i2 = (i2 == 0) ? id : i2+nCell; 
      int_8 i3 = j*n1 + i-1;     i3 = (i3 == 0) ? id : i3+nCell; 
      if (i > 0 && i < n1 && j > 0 && j < n2) {
	v->setCell({i0, i1, i2, i3}); 
      } else if (j == 0 && i > 0 && i < n1) {
	auto v0 = grid->listVertex.begin()+ind[i-1][j]; 
	v->setCell({(*v0)->cell[1], (*v0)->cell[1], i2, i3}); 
      } else if (j == n2 && i > 0 && i < n1) {
	auto v0 = grid->listVertex.begin() + ind[i-1][j]; 
	v->setCell({i0, i1, (*v0)->cell[2], (*v0)->cell[2]}); 
      } else if (i == 0 && j > 0 && j < n2) {
	auto v0 = grid->listVertex.begin() + ind[i][j-1]; 
	v->setCell({(*v0)->cell[3], i1, i2, (*v0)->cell[3]}); 
      } else if (i == n1 && j > 0 && j < n2) {
	auto v0 = grid->listVertex.begin() + ind[i][j-1]; 
	v->setCell({i0, (*v0)->cell[2],(*v0)->cell[2], i3});
      }
    }
  }

  // Patch for 2-1 or 1-2 splits with hanging node on opposite edges; 
  if (n1 == 2 && n2 == 1) {
    //find vertex on face 1 (east); 
    auto i = hangingVertexOnFace(1); 
    if (i >= 0) { 
      auto v = (*(grid->listVertex.begin()+i));
      v->cell[0] = grid->listCell.size(); 
      v->cell[3] = v->cell[0]; 
    }
  } else if (n1 == 1 && n2 == 2) {
    auto i = hangingVertexOnFace(2); 
    if (i >= 0) { 
      auto v = (*(grid->listVertex.begin()+i));
      v->cell[0] = grid->listCell.size(); 
      v->cell[1] = v->cell[0]; 
    }
  }

  for (auto j = 0; j < n2; ++j) {
    for (auto i = 0; i < n1; ++i) {
      if (i == 0 && j == 0) {
	reset({ind[i][j], ind[i+1][j], ind[i+1][j+1], ind[i][j+1]});	
	
	// auto x = getCoord(); 	
	// for (auto k = 0; k < grid->listVar.size(); ++k) {
	//   //auto val = 0.25*(ovar[i][j][k] + ovar[i+1][j][k] + ovar[i+1][j+1][k] + ovar[i][j+1][k]); 
	//   grid->listVar[k]->set(id, oldv[k] + (x-xcold)*oldg[k]); 	
	// }
      } else {
	grid->addCell({ind[i][j], ind[i+1][j], ind[i+1][j+1], ind[i][j+1]});
	auto x = (*grid->listCell.rbegin())->getCoord();	
	for (auto k = 0; k < grid->listVar.size(); ++k) {
	  //auto val = 0.25*(ovar[i][j][k] + ovar[i+1][j][k] + ovar[i+1][j+1][k] + ovar[i][j+1][k]); 
	  grid->listVar[k]->set(grid->listCell.size()-1,grid->listVar[k]->get(id)); // oldv[k] + (x-xcold)*oldg[k]); 
	}
	// for (auto var: grid->listVar) 
	//   var->set(grid->listCell.size()-1, var->get(id)); 
	
	//   }
      }	
    }
  }

  // for (auto j = 0; j < n2+1; ++j) {
  //   for (auto i = 0; i < n1+1; ++i) { 
  //     grid->listVertex[ind[i][j]]->setInterpCoef(); 
  //   }
  // }

  return;
}


void Quad::refine(int dir) {
  if (adapt[dir] < 1) return;
  if (level[dir] == grid->levelHighBound[dir]) return; 
  if (dir == 0) 
    convertToSimpleBlock({2, 1}); 
  else if (dir == 1)
    convertToSimpleBlock({1, 2});
  else 
    return; 
  
  adapt[dir]--; level[dir]++;
  (*(grid->listCell.rbegin()))->masterx.assign(masterx.begin(), masterx.end()); 
  (*(grid->listCell.rbegin()))->mastery.assign(mastery.begin(), mastery.end()); 
  (*(grid->listCell.rbegin()))->masterz.assign(masterz.begin(), masterz.end()); 

  if (dir == 0) masterx[level[dir]] = true; 
  if (dir == 1) mastery[level[dir]] = true; 
  if (dir == 2) masterz[level[dir]] = true; 

  (*(grid->listCell.rbegin()))->adapt.assign(adapt.begin(), adapt.end()); 
  (*(grid->listCell.rbegin()))->level.assign(level.begin(), level.end()); 
}

// bool Quad::coarsen() {
//   int check = 0; 
//   if (masterx[level[0]] && adapt[0] < 0) check = c+1; // x-direction not valid
//   if (mastery[level[1]] && adapt[1] < 0) check = c+2; 
//   if (c == 0) return false; 
//   else if (c == 1) { // x only
    
//   } else if (c == 2) { // y only
//   } else { //try both x and y

  
//   auto v = (*getVertex(3)); 
//   // following may not hold for complex geo grids. 
//   if (v->cell[1] < 0 || v->cell[3] < 0) { cout << "BNDR cell should not be a master" << endl; exit(1);}
//   auto c1 = *c->getCell(1); 
//   auto c2 = *c->getCell(2); 
//   auto c3 = *c->getCell(3); 
  
//   // x and y together (master); 
//   if (masterx[level[0]] && mastery[level[1]]) {

//     if (c1 && (c1->level[0] != level[0])) return false; // xlevels are the same
//     if (c1 && (c1->masterx[c1->level[0]])) return false; // not a master; 
//     if (c3 && (c3->level[1] != level[1])) return false; 
//     if (c3 && (c3->mastery[c3->level[1]])) return false; 


//   } 
//   if (masterx[level[0]]) {   // x only master



//   }
//   if (mastery[level[1]]) { // y only master
    
//   }

// }


bool Quad::coarsen(int dir) {
  if (adapt[dir] > -1) return false;
  if (level[dir] == grid->levelLowBound[dir]) return false; 
  if (dir == 0)
    if (!masterx[level[dir]]) return false; //cout << " should not see me -x" << endl; 
  if (dir == 1)
    if (!mastery[level[dir]]) return false; //cout << " should not see me -y" << endl; 
  // if (dir == 0) { 
  //   if (!masterx[level[0]]) return false;  
  // }
  // if (dir == 1) {
  //   if (!mastery[level[1]]) return false; 
  // }
  // Check +x and +y direction//
  auto i0 = 1; 
  if (dir == 1) i0 = 2; 
  if (node[i0] < 0 || node[i0+1] < 0) return false; 
  vector<shared_ptr<Vertex> > v(node.size()); 
  for (auto i = 0; i < node.size(); ++i) v[i] = (*getVertex(i));
  //--------------------------------------------------------------- 
  // QUESTION: EU: Why do following block coarsening? Not clear----
  // Check if the south cell is bigger
  if (v[(i0+3)%4]->cell[(i0+3)%4] >= 0 && (v[(i0+3)%4]->cell[(i0+3)%4] == v[(i0+3)%4]->cell[i0])) return false; 
  // check if north cell is bigger; 
  if (v[(i0+2)%4]->cell[(i0+1)%4] >= 0 && (v[(i0+2)%4]->cell[(i0+1)%4] == v[(i0+2)%4]->cell[(i0+2)%4])) return false; 
  //----------------------------------------------------------------


  // BOTH c0 and c1 points to east/north cells; 
  auto c0 = v[i0]->getCell((i0+1)%4); 
  auto c1 = v[(i0+1)%4]->getCell(i0); 
  
  if (c0 && ((*c0)->id == (*c1)->id)) {                             // if same level or larger
    if ((*c0)->adapt[dir] >= 0) return false;                       // but not marked to be coarsened
    if (level[(dir+1)%2] != (*c0)->level[(dir+1)%2]) return false;  // if not same level in other direction
    if (level[dir] != (*c0)->level[dir]) return false;              // if not same level in the same direction

    // FROM this point there is only c0 (=c1) and compatible for merger
    vector<shared_ptr<Vertex> > nv(node.size()); 
    for (auto i = 0; i < (*c0)->node.size(); ++i) nv[i] = *(*c0)->getVertex(i); 
    
    // If next cell has a bigger cell in south
    if (nv[i0]->cell[(i0+3)%4]>=0 && (nv[i0]->cell[(i0+3)%4] == nv[i0]->cell[i0])) return false; 
    // if next cell has a bigger cell in north
    if (nv[(i0+1)%4]->cell[(i0+1)%4]>=0 && (nv[(i0+1)%4]->cell[(i0+1)%4] == nv[(i0+1)%4]->cell[(i0+2)%4])) return false; 
    
    // if south has smaller face (higher level)
    if ((*c0)->hangingVertexOnFace((i0+3)%4) >= 0) return false; 
    // if north has smaller faces (higher level)
    if ((*c0)->hangingVertexOnFace((i0+1)%4) >= 0) return false; 

    // 
    if (hangingVertexOnFace((i0+3)%4) >=0 ) return false; 
    if (hangingVertexOnFace((i0+1)%4) >=0 ) return false; 

    // Then coarsen; 
    
    //[0] Calculate new value of old variables first (no point after other cells removal); 
    for (auto var : grid->listVar) { 
      var->set(id, 0.5*(var->get(id) + var->get((*c0)->id))); 
    }    
    
    // [1] Correct indices;
    //     [a] Hanging nodes; 
    auto hv = (*c0)->hangingVertexOnFace(i0); 
    if (hv >= 0) {
      grid->listVertex[hv]->cell[(i0+2)%4] = id; 
      grid->listVertex[hv]->cell[(i0+3)%4] = id; 
    }

    //     [b] vertices of ngbr cells; 
    for (auto i = 0; i < nv.size(); ++i) nv[i]->cell[(i+2)%4] = id;
    
    //     [c] remove other cell; 
    (*c0)->isAlive = false; 
    (*c0)->node.assign(4, -1); 
    
    //     [d] this cell; 
    node[i0] = nv[i0]->id; 
    node[(i0+1)%4] = nv[(i0+1)%4]->id; 

    if (dir == 0) masterx[level[dir]] = false; 
    if (dir == 1) mastery[level[dir]] = false; 

    adapt[dir] = 0; 
    level[dir] -= 1; 
  }
  return true; 
}  
    


bool Quad::refine() {
  if (adapt[0] > 0 && adapt[1] > 0) {
    convertToSimpleBlock({2, 2}); 
    adapt[0]--; adapt[1]--; 
    level[0]++; level[1]++;     
    for (auto i=0; i<3; ++i) {
      (*(grid->listCell.rbegin()+i))->adapt.assign(adapt.begin(), adapt.end()); 
      (*(grid->listCell.rbegin()+i))->level.assign(level.begin(), level.end()); 
    }
  } else if (adapt[0] > 0) {
    convertToSimpleBlock({2, 1}); 
    adapt[0]--; level[0]++;
    for (auto i=0; i<1; ++i) {
      (*(grid->listCell.rbegin()+i))->adapt.assign(adapt.begin(), adapt.end()); 
      (*(grid->listCell.rbegin()+i))->level.assign(level.begin(), level.end()); 
    } 
  } else if (adapt[1] > 0) {
    convertToSimpleBlock({1, 2}); 
    adapt[1]--; level[1]++;
    for (auto i=0; i<1; ++i) {
      (*(grid->listCell.rbegin()+i))->adapt.assign(adapt.begin(), adapt.end()); 
      (*(grid->listCell.rbegin()+i))->level.assign(level.begin(), level.end()); 
    } 
  }
  return max(adapt[0], adapt[1]); 
}

vector<int_8> Quad::ngbrCellList() { 
  vector<int_8> nlist; 
  for (auto i = 0; i < node.size(); ++i) { 
    auto v = *getVertex(i); 
    auto j = (i+1)%node.size(); //forward
    if (v->cell[j] >= 0 && (std::find(nlist.begin(), nlist.end(), v->cell[j])!=nlist.end())) 
      nlist.emplace_back(v->cell[j]); 
    j = (i+3)%node.size(); //backward
    if (v->cell[j] >= 0 && (std::find(nlist.begin(), nlist.end(), v->cell[j])!=nlist.end())) 
      nlist.emplace_back(v->cell[j]); 
  }
  return nlist; 
}

// void Quad::checkCoarseLevel(int dir) {
//   int_8 f1, f2, b1, b2; 
//   if (dir == 1) {
//     f1 = (*getVertex(1))->cell[2]; 
//     f2 = (*getVertex(2))->cell[1]; 
//     b1 = (*getVertex(3))->cell[0]; 
//     b2 = (*getVertex(0))->cell[3]; 
//   } else if (dir == 0) {
//     f1 = (*getVertex(2))->cell[3]; 
//     f2 = (*getVertex(3))->cell[2]; 
//     b1 = (*getVertex(0))->cell[1]; 
//     b2 = (*getVertex(1))->cell[0]; 
//   } else {
//     return; 
//   }
// }

void Quad::checkIslandLevels(int dir) {
  auto myFinalLevel = adapt[dir] + level[dir]; 
  vector<int_8> ngbr(4); //, f2, b1, b2; 
  if (dir == 1) {
    ngbr = {(*getVertex(1))->cell[2], (*getVertex(2))->cell[1]  
	    , (*getVertex(3))->cell[0], (*getVertex(0))->cell[3]}; 
  } else if (dir == 0) {
    ngbr = {(*getVertex(2))->cell[3], (*getVertex(3))->cell[2]  
	    , (*getVertex(0))->cell[1], (*getVertex(1))->cell[0] }; 
  } else {
    return; 
  }
  // island cells should follow neighbors; 
  for (auto i = 0; i<2; ++i) {
    for (auto j = 2; j<4; ++j) {
      if (ngbr[i] >= 0 && ngbr[j] >=0 ) { //not a boundary cell
	auto l0 = grid->listCell[ngbr[i]]->level[dir]+grid->listCell[ngbr[i]]->adapt[dir];  
	auto l2 = grid->listCell[ngbr[j]]->level[dir]+grid->listCell[ngbr[j]]->adapt[dir];  
	if (l0 == l2) {
	  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  // !!!!!!!!!! Change 0 to -1 when coarse levels exist; 
	  adapt[dir] = max(-1, min(1, l0 - level[dir])); 
	  return; 
	}
      }
    }
  }
}


void Quad::checkNgbrLevel(int dir) {
  auto myFinalLevel = adapt[dir] + level[dir]; 
  vector<int_8> ngbr(8); //, f2, b1, b2; 
  //  if (dir == 1) {
    ngbr = {(*getVertex(1))->cell[2], (*getVertex(2))->cell[1]  
	    , (*getVertex(3))->cell[0], (*getVertex(0))->cell[3] 
	    , (*getVertex(2))->cell[3], (*getVertex(3))->cell[2] 
	    , (*getVertex(0))->cell[1], (*getVertex(1))->cell[0]}; 
  //  } else if (dir == 0) {
  //   ngbr = {(*getVertex(2))->cell[3], (*getVertex(3))->cell[2]  
  // 	    , (*getVertex(0))->cell[1], (*getVertex(1))->cell[0] }; 
  // } else {
  //   return; 
  // }
  // adjacent cells should not be simultaneously refined and coarsened
  short ref = 0; 
  for (auto i = 0; i<ngbr.size(); ++i) {
    if (ngbr[i] >= 0) {
      if (grid->listCell[ngbr[i]]->adapt[dir] * adapt[dir] < 0) {
	grid->listCell[ngbr[i]]->adapt[dir]  = max(short(0), grid->listCell[ngbr[i]]->adapt[dir]); 
	adapt[dir] = max(short(0), adapt[dir]); 
      }	
    }
  }
  short ngbrFinalLevel; 
  for (auto i = 0; i < ngbr.size(); ++i) {
    if (ngbr[i] >= 0) {
      ngbrFinalLevel = grid->listCell[ngbr[i]]->adapt[dir] + grid->listCell[ngbr[i]]->level[dir]; 
      if (myFinalLevel - ngbrFinalLevel > 1) {
	grid->listCell[ngbr[i]]->adapt[dir] = myFinalLevel - grid->listCell[ngbr[i]]->level[dir] - 1;
	grid->listCell[ngbr[i]]->checkNgbrLevel(dir); 
      } else if (ngbrFinalLevel - myFinalLevel > 1) {
	adapt[dir] = ngbrFinalLevel - level[dir] - 1; 
      }
    }
  }
  return; 
}

// void checkNgbrLevel(int dir) { 
//   auto ngbrList = ngbrCellList(); 
//   for (auto n : ngbrList) {
//     auto c = grid->listCell[n]; 
//     if (c) { 
//       auto cdiff = level[dir] - c->level[dir]; 
//       auto ndiff = cdiff + adapt[dir] - c->adapt[dir];
//       if (ndiff > 1) {
// 	c->adapt[dir] = level[dir]+adapt[dir]-1; 
// 	//c->checkNgbrLevel();
//       } else if (ndiff < -1) {
// 	adapt[dir] = c->level[dir]+c->adapt[dir]-1;
       
//       }
      
      

// vector<share_ptr<Cell > > getFaces() {
//   vector<share_ptr<Cell > > l; 
//   auto v0 = getVertex(0); 
//   if (auto f = (*v0)->face[0]) {
//     l.emplace_back((*v0)->face[0]); 
//     if ((*v0)->face[0]->node[1] != node[1]) {
//       v0 = v0->ngbr(1); 
//       l.emplace_back((*v0)->face[0]); 
//     }
//     // check 
//     cout << (*v0)->face[0]->node[1] << " should be equal to " << node[1] << endl; 
//   } else { 
//     cout << "there is no face!!! error??? " << endl; 
//     exit(1); 
//   }
// }


// double div() {
//   for (auto i = 0; i<node.size(); ++i) {
//     if (auto v = getVertex(i)) {
//       auto f = v->face[0]; 
      
//       j = (i+1)%4; 
//       if ((*v)->cell[j] >= 0) {
// 	fVel = 0.5*(thisVar->data[id] + thisVar->data[(*v)->cell[j]]); 
	

   
// }




