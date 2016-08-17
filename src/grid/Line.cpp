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
#include <girdap>

void Line::assignCelltoNode() {
  // check size first; if this is the first time assign only current nodes; 
  for (int_2 i = 0; i<2; ++i) {
    if (auto v = getVertex(i)) {  
      if ((*v)->cell.size() != 2) (*v)->reset(2); 
      (*v)->replaceCell(int(i+1)%2, id); 
    }
  }
}

void Line::convertToSimpleBlock(initializer_list<int> n, bool debug) {
  int n1=1; 
  // First find the existing vertex if non then proceed; 
  if (n.size() >= 1) {n1 = *(n.begin());} 
  else {  return; }
  if (debug) cout << "Converting cell: "<<id << " to block ("<< n1<< ")"<<endl; 

  vector<int_8 > ind; 
  ind.assign(n1+1, -1); 
  
  ind[0] = node[0];   
  ind[n1] = node[1];

  Vec3 dx, x0,x1;
  x0 = Vec3(**getVertex(0)); 
  dx = edge()/double(n1); 

  // Fill in non-existing Vertices
  for (auto i = 0; i < n1+1 ; ++i) {
    if (debug) cout << "Processing i = "<< i << " using index "<< ind[i] << endl; 
    if (ind[i]>=0) {
      if (debug)
	cout<< "Existing vertex: " << ind[i] << " ("<< *(*(grid->listVertex.begin()+ind[i])) <<")"<< endl; 
      //(*(grid->listVertex.begin()+ind[i][j]))->reset(4); 
      continue;
    }
    grid->addVertex(x0 + double(i)*dx); 
    ind[i] = grid->listVertex.size()-1;
    (*(grid->listVertex.rbegin()))->reset(4);
    if (debug) 	  
      cout<< "New vertex: " << ind[i] << " ("<< *(*(grid->listVertex.begin()+ind[i])) <<")"<< endl; 
    // boundary
    //      (*(grid->listVertex.rbegin()))->setCell(0, (j-1)*n1+i-1+nCell);   
  }

  auto nCell = grid->listCell.size()-1; 
  for (auto i = 0; i < n1+1; ++i) { 
    auto v = (*(grid->listVertex.begin()+ind[i]));
    int_8 i0 = i-1; i0 = (i0 == 0) ? id : i0+nCell; 
    int_8 i1 = i;   i1 = (i1 == 0) ? id : i1+nCell; 
    
    if (i > 0 && i < n1) {
      v->setCell({i0, i1}); 
    } else if (i == 0) {
      auto v0 = grid->listVertex.begin() + ind[i]; 
      v->setCell({(*v0)->cell[0], i1}); 
    } else if (i == n1) {
      auto v0 = grid->listVertex.begin() + ind[i]; 
      v->setCell({i0, (*v0)->cell[1]}); 
    }
  }

  for (auto i = 0; i < n1; ++i) {
    if (i == 0) {
      reset({ind[i], ind[i+1]}); 
    } else {
      grid->addCell({ind[i], ind[i+1]}); 
      auto x = (*grid->listCell.rbegin())->getCoord();	
      for (auto k = 0; k < grid->listVar.size(); ++k) {
	grid->listVar[k]->set(grid->listCell.size()-1,grid->listVar[k]->get(id)); // oldv[k] + (x-xcold)*oldg[k]); 
      }
    }
  }

  return;
}


void Line::refine(int dir) {
  dir =0; 
  if (adapt[0] < 1) return;
  if (level[0] == grid->levelHighBound[0]) return; 
  convertToSimpleBlock({2}); 
  
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





// Vec3 Line::grad(shared_ptr<Var > phi) {
//   auto c0 = (prev>=0) ? &(**(grid->listCell.begin() + prev)) : this; 
//   auto c1 = (next>=0) ? &(**(grid->listCell.begin() + next)) : this; 
//   Vec3 dx = 1./(c1->getCoord() - c0->getCoord());
//   if (next >= 0 && prev >= 0) {
//     return (phi->data[next] - phi->data[prev])*dx; 
//   } else {
//     double a=0, b=0; 
//     getBCPhi(phi, a, b); 
//     return (next >=0) ? (phi->data[next] - a*(phi->data[next] + b))*dx
//       : (a*(phi->data[prev] + b) - phi->data[prev])*dx; 
//   }
// }


// void Line::getBCPhi(shared_ptr<Var> phi, double &a, double &b) {
//   a = 0; b = 0; 
//   auto bndr = (prev >= 0) ? -next-1 : -prev-1; 
//   if (bndr < 0 || bndr >= phi->listBC.size() || !(phi->listBC[bndr])) { 
//     cout << "Something wrong with bndr : Value "<< bndr <<endl; 
//     return;
//   }
//   if (phi->listBC[bndr]->type == 0) {
//     a = phi->listBC[bndr]->a_val; 
//     b = phi->listBC[bndr]->b_val; 
//   } else if (phi->listBC[bndr]->type == 1) { 
//     auto row  = (prev >= 0) ? prev : next; 
//     auto c0 = *(grid->listCell.begin() + row); 
//     auto norm = vol(); norm = norm/norm.abs(); 
//     double dx = norm*(getCoord() - c0->getCoord()); 
//     a = (1.0 + dx*phi->listBC[bndr]->a_val);
//     b = dx*phi->listBC[bndr]->b_val; 
//   }
//   return;
// }
