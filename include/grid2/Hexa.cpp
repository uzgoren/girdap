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

void Hexa::split(initializer_list<int> n) {
  int n1 = 1; int n2 = 1; int n3 = 1; 
  if (n.size() >= 3) {n1 = *(n.begin()); n2 = *(n.begin()+1); n3 = *(n.begin()+2); }
  else if (n.size() == 2) {n1 = *(n.begin()); n2 = *(n.begin()+1);}
  else if (n.size() == 1) {n1 = *(n.begin());}
  else { return; }

  Vec3 dx1, dx2, dx3, dx4, dy1, dy2, dz, dx, dy, dxa, dxb; 
  vector<int_8 > kk, ll(8,-1);
  for (auto v :node) {
    kk.push_back(v->id); 
  }

  dx1 = edge(0)/double(n1); 
  dx2 = edge(2)/double(n1); 
  dx3 = edge(4)/double(n1); 
  dx4 = edge(6)/double(n1); 
  dy1 = edge(3)/double(n2); 
  dy2 = edge(7)/double(n2); 
  dz  = edge(8)/double(n3); 
  
  int_8 ind[n1+1][n2+1][n3+1];
  for (auto k = 0; k < n3+1 ; ++k) {
    for (auto j = 0; j < n2+1 ; ++j) {
      for (auto i = 0; i < n1+1 ; ++i) {
	ind[i][j][k] = -1; 
	if (i == 0 && j == 0 && k == 0) {
	  ind[i][j][k] = kk[0]; 
	} else if (i == 0 && j == n2 && k == 0) {
	  ind[i][j][k] = kk[3]; 
	} else if (i == n1 && j == n2 && k == 0) {
	  ind[i][j][k] = kk[2]; 
	} else if (i == n1 && j == 0 && k == 0) {
	  ind[i][j][k] = kk[1];
	} else if (i == 0 && j == 0 && k == n3) {
	  ind[i][j][k] = kk[4]; 
	} else if (i == 0 && j == n2 && k == n3) {
	  ind[i][j][k] = kk[7]; 
	} else if (i == n1 && j == n2 && k == n3) {
	  ind[i][j][k] = kk[6]; 
	} else if (i == n1 && j == 0 && k == n3) {
	  ind[i][j][k] = kk[5];  
	} else {
	  dy   = dy1 + (dy2 - dy1)*double(k)/double(n3); 
	  dxa  = dx1 + (dx3 - dx1)*double(k)/double(n3); 
	  dxb  = dx2 + (dx4 - dx2)*double(k)/double(n3); 
	  dx = dxa + (dxb - dxa)*double(j)/double(n2); 
	  grid->addVertex(*(node.at(0)) + double(i)*dx + double(j)*dy +double(k)*dz);
	  (*grid->vbegin())->grid = grid; 
	  ind[i][j][k] = (grid->listVertex.size())-1; 
	}
      }
    }
  }

  for (auto k = 0; k < n3; ++k) {
    for (auto j = 0; j < n2; ++j) {
      for (auto i = 0; i < n1; ++i) {
	ll[0] = ind[i][j][k]; 
	ll[1] = ind[i+1][j][k]; 
	ll[2] = ind[i+1][j+1][k]; 
	ll[3] = ind[i][j+1][k]; 
	ll[4] = ind[i][j][k+1]; 
	ll[5] = ind[i+1][j][k+1]; 
	ll[6] = ind[i+1][j+1][k+1]; 
	ll[7] = ind[i][j+1][k+1]; 
	if (i == 0 && j == 0 && k == 0) {
	  int_4 icnt = 0; 
	  for (auto v = node.begin(); v !=node.end(); ++v) {
	    (*v) = grid->listVertex.at(ll[icnt]); ++icnt; 
	  }
	  initNodeList(*(grid->listCell.begin() + id));
	} else {	
	  grid->addCell(ll);
	}
      }
    }
  }
  
}
