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

void Vertex::cellResize(int_2 size)  {cell.resize(size);} 
void Vertex::setCell(int_2 ind, int_8 value) {
  if (cell[ind] < 0) cell[ind] = value; 
}
void Vertex::replaceCell(int_2 ind, int_8 value) { 
  cell[ind] = value;
}

shared_ptr<Cell >  *Vertex::getCell(int_2 ind, bool debug) {
  if (debug) cout << " getCell: " << ind << " " << cell[ind] << " "; 
  if (ind >= 0 && ind < cell.size() && cell[ind] >= 0) {
    if (debug) cout << cell[ind] << endl; 
    return &(*(grid->cbegin() + cell[ind]));
  } else { 
    if (debug) cout << "NULL"<<endl;
    return NULL;
  } 
}

shared_ptr<Vertex > *Vertex::ngbr(int_2 d)  {
  if (cell.size() == 2) {
    if (d == 1) {
      if (cell[1] > 0) return (*getCell(1))->getVertex(1); 
    } else if (d == -1) {
      if (cell[0] > 0) return (*getCell(0))->getVertex(0); 
    } 
    return NULL; 
  }
  if (cell.size() == 4) {
    int i; int j; int i0; int j0;  
    if (d == 1) { i = 1; j = (i+1)%4; i0 = (i+3)%4; j0 = (j+1)%4;}
    else if (d == -1) {i = 0; j = (i+3)%4; i0 = (i+1)%4; j0 = (j+3)%4;}
    else if (d == 2) {i = 2; j = (i+1)%4; i0 = (i+3)%4; j0 = (j+1)%4;}
    else if (d == -2) {i = 1; j = (i+3)%4; i0 = (i+1)%4; j0 = (j+3)%4;}
    else return NULL; 
    if (cell[i] >= 0 && cell[j] >= 0) {
      if (cell[i] == cell[j]) return NULL; 
      auto v1 = (*getCell(i))->getVertex(j); //cell[i]->node[j].lock(); 
      auto v2 = (*getCell(j))->getVertex(i); //cell[j]->node[i].lock(); 
      if (v1 && v2 && *v1 && *v2) {
	if (v1 == v2) {
	  return v1; 
	} else {
	  if (cell[i0] == cell[i]) return v2; 
	  if (cell[j] == cell[j0]) return v1;
	  if (((*v2)->cell[i]>=0) && ((*v2)->cell[i0]>=0) && ((*v2)->cell[i] == (*v2)->cell[i0])) return v2; 
	  if (((*v1)->cell[j]>=0) && ((*v1)->cell[j0]>=0) && (*v1)->cell[j] == (*v1)->cell[j0]) return v1;	
	  cout << "d: " << d << " i: "<< i << " j: " << " i0: " << i0 << " j0: " << j0 << endl;   
	  cout << "ne v1 ne de v2"<< endl;
	  cout << cell[i0] << " " << cell[i] << endl; 
	  cout << cell[j] << " " << cell[j0] << endl;
	  cout << (*v2)->cell[i] << " " << (*v2)->cell[i0] << endl; 
	  cout << (*v1)->cell[j] << " " << (*v1)->cell[j0] << endl; 
	  exit(1); 
	}
      } else { 
	cout << "d: " << d << " i: "<< i << " j: " << " i0: " << i0 << " j0: " << j0 << endl; 
	cout << "Vertex:e()! Vertex of the cell doesn't exist. Something went wrong!"<< endl;
	exit(1); 
      }
    } else if (cell[i] >= 0) {
      return (*getCell(i))->getVertex(j); //cell[i]->node[j].lock(); 
    } else if (cell[j] >= 0) { 
      return (*getCell(j))->getVertex(i); //cell[j]->node[i].lock(); 
    }
    return NULL; 
  } else {
    return NULL;
  } 
}




Scheme<double> Vertex::phi(vector<shared_ptr<Boundary> > const &bc, int_2 bias) { 
  Scheme<double> sch; double sum = 0;
  vector<int_8> flag; 
  for (auto i = 0; i < cell.size(); ++i) {
    if (cell[i] >= 0) {
      if (std::find(flag.begin(), flag.end(), cell[i])!=flag.end()) continue; 
      flag.push_back(cell[i]); 
      double w = (*getCell(i))->vol().abs(); // /((*getCell(i))->getCoord() - *this).abs(); 
      sch.push_pair(cell[i], w);
      sum += w; 
    } else { 
      auto bndr = -cell[i]-1; 
      if (bc[bndr]->type == 0) {
	sch.ind.clear(); 
	sch.val.clear(); 
	sch.c = 0; 
	sch.push_constant(bc[bndr]->b_val); 
	return sch; 
      }
    }
  }
  flag.clear(); 
  sch /= sum; 
  return sch; 
}

Scheme<double> Vertex::phi(shared_ptr<Var> var, int_2 bias) { 
  return phi(var->listBC, bias);
}


