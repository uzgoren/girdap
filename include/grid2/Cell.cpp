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

shared_ptr<Vertex> *Cell::getVertex(int_2 i) {
  if (i < 0) return NULL;
  if (i > int(node.size())-1) return NULL; 
  return &(*(grid->vbegin()+node[i])); 
}

double Cell::phiVal(VecX<double> &phi, int_2 bias) {
  if (next < 0 && prev < 0) {
    return phi[id]; 
  } else {
    if (next >= 0 && prev >= 0) {
      return phiVal_iface(phi, bias); 
    } else if (next >= 0) {
      return phi[next]; 
    } else if (prev >= 0) {
      return phi[prev]; 
    }
  }
  return 0;
}

double Cell::phiVal(shared_ptr<Var> &var, int_2 bias) {
  return phiVal(var->data, var->listBC, bias); 
}

double Cell::phiVal(VecX<double> &phi, vector<shared_ptr<Boundary> > const &bc, int_2 bias) {
  if (next < 0 && prev < 0) {
    return phi[id]; 
  } else {
    if (next >= 0 && prev >= 0) {
      return phiVal_iface(phi, bias); 
    } else {
      return phiVal_bcface(phi, bc, bias); 
    }
  }
  return 0; 
}

double Cell::phiVal_iface(VecX<double> &phi, int_2 const &bias) {
  if (next < 0 || prev < 0) {
    cout << "phiVal_iface: should be used on internal faces!" << endl; 
    exit(1);
  }
  if (bias == 0) {
    auto area = vol();
    auto xf = getCoord(); 
    auto xp = grid->listCell[prev]->getCoord(); 
    auto xn = grid->listCell[next]->getCoord(); 
    double dn = abs( (xn-xf)*area ); 
    double dp = abs( (xp-xf)*area ); 
    return (phi[prev]*dn + phi[next]*dp)/(dn+dp); 
  } else if (bias < 0) {
    return phi[prev]; 
  } else {
    return phi[next];  
  }
}

double Cell::phiVal_bcface(VecX<double> &phi, vector<shared_ptr<Boundary> > const &bc, int_2 const &bias) {
  auto bndr = (prev >= 0) ? -next-1 : -prev-1; 
  auto row = (prev >= 0) ? prev : next; 
  if (bc[bndr]->type == 0) { 
    return bc[bndr]->b_val + phi[row]*bc[bndr]->a_val; 
  } else if (bc[bndr]->type == 1) {
    auto c0 = grid->listCell[row]; 
    auto norm = vol(); norm = norm/norm.abs(); 
    double dx = norm*(getCoord() - c0->getCoord()); 
    return dx*bc[bndr]->b_val + phi[row]*(1.0 + dx*bc[bndr]->a_val); 
  } else { 
    cout << "Cell (Quad) :: interp :: boundary condition type not recognized " << endl; 
    exit(1); 
  }
}






Scheme<double> Cell::phi(vector<shared_ptr<Boundary> > const &bc, int_2 bias) {
  Scheme<double> sch; 
  if (next < 0 && prev < 0) {
    sch.push_pair(id, 1.0); 
  } else {
    // use next and prev to compute phi;
    if (next >= 0 && prev >= 0) {
      if (bias == 0) {
	double dn = abs((grid->listCell[next]->getCoord() - getCoord())*vol()); 
	double dp = abs((grid->listCell[prev]->getCoord() - getCoord())*vol()); 
	sch.push_pair(prev, dn/(dn+dp)); 
	sch.push_pair(next, dp/(dn+dp)); 
      } else if (bias == -1) { 
	sch.push_pair(prev, 1.0); 
      } else if (bias == 1) { 
	sch.push_pair(next, 1.0); 
      }
    } else {
      auto bndr = (prev >= 0) ? -next-1 : -prev-1; 
      auto row = (prev >= 0) ? prev : next; 
      if (bc[bndr]->type == 0) { 
	sch.push_constant(bc[bndr]->b_val);
	sch.push_pair(row, bc[bndr]->a_val); 
      } else if (bc[bndr]->type == 1) {
	auto c0 = grid->listCell[row]; 
	auto norm = vol(); norm = norm/norm.abs(); 
	double dx = norm*(getCoord() - c0->getCoord()); 
	sch.push_constant(dx*bc[bndr]->b_val);
	sch.push_pair(row, (1.0 + dx*bc[bndr]->a_val)); 
      } else { 
	cout << "Cell (Quad) :: interp :: boundary condition type not recognized " << endl; 
	exit(1); 
      }
 
    }
  }
  return sch; 
}

Scheme<double> Cell::phi(shared_ptr<Var> var, int_2 bias) {
  return phi(var->listBC, bias); 
}


Scheme<Vec3> Cell::grad(vector<shared_ptr<Boundary> > const &bc) {
  Scheme<Vec3> sch; 
  if (next < 0 && prev < 0) {
    // control volume of the cell;
    double cvol = vol().abs();
    for (auto i = 0; i < face.size(); ++i) { 
      auto f = grid->listFace[face[i]]; 
      auto area = (f->prev == id) ? f->vol() : -f->vol(); 
      Scheme<double> tmp = f->phi(bc); 
      for (auto j =0; j<tmp.ind.size(); ++j) {
	sch.push_pair(tmp.ind[j], tmp.val[j]*area/cvol); 
      }      
    }
  } else {
    // use next and prev to compute phi;
    if (next >= 0 && prev >= 0) {
      if (grid->listCell[next]->level[orient] == grid->listCell[prev]->level[orient]) {
      	auto norm = vol(); 
      	auto area = norm.abs(); 
       	norm = norm/area;
       	auto dx = grid->listCell[next]->getCoord() - grid->listCell[prev]->getCoord(); 
	auto onebydx = norm/(dx*norm); 
       	sch.push_pair(next, onebydx);
       	sch.push_pair(prev, -onebydx); 
      }	else { // this part is cell specific // 2d-3d
	vector<Vec3> v; 
	vector<Scheme<double>> tmp; 
      	auto norm = vol(); 
       	norm = norm/norm.abs();
	v.push_back(grid->listCell[prev]->getCoord()); 
	v.push_back(*(grid->listVertex[node[0]])); 
	v.push_back(grid->listCell[next]->getCoord()); 
	v.push_back(*(grid->listVertex[node[1]])); 
	
	tmp.push_back(grid->listCell[prev]->phi(bc));
	tmp.push_back(grid->listVertex[node[0]]->phi(bc));
	tmp.push_back(grid->listCell[next]->phi(bc));
	tmp.push_back(grid->listVertex[node[1]]->phi(bc)); 
	
	tmp[3] += grid->listCell[prev]->phi(bc);
	tmp[0] += grid->listVertex[node[0]]->phi(bc);
	tmp[1] += grid->listCell[next]->phi(bc);
	tmp[2] += grid->listVertex[node[1]]->phi(bc); 
	
	auto vol = 0.5*((v[3]-v[1])^(v[2]-v[0])).abs(); 
	
	for (auto j = 0; j < 4; ++j) {
	  auto del = v[(j+1)%4] - v[j]; 
	  auto area = Vec3(-del[1], del[0], 0); 
	  for (auto i = 0; i < tmp[j].size(); ++i)
	    sch.push_pair(tmp[j].ind[i], (0.5*tmp[j].val[i]/vol)*area); //0.5 from average;
	}	
       }      
    } else {
      //      cout << "**** " << vol() << " " << type << endl; 
      //      cout << *(grid->listVertex[node[0]]) << " " << *(grid->listVertex[node[1]]) << endl; 
      auto bndr = (prev >= 0) ? -next-1 : -prev-1; 
      auto row = (prev >= 0) ? prev : next; 
      auto norm = vol(); 
      auto area = norm.abs(); 
      norm = norm/area;
      auto dx = (next >= 0) ? grid->listCell[next]->getCoord() - getCoord() : getCoord() - grid->listCell[prev]->getCoord(); 
      auto onebydx = norm/(dx*norm);

      if (bc[bndr]->type == 0) {
	if (row == next) onebydx = -onebydx; 
      	sch.push_constant(bc[bndr]->b_val * onebydx);
      	sch.push_pair(row, (bc[bndr]->a_val - 1.0) * onebydx); 
      } else if (bc[bndr]->type == 1) {
      	sch.push_constant((bc[bndr]->b_val)*norm);
      	sch.push_pair(row, (bc[bndr]->a_val)*norm); 	
      } else { 
      	cout << "Cell (Line) :: grad :: boundary condition type not recognized " << endl; 
      	exit(1); 
      }
    }
  }
  return sch; 
}

Scheme<Vec3> Cell::grad(shared_ptr<Var> var) {
  return grad(var->listBC); 
}



Scheme<double> Cell::normGrad(vector<shared_ptr<Boundary> > const &bc) {
  Scheme<double> sch; 
  if (next < 0 && prev < 0) {
    // control volume of the cell;
    cout << "You can only call form normGrad for a face!!! " << endl; 
    exit(1); 
  } else {
    // use next and prev to compute phi;
    if (next >= 0 && prev >= 0) {
      if (grid->listCell[next]->level[orient] == grid->listCell[prev]->level[orient]) {
      	auto norm = vol(); 
      	auto area = norm.abs(); 
       	norm = norm/area;
       	auto dx = grid->listCell[next]->getCoord() - grid->listCell[prev]->getCoord(); 
	auto onebydx = 1.0/(dx*norm); 
       	sch.push_pair(next, onebydx);
       	sch.push_pair(prev, -onebydx); 
      }	else { // this part is cell specific // 2d-3d
	vector<Vec3> v; 
	vector<Scheme<double>> tmp; 
      	auto norm = vol(); 
       	norm = norm/norm.abs();
	v.push_back(grid->listCell[prev]->getCoord()); 
	v.push_back(*(grid->listVertex[node[0]])); 
	v.push_back(grid->listCell[next]->getCoord()); 
	v.push_back(*(grid->listVertex[node[1]])); 
	
	tmp.push_back(grid->listCell[prev]->phi(bc));
	tmp.push_back(grid->listVertex[node[0]]->phi(bc));
	tmp.push_back(grid->listCell[next]->phi(bc));
	tmp.push_back(grid->listVertex[node[1]]->phi(bc)); 
	
	tmp[3] += grid->listCell[prev]->phi(bc);
	tmp[0] += grid->listVertex[node[0]]->phi(bc);
	tmp[1] += grid->listCell[next]->phi(bc);
	tmp[2] += grid->listVertex[node[1]]->phi(bc); 
	
	auto vol = 0.5*((v[3]-v[1])^(v[2]-v[0])).abs(); 
	
	for (auto j = 0; j < 4; ++j) {
	  auto del = v[(j+1)%4] - v[j]; 
	  auto area = Vec3(-del[1], del[0], 0); 
	  for (auto i = 0; i < tmp[j].size(); ++i)
	    sch.push_pair(tmp[j].ind[i], (0.5*tmp[j].val[i]/vol)*area*norm); //0.5 from average;
	}	
       }      
    } else {
      //      cout << "**** " << vol() << " " << type << endl; 
      //      cout << *(grid->listVertex[node[0]]) << " " << *(grid->listVertex[node[1]]) << endl; 
      auto bndr = (prev >= 0) ? -next-1 : -prev-1; 
      auto row = (prev >= 0) ? prev : next; 
      auto norm = vol(); 
      auto area = norm.abs(); 
      norm = norm/area;
      auto dx = (next >= 0) ? grid->listCell[next]->getCoord() - getCoord() : getCoord() - grid->listCell[prev]->getCoord(); 
      auto onebydx = 1.0/(dx*norm);

      if (bc[bndr]->type == 0) {
	if (row == next) onebydx = -onebydx; 
      	sch.push_constant(bc[bndr]->b_val * onebydx);
      	sch.push_pair(row, (bc[bndr]->a_val - 1.0) * onebydx); 
      } else if (bc[bndr]->type == 1) {
      	sch.push_constant((bc[bndr]->b_val));
      	sch.push_pair(row, (bc[bndr]->a_val)); 	
      } else { 
      	cout << "Cell (Line) :: grad :: boundary condition type not recognized " << endl; 
      	exit(1); 
      }
    }
  }
  return sch; 
}


Scheme<double> Cell::normGrad(shared_ptr<Var> var) {
  return normGrad(var->listBC); 
}






// double Cell:grad_iface(VecX<double> &phi) {
//   if (next >= 0 && prev >= 0) {
//     auto n = grid->listCell[next]; 
//     auto p = grid->listCell[prev]; 
    
//   }

// }
