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


shared_ptr<Cell> Vertex::getFace(int_2 order) {
  if (cell.size() == 4 && order >= 4) {
    cout << "Error, getFace: face order can be 0,1,2,3 but you asked for "<< order<<endl; 
    exit(1); 
  }
  int_8 c0, c1; 
  if (order < 4) {
    if (order < 2) { // To make it compatible with next prev notation 
      c0 = order; c1 = (order+1)%4;
    } else {
      c0 = (order+1)%4; c1 = order; 
    }
  } else if (order < 8) {
    if (order < 6) {
      c0 = order; c1 = order + 1; 
    } else {
      c0 = order + 1; c1 = order;  
      if (c0 == 8) c0 = 4;
    } 
  } else {
    c0 = order-8; c1 = c0 + 4; 
  }
  if (cell[c0] >= 0 && cell[c1] >=0) {
    // both exists
    auto tmp = (*getCell(c0)); 
    for (auto f : tmp->face) {
      if (grid->listFace[f]->prev != cell[c0]) continue; 
      if (grid->listFace[f]->next != cell[c1]) continue; 
      return grid->listFace[f]; 
    }
  } else if (cell[c0] >= 0) {
    auto tmp = (*getCell(c0)); 
    for (auto f : tmp->face) {
      if (grid->listFace[f]->prev != cell[c0]) continue;
      return grid->listFace[f]; 
    }
  } else if (cell[c1] >= 0) { 
    auto tmp = (*getCell(c1)); 
    for (auto f : tmp->face) {
      if (grid->listFace[f]->next != cell[c1]) continue;     
      return grid->listFace[f]; 
    }
  } 
  return NULL; 
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

// Vec3 Vertex::x2xhat(x) {
  
// }

// -------------------------------------------------------------------------
// @EU#9-21-2015: distance based interpolation (easiest thing possible)
// @EU#11-20-2015: Linear interpolation (not that easy)
// -------------------------------------------------------------------------
bool Vertex::setInterpCoef() {
  // save some time by ommitting valid coefs during adaptation; 
  if (!coefUpdate) return false; 

  // Form a CV and prepare for transformation; 
  VecX<double> x(8), y(8), z(8);

  bool bndr = true; 
  // 1. GET CELL COORDINATES
  for (auto i=0; i<8; ++i) {
    if (i < cell.size()) {
      if (cell[i] >= 0) {
	auto tmp = (*getCell(i))->getCoord();
	x[i] = tmp.x(); 
	y[i] = tmp.y(); 
	z[i] = tmp.z();
      } else {
	if (cell.size() == 4) {
	  // Set BC and coefficients and move on
	  if (cell[(i+1)%4] >= 0 && cell[(i+3)%4] >= 0) {
	    cout << "Vertex.cpp: Internal corner cell is not yet supported" << endl; 
	    exit(1); 
	    // internal corner cell //  Not yet supported !!
	  } else if (cell[(i+1)%4] >=0) {
	    // flat boundary on (i+1)%4
	    auto n2 = *(*getCell((i+1)%4))->getVertex(i); 
	    Vec3 tmp = 0.5*(this->getCoord() + n2->getCoord()); 
	    x[i] = tmp.x(); 
	    y[i] = tmp.y(); 
	    z[i] = tmp.z(); 
	    
	  } else if (cell[(i+3)%4] >=0) {
	    // flat boundary on (i+3)%4
	    auto n2 = *(*getCell((i+3)%4))->getVertex(i); 
	    Vec3 tmp = 0.5*(this->getCoord() + n2->getCoord()); 
	    x[i] = tmp.x(); 
	    y[i] = tmp.y(); 
	    z[i] = tmp.z(); 

	  } else {
	    // external corner cell
	    x[i] = this->x(); 
	    y[i] = this->y(); 
	    z[i] = this->z(); 
	  } 
	}
      }
    } else { 
      x[i] = x[i-4]; 
      y[i] = y[i-4]; 
      z[i] = z[i-4]; 
    }
  }
  
  // 2. USE AI to calculate its coefficients; 
  xcoef = grid->interp.getCoef(x); 
  ycoef = grid->interp.getCoef(y); 
  zcoef = grid->interp.getCoef(z);

  // 3. ALSO get xhat at its location; 

  xhat = grid->interp.findXhat(Vec3(this->x(), this->y(), this->z()), 
			      xcoef, ycoef, zcoef); 
  coefUpdate = false; 
  return true; 
}

vector<double> Vertex::getIntWeight() {
  if (coefUpdate) {
    cout << "Error: Coefficients are not yet computed! " << id << endl;
    exit(1);
  }
  vector<double> w(cell.size(), 0);
  w[0] = (1-xhat.x())*(1-xhat.y())*(1-xhat.z()); 
  w[1] = xhat.x()*(1-xhat.y())*(1-xhat.z()); 
  w[2] = xhat.x()*xhat.y()*(1-xhat.z()); 
  w[3] = (1-xhat.x())*xhat.y()*(1-xhat.z()); 
  if (cell.size() == 8) {
    w[4] = (1-xhat.x())*(1-xhat.y())*xhat.z(); 
    w[5] = xhat.x()*(1-xhat.y())*xhat.z(); 
    w[6] = xhat.x()*xhat.y()*xhat.z(); 
    w[7] = (1-xhat.x())*xhat.y()*xhat.z(); 
  } else {
    w[0] += (1-xhat.x())*(1-xhat.y())*xhat.z(); 
    w[1] += xhat.x()*(1-xhat.y())*xhat.z(); 
    w[2] += xhat.x()*xhat.y()*xhat.z(); 
    w[3] += (1-xhat.x())*xhat.y()*xhat.z(); 
  }
  return w; 
}

double Vertex::evalPhi(shared_ptr<Var> &var) {
  return evalPhi(var->data, var->listBC);  
}

double Vertex::evalPhi(VecX<double> &phi, vector<shared_ptr<Boundary> > const &bc) {
  auto w = getIntWeight();
  double a =0; 
  for (auto j=0; j<w.size(); ++j) {
    auto i = (j < cell.size()) ? j : j-4;  
    if (cell[i] >= 0) {
      a += w[i]*phi[cell[i]]; 
    } else {
      auto bndr = -cell[i]-1; 
      double bcval=0; 
      if (bc[bndr]->type == 0) {
	if (cell.size() == 4) {       
	  if (cell[(i+1)%4] >= 0 && cell[(i+3)%4] >= 0) {
	  } else if (cell[(i+1)%4] >=0) {
	    bcval = bc[bndr]->a_val*phi[cell[(i+1)%4]] + bc[bndr]->b_val; 
	  } else if (cell[(i+3)%4] >=0) {
	    bcval = bc[bndr]->a_val*phi[cell[(i+3)%4]] + bc[bndr]->b_val; 
	  } else {
	    bcval = bc[bndr]->a_val*phi[cell[(i+2)%4]] + bc[bndr]->b_val; 
	  } 
	}
      } else if (bc[bndr]->type == 1) { 
	if (cell.size() == 4) {
	  if (cell[(i+1)%4] >= 0 && cell[(i+3)%4] >= 0) {
	  } else if (cell[(i+1)%4] >=0) {
	    auto cg = (*getCell((i+1)%4)); 
	    auto n2 = *(cg->getVertex(i)); 
	    auto dn = (0.5*(this->getCoord() + n2->getCoord()) - cg->getCoord()); 
	    auto del = (bndr%2 == 0) ? -dn.abs() : dn.abs(); 	    	    
	    bcval = (del*bc[bndr]->a_val+1.0)*phi[cell[(i+1)%4]] + del*bc[bndr]->b_val; 
	  } else if (cell[(i+3)%4] >=0) {
	    auto cg = (*getCell((i+3)%4)); 
	    auto n2 = *(cg->getVertex(i)); 
	    auto dn = (0.5*(this->getCoord() + n2->getCoord()) - cg->getCoord()); 
	    auto del = (bndr%2 == 0) ? -dn.abs() : dn.abs(); 	    	    
	    bcval = (del*bc[bndr]->a_val+1)*phi[cell[(i+3)%4]] + del*bc[bndr]->b_val; 
	  } else {
	    auto cg = (*getCell((i+2)%4));
	    auto dn = this->getCoord() - cg->getCoord(); 
	    auto del = (bndr%2 == 0) ? -dn.abs() : dn.abs();
	    bcval = (del*bc[bndr]->a_val+1.0)*phi[cell[(i+2)%4]] + del*bc[bndr]->b_val; 
	  } 
	}
      }
      a += w[i]*bcval; 
    }
  } 
  return a; 
}

  // for (auto i = 0; i < cell.size(); ++i) {
  //   if (cell[i] != cell[(i+cell.size()-1)%cell.size()]) { // remove duplicates; 
  //     auto x = grid->listGrid[i]->getCoord(); 
  
  // if (cell.size() == 4) {
  //   MatX<double> AI = {{1 0 0 0}, {-1 1 0 0}, {-1 0 0 1}, {1 -1 1 -1}};
    
  // }

  /* RELEVANT MATLAB CODE
    a = A\cell(:,1);
    b = A\cell(:,2);
    c = A\cell(:,3);    
    
    gg1 = @(l, c) c(2) + c(5)*l(2) + c(6)*l(3) + c(8)*l(2)*l(3); 
    gg2 = @(l, c) c(3) + c(5)*l(1) + c(7)*l(3) + c(8)*l(1)*l(3); 
    gg3 = @(l, c) c(4) + c(6)*l(1) + c(7)*l(2) + c(8)*l(1)*l(2);     
 
    l = [0.5, 0.5, 0.5]; 
    dl = [1,1,1];
    if (a(2) == 0 && a(5) == 0 && a(6) == 0 && a(8) == 0)
        l(1) = 0; 
        F = @(l) [x(2)-f(l, a); x(3)-f(l, c)];
        J = @(l) [gg1(l, b), gg2(l, b); ...
            gg1(l, c), gg2(l, c)];
        while (sum(abs(dl)) > 1e-6)
            dl = [0; J(l)\F(l)];
            l = l+dl.';
        end
    elseif (b(2) == 0 && b(5) == 0 && b(6) == 0 && b(8) == 0)
        l(2) = 0; 
        F = @(l) [x(1)-f(l, a); x(3)-f(l, c)];
        J = @(l) [gg1(l, a), gg2(l, a); ...
            gg1(l, c), gg2(l, c)];
        while (sum(abs(dl)) > 1e-6)
            dl = J(l)\F(l);
            dl = [dl(1); 0; dl(2)]; 
            l = l+dl.';
        end        
    elseif (c(2) == 0 && c(5) == 0 && c(6) == 0 && c(8) == 0)
        l(3) = 0; 
        F = @(l) [x(1)-f(l, a); x(2)-f(l, b)];
        J = @(l) [gg1(l, a), gg2(l, a); ...
            gg1(l, b), gg2(l, b)];
        while (sum(abs(dl)) > 1e-6)
            dl = [J(l)\F(l); 0];
            if (isnan(sum(dl)) || isinf(sum(dl)))
                break; 
            end            
            l = l+dl.';            
        end
    else
        F = @(l) [x(1)-f(l, a); x(2)-f(l, b); x(3)-f(l, c)];
        J = @(l) [gg1(l, a), gg2(l, a), gg3(l, a); ...
            gg1(l, b), gg2(l, b), gg3(l, b); ...
            gg1(l, c), gg2(l, c), gg3(l, c)];
        while (sum(abs(F(l))) > 1e-6)
            dl = J(l)\F(l);
            if (isnan(sum(dl)) || isinf(sum(dl)))
                break; 
            end
            l = l+dl.';            
        end

    end

    

   */
//}


