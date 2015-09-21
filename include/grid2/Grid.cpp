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

void Grid::addVar(std::string n, int t) {
  if (!getVar(n)) {
    listVar.emplace_back(shared_ptr<Var>(new Var(n, t))); 
    if (t == 1) {
      (*listVar.rbegin())->data.resize(listVertex.size());
    } else if (t == 2) {
      (*listVar.rbegin())->data.resize(listFace.size());
    } else {
      (*listVar.rbegin())->data.resize(listCell.size());
    }
    (*listVar.rbegin())->grid = this; 
    (*listVar.rbegin())->solver = "BiCGSTAB";
    (*listVar.rbegin())->initBC(); 
  } else 
    cout << "Skipping " << n << " as it exists"<< endl; 
} 

void Grid::addVar(initializer_list<std::string> nl) {
  for (auto n: nl) addVar(n); 
}

shared_ptr<Var> Grid::getVar(std::string n) { 
  auto matching_iter = std::find_if(listVar.begin(), listVar.end(),
				    [&n](const shared_ptr<Var> p) {
				      return n == p->name;
				    });    
  return (matching_iter!=listVar.end()) ? *matching_iter : NULL; 
}
  
void Grid::lockBC(shared_ptr<Var> v) {
  thisVar = v;
  v->grad = valGrad(v); 
  //  v->grad.assign(listCell.size(), 0)
  // v->grad.clear(); 
  // v->grad.resize(listCell.size()); 
  // for (auto i = 0; i < listCell.size(); ++i) {
  //   v->grad[i] = listCell[i]->grad(v); //v->grad = valGrad(v); 
  // }
}
void Grid::unlockBC() {
  thisVar.reset(); 
}
