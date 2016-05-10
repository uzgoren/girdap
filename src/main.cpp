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
#include <linSolve/MatY.hpp>


int main() { 
  triLinSys Ab; 
  
  Ab.vt.emplace_back(Triplet(0,0,4)); 
  Ab.vt.emplace_back(Triplet(0,1,-1)); 
  Ab.vt.emplace_back(Triplet(1,2,-1)); 
  Ab.vt.emplace_back(Triplet(2,1,-1)); 
  Ab.vt.emplace_back(Triplet(2,2,4));   
  Ab.vt.emplace_back(Triplet(3,2,-1)); 
  Ab.vt.emplace_back(Triplet(3,3,4));   
  Ab.vt.emplace_back(Triplet(3,4,-1)); 
  Ab.vt.emplace_back(Triplet(4,4,4));   
  Ab.vt.emplace_back(Triplet(4,3,-1)); 
  Ab.vt.emplace_back(Triplet(1,1,4));   
  Ab.vt.emplace_back(Triplet(2,3,-1)); 
  Ab.vt.emplace_back(Triplet(1,0,-1)); 
  Ab.vt.emplace_back(Triplet(2,1,-2)); 
  Ab.vt.emplace_back(Triplet(2,2,1)); 
  Ab.vt.emplace_back(Triplet(1,4,-6)); 

  Ab.setMat(); 
  
  for (auto i = 0; i < Ab.aij.size(); ++i) { 
    cout << i << " " << Ab.aij[i] << " " << Ab.val[i] << endl; 
  }

  return 0; 
};
