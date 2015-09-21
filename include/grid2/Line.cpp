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
