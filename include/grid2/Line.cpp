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
