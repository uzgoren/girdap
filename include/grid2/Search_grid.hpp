#ifndef SEARCHGRID
#define SEARCHGRID


double phi(shared_ptr<Var> &q, Vec3 x, int_8 i0=0) {
  auto n = searchVertexbyCoords(x, i0); 
  if (n >= 0) listVertex[n]->evalPhi(q, &x); 
} 



int_8 searchVertexbyCoords(Vec3 a, int_8 i0=0) {
  //cout << " Search in " << endl; 
  if (i0 >= listVertex.size() && i0 < 0) {
    cout << "Initial vertex do not exists" << endl; 
    exit(-1); 
  }
  auto v = listVertex[i0]; 
  while (v) {            
    bool isFound = false; 
    auto rc = (*v-a).abs(); auto rmin = rc; 
    if (v->cell.size() == 4) {
      for (auto ci : v->cell) {
	if (ci < 0) continue; 
	for (auto vi : listCell[ci]->node) {
	  auto r = (*listVertex[vi] - a).abs();
	  if (r < rmin) {
	    isFound = true; v = listVertex[vi]; rmin = r; 
	  }
	}
      }
      if (!isFound) break; 
    }     
  }
  //cout << " Search out " << *v << endl; 
  if (v) return v->id; 
  else return -1; 
}

int_8 searchCellbyCoords(Vec3 a, int_8 i0=0, bool isVertex=false) {
  auto j0 = isVertex ? i0 : listCell[i0]->node[0]; 
  auto vi = searchVertexbyCoords(a, j0); 
  if (vi < 0) return -1; 
  auto v = listVertex[vi]; 
  
  int_8 vc = -1; int_8 bc = 0; 
  for (auto ic = 0; ic < v->cell.size(); ++ic) {
    if (v->cell[ic] < 0) {bc = v->cell[ic]; continue;}
    // Not correct for curvilinear //
    auto xmin = Vec3(*listVertex[listCell[v->cell[ic]]->node[0]]);
    auto xmax = Vec3(*listVertex[listCell[v->cell[ic]]->node[2]]);
    if (a[0] < xmin[0]) continue; 
    if (a[1] < xmin[1]) continue; 
    if (a[0] > xmax[0]) continue; 
    if (a[1] > xmax[1]) continue; 
    vc = v->cell[ic]; 
  }
  if (vc >= 0) return vc; 
  else if (bc < 0) return bc; 
  else {
    cout << " Cell not found! vc = " << vc; 
    cout << " vi = " << vi << " a="<< a<< endl;
    exit(-1); 
  }
  // auto xhat = int3D.findXhat(a, v->xcoef, v->ycoef, v->zcoef); 
  // if (xhat[0] <= 0.5 && xhat[1] <= 0.5) return v->cell[0];
  // else if (xhat[0] > 0.5 && xhat[1] <= 0.5) return v->cell[1];
  // else if (xhat[0] > 0.5 && xhat[1] > 0.5) return v->cell[2]; 
  // else if (xhat[0] <= 0.5 && xhat[1] > 0.5) return v->cell[3];
  // else return -1; 
}

// int_8 searchInterpVertexbyCoords(Vec3 a, int_8 i0=0) {
//   auto vi = searchVertexbyCoords(a, i0); 
//   if (vi < 0) return -1; 
//   auto v = listVertex[vi]; 
//   auto xhat = int3D.findXhat(a, v->xcoef, v->ycoef, v->zcoef); 
//   if (xhat[0] <= 0.5 && xhat[1] <= 0.5) return v->cell[0];
//   else if (xhat[0] > 0.5 && xhat[1] <= 0.5) return v->cell[1];
//   else if (xhat[0] > 0.5 && xhat[1] > 0.5) return v->cell[2]; 
//   else if (xhat[0] <= 0.5 && xhat[1] > 0.5) return v->cell[3];
//   else return v->cell[0]; 
// }


#endif
