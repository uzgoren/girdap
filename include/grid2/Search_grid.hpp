#ifndef SEARCHGRID
#define SEARCHGRID

int_8 searchVertexbyCoords(Vec3 a, int_8 i0=0) {
  if (i0 >= listVertex.size() && i0 < 0) {
    cout << "Initial vertex do not exists" << endl; 
    exit(-1); 
  }
  auto v = listVertex[i0]; 
  shared_ptr<Vertex >* tmp; 
  while (v) {            
    bool isFound = false; 
    auto rc = (*v-a).abs();  auto rmin = rc; 
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
  if (v) return v->id;     
  else return -1; 
}

int_8 searchCellbyCoords(Vec3 a, int_8 i0=0, bool isVertex=false) {
  auto j0 = isVertex ? i0 : listCell[i0]->node[0]; 
  auto vi = searchVertexbyCoords(a, j0); 
  if (vi < 0) return -1; 
  auto v = listVertex[vi]; 
  auto xhat = interp.findXhat(a, v->xcoef, v->ycoef, v->zcoef); 
  if (xhat[0] <= 0.5 && xhat[1] <= 0.5) return v->cell[0];
  else if (xhat[0] > 0.5 && xhat[1] <= 0.5) return v->cell[1];
  else if (xhat[0] > 0.5 && xhat[1] > 0.5) return v->cell[2]; 
  else if (xhat[0] <= 0.5 && xhat[1] > 0.5) return v->cell[3];
  else return v->cell[0]; 
}


#endif
