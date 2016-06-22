#ifndef SEARCHGRID
#define SEARCHGRID


double phi(shared_ptr<Var> &q, Vec3 x, int_8 i0=0) {
  auto n = searchVertexbyCoords(x, i0); 
  if (n >= 0) listVertex[n]->evalPhi(q, &x); 
} 



int_8 searchVertexbyCoords(Vec3 a, int_8 i0=0) {
  if (i0 >= listVertex.size() && i0 < 0) {
    cout << "Initial vertex do not exists" << endl; 
    exit(-1); 
  }
  // for (auto v : listVertex) {
  //   if (v->isIn(a)) return v->id; 
  // }
  // return -1; 

  //  cout << endl << " ---------- " << " Searching for " << a << " -------- " << endl; 

  auto v = listVertex[i0]; 
  vector<bool> blist(listVertex.size(), false); 
  int_8 i1 = i0; 
  // else search
  while (v.get()) {
    //cout << " + Try: "<< i0 << " @ ( " << *v << " ) "; 
    if (v->isIn(a)) { 
      //cout << "[OK] " << endl; 
      return v->id; 
    }
    //cout << "[!] blisted " << endl; 
    blist[v->id] = true; 

    // Search for next node; 
    double rmin = 100*(*v-a).abs();
    //cout << " + dist: " << rmin << endl; 
    for (auto idir = -2; idir < 3; ++idir) {
      if (idir == 0) continue; 
      //cout << "   - dir : " << idir << " " ; 
      if (auto vnptr = v->ngbr(idir)) {
	auto vn = *vnptr; 
	//cout << vn->id << " "; 
	if (blist[vn->id]) { 
	  //cout << " [!] blist "<< endl; 
	  continue;
	}
	//	if (vn->isIn(a)) { cout << " [OK] " << endl; return vn->id; }
	// else measure distance; 
	auto r = (*vn-a).abs(); 
	// cout << "dist: " << r << " ( " << rmin << ") "; 
	if (r < rmin) {
	  i1 = vn->id;
	  rmin = r; 
	  //cout << " next cand v @ " << *v << endl; 
	} else {
	  //cout << " [!] try next dir " <<endl; 
	}
      } else { 
	//cout << "[!] bndr/hanging " << endl;  
      }    
    }
    //if (i0 == i1) { cout << " No new candidate " << endl; }
    v = listVertex[i1]; i0 = i1; 
    //    cout << " + New cand: "<< i1 << " @ ( " << *v << " ) " << endl; 
    //cin.ignore().get(); 
    if (blist[v->id]) { cout << "Search failed!!!!" << endl; cin.ignore().get(); return -1;}    
  }
    
  return -1; 


  
  
  // while (v.get()) {   
    
  //   bool isMinFound = false; 
  //   auto rmin = (*v-a).abs(); 
  //   if (v->cell.size() == 4) {      
  //     for (auto ci : v->cell) {
  // 	if (ci < 0 || ci >= listCell.size()) continue;	
  // 	for (auto vi : listCell[ci]->node) {
  // 	  if (vi < 0 || vi >= listVertex.size()) continue; 	  
  // 	  if (listVertex[vi]->isIn(a)) {
  // 	    return vi; 
  // 	  } else {
  // 	  //if (listVertex[vi]->isHanging()) continue;
  // 	  auto r = (*listVertex[vi] - a).abs();
  // 	  if (r < rmin) {
  // 	    isMinFound = true; v = listVertex[vi]; rmin = r; 
  // 	  }
  // 	}
  //     }
  //     if (!isMinFound) {
  // 	if (v->isIn(a)) {
  // 	  break; 
  // 	} else { 
  // 	  rmin = 100*rmin; 
  // 	  vblacklist.emplace_back(v->id); 
  // 	  continue; 
  // 	}
  //     }
  //   }     
  // }
  
  // // minimum distance does not guarantee connectivity 
  // // another search should be addressed! 
  // // if (!isFound) {
    
  // //   // check
  // // }


  // if (v.get()) return v->id; 
  // else return -1; 
}


int_8 searchInterpPoint(Vec3 a, int_8 i0) {
  if (i0 >= listVertex.size() && i0 < 0) {
    cout << "Initial vertex do not exists" << endl; 
    exit(-1); 
  }
  auto v = listVertex[i0]; 
  cout << endl << "---------" << endl << "Search: " << a << " from nearby vertex " << i0 << " @ " << *v << endl; 
  for (auto ci : v->cell) {
    if (ci < 0 || ci >= listCell.size()) continue; 	
    for (auto vi : listCell[ci]->node) {
      if (vi < 0 || vi >= listVertex.size()) continue; 
      cout << "+ vertex " << vi << " @ "<< *listVertex[vi] << " CoefUpdate " << listVertex[vi]->coefUpdate << endl; 
      if (listVertex[vi]->coefUpdate) continue; 
      auto xhat = int3D.findXhat(a, listVertex[vi]->xcoef, listVertex[vi]->ycoef, listVertex[vi]->zcoef); 
      cout << "         xhat = " << xhat ; 
      if (int3D.isIn(xhat)) { cout << " [FOUND] " << vi << endl; return vi; }
      cout << " [NOT FOUND] continue " << endl; 
    }
  }
  cout << " @  NOT FOUND " << endl; 
  cin.ignore().get(); 
  return -1; 
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
