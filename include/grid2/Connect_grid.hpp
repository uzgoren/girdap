#ifndef CONNECTGRID
#define CONNECTGRID

void contour(Grid* &cGrid, shared_ptr<Var> &a, double d) {
  // Cell by cell construction
  // For each cell create elements for the new grid
  // Check before surface is needed; 

  double val[listVertex.size()]; 
  int_8 connect[listFace.size()]; 
  vector<bool> color(listCell.size(), false); 

  vector<int_8> connectCell;
  vector<vector<int_8> > connectCellNodeList(listCell.size(), vector<int_8>(8,-1));     
    
  // 00. construct values at vertex; not cell dependent
  for (auto i = 0; i < listVertex.size(); ++i) {
    val[i] = listVertex[i]->evalPhi(a); 
  }
  
  // 01. Create vertices at faces; assign them to a cell; cell dependent
  // using val, other grid, let it stay; 
  for (auto f : listFace) {
    auto v0 = f->node[0]; 
    auto v1 = f->node[1]; 
	
    auto dp = d;
    if (abs(val[v0] - dp) < 1e-15) dp = d - 1e-14; 
    if (abs(val[v1] - dp) < 1e-15) dp = d - 1e-14; 
      
    if ((val[v0] - dp)*(val[v1] - dp) < 0) { 
      auto y = (dp - val[v0])/(val[v1] - val[v0]); 
      auto dr = *listVertex[v1] - *listVertex[v0]; 
      auto x = *listVertex[v0] + y*dr; 
      cGrid->addVertex(x);    // add a marker on the vertex; 
      int_8 inew = cGrid->listVertex.size()-1; 
	
      auto n = f->next; 
      auto p = f->prev; 
	
      // Create a link between 
      // Divide cells' edges into two halves and name them 
      // Quad:0-7;  Line:0-1; Tri:0-5; [Hexa:0-41 ??;]
      int_8 nloc, ploc; 
      if (n >= 0) {	  
	auto cn = listCell[n];
	if (f->orient == 0) nloc = (cn->node[0] == v0) ? 0 : 1; 
	else if (f->orient == 1) nloc = (cn->node[3] == v1) ? 6 : 7; 
	else continue; 
	if (!color[n]) connectCell.push_back(n);
	connectCellNodeList[n][nloc] = inew; 
	color[n] = true; 
      }
      if (p >= 0) {	  
	auto cn = listCell[p]; 
	if (f->orient == 0) ploc = (cn->node[2] == v1) ? 4 : 5; 
	else if (f->orient == 1) ploc = (cn->node[1] == v0) ? 2 : 3; 
	else continue; 	 
	if (!color[p]) connectCell.push_back(p);
	connectCellNodeList[p][ploc] = inew; 
	color[p] = true; 
      }	
    }      
  }
  for (auto ic = 0; ic < connectCell.size(); ++ic) {
    auto i = connectCell[ic]; auto sum = 0; 
    for (auto jc = 0; jc<connectCellNodeList[i].size(); ++jc) {
      if (connectCellNodeList[i][jc] < 0) continue; 
      ++sum; 
    }
    if (sum == 2) continue; 
    cout << ic << " Warning! " << i << " " << sum << " are not handled " << endl; 
  }

  // // 02. check cells through nodes;
  //    ic = 0; 
  for (auto ic = 0; ic < connectCell.size(); ++ic) {
    auto i = connectCell[ic]; 
    if (!listCell[i]->isAlive) continue;
    int_8 j0 = 0, j1 = 0;
    //     while (j1 <connectCellNodeList[i].size()-1) {
    int_8 k; 
    for (k = j1; k < connectCellNodeList[i].size()-1; ++k) 
      if (connectCellNodeList[i][k] >= 0) break;
    if (k >= connectCellNodeList[i].size()) continue;       
	
    for (j1 = k+1; j1 < connectCellNodeList[i].size(); ++j1) {
      //	j1 = (k + j) % 8;
      //cout << k << " " << j << " " << j1 << " " << j0 << endl; 
      if (connectCellNodeList[i][j1] >= 0) {
	if (val[listCell[i]->node[j0]] >= d) 
	  cGrid->addCell({connectCellNodeList[i][k], connectCellNodeList[i][j1]}); 
	else
	  cGrid->addCell({connectCellNodeList[i][j1], connectCellNodeList[i][k]}); 
	break; 
      } 
    } 
    // if (j1 >= connectCellNodeList[i].size()) {
    //   cout << "Couldnot find proper match with " << k << "- skipping" << endl; 
    //   continue; 
    // }
    // }
    //      k = j1;
    //cin.ignore().get(); 
  }

  // 03. disregard connect information after putting them in
  //     object array otherVertex; vertex-to-vertex connectivity

  // cGrid->otherVertex.assign(cGrid->listVertex.size(), -1); 
  // otherVertex.resize(listVertex.size()); 

  // auto i0 = connect[0]; 
  // for (auto i1=0; i1 < cGrid->listVertex.size(); ++i1) {
  //   cout << i0 << " " << i1 << " " << cGrid->listVertex.size() <<endl; 
  //   cout << *(cGrid->listVertex[i1]) << " " << listCell[i0]->id << " " << listCell[i0]->node[0] << endl; 
  //   auto v0 = searchVertexbyCoords(*(cGrid->listVertex[i1]), listCell[i0]->node[0]); 
  //   if (v0 >= 0) {
  //     // cGrid->otherVertex[i1] = listVertex[v0]->id;
  //     // i0 = (listVertex[v0]->cell[0]<0) ? i0:listVertex[v0]->cell[0] ; 
  //     // cout << "Grid: Cell = "<< i0 << ", Vertex" << v0 << " -- Other: " << i1 << endl; 
  //   } else {
  //     cout << "Could not found the next vertex" << i0 << " " << *(cGrid->listVertex[i1]) << endl;
  //     exit(1); 
  //   }
  //   cGrid->otherVertex[i1] = listVertex[v0]->id; 
  //   otherVertex[listVertex[v0]->id] = cGrid->listVertex[i1]->id; 
  // }    

  for (auto v : cGrid->listVar) {
    if (v->loc == 1) v->data.resize(cGrid->listVertex.size());
    else if (v->loc == 0) v->data.resize(cGrid->listCell.size()); 
  }


  cGrid->cleanGrid(); 
  updateOtherVertex(cGrid); 
  
  return; 
}

void releaseConnections() {
  for (auto i = 0; i < otherVertex.size(); ++i)
    otherVertex[i] = -1; 
}

void updateOtherVertex(Grid* cGrid) {
  otherVertex.assign(listVertex.size(), -1); 
  cGrid->otherVertex.resize(cGrid->listVertex.size()); 
  
  for (auto i1 = 0; i1 < cGrid->listVertex.size(); ++i1) {
    int_8 i0; 
    if (cGrid->otherVertex[i1] >= 0) {
      i0 = searchVertexbyCoords(*(cGrid->listVertex[i1]), cGrid->otherVertex[i1]); 
      otherVertex[cGrid->otherVertex[i1]]= -1; 
    } else {
      i0 = searchVertexbyCoords(*(cGrid->listVertex[i1])); 
    }
    if (i0 >= 0) {
      cGrid->otherVertex[i1] = listVertex[i0]->id; 
      otherVertex[i0] = cGrid->listVertex[i1]->id; 
    }
  }
 
}


void indicator(Grid* cGrid, shared_ptr<Var> T) {
  updateOtherVertex(cGrid); 

  T->set(0); 
  auto dx = listCell[0]->edge(0).abs();
  auto dy = listCell[0]->edge(1).abs(); 
 
  auto pi = 4.0*atan(1.0); 
  for (auto i = listCell[0]->level[0]; i < levelMax[0]+1; ++i) dx /= 2; 
  for (auto i = listCell[1]->level[1]; i < levelMax[1]+1; ++i) dy /= 2; 

  cout << " ------ " << dx << ", "<< dy << " : " << " " << levelMax[0] << endl; 

  // 00. Now color a new function at the cell centers; 
  double dist[listCell.size()];
  for (auto i = 0; i<listCell.size(); ++i) dist[i] = -1; 

  for (auto ci1 = 0; ci1 < cGrid->listCell.size(); ++ci1) {
    auto c1 = cGrid->listCell[ci1]; 
    auto x1 = c1->getCoord(); 
    auto norm = c1->vol(); norm = norm/norm.abs(); 
    auto i1 = c1->node[0]; 
    if (i1 < 0 && i1 >= cGrid->listVertex.size()) {
      cout << " Cell does not have a vertex ??? " << i1 << endl; 
      exit(1); 
    }
    auto v1 = cGrid->listVertex[i1];     
    auto i0 = cGrid->otherVertex[i1]; 
    if (i0 < 0 && i0 >= listVertex.size()) {
      cout << "WARNING - no associated vertex found! now what? " << endl; 
      exit(1); 
    }
    auto v0 = listVertex[i0];

    int_8 j; 
    for (j = 0; j < v0->cell.size(); ++j) if (v0->cell[j] >= 0) break; 
    
    vector<int_8> ngbrCells(1, v0->cell[j]); 
    ngbrCellbyLayer(5, ngbrCells); 

    for (auto j = 0; j < ngbrCells.size(); ++j) {
      auto ci0 = ngbrCells[j]; 
      auto c0 = listCell[ci0];       
      auto r = (c0->getCoord() - x1).abs(); 
      if (dist[ci0] < 0) dist[ci0] = r; 
      if (dist[ci0] < r) continue; 
      dist[ci0] = r; 
      if (r > 0.04) continue; 
      r = (c0->getCoord() - x1)*norm;
      T->set(ci0, 1.0/(1+exp(-160*r))); //0.25*r/dx*pi)));      
    }    
  }


  int icnt = 0; 
  for (auto i = 0; i < listCell.size(); ++i) {
    int_8 sum = 0; 
    paintCell(T, 0.0, 0.5, 1.0, i, sum);
    // if (sum > 1) {
    //   cout << endl << " Paint cells: " << sum << " icnt: " << icnt << endl; 
    //   if (icnt > 1) {
    // 	cout << " Found more than one bodies " << listCell[i]->getCoord() << endl; 
    // 	exit(1); 
    //   }
    //   ++icnt; 
    // }
  }
       
}


void passVar(Grid* cGrid, shared_ptr<Var> from, shared_ptr<Var> to) {
  for (auto i = 0; i < cGrid->listVertex.size(); ++i) { 
    auto v1 = cGrid->listVertex[i]; 
    auto v0s = cGrid->otherVertex[i]; //searchVertexbyCoords(*v1, v0->id);
    if (v0s >= 0) {
      auto xhat = interp.findXhat(*v1, listVertex[v0s]->xcoef, listVertex[v0s]->ycoef, listVertex[v0s]->zcoef);
      to->set(i, listVertex[v0s]->evalPhi(from)); 
    } else {
      cout << " no vertex found! now what! 2"<< endl;
    }
  }
}


// vector<int_8> cellsAroundVertex(Vec3 &p, int layer, int i0=0, bool isVertex=false) {
//   auto c = searchCellbyCoords(p, i0, isVertex); 
// }


vector<int_8> cellsAroundInterface(Grid* &a, int layer) {

  vector<int_8> tmp; 
  bool color[listCell.size()]; 
  layer = 0; 

  for (auto i = 0; i<listCell.size(); ++i) color[i] = false; 

  for (auto i1 = 0; i1 < a->listVertex.size(); ++i1) {
    auto v0i = a->otherVertex[i1]; 
    if (v0i <= 0) {
      cout << "WARNING - no associated vertex found! now what? " <<endl; 
      exit(1); 
    } 
    auto v0 = listVertex[v0i];
    for (auto li = 0; li < v0->cell.size(); ++li) {
      auto ci0 = v0->cell[li]; 
      if (ci0 < 0) continue; 
      if (color[ci0]) continue; 
      tmp.emplace_back(ci0); 
      color[ci0] = true;
    }
  }
  return std::move(tmp); 
}

// Given layer, function fills in ngbrCells recursively
// layer is what stops adding
//     color : indicates those visited; 
//     cp    : is the index of the current cell sought in ngbrCells; 
//     
// The function should be called with the location of the initial cell 
// placed as an element in ngbrCells
void ngbrCellbyLayer(int layer, vector<int_8> &ngbrCells, int_8 cp=0) {  
  if (layer <= 0) return; // no change
  if (cp >= ngbrCells.size()) return; 
  if (cp < 0) return; 
  if (ngbrCells[cp] >= listCell.size()) return; 
  if (ngbrCells[cp] < 0) return; 
  auto c = listCell[ngbrCells[cp]];
  for (auto i = 0; i < c->node.size(); ++i) {
    auto n = c->node[i]; 
    auto cp2 = listVertex[n]->cell[(i+1)%4]; 
    if (cp2 >= 0) {
      auto is = std::find(ngbrCells.begin(), ngbrCells.end(), cp2); 
      if (is == ngbrCells.end()) { // new entry
  	ngbrCells.emplace_back(cp2); 
	is = ngbrCells.end()-1; 
      }
      ngbrCellbyLayer(layer-1, ngbrCells, is - ngbrCells.begin()); //ngbrCells.size()-1);       
    }
  }  
}

void ngbrCellbyDist(Vec3 r, vector<int_8> &ngbrCells, int_8 cp=0) {
  //cout << "ind = "<< cp << " r = (" << r << " ), size = " << ngbrCells.size() << endl; 
  if (r[0] < 0 || r[1] < 0  || r[2] < 0) return; 
  if (cp >= ngbrCells.size()) return; 
  if (cp < 0) return; 
  if (ngbrCells[cp] >= listCell.size()) return; 
  if (ngbrCells[cp] < 0) return; 
  auto c = listCell[ngbrCells[cp]];
  for (auto i = 0; i < c->node.size(); ++i) {
    auto n = c->node[i]; 
    auto cp2 = listVertex[n]->cell[(i+1)%4]; 
    if (cp2 >= 0) {
      auto is = std::find(ngbrCells.begin(), ngbrCells.end(), cp2); 
      if (is == ngbrCells.end()) { // new entry
  	ngbrCells.emplace_back(cp2); 
	is = ngbrCells.end()-1; 
      }      
      Vec3 dx = listCell[cp2]->getCoord() - c->getCoord(); 
      dx[0] = r[0] - abs(dx[0])+1e-15; 
      dx[1] = r[1] - abs(dx[1])+1e-15; 
      dx[2] = r[2] - abs(dx[2])+1e-15; 
      ngbrCellbyDist(dx, ngbrCells, is - ngbrCells.begin()); //ngbrCells.size()-1);       
    }
  }  
}

void paintCell(shared_ptr<Var> &v, double replace, double bndr, double value, int_8 cid, int_8 &sum) {
  //cout << "ind = "<< cp << " r = (" << r << " ), size = " << ngbrCells.size() << endl; 
  if (v->get(cid) < bndr) return; 
  // get neighbors; 
  vector<int_8> ngbrCell(1, cid); 
  ngbrCellbyLayer(1, ngbrCell); 
  for (auto i = 1; i < ngbrCell.size(); ++i) {
    if (v->get(ngbrCell[i]) != replace) continue; 
    ++sum; 
    v->set(ngbrCell[i], value);
    paintCell(v, replace, bndr, value, ngbrCell[i], sum); 
  }
}


void contour2(Grid* &cGrid, shared_ptr<Var> &a, double d) {
  if (cGrid->listVertex.size() > 0) {
    cGrid->listVertex.clear();
    cGrid->listCell.clear(); 
    for (auto v : cGrid->listVar) 
      v->resize(0);
    cGrid->otherVertex.clear(); 
  } else { 
    //contourFromPrev(cGrid, a, d); 
  }
  otherVertex.assign(listVertex.size(), -1);  
  contourFromScratch(cGrid, a, d); 
} 



void contourFromScratch(Grid* &cGrid, shared_ptr<Var> &a, double d) {

  double val[listVertex.size()]; 
  
  // 00. construct values at vertex; not cell dependent
  for (auto i = 0; i < listVertex.size(); ++i) 
    val[i] = listVertex[i]->evalPhi(a); 
  
  // 01. check near vertices if contour appears
  for (auto ia = 0; ia < listVertex.size(); ++ia) {
    if (otherVertex[ia] >= 0) continue; // already assigned
    auto va = listVertex[ia]; 
    for (auto di = 1; di <= 3; ++di) { 
      if (di == 0) continue; // no such direction
      auto vb = va->ngbr(di); 
      if (vb == nullptr) continue; 
      auto ib = (*vb)->id; 

      auto dp = d;
      if (abs(val[ia] - dp) < 1e-15) dp = d - 1e-14; 
      if (abs(val[ib] - dp) < 1e-15) dp = d - 1e-14; 
            
      if ((val[ia] - dp)*(val[ib] - dp) >= 0) continue; 

      // create a new point along xia - xib if yet not assigned. 
      if (otherVertex[ia] < 0) {
	auto y = (dp - val[ia])/(val[ib] - val[ia]); 
	auto dr = *listVertex[ib] - *listVertex[ia]; 
	auto x = *listVertex[ia] + y*dr; 
	auto i0a = searchVertexbyCoords(x, ia);
	if (i0a >= 0 && i0a < listVertex.size() && otherVertex[i0a] < 0) {
	  cGrid->addVertex(x);
	  auto i1a = cGrid->listVertex.size()-1; 	  
	  // assign this node to a; 	  
	  cGrid->otherVertex[i1a] = i0a; 
	  otherVertex[i0a] = i1a; 
	}
      }      
    }      
  }

  if (cGrid->listVertex.size() == 0) return; 

  cout << " Vertex count: " << cGrid->listVertex.size() << " "  ; 
  // initial point to start tracing interface ! and ends here; 
  auto i1a = 0; 
  auto i0a = cGrid->otherVertex[i1a];
  if (i0a < 0 || i0a > listVertex.size()) {
    cout << " Vertex cannot be placed on the grid" << endl; 
    exit(1); 
  }
  
  ofstream ot; 
  ot.open("log"+std::to_string(filecnt)+".dat"); 
  // now I have the vertices; 
  auto v0 = listVertex[i0a]; 
  int_8 i1b; 
  shared_ptr<Vertex> *vn; 
  // Calculate gradient first; 
  int di; Vec3 gradv; 
  auto i1ref = i1a;

  auto icnt = 0; 
  vector<bool> color(listVertex.size(), false); 

  while (icnt < listVertex.size()) { 
    ++icnt; 
    vector<pair<int, double> > pp; 
    for (auto di = -3; di <= 3; di++) {
      if (di == 0) continue; 
      vn = v0->ngbr(di); 
      if (vn == nullptr) continue; 
      if (color[(*vn)->id]) continue;
      pp.push_back(std::make_pair<int, double>(int(di), abs(val[(*vn)->id] - d)));
      gradv[abs(di)-1] = (val[(*vn)->id]-val[i0a])/((**vn)[abs(di)-1]-(*v0)[abs(di)-1]) ;
    }
    if (pp.size() == 0) break; // nothing to pursue; 

    std::sort(pp.begin(), pp.end(), 
	      [](const std::pair<int,double> &left, const std::pair<int,double> &right) {
		return left.second < right.second;
	      });

    vn = nullptr; 
    for (auto ii = pp.begin(); ii !=pp.end(); ++ii) {
      if (ii->first == 1 && gradv[1] <= 0) continue; 
      if (ii->first == -1 && gradv[1] >= 0) continue; 
      if (ii->first == 2 && gradv[0] >= 0) continue; 
      if (ii->first == -2 && gradv[0] <= 0) continue; 
      vn = v0->ngbr(ii->first); 
      break; 
    } 
    if (vn == nullptr) {
      cout << endl << "Couldn't get the next vertex" << endl; 
      exit(1); 
    } 
    
    ot << i0a << " " << (*vn)->id << endl;     

    i0a = (*vn)->id; 
    v0 = (*vn); 
    color[i0a] = true; 
    if (abs(val[i0a]) < 1e-5 && abs(val[i0a] - 1.0) < 1e-5) {
      cout << "trace is not successful! stopping"<< endl; 
      cGrid->writeVTK("lag"); 
      exit(1); 
    }
    
    if (otherVertex[i0a] < 0) continue; 
    cGrid->addCell({i1a, otherVertex[i0a]}); 
    if (otherVertex[i0a] == i1ref) break; 
    i1a = otherVertex[i0a]; 
  }

  cout  <<" Cell count: " << cGrid->listCell.size() <<endl; 
  if ((*(cGrid->listCell.rbegin()))->node[1] != (*(cGrid->listCell.begin()))->node[0]) {
    cout << " curve is not closed" << endl; 
        cGrid->listVar.clear(); 
    writeVTK("euler"); 
    cGrid->writeVTK("lag"); 
    exit(1); 
  }

  if (icnt >= listVertex.size()) { 
    cout << "Trace is not successful! stopping " << endl; 
    cGrid->listVar.clear(); 
    writeVTK("euler"); 
    cGrid->writeVTK("lag"); 
    exit(1); 
  }
    
  ot.close(); 
  auto nv2 = cGrid->listVertex.size(); 
  auto nc2 = cGrid->listVertex.size(); 


  for (auto v : cGrid->listVar) {
    if (v->loc == 1) v->data.resize(cGrid->listVertex.size());
    else if (v->loc == 0) v->data.resize(cGrid->listCell.size()); 
  }


  cGrid->cleanGrid(); 
  // if (nv2 != cGrid->listVertex.size()) {
  updateOtherVertex(cGrid); 
  //   // for (auto i1 = 0; i1 < cGrid->otherVertex.size(); ++i1) {
  //   //   if (i1 != cGrid->listVertex[i1]->id) {
  //   // 	cout << " Complaint; ids have problems" << i1; 
  //   // 	cout << " " << cGrid->listVertex[i1]->id << endl;
  //   // 	exit(1); 
  //   //   }
  //   //   auto i0 = cGrid->otherVertex[i1]; 
  //   //   auto is = searchVertexbyCoords(*cGrid->listVertex[i1], i0); 
  //   //   otherVertex[i0] = i1; 
  //   // }
  // }

  return; 
}
 

   
       


  

#endif
