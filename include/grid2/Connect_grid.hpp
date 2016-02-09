#ifndef CONNECTGRID
#define CONNECTGRID

Grid* contour(shared_ptr<Var> &a, double d) {
  Grid* cGrid = new Grid(); 
    
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
      int nloc, ploc; 
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
  cGrid->otherVertex.resize(cGrid->listVertex.size()); 
  otherVertex.resize(listVertex.size()); 

  addVar("b"); 
  auto b = getVar("b"); 
  b->set(0); 
  cGrid->addVar("dummy",1);
    
  auto dummy = cGrid->getVar("dummy"); 
  dummy->set(0); 

  auto i0 = connect[0]; 
  for (auto i1=0; i1 < cGrid->listVertex.size(); ++i1) {
    auto v0 = searchVertexbyCoords(*(cGrid->listVertex[i1]), listCell[i0]->node[0]); 
    if (v0 >= 0) {
      cGrid->otherVertex[i1] = listVertex[v0];
      i0 = (listVertex[v0]->cell[0]<0) ? i0:listVertex[v0]->cell[0] ; 
      // cout << "Grid: Cell = "<< i0 << ", Vertex" << v0 << " -- Other: " << i1 << endl; 
    } else {
      cout << "Could not found the next vertex" << i0 << " " << *(cGrid->listVertex[i1]) << endl;
      exit(1); 
    }
    cGrid->otherVertex[i1] = listVertex[v0]; 
    otherVertex[listVertex[v0]->id] = cGrid->listVertex[i1]; 
  }    


  auto dr = listCell[0]->edge(0).abs();
 
  auto f = 4.0*atan(1.0); 
  for (auto i = listCell[0]->level[0]; i < levelMax[0]+1; ++i) dr /= 2; 

  cout << " ------ " << dr << " : " << f << " " << levelMax[0] << endl; 
  

  cGrid->addVec("u", 1); 
  cGrid->addVec("n", 1);
  auto nx = cGrid->getVar("nx"); 
  auto ny = cGrid->getVar("ny"); 
  auto nz = cGrid->getVar("nz");   
  // 04. Now color a new function at the cell centers; 
  double dist[listCell.size()];
  for (auto i = 0; i<listCell.size(); ++i) dist[i] = -1; 

  for (auto i1 = 0; i1 < cGrid->listVertex.size(); ++i1) {
    auto v1 = cGrid->listVertex[i1]; 
    Vec3 norm; 
    if (v1->cell[0] >=0 && v1->cell[1] >= 0) {
      auto e0 = (*v1->getCell(0))->vol(); auto b0 = e0.abs(); 
      auto e1 = (*v1->getCell(1))->vol(); auto b1 = e1.abs(); 
      
      norm = (b1*e0 + b0*e1)/(b0+b1); 
    } else if (v1->cell[1] >= 0) {
      norm = (*v1->getCell(1))->vol(); 
    } else if (v1->cell[0] >= 0) {
      norm = (*v1->getCell(0))->vol(); 
    } else {
      norm = Vec3(-1, 0, 0); 
    }
    norm /= norm.abs();
    nx->set(i1, norm[0]); 
    ny->set(i1, norm[1]); 
    nz->set(i1, norm[2]);

    auto v0 = cGrid->otherVertex[i1]; 
    if (v0 < 0) {
      cout << "WARNING - no associated vertex found! now what? " <<endl; 
      exit(1); 
    }
    for (auto j = 0; j < v0->cell.size(); ++j) {
      auto ci0 = v0->cell[j]; 
      if (ci0 <= 0) continue; 
      auto c0 = listCell[ci0]; 
      auto r = (c0->getCoord() - *v1).abs(); 
      if (dist[ci0] < 0) dist[ci0] = r; 
      if (dist[ci0] < r) continue; 
      r = (c0->getCoord() - *v1)*norm/dr;       
      b->set(ci0, 1.0/(1+exp(-r*f)));
       
      // for (auto jj = 0; jj < c0->node.size(); ++jj) {
      // 	auto vv0 = listVertex[c0->node[jj]]; 
      // 	//cout << " ** " << vv0 << endl; 
      // 	for (auto ii = 0; ii < vv0->cell.size(); ++ii) {
      // 	  auto ci00 = vv0->cell[ii]; 
      // 	  if (ci00 <= 0) continue; 
      // 	  auto cc = listCell[ci00]; 	 
      // 	  r = (cc->getCoord() - *v1).abs(); 
      // 	  if (dist[ci00] < 0) dist[ci00] = r;
      // 	  if (dist[ci00] < r) continue; 
      // 	  r = (cc->getCoord() - *v1)*norm/dr; 
      // 	  b->set(ci0, 1.0/(1+exp(-r*f)));
      // 	  // val = 1.0/(1+exp(r*f));
      // 	  // cout << " ***r = "<<r << " - val = "<<val<< " curval = " << b->get(ci0); 
      // 	  // if (abs(val-0.5) < abs(b->get(ci00) - 0.5)) b->set(ci00, val); 
      // 	  // cout << " newval = " << b->get(ci0) <<endl; 
      // 	}	  	 
      // }	
    }       
  }
    // auto v1p = v1->ngbr(-1); 
    // auto v1n = v1->ngbr(1); 
    // norm1
  
  vector<int_8> aroundCells = cellsAroundInterface(cGrid, 0); 

       
  return cGrid; 
}


vector<int_8> cellsAroundVertex(Vec3 &p, int layer, int i0=0, bool isVertex=false) {
  auto c = searchCellbyCoords(&p, i0, isVertex); 
}


vector<int_8> cellsAroundInterface(Grid* &a, int layer) {

  vector<int_8> tmp; 
  bool color[listCell.size()]; 
  layer = 0; 

  for (auto i = 0; i<listCell.size(); ++i) color[i] = false; 

  for (auto i1 = 0; i1 < a->listVertex.size(); ++i1) {
    auto v0 = a->otherVertex[i1]; 
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
  

#endif
