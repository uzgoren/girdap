#ifndef IOGRID
#define IOGRID

void writeVTK(string name, initializer_list<string > v = {"all"}) {

  // Image 

  bool isAll = false; 
  vector<bool> lvar; 
  if (listVar.size() > 0) {
    isAll = (v.size() == 1 && *v.begin() == "all"); 
    lvar.resize(listVar.size(), false); 
    for (auto i = 0; i < listVar.size(); ++i) {
      if (!isAll) {
	for (auto j = v.begin(); j != v.end(); ++j) 
	  if (*j == listVar[i]->name || (*j == "vel" && listVar[i]->name == "u") 
	      || (*j+"x" == listVar[i]->name)) {
	    lvar[i] = true; 
	    break;
	  }
      } else {
	lvar[i] = true;
	if (listVar[i]->isVec) {
	  if (listVar[i]->name == "v" || listVar[i]->name == "w" 
	      || listVar[i]->name.substr(listVar[i]->name.length()-1) == "y"
	      || listVar[i]->name.substr(listVar[i]->name.length()-1) == "z")
	    lvar[i] = false; 
	}
      }	
    }
  }
       
  ofstream out;
  out.open(name+std::to_string(filecnt++)+".vtk"); 


  int_8 nCell = listCell.size();
  int t; // = listCell[0]->node.size()+1;
  out << "# vtk DataFile Version 2.0" << endl; 
    out << "Unstructure Grid" << endl; 
    out << "ASCII"<< endl; 
    out << endl << "DATASET UNSTRUCTURED_GRID"<< endl; 
    out << "POINTS " << listVertex.size() << " float" << endl; 
    for (auto p = vbegin(); p != vend(); ++p) {
      out << **p <<  endl ; 
    }
    if (nCell > 0) {
      t = listCell[0]->node.size()+1;
      out << endl << "CELLS "<< nCell << " " << nCell*t << endl; 
      for (auto c = cbegin(); c != cend(); ++c) {
	if ((*c)->isAlive) out << *c ; 
      }
      out << endl << "CELL_TYPES " << nCell <<endl; 
      for (auto c = cbegin(); c != cend(); ++c) {
	if ((*c)->isAlive) out << (*c)->getType() << endl; 
      }
    } else {
      out << endl << "CELLS "<< listVertex.size() << " " << listVertex.size()*2 << endl; 
      for (auto i = 0; i < listVertex.size(); ++i) 
        out << "1 "<< i<< endl;     
      out << endl << "CELL_TYPES " << listVertex.size() <<endl; 
      for (auto i = 0; i < listVertex.size(); ++i) out << "1 "<< endl; 
    }
    out << endl << "POINT_DATA " << listVertex.size() << endl;
    for (auto i = 0; i < listVar.size(); ++i) {
      if (!lvar[i]) continue; 
      auto var = listVar[i]; 
      if (listVar[i]->isVec) continue; 
      
      out << "SCALARS " << var->name << " float 1"<<endl; 
      out << "LOOKUP_TABLE default"<<endl; auto icnt = 0; 
      for (auto v: listVertex) {
	//auto val = v->phi(var).eval(var); 
	double val = 0; 
	if (var->loc == 0) 
	  val=v->evalPhi(var); 
	else if (var->loc == 1) 
	  val=var->get(v->id); 
	out << ((abs(val) < 1e-10) ? 0 : val) << endl; 
      }
      out << endl; 
    }
    if (listCell.size() > 0) {
      for (auto i = 0; i < listVar.size(); ++i) {
	if (!lvar[i]) continue; 
	auto var = listVar[i]; 
	if (!var->isVec) {continue;}
	std::string name = var->name; 
	if (name != "u") {
	  if (name.substr(name.length()-1) != "x") {continue;}
	  else { name.pop_back(); }
	}
	auto v = getVecVertex(name);
	if (name == "u") name = "vel"; 
	out << "VECTORS "+name+" float" << endl; 
	for (auto i = 0; i<v.size(); ++i) out << float(v[i][0]) << " " << float(v[i][1]) << " " << float(v[i][2]) << endl;
      }
    }
    out<<endl; 	  
    out.close(); 
  }


void writePast(string name, shared_ptr<Var > var) {
  auto vel = getVecVertex("u");
  // Image        
  ofstream out;
  out.open(name+std::to_string(filecnt++)+".vtk"); 

  int_8 nCell = listCell.size();
  int t; // = listCell[0]->node.size()+1;
  out << "# vtk DataFile Version 2.0" << endl; 
  out << "Unstructure Grid" << endl; 
  out << "ASCII"<< endl; 
  out << endl << "DATASET UNSTRUCTURED_GRID"<< endl; 
  out << "POINTS " << listVertex.size() << " float" << endl; 
  for (auto p = vbegin(); p != vend(); ++p) {
    out << **p-vel[(*p)->id]*0.05 <<  endl ; 
  }
  if (nCell > 0) {
    t = listCell[0]->node.size()+1;
    out << endl << "CELLS "<< nCell << " " << nCell*t << endl; 
    for (auto c = cbegin(); c != cend(); ++c) {
      if ((*c)->isAlive) out << *c ; 
    }
    out << endl << "CELL_TYPES " << nCell <<endl; 
    for (auto c = cbegin(); c != cend(); ++c) {
      if ((*c)->isAlive) out << (*c)->getType() << endl; 
    }
  } else {
    out << endl << "CELLS "<< listVertex.size() << " " << listVertex.size()*2 << endl; 
    for (auto i = 0; i < listVertex.size(); ++i) 
      out << "1 "<< i<< endl;     
    out << endl << "CELL_TYPES " << listVertex.size() <<endl; 
      for (auto i = 0; i < listVertex.size(); ++i) out << "1 "<< endl; 
  }
  out << endl << "POINT_DATA " << listVertex.size() << endl;
  
  out << "SCALARS " << var->name << " float 1"<<endl; 
  out << "LOOKUP_TABLE default"<<endl; auto icnt = 0; 
  for (auto v: listVertex) {
    //auto val = v->phi(var).eval(var); 
    double val = 0; 
    if (var->loc == 0) {
      auto x = Vec3(*v - vel[v->id]*dt*0.5); 
      val=interp(var, x, v->id); //v->evalPhi(var, &x); 
    } else if (var->loc == 1) 
      val=var->get(v->id); 
    out << ((abs(val) < 1e-10) ? 0 : val) << endl; 
  }
  out << endl; 
  
  out<<endl; 	  
  out.close(); 
}


void writeInterp(string name) {
  ofstream out;
  out.open(name+std::to_string(filecnt)+".vtk"); 
  int_8 nCell = listVertex.size();
 
  out << "# vtk DataFile Version 2.0" << endl; 
  out << "Unstructure Grid" << endl; 
  out << "ASCII"<< endl; 
  out << endl << "DATASET UNSTRUCTURED_GRID"<< endl; 
  out << "POINTS " << 4*nCell << " float" << endl;
  for (auto v : listVertex) {
    out << int3D.linFunc(Vec3(0,0,0), v->xcoef) << " "; 
    out << int3D.linFunc(Vec3(0,0,0), v->ycoef) << " "; 
    out << int3D.linFunc(Vec3(0,0,0), v->zcoef) << endl; 
    out << int3D.linFunc(Vec3(1,0,0), v->xcoef) << " "; 
    out << int3D.linFunc(Vec3(1,0,0), v->ycoef) << " "; 
    out << int3D.linFunc(Vec3(1,0,0), v->zcoef) << endl; 
    out << int3D.linFunc(Vec3(1,1,0), v->xcoef) << " "; 
    out << int3D.linFunc(Vec3(1,1,0), v->ycoef) << " "; 
    out << int3D.linFunc(Vec3(1,1,0), v->zcoef) << endl; 
    out << int3D.linFunc(Vec3(0,1,0), v->xcoef) << " "; 
    out << int3D.linFunc(Vec3(0,1,0), v->ycoef) << " "; 
    out << int3D.linFunc(Vec3(0,1,0), v->zcoef) << endl;     
  }
  auto icnt = 0; 
  out << endl << "CELLS "<< nCell << " " << nCell*5 << endl; 
  for (auto v : listVertex) { 
    out << "4 " << icnt << " " << icnt+1 << " " << icnt+2 << " " << icnt+3 << endl; 
    icnt += 4; 
  }
  out << endl << "CELL_TYPES " << nCell <<endl; 
  for (auto v : listVertex) {
    out << "9 "<< endl; 
  }
  out.close();  

}


  friend ostream &operator<<(ostream &out, Grid* const &a) {
    if (a == NULL) {out << " "; return out; }
    int_8 nCell = a->listCell.size(); //a->nCellSize = 0; 
    // for (auto i =0; i < a->listCell.size(); ++i) {
    //   if (a->listCell[i]->isAlive) {
    // 	nCell++; 
    // 	a->nCellSize += a->listCell[i]->node.size() + 1; 
    //   }
    // }
   int t = a->listCell[0]->node.size()+1;
    out << "# vtk DataFile Version 2.0" << endl; 
    out << "Unstructure Grid" << endl; 
    out << "ASCII"<< endl; 
    out << endl << "DATASET UNSTRUCTURED_GRID"<< endl; 
    out << "POINTS " << a->listVertex.size() << " float" << endl; 
    for (auto p = a->vbegin(); p != a->vend(); ++p) {
      out << **p <<  endl ; 
    }
    if (nCell > 0) {
      out << endl << "CELLS "<< nCell << " " << nCell*t << endl; 
      for (auto c = a->cbegin(); c != a->cend(); ++c) {
	if ((*c)->isAlive) out << *c ; 
      }
      out << endl << "CELL_TYPES " << nCell <<endl; 
      for (auto c = a->cbegin(); c != a->cend(); ++c) {
	if ((*c)->isAlive) out << (*c)->getType() << endl; 
      }
    } else {
      out << endl << "CELLS "<< a->listVertex.size() << " " << a->listVertex.size()*2 << endl; 
      for (auto i = 0; i < a->listVertex.size(); ++i) 
        out << "1 "<< i<< endl; 
    
      out << endl << "CELL_TYPES " << a->listVertex.size() <<endl; 
      for (auto i = 0; i < a->listVertex.size(); ++i) out << "1 "<< endl; 
    }
    out << endl << "POINT_DATA " << a->listVertex.size() << endl;
    for (auto i = 0; i < a->listVar.size(); ++i) {
      auto var = a->listVar[i]; 
      if (a->listVar[i]->isVec) continue; 
      // cout << var->name << " " << var->loc << " "<< var->data.size() << " " << a->listCell.size() << endl; 
      // if (var->data.size() != a->listCell.size()) continue; 
      out << "SCALARS " << var->name << " float 1"<<endl; 
      out << "LOOKUP_TABLE default"<<endl; auto icnt = 0; 
      for (auto v: a->listVertex) {
	//auto val = v->phi(var).eval(var); 
	double val = 0; 
	if (var->loc == 0) 
	  val=v->evalPhi(var); 
	else if (var->loc == 1) 
	  val=var->get(v->id); 
	out << ((abs(val) < 1e-10) ? 0 : val) << endl; 
      }
      //auto d = a->getPhiVertex(v);
      //for (auto i = 0; i<d.size(); ++i) out << d[i] << endl; 
      //      for (auto s = v->data.begin(); s != v->data.end(); ++s) {
      //	out << *s << " "; ++icnt; }
      out << endl; 
    }
    if (a->listCell.size() > 0) {
      for (auto i = 0; i < a->listVar.size(); ++i) {
	auto var = a->listVar[i]; 
	if (!var->isVec) {continue;}
	std::string name = var->name; 
	if (name != "u") {
	  if (name.substr(name.length()-1) != "x") {continue;}
	  else {name.pop_back(); }
	}
	
	auto v = a->getVecVertex(name);
	if (name == "u") name = "vel"; 
	out << "VECTORS "+name+" float" << endl; 
	for (auto i = 0; i<v.size(); ++i) out << float(v[i][0]) << " " << float(v[i][1]) << " " << float(v[i][2]) << endl;
      }
    }
    // extra -remove
    // out << endl << "CELL_DATA " << nCell << endl;
    // for (auto j = 0; j < 3; ++j) {
    //   out << "SCALARS adapt"<<j<<" int 1"<<endl; 
    //   out << "LOOKUP_TABLE default"<<endl; auto icnt = 0; 
    //   for (auto i = 0; i < a->listCell.size(); ++i)
    // 	if (a->listCell[i]->isAlive) out << a->listCell[i]->adapt[j] << endl;
    // }
    // for (auto j = 0; j < 3; ++j) {
    //   out << "SCALARS level"<<j<<" int 1"<<endl; 
    //   out << "LOOKUP_TABLE default"<<endl; auto icnt = 0; 
    //   for (auto i = 0; i < a->listCell.size(); ++i)
    // 	if(a->listCell[i]->isAlive) out << a->listCell[i]->level[j] << endl;
    // } 
    // out << "SCALARS masterx int 1"<<endl; 
    // out << "LOOKUP_TABLE default"<<endl; 
    // for (auto i = 0; i < a->listCell.size(); ++i)
    //   if(a->listCell[i]->isAlive) out << a->listCell[i]->masterx[a->listCell[i]->level[0]] << endl;
    // out << "SCALARS mastery int 1"<<endl; 
    // out << "LOOKUP_TABLE default"<<endl; 
    // for (auto i = 0; i < a->listCell.size(); ++i)
    //   if(a->listCell[i]->isAlive) out << a->listCell[i]->mastery[a->listCell[i]->level[1]] << endl;


    out<<endl; 
    return out;
  }
 


  void writeAMRdata(std::string a = "amr.vtk") {
    auto nCell = listCell.size(); 
    ofstream out; 
    out.open(a); 
    out << "# vtk DataFile Version 2.0" << endl; 
    out << "Unstructure Grid" << endl; 
    out << "ASCII"<< endl; 
    out << endl << "DATASET UNSTRUCTURED_GRID"<< endl; 
    out << "POINTS " << listVertex.size() << " float" << endl; 
    for (auto p = vbegin(); p != vend(); ++p) {
      out << **p <<  endl ; 
    }
    out << endl << "CELLS "<< nCell << " " << nCell*5 << endl; 
    for (auto c: listCell) 
      if (c->isAlive) out << c; 
  
    out << endl << "CELL_TYPES " << nCell <<endl; 
    for (auto c: listCell) 
      if (c->isAlive) out << c->getType() << endl; 

    out << endl << "CELL_DATA " << nCell << endl; 
    out << "SCALARS master_x float 1"<<endl; 
    out << "LOOKUP_TABLE default"<<endl; 
    for (auto c: listCell) 
      if (c->isAlive) out << (c->masterx[c->level[0]]?1:0) << endl; 
    out << endl; 

    out << "SCALARS master_y float 1"<<endl; 
    out << "LOOKUP_TABLE default"<<endl; 
    for (auto c: listCell) 
      if (c->isAlive) out << (c->mastery[c->level[1]]?1:0) << endl; 

    out <<endl;
    out.close(); 
  }

  
  void writeFace(std::string a = "face.vtk") { 
    ofstream out; 
    out.open(a);
    out << "# vtk DataFile Version 2.0" << endl; 
    out << "Unstructure Grid" << endl; 
    out << "ASCII"<< endl; 
    out << endl << "DATASET UNSTRUCTURED_GRID"<< endl; 
    out << "POINTS " << listVertex.size() << " float" << endl; 
    for (auto p = vbegin(); p != vend(); ++p) {
      out << **p <<  endl ; 
    }
    out << endl << "CELLS "<< nFace << " " << 3*(nFace) << endl; 
    for (auto f: listFace) 
      if (f) out << f; 
  
    out << endl << "CELL_TYPES " << nFace <<endl; 
    for (auto f: listFace) 
      if (f) out << f->getType() << endl; 

    out << endl << "CELL_DATA " << nFace << endl; 
    for (auto i = 3; i < listVar.size(); ++i) {
      auto var = listVar[i]; 
      if (var->data.data.size() != nFace) continue; 
      out << "SCALARS " << var->name << " float 1"<<endl; 
      out << "LOOKUP_TABLE default"<<endl; 
      out << var->data << endl; 
      out << endl; 
    }
    
    // out << "SCALARS next" << " int 1"<<endl; 
   
    // for (auto f: listFace) 
    //   if (f) out << f->next << endl; 

    // out << "SCALARS prev" << " int 1"<<endl; 
    // out << "LOOKUP_TABLE default"<<endl; 
    // for (auto f: listFace) 
    //   if (f) out << f->prev << endl; 
  }


#endif
