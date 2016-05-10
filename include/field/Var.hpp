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
#ifndef VAR
#define VAR

#include <memory>
#include <linSolve/LinSys> 
#include <field/Boundary>

class Grid;
//template<typename T> class Scheme;

class Var {
public:
  std::string name;
  VecX<double> data;
  VecX<double> laplace, div; 
  //vector<Scheme<Vec3> > grad;
  VecX<Vec3> grad; 
  VecX<double> error; 
  //  LinSys mat; 
  Grid *grid; 
  int_2 loc;  // 0 for cell center; 1 for vertex; 2 for faces [faces are difficult]
  bool isVec; 
  double dt; 
  vector<shared_ptr<Boundary> > listBC; 


  string solver; 
  int itmax; 
  double tol; 

  
  Var(std::string a): name(a), isVec(false) { dt = 1.0; loc = 0; solver = "BiCGSTAB"; itmax = -1; tol = -1; }
  Var(std::string a, int b) :Var(a) { if (b<0 && b>2) b=0; loc = b; }

  VecX<double> get() {return data;}
  double get(int_8 i); 
  double* operator[](int_8 const &index); 
  
  void changeName(std::string a) { name = a;}
  void set(double a) { data.data.assign(data.size(), a); }
  void set(VecX<double > a) { data.data.assign(a.data.begin(), a.data.end()); }
  void set(int_8 i, double a) { 
    if (i >=0 && i < data.data.size()) data.data[i] = a; 
    else {
      cout << "Var " << name << " failed to set " << i << " to " << a << endl; 
      exit(1); 
    }
  }
  void set(initializer_list<std::string> vname, std::string exp); 

  void resize(int_8 n) {
    if (n <=  0) data.data.clear(); 
    else data.data.resize(n); 
  }

  void initBC() {
    listBC.resize(6); 
    auto icnt =0; 
    for (auto i = 0; i<listBC.size(); ++i) {
      listBC[i] = shared_ptr<Boundary>(new Boundary()); 
     }
   }

  void setBC(string id, string t, double b, double a=0) {    
    int i; int type; 
    if (id.compare("south") == 0) i = 0; 
    else if (id.compare("north") == 0) i = 1; 
    else if (id.compare("west") == 0) i = 2; 
    else if (id.compare("east") == 0) i = 3; 
    else return; 
    if (t.compare("val") == 0) type = 0; 
    else if (t.compare("grad") == 0) type = 1; 
    else return; 

    if (listBC[i]) {
      listBC[i]->setBC(type, b, a); 
    }
  }
    

  void solve(LinSys a);
};


#endif
