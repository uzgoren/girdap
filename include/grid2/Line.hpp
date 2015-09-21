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
#ifndef LINE
#define LINE

class Line: public Cell {
protected:
  const int_2 map[2] = {1, 0};  
public:
  Line(initializer_list<int_8> l):Cell(l) {
    setType(3); 
  }
  Vec3 vol() { 
    Vec3 a = edge(); 
    return Vec3(-a.data[1], a.data[0], 0); 
  }
  Vec3 edge(int_2 i=0) { return (**getVertex(1)) - (**getVertex(0));} 
  //  Vec3 grad(shared_ptr<Var > phi);
  //void getBCPhi(shared_ptr<Var> phi, double &a, double &b);  
};


#endif
