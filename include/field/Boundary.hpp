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
#ifndef BOUNDARY
#define BOUNDARY


class Boundary {
public:
  int type; // 0: Dirichlet and 1: Neumann
  double a_val, b_val; // Linearized a_val*phi_p + b_val; 
 
  Boundary():type(1), b_val(0), a_val(0) {}
  Boundary(int t) :type(t), b_val(0), a_val(0) {}
  Boundary(int t, double b) : type(t), b_val(b), a_val(0) {}
  Boundary(int t, double b, double a) : type(t), b_val(b), a_val(a){}

  void setBC() {
    type = 1; a_val = 0; b_val = 0; 
  };

  void setBC(int t) {
    type = t; b_val = 0; a_val = 0; 
  }; 
 
  void setBC(int t, double b, double a=0) {
    type = t; b_val = b; a_val = a; 
  };
     
};


#endif
