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
#ifndef MATY
#define MATY

//-------------------------------------------
// Matrix storage scheme
//    Defines: basic matrix operations (not all)
//             Vec*Matrix addition
//-------------------------------------------

#include <VecX>

class Triplet {
public:
  int_8 row; 
  int_8 col; 
  double val; 
  Triplet(int i, int j, double a) { row = i; col = j; val = a;}

  // Operators
  // sum; 
}; 

  
class triLinSys {
public:
  vector<Triplet> vt; 
  vector<int_8> aij; 
  vector<double> val; 
  int_8 rank; 

  VecX<double> b;
  
  triLinSys() {rank = 0;}; 
  triLinSys(int_8 const &N) { rank = N;};
  
  void setMat() {
    if (vt.empty()) {
      if (rank > 0) {
	aij.assign(rank+1, rank+1);
	val.assign(rank+1, 0); 
      }
      return; 
    }
    // first sort a (1st col then row); 
    sort(vt.begin(), vt.end(), 
	 [](const Triplet & a, const Triplet & b) -> bool
	 { 
	   return (a.row < b.row) || ((a.row == b.row) && (a.col < b.col)); 
	 });
    rank = max(rank, vt.rbegin()->row + 1); 
    
    cout << "Rank: " << rank << endl; 
    cout << " ------------------ " <<endl; 
    
    for (auto c : vt) { 
      cout << c.row << " " << c.col << " " << c.val <<endl; 
    }
    cout << " ------------------ " <<endl; 

    aij.reserve(30*rank); 
    val.reserve(30*rank); 

    aij.assign(rank+1, rank+1);
    val.assign(rank+1, 0); 

    cout << " Size: " << aij.size() << " " << val.size() <<endl; 

    int_8 oldrow = vt.begin()->row;

    for (auto c : vt) {
      // c.row == c.col
      if (c.row == c.col) {
	val[c.row] += c.val;
	continue; 
      }
      if (c.row != oldrow) {
	aij[c.row+1] = aij[c.row]; 
	oldrow = c.row; 
      }

      // check if col exist at the last location
      if (*(aij.rbegin()) == c.col) {
	(*val.rbegin()) += c.val; 
	continue; 
      }

      aij.emplace_back(c.col); 
      val.emplace_back(c.val); 
      aij[c.row+1]++; 

    }
    vt.clear(); 

  } 


};



#endif
