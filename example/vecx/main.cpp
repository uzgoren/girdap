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

#include <iostream>
#include <fstream>
#include <girdap>


int main() {
  cout << "VecX tests: "<< endl;
  cout << "  + Constructor: empty VecX constructor: " << endl;
  VecX<double> vempty;
  cout << "    + vempty is with size " << vempty.size() << endl << endl;

  cout << "  + Resize: vempty to 7: " << endl;
  vempty.resize(7);
  cout << "    + vempty is with size " << vempty.size() << endl;
  cout << "    + v[0] = " << vempty[0] << ", v[5] = " << vempty[5] << endl<<endl;

  cout << "  + Assign: vempty to ones(6,1): " << endl;
  vempty.assign(6, 1);
  cout << "    + vempty is with size " << vempty.size() << endl;
  cout << "    + v[0] = " << vempty[0] << ", v[5] = " << vempty[5] << endl<<endl;
  
  cout << "  + Constructor: init 10 elements uninit " << endl;
  VecX<double> v10(10);
  cout << "    + v10 is with size " << v10.size() << endl;
  cout << "    + v[0] = " << v10[0] << ", v[5] = " << v10[5] << endl<<endl;

  cout << "  + Constructor: init 5 elements with int(0) " << endl;
  VecX<int_8> v5(5, 0);
  cout << "    + v5 is with size " << v5.size() << endl;
  cout << "    + v[0] = " << v5[0] << ", v[4] = " << v5[4] << endl<<endl;

  cout << "  + Constructor: init 3 elements with int(0) " << endl;
  VecX<double> v3(3, 3.14);
  cout << "    + v3 is with size " << v3.size() << endl;
  cout << "    + v[0] = " << v3[0] << ", v[4] = " << v3[2] << endl<<endl;

  cout << "  + Constructor: init 3 elements with list " << endl;
  VecX<double> vlist01{1.1, 3.14, -2.1};
  cout << "    + vlist01 is with size " << vlist01.size() << endl;
  cout << "    + v = {" << vlist01[0] << ", " << vlist01[1];
  cout << ", " << vlist01[2] << "}" << endl<<endl;

  cout << "  + Assign: vempty = vlist01: " << endl;
  vempty.assign(vlist01);
  cout << "    + vempty is with size " << vempty.size() << endl;
  cout << "    + v = {" << vempty[0] << ", " << vempty[1];
  cout << ", " << vempty[2] << "}" << endl;
  cout << "    + vlist01 is with size " << vlist01.size() << endl;
  cout << "    + v = {" << vlist01[0] << ", " << vlist01[1];
  cout << ", " << vlist01[2] << "}" << endl<<endl;

  vlist01 = { 9.1, 7.5, 7.2, -6.3, 4.1, 0, -1 };
  
  cout << "  + Swap: vempty = vlist01: " << endl;
  vempty.swap(vlist01);
  cout << "    + vempty = {"; for (auto v : vempty) cout << v << ", " ; cout << "}" << endl;
  cout << "    + vlist01 = {"; for (auto v : vlist01) cout << v << ", " ; cout << "}" << endl << endl;

  vlist01.assign({1.0, 2.0, 3.0}, 2);
  cout << "  + Assign: vlist01 @ 2 " << endl;
  cout << "    + vlist01 = {"; for (auto v : vlist01) cout << v << ", " ; cout << "}" << endl << endl;

  cout << "  + Min/Max " << endl;
  cout << "    + vlist01 = <" << vlist01.min() << ", " << vlist01.max() <<">" <<endl;
  cout << "    + vempty = <" << vempty.min() << ", " << vempty.max() <<">"<<endl<<endl;

  cout << " + Empty: "<< endl;
  vempty.empty();
  cout << "   + vempty, size = "<< vempty.size() << " capacity =" << vempty.capacity() << endl<<endl;

  cout << " + Clear: "<< endl;
  vlist01.clear();
  cout << "   + vlist01, size = "<< vlist01.size()  << " capacity =" << vlist01.capacity() << endl<<endl;

  cout << " + push_pair: " << endl;
  for (auto i = 0; i<10; ++i)
    vempty.push_pair(i, (double)2*i);
  cout << "   + Vempty: " << vempty <<endl<<endl; 

  cout << " + zeros: " << endl;
  vempty.zeros(10);
  cout << "   + vempty: " << vempty <<endl<<endl;; 

  cout << " + ones: " << endl;
  vempty.ones(10);
  cout << "   + vempty: " << vempty <<endl<<endl;; 

  cout << " + linear: with N" << endl;
  vempty.linear(0, 1, (int_8)10);
  cout << "   + vempty: " << vempty <<endl<<endl;;
  
  cout << " + linear: with delta" << endl;
  vempty.linear(1.0, 0.1, 2.0);
  cout << "   + vempty: " << vempty <<endl<<endl;; 

  cout << " + push_pair: vectors " << endl;
  vector<int_8> a = {1,3,5,7, 12}; vector<double> b = {2, 4, 6, 8, 10};
  vempty.push_pair(a, b); 
  cout << "   + Vempty: " << vempty <<endl<<endl;

  cout << " + push_pair: initlist " << endl;
  vempty.push_pair({2,4,6,8, 14}, {-1, -3, -5, -7, -10}); 
  cout << "   + Vempty: " << vempty <<endl<<endl;

  cout << " + abs: friend - element by element " << endl;
  vempty.push_pair({2,4,6,8, 14}, {-1, -3, -5, -7, -10}); 
  cout << "   + Vempty: " << abs(vempty) <<endl<<endl; 

  cout << " + abs: " <<endl;
  cout << "   + |vempty|: " << vempty.abs() <<endl<<endl; 

  cout << " + Vector of Vec3: " << endl;
  VecX<VecX<double> > vv0 = { {0,1,1}, {0,2,2} };
  cout << "   + vv0 = " << vv0 << endl << endl;

  cout << " + Assignment =  " << endl;
  VecX<VecX<double> > vv1;
  vv1 = vv0;
  cout << "   + vv1 = " << vv1 << endl;
  cout << "   + vv0 = " << vv0 << endl<<endl;

  cout << " + Assignment = initlist" << endl; 
  VecX<double> vv2;
  vv2 = {3, 4, 5, -3, -2, 9.8};
  cout << "   + vv2 = " << vv2 << endl<<endl; 

  cout << " + Dot product " << endl;
  VecX<double> vv3; vv3.linear(2, 3, (int_8)6); 
  cout << "   + vv2*vv3 = "<< vv2*vv3 <<endl;

  cout << " + arithmetic " << endl; 
  
  return 0; 
};
