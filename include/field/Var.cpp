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
#include <string>
#include <grid2/Grid>
#include <memory>
#include <iomanip>
//#include <MathExp>
#include <linSolve/LinSys>
#include <field/Var>

void Var::solve(LinSys a) {

  // cout << "Size: " << a.b.size() << endl; //listBC.size() <<endl; // << " " << grid->listBNDR.size() << endl; 
  // Set BC; 
  a.x = &data;
  a.error = &error; 
  
  auto t = clock(); 
  // double d; 
  // for (auto row = 0; row < a.b.size(); ++row) {
  //   d = 1.0/a.A[row][row]; 
  //   a.A[row] = a.A[row]*d; 
  //   a.b[row] = a.b[row]*d; 
  // }
  //t = clock()-t; 
  // cout << " --> Normalization took "<< t/(double) CLOCKS_PER_SEC << " secs"<< endl; 
  // cout << a.A << endl; 
  //cout << a.b << endl; 
  a.setLimits(tol, itmax); 
  t = clock(); 
  if (solver.compare("BiCGSTAB") == 0) {
    a.BiCGSTAB(); 
  } else if (solver.compare("CG") == 0) {
    a.CG(); 
  } else {
    a.Gauss_Seidel(); 
  }
  t = clock() - t; 
  cout << " + " << name << " " << setw(10)<< data.min() << setw(1); 
  cout << ":" << setw(10) << data.max() ;
  cout <<setw(3)  <<" | "<< setw(12) << solver; 
  cout << setw(5) << " it: " << setw(6) << a.it << " tol: " << a.err; 
  cout << " t: "<< t/(double) CLOCKS_PER_SEC << " s"<< endl; 
}

// void Var::set(initializer_list<std::string> vname, std::string exp) {
//   typedef exprtk::symbol_table< double > symbol_table_t;
//   typedef exprtk::expression< double > expression_t;
//   typedef exprtk::parser< double > parser_t;
//   typedef exprtk::parser_error::type error_t;
  
//   int nVar = vname.size(); // size of the variables
//   vector< VecX<double > > phi(nVar); 
//   vector< double > expval(nVar); 

//   symbol_table_t symbol_table; 
//   symbol_table.add_constants(); 
  
//   VecX<double > d; 

//   int i = 0; 
//   for (auto v : vname) {
//     if (v.compare("x")) {
//       d = grid->getCoord(0); 
//     } else if (v.compare("y")) {
//       d = grid->getCoord(1); 
//     } else if (v.compare("z")) {
//       d = grid->getCoord(2); 
//     } else if (auto vr = grid->getVar(v)) {
//       d = vr->get(); 
//     } else {
//       cout << "Var::set. Variable name not recognized!"<<endl; 
//       exit(1); 
//     }
//     phi[i].data.assign(d.data.begin(), d.data.end());
//     d.clear(); 
//     symbol_table.add_variable(v, expval[i]); 
//     ++i;
//   }

//   expression_t expression; 
//   expression.register_symbol_table(symbol_table);       
//   parser_t parser; 

//   if (!parser.compile(exp, expression)) {
//     printf("Error: %s\tExpression: %s\n",
// 	   parser.error().c_str(),
// 	   exp.c_str());
      
//     for (std::size_t j = 0; j < parser.error_count(); ++j) {
//       error_t error = parser.get_error(j);
//       printf("Error: %02d Position: %02d Type: [%s] Msg: %s Expr: %s\n",
// 	     static_cast<int>(j),
// 	     static_cast<int>(error.token.position),
// 	     exprtk::parser_error::to_str(error.mode).c_str(),
// 	     error.diagnostic.c_str(),
// 	     exp.c_str());
//     }
//     exit(1);
//   }

//   for (auto i = 0; i < data.size(); ++i) {
//     for (auto j = 0; j < nVar; ++j) { 
//       expval[j] = phi[j].data.at(i);
//     }
//     data[i] =  expression.value();
//   } 

//   //data.data.assign(phi[i].begin(), phi[i].end());
  
// }
