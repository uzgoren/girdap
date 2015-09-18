#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>
#include <time.h>
#include <iomanip>
#include <random>
#include <field/Var>
#include <grid2/Grid>


int main() {
  ofstream myfile; 
  ifstream input;
  int filecnt = 0; 

  std::default_random_engine gen;
  std::bernoulli_distribution incx(0.3);
  std::bernoulli_distribution incy(0.5);
  std::bernoulli_distribution decx(0.7);
  std::bernoulli_distribution decy(0.7);
  bool prod = false; 
  // GRID 
  Block2* grid = new Block2({0, 0, 0}, {1.0, 1.0, 0}, 5, 5);
  //  Block2* grid2 = new Block2({0, 0, 0}, {1.0, 1.0, 0}, 3, 3);
  // grid->adaptCriteria(); 

  auto u = grid->getVar("u"); 
  auto v = grid->getVar("v"); 

  for (auto pass = 0; pass < 4; ++pass) {
    std::string flname="op"+std::to_string(pass)+".dat"; 
    if (prod) {
      myfile.open(flname); 
      for (auto c : grid->listCell) {
	c->adapt[0] = 0; 
	c->adapt[1] = 0; 
	if (incx(gen) && c->level[0] != grid->levelHighBound[0]) { c->adapt[0]=1;}
	if (incy(gen) && c->level[1] != grid->levelHighBound[1]) { c->adapt[1]=1;}
	if (c->adapt[0] == 0 && c->adapt[1] == 0) {
	  if (decx(gen) && c->level[0] != grid->levelLowBound[0]) { c->adapt[0]=-1;}
	  if (decy(gen) && c->level[1] != grid->levelLowBound[1]) { c->adapt[1]=-1;}
	}
	myfile << c->adapt[0] << " " << c->adapt[1] << endl; 
      }
      myfile.close();  
    } else {
      input.open(flname); 
      if (input.is_open()) {
	for (auto c : grid->listCell) {
	  input >> c->adapt[0];
	  input >> c->adapt[1];
	  // c->adapt[0] = max(short(0), c->adapt[0]); 
	  // c->adapt[1] = max(short(0), c->adapt[1]); 
	}
      input.close();
      } else {
	cout << "cannot read data"<< endl; 
	exit(1); 
      }
    }

    cout << filecnt << " pass: " << pass << " random " << endl; 
    flname="grid"+std::to_string(filecnt++)+".vtk"; 
    myfile.open(flname); 
    myfile << grid << endl;
    myfile.close();  


    for (auto i = 0; i < grid->listCell.size(); ++i) { 
      grid->listCell[i]->checkNgbrLevel(0); 
      grid->listCell[i]->checkNgbrLevel(1);
      grid->listCell[i]->checkNgbrLevel(2);
    } 

    cout << filecnt << " pass: " << pass << " corrected " << endl; 
    flname="grid"+std::to_string(filecnt++)+".vtk"; 
    myfile.open(flname); 
    myfile << grid << endl;
    myfile.close();  
    

    grid->refine2();
    // auto lmin = **(std::min_element(&(grid->levelMin), &(grid->levelMin)+3)); 
    // auto lmax = **(std::max_element(&(grid->levelMax), &(grid->levelMax)+3)); 
    // for (auto j = 0; j < 2; ++j) {
    //   for (auto l = grid->levelMin[j]; l < grid->levelMax[j]+1; ++l) {
    // 	auto nCell = grid->listCell.size();
    // 	for (auto i = 0; i < nCell; ++i) {
    // 	  auto c = grid->listCell[i]; 
    // 	  if (c->level[j] != l) continue;
    // 	  if (c->adapt[j] == 0) continue;
    // 	  c->refine(j); 
    // 	  // if (j == 0) 
    // 	  //   c->convertToSimpleBlock({2,1});
    // 	  // else 
    // 	  //   c->convertToSimpleBlock({1,2});
    // 	  // c->adapt[j]--; c->level[j]++;
    // 	  // for (auto k=0; k< 2; ++k) {
    // 	  //   grid->listCell[k]->adapt.assign(c->adapt.begin(), c->adapt.end()); 
    // 	  //   grid->listCell[k]->level.assign(c->level.begin(), c->level.end()); 
    // 	  // } 
    // 	}
    // 	grid->setCurrentLevels();
    // 	cout << "("<<grid->levelMin[0] << "-"<< grid->levelMax[0] <<")-"; 
    // 	cout << "("<<grid->levelMin[1] << "-"<< grid->levelMax[1] <<")"; 
    // 	cout << filecnt << " pass: " << pass << " l: " << l << " dir: "<< j << endl; 
    // 	for (auto i = 0; i<grid->listCell.size(); ++i) {
    // 	  u->set(i, grid->listCell[i]->adapt[0]); 
    // 	  v->set(i, grid->listCell[i]->adapt[1]);
    // 	} 
    // 	flname="grid"+std::to_string(filecnt++)+".vtk"; 
    // 	myfile.open(flname); 
    // 	myfile << grid << endl;
    // 	myfile.close();  
    //   } 
    // }
    flname="grid"+std::to_string(filecnt++)+".vtk"; 
    myfile.open(flname); 
    myfile << grid << endl;
    myfile.close();      
  }
    
    string flname="grid"+std::to_string(filecnt++)+".vtk"; 
    myfile.open(flname); 
    myfile << grid << endl;
    myfile.close();  
    // flname="grid"+std::to_string(filecnt++)+".vtk"; 
    // myfile.open(flname); 
    // myfile << grid2 << endl; 
    // myfile.close(); 

  delete(grid); 

  return 0; 
}
