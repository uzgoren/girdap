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
#include <math.h>
#include <fstream>
#include <sstream>
#include <time.h>
#include <iomanip>

#include <field/Var>
#include <grid2/Grid>


double heavy(double x, double e) {
  if (x/e <= -1)
    return 0; 
  else if (x/e <= 0)
    return 0.5 + x/e + 0.5*pow(x/e,2); 
  else if (x/e <= 1)
    return 0.5 + x/e - 0.5*pow(x/e,2); 
  else 
   return 1; 
}

double heavy(Vec3 x, Vec3 xp, double e) {
  return heavy((x - xp).abs(), e); 
}

int main(int argc,char *argv[]) {

  string casedir = "./"; 
  if (argc == 2) casedir = argv[1]; 
  if (casedir.back() != '/') casedir = casedir + "/"; 

  int_8 Nx=20, Ny=20; 
  int_2 lx=2, ly=2; 
  int intmethod=2; 
  double intscheme[3]={0.5, 0, 0.5}; 
  int rk0 = 2; 
  bool flux[3] = {false, true, false}; 
  double cfl = 0.5;   
  string veltype = "deform"; 
  bool isVelTime = false; 
  bool isVelSpace = false; 
  double u0, v0; 
  string geotype = "circle"; 
  double xm=0.3, xp=0.5, ym=0.55, yp=0.75, rp=0.15; 
  double endTime = 8; //dt*50; 
  double freq = 1; 
  string cname = "ls";
  auto writeTime = 0.005;   

  // read file; 
  std::ifstream infile(casedir+"case.txt"); 
  string name; double val; 
  while (infile >> name) { 
    //cout << name << endl; 
    if (name.compare("Nx")==0) {infile >> Nx;  continue;}
    else if (name.compare("Ny")==0) { infile >> Ny; continue;}
    else if (name.compare("lx")==0) { infile >> lx; continue;}
    else if (name.compare("ly")==0) { infile >> ly; continue;}
    else if (name.compare("interp")==0) { infile >> intmethod; continue;}
    else if (name.compare("scheme")==0) { infile >> intscheme[0] >> intscheme[1] >> intscheme[2]; continue;}
    else if (name.compare("rk")==0) { infile >> rk0; continue;}
    else if (name.compare("flux")==0) { infile >> flux[0] >> flux[1] >> flux[2]; continue;}
    else if (name.compare("cfl")==0) {infile >> cfl; continue;} 
    else if (name.compare("freq")==0) {infile >> freq; continue;} 
    else if (name.compare("endTime")==0) {infile >> endTime; continue;} 
    else if (name.compare("writeTime")==0) {infile >> writeTime; continue;}
    else if (name.compare("case")==0) {infile >> cname; continue;}
    else if (name.compare("veltype")==0) {
      infile >> veltype; 
      if (veltype.compare("-")==0) {
	infile >> u0 >> v0; continue; 
      } else if (veltype.compare("t")==0) {
	isVelTime =true;
	infile >> u0 >> v0; 
	continue; 
      } else if (veltype.compare("x")==0) {
	isVelSpace = true; 
	infile >> u0 >> v0; 
      } else if (veltype.compare("xt")==0) {
	isVelTime = true; 
	isVelSpace = true; 
	infile >> u0 >> v0; 
	continue;
      }     
    } else if (name.compare("geo")==0) {
      infile >> geotype; 
      if (geotype.compare("bar")==0) {	
	infile >> xp >> xm; 
	continue; 
      } else if(geotype.compare("rect")==0) {
	infile >> xp >> xm; 
	infile >> yp >> ym; 
	continue; 
      } else if (geotype.compare("circle")==0) {
	infile >> xp >> yp >> rp; 
	continue; 
      }
    }    
  }

  std::ofstream outcase;
  outcase.open(casedir+"run.txt"); 
  outcase << "Name = " << cname << endl; 
  outcase << "N = " << Nx << "x" << Ny << endl; 
  outcase << "levels = " << lx <<"x"<<ly<< endl; 
  if (intmethod == 0) outcase << "FOU" << endl; 
  else if (intmethod == 1) outcase << "Trilinear" << endl;  
  else outcase << "Grad" << endl; 
  outcase << "n=" << intscheme[0] << ", n+1/2="<< intscheme[1] << ", n+1="<< intscheme[2] << endl; 
  outcase << "rk=" << rk0 << endl;
  outcase << "flux calc: " << (flux[0] ? "phi " : "") << (flux[1] ? "vel " : "") << (flux[2] ? " normal" : "")<< endl; 
  outcase << "cfl = " << cfl << endl; 
  outcase << "geo = " << geotype << " xm=" << xm << " ym=" << ym << " xp=" << xp << " yp=" << yp << " rp=" << rp << endl;  
  outcase << "vel = " << veltype << " u0=" << u0 << " v0=" << v0 << endl;
  outcase << "endTime = " << endTime << " writeTime=" <<writeTime << " freq=" << freq << endl; 
  outcase.close(); 
			
  
  auto dt = 0.5; 
  auto cput = clock(); int iter = 0; 
  Block2* grid = new Block2({0, 0, 0}, {1, 1, 0}, Nx, Ny);
  grid->levelHighBound[0] = lx; 
  grid->levelHighBound[1] = ly; 
  grid->cfl = cfl; 

  //Block2* grid = new Block2({0, 0, 0}, {0.02, 0.02, 0}, 10, 1);
  double time= 0; 

  // Field variables; 
  double rho=10000, cp=1000, k=1; 

  grid->addVar({"T", "ui", "vi"}); 
    
  auto T = grid->getVar("T");
  auto u = grid->getVar("u");
  auto v = grid->getVar("v");
  auto ui = grid->getVar("ui"); 
  auto vi = grid->getVar("vi"); 

  T->solver = "Gauss";   
  T->set(0); 
  double pi = 4.0*atan(1); 
  
  ui->set(u0); 
  vi->set(v0); 
  u->set(u0); 
  v->set(v0); 
  for (auto j = 0; j < grid->levelHighBound[0] + 1; ++j) {

    auto dx = 1/Nx/pow(2,grid->levelMax[0]-1); 
    auto e = 4*dx; 

    for (auto i = 0; i < grid->listCell.size(); ++i) {
      auto x = grid->listCell[i]->getCoord(); // - Vec3(0.5, 0.5); 
      if (veltype.find("x") != std::string::npos) {
	if (geotype.compare("bar") != 0) {
	  ui->set(i, -u0*2*sin(pi*x[1])*cos(pi*x[1])*sin(pi*x[0])*sin(pi*x[0]));
	  u->set(i, ui->get(i)); 
	  vi->set(i, -v0*2*sin(pi*x[0])*cos(pi*x[0])*sin(pi*x[1])*sin(pi*x[1])); 
	  v->set(i, vi->get(i)); 
	} else {
	  ui->set(i, u0*(1 + sin(pi*x[0]))); 
	  u->set(i, ui->get(i));
	  vi->set(i, 0); 
	  v->set(i, 0); 
	}
      }
	
      if (geotype.compare("bar") == 0) 
	T->set(i, heavy(x[0]-xp, e) - heavy(x[0]-xm, e)); 
      else if (geotype.compare("rect") == 0) {
	T->set(i, (heavy(x[0]-xp, e) - heavy(x[0]-xm, e))*(heavy(x[1]-yp, e) - heavy(x[1]-ym, e))); 
	
      } else if (geotype.compare("circle") == 0) {
	auto r = (grid->listCell[i]->getCoord() - Vec3(xp, yp)).abs(); 
	//T->set(i, heavy(0.15-r, e)); 
	T->set(i, 1.0/(1.0 + exp(-2.0*80*(rp-r)))); 
      }
    }
    if (j == grid->levelHighBound[0]) break; 
    //    auto gt = grid->valGrad(T); 
    //grid->solBasedAdapt(gt); 
    grid->solBasedAdapt2(grid->getError(T));
    grid->adapt(); 
  }
  
  double mass0=0; double mass=0; 
  for (auto i=0; i < grid->listCell.size(); ++i) {
    mass0 += grid->listCell[i]->vol().abs()*T->get(i); 
  }

  grid->writeVTK(casedir+cname);


  int filecnt = 0; int it = 0, writeInt = 1; 

  T->setBC("west", "grad", 0);
  T->setBC("south", "grad", 0); 
  T->setBC("east", "grad", 0); 
  T->setBC("north", "grad", 0); 

  T->itmax = 1000; 
  T->tol = 1e-6;  

  auto writeCnt = writeTime; 
  while (time < endTime) {
    iter++; 
    for (auto i = 0; i < grid->listCell.size(); ++i) {
      auto x = grid->listCell[i]->getCoord(); // - Vec3(0.5, 0.5); 
      if (veltype.find("t") != std::string::npos) {
	u->set(i, ui->get(i)*cos(freq*pi*time/endTime));
	v->set(i, vi->get(i)*cos(freq*pi*time/endTime));
      }
    }

    grid->setDt(2.0); 
    if (writeCnt > 0 && writeCnt < grid->dt) grid->dt = writeCnt; 
    //    grid->writePast("v_past", T); 
 
    cout << setiosflags(ios::fixed) << setprecision(6); 
    cout << "------------- Processing TIME = "<< time << " ------------------"<<endl; 
    cout << "  Mass: " << mass << " Initial: " << mass0 << endl; 

    auto vel = grid->getVel();
    

    grid->advanceDiv(T, vel, rk0, intscheme, intmethod, flux); 
    cout << "Min T: " << T->data.min() << ", Max T: "<< T->data.max() << endl; 
    //if (time == 0) cin.ignore().get(); 
    //grid->writeVTK("den"); 
    //exit(1); 
    mass=0; double part=0; 
    for (auto i=0; i < grid->listCell.size(); ++i) {
      double vol = grid->listCell[i]->vol().abs();
      auto Tval = T->get(i); 
      if (Tval > 0.995) T->set(i, 1.0); 
      if (Tval < 0.005) T->set(i, 0.0); 
      mass += vol*T->get(i);
      if (Tval > 0) part += vol; 
    }
    //exit(1); 
    //    auto gt = grid->valGrad(T); 
    if (iter == 0 || iter % 5 == 0) {

      grid->solBasedAdapt2(grid->getError(T)); 

      grid->adapt(); 

    }

    time += grid->dt; 
    writeCnt -= grid->dt; 

    if (writeCnt <= 0 || time >= endTime) {
      grid->writeVTK(casedir+cname);
      // std::string flname="heat"+std::to_string(filecnt++)+".vtk"; 
      // myfile.open(flname); 
      // myfile << grid << endl;
      // myfile.close();   
      writeCnt = writeTime; 
      //cin.ignore().get(); 
     } 
    // exit(1);
    cout << "---------------------------------------------------"<<endl; 
  }
  cput = clock()-cput; 

  outcase.open(casedir+"time.txt"); 
  outcase << ((double)(cput)) / (double)CLOCKS_PER_SEC << endl; 
  outcase.close(); 

  delete(grid); 


  // LinSys ls(6); 
  // VecX<double> x(6); x.uncompress();
  // ls.A = {
  //   { 2, -1,  0,  0,  0,  0}, 
  //   {-1,  2, -1,  0,  0,  0}, 
  //   { 0, -1,  2, -1,  0,  0}, 
  //   { 0,  0, -1,  2, -1,  0}, 
  //   { 0,  0,  0, -1,  2, -1},
  //   { 0,  0,  0,  0, -1,  2}    
  // };
  // // ls.A[5] = 
  // // //vel.A[5][4] = -0.5; 

  // ls.A.info(); 
  // ls.b.info(); 

  // x = {1, -2, 3, -4, 9, 3}; 
  // cout << ls.A << endl; 
  // cout << ls.b << endl; 
  // cout << x<<endl; 
  // ls.b = (ls.A)*(x);
  // cout << ls.b << endl;
  // x = {0, 0, 0, 0, 0, 0}; 

  // ls.x = &x;
  // ls.setLimits(1e-8); 

  // // //vel.Gauss_Seidel();
  // ls.BiCGSTAB();
  // cout <<"x: "<< x<< endl;



  return 0; 
}
