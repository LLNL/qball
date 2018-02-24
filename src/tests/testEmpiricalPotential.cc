////////////////////////////////////////////////////////////////////////////////  
// Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
// qb@ll:  Qbox at Lawrence Livermore
//
// This file is part of qb@ll.
//
// Produced at the Lawrence Livermore National Laboratory. 
// Written by Erik Draeger (draeger1@llnl.gov) and Francois Gygi (fgygi@ucdavis.edu).
// Based on the Qbox code by Francois Gygi Copyright (c) 2008 
// LLNL-CODE-635376. All rights reserved. 
//
// qb@ll is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details, in the file COPYING in the
// root directory of this distribution or <http://www.gnu.org/licenses/>.
//

#include <config.h>
#include <fstream>
#include <iostream>
#include <cstring>
#include <string>
using namespace std;
#include <qball/Context.h>
#include "D3vector.h"
#include "EmpiricalPotential.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

int main(int argc, char **argv) {
#if USE_MPI
  MPI_Init(&argc,&argv);
#endif
  {
  
  Context ctxt;

  string pottype;
  if (ctxt.oncoutpe()) {
    cout << "Enter potential type (lennard-jones, morse or file): ";
    cin >> pottype;
  }
  if (pottype == "lennard-jones") { 
    cout << "Lennard-Jones potential selected, input epsilon and sigma:  ";
    double eps,sigma;
    cin >> eps >> sigma;

    EmpiricalPotential* ep_ = new EmpiricalPotential(ctxt,pottype,eps,sigma);

    string query = "n";
    while (! (query  == "y")) {
      if (ctxt.oncoutpe()) {
        double x1,y1,z1,x2,y2,z2;
        cout << "Enter position of atom 1:  ";
        cin >> x1 >> y1 >> z1;
        cout << "Enter position of atom 2:  ";
        cin >> x2 >> y2 >> z2;

        D3vector r1(x1,y1,z1);
        D3vector r2(x2,y2,z2);

        double rval = norm(r1-r2);
        
        D3vector f1 = ep_->force(r1,r2);

        cout << " r = " << rval << ", potential = " << ep_->pot(rval) << ", force on r1 = " << f1 << endl;
        
        cout << "Exit?  (y/n):  ";
        cin >> query;
      }
    }
    delete ep_;

  }
  else if (pottype == "morse") { 
    cout << "Morse potential selected, input D_e,alpha, and r_e:  ";
    double de,alpha,r_e;
    cin >> de >> alpha >> r_e;

    EmpiricalPotential* ep_ = new EmpiricalPotential(ctxt,pottype,de,alpha,r_e);

    string query = "n";
    while (! (query  == "y")) {
      if (ctxt.oncoutpe()) {
        double x1,y1,z1,x2,y2,z2;
        cout << "Enter position of atom 1:  ";
        cin >> x1 >> y1 >> z1;
        cout << "Enter position of atom 2:  ";
        cin >> x2 >> y2 >> z2;

        D3vector r1(x1,y1,z1);
        D3vector r2(x2,y2,z2);

        double rval = norm(r1-r2);
        
        D3vector f1 = ep_->force(r1,r2);

        cout << " r = " << rval << ", potential = " << ep_->pot(rval) << ", force on r1 = " << f1 << endl;
        
        cout << "Exit?  (y/n):  ";
        cin >> query;
      }
    }
    delete ep_;

  }
  else if (pottype == "file") { 

    string filename;
    //string filename("test.dat");

    if (ctxt.oncoutpe()) {
      cout << "Enter potential filename: ";
      cin >> filename;
    }
  
    EmpiricalPotential* ep_ = new EmpiricalPotential(ctxt,filename);


    if (ep_->filename() == "") {
      cout << "File not found, exiting..." << endl;
    }
    else {
      
      ofstream os;
      string outfile = filename + ".spline";
      os.open(outfile.c_str(),ofstream::out);    // text output
      
      //os.setf(ios::fixed,ios::floatfield);
      //os << setprecision(8);
      
      int np = ep_->npts();
      double rmax = ep_->r(np-1);
      int npsp = np*5 + 11;
      double deltar = rmax/npsp;
      for (int i=1; i<npsp; i++) {
        double rval = 0.05 + deltar*i;
        double potval = 0.0;
        if (rval < rmax)
          double potval = ep_->pot(rval);
        os << rval << "  " << potval << endl;
      }
      os.close();


      outfile = "epdata.dat";
      os.open(outfile.c_str(),ofstream::out);    // text output
      for (int i=0; i<np; i++) {
        os << ep_->r(i) << "  " << ep_->pot_grid(i) << endl;
      }
      os.close();
    }

    string query = "n";
    while (! (query  == "y")) {
      if (ctxt.oncoutpe()) {
        double x1,y1,z1,x2,y2,z2;
        cout << "Enter position of atom 1:  ";
        cin >> x1 >> y1 >> z1;
        cout << "Enter position of atom 2:  ";
        cin >> x2 >> y2 >> z2;

        D3vector r1(x1,y1,z1);
        D3vector r2(x2,y2,z2);

        double rval = norm(r1-r2);
        
        D3vector f1 = ep_->force(r1,r2);

        cout << " r = " << rval << ", potential = " << ep_->pot(rval) << ", force on r1 = " << f1 << endl;
        
        cout << "Exit?  (y/n):  ";
        cin >> query;
      }
    }
    delete ep_;
  }
  
  }
#if USE_MPI
  MPI_Finalize();
#endif
  return 0;
}  
  
  
  
