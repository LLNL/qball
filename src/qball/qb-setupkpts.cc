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
//
// qb-setupkpts reads in a Qbox coordinate (.sys) file, and calculates the symmetry operations
// from the cell vectors and atomic positions.  An irreducible Monkhorst-Pack k-point grid is
// then generated, with symmetry-equivalent k-points removed.
//
// written by Erik Draeger, 1/22/2007

#include <config.h>

#include <ui/UserInterface.h>
#include "Context.h"
#include "Sample.h"
#include "AtomSet.h"
#include "SymOpSet.h"
#include "SymOp.h"

#include <ui/AtomCmd.h>
#include <ui/SpeciesCmd.h>
#include <ui/SetCmd.h>
#include "Cell.h"
#include "RefCell.h"
#include <qball/UnitCell.h>
#include "Basis.h"

#include <cassert>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <valarray>
#include <map>
using namespace std;

const double threshold = 1.E-7;
const double twopiinv = 0.5/M_PI;

int main(int argc, char **argv) {

#if USE_MPI
  MPI_Init(&argc,&argv);
#endif

  if (argc <= 5 && argc != 1) {
    cout << "Usage:  qb-setupkpts [nx] [ny] [nz] [ecut] [.sys file(s)]" << endl;
    cout << "   nx,ny,nz = size of Monkhorst-Pack k-point grid" << endl;
    cout << "   ecut = plane wave cutoff in Ry" << endl;
    exit(1);
  }

  int nx,ny,nz;
  double ecut;
  //vector<char*> sysfiles;
  vector<string> sysfiles;
  if (argc > 1) {
    nx = atoi(argv[1]);
    ny = atoi(argv[2]);
    nz = atoi(argv[3]);
    ecut = atof(argv[4]);
    for (int f=5; f<argc; f++) {
      string tmp(argv[f]);
      sysfiles.push_back(tmp);
      //sysfiles.push_back(argv[f]);
    }
  }
  else { // interactive mode
    cout << "Usage:  qb-setupkpts [nx] [ny] [nz] [ecut] [.sys file(s)]" << endl;
    cout << "   nx,ny,nz = size of Monkhorst-Pack k-point grid" << endl;
    cout << "   ecut = plane wave cutoff in Ry" << endl;
    cout << endl;
    cout << "Interactive mode:" << endl;
    cout << "  Enter order of Monkhorst-Pack grid (nx ny nz): " << endl;
    cin >> nx >> ny >> nz;
    cout << "  Enter plane wave cutoff in Rydbergs:  " << endl;
    cin >> ecut;
    string insys;
    cout << "  Enter .sys filename:  " << endl;
    cin >> insys;
    //char* tmp = (char*)insys.c_str();
    //sysfiles.push_back(tmp);
    sysfiles.push_back(insys);
  }
  cout << "Calculating symmetry and k-point for " << nx << "x" << ny << "x" << nz << " Monkhorst-Pack grid, with ecut = " << ecut << endl;

  Context ctxt;
  //for (int f=0; f<sysfiles.size(); f++) {
  for (int f=0; f<1; f++) {
    //char* infilename = sysfiles[f];
    string infilename = sysfiles[f];
    ifstream infile;

    Sample* s = new Sample(ctxt);
    UserInterface ui(ctxt);
    ui.addCmd(new AtomCmd(s));
    ui.addCmd(new SpeciesCmd(s));
    ui.addCmd(new SetCmd(s));
    ui.addVar(new Cell(s));
    ui.addVar(new RefCell(s));

    // open sys file and read cell vectors and atom positions
    if (ctxt.oncoutpe())
      cout << "Opening file " << infilename << endl;
    infile.open(infilename.c_str(),ios::in);
    bool echo = false;
    ui.processCmds(infile, "[qbox]", echo);

    if (ctxt.oncoutpe())
      cout << s->atoms;

    // from unit cell, choose set of symmetry operators
    SymOpSet fullset, sysset;
    const UnitCell& cell = s->wf.cell();
    const UnitCell& refcell = s->wf.refcell();

    // for now, always use cubic operators -- later, we'll want to generalize this
    fullset.generateOps("cubic");
    if (ctxt.oncoutpe())
      cout << "Generated full set of cubic symmetry operations:  nsym = " << fullset.nsym() << endl;

    // convert symmetry operators to crystal coordinates
    fullset.convertOpsToXtal(cell);

    // convert atom positions to crystal coordinates
    vector< vector< D3vector> > r_xtal;
    r_xtal.resize(s->atoms.atom_list.size());
    for ( int is = 0; is < s->atoms.atom_list.size(); is++ ) {
      r_xtal[is].resize(s->atoms.atom_list[is].size());
      for ( int ia = 0; ia < s->atoms.atom_list[is].size(); ia++ )
        r_xtal[is][ia] = cell.cart_to_crystal(s->atoms.atom_list[is][ia]->position());
    }


    //cout << "before calling printXtal" << endl;
    //fullset.printXtal(cout);
    //cout << endl;
    //cout << "after calling printXtal" << endl;


    // use atom positions to select allowed symmetry operators
    int symlist[fullset.nsym()];
    for (int i=0; i<fullset.nsym(); i++) 
      symlist[i] = 0;
    for (int sy = 0; sy < fullset.nsym(); sy++) {
      SymOp* symop = fullset.returnOp(sy);

      bool allowed = true;
      for ( int is = 0; is < s->atoms.atom_list.size(); is++ ) {

        // for each atom, calculate the position resulting from each symmetry operation
        // and check that there is an atom of the same species at that position
        for ( int ia = 0; ia < s->atoms.atom_list[is].size(); ia++ ) {
          D3vector r_sym = symop->applyToVector(r_xtal[is][ia],true);

          bool atommatch = false;
          // compare this position against other atoms
          for ( int ja = 0; ja < s->atoms.atom_list[is].size(); ja++ ) {
            // in crystal coordinates, agreement can be off by an integer
            double xdiff = (r_sym.x - r_xtal[is][ja].x);
            int xdiffnint = (int) ( xdiff > 0.0 ? xdiff+0.5 : xdiff-0.5);
            double ydiff = (r_sym.y - r_xtal[is][ja].y);
            int ydiffnint = (int) ( ydiff > 0.0 ? ydiff+0.5 : ydiff-0.5);
            double zdiff = (r_sym.z - r_xtal[is][ja].z);
            int zdiffnint = (int) ( zdiff > 0.0 ? zdiff+0.5 : zdiff-0.5);

            if (fabs(xdiff - (double)xdiffnint) < threshold && fabs(ydiff - (double)ydiffnint) < threshold && fabs(zdiff - (double)zdiffnint) < threshold) {
              if (atommatch) {
                if (ctxt.oncoutpe())
                  cout << "ERROR:  multiple symmetry equivalent atoms found for atom " << ia << ", check sys file!" << endl;
                return 1;
              }
              else 
                atommatch = true;
            }

          }
          if (!atommatch)
            allowed = false;
        }
        if (allowed) {
          if (ctxt.oncoutpe())
            cout << "Symmetry operation " << sy << " accepted." << endl;
          sysset.addOp(symop);
          symlist[sy] = 1;
        }

      }
    }
  
    if (ctxt.oncoutpe())
      cout << "From atom positions, " << sysset.nsym() << " sym ops allowed." << endl;

    // check if cell is a supercell (which disallows fractional translations)
    bool supercell = false;
    for ( int is = 0; is < s->atoms.atom_list.size(); is++ ) {
      for ( int ia = 0; ia < s->atoms.atom_list[is].size(); ia++ ) {
        for ( int ja = 0; ja < s->atoms.atom_list[is].size(); ja++ ) {
          if (ia != ja) {
            
            // test every possible supercell displacement:  (r_ia - r_ja)
            D3vector ft = r_xtal[is][ia] - r_xtal[is][ja];
            bool thisft = true;
            
            for ( int js = 0; js < s->atoms.atom_list.size(); js++ ) {
              for ( int ka = 0; ka < s->atoms.atom_list[js].size(); ka++ ) {
                D3vector kdisp = r_xtal[js][ka] - ft;
                bool match = false;
                for ( int la = 0; la < s->atoms.atom_list[js].size(); la++ ) {
                  D3vector diff = r_xtal[js][la] - kdisp;
                  while (diff.x < -0.5) diff.x += 1.0;
                  while (diff.y < -0.5) diff.y += 1.0;
                  while (diff.z < -0.5) diff.z += 1.0;
                  while (diff.x > 0.5) diff.x -= 1.0;
                  while (diff.y > 0.5) diff.y -= 1.0;
                  while (diff.z > 0.5) diff.z -= 1.0;
                  if (length(diff) < 2.*threshold) {
                    match = true;
                    break;
                  }
                }
                if (!match) {
                  thisft = false;
                  break;
                }
              }
            }
            if (thisft) {
              if (ctxt.oncoutpe())
                cout << "This system is a supercell -- disabling fractional translations." << endl;
              supercell = true;
            }
          }
          if (supercell) break;
        }
        if (supercell) break;
      }
      if (supercell) break; 
    }


    // need to calculate charge density grid so we can check that 
    // fractional translations are compatible with it
    D3vector tmpkpoint(0.0,0.0,0.0);
    Basis* vbasis_ = new Basis(ctxt, tmpkpoint);
    vbasis_->resize(cell,refcell,4.0*ecut);
    int np0v = vbasis_->np(0);
    int np1v = vbasis_->np(1);
    int np2v = vbasis_->np(2);
    delete vbasis_;
    
    // look for fractional translations
    if (!supercell) {
      for (int sy = 0; sy < fullset.nsym(); sy++) {
        if (symlist[sy] == 0) {
          SymOp* symop = fullset.returnOp(sy);
          bool ftfound = false;
          for ( int is = 0; is < s->atoms.atom_list.size(); is++ ) {
            for ( int ia = 0; ia < s->atoms.atom_list[is].size(); ia++ ) {
              D3vector r_sym = symop->applyToVector(r_xtal[is][ia],false);
              for ( int ja = 0; ja < s->atoms.atom_list[is].size(); ja++ ) {
                
                // test every possible displacement:  (r_ia* - r_ja)
                D3vector ft = r_sym - r_xtal[is][ja];
                while (ft.x < 0.0) ft.x += 1.0;
                while (ft.y < 0.0) ft.y += 1.0;
                while (ft.z < 0.0) ft.z += 1.0;
                while (ft.x >= 1.0) ft.x -= 1.0;
                while (ft.y >= 1.0) ft.y -= 1.0;
                while (ft.z >= 1.0) ft.z -= 1.0;
                if (fabs(ft.x) < threshold) ft.x = 0.0;
                if (fabs(ft.y) < threshold) ft.y = 0.0;
                if (fabs(ft.z) < threshold) ft.z = 0.0;
                
                if (length(ft) > threshold) { 
                  bool thisft = true;
                  for ( int js = 0; js < s->atoms.atom_list.size(); js++ ) {
                    for ( int ka = 0; ka < s->atoms.atom_list[js].size(); ka++ ) {
                      D3vector rk_sym = symop->applyToVector(r_xtal[js][ka],false);
                      D3vector rk_star = rk_sym - ft;
                      bool match = false;
                      for ( int la = 0; la < s->atoms.atom_list[js].size(); la++ ) {
                        D3vector diff = r_xtal[js][la] - rk_star;
                        while (diff.x < -0.5) diff.x += 1.0;
                        while (diff.y < -0.5) diff.y += 1.0;
                        while (diff.z < -0.5) diff.z += 1.0;
                        while (diff.x > 0.5) diff.x -= 1.0;
                        while (diff.y > 0.5) diff.y -= 1.0;
                        while (diff.z > 0.5) diff.z -= 1.0;
                        if (length(diff) < 2.*threshold) {
                          match = true;
                          break;
                        }
                      }
                      if (!match) {
                        thisft = false;
                        break;
                      }
                    }
                  }
                  if (thisft) {
                    // check against density grid
                    bool gridok = true;
                    int ft0_ = (int)(ft.x*np0v);
                    int ft1_ = (int)(ft.y*np1v);
                    int ft2_ = (int)(ft.z*np2v);

                    if (ft0_ != 0 && (ft.x*np0v)/(double)ft0_ != 1.) gridok = false;
                    if (ft1_ != 0 && (ft.y*np1v)/(double)ft1_ != 1.) gridok = false;
                    if (ft2_ != 0 && (ft.z*np2v)/(double)ft2_ != 1.) gridok = false;

                    if (!gridok) { 
                      if (ctxt.oncoutpe())
                        cout << "Sym op " << sy << ", fractional translation found: " << ft << ", but incompatible with density grid at ecut = " << ecut << ".  Skipping." << endl;
                    }
                    else {
                      if (ctxt.oncoutpe())
                        cout << "Sym op " << sy << ", fractional translation found! " << ft << endl;
                      ftfound = true;
                      symlist[sy] = 1;
                      symop->setFractionalTranslation(ft);
                      sysset.addOp(symop);
                    }                    
                  }
                }
                if (ftfound) break;
              }
              if (ftfound) break;
            }
            if (ftfound) break;
          }
        }
      }
    }
    
    if (ctxt.oncoutpe())
      cout << "Including fractional translations, " << sysset.nsym() << " sym ops allowed." << endl;

    
    // generate Monkhorst-Pack grid of k-points
    vector< D3vector> fullgrid;
    fullgrid.resize(nx*ny*nz);
    int nkpts = 0;
    for (int kz=1; kz<=nz; kz++) {
      for (int ky=1; ky<=ny; ky++) {
        for (int kx=1; kx<=nx; kx++) {
          fullgrid[nkpts].x = (double)(2.*kx-nx-1.)/(2.*nx);
          fullgrid[nkpts].y = (double)(2.*ky-ny-1.)/(2.*ny);
          fullgrid[nkpts].z = (double)(2.*kz-nz-1.)/(2.*nz);
          nkpts++;
        }
      }
    }

    // reduce grid to irreducible set
    vector< D3vector> kptsym;
    vector< double> kptwt;
    
    kptsym.push_back(fullgrid[0]);
    kptwt.push_back(1.0);
    for (int k=1; k<fullgrid.size(); k++) {

      bool symequiv = false;
      for (int sy = 0; sy < sysset.nsym(); sy++) {
        SymOp* symop = sysset.returnOp(sy);
        D3vector fullsym = symop->applyToVector(fullgrid[k],false);
        while (fullsym.x > 0.5) fullsym.x -= 1.0;
        while (fullsym.y > 0.5) fullsym.y -= 1.0;
        while (fullsym.z > 0.5) fullsym.z -= 1.0;
        while (fullsym.x < -0.5) fullsym.x += 1.0;
        while (fullsym.y < -0.5) fullsym.y += 1.0;
        while (fullsym.z < -0.5) fullsym.z += 1.0;

        for (int j=0; j<kptsym.size(); j++) {
          D3vector kdiff = fullsym-kptsym[j];
          if (length(kdiff) < threshold) {
            symequiv = true;
            kptwt[j] += 1.0;
            break;
          }
        }
        if (symequiv) break;
      }
      if (!symequiv) {
        kptsym.push_back(fullgrid[k]);
        kptwt.push_back(1.0);
      }
    }

    // print out results
    if (ctxt.oncoutpe()) {
      for (int sy = 0; sy < sysset.nsym(); sy++) {
        SymOp* symop = sysset.returnOp(sy);
        symop->print(cout);
      }
      cout << endl;
      for (int j=0; j<kptsym.size(); j++) 
        cout << "kpoint " << kptsym[j] << " " << kptwt[j] << endl;
    }

    fullset.clear(); // deletes SymOp objects used by both fullset and sysset
    //delete s;
  }
#if USE_MPI
  MPI_Finalize();
#endif
  return 0;

}
