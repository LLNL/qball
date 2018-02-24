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
////////////////////////////////////////////////////////////////////////////////
//
// ParOptCmd.C:
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#include <fstream>
#include<iostream>
#include<ctime>
#include<cassert>
#include "ParOptCmd.h"
#include <qball/ParallelOptimizer.h>
#include <qball/Timer.h>
using namespace std;


int ParOptCmd::action(int argc, char **argv) {

  if ( argc < 3 || argc > 5) {
    if ( ui->oncoutpe() )
      cout << " use: paropt filename niter_ionic" << endl;
      cout << " use: paropt filename niter_ionic niter_scf" << endl;
      cout << " use: paropt filename niter_ionic niter_scf niter_nonscf" << endl;
    return 1;
  }
  
  if ( s->wf.nst() == 0 ) {
    if ( ui->oncoutpe() )
      cout << " <!-- ParOptCmd: no states, cannot run -->" << endl;
    return 1;
  }
  if ( s->wf.ecut() == 0.0 ) {
    if ( ui->oncoutpe() )
      cout << " <!-- ParOptCmd: ecut = 0.0, cannot run -->" << endl;
    return 1;
  }
  
  char* parcmdfile = argv[1];
  
  int niter = atoi(argv[2]);
  int nite = 1;
  int nitscf = 1;
  if ( argc == 4 ) {
    nitscf = atoi(argv[3]);
  }
  else if ( argc == 5 ) {
    nitscf = atoi(argv[3]);
    nite = atoi(argv[4]);
  }

  ParallelOptimizer* popt = new ParallelOptimizer(*s);
  
  if (ui->oncoutpe() )
    cout << "<parallel_optimize niter_ionic=\"" << niter << "\" niter_scf=\"" << nitscf << "\" niter_nonscf=\"" << nite << "\">" << endl;


  popt->optimize(niter,nitscf,nite);

  if (ui->oncoutpe() )
    cout << "</parallel_optimize>" << endl;

  // write out input commands to parcmdfile
  if (ui->oncoutpe() ) {
    ofstream os;
    os.open(parcmdfile,ofstream::out);
    os << "set nrowmax " << s->wf.nrowmax() << endl;
    if (s->wf.nkp() > 1) 
      os << "set nparallelkpts " << s->wf.nparallelkpts() << endl;
    os.close();
  }

  delete popt;
  return 0;
}
