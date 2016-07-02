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
// testSample.C
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#include <iostream>
using namespace std;

#include "Context.h"
#include "SlaterDet.h"
#include "Sample.h"
#include "D3vector.h"

int main(int argc, char** argv)
{
#if USE_MPI
  MPI_Init(&argc,&argv);
#endif
  // extra scope to ensure that BlacsContext objects get destructed before
  // the MPI_Finalize call
  {
    Context ctxt;
 
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int namelen;
    PMPI_Get_processor_name(processor_name,&namelen);
    cout << " Process " << ctxt.mype() << " on " << processor_name << endl;
 
    Sample s(ctxt);
    
    D3vector cell(18,18,18);
    double ecut = 25.0;
    s.wf.resize(cell,cell,ecut);
    s.wf.set_nel(12*54);
    
    s.wf.randomize(1.e-4);
    s.wf.gram();
    cout << " ortho_error: " << s.wf.sd[0][0]->ortho_error() << endl;
  }
#if USE_MPI
  MPI_Finalize();
#endif
  return 0;
}
