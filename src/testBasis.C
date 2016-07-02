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
// testBasis.C
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#include "Basis.h"
#include "Context.h"
#include <iostream>
#include <new>
#include <cstdlib>
#include <cassert>
using namespace std;

#ifdef USE_MPI
#include <mpi.h>
#endif

int main(int argc, char **argv)
{
  // use: testBasis a0x a0y a0z a1x a1y a1z a2x a2y a2z ecut kx ky kz npr npc
#if USE_MPI
  MPI_Init(&argc,&argv);
#endif
  {
    if ( argc !=16 )
    {
      cout <<
      " use: testBasis a0x a0y a0z a1x a1y a1z a2x a2y a2z ecut kx ky kz npr npc"
      << endl;
      return 1;
    }
    const D3vector a0(atof(argv[1]),atof(argv[2]),atof(argv[3]));
    const D3vector a1(atof(argv[4]),atof(argv[5]),atof(argv[6]));
    const D3vector a2(atof(argv[7]),atof(argv[8]),atof(argv[9]));
    
    double ecut = atof(argv[10]);
    D3vector kpoint(atof(argv[11]),atof(argv[12]),atof(argv[13]));
    int npr = atoi(argv[14]);
    int npc = atoi(argv[15]);
    
    Context ctxt(npr,npc);
    UnitCell cell(a0,a1,a2);
    
    Basis basis(ctxt,kpoint);
    try
    {
      basis.resize(cell,cell,ecut);
    }
    catch ( bad_alloc )
    {
      cout << " bad_alloc caught in Basis::resize" << endl;
      throw;
    }
    
    //cout << basis;
    
    //Basis b2(basis);
    //cout << b2;
  }
#if USE_MPI
  MPI_Finalize();
#endif
}
