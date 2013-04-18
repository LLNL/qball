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
// testBasisMapping.C
//
////////////////////////////////////////////////////////////////////////////////

#include "Context.h"
#include "Basis.h"
#include "UnitCell.h"
#include "BasisMapping.h"
#include "Timer.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>
using namespace std;

#ifdef USE_MPI
#include <mpi.h>
#endif

int main(int argc, char **argv)
{
#if USE_MPI
  MPI_Init(&argc,&argv);
#endif
  {
    // use: testBasisMapping a0 a1 a2 b0 b1 b2 c0 c1 c2 ecut npr npc
    double err;
    assert(argc==13);
    D3vector a(atof(argv[1]),atof(argv[2]),atof(argv[3]));
    D3vector b(atof(argv[4]),atof(argv[5]),atof(argv[6]));
    D3vector c(atof(argv[7]),atof(argv[8]),atof(argv[9]));
    UnitCell cell(a,b,c);
    double ecut = atof(argv[10]);
    D3vector kpoint;

    int npr = atoi(argv[11]);
    int npc = atoi(argv[12]);

    Timer tm;

    Context ctxt(npr,npc);

    Basis basis(ctxt,kpoint);
    basis.resize(cell,cell,ecut);

    cout << " np0=" << basis.np(0)
         << " np1=" << basis.np(1)
         << " np2=" << basis.np(2) << endl;
    cout << " basis.size=" << basis.size() << endl;
    BasisMapping bmap(basis);
    cout << " zvec_size=" << bmap.zvec_size() << endl;
    cout << " np012loc=" << bmap.np012loc() << endl;

    vector<complex<double> > zvec(bmap.zvec_size());
    vector<complex<double> > ct(bmap.np012loc());

    vector<complex<double> > f(basis.localsize());

    for ( int i = 0; i < f.size(); i++ )
    {
      f[i] = exp(-basis.g2(i));
      cout << ctxt.mype() << ": "
           << i << " " << basis.idx(3*i) << " "
           << basis.idx(3*i+1) << " " << basis.idx(3*i+2) << " "
           << f[i] << endl;
    }
    bmap.vector_to_zvec(&f[0],&zvec[0]);
    bmap.transpose_fwd(&zvec[0],&ct[0]);

    for ( int k = 0; k < bmap.np2loc(); k++ )
      for ( int j = 0; j < bmap.np1(); j++ )
        for ( int i = 0; i < bmap.np0(); i++ )
        {
          int index = i + bmap.np0() * ( j + bmap.np1() * k );
          cout << ctxt.mype() << ": "
               << i << " " << j << " " << k << " " << ct[index] << endl;
        }

    // transpose back to zvec
    for ( int i = 0; i < zvec.size(); i++ )
      zvec[i] = 0.0;
    bmap.transpose_bwd(&ct[0],&zvec[0]);

    // transpose back to array f2
    vector<complex<double> > f2(basis.localsize());
    bmap.zvec_to_vector(&zvec[0],&f2[0]);

    double sum = 0.0;
    for ( int i = 0; i < f.size(); i++ )
    {
      sum += abs(f[i]-f2[i]);
    }
    cout << " total error: " << sum << endl;
  }
#if USE_MPI
  MPI_Finalize();
#endif
}
