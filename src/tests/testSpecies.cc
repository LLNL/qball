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
// testSpecies.C
//
// use: testSpecies uri
//

#include <config.h>

#include "Species.h"
#include <qball/Context.h>
#include "SpeciesReader.h"
#include <iostream>
#include <cassert>
#include <string>
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
  
  Context ctxt;
  if ( argc != 2 )
  {
    cerr << "use: testSpecies uri" << endl;
    return 1;
  }
  
  Species s(ctxt,"unknown_name");
  
  SpeciesReader rdr(ctxt);
  
  string uri(argv[1]);
  
  try
  {
    cout << " s.uri() = " << s.uri() << endl;
    cout << " testSpecies: invoking readSpecies: uri=" << uri << endl;
    rdr.readSpecies(s,uri);
  }
  catch ( const SpeciesReaderException& e )
  {
    cout << " SpeciesReaderException caught in testSpecies" << endl;
    throw;
  }
  cerr << " SpeciesReader::readSpecies done" << endl;
  
  const double rcps = 1.0;
  
  try
  {
    s.initialize(rcps);
  }
  catch ( SpeciesInitException& e )
  {
    cout << " Exception in Species initialization: " << e.msg << endl;
    throw;
  }
  cerr << s;
  
#if 1
  
  if ( ctxt.oncoutpe() )
  {
  
  double dr = 0.01;
  double dg = 0.02;
  double v,dv;

  int n = 1000;
  
  for ( int l = 0; l <= s.lmax(); l++ )
  {
    cout << n << " Vps(l=" << l << ",r) " << endl;
    for ( int i = 0; i < n; i++ )
    {
      double r = i * dr;
      s.vpsr(l,r,v);
      cout << r << " " << v << endl;
    }

    cout << n << " dVps(l=" << l << ",r)/dr " << endl;
    for ( int i = 0; i < n; i++ )
    {
      double r = i * dr;
      s.dvpsr(l,r,v,dv);
      cout << r << " " << dv << endl;
    }
  }
  
  cout << n << " Vloc(g) " << endl;
  for ( int i = 0; i < n; i++ )
  {
    double g = i * dg;
    s.vlocg(g,v);
    cout << g << " " << v << endl;
  }
  
  cout << n << " dVloc(g)/dg " << endl;
  for ( int i = 0; i < n; i++ )
  {
    double g = i * dg;
    s.dvlocg(g,v,dv);
    cout << g << " " << dv << endl;
  }
  
  for ( int l = 0; l <= s.lmax(); l++ )
  {
    if ( l != s.llocal() )
    {
      cout << n << " Vnl(l=" << l << ",g) " << endl;
      for ( int i = 0; i < n; i++ )
      {
        double g = i * dg;
        s.vnlg(l,g,v);
        cout << g << " " << v << endl;
      }
   
      cout << n << " dVnl(l=" << l << ",g)/dg " << endl;
      for ( int i = 0; i < n; i++ )
      {
        double g = i * dg;
        s.dvnlg(l,g,v,dv);
        cout << g << " " << dv << endl;
      }
    }
  }
  
  } // oncoutpe
#endif

  }
#if USE_MPI
  MPI_Finalize();
#endif
  return 0;
}
