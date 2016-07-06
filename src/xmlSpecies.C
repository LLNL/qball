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
// xmlSpecies.C: transform a GP pseudopotential file into a fpmd xml species
//
// use: xmlSpecies < psfile > xmlfile
//

#include <config.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include "qbox_xmlns.h"
using namespace std;

int main(int argc, char **argv)
{
  int np, zion,lloc,nquad,lmax1,lmax;
  double rcut,pmass,radius,covradius;

  char buf[256];
  char psname_buf[256];
  char color_buf[256];
  int ntokens;
  // skip # comment lines
  while ( ( (char) cin.peek() ) == '#' )
  {
    while ( cin.get() != '\n' );
  }

  cin.getline(&buf[0],256);

  ntokens =
    sscanf(&buf[0]," %d %d %lf %lf %d %s %s %lf %lf %d %d",
    &np, &zion, &rcut, &pmass, &lmax1, &psname_buf[0],
    &color_buf[0], &radius, &covradius, &lloc, &nquad );

  lmax = lmax1 - 1;

  if ( ntokens == 9 )
  {
    // lloc and nquad were not read: use default values
    // default value of lloc is lmax
    lloc = lmax;
    // default value of nquad is 0 (Kleinman-Bylander)
    nquad = 0;
  }

  cout << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
  cout << "<fpmd:species xmlns:fpmd=\""
       << qbox_xmlns()
       << "\"" << endl;
  cout << "  xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" " << endl;
  cout << "  xsi:schemaLocation=\"";
  cout << qbox_xmlns();
  cout << "  species.xsd\">" << endl;
  cout << "<description>" << psname_buf << "</description>" << endl;
  cout << "<symbol>X</symbol>" << endl;
  cout << "<atomic_number> A </atomic_number>" << endl;
  cout << "<mass>" << pmass << "</mass>" << endl;
  cout << "<norm_conserving_pseudopotential>" << endl;
  cout << "<valence_charge>" << zion << "</valence_charge>" << endl;
  cout << "<lmax>" << lmax << "</lmax>" << endl;
  cout << "<llocal>" << lloc << "</llocal>" << endl;
  cout << "<nquad>" << nquad << "</nquad>" << endl;
  cout << "<rquad> RQUAD </rquad>" << endl;
  cout << "<mesh_spacing> DR </mesh_spacing>" << endl;
  cout.setf(ios::scientific,ios::floatfield);
  for ( int l = 0; l < lmax1; l++ )
  {
    cout << "<projector l=\"" << l << "\" size=\"" << np << "\">" << endl;
    cout << "<radial_potential>" << endl;
    int npdum;
    double rdum;
    if ( l > 0 )
    {
      // skip # comment lines
      while ( ( (char) cin.peek() ) == '#' )
      {
        while ( cin.get() != '\n' );
      }
      cin >> npdum;
      while ( cin.get() != '\n' );
    }

    double r,v,phi;
    for ( int i = 0; i < np; i++ )
    {
      cin >> r >> v;
      cout << setw(13) << setprecision(6) << v << endl;
      while ( cin.get() != '\n' );
    }
    cout << "</radial_potential>" << endl;

    // read phi except if l == lmax && lloc == lmax_
    if ( !( l == lmax && lloc == lmax ) )
    {
      cout << "<radial_function>" << endl;
      // skip # comment lines
      while ( ( (char) cin.peek() ) == '#' )
      {
        while ( cin.get() != '\n' );
      }
      cin >> npdum;
      while ( cin.get() != '\n' );
      for ( int i = 0; i < np; i++ )
      {
        cin >> r >> phi; while ( cin.get() != '\n' );
        cout << setw(13) << setprecision(6) << phi << endl;
      }
      cout << "</radial_function>" << endl;
    }
    cout << "</projector>" << endl;
  }

  cout << "</norm_conserving_pseudopotential>" << endl;
  cout << "</fpmd:species>" << endl;

  return 0;
}
