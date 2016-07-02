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
// StrainCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef STRAINCMD_H
#define STRAINCMD_H

#include <string>
#include <cstdlib>
#include <iostream>
using namespace std;

#include "UserInterface.h"
#include "Sample.h"

class StrainCmd : public Cmd
{
  public:

  Sample *s;

  StrainCmd(Sample *sample) : s(sample) {};

  char *name(void) const { return "strain"; }
  char *help_msg(void) const
  {
    return
    "\n strain\n\n"
    " syntax: strain [-atomsonly] [-inverse] uxx uyy uzz uxy uyz uxz \n\n"
    "   The strain command deforms the unit cell according to\n"
    "   the given strain tensor u. Atomic positions are modified\n"
    "   correspondingly. If -atomsonly is specified, only atomic\n"
    "   positions are modified and the cell is unchanged. If the\n"
    "   -inverse flag is used, the inverse transformation is applied.\n\n";
  }

  int action(int argc, char **argv)
  {
    const string usage =
    " use: strain [-atomsonly] [-inverse] uxx uyy uzz uxy uyz uxz";
    // strain must have 7, 8 or 9 arguments including the command name
    if ( argc < 7 || argc > 9 )
    {
      if ( ui->oncoutpe() )
        cout << usage << endl;
      return 1;
    }

    bool atomsonly = false;
    bool inverse = false;
    int iu = 1;
    for ( int iarg = 0; iarg < 2; iarg++ )
    {
      string arg = argv[iarg];
      if ( arg == "-atomsonly" )
      {
        atomsonly = true;
        iu++;
      }
      if ( arg == "-inverse" )
      {
        inverse = true;
        iu++;
      }
    }
    if ( argc != iu+6 )
    {
      if ( ui->oncoutpe() )
        cout << usage << endl;
      return 1;
    }

    vector<double> u(6);
    for ( int i = 0; i < 6; i++ )
      u[i] = atof(argv[iu+i]);

#ifdef DEBUG
    if ( ui->oncoutpe() )
    {
      cout << " u: input: " << endl;
      cout << u[0] << " " << u[3] << " " << u[5] << endl;
      cout << u[3] << " " << u[1] << " " << u[4] << endl;
      cout << u[5] << " " << u[4] << " " << u[2] << endl;
    }
#endif

    if ( inverse )
    {
      // replace u by inv(I+u)
      const double a0 = u[0] + 1.0;
      const double a1 = u[1] + 1.0;
      const double a2 = u[2] + 1.0;
      const double a3 = u[3];
      const double a4 = u[4];
      const double a5 = u[5];

      const double det = a0*(a1*a2-a4*a4)-a3*(a3*a2-a5*a4)+a5*(a3*a4-a5*a1);
      if ( fabs(det) < 1.e-8 )
      {
        if ( ui->oncoutpe() )
        {
          cout << " transformation is near-singular: det(I+u) < 1.e-8" << endl;
        }
        return 1;
      }
      const double detinv = 1.0 / det;
      // replace u with (a^-1 - I)
      u[0] = detinv * ( a1*a2-a4*a4 ) - 1.0;
      u[1] = detinv * ( a0*a2-a5*a5 ) - 1.0;
      u[2] = detinv * ( a0*a1-a3*a3 ) - 1.0;
      u[3] = detinv * ( a4*a5-a3*a2 );
      u[4] = detinv * ( a5*a3-a0*a4 );
      u[5] = detinv * ( a3*a4-a5*a1 );

    }
#ifdef DEBUG
    if ( ui->oncoutpe() )
    {
      cout << " u: used: " << endl;
      cout << u[0] << " " << u[3] << " " << u[5] << endl;
      cout << u[3] << " " << u[1] << " " << u[4] << endl;
      cout << u[5] << " " << u[4] << " " << u[2] << endl;
    }
#endif

    UnitCell cell = s->atoms.cell();
    if ( !atomsonly )
    {
      // compute cell deformation: damat = u * amat
      vector<double> damat(9);

      cell.smatmult3x3(&u[0], cell.amat(), &damat[0]);

      D3vector a0new(cell.a(0)+D3vector(damat[0],damat[1],damat[2]));
      D3vector a1new(cell.a(1)+D3vector(damat[3],damat[4],damat[5]));
      D3vector a2new(cell.a(2)+D3vector(damat[6],damat[7],damat[8]));
      cell.set(a0new,a1new,a2new);

      s->atoms.set_cell(a0new,a1new,a2new);

      if (s->wf.hasdata()) {
        s->wf.resize(cell,s->wf.refcell(),s->wf.ecut());
        if ( s->wfv != 0 )
        {
          s->wfv->resize(cell,s->wf.refcell(),s->wf.ecut());
          s->wfv->clear();
        }  
      }
      else {
        s->wf.set_cell(cell);
        if ( s->wfv != 0 )
          s->wfv->set_cell(cell);
      }
      
      if ( ui->oncoutpe() )
      {
        cout << s->atoms.cell();
      }
    }

    // atomic displacements
    vector<vector<double> > tau;
    s->atoms.get_positions(tau);

    vector<vector<double> > dtau;
    dtau.resize(tau.size());
    for ( int is = 0; is < dtau.size(); is++ )
      dtau[is].resize(tau[is].size());

    for ( int is = 0; is < tau.size(); is++ )
      for ( int ia = 0; ia < s->atoms.na(is); ia++ )
        cell.vecsmult3x3(&u[0], &tau[is][3*ia], &dtau[is][3*ia]);
    for ( int is = 0; is < tau.size(); is++ )
      for ( int i = 0; i < tau[is].size(); i++ )
        tau[is][i] += dtau[is][i];

    s->atoms.set_positions(tau);

    return 0;
  }
};
#endif
