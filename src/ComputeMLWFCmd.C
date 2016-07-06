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
// ComputeMLWFCmd.C:
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#include "ComputeMLWFCmd.h"
#include<iostream>
#include "Context.h"
#include "SlaterDet.h"
using namespace std;

int ComputeMLWFCmd::action(int argc, char **argv)
{
  Wavefunction& wf = s->wf;
  SlaterDet& sd = *(wf.sd(0,0));

  if (wf.nkp() > 1 || wf.kpoint(0) != D3vector(0,0,0)) {
    if ( ui->oncoutpe() )
      cout << "<ERROR> compute_mlwf command only works for gamma-point calculations! </ERROR>" << endl;
    return 1;
  }


  mlwft = new MLWFTransform(sd);

  mlwft->compute_transform();
  mlwft->apply_transform(sd);

  if ( ui->oncoutpe() )
  {
    cout << " <mlwf_set size=\"" << sd.nst() << "\">" << endl;
    for ( int i = 0; i < sd.nst(); i++ )
    {
      D3vector ctr = mlwft->center(i);
      double sp = mlwft->spread(i);
      cout.setf(ios::fixed, ios::floatfield);
      cout.setf(ios::right, ios::adjustfield);
      cout << "   <mlwf center=\"" << setprecision(6)
           << setw(12) << ctr.x
           << setw(12) << ctr.y
           << setw(12) << ctr.z
           << " \" spread=\" " << sp << " \"/>"
           << endl;
    }
    cout << " </mlwf_set>" << endl;
    D3vector edipole = mlwft->dipole();
    cout << " <electronic_dipole> " << edipole
         << " </electronic_dipole>" << endl;
    D3vector idipole = s->atoms.dipole();
    cout << " <ionic_dipole> " << idipole
         << " </ionic_dipole>" << endl;
    cout << " <total_dipole> " << idipole + edipole
         << " </total_dipole>" << endl;
    cout << " <total_dipole_length> " << length(idipole + edipole)
         << " </total_dipole_length>" << endl;
  }
  delete mlwft;
  return 0;
}
