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
// PlotCmd.C:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: PlotCmd.C,v 1.4 2009-08-26 15:03:50 fgygi Exp $

#include <config.h>

#include "PlotCmd.h"
#include "isodate.h"
#include "release.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>

#include "Context.h"
#include "Sample.h"
#include "SampleReader.h"
#include "Basis.h"
#include "FourierTransform.h"
#include "SlaterDet.h"
#include <math/Matrix.h>
#include "Species.h"
#include "Atom.h"
#include "ChargeDensity.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
int PlotCmd::action(int argc, char **argv)
{
  string usage("  Use: plot [-density|-wf <nmin> [<nmax>]] filename");

  // parse arguments
  // plot filename               : plot atoms in xyz format
  // plot -density filename      : plot atoms and density in cube format
  // plot -wf <n> filename       : plot atoms and wf <n> in cube format
  // plot -wf <n1> <n2> filename : plot atoms and wfs <n1> to <n2> in cube fmt
  if ( (argc < 2) || (argc > 5) )
  {
    if ( ui->oncoutpe() )
      cout << usage << endl;
    return 1;
  }

  bool plot_atoms = false;
  bool plot_density = false;
  bool plot_wf = false;
  bool xyz = false;
  int nmin,nmax,nwf;
  string filename;

  if ( argc == 2 )
  {
    // plot filename : plot atoms in xyz format
    plot_atoms = true;
    xyz = true;
    filename = argv[1];
  }
  else if ( argc == 3 )
  {
    // plot -density filename  : plot atoms and density in cube format
    if ( strcmp(argv[1],"-density") )
    {
      if ( ui->oncoutpe() )
        cout << usage << endl;
      return 1;
    }
    filename = argv[2];
    plot_atoms = true;
    plot_density = true;
  }
  else if ( argc == 4 )
  {
    // plot -wf <n> filename : plot wavefunction <n>
    if ( strcmp(argv[1],"-wf") )
    {
      if ( ui->oncoutpe() )
        cout << usage << endl;
      return 1;
    }
    filename = argv[3];
    plot_atoms = true;
    plot_wf = true;
    nmin = atoi(argv[2]) - 1;
    nmax = nmin;
    nwf = 1;

    if ( nmin < 0 || nmax >= s->wf.nst() || nmin > nmax )
    {
      if ( ui->oncoutpe() )
        cout << " nmin or nmax incompatible with nst="
             << s->wf.nst() << endl;
      return 1;
    }
  }
  else if ( argc == 5 )
  {
    // plot -wf <nmin> <nmin> filename :
    // plot density of wfs <nmin> to <nmax>
    if ( strcmp(argv[1],"-wf") )
    {
      if ( ui->oncoutpe() )
        cout << usage << endl;
      return 1;
    }
    filename = argv[4];
    plot_atoms = true;
    plot_wf = true;
    nmin = atoi(argv[2]) - 1;
    nmax = atoi(argv[3]) - 1;
    nwf = nmax-nmin+1;

    if ( nmin < 0 || nmax >= s->wf.nst() || nmin > nmax )
    {
      if ( ui->oncoutpe() )
        cout << " nmin or nmax incompatible with nst="
             << s->wf.nst() << endl;
      return 1;
    }
  }

  ofstream os;
  int np0=1, np1=1, np2=1;
  vector<double> tmpr;

  const Context& ctxt = *s->wf.spincontext(0);
  if ( plot_density )
  {
    ChargeDensity cd(*s);
    cd.update_density();
    cd.update_rhor();
    np0 = cd.vft()->np0();
    np1 = cd.vft()->np1();
    np2 = cd.vft()->np2();
    const int np012 = cd.vft()->np012();
    tmpr.resize(np012);
    for ( int i = 0; i < cd.vft()->np012loc(); i++ )
      tmpr[i] = cd.rhor[0][i];

    // send blocks of tmpr to pe0
    // send from first context column only
    for ( int i = 0; i < ctxt.nprow(); i++ )
    {
      bool iamsending = ctxt.mycol() == 0 && i == ctxt.myrow();

      // send size of tmpr block
      int size=-1;
      if ( ctxt.oncoutpe() )
      {
        if ( iamsending )
        {
          // sending to self, size not needed
        }
        else
          ctxt.irecv(1,1,&size,1,i,0);
      }
      else
      {
        if ( iamsending )
        {
          size = cd.vft()->np012loc();
          ctxt.isend(1,1,&size,1,0,0);
        }
      }

      // send tmpr block
      if ( ctxt.oncoutpe() )
      {
        if ( iamsending )
        {
          // do nothing, data is already in place
        }
        else
        {
          int istart = cd.vft()->np0() * cd.vft()->np1() *
                       cd.vft()->np2_first(i);
          ctxt.drecv(size,1,&tmpr[istart],1,i,0);
        }
      }
      else
      {
        if ( iamsending )
        {
          ctxt.dsend(size,1,&tmpr[0],1,0,0);
        }
      }
    }
  } // plot_density
  else if ( plot_wf )
  {
    // compute wf and store in tmpr
    if ( ctxt.oncoutpe() )
    {
      ctxt.ibcast_send(1,1,&nwf,1);
      ctxt.ibcast_send(1,1,&nmin,1);
      ctxt.ibcast_send(1,1,&nmax,1);
    }
    else
    {
      ctxt.ibcast_recv(1,1,&nwf,1,0,0);
      ctxt.ibcast_recv(1,1,&nmin,1,0,0);
      ctxt.ibcast_recv(1,1,&nmax,1,0,0);
    }

    if ( nwf > 0 && s->wf.nst() == 0 )
    {
      cout << " no states in sample" << endl;
      return 1;
    }

    SlaterDet *sdp = s->wf.sd(0,0);
    const Basis& basis = sdp->basis();
    np0 = basis.np(0);
    np1 = basis.np(1);
    np2 = basis.np(2);
    FourierTransform ft(basis,np0,np1,np2);
    const ComplexMatrix& c = sdp->c();

    vector<complex<double> > wftmp(ft.np012loc());
    vector<double> wftmpr(ft.np012());
    tmpr.resize(ft.np012());
    for ( int n = nmin; n <= nmax; n++ )
    {
      assert(n < s->wf.nst());

      // compute real-space wavefunction

      // transform wf on ctxt.mycol() hosting state n
      if ( c.pc(n) == c.context().mycol() )
      {
        //os << " state " << n << " is stored on column "
        //     << ctxt_.mycol() << " local index: " << c_.y(n) << endl;
        int nloc = c.y(n); // local index
        ft.backward(c.cvalptr(c.mloc()*nloc),&wftmp[0]);

        double *a = (double*) &wftmp[0];
        if ( basis.real() )
        {
          // real function: plot wf
          for ( int i = 0; i < ft.np012loc(); i++ )
            wftmpr[i] = a[2*i];
        }
        else
        {
          // complex function: plot modulus
          for ( int i = 0; i < ft.np012loc(); i++ ) {
            wftmpr[i] = sqrt(a[2*i]*a[2*i] + a[2*i+1]*a[2*i+1]);
            cout.precision(15);
            cout << "AS: WF " << a[2*i] << "    " << a[2*i+1] << endl;
          }
        }
      }

      // send blocks of wftmpr to pe0
      for ( int i = 0; i < c.context().nprow(); i++ )
      {
        bool iamsending = c.pc(n) == c.context().mycol() &&
                          i == c.context().myrow();

        // send size of wftmpr block
        int size=-1;
        if ( c.context().oncoutpe() )
        {
          if ( iamsending )
          {
            // sending to self, size not needed
          }
          else
            c.context().irecv(1,1,&size,1,i,c.pc(n));
        }
        else
        {
          if ( iamsending )
          {
            size = ft.np012loc();
            c.context().isend(1,1,&size,1,0,0);
          }
        }

        // send wftmpr block
        if ( c.context().oncoutpe() )
        {
          if ( iamsending )
          {
            // do nothing, data is already in place
          }
          else
          {
            int istart = ft.np0() * ft.np1() * ft.np2_first(i);
            c.context().drecv(size,1,&wftmpr[istart],1,i,c.pc(n));
          }
        }
        else
        {
          if ( iamsending )
          {
            c.context().dsend(size,1,&wftmpr[0],1,0,0);
          }
        }
      }

      // process the data on task 0
      if ( c.context().oncoutpe() )
      {
        // wftmpr is now complete on task 0
        if ( nwf == 1 )
        {
          // only one wf
          for ( int i = 0; i < ft.np012(); i++ )
          {
            tmpr[i] = wftmpr[i];
          }
        }
        else
        {
          // multiple wfs, accumulate square
          for ( int i = 0; i < ft.np012(); i++ )
          {
            tmpr[i] += wftmpr[i]*wftmpr[i];
          }
        }
      }
    } // for n
  } // if plot_wf

  // tmpr now contains the function to plot on task 0

  if ( ctxt.oncoutpe() )
    os.open(filename.c_str());

  if ( plot_atoms )
  {
    if ( ctxt.oncoutpe() )
    {
      if ( xyz )
      {
        const double a0 = 0.529177;
        int natoms = s->atoms.size();
        os << natoms << endl;
        os << "Created " << isodate() << " by " << release() << endl;
        const int nsp = s->atoms.nsp();
        for ( int is = 0; is < nsp; is++ )
        {
          Species* sp = s->atoms.species_list[is];
          string symbol = sp->symbol();
          const int na = s->atoms.na(is);
          for ( int ia = 0; ia < na; ia++ )
          {
            Atom *ap = s->atoms.atom_list[is][ia];
            os << setprecision(5);
            os << symbol << " " << a0*ap->position() << endl;
          }
        }
      }
      else
      {
        // write header and atoms
        os << "Created " << isodate() << " by " << release() << endl;
        os << endl;

        int natoms = s->atoms.size();
        D3vector a0 = s->atoms.cell().a(0);
        D3vector a1 = s->atoms.cell().a(1);
        D3vector a2 = s->atoms.cell().a(2);
        os << natoms << " " << -0.5*(a0+a1+a2) << endl;

        // write unit cell
        os << np0 << " " << a0/np0 << endl;
        os << np1 << " " << a1/np1 << endl;
        os << np2 << " " << a2/np2 << endl;
        const int nsp = s->atoms.nsp();
        for ( int is = 0; is < nsp; is++ )
        {
          Species* sp = s->atoms.species_list[is];
          const int z = sp->atomic_number();
          const int na = s->atoms.na(is);
          for ( int ia = 0; ia < na; ia++ )
          {
            Atom *ap = s->atoms.atom_list[is][ia];
            os << setprecision(5);
            os << z << " " << ((double) z) << " " << ap->position() << endl;
          }
        }
      }
    }
  } // if plot_atoms

  if ( plot_density || plot_wf )
  {
    // process the function in tmpr
    if ( ctxt.oncoutpe() )
    {
      os.setf(ios::scientific,ios::floatfield);
      os << setprecision(5);
      for ( int i = 0; i < np0; i++ )
      {
        const int ip = (i + np0/2 ) % np0;
        for ( int j = 0; j < np1; j++ )
        {
          const int jp = (j + np1/2 ) % np1;
          for ( int k = 0; k < np2; k++ )
          {
            const int kp = (k + np2/2 ) % np2;
            os << setw(13) << tmpr[ip+np0*(jp+np1*kp)];
            if ( ( k % 6 ) == 5 )
              os << endl;
          }
          if ( ( np2 % 6 ) != 0 )
            os << endl;
        }
      }
    }
  } // if plot_density || plot_wf

  os.close();

  return 0;
}
