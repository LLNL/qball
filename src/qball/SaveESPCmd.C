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
// SaveESPCmd.C:
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#include "SaveESPCmd.h"
#include "EnergyFunctional.h"
#include "Wavefunction.h"
#include "Basis.h"
#include "SlaterDet.h"
#include "ChargeDensity.h"
#include "FourierTransform.h"
#include "Context.h"
#include "AtomSet.h"
#include "Species.h"
#include<fstream>
#include<iostream>
using namespace std;

int SaveESPCmd::action(int argc, char **argv)
{

  if ( !(argc == 2 || argc == 3) ) {
    if ( ui->oncoutpe() )
      cout << "  <!-- use: saveesp filename -->" << endl;
    return 1;
  }

  char* filename = argv[1];
  string format = "vmd";
  for ( int i = 1; i < argc; i++ ) {
    string arg(argv[i]);
    if ( arg=="-molmol" ) 
      format = "molmol";
    else if ( arg=="-gopenmol" ) 
      format = "gopenmol";
    else if ( arg=="-vmd" ) 
      format = "vmd";
    else if ( arg[0] != '-' && i == argc-1 )
      filename = argv[i];
    else {
      if ( ui->oncoutpe() )
        cout << "  <!-- use: saveden [-vmd|-molmol|-gopenmol] filename -->" 
             << endl;
      return 1;
    }
  }

  int mype;
#if USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&mype);
#else
  mype = 0;
#endif

  if (s->wf.nspin() == 2 && mype == 0 )
    cout << "<WARNING> saveesp command only prints spin = 0 potential </WARNING>" << endl;
  
  ofstream os;
  os.open(filename,ofstream::out);    // text output
  os.setf(ios::fixed,ios::floatfield);
  os << setprecision(8);

  StructureFactor sf;
  ChargeDensity cd_(*s);
  // if ultrasoft, calculate position-dependent functions
  if (s->ctrl.ultrasoft)
     cd_.update_usfns();
  cd_.update_density();
  EnergyFunctional ef_(*s,s->wf,cd_);
  const Context* ctxt_ = s->wf.spincontext(0);
  FourierTransform* ft_ = cd_.vft();
  Basis* vbasis_ = cd_.vbasis();
  const UnitCell& cell = s->wf.cell();
  const AtomSet& atoms = s->atoms;

  // compute electrostatic potential directly
  vector<vector<double> > tau0, rhops;
  vector<complex<double> > rhopst, rhog, v_g, v_r;
  const int nsp_  = atoms.nsp();
  const int ngloc = vbasis_->localsize();
  const double omega = cell.volume();
  const double omega_inv = 1.0 / omega;
  const double *const g2i = vbasis_->g2i_ptr();
  const double fpi = 4.0 * M_PI;

  tau0.resize(nsp_);
  for ( int is = 0; is < nsp_; is++ ) {
    const int na = atoms.na(is);
    tau0[is].resize(3*na);
  }
  atoms.get_positions(tau0,true);
  sf.init(tau0,*vbasis_);
  sf.update(tau0,*vbasis_);
  
  rhops.resize(nsp_);
  for ( int is = 0; is < nsp_; is++ )
    rhops[is].resize(ngloc);

  int natoms_total = 0;
  for ( int is = 0; is < nsp_; is++ ) {
     natoms_total += atoms.na(is);
     Species *s = atoms.species_list[is];
     const double * const g = vbasis_->g_ptr();
     complex<double> *sfg = &sf.sfac[is][0];
     for ( int ig = 0; ig < ngloc; ig++ ) {
        const complex<double> sg = sfg[ig];
        rhops[is][ig] = s->rhopsg(g[ig]) * omega_inv;
     }    
  }
  rhopst.resize(ngloc);
  memset( (void*)&rhopst[0], 0, 2*ngloc*sizeof(double) );
  for ( int is = 0; is < atoms.nsp(); is++ )
  {
     complex<double> *s = &sf.sfac[is][0];
     for ( int ig = 0; ig < ngloc; ig++ )
     {
        const complex<double> sg = s[ig];
        rhopst[ig] += sg * rhops[is][ig];
     }
  }
  // compute total charge density:  electronic + ionic cores
  rhog.resize(ngloc);
  if ( s->wf.nspin() == 1 )
     for ( int ig = 0; ig < ngloc; ig++ ) 
        rhog[ig] = omega_inv * cd_.rhog[0][ig] + rhopst[ig];
  else 
     for ( int ig = 0; ig < ngloc; ig++ ) 
        rhog[ig] = omega_inv * ( cd_.rhog[0][ig] + cd_.rhog[1][ig] ) + rhopst[ig];

  // vhart_g = 4 * pi * (rhoel + rhops) * g2i
  v_g.resize(ngloc);
  for ( int ig = 0; ig < ngloc; ig++ ) 
     v_g[ig] = fpi * rhog[ig] * g2i[ig];
 
  // FT to v_r
  v_r.resize(ft_->np012loc());
  ft_->backward(&v_g[0],&v_r[0]);

  vector<double> vrtmp(ft_->np012loc());
  for (int j = 0; j < ft_->np012loc(); j++)
     vrtmp[j] = real(v_r[j]);
  

  if (format == "vmd") {
     const Context* wfctxt = s->wf.spincontext(0);
     const Context* vctxt = &cd_.vcontext();

     //ewd DEBUG
     assert(wfctxt->nprow() == vctxt->nprow());

     if (wfctxt->mycol() == 0) {

        for ( int i = 0; i < wfctxt->nprow(); i++ ) {
           if ( i == wfctxt->myrow() ) {
              int size = ft_->np012loc();
              wfctxt->isend(1,1,&size,1,0,0);
              if (size > 0)
                 wfctxt->dsend(size,1,&vrtmp[0],1,0,0);
           }
        }
        if ( wfctxt->mype() == 0 ) {
           vector<double> tmprecv(ft_->np012());
           int recvoffset = 0;

           D3vector a0 = s->wf.cell().a(0);
           D3vector a1 = s->wf.cell().a(1);
           D3vector a2 = s->wf.cell().a(2);
           const int np0 = ft_->np0();
           const int np1 = ft_->np1();
           const int np2 = ft_->np2();
           D3vector dft0 = a0/(double)np0;
           D3vector dft1 = a1/(double)np1;
           D3vector dft2 = a2/(double)np2;

           for ( int i = 0; i < wfctxt->nprow(); i++ ) {
              int size = 0;
              wfctxt->irecv(1,1,&size,1,i,0);
              if (size > 0)
                 wfctxt->drecv(size,1,&tmprecv[recvoffset],1,i,0);
              recvoffset += size;

              if (i==0) {
                 // write out VMD CUBE format header
                 os << "Qbox output in VMD CUBE format" << endl;
                 os << "  electrostatic potential" << endl;
                 
                 D3vector origin(0.0,0.0,0.0);
                 os << natoms_total << " " << origin << endl;
                 
                 // print FFT grid info
                 os << np0 << " " << dft0 << endl;
                 os << np1 << " " << dft1 << endl;
                 os << np2 << " " << dft2 << endl;

                 // print atom coordinates
                 for ( int is = 0; is < atoms.nsp(); is++ ) {
                    const int atnum = atoms.atomic_number(is);
                    double atnumd = (double)atnum;
                    for ( int ia = 0; ia < atoms.na(is); ia++ ) 
                       os << atnum << " " << atnumd << " " << tau0[is][3*ia] << " " << tau0[is][3*ia+1] << " " << tau0[is][3*ia+2] << endl;
                 }
              }
           }

           // write potential data to file
           int cnt = 0;
           for (int ii = 0; ii < np0; ii++) {
              ostringstream oss;
              oss.setf(ios::scientific,ios::floatfield);
              oss << setprecision(5);
              for (int jj = 0; jj < np1; jj++) {
                 for (int kk = 0; kk < np2; kk++) {
                    int index = ii + jj*np0 + kk*np0*np1;
                    oss << tmprecv[index] << " ";
                    cnt++;
                    if (cnt >= 6) {
                       cnt = 0;
                       oss << endl;
                    }
                 }
              }
              string tos = oss.str();
              os.write(tos.c_str(),tos.length());
           }
        }
     }
  }
  else // if format != vmd
  {
     if (ctxt_->mycol() == 0) {
        
        vector<double> rhortmp(ft_->np012loc());
        for (int j = 0; j < ft_->np012loc(); j++)
           rhortmp[j] = cd_.rhor[0][j];
        
        for ( int i = 0; i < ctxt_->nprow(); i++ ) {
           if ( i == ctxt_->myrow() ) {
              int size = ft_->np012loc();
              //cout << " process " << ctxt_->mype() << " sending block " << i
              //     << " of density to task 0, size = " << size << endl;
              ctxt_->isend(1,1,&size,1,0,0);
              if (size > 0)
                 ctxt_->dsend(size,1,&vrtmp[0],1,0,0);
           }
        }
        if ( ctxt_->oncoutpe() ) {
           for ( int i = 0; i < ctxt_->nprow(); i++ ) {
              int size = 0;
              ctxt_->irecv(1,1,&size,1,i,0);
              //int istart = cd_.vft.np0() * cd_.vft.np1() * cd_.vft.np2_first(i);
              //cout << " process " << ctxt_->mype() << " receiving block " << i
              //     << " of density on task 0, size = " << size << endl;
              if (size > 0)
                 ctxt_->drecv(size,1,&vrtmp[0],1,i,0);

              // headers for visualization formats
              if (i==0) {
                 if (format == "molmol") {
                    D3vector a0 = s->wf.cell().a(0);
                    D3vector a1 = s->wf.cell().a(1);
                    D3vector a2 = s->wf.cell().a(2);
                    if( a0.y != 0.0 || a0.z != 0.0 || a1.x != 0.0 || a1.z != 0.0 ||
                        a2.x != 0.0 || a2.y != 0.0 )
                       os << "Error writing header:  Molmol isosurface requires rectangular box!" << endl;
                    double dx = a0.x/(double)ft_->np0();
                    double dy = a1.y/(double)ft_->np1();
                    double dz = a2.z/(double)ft_->np2();
                    //    origin    npoints       grid spacing
                    os << "0.0 " << ft_->np0() << " " << dx << endl;
                    os << "0.0 " << ft_->np1() << " " << dy << endl;
                    os << "0.0 " << ft_->np2() << " " << dz << endl;
                 }
                 else if (format == "gopenmol") {
                    D3vector a0 = s->wf.cell().a(0);
                    D3vector a1 = s->wf.cell().a(1);
                    D3vector a2 = s->wf.cell().a(2);
                    if( a0.y != 0.0 || a0.z != 0.0 || a1.x != 0.0 || a1.z != 0.0 ||
                        a2.x != 0.0 || a2.y != 0.0 )
                       os << "Error writing header:  gOpenMol isosurface requires rectangular box!" << endl;
                    os << "3 200" << endl;  //ewd copying this from another utility, not sure what it means
                    os << ft_->np0() << " " << ft_->np1() << " " << ft_->np2() << endl;
                    os << "0.0 " << a2.z << endl;
                    os << "0.0 " << a1.y << endl;
                    os << "0.0 " << a0.x << endl;
                 }
              }
              ostringstream oss;
              for (int j=0; j<size; j++) 
                 oss << vrtmp[j] << endl;
              string tos = oss.str();
              os.write(tos.c_str(),tos.length());
           }
        }
     }
  }
  os.close();

  return 0;
}
