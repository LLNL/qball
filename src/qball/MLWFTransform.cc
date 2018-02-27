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
// MLWFTransform.C
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#include <iostream>
#include <iomanip>
#include <complex>
#include <cassert>

#include "MLWFTransform.h"
#include "D3vector.h"
#include "Basis.h"
#include "SlaterDet.h"
#include <qball/UnitCell.h>
#include "jade.h"
#include <math/blas.h>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
MLWFTransform::MLWFTransform(const SlaterDet& sd) : sd_(sd),
cell_(sd.basis().cell()), ctxt_(sd.context()),  bm_(BasisMapping(sd.basis()))
{
  a_.resize(6);
  adiag_.resize(6);
  const int n = sd.c().n();
  const int nb = sd.c().nb();
  for ( int k = 0; k < 6; k++ )
  {
    a_[k] = new DoubleMatrix(ctxt_,n,n,nb,nb);
    adiag_[k].resize(n);
  }
  u_ = new DoubleMatrix(ctxt_,n,n,nb,nb);
}

////////////////////////////////////////////////////////////////////////////////
MLWFTransform::~MLWFTransform(void)
{
  for ( int k = 0; k < 6; k++ )
    delete a_[k];
  delete u_;
}

////////////////////////////////////////////////////////////////////////////////
void MLWFTransform::compute_transform(void)
{
  const int maxsweep = 50;
  const double tol = 1.e-8;

  SlaterDet sdcosx(sd_), sdsinx(sd_),
            sdcosy(sd_), sdsiny(sd_),
            sdcosz(sd_), sdsinz(sd_);
  const ComplexMatrix& c = sd_.c();
  ComplexMatrix& ccosx = sdcosx.c();
  ComplexMatrix& csinx = sdsinx.c();
  ComplexMatrix& ccosy = sdcosy.c();
  ComplexMatrix& csiny = sdsiny.c();
  ComplexMatrix& ccosz = sdcosz.c();
  ComplexMatrix& csinz = sdsinz.c();
  // proxy real matrices cr, cc, cs
  DoubleMatrix cr(c);
  DoubleMatrix ccx(ccosx);
  DoubleMatrix csx(csinx);
  DoubleMatrix ccy(ccosy);
  DoubleMatrix csy(csiny);
  DoubleMatrix ccz(ccosz);
  DoubleMatrix csz(csinz);
  vector<complex<double> > zvec(bm_.zvec_size()),
    zvec_cos(bm_.zvec_size()), zvec_sin(bm_.zvec_size()),
    ct(bm_.np012loc()), ct_cos(bm_.np012loc()), ct_sin(bm_.np012loc());

  for ( int i = 0; i < 6; i++ )
  {
    a_[i]->resize(c.n(), c.n(), c.nb(), c.nb());
    adiag_[i].resize(c.n());
  }
  u_->resize(c.n(), c.n(), c.nb(), c.nb());

  // loop over all local states
  const int np0 = bm_.np0();
  const int np1 = bm_.np1();
  const int np2 = bm_.np2();
  const int np01 = np0 * np1;
  const int np2loc = bm_.np2loc();
  const int nvec = bm_.nvec();
  for ( int n = 0; n < c.nloc(); n++ )
  {
    const complex<double>* f = c.cvalptr(n*c.mloc());
    complex<double>* fcx = ccosx.valptr(n*c.mloc());
    complex<double>* fsx = csinx.valptr(n*c.mloc());
    complex<double>* fcy = ccosy.valptr(n*c.mloc());
    complex<double>* fsy = csiny.valptr(n*c.mloc());
    complex<double>* fcz = ccosz.valptr(n*c.mloc());
    complex<double>* fsz = csinz.valptr(n*c.mloc());

    // direction z
    // map state to array zvec_
    bm_.vector_to_zvec(&f[0],&zvec[0]);

    for ( int ivec = 0; ivec < nvec; ivec++ )
    {
      const int ibase = ivec * np2;
      compute_sincos(np2,&zvec[ibase],&zvec_cos[ibase],&zvec_sin[ibase]);
    }
    // map back zvec_cos to sdcos and zvec_sin to sdsin
    bm_.zvec_to_vector(&zvec_cos[0],&fcz[0]);
    bm_.zvec_to_vector(&zvec_sin[0],&fsz[0]);

    // x direction
    // map zvec to ct
    bm_.transpose_fwd(&zvec[0],&ct[0]);

    for ( int iz = 0; iz < np2loc; iz++ )
    {
      for ( int iy = 0; iy < np1; iy++ )
      {
        const int ibase = iz * np01 + iy * np0;
        compute_sincos(np0,&ct[ibase],&ct_cos[ibase],&ct_sin[ibase]);
      }
    }
    // transpose back ct_cos to zvec_cos
    bm_.transpose_bwd(&ct_cos[0],&zvec_cos[0]);
    // transpose back ct_sin to zvec_sin
    bm_.transpose_bwd(&ct_sin[0],&zvec_sin[0]);

    // map back zvec_cos to sdcos and zvec_sin to sdsin
    bm_.zvec_to_vector(&zvec_cos[0],&fcx[0]);
    bm_.zvec_to_vector(&zvec_sin[0],&fsx[0]);

    // y direction
    vector<complex<double> > c_tmp(np1),ccos_tmp(np1),csin_tmp(np1);
    int one = 1;
    int len = np1;
    int stride = np0;
    for ( int iz = 0; iz < np2loc; iz++ )
    {
      for ( int ix = 0; ix < np0; ix++ )
      {
        const int ibase = iz * np01 + ix;
        zcopy(&len,&ct[ibase],&stride,&c_tmp[0],&one);
        compute_sincos(np1,&c_tmp[0],&ccos_tmp[0],&csin_tmp[0]);
        zcopy(&len,&ccos_tmp[0],&one,&ct_cos[ibase],&stride);
        zcopy(&len,&csin_tmp[0],&one,&ct_sin[ibase],&stride);
      }
    }
    // transpose back ct_cos to zvec_cos
    bm_.transpose_bwd(&ct_cos[0],&zvec_cos[0]);
    // transpose back ct_sin to zvec_sin
    bm_.transpose_bwd(&ct_sin[0],&zvec_sin[0]);

    // map back zvec_cos and zvec_sin
    bm_.zvec_to_vector(&zvec_cos[0],&fcy[0]);
    bm_.zvec_to_vector(&zvec_sin[0],&fsy[0]);
  }

  // dot products a_[0] = <cos x>, a_[1] = <sin x>
  a_[0]->gemm('t','n',2.0,cr,ccx,0.0);
  a_[0]->ger(-1.0,cr,0,ccx,0);
  a_[1]->gemm('t','n',2.0,cr,csx,0.0);
  a_[1]->ger(-1.0,cr,0,csx,0);

  // dot products a_[2] = <cos y>, a_[3] = <sin y>
  a_[2]->gemm('t','n',2.0,cr,ccy,0.0);
  a_[2]->ger(-1.0,cr,0,ccy,0);
  a_[3]->gemm('t','n',2.0,cr,csy,0.0);
  a_[3]->ger(-1.0,cr,0,csy,0);

  // dot products a_[4] = <cos z>, a_[5] = <sin z>
  a_[4]->gemm('t','n',2.0,cr,ccz,0.0);
  a_[4]->ger(-1.0,cr,0,ccz,0);
  a_[5]->gemm('t','n',2.0,cr,csz,0.0);
  a_[5]->ger(-1.0,cr,0,csz,0);

  int nsweep = jade(maxsweep,tol,a_,*u_,adiag_);
}

////////////////////////////////////////////////////////////////////////////////
void MLWFTransform::compute_sincos(const int n, const complex<double>* f,
  complex<double>* fc, complex<double>* fs)
{
  // fc[i] =     0.5 * ( f[i-1] + f[i+1] )
  // fs[i] = (0.5/i) * ( f[i-1] - f[i+1] )

  // i = 0
  complex<double> zp = f[n-1];
  complex<double> zm = f[1];
  fc[0] = 0.5 * ( zp + zm );
  complex<double> zdiff = zp - zm;
  fs[0] = 0.5 * complex<double>(imag(zdiff),-real(zdiff));
  for ( int i = 1; i < n-1; i++ )
  {
    const complex<double> zzp = f[i-1];
    const complex<double> zzm = f[i+1];
    fc[i] = 0.5 * ( zzp + zzm );
    const complex<double> zzdiff = zzp - zzm;
    fs[i] = 0.5 * complex<double>(imag(zzdiff),-real(zzdiff));
  }
  // i = n-1
  zp = f[n-2];
  zm = f[0];
  fc[n-1] = 0.5 * ( zp + zm );
  zdiff = zp - zm;
  fs[n-1] = 0.5 * complex<double>(imag(zdiff),-real(zdiff));
}

////////////////////////////////////////////////////////////////////////////////
D3vector MLWFTransform::center(int i)
{
  assert(i>=0 && i<sd_.nst());
  const double cx = adiag_[0][i];
  const double sx = adiag_[1][i];
  const double cy = adiag_[2][i];
  const double sy = adiag_[3][i];
  const double cz = adiag_[4][i];
  const double sz = adiag_[5][i];
  // Next lines: M_1_PI = 1.0/pi
  const double itwopi = 1.0 / ( 2.0 * M_PI );
  const double t0 = itwopi * atan2(sx,cx);
  const double t1 = itwopi * atan2(sy,cy);
  const double t2 = itwopi * atan2(sz,cz);
  const double x = t0*cell_.a(0).x + t1*cell_.a(1).x + t2*cell_.a(2).x;
  const double y = t0*cell_.a(0).y + t1*cell_.a(1).y + t2*cell_.a(2).y;
  const double z = t0*cell_.a(0).z + t1*cell_.a(1).z + t2*cell_.a(2).z;

  return D3vector(x,y,z);
}

////////////////////////////////////////////////////////////////////////////////
double MLWFTransform::spread2(int i, int j)
{
  assert(i>=0 && i<sd_.nst());
  assert(j>=0 && j<3);
  const double c = adiag_[2*j][i];
  const double s = adiag_[2*j+1][i];
  // Next line: M_1_PI = 1.0/pi
  const double fac = 1.0 / length(cell_.b(j));
  return fac*fac * ( 1.0 - c*c - s*s );
}

////////////////////////////////////////////////////////////////////////////////
double MLWFTransform::spread2(int i)
{
  assert(i>=0 & i<sd_.nst());
  return spread2(i,0) + spread2(i,1) + spread2(i,2);
}

////////////////////////////////////////////////////////////////////////////////
double MLWFTransform::spread(int i)
{
  return sqrt(spread2(i));
}

////////////////////////////////////////////////////////////////////////////////
double MLWFTransform::spread2(void)
{
  double sum = 0.0;
  for ( int i = 0; i < sd_.nst(); i++ )
    sum += spread2(i);
  return sum;
}

////////////////////////////////////////////////////////////////////////////////
double MLWFTransform::spread(void)
{
  return sqrt(spread2());
}

////////////////////////////////////////////////////////////////////////////////
D3vector MLWFTransform::dipole(void)
{
  // total electronic dipole
  D3vector sum(0.0,0.0,0.0);
  for ( int i = 0; i < sd_.nst(); i++ )
    sum -= sd_.occ(i) * center(i);
  return sum;
}

////////////////////////////////////////////////////////////////////////////////
void MLWFTransform::apply_transform(SlaterDet& sd)
{
  // proxy double matrix c
  DoubleMatrix c(sd.c());
  DoubleMatrix cp(c);
  c.gemm('n','n',1.0,cp,*u_,0.0);
}
