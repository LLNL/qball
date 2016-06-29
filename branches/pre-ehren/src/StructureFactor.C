////////////////////////////////////////////////////////////////////////////////
//
// StructureFactor.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: StructureFactor.C,v 1.3 2009/04/05 15:27:36 draeger1 Exp $

#include "StructureFactor.h"
#include "Basis.h"
#include "UnitCell.h"
#include <cstring>
#include <cassert>

////////////////////////////////////////////////////////////////////////////////
void StructureFactor::init(const vector<vector<double> >& tau, 
  const Basis& basis)
{
  _k0min = basis.idxmin(0);
  _k1min = basis.idxmin(1);
  _k2min = basis.idxmin(2);
  _k0max = basis.idxmax(0);
  _k1max = basis.idxmax(1);
  _k2max = basis.idxmax(2);
  _k0range = _k0max - _k0min + 1;
  _k1range = _k1max - _k1min + 1;
  _k2range = _k2max - _k2min + 1;
  
  // get dimensions of tau[nsp][3*na[is]]
  
  _nsp = tau.size();
  _na.resize(_nsp);

  cos0.resize(_nsp);
  cos1.resize(_nsp);
  cos2.resize(_nsp);
  sin0.resize(_nsp);
  sin1.resize(_nsp);
  sin2.resize(_nsp);
  sfac.resize(_nsp);
  
  _ng = basis.localsize();
  
  for ( int is = 0; is < _nsp; is++ )
  {
    assert( tau[is].size() % 3 == 0 );
    _na[is] = tau[is].size() / 3; // tau[nsp][3*na[is]]
    cos0[is].resize(_na[is]*_k0range);
    cos1[is].resize(_na[is]*_k1range);
    cos2[is].resize(_na[is]*_k2range);
    sin0[is].resize(_na[is]*_k0range);
    sin1[is].resize(_na[is]*_k1range);
    sin2[is].resize(_na[is]*_k2range);
    sfac[is].resize(_ng);
  }
}

////////////////////////////////////////////////////////////////////////////////
void StructureFactor::update(const vector<vector<double> >& tau, 
  const Basis& basis)
{
  // it is assumed that the dimensions of tau and the basis have 
  // not changed since the last call to StructureFactor::init

  // check that number of species has not changed  
  assert(tau.size() == _nsp);  
  assert(basis.localsize() == _ng);
  
  const int * const idx = basis.idx_ptr();
  
  const UnitCell& cell = basis.cell();
  const D3vector b0 = cell.b(0);
  const D3vector b1 = cell.b(1);
  const D3vector b2 = cell.b(2);

  for ( int is = 0; is < _nsp; is++ )
  {
    assert( 3 * _na[is] == tau[is].size() );
    memset( (void*)&sfac[is][0], 0, 2*_ng*sizeof(double) );
    
    for ( int ia = 0; ia < _na[is]; ia++ )
    {
      double *c0 = cos0_ptr(is,ia);
      double *c1 = cos1_ptr(is,ia);
      double *c2 = cos2_ptr(is,ia);
      double *s0 = sin0_ptr(is,ia);
      double *s1 = sin1_ptr(is,ia);
      double *s2 = sin2_ptr(is,ia);
      
      const double * const tauptr = &tau[is][3*ia];
      
      const D3vector t(tauptr[0],tauptr[1],tauptr[2]);
      
      /* x direction */
      const double fac0 = b0 * t;

      for ( int i = _k0min; i < _k0max+1; i++ )
      {
        const double arg = i * fac0;
        c0[i] = cos(arg);
        s0[i] = sin(arg);
      }

      /* y direction */
      const double fac1 = b1 * t;

      for ( int i = _k1min; i < _k1max+1; i++ )
      {
        const double arg = i * fac1;
        c1[i] = cos(arg);
        s1[i] = sin(arg);
      }

      /* z direction */
      const double fac2 = b2 * t;

      for ( int i = _k2min; i < _k2max+1; i++ )
      {
        const double arg = i * fac2;
        c2[i] = cos(arg);
        s2[i] = sin(arg);
      }
      
      // compute sfac[is][i] 
 
      for ( int i = 0; i < _ng; i++ )
      {
        const int iii = i+i+i;
        const int kx = idx[iii];
        const int ky = idx[iii+1];
        const int kz = idx[iii+2];
        
        const double cos_a = c0[kx];
        const double cos_b = c1[ky];
        const double cos_c = c2[kz];
        
        const double sin_a = s0[kx];
        const double sin_b = s1[ky];
        const double sin_c = s2[kz];
        
        // Next line: exp(-i*gr) = 
        // (cos_a - I sin_a)*(cos_b - I sin_b)*(cos_c - I sin_c)
        sfac[is][i] += complex<double>(
        cos_a*cos_b*cos_c - sin_a*sin_b*cos_c - 
        sin_a*cos_b*sin_c - cos_a*sin_b*sin_c,
        sin_a*sin_b*sin_c - sin_a*cos_b*cos_c - 
        cos_a*sin_b*cos_c - cos_a*cos_b*sin_c );         
      }
    }
  }
}
