////////////////////////////////////////////////////////////////////////////////
//
// BLYPFunctional.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: BLYPFunctional.C,v 1.3 2005/11/29 18:59:47 draeger1 Exp $

#include <cmath>
#include <cassert>
#include "BLYPFunctional.h"
using namespace std;

////////////////////////////////////////////////////////////////////////////////
BLYPFunctional::BLYPFunctional(const vector<vector<double> > &rhoe) {
  _nspin = rhoe.size();
  if ( _nspin > 1 ) assert(rhoe[0].size() == rhoe[1].size());
  _np = rhoe[0].size();
 
  if ( _nspin == 1 )
  {
    _exc.resize(_np);
    _vxc1.resize(_np);
    _vxc2.resize(_np);
    _grad_rho[0].resize(_np);
    _grad_rho[1].resize(_np);
    _grad_rho[2].resize(_np);
    rho = &rhoe[0][0];
    grad_rho[0] = &_grad_rho[0][0];
    grad_rho[1] = &_grad_rho[1][0];
    grad_rho[2] = &_grad_rho[2][0];
    exc = &_exc[0];
    vxc1 = &_vxc1[0];
    vxc2 = &_vxc2[0];
  }
  else
  {
    // not implemented
    assert(false);
  }
}

////////////////////////////////////////////////////////////////////////////////
void BLYPFunctional::setxc(void) 
{
  if ( _np == 0 ) return;
  if ( _nspin == 1 )
  {
    assert( rho != 0 );
    assert( grad_rho[0] != 0 && grad_rho[1] != 0 && grad_rho[2] != 0 );
    assert( exc != 0 );
    assert( vxc1 != 0 );
    assert( vxc2 != 0 );
    for ( int i = 0; i < _np; i++ )
    {
      double grad = sqrt(grad_rho[0][i]*grad_rho[0][i] +
                         grad_rho[1][i]*grad_rho[1][i] +
                         grad_rho[2][i]*grad_rho[2][i] );
      excblyp(rho[i],grad,&exc[i],&vxc1[i],&vxc2[i]);
    }
  }
  else
  {
#if 0 // not implemented
    assert( rho_up != 0 );
    assert( rho_dn != 0 );
    assert( grad_rho_up[0] != 0 && grad_rho_up[1] != 0 && grad_rho_up[2] != 0 );
    assert( grad_rho_dn[0] != 0 && grad_rho_dn[1] != 0 && grad_rho_dn[2] != 0 );
    assert( exc_up != 0 );
    assert( exc_dn != 0 );
    assert( vxc1_up != 0 );
    assert( vxc1_dn != 0 );
    assert( vxc2_upup != 0 );
    assert( vxc2_updn != 0 );
    assert( vxc2_dnup != 0 );
    assert( vxc2_dndn != 0 );

    for ( int i = 0; i < _np; i++ )
    {
      double grx_up = grad_rho_up[0][i];
      double gry_up = grad_rho_up[1][i];
      double grz_up = grad_rho_up[2][i];
      double grx_dn = grad_rho_dn[0][i];
      double gry_dn = grad_rho_dn[1][i];
      double grz_dn = grad_rho_dn[2][i];
      double grx = grx_up + grx_dn;
      double gry = gry_up + gry_dn;
      double grz = grz_up + grz_dn;
      double grad_up = sqrt(grx_up*grx_up + gry_up*gry_up + grz_up*grz_up);
      double grad_dn = sqrt(grx_dn*grx_dn + gry_dn*gry_dn + grz_dn*grz_dn);
      double grad    = sqrt(grx*grx + gry*gry + grz*grz);
      excpbe_sp(rho_up[i],rho_dn[i],grad_up,grad_dn,grad,&exc_up[i],&exc_dn[i],
                &vxc1_up[i],&vxc1_dn[i],&vxc2_upup[i],&vxc2_dndn[i],
                &vxc2_updn[i], &vxc2_dnup[i]);
    }
#endif
  }
}
////////////////////////////////////////////////////////////////////////////////
void BLYPFunctional::excblyp(double rho, double grad, 
  double *exc, double *vxc1, double *vxc2)
{
  /* Becke exchange constants */
  const double third  = 1.0 / 3.0;
  const double fourthirds = 4.0 / 3.0;
  const double fivethirds = 5.0 / 3.0;
  const double beta=0.0042;
  const double ax = -0.7385587663820224058; /* -0.75*pow(3.0/pi,third) */
  const double axa = -0.9305257363490999; /* -1.5*pow(3.0/(4*pi),third) */

  /* LYP constants */
  const double a = 0.04918;
  const double b = 0.132;
  const double ab36 = a * b / 36.0;
  const double c = 0.2533;
  const double c_third = c / 3.0;
  const double d = 0.349;
  const double d_third = d / 3.0;
  const double cf = 2.87123400018819; /* (3/10)*pow(3*pi*pi,2/3) */
  const double cfb = cf * b;
  
  *exc = 0.0;
  *vxc1 = 0.0;
  *vxc2 = 0.0;

  if ( rho < 1.e-18  ) 
  {
    return;
  }

  /*
   * Becke's exchange
   * A.D.Becke, Phys.Rev. B38, 3098 (1988)
   */
   
  const double rha = 0.5 * rho;
  const double grada = 0.5 * grad;
  
  const double rha13 = pow ( rha, third );
  const double rha43 = rha * rha13;
  const double xa = grada / rha43;
  const double xa2 = xa*xa;
  const double asinhxa = asinh(xa);
  const double frac = 1.0 / ( 1.0 + 6.0 * beta * xa * asinhxa );
  const double ga = axa - beta * xa2 * frac;
  /* N.B. in next line, ex is the energy density, hence rh13 */
  const double ex = rha13 * ga;
  
  /* energy done, now the potential */
  const double gpa = ( 6.0*beta*beta*xa2 * ( xa/sqrt(xa2+1.0) - asinhxa ) - 2.0*beta*xa ) *
        frac*frac;
  const double vx1 = rha13 * fourthirds * ( ga - xa * gpa );
  const double vx2 = - 0.5 * gpa / grada;
  
  /*------------------------------------------------------------*/
  /* LYP correlation */
  /* Phys. Rev. B 37, 785 (1988). */
  /* next lines specialized to the unpolarized case */
  const double rh13 = pow ( rho, third );
  const double rhm13 = 1.0 / rh13;
  const double rhm43 = rhm13 / rho;
  const double e = exp ( - c * rhm13 );
  const double num = 1.0 + cfb * e;
  const double den = 1.0 + d * rhm13;
  const double deninv = 1.0 / den;
  const double cfrac = num * deninv;
  
  const double delta = rhm13 * ( c + d * deninv );  
  const double rhm53 = rhm43 * rhm13;
  const double t1 = e * deninv;
  const double t2 = rhm53;
  const double t3 = 6.0 + 14.0 * delta;
  
  const double g = ab36 * t1 * t2 * t3;
  
  /* next line, ec is the energy density, hence divide the energy by rho */
  const double ec = - a * cfrac + 0.25 * g * grad * grad / rho;
  
  /* energy done, now the potential */
  const double de = c_third * rhm43 * e;
  const double dnum = cfb * de;
  const double dden = - d_third * rhm43;
  const double dfrac = ( dnum * den - dden * num ) * deninv * deninv;
  
  const double ddelta = - third * rhm43 * ( c + d * deninv ) -
           rhm13 * d * dden * deninv * deninv;
  const double dt1 = de * deninv - e * dden * deninv * deninv;
  const double dt2 = - fivethirds * rhm53/rho;
  const double dt3 = 14.0 * ddelta;
  
  const double dg = ab36 * ( dt1 * t2 * t3 + t1 * dt2 * t3 + t1 * t2 * dt3 );
  
  const double vc1 = - a * ( cfrac + rho * dfrac ) + 0.25 * dg * grad * grad; 
  const double vc2 = -0.5 * g;
  
  *exc = ex + ec;
  *vxc1 = vx1 + vc1;
  *vxc2 = vx2 + vc2;
}

