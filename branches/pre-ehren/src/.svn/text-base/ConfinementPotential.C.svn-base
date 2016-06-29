////////////////////////////////////////////////////////////////////////////////
//
// ConfinementPotential.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ConfinementPotential.C,v 1.2 2007/07/24 16:40:14 draeger1 Exp $

#include "ConfinementPotential.h"
#include "Basis.h"

////////////////////////////////////////////////////////////////////////////////
ConfinementPotential::ConfinementPotential(double ecuts, double facs, 
  double sigmas, const Basis& basis): ecuts_(ecuts), facs_(facs),
  sigmas_(sigmas), basis_(basis)
{
  const int ngwloc = basis_.localsize();
  fstress_.resize(ngwloc);
  dfstress_.resize(ngwloc);
  update();
}

////////////////////////////////////////////////////////////////////////////////
void ConfinementPotential::update(void)
{
  // recompute confinement potential and its derivative
  const double sigmas_inv = 1.0 / sigmas_;
  const int ngwloc = basis_.localsize();
  
  for ( int ig = 0; ig < ngwloc; ig++ )
  {
    const double gsq = basis_.kpg2(ig);
    // Next line: 0.5 from 1/2m
    const double arg = ( 0.5 * gsq - ecuts_ ) * sigmas_inv;
    // Next lines: fp = fermi(arg);
    // fm = fermi(-arg) = 1 - fp;
    double fm,fp;
    if ( arg > 50.0 )
    {
      fm = 1.0;
      fp = 0.0;
    }
    else if ( arg < -50.0 )
    {
      fm = 0.0;
      fp = 1.0;
    }
    else
    {
      fp = 1.0 / ( 1.0 + exp ( arg ) );
      fm = 1.0 - fp;
    }

    // f(G) = facs * ( 1 - fermi( (G^2-ecuts)/sigmas ) )
    // fg = f(G)
    const double fg = facs_ * fm;
    // gfgp = G f'(G)
    const double gfgp = gsq * fg * fp * sigmas_inv;
    
    // fstress[ig] = (k+G)^2 * f(G)
    fstress_[ig] = gsq * fg;
    
    // dfstress =  2 f(G) + G * f'(G)
    dfstress_[ig] = 2.0 * fg + gfgp;

    // ekin = sum_G |c_G|^2  G^2
    // econf = sum_G |c_G|^2 fstress[G]
    // stress_ekin_ij = (1/Omega) sum_G |c_G|^2 * 2 * G_i * G_j
    // stress_econf_ij = (1/Omega) sum_G |c_G|^2 * dfstress[G] * G_i * G_j
  }
}
