////////////////////////////////////////////////////////////////////////////////
//
// Preconditioner.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Preconditioner.C,v 1.7 2010/01/07 18:01:48 draeger1 Exp $

#include "Preconditioner.h"
#include "EnergyFunctional.h"
#include "Wavefunction.h"
#include "Sample.h"
#include "Basis.h"
#include "SlaterDet.h"
#include "ConfinementPotential.h"

////////////////////////////////////////////////////////////////////////////////
Preconditioner::Preconditioner(const Sample& s, const Wavefunction& wf, const EnergyFunctional& ef) : 
    s_(s), wf_(wf), ef_(ef)
{
  update();
}

////////////////////////////////////////////////////////////////////////////////
void Preconditioner::update(void)
{
  // reinitialize preconditioner
  bool use_confinement = s_.ctrl.ecuts > 0.0;
  // If ecutprec is zero, use ecut
  const double ecutpr = s_.ctrl.ecutprec> 0.0 ? s_.ctrl.ecutprec : wf_.ecut();
  
  diag_.resize(wf_.nspin());
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ ) {
    if (wf_.spinactive(ispin)) {
      diag_[ispin].resize(wf_.nkp());
      for ( int ikp=0; ikp<wf_.nkp(); ikp++) {
        if (wf_.kptactive(ikp)) {
          assert(wf_.sd(ispin,ikp) != 0);
          // Only resize and initialize diag_ if ikp is active on this task
          const Basis& basis = wf_.sd(ispin,ikp)->basis();
          const int ngwloc = basis.localsize();
          diag_[ispin][ikp].resize(ngwloc);
          const double *kpg2_ptr = basis.kpg2_ptr();
          
          if ( use_confinement ) {
            const valarray<double>& fstress = ef_.confpot(ispin,ikp)->fstress();
            for ( int ig = 0; ig < ngwloc; ig++ ) {
              double e = 0.5 * ( kpg2_ptr[ig] + fstress[ig] );
              diag_[ispin][ikp][ig] = ( e < ecutpr ) ? 0.5 / ecutpr : 0.5 / e;
            }
          }
          else {
            for ( int ig = 0; ig < ngwloc; ig++ ) {
              double e = 0.5 * kpg2_ptr[ig];
              diag_[ispin][ikp][ig] = ( e < ecutpr ) ? 0.5 / ecutpr : 0.5 / e;
            }
          }
        }
      }
    }
  }
}
