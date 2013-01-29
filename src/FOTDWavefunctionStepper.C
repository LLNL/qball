////////////////////////////////////////////////////////////////////////////////
//
// FOTDWavefunctionStepper.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: FOTDWavefunctionStepper.C,v 1.8 2011-06-02 15:56:19 schleife Exp $

#include "FOTDWavefunctionStepper.h"
#include "Wavefunction.h"
#include "SlaterDet.h"
#include "Sample.h"
#include <iostream>
#include <deque>
using namespace std;

// AS: this class provides the FOTD wave function stepper
// AS: FOTD enables fourth-order wave-function propagation according to
// AS: |psi(t+2tddt)> = 8|psi(t+tddt)> + 12*i*tddt*|H psi(t)> - 8|psi(t-tddt)> + |psi(t-2tddt)> 
// AS:                = 8*wf2_         + 12*i*tddt*dwf        - 8*next_wf_     + wf1_

////////////////////////////////////////////////////////////////////////////////
FOTDWavefunctionStepper::FOTDWavefunctionStepper(Wavefunction& wf, double tddt,
 TimerMap& tmap, deque<Wavefunction> *wfdeque) : tddt_(tddt), WavefunctionStepper(wf,tmap), wfdeque_(wfdeque)
{ }

// AS: the following block worked without dequeue
// FOTDWavefunctionStepper::FOTDWavefunctionStepper(Wavefunction& wf, double tddt,
// TimerMap& tmap, Wavefunction* next_wf, Wavefunction* wf1, Wavefunction* wf2) :
//  tddt_(tddt), WavefunctionStepper(wf,tmap), next_wf_(next_wf), wf1_(wf1), wf2_(wf2)
// { }

////////////////////////////////////////////////////////////////////////////////
void FOTDWavefunctionStepper::update(Wavefunction& dwf)
{
  // AS: create some space for the new wave function
  (*wfdeque_).push_back( wf_ );

  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
  {
    for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
    {
      tmap_["sd_update_wf"].start();

      ((*wfdeque_)[4]).sd(ispin,ikp)->c() = ((*wfdeque_)[0]).sd(ispin,ikp)->c();
      (((*wfdeque_)[4]).sd(ispin,ikp)->c()).axpy( 8.0, ((*wfdeque_)[3]).sd(ispin,ikp)->c() );
      (((*wfdeque_)[4]).sd(ispin,ikp)->c()).axpy(-8.0, ((*wfdeque_)[1]).sd(ispin,ikp)->c() );
      (((*wfdeque_)[4]).sd(ispin,ikp)->c()).axpy(12.0*tddt_*(complex<double>(0,1)),dwf.sd(ispin,ikp)->c());

      tmap_["sd_update_wf"].stop();

      // AS: the following block worked without dequeue
      // tmap_["sd_update_wf"].start();

      // ComplexMatrix ctmp_(wf_.sd(ispin,ikp)->c());
      // ctmp_ = (*wf1_).sd(ispin,ikp)->c();
      // ctmp_.axpy(8.0,(*wf2_).sd(ispin,ikp)->c());
      // ctmp_.axpy(-8.0,(*next_wf_).sd(ispin,ikp)->c());
      // ctmp_.axpy(12.0*tddt_*(complex<double>(0,1)),dwf.sd(ispin,ikp)->c());

      // (*wf1_).sd(ispin,ikp)->c()=(*next_wf_).sd(ispin,ikp)->c();
      // (*next_wf_).sd(ispin,ikp)->c()=wf_.sd(ispin,ikp)->c();
      // wf_.sd(ispin,ikp)->c()=(*wf2_).sd(ispin,ikp)->c();
      // (*wf2_).sd(ispin,ikp)->c()=ctmp_;

      // tmap_["sd_update_wf"].stop();
    }
  }

  // AS: trash the oldest wave function
  (*wfdeque_).pop_front();

  // AS: copy back to wf, because EnergyFunctional is applied based on wf
  wf_=(*wfdeque_)[2];

}
