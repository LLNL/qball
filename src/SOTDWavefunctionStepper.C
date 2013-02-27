////////////////////////////////////////////////////////////////////////////////
//
// SOTDWavefunctionStepper.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SOTDWavefunctionStepper.C,v 1.8 2011-05-03 15:56:19 schleife Exp $

#include "SOTDWavefunctionStepper.h"
#include "Wavefunction.h"
#include "SlaterDet.h"
#include "Sample.h"
#include <iostream>
#include <deque>
using namespace std;

// AS: this class provides the SOTD wave function stepper
// AS: SOTD enables second-order wave-function propagation according to |psi(t+tddt)> = |psi(t-tddt)> - 2*i*tddt*|H psi(t)>

////////////////////////////////////////////////////////////////////////////////
SOTDWavefunctionStepper::SOTDWavefunctionStepper(Wavefunction& wf, double tddt,
 TimerMap& tmap, deque<Wavefunction*> *wfdeque) : tddt_(tddt), WavefunctionStepper(wf,tmap), wfdeque_(wfdeque)
{ }

////////////////////////////////////////////////////////////////////////////////
void SOTDWavefunctionStepper::update(Wavefunction& dwf)
{
  // AS: create some space for the new wave function
  (*wfdeque_).push_back( &wf_ );

  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
  {
    for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
    {
      tmap_["sd_update_wf"].start();

      (*(*wfdeque_)[2]).sd(ispin,ikp)->c() = (*(*wfdeque_)[0]).sd(ispin,ikp)->c();
      ((*(*wfdeque_)[2]).sd(ispin,ikp)->c()).axpy(-2.0*tddt_*(complex<double>(0,1)),dwf.sd(ispin,ikp)->c());

      tmap_["sd_update_wf"].stop();
    }
  }

  // AS: trash the oldest wave function
  (*wfdeque_).pop_front();

  // AS: copy back to wf, because EnergyFunctional is applied based on wf
  wf_=*(*wfdeque_)[1];
}
