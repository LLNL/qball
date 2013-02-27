////////////////////////////////////////////////////////////////////////////////
//
// TDEULERWavefunctionStepper.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: TDEULERWavefunctionStepper.C,v 1.8 2011-03-31 15:56:19 schleife Exp $

#include "TDEULERWavefunctionStepper.h"
#include "Wavefunction.h"
#include "SlaterDet.h"
#include "Sample.h"
#include <iostream>
using namespace std;

// AS: this class provides the TDEULER wave function stepper

////////////////////////////////////////////////////////////////////////////////
TDEULERWavefunctionStepper::TDEULERWavefunctionStepper(Wavefunction& wf, double tddt,
  TimerMap& tmap) :
  tddt_(tddt), WavefunctionStepper(wf,tmap)
{}

////////////////////////////////////////////////////////////////////////////////
void TDEULERWavefunctionStepper::update(Wavefunction& dwf)
{
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
  {
    for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
    {
      tmap_["sd_update_wf"].start();
      // AS: TDEULER enables first-order wave-function propagation according to |psi(t+tddt)> = |psi(t)> - i*tddt*|H psi(t)>
      wf_.sd(ispin,ikp)->c().axpy(-tddt_*(complex<double>(0,1)),dwf.sd(ispin,ikp)->c());
      // AS: if there was an energy renormalization to prevent time propagation from blowing up
      // AS: it could be undone here; this is just an example:
      // wf_.sd(ispin,ikp)->c() *= exp( complex<double>(0,1) * tddt_ * (2.263) );
      // AS: for the Euler integrator it does not change anything, possibly we need it later
      tmap_["sd_update_wf"].stop();
    }
  }
}
