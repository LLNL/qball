////////////////////////////////////////////////////////////////////////////////
//
// SDWavefunctionStepper.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SDWavefunctionStepper.C,v 1.8 2009/09/25 23:18:11 draeger1 Exp $

#include "SDWavefunctionStepper.h"
#include "Wavefunction.h"
#include "SlaterDet.h"
#include <iostream>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
SDWavefunctionStepper::SDWavefunctionStepper(Wavefunction& wf, double alpha, TimerMap& tmap) : 
    WavefunctionStepper(wf,tmap),alpha_(alpha)
{}

////////////////////////////////////////////////////////////////////////////////
void SDWavefunctionStepper::update(Wavefunction& dwf) {
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ ) {
    if (wf_.spinactive(ispin)) {
      for ( int ikp=0; ikp<wf_.nkp(); ikp++) {
        if (wf_.kptactive(ikp)) {
          assert(wf_.sd(ispin,ikp) != 0);
          // c = c - dt2bye * hpsi
          tmap_["sd_update_wf"].start();
          wf_.sd(ispin,ikp)->c().axpy(-alpha_,dwf.sd(ispin,ikp)->c());
          tmap_["sd_update_wf"].stop();
        }
      }
    }
  }
}
