////////////////////////////////////////////////////////////////////////////////
//
// RunCmd.C:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: RunCmd.C,v 1.6 2010/05/12 20:05:25 draeger1 Exp $

#include "RunCmd.h"
#include<iostream>
using namespace std;
#include "BOSampleStepper.h"
#include "CPSampleStepper.h"

#include<ctime>
#include<cassert>

int RunCmd::action(int argc, char **argv)
{

  if ( argc < 2 || argc > 5)
  {
    if ( ui->oncoutpe() )
    {
      cout << " use: run [-atomic_density] niter" << endl;
      cout << "      run [-atomic_density] niter nitscf" << endl;
      cout << "      run [-atomic_density] niter nitscf nite" << endl;
      return 1;
    }
  }

  if (s->ctrl.timer_hit) {
    if ( ui->oncoutpe() )
      cout << " <!-- RunCmd: run_timer exceeded in previous run, skipping current run command. -->" << endl;
    return 0;
  }
  
  if ( s->wf.nst() == 0 )
  {
    if ( ui->oncoutpe() )
      cout << " <!-- RunCmd: no states, cannot run -->" << endl;
    return 1;
  }
  if ( s->wf.ecut() == 0.0 )
  {
    if ( ui->oncoutpe() )
      cout << " <!-- RunCmd: ecut = 0.0, cannot run -->" << endl;
    return 1;
  }
  
  SampleStepper* stepper;
  
  int iarg = 1;
  bool atomic_density = false;
  if ( !strcmp(argv[iarg],"-atomic_density") )
  {
    atomic_density = true;
    iarg++;
    argc--;
  }

  int niter = atoi(argv[iarg]);
  int nite = 1;
  int nitscf = 1;
  if ( argc == 3 )
  {
    // run niter nitscf
    nitscf = atoi(argv[iarg+1]);
  }
  else if ( argc == 4 )
  {
    // run niter nitscf nite
    nitscf = atoi(argv[iarg+1]);
    nite = atoi(argv[iarg+2]);
  }

  if (!s->wf.hasdata()) {
    s->wf.set_hasdata(true);
    if ( !atomic_density && ui->oncoutpe()) {
      cout << "<WARNING> Wavefunction has been neither loaded nor initialized! </WARNING>" << endl;
      cout << "<WARNING> Use randomize_wf or run -atomic_density to initialize.  </WARNING>" << endl;
    }
  }
  
  if ( s->ctrl.wf_dyn == "MD" )
    stepper = new CPSampleStepper(*s);
  else
    stepper = new BOSampleStepper(*s,nitscf,nite);
  
  assert(stepper!=0);
  
  if (ui->oncoutpe() )
    cout << "<run niter_ionic=\"" << niter << "\" niter_scf=\"" << nitscf << "\" niter_nonscf=\"" << nite << "\">" << endl;

  if ( atomic_density )
    stepper->initialize_density();

  s->wf.info(cout,"wavefunction");
  stepper->step(niter);
  if (ui->oncoutpe() )
    cout << "</run>" << endl;
  
  delete stepper;
  
  return 0;
}
