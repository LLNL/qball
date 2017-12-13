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
// qb.C
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
using namespace std;

#include <config.h>

#include <sys/utsname.h>
#include <unistd.h>
#include <cstdlib>
#include <fstream>
#if AIX 
#include<filehdr.h>
#endif
#ifdef USE_APC
#include "apc.h"
#endif

#include "omp.h"
#include "isodate.h"
#include "release.h"
#include "qbox_xmlns.h"
#include "Context.h"
#include "UserInterface.h"
#include "Sample.h"
#include "Timer.h"

#include "AtomCmd.h"
#include "MMAtomCmd.h"
#include "HelpCmd.h"
#include "ListAtomsCmd.h"
#include "ListSpeciesCmd.h"
#include "LoadCmd.h"
#include "PrintCmd.h"
#include "PromoteOccCmd.h"
#include "QuitCmd.h"
#include "RandomizeWfCmd.h"
#include "RandomizeRealWfCmd.h"
#include "RandomizeVelCmd.h"
#include "RunCmd.h"
#include "MDSaveCmd.h"
#include "SaveCmd.h"
#include "SavesysCmd.h"
#include "SavedenCmd.h"
#include "SaveESPCmd.h"
#include "SetCmd.h"
#include "ShiftWFCmd.h"
#include "WFPhaseRealCmd.h"
#include "SpeciesCmd.h"
#include "MMSpeciesCmd.h"
#include "EmpiricalPotentialCmd.h"
#include "StatusCmd.h"
#include "KpointCmd.h"
#include "SymmetryCmd.h"
#include "ParOptCmd.h"
#include "LockCmd.h"
#include "UnlockCmd.h"
#include "SetVelCmd.h"
#include "ComputeMLWFCmd.h"
#include "AngleCmd.h"
#include "ConstraintCmd.h"
#include "DistanceCmd.h"
#include "FoldInWsCmd.h"
#include "StrainCmd.h"
#include "TorsionCmd.h"
#include "ResetVcmCmd.h"
#include "ListConstraintsCmd.h"
#include "PlotCmd.h"
#include "CoordinatesCmd.h"

#include "AtomsDyn.h"
#include "Cell.h"
#include "CellDyn.h"
#include "CellLock.h"
#include "CellMass.h"
#include "ChargeMixing.h"
#include "ChargeMixCoeff.h"
#include "ChargeMixRcut.h"
#include "ChargeMixNdim.h"
#include "Debug.h"
#include "Ecut.h"
#include "Ecutprec.h"
#include "Ecutden.h"
#include "Ecuts.h"
#include "Emass.h"
#include "ExtStress.h"
#include "Smearing.h"
#include "SmearingWidth.h"
#include "FermiTemp.h"
#include "Force_Complex_WF.h"
#include "Non_Selfconsistent_Energy_Output.h"
#include "TDDt.h"
#include "NA_overlaps.h"
#include "Dt.h"
#include "Nempty.h"
#include "Nrowmax.h"
#include "RefCell.h"
#include "Spin.h"
#include "Stress.h"
#include "Thermostat.h"
#include "ThresholdScf.h"
#include "ThresholdForce.h"
#include "ThresholdStress.h"
#include "ThTemp.h"
#include "ThTime.h"
#include "ThWidth.h"
#include "CenterOfMass.h"
#include "WfDiag.h"
#include "WfDyn.h"
#include "WfExtrap.h"
#include "WFPhaseRealVar.h"
#include "Xc.h"
#include "Nparallelkpts.h"
#include "Nkpoints.h"
#include "IPrint.h"
#include "CellStepFreq.h"
#include "EnthalpyPressure.h"
#include "EnthalpyThreshold.h"
#include "HugoniostatVar.h"
#include "HugDeltaTemp.h"
#include "HugFreq.h"
#include "RunTimer.h"
#include "HubbardU.h"
#include "Memory.h"
#include "MDIter.h"
#include "profile.h"
#include "MatrixLoc.h"
#include "Pblock.h"
#include "SaveFreq.h"
#include "SaveDenFreq.h"
#include "SaveWfFreq.h"
#include "NetCharge.h"
#include "EsmBC.h"
#include "EsmW.h"
#include "FcpThermostat.h"
#include "FcpThTemp.h"
#include "FcpThTime.h"
#include "FcpThWidth.h"
#include "FcpPmass.h"
#include "FcpMu.h"
#include "VdW.h"

#ifdef HAVE_BGQLIBS
#include <spi/include/kernel/process.h>
#include <spi/include/kernel/location.h>
#endif

#ifdef USE_OLD_CTF
#include "cyclopstf.h"
#endif

#ifdef USE_JAGGEMM
extern "C" int setup_grid();
#endif

int main(int argc, char **argv, char **envp)
{
  
  //Check command line arguments for -c (print configuration options)
  if ( argc == 2 && argv[1] == string("-c") ){
#ifdef USE_MPI
    cout << "mpi ";
#endif
#ifdef HAVE_OPENMP
    cout << "openmp ";
#endif
    cout << endl;
    return 0;
  }

  Timer tm;
  tm.start();

#if USE_MPI
  MPI_Init(&argc,&argv);
#endif
#if USE_APC
  ApcInit();
#endif
#ifdef USE_MPIP
  MPI_Pcontrol(0);
#endif
#ifdef USE_OLD_CTF
  {
    int myRank,numPes;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numPes);
    CTF_init(MPI_COMM_WORLD, MACHINE_BGQ, myRank, numPes); 
    CTF_init_complex(MPI_COMM_WORLD, MACHINE_BGQ, myRank, numPes); 
    //CTF_init(MPI_COMM_WORLD, myRank, numPes); 
    //CTF_init_complex(MPI_COMM_WORLD, myRank, numPes); 
  }    
#endif
#ifdef TAU
  int myRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  TAU_PROFILE_TIMER(timer, "main", "int (int, char**)", TAU_USER);
  TAU_PROFILE_START(timer);
  TAU_PROFILE_INIT(argc, argv);
  TAU_PROFILE_SET_NODE(myRank);
  TAU_PROFILE_SET_CONTEXT(0);
#endif

  {

  Context ctxt;
  ctxt.set_coutpe(0);
  
  if ( ctxt.oncoutpe() )
  {
    cout << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
    cout << "<!--\n\n";
    cout << "                   ___________________________\n";
    cout << "                   |                         |\n";
    cout << "                   | " << setw(22) << left << release() << "  |\n";
    cout << "                   |                         |\n";
    cout << "                   |                         |\n";
    cout << "                   |                         |\n";
    cout << "                   |                         |\n";
    cout << "                   |                         |\n";
    cout << "                   |                         |\n";
    cout << "                   |                         |\n";
    cout << "                   |                         |\n";
    cout << "                   |      Lawrence Livermore |\n";
    cout << "                   |     National Laboratory |\n";
    cout << "                   |                         |\n";
    cout << "                   | Copyright (c) 2003-2016 |\n";     
    cout << "                   |_________________________|\n\n";
    cout << "-->\n";
    cout << "<qbox:simulation xmlns:qbox=\"" << qbox_xmlns() << "\">" << endl;
    cout << "<release> " << release() << " </release>" << endl;
    cout << "<npes> " << ctxt.size() << " </npes>" << endl;
    int nthreads = omp_get_max_threads();
    cout << "<nthreads> " << nthreads << " </nthreads>" << endl;

    // Identify executable name, checksum, size and link date
    if ( getlogin() != 0 ) 
      cout << "<user> " << getlogin() << " </user>" << endl;
#if AIX || OSF1
    // read filehdr for link time
    filehdr hdr;
    FILE *fx = fopen(argv[0],"r");
    if ( fx != 0 ) {
      size_t sz = fread((void*)&hdr,sizeof(filehdr),1,fx);
      fclose(fx);
      string s = ctime((time_t*)&hdr.f_timdat);
      cout << "<linktime> " << s << " </linktime>" << endl;
    }
#endif
  
    // Identify platform
    {
       struct utsname un;
       uname (&un);
       cout << "<sysname> " << un.sysname << " </sysname>" << endl;
       cout << "<nodename> " << un.nodename << " </nodename>" << endl;
    }
    
    cout << "<start_time> " << isodate() << " </start_time>" << endl;

  }

#ifdef HAVE_BGQLIBS
  Personality_t pers;
  Kernel_GetPersonality(&pers, sizeof(pers));
  const int nTDim = 5;
  vector<int> torusdim(nTDim);
  torusdim[0] = pers.Network_Config.Anodes;
  torusdim[1] = pers.Network_Config.Bnodes;
  torusdim[2] = pers.Network_Config.Cnodes;
  torusdim[3] = pers.Network_Config.Dnodes;
  torusdim[4] = pers.Network_Config.Enodes;
  if ( ctxt.oncoutpe() )
     cout << "BG/Q torus dimensions:  " << torusdim[0] << " x " << torusdim[1] << " x " << torusdim[2]
          << " x " << torusdim[3] << " x " << torusdim[4] << endl;
#endif
  
  Sample* s = new Sample(ctxt);

  //store timing for run_timer
  s->ctrl.time_init = MPI_Wtime();
  s->ctrl.timer_hit = false;
  s->ctrl.timer_mdsavecmd = false;
  s->ctrl.timer_savecmd = false;
  s->ctrl.timer_savesyscmd = false;
  
  UserInterface* ui = new UserInterface(ctxt);

  ui->addCmd(new AtomCmd(s));
  ui->addCmd(new MMAtomCmd(s));
  ui->addCmd(new HelpCmd(s));
  ui->addCmd(new ListAtomsCmd(s));
  ui->addCmd(new ListSpeciesCmd(s));
  ui->addCmd(new LoadCmd(s));
  ui->addCmd(new PrintCmd(s));
  ui->addCmd(new PromoteOccCmd(s));
  ui->addCmd(new QuitCmd(s));
  ui->addCmd(new RandomizeWfCmd(s));
  ui->addCmd(new RandomizeRealWfCmd(s));
  ui->addCmd(new RandomizeVelCmd(s));
  ui->addCmd(new RunCmd(s));
  ui->addCmd(new MDSaveCmd(s));
  ui->addCmd(new SaveCmd(s));
  ui->addCmd(new SavesysCmd(s));
  ui->addCmd(new SavedenCmd(s));
  ui->addCmd(new SaveESPCmd(s));
  ui->addCmd(new SetCmd(s));
  ui->addCmd(new SpeciesCmd(s));
  ui->addCmd(new MMSpeciesCmd(s));
  ui->addCmd(new EmpiricalPotentialCmd(s));
  ui->addCmd(new StatusCmd(s));
  ui->addCmd(new KpointCmd(s));
  ui->addCmd(new SymmetryCmd(s));
  ui->addCmd(new ParOptCmd(s));
  ui->addCmd(new LockCmd(s));
  ui->addCmd(new UnlockCmd(s));
  ui->addCmd(new SetVelCmd(s));
  ui->addCmd(new ComputeMLWFCmd(s));
  ui->addCmd(new ConstraintCmd(s));
  ui->addCmd(new ShiftWFCmd(s));
  ui->addCmd(new WFPhaseRealCmd(s));
  ui->addCmd(new PlotCmd(s));
  ui->addCmd(new ResetVcmCmd(s));
  ui->addCmd(new CoordinatesCmd(s));

  ui->addVar(new AtomsDyn(s));
  ui->addVar(new Cell(s));
  ui->addVar(new CellDyn(s));
  ui->addVar(new CellLock(s));
  ui->addVar(new ChargeMixing(s));
  ui->addVar(new ChargeMixCoeff(s));
  ui->addVar(new ChargeMixRcut(s));
  ui->addVar(new ChargeMixNdim(s));
  ui->addVar(new CellMass(s));
  ui->addVar(new Debug(s));
  ui->addVar(new Ecut(s));
  ui->addVar(new Ecutprec(s));
  ui->addVar(new Ecutden(s));
  ui->addVar(new Ecuts(s));
  ui->addVar(new Emass(s));
  ui->addVar(new ExtStress(s));
  ui->addVar(new Smearing(s));
  ui->addVar(new SmearingWidth(s));
  ui->addVar(new FermiTemp(s));
  ui->addVar(new Dt(s));
  ui->addVar(new Nempty(s));
  ui->addVar(new Nrowmax(s));
  ui->addVar(new RefCell(s));
  ui->addVar(new Spin(s));
  ui->addVar(new Stress(s));
  ui->addVar(new Thermostat(s));
  ui->addVar(new ThresholdScf(s));
  ui->addVar(new ThresholdForce(s));
  ui->addVar(new ThresholdStress(s));
  ui->addVar(new ThTemp(s));
  ui->addVar(new ThTime(s));
  ui->addVar(new ThWidth(s));
  ui->addVar(new CenterOfMass(s));
  ui->addVar(new WfDiag(s));
  ui->addVar(new WfDyn(s));
  ui->addVar(new WfExtrap(s));
  ui->addVar(new Xc(s));
  ui->addVar(new Nparallelkpts(s));
  ui->addVar(new Nkpoints(s));
  ui->addVar(new IPrint(s));
  ui->addVar(new CellStepFreq(s));
  ui->addVar(new EnthalpyPressure(s));
  ui->addVar(new EnthalpyThreshold(s));
  ui->addVar(new HugoniostatVar(s));
  ui->addVar(new HugDeltaTemp(s));
  ui->addVar(new HugFreq(s));
  ui->addVar(new RunTimer(s));
  ui->addVar(new HubbardU(s));
  ui->addVar(new Memory(s));
  ui->addVar(new MatrixLoc(s));
  ui->addVar(new Pblock(s));
  ui->addVar(new MDIter(s));
  ui->addVar(new Force_Complex_WF(s));
  ui->addVar(new Non_Selfconsistent_Energy_Output(s));
  ui->addVar(new TDDt(s));
  ui->addVar(new NA_overlaps(s));
  ui->addVar(new WF_Phase_RealVar(s));
  ui->addVar(new SaveFreq(s));
  ui->addVar(new SaveDenFreq(s));
  ui->addVar(new SaveWfFreq(s));
  ui->addVar(new NetCharge(s));
  ui->addVar(new EsmBC(s));
  ui->addVar(new EsmW(s));
  ui->addVar(new FcpThermostat(s));
  ui->addVar(new FcpThTemp(s));
  ui->addVar(new FcpThTime(s));
  ui->addVar(new FcpThWidth(s));
  ui->addVar(new FcpPmass(s));
  ui->addVar(new FcpMu(s));
  ui->addVar(new VdW(s));
  
#ifdef USE_JAGGEMM
  setup_grid();
#endif  

  if ( argc == 2 )
  {
    // input file was given as a command line argument
    bool echo = true;
    ifstream in;
    if ( ctxt.oncoutpe() ) {
      in.open(argv[1], ios::in);
      if(!in.is_open()) {
	ui->error(string("Cannot open input file: ") + argv[1]);
	exit(1);
      }
    }
    ui->processCmds(in, "[qball]", echo, /*interactive =*/ false);
  }
  else
  {
    // use standard input
    bool echo = !isatty(0);
    ui->processCmds(cin, "[qball]", echo);
  }

  // exit using the quit command when a encountering EOF in a script
  Cmd *c = ui->findCmd("quit");
  c->action(1,NULL);

  delete ui;
  delete s;
  
  if ( ctxt.oncoutpe() )
  {
    cout << "<real_time> " << tm.real() << " </real_time>" << endl;
    cout << "<end_time> " << isodate() << " </end_time>" << endl;
    cout << "</qbox:simulation>" << endl;
  }

  } // end of Context scope
  TAU_PROFILE_STOP(timer);
#if USE_APC
  ApcFinalize();
#endif
#ifdef USE_OLD_CTF
  CTF_exit();
#endif
#if USE_MPI
  MPI_Finalize();
#endif

  return 0;
}
