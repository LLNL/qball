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
// qblink.C
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
using namespace std;

#include <sys/utsname.h>
#include <unistd.h>
#include <cstdlib>
#include <fstream>
#if AIX 
#include<filehdr.h>
#endif

#include "isodate.h"
#include "release.h"
#include "qbox_xmlns.h"

#include "Context.h"
#include "UserInterface.h"
#include "Sample.h"
#include "Timer.h"
#include "BOSampleStepper.h"
#include "CPSampleStepper.h"
#include "qbLink.h"

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

#ifdef USE_JAGGEMM
extern "C" int setup_grid();
#endif

qbLink::qbLink() {
  ctxt = new Context();
  active_ = true;
  // redirect cout
  qboxlog = new ofstream("qboxoutput.log");
  cout.flush();
  saved_cout = cout.rdbuf();
  cout_to_qboxlog();
  init();
  restore_cout();
}
qbLink::qbLink(string logfilename) {
  ctxt = new Context();
  active_ = true;
  // redirect cout
  qboxlog = new ofstream(logfilename.c_str());
  cout.flush();
  saved_cout = cout.rdbuf();
  cout_to_qboxlog();
  init();
  restore_cout();
}
qbLink::qbLink(string logfilename, int firstpe, int lastpe) {

  int npes;
#ifdef USE_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &npes);
#else
  npes = 1;
#endif
  int nrowmax = 32;

  int npes_sub = lastpe - firstpe + 1;
  if (npes_sub > npes) {
    npes_sub = npes;
    firstpe = 0;
    lastpe = npes-1;
  }

  int npr = nrowmax;
  while ( npes_sub%npr != 0 ) npr--;
  int npc = npes_sub/npr;
  int irstart = 0;
  int icstart = (int)(firstpe/npr);
  int npcfull = (int)(npes/npr);

  Context* fullctxt = new Context();
  active_ = false;
  onfirstpe_ = false;
  if (fullctxt->mype() >= firstpe && fullctxt->mype() <= lastpe){
    active_ = true; 
    if (fullctxt->mype() == firstpe) 
      onfirstpe_ = true;
  }

  // create context on processors we want to run on
  Context* tmpctxt = new Context((*fullctxt),1,npes_sub,0,firstpe);
  ctxt = 0;
  if (active_) {
    assert(tmpctxt->active());
    // copy context so process numbering is relative (i.e. mype = 0-npes, not firstpe-lastpe)
    ctxt = new Context(tmpctxt->comm(),1,npes_sub);
    ctxt->set_coutpe(0);
  }
  else {
    assert(!tmpctxt->active());
  }
  delete fullctxt;
  delete tmpctxt;

  /*
  //ewd:  this context creation caused an error in force calculation, since EnergyFunctional performs dwf.dsum to compute averaged forces
  //ewd:  and qbLink s_.wf.ctxt is nprow x npcol, while qb s_.wf.ctxt is single-row.   
  Context* fullctxt = new Context(npr,npcfull);
  active_ = false;

  onfirstpe_ = false;
  if (fullctxt->mype() >= firstpe && fullctxt->mype() <= lastpe){
    active_ = true; 
    if (fullctxt->mype() == firstpe) 
      onfirstpe_ = true;
  }

  // create context on processors we want to run on
  Context* tmpctxt = new Context((*fullctxt),npr,npc,irstart,icstart);
  ctxt = 0;
  if (active_) {
    assert(tmpctxt->active());
    // copy context so process numbering is relative (i.e. mype = 0-npes, not firstpe-lastpe)
    ctxt = new Context(tmpctxt->comm(),npr,npc);
    ctxt->set_coutpe(0);
  }
  else {
    assert(!tmpctxt->active());
  }
  delete fullctxt;
  delete tmpctxt;
  */
  
  saved_cout = cout.rdbuf();
  if (active_) {
    // redirect cout
    qboxlog = new ofstream(logfilename.c_str());
    cout.flush();
    cout_to_qboxlog();
    init();
    restore_cout();
  }
  else {
    //cout.rdbuf(empty_cout);
  }
}
qbLink::~qbLink() {
  if (active_) {
    cout_to_qboxlog();
    if (stepper != 0) 
      delete stepper;
    delete ui;
    delete s;
    if (ctxt->oncoutpe())
      cout << "</qbox:simulation>" << endl;
    restore_cout();
  }
  if (active_) {
    qboxlog->close();
    delete qboxlog;
  }
  if (ctxt != 0)
    delete ctxt;
}
void qbLink::init(void) {
  s = new Sample(*ctxt);
  ui = new UserInterface(*ctxt);
  ui->set_coutpe(0);
  stepper = 0;
  if ( ctxt->oncoutpe() ) {
    cout << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
    cout << "<!--\n\n";
    cout << "                   ===========================\n";
    cout << "                   I qbox " 
       << setw(17) << left << release() << "  I\n";
    cout << "                   I  link interface         I\n";
    cout << "                   I                         I\n";
    cout << "                   I                         I\n";
    cout << "                   I                         I\n";
    cout << "                   I                         I\n";
    cout << "                   I                         I\n";
    cout << "                   I                         I\n";
    cout << "                   I                         I\n";
    cout << "                   I                         I\n";
    cout << "                   I                         I\n";
    cout << "                   I                         I\n";
    cout << "                   I  F.Gygi and E. Draeger  I\n";
    cout << "                   I                  LLNL   I\n";
    cout << "                   I Copyright (c) 2003-2008 I\n";     
    cout << "                   ===========================\n\n";
    cout << "-->\n";
    cout << "<qbox:simulation xmlns:qbox=\"" << qbox_xmlns() << "\">" << endl;
    cout << "<npes> " << ctxt->size() << " </npes>" << endl;
  }
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
  ui->addCmd(new RandomizeVelCmd(s));
  ui->addCmd(new RunCmd(s));
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
  ui->addCmd(new ComputeMLWFCmd(s));
  ui->addCmd(new ConstraintCmd(s));
  ui->addCmd(new ShiftWFCmd(s));
  ui->addCmd(new WFPhaseRealCmd(s));
  ui->addCmd(new PlotCmd(s));
  ui->addCmd(new ResetVcmCmd(s));
  
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

#ifdef USE_JAGGEMM
  setup_grid();
#endif  

}
void qbLink::cout_to_qboxlog(void) {
  cout.flush();
  if (ctxt != 0) 
    if (ctxt->oncoutpe()) 
      cout.rdbuf(qboxlog->rdbuf());

 //ewd: original version:
  /*
  if (ctxt != 0) {
    if (ctxt->oncoutpe()) 
      cout.rdbuf(qboxlog->rdbuf());
    //else 
      //cout.rdbuf(empty_cout);
  }
  //else 
    //cout.rdbuf(empty_cout);
  */
}
void qbLink::restore_cout(void) {
  cout.flush();
  cout.rdbuf(saved_cout);
}
void qbLink::processInputLine(string inputline) {
  if (active_) {
    cout_to_qboxlog();
    stringstream is; 
    is << inputline << endl;
    ui->processCmds(is,"[qbLink]",true);
    restore_cout();
  }
}
void qbLink::processInputFile(string filename) {
  if (active_) {
    cout_to_qboxlog();
    ifstream infile;
    if ( ctxt->oncoutpe() )
      infile.open(filename.c_str(),ios::in);
    ui->processCmds(infile, "[qbLink]", true);

    restore_cout();
  }
}
void qbLink::set_ecut(double ecut) {
  if (active_) {
    cout_to_qboxlog();
    const int nargs = 1;
    Var*  ecutvar = new Ecut(s);
    ecutvar->ui = ui;
    char** args = new char*[nargs+1];

    args[0] = (char*)"ecut";

    double a[nargs] = {ecut};
    int sflag;
    char tmpbuf[20][nargs];
    for (int i=0; i<nargs; i++) {
      sflag = sprintf(tmpbuf[i],"%0.9f",a[i]);
      args[i+1] = tmpbuf[i];
    }

    ecutvar->set(nargs+1,args);
    delete ecutvar;
    delete args;
    restore_cout();
  }
}
void qbLink::set_cell(double a0x, double a0y, double a0z, double a1x, double a1y, double a1z, double a2x, double a2y, double a2z) {

  if (active_) {
    cout_to_qboxlog();
    const int nargs = 9;
    Var*  cellvar = new Cell(s);
    cellvar->ui = ui;
    char** args = new char*[nargs+1];
    
    args[0] = (char*)"cell";

    double a[nargs] = {a0x,a0y,a0z,a1x,a1y,a1z,a2x,a2y,a2z};
    int sflag;
    char tmpbuf[20][nargs];
    for (int i=0; i<nargs; i++) {
      sflag = sprintf(tmpbuf[i],"%0.9f",a[i]);
      args[i+1] = tmpbuf[i];
    }

    cellvar->set(nargs+1,args);
    delete cellvar;
    delete args;
    restore_cout();
  }
}
void qbLink::get_cell(double &a0x, double &a0y, double &a0z, double &a1x, double &a1y, double &a1z, double &a2x, double &a2y, double &a2z) {
  
  D3vector a0 = s->wf.cell().a(0);
  D3vector a1 = s->wf.cell().a(1);
  D3vector a2 = s->wf.cell().a(2);

  a0x = a0.x;
  a0y = a0.y;
  a0z = a0.z;
  a1x = a1.x;
  a1y = a1.y;
  a1z = a1.z;
  a2x = a2.x;
  a2y = a2.y;
  a2z = a2.z;
}
void qbLink::set_stress(string onoff) {
  if (active_) {
    cout_to_qboxlog();
    const int nargs = 1;
    Var*  stressvar = new Stress(s);
    stressvar->ui = ui;
    char** args = new char*[nargs+1];

    args[0] = (char*)"stress";
    args[1] = (char*)onoff.c_str();

    stressvar->set(nargs+1,args);
    delete stressvar;
    delete args;
  }
}
void qbLink::set_th_time(double time) {
  if (active_) {
    cout_to_qboxlog();
    const int nargs = 1;
    Var*  thtimevar = new ThTime(s);
    thtimevar->ui = ui;
    char** args = new char*[nargs+1];

    args[0] = (char*)"th_time";
    int sflag;
    char tmpbuf[20];
    sflag = sprintf(tmpbuf,"%0.9f",time);
    args[1] = tmpbuf;
    thtimevar->set(nargs+1,args);
    delete thtimevar;
    delete args;
  }
}
void qbLink::set_th_temp(double temp) {
  if (active_) {
    cout_to_qboxlog();
    const int nargs = 1;
    Var*  thtempvar = new ThTemp(s);
    thtempvar->ui = ui;
    char** args = new char*[nargs+1];

    args[0] = (char*)"th_temp";
    int sflag;
    char tmpbuf[20];
    sflag = sprintf(tmpbuf,"%0.9f",temp);
    args[1] = tmpbuf;
    thtempvar->set(nargs+1,args);
    delete thtempvar;
    delete args;
  }
}
void qbLink::set_dt(double dt) {
  if (active_) {
    cout_to_qboxlog();
    const int nargs = 1;
    Var*  dtvar = new Dt(s);
    dtvar->ui = ui;
    char** args = new char*[nargs+1];

    args[0] = (char*)"dt";
    int sflag;
    char tmpbuf[20];
    sflag = sprintf(tmpbuf,"%0.9f",dt);
    args[1] = tmpbuf;
    dtvar->set(nargs+1,args);
    delete dtvar;
    delete args;
  }
}
void qbLink::set_fermi_temp(double temp) {
  if (active_) {
    cout_to_qboxlog();
    const int nargs = 1;
    Var*  fermitempvar = new FermiTemp(s);
    fermitempvar->ui = ui;
    char** args = new char*[nargs+1];

    args[0] = (char*)"fermi_temp";
    int sflag;
    char tmpbuf[20];
    sflag = sprintf(tmpbuf,"%0.9f",temp);
    args[1] = tmpbuf;
    fermitempvar->set(nargs+1,args);
    delete fermitempvar;
    delete args;
  }
}
void qbLink::runBOSteps(int niter, int nitscf, int nite) {
  if (active_) {
    cout_to_qboxlog();
    if (stepper != 0)
      delete stepper;

    stepper = new BOSampleStepper(*s,nitscf,nite);
    assert(stepper!=0);
    if (ctxt->oncoutpe() ) 
      cout << "<run niter_ionic=\"" << niter << "\" niter_scf=\"" << nitscf << "\" niter_nonscf=\"" << nite << "\">" << endl;
    s->wf.info(cout,"wavefunction");
    stepper->step(niter);
    if (ctxt->oncoutpe() ) 
      cout << "</run>" << endl;

    restore_cout();
  }
}
void qbLink::runCPSteps(int niter) {
  if (active_) {
    cout_to_qboxlog();
    if (stepper != 0)
      delete stepper;

    stepper = new CPSampleStepper(*s);
    assert(stepper!=0);
    if (ctxt->oncoutpe() ) 
      cout << "<run niter_ionic=\"" << niter << "\">" << endl;
    s->wf.info(cout,"wavefunction");
    stepper->step(niter);
    if (ctxt->oncoutpe() ) 
      cout << "</run>" << endl;
    restore_cout();
  }
}
int qbLink::nsp(void) {
  if (active_) 
    return s->atoms.nsp();
  else 
    return 0;
}
int qbLink::nsp_mm(void) {
  if (active_) 
    return s->atoms.nsp_mm();
  else 
    return 0;
}
int qbLink::na(int isp) {
  if (active_) 
    return s->atoms.na(isp);
  else 
    return 0;
}
int qbLink::na_mm(int isp) {
  if (active_) 
    return s->atoms.na_mm(isp);
  else 
    return 0;
}
double qbLink::mass(int isp) {
  if (active_) 
    return s->atoms.mass(isp);
  else 
    return 0.0;
}
double qbLink::mass_mm(int isp) {
  if (active_) 
    return s->atoms.mass_mm(isp);
  else 
    return 0.0;
}
int qbLink::atomic_number(int isp) {
  if (active_) 
    return s->atoms.atomic_number(isp);
  else 
    return 0;
}
void qbLink::get_positions(vector<vector<double> > &r) {
  if (active_) 
    s->atoms.get_positions(r,false);
}
void qbLink::set_positions(vector<vector<double> > &r) {
  if (active_) 
    s->atoms.set_positions(r,true);
}
void qbLink::get_velocities(vector<vector<double> > &v) {
  if (active_) 
    s->atoms.get_velocities(v,false);
}
void qbLink::set_velocities(vector<vector<double> > &v) {
  if (active_) 
    s->atoms.set_velocities(v);
}
void qbLink::get_forces(vector<vector<double> > &f) {
  if (active_) {
    assert(stepper != 0);
    stepper->get_forces(f);
  }
}
void qbLink::get_fion_ext(vector<vector<double> > &f) {
  if (active_) 
    s->atoms.get_fion_ext(f);
}
void qbLink::set_fion_ext(vector<vector<double> > &f) {
  if (active_) 
    s->atoms.set_fion_ext(f);
}
double qbLink::get_etotal(void) { 
  if (stepper != 0)
    return stepper->get_energy("etotal"); 
  else 
    return 0.0;
}
double qbLink::get_ekin(void) { 
  if (stepper != 0)
    return stepper->get_energy("ekin"); 
  else 
    return 0.0;
}
double qbLink::get_econf(void) { 
  if (stepper != 0)
    return stepper->get_energy("econf"); 
  else 
    return 0.0;
}
double qbLink::get_eps(void) { 
  if (stepper != 0)
    return stepper->get_energy("eps"); 
  else 
    return 0.0;
}
double qbLink::get_enl(void) { 
  if (stepper != 0)
    return stepper->get_energy("enl"); 
  else 
    return 0.0;
}
double qbLink::get_ehart(void) { 
  if (stepper != 0)
    return stepper->get_energy("ehart"); 
  else 
    return 0.0;
}
double qbLink::get_ecoul(void) { 
  if (stepper != 0)
    return stepper->get_energy("ecoul"); 
  else 
    return 0.0;
}
double qbLink::get_exc(void) { 
  if (stepper != 0)
    return stepper->get_energy("exc"); 
  else 
    return 0.0;
}
double qbLink::get_esr(void) { 
  if (stepper != 0)
    return stepper->get_energy("esr"); 
  else 
    return 0.0;
}
double qbLink::get_eself(void) { 
  if (stepper != 0)
    return stepper->get_energy("eself"); 
  else 
    return 0.0;
}
double qbLink::get_ets(void) { 
  if (stepper != 0)
    return stepper->get_energy("ets"); 
  else 
    return 0.0;
}
valarray<double> qbLink::get_stress_tot(void) { 
  const double gpa = 29421.0120;
  if (stepper != 0)
    return gpa*stepper->get_stress("total"); 
  else {
    valarray<double> nullv(6);
    for (int i=0; i<6; i++)
      nullv[i] = 0.0;
    return nullv; 
  }
}
valarray<double> qbLink::get_stress_kin(void) { 
  const double gpa = 29421.0120;
  if (stepper != 0)
    return gpa*stepper->get_stress("kin"); 
  else {
    valarray<double> nullv(6);
    for (int i=0; i<6; i++)
      nullv[i] = 0.0;
    return nullv; 
  }
}
valarray<double> qbLink::get_stress_ext(void) { 
  const double gpa = 29421.0120;
  if (stepper != 0)
    return gpa*stepper->get_stress("ext"); 
  else {
    valarray<double> nullv(6);
    for (int i=0; i<6; i++)
      nullv[i] = 0.0;
    return nullv; 
  }
}
valarray<double> qbLink::get_stress_eks(void) { 
  const double gpa = 29421.0120;
  if (stepper != 0)
    return gpa*stepper->get_stress("eks"); 
  else {
    valarray<double> nullv(6);
    for (int i=0; i<6; i++)
      nullv[i] = 0.0;
    return nullv; 
  }
}
