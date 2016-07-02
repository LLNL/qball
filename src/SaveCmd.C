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
// SaveCmd.C:
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>


#include "SaveCmd.h"
#include "SampleWriter.h"
#include "fstream"
#include "isodate.h"
#include "release.h"
#include "qbox_xmlns.h"
#include "Timer.h"
#include "Context.h"
#include "ChargeDensity.h"
#include "FourierTransform.h"
#include "Wavefunction.h"
#include "EnergyFunctional.h"
#include "AtomSet.h"
#include "release.h"

#ifdef USE_CSTDIO_LFS
#include <cstdio>
#include <cstdlib>
#include <sstream>
#endif
using namespace std;

////////////////////////////////////////////////////////////////////////////////
int SaveCmd::action(int argc, char **argv) {

  if ( !(argc>=2 && argc<=4 ) ) {
    if ( ui->oncoutpe() ) {
      cout << "  <!-- use: save [encoding options] [state format options] filename   -->" << endl;
      cout << "  <!--     encoding options:                                          -->" << endl;
      cout << "  <!--       -dump:  one binary file for each process                 -->" << endl;
      cout << "  <!--       -fast:  one binary file for each I/O node                 -->" << endl;
      cout << "  <!--       -states:  one file for each state                        -->" << endl;
      cout << "  <!--       -xml:  entire system in one xml file (Gamma-point only)  -->" << endl;
      cout << "  <!--       -casino:  wavefunction in CASINO trial function format  -->" << endl;
      cout << "  <!--       -vmd:  wavefunction and density in VMD cube visualization format  -->" << endl;

      cout << "  <!--     state format options:                                      -->" << endl;
      cout << "  <!--       -binary:  save each state as binary file (default)       -->" << endl;
      cout << "  <!--       -molmol:  Molmol isosurface format for visualization     -->" << endl;
      cout << "  <!--       -gopenmol:  GOpenMol isosurface format for visualization -->" << endl;
    } 
    return 1;
  }
  
  // set default encoding and format
  string encoding = "dump";
  string format = "binary";
  char* filename = 0;
  bool atomsonly = false;
  bool writevel = false;
  bool serial = false; 
  bool base64 = true;

  // parse arguments
  for ( int i = 1; i < argc; i++ ) {
    string arg(argv[i]);
    
    if ( arg=="-dump" )
      encoding = "dump";
    else if ( arg=="-fast" )
      encoding = "fast";
    else if ( arg=="-states" )
      encoding = "states";
    else if ( arg=="-xml" )
      encoding = "xml";
    else if ( arg=="-casino" )
      encoding = "casino";
    else if ( arg=="-vel" )
      writevel = true;
    else if ( arg=="-binary" )
      format = "binary";
    else if ( arg=="-molmol" ) {
      encoding = "states";
      format = "molmol";
    }
    else if ( arg=="-vmd" ) {
      encoding = "vmd";
    }
    else if ( arg=="-text" ) {
      encoding = "states";
      format = "text";
    }
    else if ( arg=="-gopenmol" ) {
      encoding = "states";
      format = "gopenmol";
    }
    else if ( arg=="-atomsonly" )
      atomsonly = true;
    else if ( arg[0] != '-' && i == argc-1 )
      filename = argv[i];
    else {
      if ( ui->oncoutpe() ) {
        cout << "  <!-- use: save [encoding options] [state format options] filename   -->" << endl;
        cout << "  <!--     encoding options:                                          -->" << endl;
        cout << "  <!--       -dump:  one binary file for each process                 -->" << endl;
        cout << "  <!--       -states:  one file for each state                        -->" << endl;
        cout << "  <!--       -xml:  entire system in one xml file (Gamma-point only)  -->" << endl;
        cout << "  <!--       -vmd:  wavefunction and density in VMD cube visualization format  -->" << endl;
        cout << "  <!--     state format options:                                      -->" << endl;
        cout << "  <!--       -binary:  save each state as binary file (default)       -->" << endl;
        cout << "  <!--       -molmol:  Molmol isosurface format for visualization     -->" << endl;
        cout << "  <!--       -gopenmol:  GOpenMol isosurface format for visualization -->" << endl;
      } 
      return 1;
    }
  }
  
  if ( filename == 0 ) {
    if ( ui->oncoutpe() ) {
      cout << "  <!-- use: save [encoding options] [state format options] filename   -->" << endl;
      cout << "  <!--     encoding options:                                          -->" << endl;
      cout << "  <!--       -dump:  one binary file for each process                 -->" << endl;
      cout << "  <!--       -states:  one file for each state                        -->" << endl;
      cout << "  <!--       -xml:  entire system in one xml file (Gamma-point only)  -->" << endl;
      cout << "  <!--       -vmd:  wavefunction and density in VMD cube visualization format  -->" << endl;
      cout << "  <!--     state format options:                                      -->" << endl;
      cout << "  <!--       -binary:  save each state as binary file (default)       -->" << endl;
      cout << "  <!--       -molmol:  Molmol isosurface format for visualization     -->" << endl;
      cout << "  <!--       -gopenmol:  GOpenMol isosurface format for visualization -->" << endl;
    }
    return 1;
  }

  if (s->ctrl.timer_hit) {
    if (s->ctrl.timer_savecmd) {
      if ( ui->oncoutpe() )
        cout << " <!-- SaveCmd: run_timer exceeded: all save commands beyond first will be ignored. -->" << endl;
      return 0;
    }
    else
      s->ctrl.timer_savecmd = true;
  }
  
  Timer savetm;
  savetm.start();

  string filestr(filename);

  // binary output
  if (encoding == "dump" ) {
     if ( ui->oncoutpe() )
        cout << "<!-- SaveCmd:  writing wf " << filestr << "... -->" << endl;
     s->wf.write_dump(filestr);
     s->wf.write_mditer(filestr,s->ctrl.mditer);
     if (s->ctrl.tddft_involved)
     {
       // write s->hamil_wf
       string hamwffile = filestr + "hamwf";
       if ( ui->oncoutpe() )
          cout << "<!-- SaveCmd:  wf write finished, writing hamil_wf " << hamwffile << "... -->" << endl;
       s->hamil_wf->write_dump(hamwffile);
    }
    else
    {
       // ewd:  write rhor_last to file
       ChargeDensity cd_(*s);
       if (s->ctrl.ultrasoft)    
          cd_.update_usfns();

       const int nspin = s->wf.nspin();
       // load mixed charge density into cd_
       for (int ispin = 0; ispin < nspin; ispin++) {
          for ( int i=0; i < s->rhog_last[ispin].size(); i++ )
             cd_.rhog[ispin][i] = s->rhog_last[ispin][i];
       }
       cd_.update_rhor();
       
       const Context* wfctxt = s->wf.wfcontext();
       const Context* vctxt = &cd_.vcontext();
       FourierTransform* ft_ = cd_.vft();
       for (int ispin = 0; ispin < nspin; ispin++) {
          if (wfctxt->mycol() == 0) {
             vector<double> rhortmp(ft_->np012loc());
             for (int j = 0; j < ft_->np012loc(); j++)
                rhortmp[j] = cd_.rhor[ispin][j];

             ofstream os;
             string rhorfile;
             if (nspin == 1)
                rhorfile = filestr + ".lastrhor";
             else {
                ostringstream oss;
                oss.width(1);  oss.fill('0');  oss << ispin;
                rhorfile = filestr + ".s" + oss.str() + ".lastrhor";
             }
             
             if (wfctxt->onpe0()) {
                os.open(rhorfile.c_str(),ofstream::binary);
                // hack to make checkpointing work w. BlueGene compilers
#ifdef BGQ
                os.write(rhorfile.c_str(),sizeof(char)*rhorfile.length());
#endif
                
             }
             for ( int i = 0; i < wfctxt->nprow(); i++ ) {
                if ( i == wfctxt->myrow() ) {
                   int size = ft_->np012loc();
                   wfctxt->isend(1,1,&size,1,0,0);
                   if (size > 0)
                      wfctxt->dsend(size,1,&rhortmp[0],1,0,0);
                }
             }
             if ( wfctxt->onpe0() ) {
                for ( int i = 0; i < wfctxt->nprow(); i++ ) {
                   int size = 0;
                   wfctxt->irecv(1,1,&size,1,i,0);
                   if (size > 0)
                      wfctxt->drecv(size,1,&rhortmp[0],1,i,0);
                   os.write((char*)&rhortmp[0],sizeof(double)*size);
                }
                os.close();
             }
          }
       }
    }
    
    if (writevel) { 
      const string atoms_dyn = s->ctrl.atoms_dyn;
      const bool compute_forces = ( atoms_dyn != "LOCKED" );
      if (compute_forces) {
        string wfvfile = filestr + "wfv";
        s->wfv->write_dump(wfvfile);
      }
      else {
        if ( ui->oncoutpe() )
          cout << "<WARNING>SaveCmd:  wavefunction velocities only available when atoms_dyn set, not saving.</WARNING>" << endl;
      }
    }

    if (atomsonly)
      if ( ui->oncoutpe() )
        cout << "<!-- SaveCmd:  atomsonly flag only used with xml output, ignoring. -->" << endl;

    if (format != "binary" )
      if ( ui->oncoutpe() )
        cout << "<!-- SaveCmd:  " << format << " flag only used with -states output, ignoring. -->" << endl;  }

  // binary output
  else if (encoding == "fast" ) {
     if ( ui->oncoutpe() )
        cout << "<!-- SaveCmd:  writing wf " << filestr << "... -->" << endl;
     s->wf.write_fast(filestr);
     s->wf.write_mditer(filestr,s->ctrl.mditer);
     if (s->ctrl.tddft_involved)
     {
       // write s->hamil_wf
       string hamwffile = filestr + "hamwf";
       if ( ui->oncoutpe() )
          cout << "<!-- SaveCmd:  wf write finished, writing hamil_wf " << hamwffile << "... -->" << endl;
       s->hamil_wf->write_fast(hamwffile);
    }
    else
    {
       // ewd:  write rhor_last to file
       ChargeDensity cd_(*s);
       if (s->ctrl.ultrasoft)    
          cd_.update_usfns();

       const int nspin = s->wf.nspin();
       // load mixed charge density into cd_
       for (int ispin = 0; ispin < nspin; ispin++) {
          for ( int i=0; i < s->rhog_last[ispin].size(); i++ )
             cd_.rhog[ispin][i] = s->rhog_last[ispin][i];
       }
       cd_.update_rhor();
       
       const Context* wfctxt = s->wf.wfcontext();
       const Context* vctxt = &cd_.vcontext();
       FourierTransform* ft_ = cd_.vft();
       for (int ispin = 0; ispin < nspin; ispin++) {
          if (wfctxt->mycol() == 0) {
             vector<double> rhortmp(ft_->np012loc());
             for (int j = 0; j < ft_->np012loc(); j++)
                rhortmp[j] = cd_.rhor[ispin][j];

             ofstream os;
             string rhorfile;
             if (nspin == 1)
                rhorfile = filestr + ".lastrhor";
             else {
                ostringstream oss;
                oss.width(1);  oss.fill('0');  oss << ispin;
                rhorfile = filestr + ".s" + oss.str() + ".lastrhor";
             }
             
             if (wfctxt->onpe0()) {
                os.open(rhorfile.c_str(),ofstream::binary);
                // hack to make checkpointing work w. BlueGene compilers
#ifdef BGQ
                os.write(rhorfile.c_str(),sizeof(char)*rhorfile.length());
#endif
                
             }
             for ( int i = 0; i < wfctxt->nprow(); i++ ) {
                if ( i == wfctxt->myrow() ) {
                   int size = ft_->np012loc();
                   wfctxt->isend(1,1,&size,1,0,0);
                   if (size > 0)
                      wfctxt->dsend(size,1,&rhortmp[0],1,0,0);
                }
             }
             if ( wfctxt->onpe0() ) {
                for ( int i = 0; i < wfctxt->nprow(); i++ ) {
                   int size = 0;
                   wfctxt->irecv(1,1,&size,1,i,0);
                   if (size > 0)
                      wfctxt->drecv(size,1,&rhortmp[0],1,i,0);
                   os.write((char*)&rhortmp[0],sizeof(double)*size);
                }
                os.close();
             }
          }
       }
    }
    
    if (writevel) { 
      const string atoms_dyn = s->ctrl.atoms_dyn;
      const bool compute_forces = ( atoms_dyn != "LOCKED" );
      if (compute_forces) {
        string wfvfile = filestr + "wfv";
        s->wfv->write_fast(wfvfile);
      }
      else {
        if ( ui->oncoutpe() )
          cout << "<WARNING>SaveCmd:  wavefunction velocities only available when atoms_dyn set, not saving.</WARNING>" << endl;
      }
    }

    if (atomsonly)
      if ( ui->oncoutpe() )
        cout << "<!-- SaveCmd:  atomsonly flag only used with xml output, ignoring. -->" << endl;

    if (format != "binary" )
      if ( ui->oncoutpe() )
        cout << "<!-- SaveCmd:  " << format << " flag only used with -states output, ignoring. -->" << endl;
    
  }

  else if (encoding == "states" ) {
     if ( ui->oncoutpe() )
        cout << "<!-- SaveCmd:  writing wf " << filestr << "... -->" << endl;
     s->wf.write_states(filestr,format);
     s->wf.write_mditer(filestr,s->ctrl.mditer);

    if (s->ctrl.tddft_involved)
    {
       // write s->hamil_wf
       string hamwffile = filestr + "hamwf";
       if ( ui->oncoutpe() )
          cout << "<!-- SaveCmd:  wf write finished, writing hamil_wf " << hamwffile << "... -->" << endl;
       s->hamil_wf->write_states(hamwffile,format);
    }
    else
    {
       ChargeDensity cd_(*s);
       if (s->ctrl.ultrasoft)    
          cd_.update_usfns();

       const int nspin = s->wf.nspin();
       // load mixed charge density into cd_
       for (int ispin = 0; ispin < nspin; ispin++) {
          for ( int i=0; i < s->rhog_last[ispin].size(); i++ )
             cd_.rhog[ispin][i] = s->rhog_last[ispin][i];
       }
       cd_.update_rhor();

       const Context* wfctxt = s->wf.wfcontext();
       const Context* vctxt = &cd_.vcontext();
       FourierTransform* ft_ = cd_.vft();
       for (int ispin = 0; ispin < nspin; ispin++) {
          if (wfctxt->mycol() == 0) {
             vector<double> rhortmp(ft_->np012loc());
             for (int j = 0; j < ft_->np012loc(); j++)
                rhortmp[j] = cd_.rhor[ispin][j];

             ofstream os;
             string rhorfile;
             if (nspin == 1)
                rhorfile = filestr + ".lastrhor";
             else {
                ostringstream oss;
                oss.width(1);  oss.fill('0');  oss << ispin;
                rhorfile = filestr + ".s" + oss.str() + ".lastrhor";
             }
             
             if (wfctxt->onpe0()) {
                os.open(rhorfile.c_str(),ofstream::binary);
                // hack to make checkpointing work w. BlueGene compilers
#ifdef BGQ
                os.write(rhorfile.c_str(),sizeof(char)*rhorfile.length());
#endif          
             }
             for ( int i = 0; i < wfctxt->nprow(); i++ ) {
                if ( i == wfctxt->myrow() ) {
                   int size = ft_->np012loc();
                   wfctxt->isend(1,1,&size,1,0,0);
                   if (size > 0)
                      wfctxt->dsend(size,1,&rhortmp[0],1,0,0);
                }
             }
             if ( wfctxt->onpe0() ) {
                for ( int i = 0; i < wfctxt->nprow(); i++ ) {
                   int size = 0;
                   wfctxt->irecv(1,1,&size,1,i,0);
                   if (size > 0)
                      wfctxt->drecv(size,1,&rhortmp[0],1,i,0);
                   os.write((char*)&rhortmp[0],sizeof(double)*size);
                }
                os.close();
             }
          }
       }
    }
    
    if (writevel) { 
      const string atoms_dyn = s->ctrl.atoms_dyn;
      const bool compute_forces = ( atoms_dyn != "LOCKED" );
      if (compute_forces) {
        string wfvfile = filestr + "wfv";
        s->wfv->write_states(wfvfile,format);
      }
      else {
        if ( ui->oncoutpe() )
          cout << "<WARNING>SaveCmd:  wavefunction velocities only available when atoms_dyn set, not saving.</WARNING>" << endl;
      }
    }

    if (atomsonly)
      if ( ui->oncoutpe() )
        cout << "<!-- SaveCmd:  atomsonly flag only used with xml output, ignoring. -->" << endl;
  }

  // xml output
  else if (encoding == "xml" ) {
    SampleWriter swriter(s->ctxt_);
    string description = string(" Created ") + isodate() +
      string(" by qbox-") + release() + string(" ");
    swriter.writeSample(*s, filename, description, base64, atomsonly, serial);
  }
  else if (encoding == "casino" ) {
    // need to recalculate energy from wavefunction
    ChargeDensity cd_(*s);
    cd_.update_density();
    EnergyFunctional ef_(*s,s->wf,cd_);
    Wavefunction dwf(s->wf);
    ef_.update_vhxc();
    const bool compute_forces = false;
    const bool compute_stress = false;
    valarray<double> sigma_eks(6);
    vector<vector<double> > fion;
    fion.resize(s->atoms.nsp()+s->atoms.nsp_mm());
    for ( int is = 0; is < fion.size(); is++ )
      fion[is].resize(3*s->atoms.na(is));

    ef_.energy(false,dwf,false,fion,false,sigma_eks);
    double eewald = ef_.casino_ewald();
    double evloc = ef_.casino_vloc();

    for (int ispin=0; ispin<s->wf.nspin(); ispin++)
    {
       if (s->wf.spinactive(ispin))
       {
          for (int kk = 0; kk<s->wf.nkp(); kk++)
          {
             if (s->wf.kptactive(kk))
             {

                // write out separate casino file for each k-point
                ostringstream oss1,oss2;
                oss1.width(4);  oss1.fill('0');  oss1 << kk;
                oss2.width(1);  oss2.fill('0');  oss2 << ispin;
                string casinofile = filestr + ".spin" + oss2.str() + "kpt" + oss1.str();

                int kpc = s->wf.mysdctxt(kk);
                Context* sdctxt = 0;
                if (kpc >=0 )
                   sdctxt = (Context*)s->wf.sdcontext(ispin,kpc);
          
                ofstream os;
                if (sdctxt != 0)
                {
                   if ( sdctxt->active() )
                   {
                      if ( sdctxt->myproc() == 0 )
                      {
                         os.open(casinofile.c_str(),ofstream::out);
                         os.setf(ios::fixed,ios::floatfield);
                         os << setprecision(10);
                         os << "k-point:  " << s->wf.kpoint(kk) << "   weight:  " << s->wf.weight(kk) << endl;
                         os << "" << endl;
                         os << "BASIC INFO" << endl;
                         os << "----------" << endl;
                         os << "Generated by:" << endl;
                         os << "Qbox " << release() << endl;
                         os << "Method:" << endl;
                         os << "DFT" << endl;
                         os << "DFT Functional" << endl;
                         os << "unknown" << endl;
                         os << "Pseudopotential" << endl;
                         os << "unknown" << endl;
                         os << "Plane wave cutoff (au)" << endl;
                         os << s->wf.ecut() << endl;
                         os << "Spin polarized:" << endl;
                         if (s->wf.nspin() > 1)
                            os << "T" << endl;
                         else
                            os << "F" << endl;
                         os << "Total energy (au per primitive cell)" << endl;
                         os << ef_.etotal() << endl;
                         os << "Kinetic energy (au per primitive cell)" << endl;
                         os << ef_.ekin() << endl;
                         os << "Local potential energy (au per primitive cell)" << endl;
                         os << ef_.eps() << endl;
                         os << "Non-local potential energy (au per primitive cell)" << endl;
                         os << ef_.enl() << endl;
                         os << "Electron-electron energy (au per primitive cell)" << endl;
                         os << ef_.ehart() << endl;
                         os << "Ion-ion energy (au per primitive cell)" << endl;
                         //os << ef_.esr()-ef_.eself() << endl;
                         os << eewald << endl;
                         os << "Number of electrons per primitive cell" << endl;
                         os << s->wf.nel() << endl;
                         os << endl;
                         s->atoms.print_casino(os);
                      }
                      sdctxt->barrier();

                      s->wf.print_casino(os,ispin,kk);

                      if ( sdctxt->myproc() == 0 ) 
                         os.close();
                      sdctxt->barrier();
                   }
                }
             }
          }
       }
    }
    
    /* old method:  write out each k-point one at a time
    for (int ispin=0; ispin<s->wf.nspin(); ispin++)
    {
       for (int kk = 0; kk<s->wf.nkp(); kk++)
       {
          // write out separate casino file for each k-point
          ostringstream oss1,oss2;
          oss1.width(4);  oss1.fill('0');  oss1 << kk;
          oss2.width(1);  oss2.fill('0');  oss2 << ispin;
          string casinofile = filestr + ".spin" + oss2.str() + "kpt" + oss1.str();

          int kpc = s->wf.mysdctxt(kk);
          Context* sdctxt = 0;
          if (kpc >=0 )
             sdctxt = (Context*)s->wf.sdcontext(ispin,kpc);
          
          ofstream os;
          if (sdctxt != 0)
          {
             if ( sdctxt->active() )
             {
                if ( sdctxt->myproc() == 0 )
                {
                   os.open(casinofile.c_str(),ofstream::out);
                   os.setf(ios::fixed,ios::floatfield);
                   os << setprecision(10);
                   os << "k-point:  " << s->wf.kpoint(kk) << "   weight:  " << s->wf.weight(kk) << endl;
                   os << "" << endl;
                   os << "BASIC INFO" << endl;
                   os << "----------" << endl;
                   os << "Generated by:" << endl;
                   os << "Qbox " << release() << endl;
                   os << "Method:" << endl;
                   os << "DFT" << endl;
                   os << "DFT Functional" << endl;
                   os << "unknown" << endl;
                   os << "Pseudopotential" << endl;
                   os << "unknown" << endl;
                   os << "Plane wave cutoff (au)" << endl;
                   os << s->wf.ecut() << endl;
                   os << "Spin polarized:" << endl;
                   if (s->wf.nspin() > 1)
                      os << "T" << endl;
                   else
                      os << "F" << endl;
                   os << "Total energy (au per primitive cell)" << endl;
                   os << ef_.etotal() << endl;
                   os << "Kinetic energy (au per primitive cell)" << endl;
                   os << ef_.ekin() << endl;
                   os << "Local potential energy (au per primitive cell)" << endl;
                   os << ef_.eps() << endl;
                   os << "Non-local potential energy (au per primitive cell)" << endl;
                   os << ef_.enl() << endl;
                   os << "Electron-electron energy (au per primitive cell)" << endl;
                   os << ef_.ehart() << endl;
                   os << "Ion-ion energy (au per primitive cell)" << endl;
                   //os << ef_.esr()-ef_.eself() << endl;
                   os << eewald << endl;
                   os << "Number of electrons per primitive cell" << endl;
                   os << s->wf.nel() << endl;
                   os << endl;
                   s->atoms.print_casino(os);
                }
             }
          }
          s->ctxt_.barrier();
          s->wf.print_casino(os,ispin,kk);
          if (sdctxt != 0)
             if ( sdctxt->active() )
                if ( sdctxt->myproc() == 0 ) 
                   os.close();
          s->ctxt_.barrier();
       }
    }
*/

  }
  else if (encoding == "vmd" ) {
     s->wf.print_vmd(filename,s->atoms);

    // print density
    ChargeDensity cd_(*s);

    /*
    // load mixed charge density into cd_
    for (int ispin = 0; ispin < s->wf.nspin(); ispin++) 
      for ( int i=0; i < s->rhog_last[ispin].size(); i++ )
        cd_.rhog[ispin][i] = s->rhog_last[ispin][i];
    cd_.update_rhor();
    */
    // if ultrasoft, calculate position-dependent functions
    if (s->ctrl.ultrasoft)
       cd_.update_usfns();
    cd_.update_density();

    const Context* wfctxt = s->wf.spincontext(0);
    const Context* vctxt = &cd_.vcontext();

    //ewd DEBUG
    assert(wfctxt->nprow() == vctxt->nprow());

    FourierTransform* ft_ = cd_.vft();

    if (wfctxt->mycol() == 0) {

      vector<double> rhortmp(ft_->np012loc());
      for (int j = 0; j < ft_->np012loc(); j++)
        rhortmp[j] = cd_.rhor[0][j];
    
      for ( int i = 0; i < wfctxt->nprow(); i++ ) {
        if ( i == wfctxt->myrow() ) {
          int size = ft_->np012loc();
          wfctxt->isend(1,1,&size,1,0,0);
          if (size > 0)
             wfctxt->dsend(size,1,&rhortmp[0],1,0,0);
        }
      }
      if ( wfctxt->oncoutpe() ) {
        vector<double> tmprecv(ft_->np012());
        int recvoffset = 0;

        D3vector a0 = s->wf.cell().a(0);
        D3vector a1 = s->wf.cell().a(1);
        D3vector a2 = s->wf.cell().a(2);
        const int np0 = ft_->np0();
        const int np1 = ft_->np1();
        const int np2 = ft_->np2();
        D3vector dft0 = a0/(double)np0;
        D3vector dft1 = a1/(double)np1;
        D3vector dft2 = a2/(double)np2;

        string denfile = filestr + ".density.cube";
        ofstream os;
        os.open(denfile.c_str(),ofstream::out);
        os.setf(ios::scientific,ios::floatfield);
        os << setprecision(8);
        
        for ( int i = 0; i < wfctxt->nprow(); i++ ) {
          int size = 0;
          wfctxt->irecv(1,1,&size,1,i,0);
          if (size > 0)
             wfctxt->drecv(size,1,&tmprecv[recvoffset],1,i,0);
          recvoffset += size;

          if (i==0) {
            // write out VMD CUBE format header
            os << "Qbox wavefunction in VMD CUBE format" << endl;
            os << "  electron density" << endl;

            // get atom positions
            AtomSet& as = s->atoms;
            vector<vector<double> > rion;
            rion.resize(as.nsp());
            int natoms_total = 0;
            for ( int is = 0; is < as.nsp(); is++ ) {
              rion[is].resize(3*as.na(is));
              natoms_total += as.na(is);
            }
            as.get_positions(rion,true);
            D3vector origin(0.0,0.0,0.0);
            os << natoms_total << " " << origin << endl;

            // print FFT grid info
            os << np0 << " " << dft0 << endl;
            os << np1 << " " << dft1 << endl;
            os << np2 << " " << dft2 << endl;

            // print atom coordinates
            for ( int is = 0; is < as.nsp(); is++ ) {
              const int atnum = as.atomic_number(is);
              double atnumd = (double)atnum;
              for ( int ia = 0; ia < as.na(is); ia++ ) 
                os << atnum << " " << atnumd << " " << rion[is][3*ia] << " " << rion[is][3*ia+1] << " " << rion[is][3*ia+2] << endl;
            }
          }
        }

        // write density data to file
        int cnt = 0;
        for (int ii = 0; ii < np0; ii++) {
          ostringstream oss;
          oss.setf(ios::scientific,ios::floatfield);
          oss << setprecision(5);
          for (int jj = 0; jj < np1; jj++) {
            for (int kk = 0; kk < np2; kk++) {
              int index = ii + jj*np0 + kk*np0*np1;
              oss << tmprecv[index] << " ";
              cnt++;
              if (cnt >= 6) {
                cnt = 0;
                oss << endl;
              }
            }
          }
          string tos = oss.str();
          os.write(tos.c_str(),tos.length());
        }
        os.close();
      }
    }
  }

  savetm.stop();

  double time = savetm.cpu();
  double tmin = time;
  double tmax = time;
    
  s->ctxt_.dmin(1,1,&tmin,1);
  s->ctxt_.dmax(1,1,&tmax,1);
  if ( ui->oncoutpe() ) {
    cout << "<!--  save timing : " << setprecision(4) << setw(9) << tmin
         << " " << setprecision(4) << setw(9) << tmax << " -->" << endl;
  }

  return 0;
}
