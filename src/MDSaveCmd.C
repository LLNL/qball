////////////////////////////////////////////////////////////////////////////////
//
// MDSaveCmd.C:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: MDSaveCmd.C,v 1.23 2010/05/12 20:05:25 draeger1 Exp $


#include "MDSaveCmd.h"
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
#include <sys/stat.h>

#ifdef USE_CSTDIO_LFS
#include <cstdio>
#include <cstdlib>
#include <sstream>
#endif
using namespace std;

////////////////////////////////////////////////////////////////////////////////
int MDSaveCmd::action(int argc, char **argv) {

  if ( !(argc>=1 && argc<=3 ) ) {
    if ( ui->oncoutpe() ) {
      cout << "  <!-- use: mdsave [encoding options] [filebase] -->" << endl;
      cout << "  <!--     encoding options:                                          -->" << endl;
      cout << "  <!--       -dump:  one binary file for each process                 -->" << endl;
      cout << "  <!--       -states:  one file for each state                        -->" << endl;
    } 
    return 1;
  }
  
  // set default encoding and format
  string format = "binary";
  string encoding = "dump";
  string dirbase = "md.";
  char* filename = "mdchk";

  // parse arguments
  for ( int i = 1; i < argc; i++ ) {
    string arg(argv[i]);
    
    if ( arg=="-dump" )
      encoding = "dump";
    else if ( arg=="-states" )
      encoding = "states";
    else if ( arg=="-binary" )
      format = "binary";
    else if ( arg[0] != '-' && i == argc-1 )
      filename = argv[i];
    else {
      if ( ui->oncoutpe() ) {
         cout << "  <!-- use: mdsave [encoding options] [filebase]   -->" << endl;
         cout << "  <!--     encoding options:                                          -->" << endl;
         cout << "  <!--       -dump:  one binary file for each process                 -->" << endl;
         cout << "  <!--       -states:  one file for each state                        -->" << endl;
      } 
      return 1;
    }
  }
  
  if ( filename == 0 ) {
     if ( ui->oncoutpe() ) {
        cout << "  <!-- use: mdsave [encoding options] [filebase]   -->" << endl;
        cout << "  <!--     encoding options:                                          -->" << endl;
        cout << "  <!--       -dump:  one binary file for each process                 -->" << endl;
        cout << "  <!--       -states:  one file for each state                        -->" << endl;
     }
     return 1;
  }

  if (s->ctrl.timer_hit) {
     if (s->ctrl.timer_mdsavecmd) {
        if ( ui->oncoutpe() )
           cout << " <!-- MDSaveCmd: run_timer exceeded: all save commands beyond first will be ignored. -->" << endl;
        return 0;
     }
     else
        s->ctrl.timer_mdsavecmd = true;
  }
  
  Timer savetm;
  savetm.start();

  // create output directory if it doesn't exist
  ostringstream oss;
  oss.width(7);  oss.fill('0');  oss << s->ctrl.mditer;
  string dirstr = dirbase + oss.str();
  if ( ui->oncoutpe() )
  {
     int mode = 0775;
     struct stat statbuf;
     int rc = stat(dirstr.c_str(), &statbuf);
     if (rc == -1)
     {
        cout << "Creating directory: " << dirstr << endl;
        rc = mkdir(dirstr.c_str(), mode);
        rc = stat(dirstr.c_str(), &statbuf);
     }
     if (rc != 0 || !(statbuf.st_mode))
     {
        cout << "<ERROR> Can't stat directory " << dirstr << " </ERROR> " << endl;
        MPI_Abort(MPI_COMM_WORLD,2);
     }
  }
  string filebase(filename);
  string filestr = dirstr + "/" + filebase;

  // binary output
  if (encoding == "dump" ) {
     s->wf.write_dump(filestr,s->ctrl.mditer);
     if (s->ctrl.tddft_involved)
     {
        if ( ui->oncoutpe() )
           cout << "<!-- MDSaveCmd:  wf write finished, writing hamil_wf... -->" << endl;
        // write s->hamil_wf
        string hamwffile = filestr + "hamwf";
        s->hamil_wf->write_dump(hamwffile,-1);
     }
     else
     {
        // ewd:  write rhor_last to file
        ChargeDensity cd_(*s);

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
                    wfctxt->dsend(size,1,&rhortmp[0],1,0,0);
                 }
              }
              if ( wfctxt->onpe0() ) {
                 for ( int i = 0; i < wfctxt->nprow(); i++ ) {
                    int size = 0;
                    wfctxt->irecv(1,1,&size,1,i,0);
                    wfctxt->drecv(size,1,&rhortmp[0],1,i,0);
                    os.write((char*)&rhortmp[0],sizeof(double)*size);
                 }
                 os.close();
              }
           }
        }
     }
  }    
  else if (encoding == "states" ) {
     s->wf.write_states(filename,format,s->ctrl.mditer);
     
     if (s->ctrl.tddft_involved)
     {
        // write s->hamil_wf
        string hamwffile = filestr + "hamwf";
        s->hamil_wf->write_states(hamwffile,format,-1);
     }
     else
     {
        ChargeDensity cd_(*s);
        
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
                    wfctxt->dsend(size,1,&rhortmp[0],1,0,0);
                 }
              }
              if ( wfctxt->onpe0() ) {
                 for ( int i = 0; i < wfctxt->nprow(); i++ ) {
                    int size = 0;
                    wfctxt->irecv(1,1,&size,1,i,0);
                    wfctxt->drecv(size,1,&rhortmp[0],1,i,0);
                    os.write((char*)&rhortmp[0],sizeof(double)*size);
                 }
                 os.close();
              }
           }
        }
     }
  }

  // write .sys file
  if (ui->oncoutpe() ){
    
     string sysfilename = dirstr + "/" + "mdsave.sys";
     ofstream os;
     os.open(sysfilename.c_str(),ofstream::out);

    // cell info
    string cmd("set cell ");
    s->wf.cell().printsys(os,cmd);

    // ref cell info, if necessary
    if ( s->wf.refcell().volume() != 0.0 ) {
      string refcmd("set ref_cell ");
      s->wf.refcell().printsys(os,refcmd);
    }

    // species info
    const int nspqm_ = s->atoms.nsp();
    for (int i=0; i<nspqm_; i++)
      s->atoms.species_list[i]->printsys(os);

    const int nspmm_ = s->atoms.nsp_mm();
    for (int i=0; i<nspmm_; i++)
      s->atoms.mmspecies_list[i]->printsys(os);

    // atom coordinates and info
    s->atoms.printsys(os);
    
    os.close();
  }
  
  savetm.stop();

  double time = savetm.cpu();
  double tmin = time;
  double tmax = time;
    
  s->ctxt_.dmin(1,1,&tmin,1);
  s->ctxt_.dmax(1,1,&tmax,1);
  if ( ui->oncoutpe() ) {
    cout << "<!--  mdsave timing : " << setprecision(4) << setw(9) << tmin
         << " " << setprecision(4) << setw(9) << tmax << " -->" << endl;
  }

  return 0;
}
