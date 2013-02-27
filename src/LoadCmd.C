////////////////////////////////////////////////////////////////////////////////
//
// LoadCmd.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: LoadCmd.C,v 1.21 2010/01/26 03:29:45 draeger1 Exp $

#include "LoadCmd.h"
#include "SampleReader.h"
#include "Sample.h"
#include "Timer.h"
#include "Context.h"
#include "ChargeDensity.h"
#include "FourierTransform.h"
#include "Basis.h"
#include "fstream"
#include <iostream>
#include <iomanip>
//ewd DEBUG
#include "SlaterDet.h"
#include "Wavefunction.h"
//ewd DEBUG
using namespace std;

////////////////////////////////////////////////////////////////////////////////
int LoadCmd::action(int argc, char **argv) {

  if ( !(argc>=2 && argc<=4 ) ) {
    if ( ui->oncoutpe() )
      cout << "  <!-- use: load [-dump|-fast|-states|-text|-xml] [-serial] filename -->" 
           << endl;
    return 1;
  }
  
  // set default encoding
  string encoding = "dump";
  char* filename = 0;
  bool serial = false;
  bool readvel = false;
  
  // parse arguments
  for ( int i = 1; i < argc; i++ ) {
    string arg(argv[i]);
    
    if ( arg=="-text" )
      encoding = "text";
    else if ( arg=="-dump" )
      encoding = "dump";
    else if ( arg=="-fast" )
      encoding = "fast";
    else if ( arg=="-states" )
      encoding = "states";
    else if ( arg=="-xml" )
      encoding = "base64";
    else if ( arg=="-vel" )
      readvel = true;
    else if ( arg=="-serial" )
      serial = true;
    else if ( arg[0] != '-' && i == argc-1 )
      filename = argv[i];
    else {
      if ( ui->oncoutpe() )
        cout << "  <!-- use: load [-dump|-states|-text|-xml] [-serial] filename -->" 
             << endl;
      return 1;
    }
  }
  
  if ( filename == 0 ) {
    if ( ui->oncoutpe() )
      cout << "  <!-- use: load [-dump|-states|-text|-xml] [-serial] filename -->" 
           << endl;
    return 1;
  }

  Timer loadtm;
  loadtm.start();

  string filestr(filename);

  /////  DUMP CHECKPOINTING  /////
  if (encoding == "dump" ) {
     s->wf.read_dump(filestr);
     s->wf.read_mditer(filestr,s->ctrl.mditer);
     if ( ui->oncoutpe())
        cout << "<!-- LoadCmd:  setting MD iteration count to " << s->ctrl.mditer << ". -->" << endl;       

    if (s->ctrl.extra_memory >= 3)
      s->wf.set_highmem();    
    if (s->ctrl.ultrasoft)
      s->wf.init_usfns(&s->atoms);
      
    if (s->ctrl.tddft_involved)
    {
        string hamwffile = filestr + "hamwf";
        if ( s->hamil_wf == 0 ) {
          s->hamil_wf = new Wavefunction(s->wf);
          (*s->hamil_wf) = s->wf;
          (*s->hamil_wf).update_occ(0.0,0);
          //s->hamil_wf->clear();
        }
        s->hamil_wf->read_dump(hamwffile);
    }
    else
    {
       // ewd:  read rhor_last from file, send to appropriate pes
       ChargeDensity cd_(*s);
       cd_.update_density();

       //ewd DEBUG
       if (s->ctrl.ultrasoft)    
          cd_.update_usfns();
    
       const Context* wfctxt = s->wf.wfcontext();
       const Context* vctxt = &cd_.vcontext();
       FourierTransform* ft_ = cd_.vft();
       const double omega = cd_.vbasis()->cell().volume();
       const int nspin = s->wf.nspin();
       s->rhog_last.resize(nspin);
       for (int ispin = 0; ispin < nspin; ispin++) {
          valarray<complex<double> > rhortmp(ft_->np012loc());
          int rhorsize = cd_.rhor[ispin].size();
          ifstream is;
          string rhorfile;

          if (nspin == 1)
             rhorfile = filestr + ".lastrhor";
          else {
             ostringstream oss;
             oss.width(1);  oss.fill('0');  oss << ispin;
             rhorfile = filestr + ".s" + oss.str() + ".lastrhor";
          }

          int file_exists = 0;
          if (wfctxt->myrow() == 0) {
             is.open(rhorfile.c_str(),ifstream::binary);
             if (is.is_open()) 
                file_exists = 1;
             else 
                file_exists = -1;
             
             // send file_exists flag to other pes in column
             for ( int i = 0; i < wfctxt->nprow(); i++ ) {
                wfctxt->isend(1,1,&file_exists,1,i,wfctxt->mycol());
             }
          }
          wfctxt->irecv(1,1,&file_exists,1,0,wfctxt->mycol());
          if (file_exists == 1) { 
             // send local charge density size from each pe in column to row 0
             if (wfctxt->mycol() == 0) {
                wfctxt->isend(1,1,&rhorsize,1,0,0);
             }
             
             if (wfctxt->onpe0()) {
                cout << "<!-- LoadCmd:  loading mixed charge density from file. -->" << endl;
                // hack to make checkpointing work w. BlueGene compilers
#ifdef BGQ
                char* tmpfilename = new char[256];
                is.read(tmpfilename,sizeof(char)*rhorfile.length());
#endif          
                for ( int i = 0; i < wfctxt->nprow(); i++ ) {
                   int tmpsize;
                   wfctxt->irecv(1,1,&tmpsize,1,i,0);
                   vector<double> tmpr(tmpsize);
                   // read this portion of charge density, send to all pes in row i
                   is.read((char*)&tmpr[0],sizeof(double)*tmpsize);
                   if (tmpsize > 0) 
                      for ( int j = 0; j < wfctxt->npcol(); j++ ) 
                         wfctxt->dsend(tmpsize,1,&tmpr[0],1,i,j);
                }
             }
             if (rhorsize > 0) {
                vector<double> rhorrecv(rhorsize);
                wfctxt->drecv(rhorsize,1,&rhorrecv[0],1,0,0);
                // copy rhor to complex<double> array for Fourier transform
                for (int j = 0; j < rhorsize; j++)
                   rhortmp[j] = complex<double>(omega*rhorrecv[j],0.0);
             }
             int rhogsize = cd_.rhog[ispin].size();
             vector<complex<double> > rhogtmp(rhogsize);
             ft_->forward(&rhortmp[0],&rhogtmp[0]);
             s->rhog_last[ispin].resize(rhogsize);
             //complex<double> *rhogp = &s->rhog_last[ispin];
             //for (int j = 0; j < rhogsize; j++)
             //  rhogp[j] = rhogtmp[j];
             for (int j = 0; j < rhogsize; j++)
                s->rhog_last[ispin][j] = rhogtmp[j];
          }
          else {
             if ( ui->oncoutpe() )
                cout << "<!-- LoadCmd: mixed charge density checkpoint file not found. -->" << endl;
          }
          if (wfctxt->myrow() == 0) 
             is.close();
       }
    }
    
    if (readvel) { 
      const string atoms_dyn = s->ctrl.atoms_dyn;
      const bool compute_forces = ( atoms_dyn != "LOCKED" );
      if (compute_forces) {
        string wfvfile = filestr + "wfv";
        if ( s->wfv == 0 ) {
          s->wfv = new Wavefunction(s->wf);
          s->wfv->clear();
        }
        s->wfv->read_dump(wfvfile);
      }
      else {
        if ( ui->oncoutpe() )
          cout << "<WARNING>LoadCmd:  wavefunction velocities only available when atoms_dyn set, can't load.</WARNING>" << endl;
      }
    }

    if (serial)
      if ( ui->oncoutpe() )
        cout << "<!-- LoadCmd:  serial flag only used with xml input, ignoring. -->" << endl;
  }
  /////  FAST CHECKPOINTING  /////
  if (encoding == "fast" ) {
     s->wf.read_fast(filestr);
     s->wf.read_mditer(filestr,s->ctrl.mditer);
     if ( ui->oncoutpe())
        cout << "<!-- LoadCmd:  setting MD iteration count to " << s->ctrl.mditer << ". -->" << endl;       

    if (s->ctrl.extra_memory >= 3)
      s->wf.set_highmem();    
    if (s->ctrl.ultrasoft)
      s->wf.init_usfns(&s->atoms);
      
    if (s->ctrl.tddft_involved)
    {
        string hamwffile = filestr + "hamwf";
        if ( s->hamil_wf == 0 ) {
          s->hamil_wf = new Wavefunction(s->wf);
          (*s->hamil_wf) = s->wf;
          (*s->hamil_wf).update_occ(0.0,0);
          //s->hamil_wf->clear();
        }
        s->hamil_wf->read_fast(hamwffile);
    }
    else
    {
       // ewd:  read rhor_last from file, send to appropriate pes
       ChargeDensity cd_(*s);

       //ewd DEBUG
       if (s->ctrl.ultrasoft)    
          cd_.update_usfns();

       cd_.update_density();
    
       const Context* wfctxt = s->wf.wfcontext();
       const Context* vctxt = &cd_.vcontext();
       FourierTransform* ft_ = cd_.vft();
       const double omega = cd_.vbasis()->cell().volume();
       const int nspin = s->wf.nspin();
       s->rhog_last.resize(nspin);
       for (int ispin = 0; ispin < nspin; ispin++) {
          valarray<complex<double> > rhortmp(ft_->np012loc());
          int rhorsize = cd_.rhor[ispin].size();
          ifstream is;
          string rhorfile;

          if (nspin == 1)
             rhorfile = filestr + ".lastrhor";
          else {
             ostringstream oss;
             oss.width(1);  oss.fill('0');  oss << ispin;
             rhorfile = filestr + ".s" + oss.str() + ".lastrhor";
          }

          int file_exists = 0;
          if (wfctxt->myrow() == 0) {
             is.open(rhorfile.c_str(),ifstream::binary);
             if (is.is_open()) 
                file_exists = 1;
             else 
                file_exists = -1;
             
             // send file_exists flag to other pes in column
             for ( int i = 0; i < wfctxt->nprow(); i++ ) {
                wfctxt->isend(1,1,&file_exists,1,i,wfctxt->mycol());
             }
          }
          wfctxt->irecv(1,1,&file_exists,1,0,wfctxt->mycol());
          if (file_exists == 1) { 
             // send local charge density size from each pe in column to row 0
             if (wfctxt->mycol() == 0) {
                wfctxt->isend(1,1,&rhorsize,1,0,0);
             }
             
             if (wfctxt->onpe0()) {
                cout << "<!-- LoadCmd:  loading mixed charge density from file. -->" << endl;
                // hack to make checkpointing work w. BlueGene compilers
#ifdef BGQ
                char* tmpfilename = new char[256];
                is.read(tmpfilename,sizeof(char)*rhorfile.length());
#endif          
                for ( int i = 0; i < wfctxt->nprow(); i++ ) {
                   int tmpsize;
                   wfctxt->irecv(1,1,&tmpsize,1,i,0);
                   vector<double> tmpr(tmpsize);
                   // read this portion of charge density, send to all pes in row i
                   is.read((char*)&tmpr[0],sizeof(double)*tmpsize);
                   if (tmpsize > 0) 
                      for ( int j = 0; j < wfctxt->npcol(); j++ ) 
                         wfctxt->dsend(tmpsize,1,&tmpr[0],1,i,j);
                }
             }
             if (rhorsize > 0) {
                vector<double> rhorrecv(rhorsize);
                wfctxt->drecv(rhorsize,1,&rhorrecv[0],1,0,0);
                // copy rhor to complex<double> array for Fourier transform
                for (int j = 0; j < rhorsize; j++)
                   rhortmp[j] = complex<double>(omega*rhorrecv[j],0.0);
             }
             int rhogsize = cd_.rhog[ispin].size();
             vector<complex<double> > rhogtmp(rhogsize);
             ft_->forward(&rhortmp[0],&rhogtmp[0]);
             s->rhog_last[ispin].resize(rhogsize);
             //complex<double> *rhogp = &s->rhog_last[ispin];
             //for (int j = 0; j < rhogsize; j++)
             //  rhogp[j] = rhogtmp[j];
             for (int j = 0; j < rhogsize; j++)
                s->rhog_last[ispin][j] = rhogtmp[j];
          }
          else {
             if ( ui->oncoutpe() )
                cout << "<!-- LoadCmd: mixed charge density checkpoint file not found. -->" << endl;
          }
          if (wfctxt->myrow() == 0) 
             is.close();
       }
    }
    
    if (readvel) { 
      const string atoms_dyn = s->ctrl.atoms_dyn;
      const bool compute_forces = ( atoms_dyn != "LOCKED" );
      if (compute_forces) {
        string wfvfile = filestr + "wfv";
        if ( s->wfv == 0 ) {
          s->wfv = new Wavefunction(s->wf);
          s->wfv->clear();
        }
        s->wfv->read_fast(wfvfile);
      }
      else {
        if ( ui->oncoutpe() )
          cout << "<WARNING>LoadCmd:  wavefunction velocities only available when atoms_dyn set, can't load.</WARNING>" << endl;
      }
    }

    if (serial)
      if ( ui->oncoutpe() )
        cout << "<!-- LoadCmd:  serial flag only used with xml input, ignoring. -->" << endl;
  }
  /////  STATES CHECKPOINTING  /////
  else if (encoding == "states" ) {
     s->wf.read_states(filestr);
     s->wf.read_mditer(filestr,s->ctrl.mditer);
     if ( ui->oncoutpe())
        cout << "<!-- LoadCmd:  setting MD iteration count to " << s->ctrl.mditer << ". -->" << endl;       

    if (s->ctrl.extra_memory >= 3)
      s->wf.set_highmem();    
    if (s->ctrl.ultrasoft)
      s->wf.init_usfns(&s->atoms);

    if (s->ctrl.tddft_involved)
    {
        string hamwffile = filestr + "hamwf";
        if ( s->hamil_wf == 0 ) {
          s->hamil_wf = new Wavefunction(s->wf);
          (*s->hamil_wf) = s->wf;
          (*s->hamil_wf).update_occ(0.0,0);
          //s->hamil_wf->clear();
        }
        s->hamil_wf->read_states(hamwffile);
    }
    else
    {
       // ewd:  read rhor_last from file, send to appropriate pes
       ChargeDensity cd_(*s);

       //ewd DEBUG
       if (s->ctrl.ultrasoft)    
          cd_.update_usfns();
       
       cd_.update_density();

       const Context* wfctxt = s->wf.wfcontext();
       const Context* vctxt = &cd_.vcontext();
       FourierTransform* ft_ = cd_.vft();
       const double omega = cd_.vbasis()->cell().volume();
       const int nspin = s->wf.nspin();
       s->rhog_last.resize(nspin);
       for (int ispin = 0; ispin < nspin; ispin++) {
          valarray<complex<double> > rhortmp(ft_->np012loc());
          int rhorsize = cd_.rhor[ispin].size();
          ifstream is;
          string rhorfile;

          if (nspin == 1)
             rhorfile = filestr + ".lastrhor";
          else {
             ostringstream oss;
             oss.width(1);  oss.fill('0');  oss << ispin;
             rhorfile = filestr + ".s" + oss.str() + ".lastrhor";
          }

          int file_exists = 0;
          if (wfctxt->myrow() == 0) {
             is.open(rhorfile.c_str(),ifstream::binary);
             if (is.is_open()) 
                file_exists = 1;
             else 
                file_exists = -1;

             // send file_exists flag to other pes in column
             for ( int i = 0; i < wfctxt->nprow(); i++ ) {
                wfctxt->isend(1,1,&file_exists,1,i,wfctxt->mycol());
             }
          }
          wfctxt->irecv(1,1,&file_exists,1,0,wfctxt->mycol());
          if (file_exists == 1) { 
             // send local charge density size from each pe in column to row 0
             if (wfctxt->mycol() == 0) {
                wfctxt->isend(1,1,&rhorsize,1,0,0);
             }
             
             if (wfctxt->onpe0()) {
                cout << "<!-- LoadCmd:  loading mixed charge density from file. -->" << endl;
                // hack to make checkpointing work w. BlueGene compilers
#ifdef BGQ
                char* tmpfilename = new char[256];
                is.read(tmpfilename,sizeof(char)*rhorfile.length());
#endif
                
                for ( int i = 0; i < wfctxt->nprow(); i++ ) {
                   int tmpsize;
                   wfctxt->irecv(1,1,&tmpsize,1,i,0);
                   vector<double> tmpr(tmpsize);
                   // read this portion of charge density, send to all pes in row i
                   is.read((char*)&tmpr[0],sizeof(double)*tmpsize);
                   if (tmpsize > 0) 
                      for ( int j = 0; j < wfctxt->npcol(); j++ ) 
                         wfctxt->dsend(tmpsize,1,&tmpr[0],1,i,j);
                }
             }
             if (rhorsize > 0) {
                vector<double> rhorrecv(rhorsize);
                wfctxt->drecv(rhorsize,1,&rhorrecv[0],1,0,0);
                // copy rhor to complex<double> array for Fourier transform
                for (int j = 0; j < rhorsize; j++)
                   rhortmp[j] = complex<double>(omega*rhorrecv[j],0.0);
             }
             int rhogsize = cd_.rhog[ispin].size();
             vector<complex<double> > rhogtmp(rhogsize);
             ft_->forward(&rhortmp[0],&rhogtmp[0]);
             s->rhog_last[ispin].resize(rhogsize);
             //complex<double> *rhogp = &s->rhog_last[ispin];
             //for (int j = 0; j < rhogsize; j++)
             //  rhogp[j] = rhogtmp[j];
             for (int j = 0; j < rhogsize; j++)
                s->rhog_last[ispin][j] = rhogtmp[j];
          }
          else {
             if ( ui->oncoutpe() )
                cout << "<!-- LoadCmd: mixed charge density checkpoint file not found. -->" << endl;
          }
          if (wfctxt->myrow() == 0) 
             is.close();
       }
    }
    
    if (readvel) { 
      const string atoms_dyn = s->ctrl.atoms_dyn;
      const bool compute_forces = ( atoms_dyn != "LOCKED" );
      if (compute_forces) {
        string wfvfile = filestr + "wfv";
        if ( s->wfv == 0 ) {
          s->wfv = new Wavefunction(s->wf);
          s->wfv->clear();
        }
        s->wfv->read_states(wfvfile);
      }
      else {
        if ( ui->oncoutpe() )
          cout << "<WARNING>LoadCmd:  wavefunction velocities only available when atoms_dyn set, can't load.</WARNING>" << endl;
      }
    }

    if (serial)
      if ( ui->oncoutpe() )
        cout << "<!-- LoadCmd:  serial flag only used with xml input, ignoring. -->" << endl;

    }
  /////  XML CHECKPOINTING  /////
  else {

    //ewd DEBUG:  try reshaping context so loading works for nspin > 1 and/or nkpts > 1
    const int nrowmax_orig = s->wf.nrowmax();
    const int npes = s->ctxt_.size();

    s->wf.set_nrowmax(npes);

    if ( ui->oncoutpe() )
      cout << " <!-- LoadCmd: loading from " << filename << " -->" << endl;
    
    if (!s->wf.hasdata())
      s->wf.set_hasdata(true);

    SampleReader s_reader(s->ctxt_);
  
    if ( ui->oncoutpe() )
      cout << " <!--" << endl;
    
    try {
      s_reader.readSample(*s,filename,serial);
    }
    catch ( const SampleReaderException& e ) {
      cout << " SampleReaderException caught in LoadCmd:" << endl;
      cout << e.msg << endl;
    }
    catch (...) {
      cout << " LoadCmd: cannot load Sample" << endl;
    }
    s->ctxt_.barrier();
    //ewd DEBUG
    s->wf.set_nrowmax(nrowmax_orig);

    // If only <atomset> was read, set nel for the wavefunction
    //cout << " LoadCmd: atoms.nel() = " << s->atoms.nel() << endl;
    //cout << " LoadCmd: wf.nel() =    " << s->wf.nel() << endl;
    if ( s->wf.nel() != s->atoms.nel() ) {
      s->wf.set_nel(s->atoms.nel());
      s->wf.update_occ(0.0,0);
    }
    //cout << " LoadCmd: atoms.nel() = " << s->atoms.nel() << endl;
    //cout << " LoadCmd: wf.nel() =    " << s->wf.nel() << endl;

    if ( ui->oncoutpe() )
      cout << " -->" << endl;
  }

  if (s->ctrl.extra_memory >= 3)
    s->wf.set_highmem();    
  if (s->ctrl.ultrasoft)
    s->wf.init_usfns(&s->atoms);
        
  loadtm.stop();

  double time = loadtm.cpu();
  double tmin = time;
  double tmax = time;
    
  s->ctxt_.dmin(1,1,&tmin,1);
  s->ctxt_.dmax(1,1,&tmax,1);
  if ( ui->oncoutpe() ) {
    cout << "<!--  load timing : " << setprecision(4) << setw(9) << tmin
         << " " << setprecision(4) << setw(9) << tmax << " -->" << endl;
  }

  return 0;
}
