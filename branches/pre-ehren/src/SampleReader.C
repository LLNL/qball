////////////////////////////////////////////////////////////////////////////////
//
// SampleReader.C:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SampleReader.C,v 1.12 2010/01/07 18:01:48 draeger1 Exp $

#define DEBUG 0

#include "Sample.h"
#include "SampleReader.h"
#include "SpeciesReader.h"
#include "Species.h"
#include "Basis.h"
#include "FourierTransform.h"
#include "SlaterDet.h"

#include "XMLGFPreprocessor.h"

#include "Timer.h"

#include <cassert>
#include <string>
#include <iostream>
#include <fstream>
#include <sys/stat.h>

#if USE_XERCES
#include "SampleHandler.h"
#include "StructuredDocumentHandler.h"
#include <xercesc/util/XMLUniDefs.hpp>
#include <xercesc/sax2/Attributes.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/framework/LocalFileInputSource.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>
#include <xercesc/framework/MemBufInputSource.hpp>
using namespace xercesc;
#endif
using namespace std;

////////////////////////////////////////////////////////////////////////////////
SampleReader::SampleReader(const Context& ctxt) : ctxt_(ctxt) {}

////////////////////////////////////////////////////////////////////////////////
void SampleReader::readSample (Sample& s, const string uri, bool serial)
{
  Timer tm;
  tm.start();
#if USE_XERCES
  const char* encodingName = "UTF-8";
  //SAX2XMLReader::ValSchemes valScheme = SAX2XMLReader::Val_Auto;
  //SAX2XMLReader::ValSchemes valScheme = SAX2XMLReader::Val_Always;
  SAX2XMLReader::ValSchemes valScheme = SAX2XMLReader::Val_Never;
  bool expandNamespaces = false;
  bool doNamespaces = true;
  bool doSchema = true;
  bool schemaFullChecking = true;
  bool namespacePrefixes = false;
  SAX2XMLReader* parser = 0;

  Wavefunction* wfvtmp = new Wavefunction(s.wf);
  Wavefunction* current_wf = &s.wf;
  vector<vector<vector<double> > > dmat;
  int nx, ny, nz; // size of <grid> in wavefunction
  int nspin = s.wf.nspin();
  int current_ispin;
  int current_ikp;
  //ewd DEBUG:  size this for nspin=2 case
  //dmat.resize(nspin);
  dmat.resize(2);

  MemBufInputSource* memBufIS = 0;

  // XMLPlatformUtils initialization on task 0 only
  int ierr = 0;
  if ( ctxt_.onpe0() )
  {
    try
    {
      XMLPlatformUtils::Initialize();
    }

    catch (const XMLException& toCatch)
    {
      cout << "  Sample::readSample: Error during XML initialization :\n"
           << StrX(toCatch.getMessage()) << endl;
      ierr = 1;
    }
    ctxt_.ibcast_send(1,1,&ierr,1);
  }
  else
  {
    ctxt_.ibcast_recv(1,1,&ierr,1,0,0);
  }
  //cout << ctxt_.mype() <<": SampleReader: ierr=" << ierr << endl;
  if ( ierr > 0 )
    throw SampleReaderException("error in XMLPlatformUtils::Initialize");

  // Determine if uri is a local file
  struct stat statbuf;
  bool read_from_file = !stat(uri.c_str(),&statbuf);

  // check for serial override
  read_from_file &= !serial;

  string xmlcontent;
  DoubleMatrix gfdata(ctxt_);
  if ( read_from_file )
  {
    const char* const filename = uri.c_str();
    if ( ctxt_.onpe0() )
      cout << " SampleReader: reading from file "
           << filename << " size: "
           << statbuf.st_size << endl;

    XMLGFPreprocessor xmlgfp;

    xmlgfp.process(filename,gfdata,xmlcontent);

    if ( ctxt_.onpe0() )
    {
      cout << " xmlcontent.size(): " << xmlcontent.size()
           << endl;
#if DEBUG
      cout << ctxt_.mype() << ": xmlcontent: " << xmlcontent << endl;
#endif
      memBufIS = new MemBufInputSource
        ( (const XMLByte*) &xmlcontent[0], xmlcontent.size(), "buf_id", false);
    }
  }

  bool read_wf = false;
  bool read_wfv = false;
  // initialize wavefunction_velocity in case it appears in the sample file

  assert(sizeof(event_type)==sizeof(int));
  // remove default kpoint in wf
  //ewd:  my version doesn't need this
  //ewds.wf.del_kpoint(D3vector(0.0,0.0,0.0));
  //ewdwfvtmp->del_kpoint(D3vector(0.0,0.0,0.0));

  if ( ctxt_.onpe0() )
  {
    parser = XMLReaderFactory::createXMLReader();
    if (valScheme == SAX2XMLReader::Val_Auto)
    {
        parser->setFeature(XMLUni::fgSAX2CoreValidation, true);
        parser->setFeature(XMLUni::fgXercesDynamic, true);
    }

    if (valScheme == SAX2XMLReader::Val_Never)
    {
        parser->setFeature(XMLUni::fgSAX2CoreValidation, false);
    }

    if (valScheme == SAX2XMLReader::Val_Always)
    {
        parser->setFeature(XMLUni::fgSAX2CoreValidation, true);
        parser->setFeature(XMLUni::fgXercesDynamic, false);
    }

    parser->setFeature(XMLUni::fgSAX2CoreNameSpaces, doNamespaces);
    parser->setFeature(XMLUni::fgXercesSchema, doSchema);
    parser->setFeature(XMLUni::fgXercesSchemaFullChecking, schemaFullChecking);
    parser->setFeature(XMLUni::fgSAX2CoreNameSpacePrefixes, namespacePrefixes);

    int errorCount = 0;
    SampleHandler* s_handler =
      new SampleHandler(s,gfdata,nx,ny,nz,dmat,*wfvtmp);

    try
    {
      StructuredDocumentHandler handler(s_handler);
      parser->setContentHandler(&handler);
      parser->setErrorHandler(&handler);

      cout << " Starting XML parsing" << endl;
      if ( read_from_file )
        parser->parse(*memBufIS);
      else
        parser->parse(uri.c_str());
      cout << " XML parsing done" << endl;

      errorCount = parser->getErrorCount();
    }

    catch (const XMLException& toCatch)
    {
        cout << "\nAn error occurred\n  Error: "
             << StrX(toCatch.getMessage())
             << "\n" << endl;
        //XMLPlatformUtils::Terminate();
        //delete parser;
        //throw;
    }

    catch (const SAXParseException& e)
    {
        cout << "\na SAXParseException occurred in file "
             << StrX(e.getSystemId())
             << ", line " << e.getLineNumber()
             << ", char " << e.getColumnNumber()
             << "\n  Message: " << StrX(e.getMessage()) << endl;
        //XMLPlatformUtils::Terminate();
        //delete parser;
        //throw;
    }
    read_wf = s_handler->read_wf;
    if ( read_wf )
      cout << " wavefunction was read" << endl;
    read_wfv = s_handler->read_wfv;
    if ( read_wfv )
      cout << " wavefunction velocity was read" << endl;
    cout << " SampleReader::readSample: grid nx,ny,nz="
         << nx << " " << ny << " " << nz << endl;
    
    delete s_handler;
    delete parser;
    XMLPlatformUtils::Terminate();
    // parsing of sample is complete, send end of sample tag to tasks > 0
    event_type event = end;
    ctxt_.ibcast_send(1,1,(int*)&event,1);
  } // onpe0
  else
  {
    // tasks > 0
    // listen for Sample events
    // cout << ctxt_.mype() << ": listening..." << endl;
    bool done = false;
    event_type event = invalid;
    while ( !done )
    {
      ctxt_.ibcast_recv(1,1,(int*)&event,1,0,0);
      if ( event == end )
        done = true;
      else if ( event == unit_cell )
      {
        // unit_cell
        double buf[9];
        ctxt_.dbcast_recv(9,1,buf,9,0,0);
        D3vector a(buf[0],buf[1],buf[2]);
        D3vector b(buf[3],buf[4],buf[5]);
        D3vector c(buf[6],buf[7],buf[8]);
        s.atoms.set_cell(a,b,c);
      }
      else if ( event == species )
      {
        // Species
        string curr_name;
        ctxt_.string_bcast(curr_name,0);
        //cout << ctxt_.mype() << ": receiving species " << curr_name << endl;
        Species* sp = new Species(ctxt_,curr_name);
        SpeciesReader sp_reader(ctxt_);
        sp_reader.bcastSpecies(*sp);
        s.atoms.addSpecies(sp,curr_name);
      }
      else if ( event == atom )
      {
        // Atom
        string curr_atom_name,curr_atom_species;
        D3vector curr_position, curr_velocity;
        ctxt_.string_bcast(curr_atom_name,0);
        ctxt_.string_bcast(curr_atom_species,0);
        double buf[3];
        ctxt_.dbcast_recv(3,1,buf,3,0,0);
        curr_position = D3vector(buf[0],buf[1],buf[2]);
        ctxt_.dbcast_recv(3,1,buf,3,0,0);
        curr_velocity = D3vector(buf[0],buf[1],buf[2]);
        //cout << ctxt_.mype() << ": receiving atom " << curr_atom_name << endl;
        Atom* a = new Atom(curr_atom_name, curr_atom_species,
                           curr_position, curr_velocity);
        s.atoms.addAtom(a);
      }
      else if ( event == wavefunction )
      {
        current_ikp = 0;
        current_ispin = 0;
        // Wavefunction
        read_wf = true;
        int nel,nspin,delta_spin,nempty;
        ctxt_.ibcast_recv(1,1,&nel,1,0,0);
        ctxt_.ibcast_recv(1,1,&nspin,1,0,0);
        ctxt_.ibcast_recv(1,1,&delta_spin,1,0,0);
        ctxt_.ibcast_recv(1,1,&nempty,1,0,0);
        //cout << ctxt_.mype() << ": SampleReader: receiving wf nel=" << nel
        //     << " nspin=" << nspin << " nempty=" << nempty << endl;
        s.wf.set_nel(nel);
        s.wf.set_nspin(nspin);
        if (nspin > 1)
          s.wf.set_deltaspin(delta_spin);
        s.wf.set_nempty(nempty);
        //cout << ctxt_.mype() << ": SampleReader: nel=" << s.wf.nel()
        //     << " nst=" << s.wf.nst() << endl;

        // domain
        double buf[9];
        ctxt_.dbcast_recv(9,1,buf,9,0,0);
        D3vector a(buf[0],buf[1],buf[2]);
        D3vector b(buf[3],buf[4],buf[5]);
        D3vector c(buf[6],buf[7],buf[8]);
        UnitCell uc(a,b,c);

        // grid
        int ibuf[3];
        ctxt_.ibcast_recv(3,1,ibuf,3,0,0);
        nx = ibuf[0];
        ny = ibuf[1];
        nz = ibuf[2];
        // receive only computed ecut
        double ecut;
        ctxt_.dbcast_recv(1,1,&ecut,1,0,0);

        // reference_domain
        ctxt_.dbcast_recv(9,1,buf,9,0,0);
        D3vector ar(buf[0],buf[1],buf[2]);
        D3vector br(buf[3],buf[4],buf[5]);
        D3vector cr(buf[6],buf[7],buf[8]);
        UnitCell ruc(ar,br,cr);
        s.wf.resize(uc,ruc,ecut);
        //ewd DEBUG
        //cout << ctxt_.mype() << ": wf resized, ecut=" << ecut << endl;
      }
      else if ( event == wavefunction_velocity )
      {
        current_wf = wfvtmp; // set ptr used by the slater det section
        current_ikp = 0;
        //cout << ctxt_.mype()
        //     << ": SampleReader received wavefunction_velocity" << endl;
        read_wfv = true;
        // Wavefunction velocity
        int nel,nspin,delta_spin,nempty;
        ctxt_.ibcast_recv(1,1,&nel,1,0,0);
        ctxt_.ibcast_recv(1,1,&nspin,1,0,0);
        ctxt_.ibcast_recv(1,1,&delta_spin,1,0,0);
        ctxt_.ibcast_recv(1,1,&nempty,1,0,0);
        //cout << ctxt_.mype() << ": receiving wf nel=" << nel
        //     << " nspin=" << nspin << " nempty=" << nempty << endl;
        wfvtmp->set_nel(nel);
        wfvtmp->set_nspin(nspin);
        wfvtmp->set_deltaspin(delta_spin);
        wfvtmp->set_nempty(nempty);

        // domain
        double buf[9];
        ctxt_.dbcast_recv(9,1,buf,9,0,0);
        D3vector a(buf[0],buf[1],buf[2]);
        D3vector b(buf[3],buf[4],buf[5]);
        D3vector c(buf[6],buf[7],buf[8]);
        UnitCell uc(a,b,c);

        // grid
        int ibuf[3];
        ctxt_.ibcast_recv(3,1,ibuf,3,0,0);
        assert(nx == ibuf[0]);
        assert(ny == ibuf[1]);
        assert(nz == ibuf[2]);

        // receive only computed ecut
        double ecut;
        ctxt_.dbcast_recv(1,1,&ecut,1,0,0);

        // reference_domain
        ctxt_.dbcast_recv(9,1,buf,9,0,0);
        D3vector ar(buf[0],buf[1],buf[2]);
        D3vector br(buf[3],buf[4],buf[5]);
        D3vector cr(buf[6],buf[7],buf[8]);
        UnitCell ruc(ar,br,cr);

        current_wf->resize(uc,ruc,ecut);

        //cout << ctxt_.mype() << ": wfv resized, ecut=" << ecut << endl;

      }
      else if ( event == slater_determinant )
      {
        // process SlaterDet
        // receive current_ispin, current_ikp
        int tind[2];
        ctxt_.ibcast_recv(2,1,tind,2,0,0);
        current_ispin = tind[0];
        current_ikp = tind[1];
        
        // receive kpoint and weight
        double buf[4];
        ctxt_.dbcast_recv(4,1,buf,4,0,0);
        //ewd need to avoid adding each k-point twice when nspin = 2
        if (current_ispin == 0)
          current_wf->add_kpoint(D3vector(buf[0],buf[1],buf[2]),buf[3]);

        // receive density_matrix
        int tsize;
        s.wf.context().ibcast_recv(1,1,&tsize,1,0,0);
        vector<double> dmat_tmp(tsize);
        s.wf.context().dbcast_recv(tsize,1,&dmat_tmp[0],tsize,0,0);
        dmat[current_ispin].push_back(dmat_tmp);
        if ( !read_from_file )
        {
          // receive grid_functions
          SlaterDet* sd = current_wf->sd(current_ispin,current_ikp);
          const Basis& basis = sd->basis();
          FourierTransform ft(basis,nx,ny,nz);
          const int wftmpr_size = basis.real() ? ft.np012loc() :
            2*ft.np012loc();
          valarray<double> wftmpr(wftmpr_size);
          vector<complex<double> > wftmp(ft.np012loc());
          int size = -1;
          for ( int nloc = 0; nloc < sd->nstloc(); nloc++ )
          {
            // read grid_function nloc fragment
            //cout << sd->context();
            //cout << sd->context().mype()
            //     << ": wf receiving nloc=" << nloc << endl;
            sd->context().irecv(1,1,&size,1,0,0);
            //cout << sd->context().mype()
            //     << ": received size=" << size << endl;
            assert(size==wftmpr_size);
            sd->context().drecv(size,1,&wftmpr[0],size,0,0);
            //cout << sd->context().mype()
            //     << ": grid_function nloc=" << nloc
            //     << "received" << endl;

            // copy to complex array
            if ( basis.real() )
            {
              for ( int i = 0; i < size; i++ )
              {
                wftmp[i] = wftmpr[i];
              }
            }
            else
            {
              memcpy((void*)&wftmp[0],(void*)&wftmpr[0],wftmpr_size*
                sizeof(double));
            }
            ComplexMatrix& c = sd->c();
            ft.forward(&wftmp[0],c.valptr(c.mloc()*nloc));
            //cout << sd->context().mype()
            //     << ": grid_function read for nloc=" << nloc << endl;
          }
        }
      }
    } // while !done receiving events from node 0
  } // if-else onpe0

  // This part is now executing on all tasks
  if ( read_from_file )
  {
    // dmat may contain two sets of density matrices if wfv was read
    // Only the first set of density matrices is used

    //ewd: comment this out for now
    //if ( read_wfv )
    //  assert(dmat[0].size() == 2 * s.wf.nkp() );
    //else
    //  assert(dmat[0].size() == s.wf.nkp() );

#if DEBUG
    //ofstream gffile("gf.dat");
    //gffile << gfdata;
#endif
    if ( read_wf )
    {
      // transfer data from the gfdata matrix to the SlaterDets
#if DEBUG
      //cout << ctxt_.mype() << ": mapping gfdata matrix..."
      //     << endl;
      //cout << ctxt_.mype() << ": gfdata: (" << gfdata.m() << "x" << gfdata.n()
      //<< ") (" << gfdata.mb() << "x" << gfdata.nb() << ") blocks" << endl;
      //cout << ctxt_.mype() << ": gfdata.mloc()=" << gfdata.mloc()
      //<< " gfdata.nloc()=" << gfdata.nloc() << endl;
#endif
      // reshape wavefunction for reading from file
      int nparkpts_save = s.wf.nparallelkpts();
      if (nparkpts_save > 1)
        s.wf.set_nparallelkpts(1);      
      for ( int ispin = 0; ispin < s.wf.nspin(); ispin++ ) {
        if (s.wf.spinactive(ispin)) {
          for ( int ikp=0; ikp<s.wf.nkp(); ikp++) {
            if (s.wf.kptactive(ikp)) {
              assert(s.wf.sd(ispin,ikp) != 0);
              SlaterDet* sd = s.wf.sd(ispin,ikp);
#if DEBUG
              cout << "SampleReader: ikp=" << ikp << " nst = " << sd->nst() << endl;
              cout << "SampleReader: ikp=" << ikp << " dmat size = "
                   << dmat[ispin][ikp].size() << endl;
#endif
              // copy density matrix information
              sd->set_occ(dmat[ispin][ikp]);
              const Basis& basis = sd->basis();
              FourierTransform ft(basis,nx,ny,nz);
#if DEBUG
              cout << ctxt_.mype() << ": ft.np012loc()=" << ft.np012loc() << endl;
              cout << ctxt_.mype() << ": ft.context(): " << ft.context();
#endif

              ComplexMatrix& c = sd->c();
              // copy wf values
              // Determine the size of the temporary real matrix wftmpr
              int wftmpr_size, wftmpr_block_size;
              if ( basis.real() )
              {
                wftmpr_size = ft.np012();
                wftmpr_block_size = ft.np012loc(0);
              }
              else
              {
                wftmpr_size = 2*ft.np012();
                wftmpr_block_size = 2*ft.np012loc(0);
              }
#if DEBUG
              cout << ctxt_.mype() << ": wftmpr_size: " << wftmpr_size << endl;
              cout << ctxt_.mype() << ": wftmpr_block_size: "
                   << wftmpr_block_size << endl;
#endif
              DoubleMatrix wftmpr(sd->context(),wftmpr_size,sd->nst(),
                                  wftmpr_block_size,c.nb());

              wftmpr.getsub(gfdata,wftmpr_size,sd->nst(),0,ikp*sd->nst());

#if DEBUG
              // Check orthogonality by computing overlap matrix
              //DoubleMatrix smat(sd->context(),sd->nst(),sd->nst(),c.nb(),c.nb());
              //smat.syrk('l','t',1.0,wftmpr,0.0);
              //cout << smat;
#endif
              
              vector<complex<double> > wftmp(ft.np012loc());
              for ( int nloc = 0; nloc < sd->nstloc(); nloc++ )
              {
                // copy column of wftmpr to complex array wftmp
                if ( wftmpr_size == ft.np012() )
                {
                  // function is real and must be expanded
                  double* p = wftmpr.valptr(wftmpr.mloc()*nloc);
                  for ( int i = 0; i < ft.np012loc(); i++ )
                    wftmp[i] = p[i];
                }
                else
                {
                  // function is complex
                  double* p = wftmpr.valptr(wftmpr.mloc()*nloc);
                  for ( int i = 0; i < ft.np012loc(); i++ )
                    wftmp[i] = complex<double>(p[2*i],p[2*i+1]);
                }
                ft.forward(&wftmp[0],c.valptr(c.mloc()*nloc));
              }
              // cout << " c matrix: ikp=" << ikp << endl;
              // cout << c;
            }
          }
        }
      }
      if (nparkpts_save > 1)
        s.wf.set_nparallelkpts(nparkpts_save);
    } // if read_wf
            
    if ( read_wfv )
      {
      // transfer wfv data from the gfdata matrix to the SlaterDets
      //cout << ctxt_.mype() << ": mapping gfdata matrix..."
      //     << endl;
      //cout << ctxt_.mype() << ": gfdata: (" << gfdata.m() << "x" << gfdata.n()
      //<< ") (" << gfdata.mb() << "x" << gfdata.nb() << ") blocks" << endl;
      //cout << ctxt_.mype() << ": gfdata.mloc()=" << gfdata.mloc()
      //<< " gfdata.nloc()=" << gfdata.nloc() << endl;

      // reshape wavefunction for reading from file
      wfvtmp->set_nparallelkpts(1);      
      for ( int ispin = 0; ispin < s.wf.nspin(); ispin++ ) {
        if (s.wf.spinactive(ispin)) {
          for ( int ikp=0; ikp<s.wf.nkp(); ikp++) {
            if (s.wf.kptactive(ikp)) {
              assert(s.wf.sd(ispin,ikp) != 0);
              SlaterDet* sd = wfvtmp->sd(ispin,ikp);
              const Basis& basis = sd->basis();
              FourierTransform ft(basis,nx,ny,nz);
              //cout << ctxt_.mype() << ": ft.np012loc()=" << ft.np012loc() << endl;
              //cout << ctxt_.mype() << ": ft.context(): " << ft.context();
            
              ComplexMatrix& c = sd->c();
              // copy wf values
              // Determine the size of the temporary real matrix wftmpr
              int wftmpr_size, wftmpr_block_size;
              if ( basis.real() && ( ft.np012() == gfdata.m() ) )
              {
                // all functions are real
                wftmpr_size = ft.np012();
                wftmpr_block_size = ft.np012loc(0);
              }
              else
              {
                // there is either 1) a mix of real and complex functions
                // or 2) only complex functions
                wftmpr_size = 2*ft.np012();
                wftmpr_block_size = 2*ft.np012loc(0);
              }
              DoubleMatrix wftmpr(sd->context(),wftmpr_size,sd->nst(),
                                  wftmpr_block_size,c.nb());
              
              wftmpr.getsub(gfdata,wftmpr_size,sd->nst(),0,
                            (ikp+s.wf.nkp())*sd->nst());
              
              vector<complex<double> > wftmp(ft.np012loc());
              for ( int nloc = 0; nloc < sd->nstloc(); nloc++ )
              {
                // copy column of wftmpr to complex array wftmp
                if ( wftmpr_size == ft.np012() )
                {
                  // function is real and must be expanded
                  double* p = wftmpr.valptr(wftmpr.mloc()*nloc);
                  for ( int i = 0; i < ft.np012loc(); i++ )
                    wftmp[i] = p[i];
                }
                else
                {
                  // function is complex
                  double* p = wftmpr.valptr(wftmpr.mloc()*nloc);
                  for ( int i = 0; i < ft.np012loc(); i++ )
                    wftmp[i] = complex<double>(p[2*i],p[2*i+1]);
                }
                ft.forward(&wftmp[0],c.valptr(c.mloc()*nloc));
              }
            }
          }
        }
      }
      int nparkpts_wf = s.wf.nparallelkpts();
      wfvtmp->set_nparallelkpts(nparkpts_wf);
    }
  }
  // force consistency of unit cell
  // copy wavefunction domain on atomset unit_cell
  s.atoms.set_cell(s.wf.cell());

  // check if wavefunction_velocity element was read, if not, delete wfvtmp
  if ( s.wfv != 0 )
  {
    delete s.wfv;
    s.wfv = 0;
  }
  if ( read_wfv )
  {
    s.wfv = wfvtmp;
  }
#else
  // USE_XERCES was not defined
  if ( ctxt_.onpe0() )
  {
    cout << "  SampleReader: could not read (parser not defined)"
         << endl;
  }
#endif
  tm.stop();
  if ( ctxt_.onpe0() )
    cout << " SampleReader: read time: " << tm.real() << " s" << endl;
}
