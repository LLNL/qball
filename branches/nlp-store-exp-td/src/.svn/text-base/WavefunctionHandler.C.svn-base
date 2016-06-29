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
// WavefunctionHandler.C
//
////////////////////////////////////////////////////////////////////////////////

#if USE_XERCES

#include "WavefunctionHandler.h"
#include "Wavefunction.h"
#include "FourierTransform.h"
#include "Timer.h"
#include "SampleReader.h"

#include "StrX.h"
// XML transcoding for loading grid_functions
#include <xercesc/util/Base64.hpp>
#include <xercesc/util/XMLString.hpp>
using namespace xercesc;
#include <iostream>
#include <cassert>
#include <sstream>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
WavefunctionHandler::WavefunctionHandler(Wavefunction& wf,
  DoubleMatrix& gfdata, int& nx, int& ny, int& nz,
  vector<vector<vector<double> > > &dmat) :
  wf_(wf), gfdata_(gfdata), dmat_(dmat), ecut(0.0), nx_(nx), ny_(ny), nz_(nz)
{
  // if the data is read from a file, gfdata has a finite size
  // since the grid functions were processed by XMLGFPreprocessor
  // if the gfdata matrix has finite dimensions, set read_from_gfdata flag
  read_from_gfdata = ( gfdata.n() > 0 );
  current_igf = 0;
}

////////////////////////////////////////////////////////////////////////////////
WavefunctionHandler::~WavefunctionHandler() {}

////////////////////////////////////////////////////////////////////////////////
void WavefunctionHandler::byteswap_double(size_t n, double* x)
{
  if (n==0) return;
  unsigned char* c = (unsigned char*) x;
  while ( n-- > 0 )
  {
    unsigned char tmp;
    tmp = c[7]; c[7] = c[0]; c[0] = tmp;
    tmp = c[6]; c[6] = c[1]; c[1] = tmp;
    tmp = c[5]; c[5] = c[2]; c[2] = tmp;
    tmp = c[4]; c[4] = c[3]; c[3] = tmp;

    c+=8;
  }
}

////////////////////////////////////////////////////////////////////////////////
void WavefunctionHandler::startElement(const XMLCh* const uri,
  const XMLCh* const localname, const XMLCh* const qname,
  const Attributes& attributes)
{
  // cout << " WavefunctionHandler::startElement " << StrX(qname) << endl;
  string locname(XMLString::transcode(localname));

  int nspin=1, nel=0, nempty=0, delta_spin=0;
  event_type event = invalid;
  // consider only elements that are dealt with directly by WavefunctionHandler

  if ( locname == "wavefunction" || locname == "wavefunction_velocity" )
  {

    unsigned int len = attributes.getLength();
    for (unsigned int index = 0; index < len; index++)
    {
      string attrname(XMLString::transcode(attributes.getLocalName(index)));
      if ( attrname == "ecut")
      {
        ecut = atof(StrX(attributes.getValue(index)).localForm());
      }
      else if ( attrname == "nspin")
      {
        nspin = atoi(StrX(attributes.getValue(index)).localForm());
      }
      else if ( attrname == "delta_spin")
      {
        delta_spin = atoi(StrX(attributes.getValue(index)).localForm());
      }
      else if ( attrname == "nempty" )
      {
        nempty = atoi(StrX(attributes.getValue(index)).localForm());
      }
      else if ( attrname == "nel" )
      {
        nel = atoi(StrX(attributes.getValue(index)).localForm());
      }
    }

    cout << " WavefunctionHandler::startElement: " << locname
         << " nspin=" << nspin << " delta_spin=" << delta_spin << " nel=" << nel
         << " nempty=" << nempty << endl;

    current_ispin = 0;
    current_ikp = 0;
    current_n = 0;

    // notify listening nodes
    if ( locname == "wavefunction" )
      event = wavefunction;
    else
      event = wavefunction_velocity;

    wf_.context().ibcast_send(1,1,(int*)&event,1);
    wf_.context().ibcast_send(1,1,&nel,1);
    wf_.context().ibcast_send(1,1,&nspin,1);
    wf_.context().ibcast_send(1,1,&delta_spin,1);
    wf_.context().ibcast_send(1,1,&nempty,1);

    wf_.set_nel(nel);
    wf_.set_nspin(nspin);
    if (nspin > 1)
      wf_.set_deltaspin(delta_spin);
    wf_.set_nempty(nempty);
  }
  else if ( locname == "domain")
  {
    D3vector a,b,c;
    unsigned int len = attributes.getLength();
    for (unsigned int index = 0; index < len; index++)
    {
      string attrname(XMLString::transcode(attributes.getLocalName(index)));
      string attrval(XMLString::transcode(attributes.getValue(index)));
      istringstream stst(attrval);
      if ( attrname == "a")
      {
        stst >> a;
      }
      else if ( attrname == "b" )
      {
        stst >> b;
      }
      else if ( attrname == "c" )
      {
        stst >> c;
      }
    }

    //cout << " WavefunctionHandler::startElement: domain" << endl;
    uc.set(a,b,c);
    //cout << uc;

    // notify listening nodes
    double buf[9];
    buf[0] = uc.a(0).x; buf[1] = uc.a(0).y; buf[2] = uc.a(0).z;
    buf[3] = uc.a(1).x; buf[4] = uc.a(1).y; buf[5] = uc.a(1).z;
    buf[6] = uc.a(2).x; buf[7] = uc.a(2).y; buf[8] = uc.a(2).z;
    wf_.context().dbcast_send(9,1,buf,9);
  }
  else if ( locname == "reference_domain")
  {
    D3vector a,b,c;
    unsigned int len = attributes.getLength();
    for (unsigned int index = 0; index < len; index++)
    {
      string attrname(XMLString::transcode(attributes.getLocalName(index)));
      string attrval(XMLString::transcode(attributes.getValue(index)));
      istringstream stst(attrval);
      if ( attrname == "a")
      {
        stst >> a;
      }
      else if ( attrname == "b" )
      {
        stst >> b;
      }
      else if ( attrname == "c" )
      {
        stst >> c;
      }
    }

    //cout << " WavefunctionHandler::startElement: reference_domain" << endl;
    ruc.set(a,b,c);
    //cout << ruc;
  }
  else if ( locname == "density_matrix")
  {
    unsigned int len = attributes.getLength();
    int dmat_size = 0;
    string dmat_form;
    for (unsigned int index = 0; index < len; index++)
    {
      string attrname(XMLString::transcode(attributes.getLocalName(index)));
      string attrval(XMLString::transcode(attributes.getValue(index)));
      istringstream stst(attrval);
      if ( attrname == "form")
      {
        stst >> dmat_form;
      }
      else if ( attrname == "size" )
      {
        stst >> dmat_size;
      }
    }
    if ( dmat_form != "diagonal" )
    {
      cout << "WavefunctionHandler: density_matrix must be diagonal" << endl;
      wf_.context().abort(1);
    }
    dmat_tmp.resize(dmat_size);
  }
  else if ( locname == "grid")
  {
    unsigned int len = attributes.getLength();
    for (unsigned int index = 0; index < len; index++)
    {
      string attrname(XMLString::transcode(attributes.getLocalName(index)));
      string attrval(XMLString::transcode(attributes.getValue(index)));
      istringstream stst(attrval);
      if ( attrname == "nx")
      {
        stst >> nx_;
      }
      else if ( attrname == "ny" )
      {
        stst >> ny_;
      }
      else if ( attrname == "nz" )
      {
        stst >> nz_;
      }
    }

    // notify listening nodes
    int ibuf[3];
    ibuf[0] = nx_;
    ibuf[1] = ny_;
    ibuf[2] = nz_;
    wf_.context().ibcast_send(3,1,ibuf,3);

    if ( ecut == 0.0 )
    {
      // ecut attribute was not specified. Infer from grid size
      // Ecut = max(ecut_x,ecut_y,ecut_z)

      // When importing grids with Dirichlet b.c. grid sizes can be odd
      // round nx,ny,nz to next even number to compute ecut
      // use nx+nx%2 instead of nx
      double g0_max = ((2*(nx_+nx_%2)-1.0)/2.0) * 2.0 * M_PI / length(uc.a(0));
      double g1_max = ((2*(ny_+ny_%2)-1.0)/2.0) * 2.0 * M_PI / length(uc.a(1));
      double g2_max = ((2*(nz_+nz_%2)-1.0)/2.0) * 2.0 * M_PI / length(uc.a(2));
      double ecut0 = 0.125 * g0_max * g0_max;
      double ecut1 = 0.125 * g1_max * g1_max;
      double ecut2 = 0.125 * g2_max * g2_max;

      ecut = max(max(ecut0,ecut1),ecut2);
      cout << " ecut=" << 2*ecut << " Ry" << endl;
    }
    // notify listening nodes of ecut
    wf_.context().dbcast_send(1,1,&ecut,1);

    // notify listening nodes of the reference_domain
    // note: the reference_domain is optional in the sample file
    // notify listening nodes
    double buf[9];
    buf[0] = ruc.a(0).x; buf[1] = ruc.a(0).y; buf[2] = ruc.a(0).z;
    buf[3] = ruc.a(1).x; buf[4] = ruc.a(1).y; buf[5] = ruc.a(1).z;
    buf[6] = ruc.a(2).x; buf[7] = ruc.a(2).y; buf[8] = ruc.a(2).z;
    wf_.context().dbcast_send(9,1,buf,9);

    wf_.resize(uc,ruc,ecut);
  }
  else if ( locname == "slater_determinant")
  {
    unsigned int len = attributes.getLength();
    for (unsigned int index = 0; index < len; index++)
    {
      string attrname(XMLString::transcode(attributes.getLocalName(index)));
      string attrval(XMLString::transcode(attributes.getValue(index)));
      istringstream stst(attrval);
      if ( attrname == "kpoint")
      {
        stst >> current_kx >> current_ky >> current_kz;
        cout << " read slater_determinant kpoint=" << current_kx
             << " " << current_ky << " " << current_kz;
      }
      else if ( attrname == "weight" )
      {
        stst >> current_weight;
        cout << " weight=" << current_weight;
      }
      else if ( attrname == "size" )
      {
        stst >> current_size;
        cout << " size=" << current_size << endl;
      }
    }
    // notify listening nodes: slater_determinant
    event = slater_determinant;
    wf_.context().ibcast_send(1,1,(int*)&event,1);

    // send current spin and kpoint index to listening nodes
    int tind[2];
    tind[0] = current_ispin;  tind[1] = current_ikp;
    wf_.context().ibcast_send(2,1,tind,2);
    
    // send kpoint and weight to listening nodes
    double buf[4];
    buf[0] = current_kx; buf[1] = current_ky; buf[2] = current_kz;
    buf[3] = current_weight;
    wf_.context().dbcast_send(4,1,buf,4);

    //ewd need to avoid adding each k-point twice when nspin = 2
    if (current_ispin == 0)
      wf_.add_kpoint(D3vector(current_kx,current_ky,current_kz),current_weight);

    if (wf_.kptactive(current_ikp) && wf_.spinactive(current_ispin)) {
      assert(current_size==wf_.sd(current_ispin,current_ikp)->nst());
      const Basis& basis = wf_.sd(current_ispin,current_ikp)->basis();
      ft = new FourierTransform(basis,nx_,ny_,nz_);
      wftmp.resize((ft->np012loc()));
    }
  }
  else if ( locname == "grid_function")
  {
    unsigned int len = attributes.getLength();
    for (unsigned int index = 0; index < len; index++)
    {
      string attrname(XMLString::transcode(attributes.getLocalName(index)));
      string attrval(XMLString::transcode(attributes.getValue(index)));
      istringstream stst(attrval);
      if ( attrname == "nx")
      {
        stst >> current_gf_nx;
      }
      else if ( attrname == "ny" )
      {
        stst >> current_gf_ny;
      }
      else if ( attrname == "nz" )
      {
        stst >> current_gf_nz;
      }
      else if ( attrname == "encoding" )
      {
        stst >> current_gf_encoding;
      }
    }

    //cout << " WavefunctionHandler::startElement: grid_function"
    //     << " nx=" << current_gf_nx
    //     << " ny=" << current_gf_ny
    //     << " nz=" << current_gf_nz
    //     << "\n encoding=" << current_gf_encoding
    //     << endl;
  }
}

////////////////////////////////////////////////////////////////////////////////
void WavefunctionHandler::endElement(const XMLCh* const uri,
  const XMLCh* const localname, const XMLCh* const qname, string& content)
{
  string locname(XMLString::transcode(localname));
  //cout << " WavefunctionHandler::endElement " << locname << endl;
  if ( locname == "density_matrix")
  {
    istringstream stst(content);
    for ( int i = 0; i < dmat_tmp.size(); i++ )
      stst >> dmat_tmp[i];

    // send dmat to listening nodes
    dmat_[current_ispin].push_back(dmat_tmp);
    int tsize = dmat_tmp.size();
    wf_.context().ibcast_send(1,1,&tsize,1);
    wf_.context().dbcast_send(dmat_tmp.size(),1,&dmat_tmp[0],dmat_tmp.size());
  }
  else if ( locname == "grid_function")
  {
    // current implementation accepts only full grids as declared in
    // the wavefunction <grid> element
    assert(current_gf_nx==nx_ &&
           current_gf_ny==ny_ &&
           current_gf_nz==nz_ );

    if ( read_from_gfdata ) {
      // do nothing
      //cout << "WavefunctionHandler::endElement: current_igf=" << current_igf
      //     << endl;
      current_igf++;
    }
    else {
      const Basis& basis = wf_.sd(current_ispin,current_ikp)->basis();
      int wftmpr_size = current_gf_nx*current_gf_ny*current_gf_nz;
      // if basis is complex, double the size of wftmpr
      if ( !basis.real() )
        wftmpr_size *= 2;

      valarray<double> wftmpr(wftmpr_size);

      if ( current_gf_encoding == "text" ) {
        istringstream stst(content);
        for ( int i = 0; i < wftmpr_size; i++ )
          stst >> wftmpr[i];
      }
      else if ( current_gf_encoding == "base64" ) {
        // base64 encoding
        unsigned int length;

#ifdef XERCESC_3
        XMLByte* b = Base64::decode((XMLByte*)content.c_str(),
                                    (XMLSize_t*) &length);
#else
        XMLByte* b = Base64::decode((XMLByte*)content.c_str(), &length);
#endif        
        assert(b!=0);
        // use data in b
        assert(length/sizeof(double)==wftmpr_size);
#if PLT_BIG_ENDIAN
        byteswap_double(wftmpr_size,(double*)b);
#endif
        memcpy(&wftmpr[0],b,wftmpr_size*sizeof(double));
        XMLString::release((char**)&b);
      }
      else
      {
        cout << "WavefunctionHandler: unknown grid_function encoding" << endl;
        return;
      }

      // wftmpr now contains the grid function

      // send subgrids to listening nodes

      if (wf_.kptactive(current_ikp) && wf_.spinactive(current_ispin))
      {
        SlaterDet* sd = wf_.sd(current_ispin,current_ikp);
        assert(sd != 0);
        // pcol = process column destination
        int pcol = sd->c().pc(current_n);
        for ( int prow = 0; prow < sd->context().nprow(); prow++ )
        {
          int size = ft->np2_loc(prow) * ft->np0() * ft->np1();
          int istart = ft->np2_first(prow) * ft->np0() * ft->np1();
          if ( !basis.real() )
          {
            size *= 2;
            istart *= 2;
          }
          // send subgrid to node (prow,pcol)
          if ( !(prow==0 && pcol==0) )
          {
            sd->context().isend(1,1,&size,1,prow,pcol);
            sd->context().dsend(size,1,&wftmpr[istart],size,prow,pcol);
          }
        }

        // if destination column is pcol=0, copy to complex array on node 0
        // and process
        if ( pcol == 0 )
        {
          if ( basis.real() )
          {
            for ( int i = 0; i < ft->np012loc(); i++ )
              wftmp[i] = wftmpr[i];
          }
          else
          {
            for ( int i = 0; i < ft->np012loc(); i++ )
              wftmp[i] = complex<double>(wftmpr[2*i],wftmpr[2*i+1]);
          }
          ComplexMatrix& c = wf_.sd(current_ispin,current_ikp)->c();
          ft->forward(&wftmp[0],c.valptr(c.mloc()*current_n));
        }
      }
    }
    //cout << " grid_function read n=" << current_n << endl;
    current_n++;
  }
  else if ( locname == "slater_determinant")
  {
    int nspin = wf_.nspin();
    if (nspin == 1) 
    {
      current_ikp++;
    }
    else
    {
      if (current_ispin == 0)
      {
        current_ispin = 1;
      }
      else if (current_ispin == 1) 
      {
        current_ispin = 0;
        current_ikp++;
      }
    }
    delete ft;
  }
}

////////////////////////////////////////////////////////////////////////////////
StructureHandler* WavefunctionHandler::startSubHandler(const XMLCh* const uri,
    const XMLCh* const localname, const XMLCh* const qname,
    const Attributes& attributes)
{
  // check if element qname can be processed by another StructureHandler
  // If it can, return a pointer to the StructureHandler, otherwise return 0
  // cout << " WavefunctionHandler::startSubHandler " << StrX(qname) << endl;
  return 0;
}

////////////////////////////////////////////////////////////////////////////////
void WavefunctionHandler::endSubHandler(const XMLCh* const uri,
    const XMLCh* const localname, const XMLCh* const qname,
    const StructureHandler* const last)
{
  string locname(XMLString::transcode(localname));
  //cout << " WavefunctionHandler::endSubHandler " << locname << endl;
  delete last;
}

#endif
