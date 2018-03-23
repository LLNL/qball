////////////////////////////////////////////////////////////////////////////////  
// Copyright (c) 2013-2016, Lawrence Livermore National Security, LLC. 
// qb@ll:  Qbox at Lawrence Livermore
//
// This file is part of qb@ll.
//
// Produced at the Lawrence Livermore National Laboratory. 
// Written by Erik Draeger (draeger1@llnl.gov) and Francois Gygi (fgygi@ucdavis.edu)
// and Xavier Andrade (xavier@tddft.org).
//
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
// SpeciesReader.C:
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#include "Species.h"
#include <qball/SpeciesReader.h>
#include "Messages.h"
#include <cassert>
#include <string>
#include <iostream>
#include <vector>
using namespace std;

#include <algorithm>
#include <string>
#include <sstream>
#include <cstdio>
#include <sys/stat.h>
#include <pseudo/psml.hpp>
#include <pseudo/qso.hpp>
#include <pseudo/upf1.hpp>
#include <pseudo/upf2.hpp>
#include <pseudo/detect_format.hpp>

////////////////////////////////////////////////////////////////////////////////
SpeciesReader::SpeciesReader(const Context& ctxt) : ctxt_(ctxt) {}

////////////////////////////////////////////////////////////////////////////////
void SpeciesReader::readSpecies(Species& sp, const string uri){

  if(ctxt_.oncoutpe()){

    pseudopotential::format format = pseudopotential::detect_format(uri);

    if(format == pseudopotential::format::FILE_NOT_FOUND) Messages::fatal("cannot find pseudopotential file '" + uri + "'");
    if(format == pseudopotential::format::UNKNOWN) Messages::fatal("cannot determine format for pseudopotential file '" + uri + "'");
    
    sp.uri_ = uri;

    cout << "  <!-- SpeciesReader opening file " << uri << " -->" << endl;

    pseudopotential::base * pseudo = NULL;
    
    switch(format){
    case pseudopotential::format::QSO:
      cout << "  <!--   format: QSO -->" << endl;
      pseudo = new pseudopotential::qso(uri);
      break;
    case pseudopotential::format::UPF1:
      cout << "  <!--   format: UPF1 -->" << endl;
      pseudo = new pseudopotential::upf1(uri, /*uniform_grid = */ true);
      break;
    case pseudopotential::format::UPF2:
      cout << "  <!--   format: UPF2 -->" << endl;
      pseudo = new pseudopotential::upf2(uri, /*uniform_grid = */ true);
      break;
    case pseudopotential::format::PSML:
      cout << "  <!--   format: PSML -->" << endl;
      pseudo = new pseudopotential::psml(uri, /*uniform_grid = */ true);
      break;
    default:
      Messages::fatal("unsupported format for pseudopotential file '" + uri + "'");
    }

    cout << "  <!--   size:   " << pseudo->size() << " -->" << endl;

    fill_species(sp, *pseudo);

    delete pseudo;
    
  }

}

void  SpeciesReader::fill_species(Species& sp, pseudopotential::base & pseudo){

  sp.nlcc_ = false;
  sp.usoft_ = false;
  sp.oncv_ = false;

  switch(pseudo.type()) {
  case pseudopotential::type::ULTRASOFT :
    sp.usoft_ = true;
    cout << "  <!-- SpeciesReader::readSpecies: potential type:  ultrasoft -->" << endl;
    break;
  case pseudopotential::type::KLEINMAN_BYLANDER :
    sp.oncv_ = true;
    cout << "  <!-- SpeciesReader::readSpecies: potential type:  Kleinman-Bylander norm-conserving -->" << endl;
    break;
  case pseudopotential::type::SEMILOCAL :
    cout << "  <!-- SpeciesReader::readSpecies: potential type:  norm-conserving -->" << endl;
    break;
  }
    
  sp.description_ = pseudo.description();
  cout << "  <!-- SpeciesReader::readSpecies: read description " << sp.description_ << " -->" << endl;
    
  sp.symbol_ = pseudo.symbol();
  cout << "  <!-- SpeciesReader::readSpecies: read symbol " << sp.symbol_ << " -->" << endl;

  sp.atomic_number_ = pseudo.atomic_number();
  cout << "  <!-- SpeciesReader::readSpecies: read atomic_number " << sp.atomic_number_ << " -->" << endl;

  sp.mass_ = pseudo.mass();
  cout << "  <!-- SpeciesReader::readSpecies: read mass " << sp.mass_ << " -->" << endl;

  sp.zval_ = pseudo.valence_charge();
  cout << "  <!-- SpeciesReader::readSpecies: read valence_charge " << sp.zval_ << " -->" << endl;

  sp.lmax_ = pseudo.lmax();
  cout << "  <!-- SpeciesReader::readSpecies: read lmax " << sp.lmax_ << " -->" << endl;

  sp.llocal_ = pseudo.llocal();
  cout << "  <!-- SpeciesReader::readSpecies: read llocal " << sp.llocal_ << " -->" << endl;
  
  if(pseudo.type() == pseudopotential::type::SEMILOCAL) { 
    sp.nquad_ = pseudo.nquad();
    cout << "  <!-- SpeciesReader::readSpecies: read nquad " << sp.nquad_ << " -->" << endl;
    sp.rquad_ = pseudo.rquad();
    cout << "  <!-- SpeciesReader::readSpecies: read rquad " << sp.rquad_ << " -->" << endl;
  } else {
    sp.nquad_ = 0;
    sp.rquad_ = 0.0;
  }

  sp.deltar_ = pseudo.mesh_spacing();
  cout << "  <!-- SpeciesReader::readSpecies: read mesh_spacing " << sp.deltar_ << " -->" << endl;
  
  sp.nchannels_ = pseudo.nchannels();

  // read the local potential
  if(pseudo.type() == pseudopotential::type::KLEINMAN_BYLANDER){
    pseudo.local_potential(sp.vloc_);
    cout << "  <!-- SpeciesReader::readSpecies: read local_potential size=" << sp.vloc_.size() << " -->" << endl;
      
    sp.projectors_.resize(sp.lmax_ + 1);

    for(int l = 0; l < sp.lmax_ + 1; l++ ) {
      sp.projectors_[l].resize(sp.nchannels_);
	  
      for ( int i = 0; i < sp.nchannels_; i++){
	pseudo.projector(l, i, sp.projectors_[l][i]);
	cout << "  <!-- SpeciesReader::readSpecies: read projector l="
	     << l << " i=" << i << " size=" << sp.projectors_[l][i].size() << " -->" << endl;
      }
    }

    /*
    std::ofstream pseudo_debug("projectors.dat");
    for(unsigned ip = 0; ip < sp.projectors_[0][0].size(); ip++){
      pseudo_debug << ip*sp.deltar_;
      for(int l = 0; l < sp.lmax_ + 1; l++ ) {
	for ( int i = 0; i < sp.nchannels_; i++){
	  pseudo_debug << '\t' << sp.projectors_[l][i][ip];
	}
      }
      pseudo_debug << endl;
    }
    pseudo_debug.close();
    */

    // the oncv weights
    sp.dij_.resize(sp.lmax_ + 1);
    for(int l = 0; l < sp.lmax_ + 1; l++){
      sp.dij_[l].resize(sp.nchannels_);
      for(int ii = 0; ii < sp.nchannels_; ii++){
	sp.dij_[l][ii].resize(sp.nchannels_);
	for(int jj = 0; jj < sp.nchannels_; jj++){
	  sp.dij_[l][ii][jj] = pseudo.d_ij(l, ii, jj);
	      
	  if(fabs(sp.dij_[l][ii][jj]) > 1.0e-6 && ii != jj){
	    cerr << endl << "Error:" << endl;
	    cerr << "       Unsupported pseudopotential file (non-diagonal ONCV)." << endl << endl;
	    exit(1);
	  }
	}
      }
      cout << "  <!-- SpeciesReader::readSpecies: read d_ij l=" << l << " -->" << endl;
    }
  }

  if(pseudo.type() == pseudopotential::type::SEMILOCAL){
    sp.vps_.resize(sp.lmax_ + 1);
    sp.phi_.resize(sp.lmax_ + 1);

    for(int l = 0; l < sp.lmax_ + 1; l++){

      pseudo.radial_potential(l, sp.vps_[l]);

      cout << "  <!-- SpeciesReader::readSpecies: read radial_potential l="
	   << l << " size=" << sp.vps_[l].size() << " -->" << endl;

      if(pseudo.has_radial_function(l)){
	pseudo.radial_function(l, sp.phi_[l]);

	cout << "  <!-- SpeciesReader::readSpecies: read radial_function  l="
	     << l << " size=" << sp.phi_[l].size() << " -->" << endl;
      }
    }
  }
    
  if(pseudo.type() == pseudopotential::type::ULTRASOFT){

    sp.nbeta_ = pseudo.nprojectors();
    cout << "  <!-- SpeciesReader::readSpecies: read nbeta " << sp.nbeta_ << " -->" << endl;
      
    sp.llocal_ = 0;
    sp.vps_.resize(1);
    pseudo.local_potential(sp.vps_[0]);
    cout << "  <!-- SpeciesReader::readSpecies: read vlocal size=" << sp.vps_[0].size() << " -->" << endl;

    sp.betar_.resize(sp.nbeta_);
    sp.betal_.resize(sp.nbeta_);
    for ( int b = 0; b < sp.nbeta_; b++ ){
      pseudo.beta(b, sp.betal_[b], sp.betar_[b]);
      cout << "  <!-- SpeciesReader::readSpecies: read beta b=" << b << " size=" << sp.betar_[b].size() << " -->" << endl;
    }

    // read dnm_zero (D_nm^0)
    pseudo.dnm_zero(sp.nbeta_, sp.dzero_);
    cout << "  <!-- SpeciesReader::readSpecies: read dnm_zero -->" << endl;

    // read rinner:  radii inside which Qnm^L are replaced by polynomial expansions
    //   (see Laasonen, Phys. Rev. B 47, 10142 (1993).)
    if(pseudo.has_rinner()){
      pseudo.rinner(sp.rinner_);
      cout << "  <!-- SpeciesReader::readSpecies: read rinner size=" << sp.rinner_.size() << " -->" << endl;
    }
      
    // read Q_nm and qfcoeff
    int nqnm = 0;
    for (int i = 1; i <= sp.nbeta_; i++) nqnm += i;
    sp.nqfun_ = nqnm;
    sp.qfunr_.resize(nqnm);
    sp.qfunl1_.resize(nqnm);
    sp.qfunl2_.resize(nqnm);
    sp.qfunb1_.resize(nqnm);
    sp.qfunb2_.resize(nqnm);
    sp.qfcoeff_.resize(nqnm);
      
    for(int q = 0; q < nqnm; q++){
      // qnm
      pseudo.qnm(q, sp.qfunl1_[q], sp.qfunl2_[q], sp.qfunb1_[q], sp.qfunb2_[q], sp.qfunr_[q]);
      cout << "  <!-- SpeciesReader::readSpecies: read qnm q="
	   << q << " l1=" << sp.qfunl1_[q] << " l2=" << sp.qfunl2_[q] << " size="
	   << sp.qfunr_[q].size() << " -->" << endl;
 
      // qcoeff
      if (sp.rinner_.size() > 0) {
	sp.qfcoeff_[q].resize(2*sp.lmax_+1);
	for (int ltot=0; ltot<2*sp.lmax_+1; ltot++) {
	  pseudo.qfcoeff(q, ltot, sp.qfcoeff_[q][ltot]);
	  cout << "  <!-- SpeciesReader::readSpecies: read qfcoeff q="
	       << q << " ltot=" << ltot << " size=" << sp.qfcoeff_[q][ltot].size() << " -->" << endl;
	}
      }
    }

  }

  if(pseudo.has_nlcc()){
    sp.nlcc_ = true;
    cout << "  <!-- SpeciesReader::readSpecies: nlcc found. -->" << endl;
    pseudo.nlcc_density(sp.rhor_nlcc_);
    cout << "  <!-- SpeciesReader::readSpecies: read rho_nlcc size=" << sp.rhor_nlcc_.size() << " -->" << endl;
  }

}

////////////////////////////////////////////////////////////////////////////////
void SpeciesReader::bcastSpecies(Species& sp)
{
  //cout << ctxt_.mype() << ": starting bcastSpecies" << endl;
  if ( ctxt_.oncoutpe() )
  {
    ctxt_.ibcast_send(1,1,&sp.atomic_number_,1);
    ctxt_.dbcast_send(1,1,&sp.mass_,1);
    ctxt_.ibcast_send(1,1,&sp.zval_,1);
    ctxt_.ibcast_send(1,1,&sp.lmax_,1);
    ctxt_.ibcast_send(1,1,&sp.llocal_,1);
    ctxt_.dbcast_send(1,1,&sp.deltar_,1);
    ctxt_.ibcast_send(1,1,&sp.nchannels_,1);
    int iusoft = (sp.ultrasoft() ? 1 : 0 );
    ctxt_.ibcast_send(1,1,&iusoft,1);
    int oncv = (sp.oncv_ ? 1:0);
    ctxt_.ibcast_send(1,1,&oncv,1);
    if (!sp.ultrasoft()) {
      ctxt_.ibcast_send(1,1,&sp.nquad_,1);
      ctxt_.dbcast_send(1,1,&sp.rquad_,1);
    }
    else {
      ctxt_.ibcast_send(1,1,&sp.nbeta_,1);
      ctxt_.ibcast_send(1,1,&sp.nqfun_,1);
    }
    int inlcc = (sp.nlcc() ? 1 : 0 );
    ctxt_.ibcast_send(1,1,&inlcc,1);
  }
  else
  {
    // calculate row and col indices of process oncoutpe
    int irow = ctxt_.coutpe();
    while (irow >= ctxt_.nprow())
      irow -= ctxt_.nprow();
    int icol = int (ctxt_.coutpe()/ctxt_.nprow());
    assert(ctxt_.pmap(irow,icol) == ctxt_.coutpe());
            
    ctxt_.ibcast_recv(1,1,&sp.atomic_number_,1,irow,icol);
    ctxt_.dbcast_recv(1,1,&sp.mass_,1,irow,icol);
    ctxt_.ibcast_recv(1,1,&sp.zval_,1,irow,icol);
    ctxt_.ibcast_recv(1,1,&sp.lmax_,1,irow,icol);
    ctxt_.ibcast_recv(1,1,&sp.llocal_,1,irow,icol);
    ctxt_.dbcast_recv(1,1,&sp.deltar_,1,irow,icol);
    ctxt_.ibcast_recv(1,1,&sp.nchannels_,1,irow,icol);

    int iusoft;
    ctxt_.ibcast_recv(1,1,&iusoft,1,irow,icol);
    sp.usoft_ = ( iusoft == 1 ? true : false );

    int oncv;
    ctxt_.ibcast_recv(1,1,&oncv,1,irow,icol);
    sp.oncv_ = ( oncv == 1);
    
    if (!sp.ultrasoft()) {
      ctxt_.ibcast_recv(1,1,&sp.nquad_,1,irow,icol);
      ctxt_.dbcast_recv(1,1,&sp.rquad_,1,irow,icol);
      sp.vps_.resize(sp.lmax_+1);
      sp.phi_.resize(sp.lmax_+1);
    }
    else {
      ctxt_.ibcast_recv(1,1,&sp.nbeta_,1,irow,icol);
      ctxt_.ibcast_recv(1,1,&sp.nqfun_,1,irow,icol);
      sp.vps_.resize(1);
      sp.betar_.resize(sp.nbeta_);
      sp.betal_.resize(sp.nbeta_);
      sp.dzero_.resize(sp.nbeta_);
      sp.qfunr_.resize(sp.nqfun_);
      sp.qfunl1_.resize(sp.nqfun_);
      sp.qfunl2_.resize(sp.nqfun_);
      sp.qfunb1_.resize(sp.nqfun_);
      sp.qfunb2_.resize(sp.nqfun_);
      sp.qfcoeff_.resize(sp.nqfun_);
    }
    int inlcc;
    ctxt_.ibcast_recv(1,1,&inlcc,1,irow,icol);
    sp.nlcc_ = ( inlcc == 1 ? true : false );
  }

  ctxt_.string_bcast(sp.symbol_,ctxt_.coutpe());
  ctxt_.string_bcast(sp.description_,ctxt_.coutpe());
  ctxt_.string_bcast(sp.uri_,ctxt_.coutpe());

  if (sp.oncv_) {
    // calculate row and col indices of process oncoutpe
    int irow = ctxt_.coutpe();
    while (irow >= ctxt_.nprow()) irow -= ctxt_.nprow();
    int icol = int (ctxt_.coutpe()/ctxt_.nprow());
    assert(ctxt_.pmap(irow,icol) == ctxt_.coutpe());

    int np;

    // the local potential
    if ( ctxt_.oncoutpe() ) {
      np = sp.vloc_.size();
      ctxt_.ibcast_send(1, 1, &np, 1);
      ctxt_.dbcast_send(np, 1, &sp.vloc_[0], np);
    } else {
      ctxt_.ibcast_recv(1, 1, &np, 1, irow, icol);
      sp.vloc_.resize(np);
      ctxt_.dbcast_recv(np, 1, &sp.vloc_[0], np, irow, icol);
    }

    if(!ctxt_.oncoutpe()){
      sp.projectors_.resize(sp.lmax_ + 1);
      sp.dij_.resize(sp.lmax_ + 1);
    }
    
    for(int ll = 0; ll <= sp.lmax_; ll++ ){

      if(!ctxt_.oncoutpe()){
	sp.projectors_[ll].resize(sp.nchannels_);
	sp.dij_[ll].resize(sp.nchannels_);
      }
    
      for(int ii = 0; ii < sp.nchannels_; ii++){

	// the projectors
	if ( ctxt_.oncoutpe() ) {
	  np = sp.projectors_[ll][ii].size();
	  ctxt_.ibcast_send(1, 1, &np, 1);
	  ctxt_.dbcast_send(np, 1, &sp.projectors_[ll][ii][0], np);
	} else {
	  ctxt_.ibcast_recv(1, 1, &np, 1, irow, icol);
	  sp.projectors_[ll][ii].resize(np);
	  ctxt_.dbcast_recv(np, 1, &sp.projectors_[ll][ii][0], np, irow, icol);
	}

	// the weights
	if ( ctxt_.oncoutpe() ) {
	  ctxt_.dbcast_send(sp.nchannels_, 1, &sp.dij_[ll][ii][0], sp.nchannels_);
	} else {
	  sp.dij_[ll][ii].resize(sp.nchannels_);
	  ctxt_.dbcast_recv(sp.nchannels_, 1, &sp.dij_[ll][ii][0], sp.nchannels_, irow, icol);
	}

      }
    }
    
  } else if (!sp.ultrasoft()) {

    for ( int l = 0; l <= sp.lmax_; l++ )
    {
      int np_vps;
      if ( ctxt_.oncoutpe() )
      {
        np_vps = sp.vps_[l].size();
        ctxt_.ibcast_send(1,1,&np_vps,1);
        ctxt_.dbcast_send(np_vps,1,&sp.vps_[l][0],np_vps);
      }
      else
      {
        // calculate row and col indices of process oncoutpe
        int irow = ctxt_.coutpe();
        while (irow >= ctxt_.nprow())
          irow -= ctxt_.nprow();
        int icol = int (ctxt_.coutpe()/ctxt_.nprow());
        assert(ctxt_.pmap(irow,icol) == ctxt_.coutpe());
            
        ctxt_.ibcast_recv(1,1,&np_vps,1,irow,icol);
        sp.vps_[l].resize(np_vps);
        ctxt_.dbcast_recv(np_vps,1,&sp.vps_[l][0],np_vps,irow,icol);
      }
    }

    // broadcast atomic orbitals
    for ( int l = 0; l <= sp.lmax_; l++ )
    {
      int np_phi;
      if ( ctxt_.oncoutpe() )
      {
        np_phi = sp.phi_[l].size();
        ctxt_.ibcast_send(1,1,&np_phi,1);
        ctxt_.dbcast_send(np_phi,1,&sp.phi_[l][0],np_phi);
      }
      else
      {
        // calculate row and col indices of process oncoutpe
        int irow = ctxt_.coutpe();
        while (irow >= ctxt_.nprow())
          irow -= ctxt_.nprow();
        int icol = int (ctxt_.coutpe()/ctxt_.nprow());
        assert(ctxt_.pmap(irow,icol) == ctxt_.coutpe());
      
        ctxt_.ibcast_recv(1,1,&np_phi,1,irow,icol);
        sp.phi_[l].resize(np_phi);
        ctxt_.dbcast_recv(np_phi,1,&sp.phi_[l][0],np_phi,irow,icol);
      }
    }

    // broadcast nlcc density, if any
    int np_nlcc;
    if ( ctxt_.oncoutpe() )
    {
      np_nlcc = sp.rhor_nlcc_.size();
      ctxt_.ibcast_send(1,1,&np_nlcc,1);
      if (np_nlcc > 0)
         ctxt_.dbcast_send(np_nlcc,1,&sp.rhor_nlcc_[0],np_nlcc);
    }
    else
    {
      // calculate row and col indices of process oncoutpe
      int irow = ctxt_.coutpe();
      while (irow >= ctxt_.nprow())
        irow -= ctxt_.nprow();
      int icol = int (ctxt_.coutpe()/ctxt_.nprow());
      assert(ctxt_.pmap(irow,icol) == ctxt_.coutpe());
      
      ctxt_.ibcast_recv(1,1,&np_nlcc,1,irow,icol);
         
      if (np_nlcc > 0) {
         sp.rhor_nlcc_.resize(np_nlcc);
         ctxt_.dbcast_recv(np_nlcc,1,&sp.rhor_nlcc_[0],np_nlcc,irow,icol);         
      }
    }

  }
  else {  // ultrasoft

    // broadcast local potential, stored on coutpe in vps_[0]
    int np_vps;
    if ( ctxt_.oncoutpe() )
    {
      np_vps = sp.vps_[0].size();
      ctxt_.ibcast_send(1,1,&np_vps,1);
      ctxt_.dbcast_send(np_vps,1,&sp.vps_[0][0],np_vps);
    }
    else
    {
      // calculate row and col indices of process oncoutpe
      int irow = ctxt_.coutpe();
      while (irow >= ctxt_.nprow())
        irow -= ctxt_.nprow();
      int icol = int (ctxt_.coutpe()/ctxt_.nprow());
      assert(ctxt_.pmap(irow,icol) == ctxt_.coutpe());
      
      ctxt_.ibcast_recv(1,1,&np_vps,1,irow,icol);
      sp.vps_[0].resize(np_vps);
      ctxt_.dbcast_recv(np_vps,1,&sp.vps_[0][0],np_vps,irow,icol);
    }

      // broadcast beta fns
    for ( int b = 0; b < sp.nbeta_; b++ )
    {
      int np_beta;
      if ( ctxt_.oncoutpe() )
      {
        np_beta = sp.betar_[b].size();
        ctxt_.ibcast_send(1,1,&np_beta,1);
        ctxt_.dbcast_send(np_beta,1,&sp.betar_[b][0],np_beta);
        ctxt_.ibcast_send(1,1,&sp.betal_[b],1);
      }
      else
      {
        // calculate row and col indices of process oncoutpe
        int irow = ctxt_.coutpe();
        while (irow >= ctxt_.nprow())
          irow -= ctxt_.nprow();
        int icol = int (ctxt_.coutpe()/ctxt_.nprow());
        assert(ctxt_.pmap(irow,icol) == ctxt_.coutpe());            
        ctxt_.ibcast_recv(1,1,&np_beta,1,irow,icol);
        sp.betar_[b].resize(np_beta);
        ctxt_.dbcast_recv(np_beta,1,&sp.betar_[b][0],np_beta,irow,icol);
        ctxt_.ibcast_recv(1,1,&sp.betal_[b],1,irow,icol);
      }
    }

    // broadcast dzero
    for ( int b = 0; b < sp.nbeta_; b++ )
    {
      int np_dzero;
      if ( ctxt_.oncoutpe() )
      {
        np_dzero = sp.dzero_[b].size();
        ctxt_.ibcast_send(1,1,&np_dzero,1);
        ctxt_.dbcast_send(np_dzero,1,&sp.dzero_[b][0],np_dzero);
      }
      else
      {
        // calculate row and col indices of process oncoutpe
        int irow = ctxt_.coutpe();
        while (irow >= ctxt_.nprow())
          irow -= ctxt_.nprow();
        int icol = int (ctxt_.coutpe()/ctxt_.nprow());
        assert(ctxt_.pmap(irow,icol) == ctxt_.coutpe());            
        ctxt_.ibcast_recv(1,1,&np_dzero,1,irow,icol);
        sp.dzero_[b].resize(np_dzero);
        ctxt_.dbcast_recv(np_dzero,1,&sp.dzero_[b][0],np_dzero,irow,icol);
      }
    }

    // broadcast rinner
    int np_rinner;
    if ( ctxt_.oncoutpe() )
    {
      np_rinner = sp.rinner_.size();
      ctxt_.ibcast_send(1,1,&np_rinner,1);
      if (np_rinner > 0) 
        ctxt_.dbcast_send(np_rinner,1,&sp.rinner_[0],np_rinner);
    }
    else
    {
      // calculate row and col indices of process oncoutpe
      int irow = ctxt_.coutpe();
      while (irow >= ctxt_.nprow())
        irow -= ctxt_.nprow();
      int icol = int (ctxt_.coutpe()/ctxt_.nprow());
      assert(ctxt_.pmap(irow,icol) == ctxt_.coutpe());            
      ctxt_.ibcast_recv(1,1,&np_rinner,1,irow,icol);
      if (np_rinner > 0) {
        sp.rinner_.resize(np_rinner);
        ctxt_.dbcast_recv(np_rinner,1,&sp.rinner_[0],np_rinner,irow,icol);
      }
    }

    // broadcast qfun, qfcoeff
    for ( int q = 0; q < sp.nqfun_; q++ )
    {
      int np_qfun, np_qfcoeff;
      if ( ctxt_.oncoutpe() )
      {
        np_qfun = sp.qfunr_[q].size();
        ctxt_.ibcast_send(1,1,&np_qfun,1);
        ctxt_.dbcast_send(np_qfun,1,&sp.qfunr_[q][0],np_qfun);
        ctxt_.ibcast_send(1,1,&sp.qfunl1_[q],1);
        ctxt_.ibcast_send(1,1,&sp.qfunl2_[q],1);
        ctxt_.ibcast_send(1,1,&sp.qfunb1_[q],1);
        ctxt_.ibcast_send(1,1,&sp.qfunb2_[q],1);
        int qfcosize = sp.qfcoeff_[q].size();
        ctxt_.ibcast_send(1,1,&qfcosize,1);
        if (qfcosize > 0) {
          for (int ltot=0; ltot<2*sp.lmax_+1; ltot++) {
            np_qfcoeff = sp.qfcoeff_[q][ltot].size();
            ctxt_.ibcast_send(1,1,&np_qfcoeff,1);
            ctxt_.dbcast_send(np_qfcoeff,1,&sp.qfcoeff_[q][ltot][0],np_qfcoeff);
          }
        }
      }
      else
      {
        // calculate row and col indices of process oncoutpe
        int irow = ctxt_.coutpe();
        while (irow >= ctxt_.nprow())
          irow -= ctxt_.nprow();
        int icol = int (ctxt_.coutpe()/ctxt_.nprow());
        assert(ctxt_.pmap(irow,icol) == ctxt_.coutpe());            
        ctxt_.ibcast_recv(1,1,&np_qfun,1,irow,icol);
        sp.qfunr_[q].resize(np_qfun);
        ctxt_.dbcast_recv(np_qfun,1,&sp.qfunr_[q][0],np_qfun,irow,icol);
        ctxt_.ibcast_recv(1,1,&sp.qfunl1_[q],1,irow,icol);
        ctxt_.ibcast_recv(1,1,&sp.qfunl2_[q],1,irow,icol);
        ctxt_.ibcast_recv(1,1,&sp.qfunb1_[q],1,irow,icol);
        ctxt_.ibcast_recv(1,1,&sp.qfunb2_[q],1,irow,icol);
        int qfcosize;
        ctxt_.ibcast_recv(1,1,&qfcosize,1,irow,icol);
        if (qfcosize > 0) { 
          assert(qfcosize == 2*sp.lmax_+1);
          sp.qfcoeff_[q].resize(2*sp.lmax_+1);
          for (int ltot=0; ltot<2*sp.lmax_+1; ltot++) {
            ctxt_.ibcast_recv(1,1,&np_qfcoeff,1,irow,icol);
            sp.qfcoeff_[q][ltot].resize(np_qfcoeff);
            ctxt_.dbcast_recv(np_qfcoeff,1,&sp.qfcoeff_[q][ltot][0],np_qfcoeff,irow,icol);
          }
        }
      }
    }

    int np_nlcc;
    if ( ctxt_.oncoutpe() )
    {
      np_nlcc = sp.rhor_nlcc_.size();
      ctxt_.ibcast_send(1,1,&np_nlcc,1);
      if (np_nlcc > 0)
         ctxt_.dbcast_send(np_nlcc,1,&sp.rhor_nlcc_[0],np_nlcc);
    }
    else
    {
      // calculate row and col indices of process oncoutpe
      int irow = ctxt_.coutpe();
      while (irow >= ctxt_.nprow())
        irow -= ctxt_.nprow();
      int icol = int (ctxt_.coutpe()/ctxt_.nprow());
      assert(ctxt_.pmap(irow,icol) == ctxt_.coutpe());
      
      ctxt_.ibcast_recv(1,1,&np_nlcc,1,irow,icol);
         
      if (np_nlcc > 0) {
         sp.rhor_nlcc_.resize(np_nlcc);
         ctxt_.dbcast_recv(np_nlcc,1,&sp.rhor_nlcc_[0],np_nlcc,irow,icol);         
      }
    }
  }
}
