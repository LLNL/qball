////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2014 The Regents of the University of California
//
// This file is part of Qbox
//
// Qbox is distributed under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 2 of
// the License, or (at your option) any later version.
// See the file COPYING in the root directory of this distribution
// or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
// CurrentDensity.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: CurrentDensity.h,v 1.13 2008-09-08 15:56:18 fgygi Exp $

#include <vector>
#include <valarray>
#include "CurrentDensity.h"
#include <math/d3vector.h>
#include "Basis.h"
#include <qball/UnitCell.h>
#include "SlaterDet.h"
#include <fstream>
#include "isodate.h"
#include "release.h"
#include "Species.h"
#include <iomanip>
#include "Base64Transcoder.h"

CurrentDensity::CurrentDensity(const Sample& s, const Wavefunction & wf):
  ChargeDensity(s), wf_(wf){
  
}

void CurrentDensity::update_current(EnergyFunctional & energy_functional, const Wavefunction & dwf,bool output){

  Wavefunction rwf(wf_);
  Wavefunction rdwf(wf_);
  Wavefunction drwf(wf_);

  std::vector<std::vector<double> > fion;
  std::valarray<double> sigma;
  std::vector<std::complex<double> > tmp(vft()->np012loc(), 0.0);

  current.resize(3);
  double volume_element = vbasis()->cell().volume()/vft()->np012();

  for(int idir = 0; idir < 3; idir++){

    current[idir].resize(wf_.nspin());
    total_current[idir] = 0.0;

    for ( int ispin = 0; ispin < wf_.nspin(); ispin++){
      current[idir][ispin].resize(vft()->np012loc());

      for(int ip = 0; ip < vft()->np012loc(); ip++){
	current[idir][ispin][ip] = 0.0;
      }

      for ( int ikp = 0; ikp < wf_.nkp(); ikp++ ){

	const int ngwloc = wf_.sd(ispin, ikp)->basis().localsize();
	const int mloc = wf_.sd(ispin, ikp)->c().mloc();

	for ( int n = 0; n < wf_.sd(ispin, ikp)->nstloc(); n++ ){
	  for ( int ig = 0; ig < ngwloc; ig++ ){
	    rwf.sd(ispin, ikp)->c()[ig + mloc*n] = std::complex<double>(0.0, 1.0)*wf_.sd(ispin, ikp)->basis().kpgx_ptr(idir)[ig]*wf_.sd(ispin, ikp)->c()[ig + mloc*n];
	  }
	}

	for(int ip = 0; ip < vft()->np012loc(); ip++){
	  tmp[ip] = 0.0;
	}

	wf_.sd(ispin, ikp)->compute_density(*ft(ispin, ikp), wf_.weight(ikp), &tmp[0], *rwf.sd(ispin, ikp));
	
	for(int ip = 0; ip < vft()->np012loc(); ip++){
	  current[idir][ispin][ip] += -std::imag(tmp[ip]);
	}
      }

      wf_.wfcontext()->dsum('r', vft()->np012loc(), 1, &current[idir][ispin][0], vft()->np012loc());
      
      for(int ip = 0; ip < vft()->np012loc(); ip++){
	total_current[idir] += volume_element*current[idir][ispin][ip];
      }

      wf_.spincontext(ispin)->dsum('c', 1, 1, &total_current[idir], 1);

      if(energy_functional.vp) total_current[idir] += energy_functional.vp->value()[idir]*energy_functional.hamil_cd()->nelectrons();
      
    }
    
  }

  // TODO: Reduce total current over spin
  assert(wf_.nspin() == 1);
	 
  if (output){
    if ( wf_.context().onpe0() ){
        std::cout << "  total_electronic_current:\t" << std::fixed << std::setw( 20 ) << std::setprecision( 12 ) << total_current[0] << '\t' << total_current[1] << '\t' << total_current[2] << std::endl;
    }
  }
}

void CurrentDensity::plot(const Sample * s, const std::string & filename){
  using namespace std;

  std::vector<std::vector<double> > global_current(3);

  for(int idir = 0; idir < 3; idir++) {
    vft()->gather(*s->wf.spincontext(0), current[idir][0], global_current[idir]);
  }

  if ( s->ctxt_.onpe0() )
    {

      const int np0 = vft()->np0();
      const int np1 = vft()->np1();
      const int np2 = vft()->np2();
      
      for(int idir = 0; idir < 3; idir++){

	std::ofstream os;

	switch(idir){
	case 0:
	  os.open(("x" + filename + ".cube").c_str());
	  break;
	case 1:
	  os.open(("y" + filename + ".cube").c_str());
	  break;
	case 2:
	  os.open(("z" + filename + ".cube").c_str());
	  break;
	}

	// write header and atoms
	os << "Created " << isodate() << " by qbox-" << release() << endl;
	os << endl;
      
	int natoms = s->atoms.size();
	D3vector a0 = s->atoms.cell().a(0);
	D3vector a1 = s->atoms.cell().a(1);
	D3vector a2 = s->atoms.cell().a(2);
	os << natoms << " " << -0.5*(a0+a1+a2) << endl;
      
	//write unit cell
	os << np0 << " " << a0/np0 << endl;
	os << np1 << " " << a1/np1 << endl;
	os << np2 << " " << a2/np2 << endl;
	const int nsp = s->atoms.nsp();
	for ( int is = 0; is < nsp; is++ )
	  {
	    Species* sp = s->atoms.species_list[is];
	    const int z = sp->atomic_number();
	    const int na = s->atoms.na(is);
	    for ( int ia = 0; ia < na; ia++ )
	      {
		Atom *ap = s->atoms.atom_list[is][ia];
		D3vector pos =  ap->position();
		while((pos-(a0 + a1 + a2)/2.0)*a0/length(a0) > 1) pos -= a0;
		while((pos-(a0 + a1 + a2)/2.0)*a0/length(a0) < 0) pos += a0;
		while((pos-(a0 + a1 + a2)/2.0)*a1/length(a1) > 1) pos -= a1;
		while((pos-(a0 + a1 + a2)/2.0)*a1/length(a1) < 0) pos += a1;
		while((pos-(a0 + a1 + a2)/2.0)*a2/length(a2) > 1) pos -= a2;
		while((pos-(a0 + a1 + a2)/2.0)*a2/length(a2) < 0) pos += a2;
		
		os << setprecision(5);
		os << z << " " << ((double) z) << " " << pos << endl;
	      }
	  }


	os.setf(ios::scientific,ios::floatfield);
	os << setprecision(5);
	for ( int i = 0; i < np0; i++ )
	  {
	    const int ip = (i + np0/2 ) % np0;
	    for ( int j = 0; j < np1; j++ )
	      {
		const int jp = (j + np1/2 ) % np1;
		for ( int k = 0; k < np2; k++ )
		  {
		    const int kp = (k + np2/2 ) % np2;
		    os << setw(13) << global_current[idir][ip+np0*(jp+np1*kp)];
		    if ( ( k % 6 ) == 5 )
		      os << '\n';
		  }
		if ( ( np2 % 6 ) != 0 )
		  os << '\n';
	      }
	  }

	os.close();

      }

    }
}

void CurrentDensity::plot_vtk(const Sample * s, const std::string & filename){
  using namespace std;
  Base64Transcoder xcdr;
  
  std::vector<std::vector<double> > global_current(3);
  
  for(int idir = 0; idir < 3; idir++) {
    vft()->gather(*s->wf.spincontext(0), current[idir][0], global_current[idir]);
  }

  if ( s->ctxt_.onpe0() ) {

    const int np0 = vft()->np0();
    const int np1 = vft()->np1();
    const int np2 = vft()->np2();
    
    D3vector a0 = s->atoms.cell().a(0);
    D3vector a1 = s->atoms.cell().a(1);
    D3vector a2 = s->atoms.cell().a(2);
    
    std::ofstream os;
    
    os.open((filename + ".vtk").c_str(), ios::binary);
    
    // write header and atoms
    os << "# vtk DataFile Version 2.0" << endl;
    os << "Created " << isodate() << " by " << release() << endl;
    os << "BINARY" << endl;
    os << "DATASET STRUCTURED_POINTS" << endl;
    os << "DIMENSIONS\t" << np0 << '\t' << np1 << '\t' << np2 << endl;
    os << "ORIGIN\t" << -a0[0]/2.0 << "\t" << -a1[1]/2.0 << "\t" << -a2[2]/2.0 << "\t" << endl;
    os << "SPACING\t" << a0[0]/np0 << '\t' << a1[1]/np1 << '\t' << a2[2]/np2 << endl;
    os << "POINT_DATA\t" << np0*np1*np2 << endl;
    os << "SCALARS current double 3" << endl;
    os << "LOOKUP_TABLE default" << endl;

    for ( int k = 0; k < np2; k++ ) {
      const int kp = (k + np2/2 ) % np2;
      
      for ( int j = 0; j < np1; j++ ) {
	const int jp = (j + np1/2 ) % np1;

	for ( int i = 0; i < np0; i++ ) {
	  const int ip = (i + np0/2 ) % np0;
	  
	  for(int idir = 0; idir < 3; idir++) {
	    double value = global_current[idir][ip + np0*(jp + np1*kp)];
#ifndef WORDS_BIGENDIAN
	    //Convert to big endian	    
	    xcdr.byteswap_double(1, &value);
#endif
	    os.write((char *)&value, sizeof(double));
	  }
	}
      }
    }
      
    os.close();
    
  }

}
