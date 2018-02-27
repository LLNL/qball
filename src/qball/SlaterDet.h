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
// SlaterDet.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef SLATERDET_H
#define SLATERDET_H

class FourierTransform;
#include "Context.h"
#include "Basis.h"
#include <math/Matrix.h>

#include <math/D3vector.h>
#include <vector>
#include <iostream>
#include "Timer.h"
#include <string>
#include <map>
#if USE_CSTDIO_LFS
#include <cstdio>
#endif
using namespace std;

class SharedFilePtr;
class AtomSet;

typedef map<string,Timer> TimerMap;

class SlaterDet {
  private:

  Context& ctxt_;
  Context& ctxtsq_;
  const Context& col_ctxt_;
  Basis* basis_;
  ComplexMatrix c_;
  vector<double> occ_;
  vector<double> eig_;
  int nbetagloc_;                     // number of local ultrasoft betag 
  vector<ComplexMatrix*> betag_;      // ultrasoft beta function
  vector<ComplexMatrix*> betapsi_;    // <beta|phi>
  vector<ComplexMatrix*> dbetapsi_;    // d<beta|phi>/dRj_I
  ComplexMatrix spsi_;                // S*phi
  AtomSet* atoms_;
  vector<vector<double> > qaug_;      // q_nm storage, used to calc S*psi
  
  void byteswap_double(size_t n, double* x);
  double fermi(double e, double mu, double fermitemp);
  double methfessel(double e, double mu, double width, int ngauss);
  bool ultrasoft_;
  bool force_complex_;
  bool highmem_;
  int mbset_, nbset_;
  int mblks_, nblks_;
  
  public:
  
  mutable TimerMap tmap;

  SlaterDet(Context& ctxt, const Context& my_col_ctxt, Context& ctxtsq, D3vector kpoint, bool ultrasoft, bool force_complex);
  SlaterDet(const SlaterDet& rhs);
  ~SlaterDet();
  void print_timing(void);
  Context& context(void) const { return ctxt_; }
  const Context& col_ctxt(void) const { return col_ctxt_; }
  const Basis& basis(void) const { return *basis_; }
  Basis& basis(void) { return *basis_; }
  const D3vector kpoint(void) const { return basis_->kpoint(); }
  const ComplexMatrix& c(void) const { return c_; }
  ComplexMatrix& c(void) { return c_; }
  const ComplexMatrix* betag(int is) const { return betag_[is]; }
  ComplexMatrix* betag(int is) { return betag_[is]; }
  const ComplexMatrix* betapsi(int is) const { return betapsi_[is]; }
  ComplexMatrix* betapsi(int is) { return betapsi_[is]; }
  ComplexMatrix* dbetapsi(int is) { return dbetapsi_[is]; }
  const ComplexMatrix* dbetapsi(int is) const { return dbetapsi_[is]; }
  const ComplexMatrix& spsi(void) const { return spsi_; }
  ComplexMatrix& spsi(void) { return spsi_; }
  const vector<double>& occ(void) const { return occ_; }
  const vector<double>& eig(void) const { return eig_; }
  int nst(void) const { return c_.n(); }
  int nstloc(void) const { return c_.nloc(); }
  void resize(const UnitCell& cell, const UnitCell& refcell,
              double ecut, int nst);
  void reshape(const Context& newctxt, const Context& new_col_ctxt, const Context& newctxtsq, bool setnewctxt);
  void copyTo(SlaterDet* newsd);
  void compute_density(FourierTransform& ft, double weight, double* rho) const;
  void compute_density(FourierTransform& ft, double weight, std::complex<double> * rho, const SlaterDet & sd2_) const;
  void rs_mul_add(FourierTransform& ft, const double* v, SlaterDet& sdp) const;
  void randomize(double amplitude);
  void randomize_real(double amplitude);
  void randomize_us(double amplitude, AtomSet& as);
  void rescale(double factor);
  // AS: shift state n_state by the vector (shift_x, shift_y, shift_z)
  void shift_wf(double shift_x,double shift_y,double shift_z,int n_state);
  // AS: change phase of the wave function to make it real for Gamma only
  void phase_wf_real(void);
  void cleanup(void);
  void reset(void);
  void gram();
  void set_local_block(int mb, int nb);
  void set_nblocks(int mblks, int nblks);
  void riccati(SlaterDet& sd);
  void lowdin();
  void align(const SlaterDet& sd);
  void ortho_align(const SlaterDet& sd);
  double dot(const SlaterDet& sd) const;
  double sdot(const SlaterDet& sd) const;
  double total_charge(void);
  void update_occ(int nel, int nspin);
  void update_occ(int nspin, double mu, double temp, int ngauss);
  void update_index_of_minus_g(void) { basis_->update_index_of_minus_g(); }
  double eig(int i) const { return eig_[i]; };
  const double* eig_ptr(void) const { return &eig_[0]; }
  const double* eig_ptr(int i) const { return &eig_[i]; }
  double occ(int i) const { return occ_[i]; };
  const double* occ_ptr(void) const { return &occ_[0]; }
  const double* occ_ptr(int i) const { return &occ_[i]; }
  void set_occ(vector<double>& occ)
    { assert(occ_.size()==occ.size()); occ_ = occ; }
  void set_eig(vector<double>& eig)
    { assert(eig_.size()==eig.size()); eig_ = eig; }
  void set_eig(valarray<double>& eig)
    { assert(eig_.size()==eig.size()); 
      for (unsigned i = 0; i < eig.size(); i++ )
        eig_[i] = eig[i];
    }
  void promote_occ(double occ_change, int origin_level, int destination_level);
  
  bool ultrasoft(void) { return ultrasoft_; }
  bool highmem(void) { return highmem_; }
  void set_highmem(void) { highmem_ = true; }
  void init_usfns(AtomSet* atoms);
  void update_usfns();
  void calc_betag();
  void calc_betapsi(void);
  void calc_dbetapsi(int j);
  void calc_spsi();
  void calc_wfphase();
  void set_qaug(int is, vector<double>& qaug);
  double qaug(int is, int qind) { return qaug_[is][qind]; }
  double entropy(int nspin);
  double ortho_error(void);
  double memsize(void) const;
  double localmemsize(void) const;
  SlaterDet& operator=(SlaterDet& rhs);
  void print(ostream& os, string encoding);
  void write(SharedFilePtr& fh, std::string encoding, double weight, int ispin,
             int nspin) const;

  void info(ostream& os);
  void print_memory(ostream&os, int kmult, int kmultloc, double& totsum, double& locsum) const;
};
ostream& operator << ( ostream& os, SlaterDet& sd );

class SlaterDetException
{
  public:
  string msg;
  SlaterDetException(string s) : msg(s) {}
};
#endif

// Local Variables:
// mode: c++
// End:
