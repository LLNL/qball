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
//  BasisMapping.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef BASISMAPPING_H
#define BASISMAPPING_H

#include <complex>
#include <vector>

class Basis;
class Context;

class BasisMapping
{
  private:

  const Context& ctxt_;
  const Basis& basis_;
  int nprocs_, myproc_;

  int np0_, np1_, np2_, np012_, np012loc_;
  int nvec_;

  std::vector<int> np2_loc_;   // np2_loc_[iproc], iproc=0, nprocs_-1
  std::vector<int> np2_first_; // np2_first_[iproc], iproc=0, nprocs_-1

  std::vector<int> scounts, sdispl, rcounts, rdispl;
  std::vector<std::complex<double> > sbuf, rbuf;

  std::vector<int> ip_, im_;
  std::vector<int> ipack_, iunpack_;

  public:

  BasisMapping (const Basis &basis);
  int np0(void) const { return np0_; }
  int np1(void) const { return np1_; }
  int np2(void) const { return np2_; }
  int np2loc(void) const { return np2_loc_[myproc_]; }
  int np012(void) const { return np012_; }
  int np012loc(void) const { return np012loc_; }
  int nvec(void) const { return nvec_; }
  int zvec_size(void) const { return nvec_ * np2_; }

  const Context& context(void) const { return ctxt_; }

  // map a function c(G) to zvec_
  void vector_to_zvec(const std::complex<double> *c,
                      std::complex<double> *zvec);
  // map zvec_ to a function c(G)
  void zvec_to_vector(const std::complex<double> *zvec,
                      std::complex<double> *c);

  void transpose_fwd(const std::complex<double> *zvec,
                     std::complex<double> *ct);
  void transpose_bwd(const std::complex<double> *ct,
                     std::complex<double> *zvec);
};
#endif
