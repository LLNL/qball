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
// WavefunctionHandler.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef WavefunctionHANDLER_H
#define WavefunctionHANDLER_H

#include "StructureHandler.h"
#include "UnitCell.h"
#include "Wavefunction.h"
#include "SlaterDet.h"

class FourierTransform;

class WavefunctionHandler : public StructureHandler
{
  private:

  Wavefunction& wf_;
  DoubleMatrix& gfdata_;
  UnitCell uc;
  UnitCell ruc;
  double ecut;
  // dmat[ispin][ikp][i]
  std::vector<std::vector<std::vector<double> > > &dmat_;
  int& nx_;
  int& ny_;
  int& nz_;
  int current_gf_nx,current_gf_ny,current_gf_nz;
  std::string current_gf_encoding;
  int current_ispin,current_ikp,current_n,current_igf;
  std::vector<double> dmat_tmp;
  double current_kx, current_ky, current_kz, current_weight;
  int current_size;
  int read_from_gfdata;
  FourierTransform* ft;
  std::vector<std::complex<double> > wftmp;

  void byteswap_double(size_t n, double* x);

  public:

  // Start of the root element in the structure being handled
  virtual void startElement(const XMLCh* const uri,const XMLCh* const localname,
      const XMLCh* const qname, const Attributes& attributes);

  // End of the root element in the structure being handled
  virtual void endElement(const XMLCh* const uri, const XMLCh* const localname,
      const XMLCh* const qname, std::string& content);

  // start a subhandler
  virtual StructureHandler* startSubHandler(const XMLCh* const uri,
    const XMLCh* const localname, const XMLCh* const qname,
    const Attributes& attributes);

  // end a subhandler
  virtual void endSubHandler(const XMLCh* const uri,
    const XMLCh* const localname, const XMLCh* const qname,
    const StructureHandler* const subHandler);

  WavefunctionHandler(Wavefunction& wf, DoubleMatrix& gfdata,
    int& nx, int& ny, int& nz,
    std::vector<std::vector<std::vector<double> > > &dmat);
  ~WavefunctionHandler();
};
#endif

// Local Variables:
// mode: c++
// End:
