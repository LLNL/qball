////////////////////////////////////////////////////////////////////////////////  
// Copyright (c) 2017, Lawrence Livermore National Security, LLC. 
// qb@ll:  Qbox at Lawrence Livermore
//
// This file is part of qb@ll.
//
// Produced at the Lawrence Livermore National Laboratory. 
// Written by Xavier Andrade (xavier@llnl.gov) and Francois Gygi (fgygi@ucdavis.edu).
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

#include <config.h>

#ifndef DFTD3_H
#define DFTD3_H

#define dftd3_init FC_FUNC_(f90_dftd3_init, F90_DFTD3_INIT) 
#define dftd3_end FC_FUNC_(f90_dftd3_end, F90_DFTD3_END)
#define dftd3_pbc_dispersion FC_FUNC_(f90_dftd3_pbc_dispersion, F90_DFTD3_PBC_DISPERSION) 
  
extern "C" {
  
  void dftd3_init();
  void dftd3_end();
  void dftd3_pbc_dispersion(const int * natoms, const double * coords, const int * izp, const double * latvecs, double * disp, double * grads, double * stress);

}


#endif

// Local Variables:
// mode: c++
// End:
