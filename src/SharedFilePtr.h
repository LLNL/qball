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
// SharedFilePtr.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef SHAREDFILEPTR_H
#define SHAREDFILEPTR_H

#include "mpi.h"

class SharedFilePtr
{
  private:

  MPI_Comm comm_;
  MPI_File& fh_;
  long long int offset_;

  public:

  MPI_File& file(void) { return fh_; }
  long long int offset(void) const { return offset_; }
  MPI_Offset mpi_offset(void) const { return (MPI_Offset) offset_; }
  void sync(void)
  {
    // set all offsets to the largest offset
    long long int s_off = offset_;
    MPI_Allreduce(&s_off,&offset_,1,MPI_LONG_LONG,MPI_MAX,comm_);
  }
  void set_offset(long long int off)
  {
    offset_ = off;
  }
  void advance(long long int dist)
  {
    offset_ += dist;
  }

  SharedFilePtr(MPI_Comm comm, MPI_File& fh) : comm_(comm),
                fh_(fh), offset_(0) {}
  ~SharedFilePtr() {}
};
#endif
