////////////////////////////////////////////////////////////////////////////////
//
// SharedFilePtr.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SharedFilePtr.h,v 1.1 2009/03/27 00:53:24 draeger1 Exp $

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
