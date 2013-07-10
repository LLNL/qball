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
//  Timer.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef TIMER_H
#define TIMER_H

#include <time.h>
#include <sys/time.h>
#include <stdint.h>

class Timer
{
  private:

  clock_t clk;
  double t,total_cpu,total_real;
  int running_;
  uint64_t nstart_;
  
  public:

  Timer() : total_cpu(0.0), total_real(0.0), running_(0), nstart_(0) {};

  void reset() { total_cpu = 0.0; total_real = 0.0; running_ = 0; };

  void start()
  {
    clk = clock();
    t = gtod();
    running_ = 1;
    nstart_++;
  };

  int running() { return running_; };

  void stop()
  {
    if ( running_ ) 
    {
      total_cpu += ((double)(clock()-clk))/CLOCKS_PER_SEC;
      total_real += gtod()-t;
      running_ = 0;
    }
  };

  double cpu()
  {
    if ( running_ ) 
    {
      return total_cpu + ((double)(clock()-clk))/CLOCKS_PER_SEC;
    }
    else
    {
      return total_cpu;
    } 
  };

  double real()
  {
    if ( running_ ) 
    {
      return total_real + gtod()-t;
    }
    else
    {
      return total_real;
    } 
  };

  double gtod(void)
  {
    static struct timeval tv;
    static struct timezone tz;
    gettimeofday(&tv,&tz);
    return tv.tv_sec + 1.e-6*tv.tv_usec;
  }

  uint64_t counts(void)
  {
     return nstart_;
  }


};
#endif
