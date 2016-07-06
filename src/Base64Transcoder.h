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
// Base64Transcoder.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef BASE64TRANSCODER_H
#define BASE64TRANSCODER_H

#include <iostream>
#include <cstdio>
#include <string>
using namespace std;
typedef unsigned char byte;

class Base64Transcoder
{
  char etable[64];  // encode table
  byte dtable[256]; // decode table

  public:
  
  Base64Transcoder();
  int encode(int nbytes, const byte* const from, char* const to);
  int decode(int nchars, const char* const from, byte* const to);
  void byteswap_double(size_t n, double* const x);
  void byteswap_int(size_t n, int* const x);
  int print(int nchars, const char* const buf, ostream& o);
  int print(const string buf, ostream& o);
  int print(int nchars, const char* const buf, FILE* outfile);
  int print(const string buf, FILE* outfile);

  // number of chars needed to encode nbytes bytes
  int nchars(int nbytes) { return 4 * ( ( nbytes + 2 ) / 3 ); }
  // number of bytes needed to decode nchars chars
  int nbytes(int nchars) { return 3 * ( nchars / 4 ); }
};

#endif
