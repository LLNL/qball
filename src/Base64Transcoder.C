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
// Base64Transcoder.C
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#include "Base64Transcoder.h"
#include <iostream>
#include <cstdio>
#include <string>
#include <cassert>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
Base64Transcoder::Base64Transcoder()
{
  // initialize encoding/decoding tables
  
  for (int i = 0; i < 26; i++)
  {
      etable[i] = 'A' + i;
      etable[26 + i] = 'a' + i;
  }
  for (int i = 0; i < 10; i++)
  {
      etable[52 + i] = '0' + i;
  }
  etable[62] = '+';
  etable[63] = '/';

  for (int i = 0; i < 255; i++)
  {
      dtable[i] = 0x80;
  }
  for (int i = 'A'; i <= 'Z'; i++)
  {
      dtable[i] = 0 + (i - 'A');
  }
  for (int i = 'a'; i <= 'z'; i++)
  {
      dtable[i] = 26 + (i - 'a');
  }
  for (int i = '0'; i <= '9'; i++)
  {
      dtable[i] = 52 + (i - '0');
  }
  dtable['+'] = 62;
  dtable['/'] = 63;
  dtable['='] = 0;
}

////////////////////////////////////////////////////////////////////////////////
int Base64Transcoder::encode(int nbytes, const byte* const from, char* const to)
{
  const byte* fptr = from;
  char* tptr = to;
  
  int n3 = nbytes / 3; // number of groups of three bytes
  
  while ( n3-- > 0 )
  {
    byte ig0 = *fptr++;
    byte ig1 = *fptr++;
    byte ig2 = *fptr++;

    *tptr++ = etable[ig0 >> 2];
    *tptr++ = etable[((ig0 & 3) << 4) | (ig1 >> 4)];
    *tptr++ = etable[((ig1 & 0xF) << 2) | (ig2 >> 6)];
    *tptr++ = etable[ig2 & 0x3F];    
  }
  
  int nr = nbytes % 3; // remaining bytes
  
  if ( nr == 2 )
  {
    byte ig0 = *fptr++;
    byte ig1 = *fptr++;
    byte ig2 = 0;
    
    *tptr++ = etable[ig0 >> 2];
    *tptr++ = etable[((ig0 & 3) << 4) | (ig1 >> 4)];
    *tptr++ = etable[((ig1 & 0xF) << 2) | (ig2 >> 6)];
    *tptr++ = '=';
  }
  else if ( nr == 1 )
  {
    byte ig0 = *fptr++;
    byte ig1 = 0;
    
    *tptr++ = etable[ig0 >> 2];
    *tptr++ = etable[((ig0 & 3) << 4) | (ig1 >> 4)];
    *tptr++ = '=';
    *tptr++ = '=';
  }
  
  return 0;  
}

////////////////////////////////////////////////////////////////////////////////
int Base64Transcoder::decode(const int nchars, const char* const from, 
  byte* const to)  
{
  // Decode Base64 chars in array "from" into bytes in array "to"
  // White space and new lines are skipped
  // extra characters at end that do not form a valid group of 4 chars are
  // ignored.
  // nchars: number of chars in array "from"
  // the number of bytes successfully translated is returned
  
  byte a0,a1,a2,a3,b0,b1,b2,b3;
  int c;
  const char* fptr = from;
  const char* const fptr_end = from+nchars+1;
  byte* tptr = to;
  
  while ( fptr < fptr_end-4 )
  {
    // get 4 valid characters from input string
    do
    {
      c = *fptr++;
    }
    while ( (c <= ' ') && (fptr < fptr_end) );
    if ( fptr >= fptr_end )
    {
#ifdef DEBUG
      cerr << " Base64Transcoder::decode: end of string reached reading c0 "
           << endl;
#endif
      break;
    }
    a0 = (byte) c;
    b0 = (byte) dtable[c];
 
    do
    {
      c = *fptr++;
    }
    while ( (c <= ' ') && (fptr < fptr_end) );
    if ( fptr >= fptr_end )
    {
#ifdef DEBUG
      cerr << " Base64Transcoder::decode: end of string reached reading c1 "
           << endl;
#endif
      break;
    }
    a1 = (byte) c;
    b1 = (byte) dtable[c];
 
    do
    {
      c = *fptr++;
    }
    while ( (c <= ' ') && (fptr < fptr_end) );
    if ( fptr >= fptr_end )
    {
#ifdef DEBUG
      cerr << " Base64Transcoder::decode: end of string reached reading c2 "
           << endl;
#endif
      break;
    }
    a2 = (byte) c;
    b2 = (byte) dtable[c];
 
    do
    {
      c = *fptr++;
    }
    while ( (c <= ' ') && (fptr < fptr_end) );
    if ( (c <= ' ') && fptr >= fptr_end )
    {
#ifdef DEBUG
      cerr << " Base64Transcoder::decode: end of string reached reading c3\n"
           << " (without reading a valid c3) " << endl;
#endif
      break;
    }
    a3 = (byte) c;
    b3 = (byte) dtable[c];
 
    if ((b0|b1|b2|b3) & 0x80)
    {
#ifdef DEBUG
      cerr << " Base64Transcoder::decode: Illegal character in input: "
           << endl;
#endif
      return tptr - to;
    }
      
    if ( a3 == '=' )
    {
      if ( a2 == '=' )
      {
        // write 1 byte only
        *tptr++ = (b0 << 2) | (b1 >> 4);
      }
      else
      {
        // write 2 bytes only
        *tptr++ = (b0 << 2) | (b1 >> 4);
        *tptr++ = (b1 << 4) | (b2 >> 2);
      }
    }
    else
    {
      // write 3 bytes
      *tptr++ = (b0 << 2) | (b1 >> 4);
      *tptr++ = (b1 << 4) | (b2 >> 2);
      *tptr++ = (b2 << 6) | b3;
    }
  
  }
#ifdef DEBUG
  if ( fptr >= fptr_end )
  {
    cerr << " Base64Transcoder::decode: end of string reached in input: "
         << endl;
  }
#endif
  
  return tptr - to;
}

////////////////////////////////////////////////////////////////////////////////
void Base64Transcoder::byteswap_double(size_t n, double* const x)
{
  if (n==0) return;
  unsigned char* c = (unsigned char*) x;
  while ( n-- > 0 )
  {
    unsigned char tmp;
    tmp = c[7]; c[7] = c[0]; c[0] = tmp;
    tmp = c[6]; c[6] = c[1]; c[1] = tmp;
    tmp = c[5]; c[5] = c[2]; c[2] = tmp;
    tmp = c[4]; c[4] = c[3]; c[3] = tmp;
    
    c+=8;
  }
}  

////////////////////////////////////////////////////////////////////////////////
void Base64Transcoder::byteswap_int(size_t n, int* const x)
{
  if (n==0) return;
  unsigned char* c = (unsigned char*) x;
  while ( n-- > 0 )
  {
    unsigned char tmp;
    tmp = c[3]; c[3] = c[0]; c[0] = tmp;
    tmp = c[2]; c[2] = c[1]; c[1] = tmp;
    
    c+=4;
  }
} 

////////////////////////////////////////////////////////////////////////////////
int Base64Transcoder::print(const string buf, ostream& o)
{
  const char* b = buf.c_str();
  int nchars = buf.size();
  int nl = nchars / 72;
  for ( int i = 0; i < nl; i++ )
  {
    o.write(b,72);
    o << endl;
    b += 72;
  }
  
  o.write(b,nchars%72);
  o << endl;
  return 0;
}

////////////////////////////////////////////////////////////////////////////////
int Base64Transcoder::print(int nchars, const char* const buf, ostream& o)
{
  const char* b = buf;
  int nl = nchars / 72;
  for ( int i = 0; i < nl; i++ )
  {
    o.write(b,72);
    o << endl;
    b += 72;
  }
  
  o.write(b,nchars%72);
  o << endl;
  return 0;
}
////////////////////////////////////////////////////////////////////////////////
int Base64Transcoder::print(const string buf, FILE* outfile)
{
  const char* b = buf.c_str();
  const char nlchar = '\n';
  int nchars = buf.size();
  int nl = nchars / 72;
  for ( int i = 0; i < nl; i++ )
  {
    fwrite(b,sizeof(char),72,outfile);
    fwrite(&nlchar,sizeof(char),1,outfile);
    b += 72;
  }
  
  if ( nchars%72 != 0 )
  {
    fwrite(b,sizeof(char),nchars%72,outfile);
    fwrite(&nlchar,sizeof(char),1,outfile);
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////////
int Base64Transcoder::print(int nchars, const char* const buf, FILE* outfile)
{
  const char* b = buf;
  const char nlchar = '\n';
  int nl = nchars / 72;
  for ( int i = 0; i < nl; i++ )
  {
    fwrite(b,sizeof(char),72,outfile);
    fwrite(&nlchar,sizeof(char),1,outfile);
    b += 72;
  }
  
  if ( nchars%72 != 0 )
  {
    fwrite(b,sizeof(char),nchars%72,outfile);
    fwrite(&nlchar,sizeof(char),1,outfile);
  }
  return 0;
}
