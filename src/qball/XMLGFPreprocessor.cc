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
// XMLGFPreprocessor.C
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <valarray>
#include <sys/stat.h>

#include "Timer.h"
#include "Context.h"
#include "Base64Transcoder.h"
#include <math/matrix.h>
#include "XMLGFPreprocessor.h"
using namespace std;

////////////////////////////////////////////////////////////////////////////////
//
// XMLGFPreprocessor class
// Used to preprocess <grid_function> elements in an XML document
// Input: filename, DoubleMatrix& gfdata, string xmlcontent
// The preprocessor reads the file in parallel, processes all <grid_function>
// elements and stores the values of the grid_functions in the matrix gfdata
// which has dimensions (ngf,maxgridsize), where ngf is the total number of
// <grid_function> elements found in the file, and maxgridsize is the size of
// the largest grid_function.
// On return, the string xmlcontent contains (on task 0) the XML file
// with <grid_function> elements reduced to empty strings.
//
////////////////////////////////////////////////////////////////////////////////
void XMLGFPreprocessor::process(const char* const filename,
    DoubleMatrix& gfdata, string& xmlcontent)
{
  cout.precision(4);

  const Context& ctxt = gfdata.context();
  // define a global single row context for segment manipulations
  Context rctxt;
#if DEBUG
  if ( rctxt.oncoutpe() )
  {
    cout << "gfdata.context(): " << ctxt;
    cout << "rctxt: " << rctxt;
  }
#endif

  Timer tm,ttm;
  ttm.start();
  string st;

  struct stat statbuf;
  bool found_file = !stat(filename,&statbuf);
  // Note: process should only be called if a file is present
  assert(found_file);
#if DEBUG
    cout << " process " << ctxt.mype() << " found file "
         << filename << " size: "
         << statbuf.st_size << endl;
#endif
  off_t sz = statbuf.st_size;

  FILE* infile;
  infile = fopen(filename,"r");
  if ( !infile )
  {
    cout << " process " << ctxt.mype() << " could not open file "
         << filename << " for reading" << endl;
    return;
  }

//  ifstream is(filename);

  // determine local size
  off_t block_size = sz / ctxt.size();
  off_t local_size = block_size;
  off_t max_local_size = local_size + sz % ctxt.size();
  // adjust local_size on last task
  if ( ctxt.mype()==ctxt.size()-1 )
  {
    local_size = max_local_size;
  }

  // use contiguous read buffer, to be copied later to a string
  char *rdbuf = new char[local_size];
#if DEBUG
  cout << ctxt.mype() << ": local_size: " << local_size << endl;
#endif

  tm.start();
  off_t offset = ctxt.mype()*block_size;
  int fseek_status = fseeko(infile,offset,SEEK_SET);
  if ( fseek_status != 0 )
  {
    cout << "fseeko failed: offset=" << offset << " file_size=" << sz << endl;
  }
  assert(fseek_status==0);
  size_t items_read = fread(rdbuf,sizeof(char),local_size,infile);
  assert(items_read==local_size);

  //std::streampos offset = ctxt.mype()*block_size;
  //is.seekg(offset);
  //std::streampos new_offset = is.tellg();
  //assert(offset == new_offset);
  //is.read(&buf[0],local_size);
  //is.close();

  string buf(rdbuf,local_size);
  delete [] rdbuf;

  tm.stop();

  if ( ctxt.oncoutpe() )
  {
    cout << " XMLGFPreprocessor: read time: " << tm.real() << endl;
    cout << " XMLGFPreprocessor: local read rate: "
         << local_size/(tm.real()*1024*1024) << " MB/s"
         << "  aggregate read rate: "
         << sz/(tm.real()*1024*1024) << " MB/s" << endl;
  }

  tm.reset();
  tm.start();

  ////////////////////////////////////////////////////////////////////////////
  // fix broken tags
  ////////////////////////////////////////////////////////////////////////////
  string::size_type first_bracket_pos = buf.find_first_of("<>");
  string::size_type last_bracket_pos = buf.find_last_of("<>");
  //cout << ctxt.mype() << ": first_bracket_pos: " << first_bracket_pos
  //     << " bracket is " << buf[first_bracket_pos] << endl;
  if ( first_bracket_pos == string::npos )
  {
    // no bracket found, nothing to fix
    // cout << " no bracket found" << endl;
  }
  else if ( first_bracket_pos == last_bracket_pos )
  {
    // only one bracket found
    char bracket = buf[first_bracket_pos];
    assert(bracket=='<' || bracket=='>');
    if ( bracket == '<' )
    {
      // receive missing data from mype+1
      //cout << " found < bracket only: receive data from mype+1" << endl;
      string right_buf;
      //cout << "receiving from task " << ctxt.mype()+1 << endl;
      rctxt.string_recv(right_buf,0,rctxt.mype()+1);
      //cout << " received: " << right_buf << endl;
      buf.append(right_buf);
    }
    else
    {
      // bracket == '>'
      // send data up to bracket (inclusive) to mype-1
      string send_buf = buf.substr(0,first_bracket_pos+1);
      rctxt.string_send(send_buf,0,rctxt.mype()-1);
      buf.erase(0,first_bracket_pos+1);
      last_bracket_pos -= first_bracket_pos+1;
    }
  }
  else
  {
    // two or more brackets found

    // process first bracket
    char bracket = buf[first_bracket_pos];
    assert(bracket=='<' || bracket=='>');
    if ( bracket == '<' )
    {
      // nothing to fix
      // cout << " found < first: nothing to fix" << endl;
    }
    else
    {
      // bracket == '>'
      // send data up to first_bracket (inclusive) to task mype-1
      // cout << " found > bracket first: sending data to mype-1" << endl;
      string send_buf = buf.substr(0,first_bracket_pos+1);
      rctxt.string_send(send_buf,0,rctxt.mype()-1);
      buf.erase(0,first_bracket_pos+1);
      last_bracket_pos -= first_bracket_pos+1;
    }

    // process last bracket
    bracket = buf[last_bracket_pos];
    assert(bracket=='<' || bracket=='>');
    if ( bracket == '<' )
    {
      // receive missing data from task mype+1
      // cout << " found < bracket last: receive data from mype+1" << endl;
      string right_buf;
      // cout << "receiving from task " << rctxt.mype()+1 << endl;
      rctxt.string_recv(right_buf,0,rctxt.mype()+1);
      // cout << " received: " << right_buf << endl;
      buf.append(right_buf);
    }
    else
    {
      // bracket == '>'
      // nothing to fix
      // cout << " found > last: nothing to fix" << endl;
    }
  }

  tm.stop();
  if ( ctxt.oncoutpe() )
    cout << " XMLGFPreprocessor: tag fixing time: " << tm.real() << endl;
  tm.reset();
  tm.start();

  ////////////////////////////////////////////////////////////////////////////
  // reduce data in <grid_function> elements
  ////////////////////////////////////////////////////////////////////////////

  // locate grid_function tags
  string::size_type first_start = buf.find("<grid_function");
  string::size_type first_end = buf.find("</grid_function>");
  string::size_type last_start = buf.rfind("<grid_function");
  string::size_type last_end = buf.rfind("</grid_function>");
  bool start_tag_found = ( first_start != string::npos ) ;
  bool end_tag_found = ( first_end != string::npos ) ;

  ////////////////////////////////////////////////////////////////////////////
  // Determine if this task starts and/or ends within a grid_function element
  ////////////////////////////////////////////////////////////////////////////
  int start_within,end_within=0;
  if ( ctxt.oncoutpe() )
  {
    start_within = 0; // task 0 starts at the beginning of the XML file
  }
  else
  {
    // on pe > 0
    rctxt.irecv(1,1,&start_within,1,0,rctxt.mype()-1);
  }

  if ( !start_within )
  {
    end_within =
    ( ( start_tag_found && !end_tag_found ) ||
      ( start_tag_found && end_tag_found && last_start > last_end ) );
  }
  else
  {
    end_within =
    ( ( !end_tag_found ) ||
      ( start_tag_found && end_tag_found && last_start > last_end ) );
  }

  if ( rctxt.mype() < ctxt.size()-1 )
    rctxt.isend(1,1,&end_within,1,0,rctxt.mype()+1);

#if DEBUG
  cout << rctxt.mype() << ": start_within=" << start_within
       << " end_within=" << end_within << endl;
#endif

  ////////////////////////////////////////////////////////////////////////////
  // Determine intervals of characters to be processed
  ////////////////////////////////////////////////////////////////////////////
  // build tables seg_start[i], seg_end[i] that define intervals lying inside
  // elements

  string::size_type pos = 0;
  bool done = false;
  vector<string::size_type> seg_start,seg_end;

  // position of start and end tags on this task
  vector<string::size_type> local_start_tag_offset, local_end_tag_offset;

  // node on which starting tag igf can be found: starting_node[igf]
  vector<string::size_type> starting_node, ending_node;

  // treat first the case where the buffer starts within an element
  if ( start_within )
  {
    // get next end tag
    string::size_type next_end = buf.find("</grid_function>",pos);
    if ( next_end == string::npos )
    {
      // did not find endin tag
      // add segment [pos,buf.size()] to list
      seg_start.push_back(pos);
      seg_end.push_back(buf.size());
      done = true;
      assert(end_within);
    }
    else
    {
      // found end tag at position next_end
      local_end_tag_offset.push_back(next_end);
      // add segment [pos,next_end] to list
      seg_start.push_back(pos);
      seg_end.push_back(next_end);
      pos = next_end+16; // size of "</grid_function>"
    }
  }

  while (!done)
  {
    string::size_type next_start = buf.find("<grid_function",pos);
    if ( next_start == string::npos )
    {
      // did not find start tag
      // no start tag found, done.
      done = true;
      assert(!end_within);
    }
    else
    {
      // found start tag
      local_start_tag_offset.push_back(next_start);
      pos = buf.find(">",next_start)+1;
      string::size_type next_end = buf.find("</grid_function>",pos);
      if ( next_end == string::npos )
      {
        // did not find end tag
        // add segment [pos,buf.size()] to list
        seg_start.push_back(pos);
        seg_end.push_back(buf.size());
        done = true;
        assert(end_within);
      }
      else
      {
        // found end tag
        local_end_tag_offset.push_back(next_end);
        // add segment [pos,next_end] to list
        seg_start.push_back(pos);
        seg_end.push_back(next_end);
        pos = next_end+16; // size of "</grid_function>"
      }
    }
  }

  // compute total number of grid functions (total number of start tags)
  const int ngf_loc = local_start_tag_offset.size();
  int ngf = ngf_loc;
  rctxt.isum(1,1,&ngf,1);

  // compute index of grid function starting at the first <grid_function> tag
  int igfmin, igfmax;
  if ( rctxt.mype() == 0 )
  {
    igfmin = 0;
  }
  else
  {
    rctxt.irecv(1,1,&igfmin,1,0,rctxt.mype()-1);
  }

  igfmax = igfmin + local_end_tag_offset.size();

  if ( rctxt.mype() < rctxt.size()-1 )
  {
    rctxt.isend(1,1,&igfmax,1,0,rctxt.mype()+1);
  }

  // associate an igf value to each segment
  vector<int> igf(seg_start.size());
  for ( int i = 0; i < seg_start.size(); i++ )
    igf[i] = i+igfmin;

  // cout << rctxt.mype() << ": igfmin=" << igfmin
  //      << " igfmax=" << igfmax << endl;

  tm.stop();
  if ( ctxt.oncoutpe() )
    cout << " XMLGFPreprocessor: segment definition time: " << tm.real() << endl;

  ////////////////////////////////////////////////////////////////////////////
  // Adjust segment boundaries to allow for transcoding
  ////////////////////////////////////////////////////////////////////////////
  // If encoding is "text", the segment boundary must fall on ' ' or '\n'
  // If encoding is "base64", the segment boundary must fall between
  // groups of 4 characters in "A-Za-z0-9+/="

  ////////////////////////////////////////////////////////////////////////////
  // Determine encoding type of the last segment and send encoding_type
  // to right neighbour

  string left_encoding = "none";
  if ( start_within && !rctxt.oncoutpe() )
  {
    rctxt.string_recv(left_encoding,0,rctxt.mype()-1);
  }

  string right_encoding = "none";
  // find position of "encoding" string after last_start position
  if ( end_within )
  {
    if ( start_tag_found )
    {
      string::size_type encoding_pos = buf.find("encoding",last_start);
      assert(encoding_pos != string::npos);
      encoding_pos = buf.find("\"",encoding_pos)+1;
      string::size_type end = buf.find("\"",encoding_pos);
      right_encoding = buf.substr(encoding_pos,end-encoding_pos);
    }
    else
    {
      right_encoding = left_encoding;
    }
  }

  if ( end_within && rctxt.mype() < ctxt.size()-1 )
  {
    rctxt.string_send(right_encoding,0,rctxt.mype()+1);
  }
#if DEBUG
  cout << rctxt.mype() << ": left_encoding = " << left_encoding << endl;
  cout << rctxt.mype() << ": right_encoding = " << right_encoding << endl;
#endif

  // build table of encodings for all segments
  vector<string> encoding(seg_start.size());
  int iseg = 0;
  if ( start_within )
    encoding[iseg++] = left_encoding;
  for ( int i = 0; i < local_start_tag_offset.size(); i++ )
  {
    // determine encoding of tag at local_start_tag_offset[i]
    string::size_type encoding_pos = buf.find("encoding",local_start_tag_offset[i]);
    assert(encoding_pos != string::npos);
    encoding_pos = buf.find("\"",encoding_pos)+1;
    string::size_type end = buf.find("\"",encoding_pos);
    encoding[iseg++] = buf.substr(encoding_pos,end-encoding_pos);
  }
  assert(iseg == seg_start.size() );

  // print out table of segments [seg_start,seg_end] and corresponding igf
  assert(seg_start.size()==seg_end.size());
#if DEBUG
  for ( int i = 0; i < seg_start.size(); i++ )
  {
    cout << rctxt.mype() << ": segment [" << seg_start[i] << ","
         << seg_end[i] << "] igf=" << igf[i] << " ("
         << seg_end[i]-seg_start[i] << " bytes)"
         << " encoding=" << encoding[i]
         << endl;
  }
#endif

  tm.reset();
  tm.start();

  ////////////////////////////////////////////////////////////////////////////
  // Fix boundaries:
  ////////////////////////////////////////////////////////////////////////////

  // first fix all text_encoded boundaries
  // on all nodes involved in a text encoded boundary
  string from_right, to_left;
  if ( left_encoding == "text" || right_encoding == "text" )
  {
    if ( right_encoding == "text" )
    {
      // post string_receive from right

      // next line: cannot have right_encoding = text on last task
      assert(rctxt.mype() < rctxt.size()-1);

      rctxt.string_recv(from_right,0,rctxt.mype()+1);
      buf.append(from_right);

      // adjust ending offset of last segment
      seg_end[seg_end.size()-1] += from_right.size();

    }
    if (left_encoding == "text")
    {
      // send string up to first separator to left
      // next line: cannot have left_encoding = text on task 0
      assert(!rctxt.oncoutpe());

      string::size_type pos = buf.find_first_of(" \n<",0);
      to_left = buf.substr(0,pos);
      rctxt.string_send(to_left,0,rctxt.mype()-1);

      // adjust starting offset of first segment
      seg_start[0] = pos;

      // remove first pos characters and shift all segment indices
      buf.erase(0,pos);
      for ( int iseg = 0; iseg < seg_start.size(); iseg++ )
      {
        seg_start[iseg] -= pos;
        seg_end[iseg] -= pos;
      }
    }
  }

  // at this point, all text_encoded boundaries are fixed

  // fix base64 boundaries

  // all nodes having a base64 encoded boundary
  if ( left_encoding == "base64" || right_encoding == "base64" )
  {
    // count valid base64 chars in last segment
    // count = number of chars in "A-z,a-z,0-9,+,/,=" in last segment
    // Note: the last segment may also be the first, if only one segment
    int count = 0;
    int last_seg = seg_start.size()-1;
    // count number of valid base64 chars in segment last_seg
    for ( int i = seg_start[last_seg]; i < seg_end[last_seg]; i++ )
    {
      int ch = buf[i];
      if ( isalnum(ch) || ch == '+' || ch == '/' || ch == '=' )
        count++;
    }
#if DEBUG
    cout << rctxt.mype() << ": valid base64 chars in last seg="
         << count << endl;
#endif

    string send_left;
    int missing_on_left = 0;
    if ( left_encoding == "base64" )
    {
      // each node except node 0 posts a receive of the number of missing
      // valid chars from its left neighbour
      if ( !rctxt.oncoutpe() )
        rctxt.irecv(1,1,&missing_on_left,1,0,rctxt.mype()-1);

      // search for missing_on_left valid characters in buf[0]
      // Note: it is assumed that there are always enough valid
      // characters in the first segment to align on a 128 char boundary
      // 3-doubles = 192 bits
      // one 4-char block encodes 24 bits
      // boundary at 8 4-char blocks = 32 chars
      assert(missing_on_left < 32);

      // Note: the number of chars missing on left can be larger than
      // the total number of chars available in the first segment.
      // In that case, all characters in the first segment are sent
      string::size_type pos = 0;
      int n = 0;
      while ( pos < seg_end[0] && n < missing_on_left )
      {
        char c = buf[pos++];
        if ( isalnum(c) || c == '+' || c == '/' || c == '=' )
        {
          send_left.push_back(c);
          n++;
        }
      }
      // adjust start of segment offset
      seg_start[0] = pos;
      // remove first pos characters and shift all segment indices
      buf.erase(0,pos);
      for ( int iseg = 0; iseg < seg_start.size(); iseg++ )
      {
        seg_start[iseg] -= pos;
        seg_end[iseg] -= pos;
      }
    }

    // if the last segment is also the first (i.e. only one segment)
    // remove chars sent to left from valid char count
    if ( seg_start.size() == 1 )
      count -= missing_on_left;

    int missing_on_right = 0;
    if ( right_encoding == "base64" )
    {
      // compute number of missing valid chars on the right
      // Request enough characters to fall on a 32-char boundary
      // number of missing *valid* chars requested is (32 - count % 32)
      // Note: there may not be enough chars in the first segment of the
      // right neighbor to return that many chars. In that case, fewer
      // chars will be received from the right neighbor
      int countmod32 = count % 32;
      if ( countmod32 != 0 )
        missing_on_right = 32 - countmod32;
#if DEBUG
      cout << rctxt.mype() << ": missing_on_right="
      << missing_on_right << endl;
#endif

      // each node except the last sends the number of missing chars to its
      // right neighbour (receive is already posted by other nodes)
      // each node except the last posts a string receive from the right node

      if ( rctxt.mype() < rctxt.size()-1 )
      {
        rctxt.isend(1,1,&missing_on_right,1,0,rctxt.mype()+1);
        // receive string from task mype+1
        string missing_chars;
        rctxt.string_recv(missing_chars,0,rctxt.mype()+1);
        // append to buf
        buf.append(missing_chars);
        seg_end[last_seg] += missing_chars.size();
      }
    }

    if ( left_encoding == "base64" )
    {
      // each node except node 0 sends string "send_left" to left node
      if ( !rctxt.oncoutpe() )
      {
        rctxt.string_send(send_left,0,rctxt.mype()-1);
#if DEBUG
        cout << rctxt.mype() << ": sent " << send_left.size()
             << " chars to left task: " << send_left << endl;
#endif
      }
    }
  }

  // print out table of segments [seg_start,seg_end] and corresponding igf
#if DEBUG
  for ( int i = 0; i < seg_start.size(); i++ )
  {
    cout << rctxt.mype() << ": segment [" << seg_start[i] << ","
         << seg_end[i] << "] igf=" << igf[i] << " ("
         << seg_end[i]-seg_start[i] << " bytes)"
         << " encoding=" << encoding[i]
         << endl;
  }
#endif

  tm.stop();
  if ( ctxt.oncoutpe() )
    cout << " XMLGFPreprocessor: boundary adjustment time: "
         << tm.real() << endl;
  tm.reset();
  tm.start();

  // Transcode segments
  // Instantiate a Base64Transcoder
  Base64Transcoder xcdr;

  vector<vector<double> > dbuf;
  dbuf.resize(seg_start.size());
  for ( int iseg = 0; iseg < seg_start.size(); iseg++ )
  {
    int nchars = seg_end[iseg]-seg_start[iseg];
    if ( encoding[iseg] == "base64" )
    {
      // Base64 case:
      // the dimension of the dbuf array is computed from the *total* number
      // of bytes in the segment, i.e. will be overestimated slightly
      // since not all bytes are valid base64 bytes
      dbuf[iseg].resize((((nchars/4)*3)/8));
#if DEBUG
      cout << rctxt.mype() << ": iseg=" << iseg
           << " dbufsize=" << dbuf[iseg].size()
           << endl;
#endif
      int nbytes = xcdr.decode(nchars,buf.data()+seg_start[iseg],
                               (byte*)&dbuf[iseg][0]);
#if DEBUG
      cout << rctxt.mype() << ": iseg=" << iseg << " nbytes=" << nbytes
           << endl;
#endif
      assert(nbytes % 8 == 0 );
      int ndoubles = nbytes / 8;
      assert(ndoubles <= dbuf[iseg].size());
      // adjust size of double array
      dbuf[iseg].resize(ndoubles);
    }
    else
    {
      // text encoding
#if DEBUG
      cout << rctxt.mype() << ": text in iseg=" << iseg << ": "
           << buf.substr(seg_start[iseg],80) << endl;
#endif
      istringstream stst(buf.substr(seg_start[iseg],nchars));
      double d;
      while ( stst >> d )
      {
        dbuf[iseg].push_back(d);
      }
#if DEBUG
      cout << rctxt.mype() << ": doubles read from text: "
          << dbuf[iseg].size() << endl;
#endif
    }
  }

  tm.stop();
  if ( ctxt.oncoutpe() )
    cout << " XMLGFPreprocessor: transcoding time: " << tm.real() << endl;
  tm.reset();
  tm.start();

  ////////////////////////////////////////////////////////////////////////////
  // byte swapping on big-endian platforms
#ifdef WORDS_BIGENDIAN
  for ( int iseg = 0; iseg < seg_start.size(); iseg++ )
    xcdr.byteswap_double(dbuf[iseg].size(),&dbuf[iseg][0]);

  tm.stop();
  if ( ctxt.oncoutpe() )
    cout << " XMLGFPreprocessor: byte swapping time: " << tm.real() << endl;
  tm.reset();
  tm.start();
#endif

#if DEBUG
  for ( int iseg = 0; iseg < seg_start.size(); iseg++ )
  {
    cout << rctxt.mype() << ": seg_data[iseg=" << iseg << "]: size=" << dbuf[iseg].size() << "  ";
    if ( dbuf[iseg].size() >=3 )
      cout << dbuf[iseg][0] << " " << dbuf[iseg][1]
           << " " << dbuf[iseg][2];
    cout << endl;
  }
#endif

  ////////////////////////////////////////////////////////////////////////////
  // Collect all data in dbuf into a DoubleMatrix of dimension maxgfdim x ngf
  // Note: the largest grid must be used to dimension the matrix.
  // Determine the largest dimension of the grid_functions
  // use dmax call

  // compute size of largest grid
  int maxgfsize = 0;
  valarray<int> gfsize(ngf);
  gfsize = 0;
  for ( int iseg = 0; iseg < seg_start.size(); iseg++ )
  {
    for ( int i = 0; i < ngf; i++ )
      if ( igf[iseg] == i )
        gfsize[i] = dbuf[iseg].size();
  }
  rctxt.isum(ngf,1,&gfsize[0],ngf);
  if ( ngf > 0 )
    maxgfsize = gfsize.max();
#if DEBUG
  cout << rctxt.mype() << ": maxgfsize=" << maxgfsize << endl;
  for ( int i = 0; i < ngf; i++ )
    cout << rctxt.mype() << ": igf=" << i << " gfsize=" << gfsize[i]
         << endl;
#endif

  // determine block sizes of matrix gfdata
  // note: use ctxt, not rctxt
  int gfmb = maxgfsize / ctxt.nprow() +
           ( maxgfsize % ctxt.nprow() != 0 ? 1 : 0 );
  int gfnb = ngf / ctxt.npcol() +
           ( ngf % ctxt.npcol() != 0 ? 1 : 0 );
  gfdata.resize(maxgfsize,ngf,gfmb,gfnb);
#if DEBUG
  cout << ctxt.mype() << ": gfdata resized: (" << maxgfsize << "x" << ngf
       << ")  (" << gfmb << "x" << gfnb << ") blocks" << endl;
  cout << ctxt.mype() << ": gfdata.context(): " << gfdata.context();
#endif

  // prepare buffer sbuf for all_to_all operation
  int sbufsize = 0;
  for ( int iseg = 0; iseg < seg_start.size(); iseg++ )
    sbufsize += dbuf[iseg].size();

  valarray<double> sbuf(sbufsize);

  // determine offset in the grid of the first element of the first
  // segment
  // Each task having start_within==true posts a receive from its left
  // task
  // Each task having end_within==true computes the offset of the last
  // element of its last segment and sends it to its right neighbour

  int left_offset = 0;
  if ( start_within )
  {
    rctxt.irecv(1,1,&left_offset,1,0,rctxt.mype()-1);
  }
  if ( end_within )
  {
    // compute offset of last element of last segment
    int right_offset = 0;
    if ( dbuf.size() == 1 )
    {
      // if there is only one segment, take left offset into account
      right_offset = left_offset + dbuf[seg_start.size()-1].size();
    }
    else
    {
      right_offset = dbuf[seg_start.size()-1].size();
    }
#if DEBUG
    cout << rctxt.mype() << ": right_offset=" << right_offset << endl;
#endif
    rctxt.isend(1,1,&right_offset,1,0,rctxt.mype()+1);
  }
#if DEBUG
  cout << rctxt.mype() << ": left_offset=" << left_offset << endl;
#endif

  // compute array sbuf
  // copy into sbuf all data to be sent
  // data going to task irow icol is in subsegment irow
  // for each task itask
  //   compute irow,icol of itask, igfmin, igfmax, min_offset, max_offset
  //   for each segment iseg
  //      if igf[iseg] falls on task itask
  //        for each block of segment iseg
  //          if offset of this block falls in [min_offset,max_offset]
  //            copy block to sbuf
  //            adjust scount[itask]
  valarray<int> scounts(0,rctxt.size()), sdispl(0,rctxt.size());
  valarray<int> rcounts(0,rctxt.size()), rdispl(0,rctxt.size());

  int sbuf_pos = 0;
  for ( int itask = 0; itask < ctxt.size(); itask++ )
  {
    // determine block sizes of the block cyclic distribution of gfdata
    const int irow = itask % ctxt.nprow();
    const int icol = itask / ctxt.nprow();

    const int igfmin = icol * gfnb;
    // nb_loc: size of block on task icol
    // nb_loc is gfnb       for icol < kn
    //           ngf-k*gfnb for icol == kn
    //           0          for icol > kn
    // where kn = int(ngf/gfnb)
    const int kn = (gfnb>0) ? ngf / gfnb : 0;
    int nb_loc;
    if ( icol < kn )
      nb_loc = gfnb;
    else if ( icol == kn )
      nb_loc = ngf - kn * gfnb;
    else
      nb_loc = 0;
    const int igfmax = igfmin + nb_loc;

    const int min_offset = irow * gfmb;
    // mb_loc: size of block on task irow
    // mb_loc is gfmb             for irow < km
    //           maxgfsize-k*gfmb for irow == km
    //           0                for irow > km
    // where km = int(maxgfsize/gfmb)
    const int km = (gfmb>0) ? maxgfsize / gfmb : 0;
    int mb_loc;
    if ( irow < km )
      mb_loc = gfmb;
    else if ( irow == km )
      mb_loc = maxgfsize - km * gfmb;
    else
      mb_loc = 0;

    const int max_offset = min_offset + mb_loc;

    sdispl[itask] = sbuf_pos;
    for ( int iseg = 0; iseg < seg_start.size(); iseg++ )
    {
      const int segsize = dbuf[iseg].size();

      if ( igf[iseg] >= igfmin && igf[iseg] < igfmax )
      {
        // igf falls within the states stored on itask

        // Find if there is an intersection of this segment
        // with the range [min_offset,max_offset[ on task itask

        int block_start, block_size = 0;

        if ( iseg == 0 )
        {
          // On segment 0, the available range is
          // [left_offset,left_offset+segsize[
          // There is overlap if min_offset < left_offset+segsize
          // and max_offset > left_offset
          // If there is overlap, the range to be sent is
          // [max(min_offset-left_offset,0) ,
          //  min(max_offset-left_offset,segsize)[
          if ( min_offset < left_offset+segsize && max_offset > left_offset )
          {
            // there is overlap
            block_start = max(min_offset-left_offset,0);
            block_size = min(max_offset-left_offset,segsize) - block_start;
          }
        }
        else
        {
          // On all segments except segment 0, the available range
          // is [0,segsize[
          // There is overlap if min_offset < segsize
          // (and if max_offset > 0, but this is always true)
          // If there is overlap, the range to be sent is
          // [min_offset, min(max_offset,segsize[
          if ( min_offset < segsize )
          {
            // there is overlap
            block_start = min_offset;
            block_size = min(max_offset,segsize) - block_start;
          }
        }

        if ( block_size > 0 )
        {
          // copy block to sbuf
          memcpy((void*)&sbuf[sbuf_pos],(void*)&dbuf[iseg][block_start],
                 block_size*sizeof(double));
          sbuf_pos += block_size;
          scounts[itask] += block_size;

#if DEBUG
          cout << rctxt.mype() << ": sbuf: iseg=" << iseg << " sending start="
          << block_start << " size=" << block_size << " to task " << itask
          << " (irow=" << irow << ",icol=" << icol << ")" << endl;
#endif
        }
      }
    }
  }

  // send scount array using all_to_all call
  valarray<int> a2a_scounts(1,ctxt.size()),a2a_rcounts(1,ctxt.size()),
                a2a_sdispl(0,ctxt.size()),a2a_rdispl(0,ctxt.size());
  for ( int i = 0; i < ctxt.size(); i++ )
  {
    a2a_rdispl[i] = a2a_sdispl[i] = i;
  }
  int status = MPI_Alltoallv(&scounts[0],&a2a_scounts[0],&a2a_sdispl[0],
     MPI_INT,&rcounts[0],&a2a_rcounts[0],&a2a_rdispl[0],MPI_INT,
     rctxt.comm());
  assert(status==0);

  rdispl[0] = 0;
  for ( int i = 1; i < ctxt.size(); i++ )
    rdispl[i] = rdispl[i-1] + rcounts[i-1];

  int rbufsize = 0;
  for ( int i = 0; i < ctxt.size(); i++ )
    rbufsize += rcounts[i];
#if DEBUG
  cout << ctxt.mype() << ": "
       << " rbufsize: " << rbufsize << endl;
#endif

  valarray<double> rbuf(rbufsize);


  status = MPI_Alltoallv(&sbuf[0],&scounts[0],&sdispl[0],MPI_DOUBLE,
           &rbuf[0],&rcounts[0],&rdispl[0],MPI_DOUBLE,rctxt.comm());
  assert(status==0);

  // copy data from rbuf to gfdata.valptr()
  // functions in rbuf can have varying length
  // determine bounds igfmin, igfmax of functions in rbuf
  // then use gfsize[igf] to copy one function at a time to gfdata
  //
  // compute block sizes for current node
  const int irow = ctxt.myrow();
  const int icol = ctxt.mycol();

  const int igfminloc = icol * gfnb;
  // nb_loc: size of block on task icol
  // nb_loc is gfnb       for icol < kn
  //           ngf-k*gfnb for icol == kn
  //           0          for icol > kn
  // where kn = int(ngf/gfnb)
  const int kn = (gfnb>0) ? ngf / gfnb : 0;
  int nb_loc;
  if ( icol < kn )
    nb_loc = gfnb;
  else if ( icol == kn )
    nb_loc = ngf - kn * gfnb;
  else
    nb_loc = 0;
  const int igfmaxloc = igfminloc + nb_loc;

#if DEBUG
  cout << ctxt.mype() << ": "
       << " igfminloc: " << igfminloc << " igfmaxloc: " << igfmaxloc << endl;
#endif
  int rbuf_pos = 0;
  int igfloc = 0;
  for ( int igf = igfminloc; igf < igfmaxloc; igf++ )
  {
    // mb_loc: size of block on task irow
    // mb_loc is gfmb             for irow < km
    //           gfsize-k*gfmb for irow == km
    //           0                for irow > km
    // where km = int(gfsize/gfmb)
    const int km = (gfmb>0) ? gfsize[igf] / gfmb : 0;
    int mb_loc;
    if ( irow < km )
      mb_loc = gfmb;
    else if ( irow == km )
      mb_loc = gfsize[igf] - km * gfmb;
    else
      mb_loc = 0;

#if DEBUG
    cout << " copying gfsize[" << igf << "] = " << gfsize[igf]
         << " block size: " << mb_loc
         << " gfdata.mloc(): " << gfdata.mloc() << endl;
#endif

    // copy rbuf[rbuf_pos] size mb_loc to
    // gfdata.valptr(icol*mloc())
    const int dest = igfloc * gfdata.mloc();
    memcpy(gfdata.valptr(dest),&rbuf[rbuf_pos],mb_loc*sizeof(double));
    rbuf_pos += mb_loc;
    igfloc++;
  }

  tm.stop();
  if ( ctxt.oncoutpe() )
    cout << " XMLGFPreprocessor: data redistribution time: "
         << tm.real() << endl;
  tm.reset();
  tm.start();

  // compact XML data:
  // erase <grid_function> contents from XML data
  // Note: when removing the data from <grid_function>, add
  // the attribute xsi:null="true"
  // Also add to the schema definition: nillable="true" in the
  // definition of element grid_function

  // delete segment data
  for ( int iseg = seg_start.size()-1; iseg >= 0; iseg--)
  {
    //cout << " erasing segment: ["
    //     << seg_start[iseg] << "," << seg_end[iseg] << "]" << endl;
    buf.erase(seg_start[iseg],seg_end[iseg]-seg_start[iseg]);
  }
  //cout << " buf.size() after erase: " << buf.size() << endl;

  // collect all XML data on task 0
  // Distribute sizes of local strings to all tasks
  // and store in array rcounts

  int xmlcontentsize = buf.size();
  status = MPI_Allgather(&xmlcontentsize,1,MPI_INT,&rcounts[0],1,
             MPI_INT,rctxt.comm());
  assert(status==0);
  // compute total size
  int totalsize = 0;
  for ( int i = 0; i < rctxt.size(); i++ )
    totalsize += rcounts[i];

  rdispl[0] = 0;
  for ( int i = 1; i < rctxt.size(); i++ )
    rdispl[i] = rdispl[i-1] + rcounts[i-1];

  xmlcontent.resize(totalsize);

  // Allgatherv xml content
  status = MPI_Allgatherv(&buf[0],xmlcontentsize,MPI_CHAR,&xmlcontent[0],
    &rcounts[0],&rdispl[0],MPI_CHAR,rctxt.comm());

  tm.stop();
  if ( ctxt.oncoutpe() )
    cout << " XMLGFPreprocessor: XML compacting time: " << tm.real() << endl;
  tm.reset();
  tm.start();

  // data is now available in the form of the DoubleMatrix gfdata
  //
  // Redistribution of the data can be done using the DoubleMatrix::getsub
  // and DoubleMatrix::operator= members
  //
  // 1) Use getsub to get only wf (not dwf) in a DoubleMatrix distributed over
  //    the global context
  // 2) Use operator= to redistribute the gridfunctions to the same
  //    context as Wavefunction

  // Fourier transforms
  // WavefunctionHandler can use the real-space data in the redistributed
  // gfdata matrix to compute Fourier coefficients, which completes the
  // load operation

  ////////////////////////////////////////////////////////////////////////////
  // end of program

  ctxt.barrier();
  ttm.stop();
  if ( ctxt.oncoutpe() )
  {
    cout << " XMLGFPreprocessor: total time: " << ttm.real() << endl;
  }
}
