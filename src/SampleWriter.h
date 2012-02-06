////////////////////////////////////////////////////////////////////////////////
//
// SampleWriter.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SampleWriter.h,v 1.1 2009/03/27 00:53:24 draeger1 Exp $

#ifndef SAMPLEWRITER_H
#define SAMPLEWRITER_H

#include "Context.h"
class Sample;

class SampleWriter
{
  private:

  const Context& ctxt_;

  public:

  SampleWriter(const Context& ctxt);
  void writeSample(const Sample& s, const std::string filename,
                   std::string description,
                   bool base64, bool atomsonly, bool serial);
};

class SampleWriterException
{
  public:
  std::string msg;
  SampleWriterException(std::string s) : msg(s) {}
};

#endif
