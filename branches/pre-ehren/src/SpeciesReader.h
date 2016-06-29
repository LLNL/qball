////////////////////////////////////////////////////////////////////////////////
//
// SpeciesReader.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SpeciesReader.h,v 1.1.1.1 2005/08/18 17:23:33 draeger1 Exp $

#ifndef SPECIESREADER_H
#define SPECIESREADER_H

#include <string>
using namespace std;
#include "Context.h"

class SpeciesReader
{
  private:
  
  const Context& ctxt_;
  
  string uri_;   // uri from which Species is read
  
  public:

  SpeciesReader(const Context& ctxt);
  void readSpecies(Species& sp, const string uri);
  void bcastSpecies(Species& sp);
};

class SpeciesReaderException {};

#endif
