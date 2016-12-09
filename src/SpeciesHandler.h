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
// SpeciesHandler.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef SPECIESHANDLER_H
#define SPECIESHANDLER_H

#include "StructureHandler.h"

class Species;

class SpeciesHandler : public StructureHandler
{
  private:

  Species& sp_;
  int current_l, current_size;
  std::string current_name, current_href;
  double current_interval;

  public:

  // Start of an element handled by SpeciesHandler
  virtual void startElement(const XMLCh* const uri,const XMLCh* const localname,
      const XMLCh* const qname, const Attributes& attributes);

  // End of an element handled by SpeciesHandler
  virtual void endElement(const XMLCh* const uri, const XMLCh* const localname,
      const XMLCh* const qname, std::string& content);

  // start a subHandler if possible
  virtual StructureHandler* startSubHandler(const XMLCh* const uri,
    const XMLCh* const localname, const XMLCh* const qname,
    const Attributes& attributes);

  // end subHandler
  virtual void endSubHandler(const XMLCh* const uri,
    const XMLCh* const localname, const XMLCh* const qname,
    const StructureHandler* const subHandler);

  SpeciesHandler(Species& sp);
  ~SpeciesHandler();
};
#endif

// Local Variables:
// mode: c++
// End:
