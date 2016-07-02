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
// TestHandler.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef TESTHANDLER_H
#define TESTHANDLER_H

#include <xercesc/sax2/DefaultHandler.hpp>
#include "StrX.h"
using namespace xercesc;

#include <string>
using namespace std;

class TestHandler : public DefaultHandler
{
  protected:
  
  public:

  TestHandler(void) {}
    
  ~TestHandler() {}
  
  // -----------------------------------------------------------------------
  //  Implementations of the SAX DocumentHandler interface
  // -----------------------------------------------------------------------
  void startDocument() {};
  void endDocument() {};

  void startElement(const XMLCh* const uri,const XMLCh* const localname,
    const XMLCh* const qname, const Attributes& attributes) {};
  void characters(const XMLCh* const chars, const unsigned int length) {};
  void endElement(const XMLCh* const uri, const XMLCh* const localname,
                  const XMLCh* const qname) {};
  void ignorableWhitespace(const XMLCh* const chars,
    const unsigned int length) {};
  void processingInstruction(const XMLCh* const target,
       const XMLCh* const data) {};
};
#endif
