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
// StructuredDocumentHandler.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef STRUCTUREDDOCUMENTHANDLER_H
#define STRUCTUREDDOCUMENTHANDLER_H

#include <xercesc/sax2/DefaultHandler.hpp>
#include "StrX.h"

#include "StructureHandler.h"

#include <stack>
#include <string>

class StructuredDocumentHandler : public DefaultHandler
{
  struct HandlerContext
  {
    int depth;
    StructureHandler* handler;
    HandlerContext(StructureHandler* handler_, int depth_) :
      handler(handler_), depth(depth_) {}
  };

  protected:

  std::stack<HandlerContext> contextStack;
  int nestingDepth;
  int contextDepth;
  StructureHandler* activeHandler;
  std::string buffer;

  public:

  StructuredDocumentHandler(StructureHandler* handler) :
    activeHandler(handler), contextDepth(0), nestingDepth(0) {}

  ~StructuredDocumentHandler() {}

  // -----------------------------------------------------------------------
  //  Implementations of the SAX DocumentHandler interface
  // -----------------------------------------------------------------------
  void startDocument();
  void endDocument();

  void startElement(const XMLCh* const uri,const XMLCh* const localname,
    const XMLCh* const qname, const Attributes& attributes);
  void characters(const XMLCh* const chars, const unsigned int length);
  void endElement(const XMLCh* const uri, const XMLCh* const localname,
                  const XMLCh* const qname);
  void ignorableWhitespace(const XMLCh* const chars,
    const unsigned int length);
  void processingInstruction(const XMLCh* const target,
       const XMLCh* const data);

  // -----------------------------------------------------------------------
  //  Implementations of the SAX ErrorHandler interface
  // -----------------------------------------------------------------------
  void warning(const SAXParseException& exception);
  void error(const SAXParseException& exception);
  void fatalError(const SAXParseException& exception);

  // -----------------------------------------------------------------------
  //  Implementation of the SAX DTDHandler interface
  // -----------------------------------------------------------------------
  void notationDecl(const XMLCh* const name, const XMLCh* const publicId,
       const XMLCh* const systemId);

  void unparsedEntityDecl(const XMLCh* const name,
    const XMLCh* const publicId, const XMLCh* const systemId,
    const XMLCh* const notationName);

};
#endif
