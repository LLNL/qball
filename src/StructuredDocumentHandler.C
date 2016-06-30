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
// StructuredDocumentHandler.C
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#if USE_XERCES

#include "StructuredDocumentHandler.h"
#include <xercesc/util/XMLUniDefs.hpp>
#include <xercesc/sax2/Attributes.hpp>
#if TIMING
#include "Timer.h"
#endif

#include "StrX.h"
#include <iostream>
#include <cassert>
using namespace xercesc;
using namespace std;

////////////////////////////////////////////////////////////////////////////////
void StructuredDocumentHandler::startElement(const XMLCh* const uri,
  const XMLCh* const localname, const XMLCh* const qname,
  const Attributes& attributes)
{
  // cout << " StructuredDocumentHandler::startElement " << StrX(qname) << endl;
  // cout << " nestingDepth before: " << nestingDepth << endl;
  buffer = "";

  // attempt to start a subhandler for this element
  StructureHandler* next = activeHandler->startSubHandler(uri,localname,qname,
    attributes);

  // yes, this element can be processed by a subhandler
  if ( next != 0 )
  {
    contextStack.push(HandlerContext(activeHandler,contextDepth));

    activeHandler = next;
    contextDepth = nestingDepth;
  }
  activeHandler->startElement(uri,localname,qname,attributes);
  nestingDepth++;
  // cout << " nestingDepth after:  " << nestingDepth << endl;
}

////////////////////////////////////////////////////////////////////////////////
#if XERCES_VERSION_MAJOR >= 3
void StructuredDocumentHandler::characters(const XMLCh* const chars,
                                           const XMLSize_t length)
#else
void StructuredDocumentHandler::characters(const XMLCh* const chars,
  const unsigned int length)
#endif
{
#if TIMING
  Timer tm;
  tm.start();
#endif
  size_t pos = buffer.size();
  // Note: buffer must be able to hold length+1 chars for '\0'
  buffer.resize(pos+length+1);
  bool status = XMLString::transcode(chars,&buffer[pos],length);
  assert(status==true);
#if TIMING
  tm.stop();
  cout << " StructuredDocumentHandler::characters: time: " << tm.real() << endl;
#endif
  //cout << "length=" << length << " buffer.size=" << buffer.size() << endl;
}

////////////////////////////////////////////////////////////////////////////////
void StructuredDocumentHandler::endElement(const XMLCh* const uri,
  const XMLCh* const localname, const XMLCh* const qname)
{
  // cout << " StructuredDocumentHandler::endElement " << StrX(qname) << endl;
  // cout << " nestingDepth before: " << nestingDepth << endl;
  nestingDepth--;
  activeHandler->endElement(uri,localname,qname,buffer);

  // Check if this element was processed by a subhandler
  // If yes, pop the previous context from the Context stack
  if ( nestingDepth != 0 && nestingDepth == contextDepth )
  {
    // cout << " popping context stack: nestingDepth=" << nestingDepth << endl;
    HandlerContext context = contextStack.top();
    StructureHandler* last = activeHandler;
    activeHandler = context.handler;
    contextDepth = context.depth;
    contextStack.pop();

    // notify activeHandler that current structure is complete
    // This allows the activeHandler to delete the last handler it has created
    activeHandler->endSubHandler(uri,localname,qname,last);
  }
  // cout << " nestingDepth after:  " << nestingDepth << endl;
  buffer = "";
}

////////////////////////////////////////////////////////////////////////////////
void StructuredDocumentHandler::startDocument() {}

////////////////////////////////////////////////////////////////////////////////
void StructuredDocumentHandler::endDocument() {}

////////////////////////////////////////////////////////////////////////////////
void StructuredDocumentHandler::ignorableWhitespace( const   XMLCh* const chars,
  const  unsigned int length) {}

////////////////////////////////////////////////////////////////////////////////
void StructuredDocumentHandler::processingInstruction(const  XMLCh* const target,
  const XMLCh* const data) {}

////////////////////////////////////////////////////////////////////////////////
void StructuredDocumentHandler::error(const SAXParseException& e)
{
    cout << "\nError at file " << StrX(e.getSystemId())
		 << ", line " << e.getLineNumber()
		 << ", char " << e.getColumnNumber()
         << "\n  Message: " << StrX(e.getMessage()) << endl;

    throw(e);
}

////////////////////////////////////////////////////////////////////////////////
void StructuredDocumentHandler::fatalError(const SAXParseException& e)
{
    cout << "\nFatal Error at file " << StrX(e.getSystemId())
		 << ", line " << e.getLineNumber()
		 << ", char " << e.getColumnNumber()
         << "\n  Message: " << StrX(e.getMessage()) << endl;
    throw(e);
}

////////////////////////////////////////////////////////////////////////////////
void StructuredDocumentHandler::warning(const SAXParseException& e)
{
    cout << "\nWarning at file " << StrX(e.getSystemId())
		 << ", line " << e.getLineNumber()
		 << ", char " << e.getColumnNumber()
         << "\n  Message: " << StrX(e.getMessage()) << endl;
}

////////////////////////////////////////////////////////////////////////////////
void StructuredDocumentHandler::unparsedEntityDecl(const     XMLCh* const name,
  const   XMLCh* const publicId, const   XMLCh* const systemId,
  const   XMLCh* const notationName)
{
    // Not used at this time
}

////////////////////////////////////////////////////////////////////////////////
void StructuredDocumentHandler::notationDecl(const   XMLCh* const name,
  const XMLCh* const publicId, const XMLCh* const systemId)
{
    // Not used at this time
}

#endif
