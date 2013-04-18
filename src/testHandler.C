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
// testHandler.C:
//
////////////////////////////////////////////////////////////////////////////////

#include <cassert>
#include <string>
#include <iostream>
#include <vector>
using namespace std;

#include <xercesc/util/XMLUniDefs.hpp>
#include <xercesc/sax2/Attributes.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>
#include "StrX.h"

#include <xercesc/sax2/DefaultHandler.hpp>
using namespace xercesc;

////////////////////////////////////////////////////////////////////////////////
class TestHandler : public DefaultHandler
{
  public:

  TestHandler(void) { cout << "TestHandler::ctor" << endl; }
    
  ~TestHandler() { cout << "TestHandler::dtor" << endl; }
  
  // -----------------------------------------------------------------------
  //  Implementations of the SAX DocumentHandler interface
  // -----------------------------------------------------------------------
  void startDocument() { cout << "TestHandler::startDocument" << endl; };
  void endDocument() { cout << "TestHandler::endDocument" << endl; };

  void startElement(const XMLCh* const uri,const XMLCh* const localname,
    const XMLCh* const qname, const Attributes& attributes)
  { cout << "TestHandler::startElement: " << StrX(qname) << endl; };
  void characters(const XMLCh* const chars, const unsigned int length)
  { 
    // cout << "TestHandler::characters" << endl;
  };
  void endElement(const XMLCh* const uri, const XMLCh* const localname,
                  const XMLCh* const qname) 
  { cout << "TestHandler::endElement:   " << StrX(qname) << endl; };
  void ignorableWhitespace(const XMLCh* const chars,
    const unsigned int length) {};
  void processingInstruction(const XMLCh* const target,
       const XMLCh* const data) {};
  // -----------------------------------------------------------------------
  //  Implementations of the SAX ErrorHandler interface
  // -----------------------------------------------------------------------
  void warning(const SAXParseException& exception) {};
  void error(const SAXParseException& exception) {};
  void fatalError(const SAXParseException& exception) {};
 
  // -----------------------------------------------------------------------
  //  Implementation of the SAX DTDHandler interface
  // -----------------------------------------------------------------------
  void notationDecl(const XMLCh* const name, const XMLCh* const publicId,
       const XMLCh* const systemId) {};

  void unparsedEntityDecl(const XMLCh* const name,
    const XMLCh* const publicId, const XMLCh* const systemId,
    const XMLCh* const notationName) {};
};

////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  const char* encodingName = "UTF-8";
  //SAX2XMLReader::ValSchemes valScheme = SAX2XMLReader::Val_Auto;
  SAX2XMLReader::ValSchemes valScheme = SAX2XMLReader::Val_Always;
  bool expandNamespaces = false;
  bool doNamespaces = true;
  bool doSchema = true;
  bool schemaFullChecking = true;
  bool namespacePrefixes = false;
  SAX2XMLReader* parser = 0;
  
  try
  {
    XMLPlatformUtils::Initialize();
  }

  catch (const XMLException& toCatch)
  {
    cout << " Error during XML initialization :\n"
         << StrX(toCatch.getMessage()) << endl;
  }
    
  parser = XMLReaderFactory::createXMLReader();
  if (valScheme == SAX2XMLReader::Val_Auto)
  {
      parser->setFeature(XMLUni::fgSAX2CoreValidation, true);
      parser->setFeature(XMLUni::fgXercesDynamic, true);
  }

  if (valScheme == SAX2XMLReader::Val_Never)
  {
      parser->setFeature(XMLUni::fgSAX2CoreValidation, false);
  }

  if (valScheme == SAX2XMLReader::Val_Always)
  {
      parser->setFeature(XMLUni::fgSAX2CoreValidation, true);
      parser->setFeature(XMLUni::fgXercesDynamic, false);
  }

  parser->setFeature(XMLUni::fgSAX2CoreNameSpaces, doNamespaces);
  parser->setFeature(XMLUni::fgXercesSchema, doSchema);
  parser->setFeature(XMLUni::fgXercesSchemaFullChecking, schemaFullChecking);
  parser->setFeature(XMLUni::fgSAX2CoreNameSpacePrefixes, namespacePrefixes);

  int errorCount = 0;
 
  try
  {
    TestHandler handler;
    parser->setContentHandler(&handler);
    parser->setErrorHandler(&handler);
    parser->parse(argv[1]);
    errorCount = parser->getErrorCount();
  }

  catch (const XMLException& toCatch)
  {
    cout << "\nAn error occurred\n  Error: "
         << StrX(toCatch.getMessage())
         << "\n" << endl;
    XMLPlatformUtils::Terminate();
    delete parser;
    throw;
  }

  catch (const SAXParseException& e)
  {
    cout << "\na SAXParseException occurred in file "
         << StrX(e.getSystemId())
         << ", line " << e.getLineNumber()
         << ", char " << e.getColumnNumber()
         << "\n  Message: " << StrX(e.getMessage()) << endl;
    XMLPlatformUtils::Terminate();
    delete parser;
    throw;
  }
 
  delete parser;
  XMLPlatformUtils::Terminate();
}
