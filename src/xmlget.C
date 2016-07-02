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
// xmlget.C:
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#include <cassert>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
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
  bool read;
  vector<string> y;
  string tag;
  public:

  string buffer;

  TestHandler(string tag_) : tag(tag_), read(false) { buffer = ""; }
    
  ~TestHandler() {}
  
  // -----------------------------------------------------------------------
  //  Implementations of the SAX DocumentHandler interface
  // -----------------------------------------------------------------------
  void startDocument() {}
  void endDocument() {}
  void print_data()
  { 
    //cout << "TestHandler::endDocument" << endl;
    cout << y.size() << endl;
    cout.setf(ios::scientific,ios::floatfield);
    for ( int i = 0; i < y.size(); i++ )
      cout << i << " " << y[i] << endl;
  };

  void startElement(const XMLCh* const uri,const XMLCh* const localname,
    const XMLCh* const qname, const Attributes& attributes)
  { 
    //cout << "TestHandler::startElement: " << StrX(qname) << endl;
    char* buf = XMLString::transcode(qname);
    if ( buf == tag )
    {
      buffer = "";
      read = true;
    }
    XMLString::release(&buf);
  };
  void characters(const XMLCh* const chars, const unsigned int length)
  { 
    // cout << "TestHandler::characters" << endl;
    if ( read )
    {
      char* buf = XMLString::transcode(chars);
      buffer += string(buf,length);
      XMLString::release(&buf);
    }
  };
  void endElement(const XMLCh* const uri, const XMLCh* const localname,
                  const XMLCh* const qname) 
  { 
    //cout << "TestHandler::endElement:   " << StrX(qname) << endl;
    if ( XMLString::transcode(qname) == tag )
    {
      y.push_back(buffer);
      read = false;
    }
  }
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
  SAX2XMLReader::ValSchemes valScheme = SAX2XMLReader::Val_Never;
  //SAX2XMLReader::ValSchemes valScheme = SAX2XMLReader::Val_Auto;
  //SAX2XMLReader::ValSchemes valScheme = SAX2XMLReader::Val_Always;
  bool expandNamespaces = false;
//   bool doNamespaces = true;
//   bool doSchema = true;
//   bool schemaFullChecking = true;
  bool doNamespaces = false;
  bool doSchema = false;
  bool schemaFullChecking = false;
  bool namespacePrefixes = false;
  SAX2XMLReader* parser = 0;
  assert(argc==3);
  string tag(argv[1]);
  
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
  XMLPScanToken token;
 
  try
  {
    TestHandler handler(tag);
    parser->setContentHandler(&handler);
    parser->setErrorHandler(&handler);

    if ( !parser->parseFirst(argv[2],token) )
    {
      cout << "parseFirst() failed" << endl;
      XMLPlatformUtils::Terminate();
      exit(1);
    }
    //cout << "parseFirst() completed" << endl;

    bool gotMore = true;
    bool done = false;
    while ( gotMore && !done )
    {
      gotMore = parser->parseNext(token);
      if ( !gotMore )
        //cout << "parseNext() returned 0" << endl;
      //cout << "buffer = " << handler.buffer << endl;
      if ( handler.buffer.size() > 32 )
        handler.buffer = "";
    }

    errorCount = parser->getErrorCount();
    
    //cout << handler.buffer << endl;
    
    handler.print_data();
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
