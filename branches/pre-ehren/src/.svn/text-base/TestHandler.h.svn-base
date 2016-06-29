////////////////////////////////////////////////////////////////////////////////
//
// TestHandler.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: TestHandler.h,v 1.1.1.1 2005/08/18 17:23:33 draeger1 Exp $

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
