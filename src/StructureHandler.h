////////////////////////////////////////////////////////////////////////////////
//
// StructureHandler.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: StructureHandler.h,v 1.3 2009/03/27 00:53:24 draeger1 Exp $

#ifndef STRUCTUREHANDLER_H
#define STRUCTUREHANDLER_H

#include <xercesc/util/XMLUniDefs.hpp>
#include <xercesc/sax2/Attributes.hpp>
#include <string>
using namespace xercesc;

class StructureHandler
{
  private:

  public:

  // Start of an element handled by the StructureHandler
  virtual void startElement(const XMLCh* const uri,const XMLCh* const localname,
      const XMLCh* const qname, const Attributes& attributes) = 0;

  // End of an element handled by the StructureHandler
  virtual void endElement(const XMLCh* const uri, const XMLCh* const localname,
      const XMLCh* const qname, std::string& content) = 0;

  // start a subhandler
  virtual StructureHandler* startSubHandler(const XMLCh* const uri,
    const XMLCh* const localname, const XMLCh* const qname,
    const Attributes& attributes) = 0;

  // end a subhandler
  virtual void endSubHandler(const XMLCh* const uri,
    const XMLCh* const localname, const XMLCh* const qname,
    const StructureHandler* const subHandler) = 0;

  StructureHandler() {}
  virtual ~StructureHandler() {}
};
#endif
