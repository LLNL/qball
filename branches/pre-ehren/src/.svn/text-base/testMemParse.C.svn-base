////////////////////////////////////////////////////////////////////////////////
//
// testMemParse.C:
//
// test in-memory parsing using progressive parse
//
////////////////////////////////////////////////////////////////////////////////

#include <cassert>
#include <string>
#include <iostream>
#include <algorithm>
using namespace std;

#include <xercesc/util/XMLUniDefs.hpp>
#include <xercesc/sax2/Attributes.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>
#include "StrX.h"

#include <xercesc/sax2/DefaultHandler.hpp>
#include <xercesc/framework/MemBufInputSource.hpp>
using namespace xercesc;

static const char*  input_buffer[2] = {
"<?xml version='1.0' encoding='utf-8'?>\n\
<!DOCTYPE company [\n\
<!ELEMENT company     (product,category,developedAt)>\n\
<!ELEMENT product     (#PCDATA)>\n\
<!ELEMENT category    (#PCDATA)>\n\
<!ATTLIST category idea CDATA #IMPLIED>\n\
<!ELEMENT developedAt (#PCDATA)>\n\
]>\n\n\
<company>\n\
    <product>XML4C</product>\n\
\0", 


"\
    <category idea='great'>XML Parsing Tools</category>\n\
    <developedAt>\n\
</company>\
\0" };

// static const char*  gXMLInMemBuf =
// "\
// <?xml version='1.0' encoding='utf-8'?>\n\
// <!DOCTYPE company [\n\
// <!ELEMENT company     (product,category,developedAt)>\n\
// <!ELEMENT product     (#PCDATA)>\n\
// <!ELEMENT category    (#PCDATA)>\n\
// <!ATTLIST category idea CDATA #IMPLIED>\n\
// <!ELEMENT developedAt (#PCDATA)>\n\
// ]>\n\n\
// <company>\n\
//     <product>XML4C</product>\n\
//     <category idea='great'>XML Parsing Tools</category>\n\
//     <developedAt>\n\
//       IBM Center for Java Technology, Silicon Valley, Cupertino, CA\n\
//     </developedAt>\n\
// </company>\
// ";
//

static const char*  gMemBufId = "prodInfo";

////////////////////////////////////////////////////////////////////////////////
class TestHandler : public DefaultHandler
{
  public:

  string buffer;

  TestHandler(void) { buffer = ""; cout << "TestHandler::ctor" << endl; }
    
  ~TestHandler() { cout << "TestHandler::dtor" << endl; }
  
  // -----------------------------------------------------------------------
  //  Implementations of the SAX DocumentHandler interface
  // -----------------------------------------------------------------------
  void startDocument() { cout << "TestHandler::startDocument" << endl; };
  void endDocument() { cout << "TestHandler::endDocument" << endl; };

  void startElement(const XMLCh* const uri,const XMLCh* const localname,
    const XMLCh* const qname, const Attributes& attributes)
  { 
    cout << "TestHandler::startElement: " << StrX(qname) << endl;
    buffer += string(StrX(qname).localForm());
  };
  void characters(const XMLCh* const chars, const unsigned int length)
  { 
    // cout << "TestHandler::characters" << endl;
    buffer += string(XMLString::transcode(chars),length);
  };
  void endElement(const XMLCh* const uri, const XMLCh* const localname,
                  const XMLCh* const qname) 
  { 
    cout << "TestHandler::endElement:   " << StrX(qname) << endl;
    buffer += string(StrX(qname).localForm());
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
  //SAX2XMLReader::ValSchemes valScheme = SAX2XMLReader::Val_Auto;
  SAX2XMLReader::ValSchemes valScheme = SAX2XMLReader::Val_Always;
  bool expandNamespaces = false;
  bool doNamespaces = true;
  bool doSchema = true;
  bool schemaFullChecking = true;
  bool namespacePrefixes = false;
  SAX2XMLReader* parser = 0;
  
  string buffer;
  buffer.resize(1024);
  fill(buffer.begin(),buffer.end(),' ');
  
  buffer = input_buffer[0];
  int ibuf = 0;
  const int nbuf = 1;
  
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
 
  MemBufInputSource* memBufIS = new MemBufInputSource
  (
      (const XMLByte*) &buffer[0]
      , buffer.size()
      , gMemBufId
      , false
  );

  try
  {
    TestHandler handler;
    parser->setContentHandler(&handler);
    parser->setErrorHandler(&handler);

    // parser->parse(*memBufIS);
    if ( !parser->parseFirst(*memBufIS,token) )
    {
      cout << "parseFirst() failed" << endl;
      XMLPlatformUtils::Terminate();
      exit(1);
    }
    cout << "parseFirst() completed" << endl;

    bool gotMore = true;
    bool done = false;
    while ( gotMore && !done )
    {
      gotMore = parser->parseNext(token);
      if ( !gotMore )
      {
        cout << "parseNext() returned 0" << endl;
        // reload buffer
        if ( ibuf < nbuf-1 )
        {
          buffer = input_buffer[ibuf];
          ibuf++;
          gotMore = true;
        }

      }
      cout << "buffer = " << handler.buffer << endl;
      if ( handler.buffer.size() > 32 )
        handler.buffer = "";
    }

    errorCount = parser->getErrorCount();
    
    cout << handler.buffer << endl;
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
  delete memBufIS;
  XMLPlatformUtils::Terminate();
}
