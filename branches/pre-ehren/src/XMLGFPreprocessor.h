////////////////////////////////////////////////////////////////////////////////
//
// XMLGFPreprocessor.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: XMLGFPreprocessor.h,v 1.3 2009/03/27 00:53:24 draeger1 Exp $

#include <string>
#include "Matrix.h"

////////////////////////////////////////////////////////////////////////////////
//
// XMLGFPreprocessor class
// Used to preprocess <grid_function> elements in an XML document
// Input: filename, DoubleMatrix& gfdata, string xmlcontent
// The preprocessor reads the file in parallel, processes all <grid_function>
// elements and stores the values of the grid_functions in the matrix gfdata
// which has dimensions (ngf,maxgridsize), where ngf is the total number of
// <grid_function> elements found in the file, and maxgridsize is the size of
// the largest grid_function.
// On return, the string xmlcontent contains (on task 0) the XML file
// with <grid_function> elements reduced to empty strings.
//
////////////////////////////////////////////////////////////////////////////////
class XMLGFPreprocessor
{
  public:

  void process(const char* const filename,
    DoubleMatrix& gfdata, std::string& xmlcontent);
};
