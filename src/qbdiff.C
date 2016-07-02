////////////////////////////////////////////////////////////////////////////////  
// Copyright (c) 2016, Lawrence Livermore National Security, LLC. 
// qb@ll:  Qbox at Lawrence Livermore
//
// This file is part of qb@ll.
//
// Produced at the Lawrence Livermore National Laboratory. 
// Written by Xavier Andrade (xavier@llnl.gov).
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
// qbdiff.C:
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#include <iostream>
#include <string>
#include <cmath>

#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/util/PlatformUtils.hpp>

using namespace std;
using namespace xercesc;

double energy_tol = 1.0e-6;

int main(int argc, char** argv)
{

  if(argc != 3){
    cerr << "Usage: " << argv[0] << " file1.xml file2.xml" << endl;
  }

  string filename1 = string(argv[1]);
  string filename2 = string(argv[2]);

  XMLPlatformUtils::Initialize();

  XercesDOMParser * parser = new XercesDOMParser();

  parser->parse(filename1.c_str());
  DOMDocument * doc1 = parser->adoptDocument();

  parser->parse(filename2.c_str());
  DOMDocument * doc2 = parser->adoptDocument();

  XMLCh* tmp = XMLString::transcode("iteration");
  DOMNodeList* list1 = doc1->getElementsByTagName(tmp);
  DOMNodeList* list2 = doc2->getElementsByTagName(tmp);
  XMLString::release(&tmp);

  if(list1->getLength() != list2->getLength()){
    cout << "Different number of iterations: " << list1->getLength() << " vs " <<list2->getLength() << std::endl;
    exit(1);
  }
  
  bool correct = true;

  for(unsigned iter = 0; iter < list1->getLength(); iter++){
    
    bool iter_correct = true;

    cout << "Iteration\t" << iter << "\t:\t";

    DOMElement * it1 = dynamic_cast<DOMElement*>(list1->item(iter));
    DOMElement * it2 = dynamic_cast<DOMElement*>(list2->item(iter));
    
    XMLCh* tmp = XMLString::transcode("etotal");
    DOMElement* child1 = dynamic_cast<DOMElement*>(it1->getElementsByTagName(tmp)->item(0));
    DOMElement* child2 = dynamic_cast<DOMElement*>(it2->getElementsByTagName(tmp)->item(0));
    XMLString::release(&tmp);

    double etotal1 = strtod(XMLString::transcode(child1->getTextContent()), NULL);
    double etotal2 = strtod(XMLString::transcode(child2->getTextContent()), NULL);
    
    //    cout << etotal1 << " " << etotal2 << " " << fabs(etotal1 - etotal2) << endl;

    if(fabs(etotal1 - etotal2) > energy_tol){
      iter_correct = false;
    }

    correct = correct && iter_correct;

    if(iter_correct){
      cout << "[   OK   ]" << endl;
    } else {
      cout << "[  FAIL  ]" << endl;
    }

  }

  doc1->release();
  doc2->release();
  
  XMLPlatformUtils::Terminate();

  if(!correct) return 1;
}
