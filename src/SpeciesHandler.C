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
// SpeciesHandler.C
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#if HAVE_XERCES

#include "SpeciesHandler.h"
#include "Species.h"
#include "StrX.h"
using namespace xercesc;
#include <iostream>
#include <sstream>
#include <cassert>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
SpeciesHandler::SpeciesHandler(Species& sp) :
  sp_(sp) {}

////////////////////////////////////////////////////////////////////////////////
SpeciesHandler::~SpeciesHandler() {}

////////////////////////////////////////////////////////////////////////////////
void SpeciesHandler::startElement(const XMLCh* const uri,
  const XMLCh* const localname, const XMLCh* const qname,
  const Attributes& attributes)
{
  // cout << " SpeciesHandler::startElement " << StrX(qname) << endl;

  string locname(XMLString::transcode(localname));

  if ( locname == "species" )
  {
    // check for the case where the species is a link to another uri
    unsigned int len = attributes.getLength();
    for (unsigned int index = 0; index < len; index++)
    {
      string attrname(XMLString::transcode(attributes.getLocalName(index)));
      if ( attrname == "name")
      {
        current_name = StrX(attributes.getValue(index)).localForm();
        sp_.name_ = current_name;
      }
      else if ( attrname == "href" )
      {
        current_href = StrX(attributes.getValue(index)).localForm();
        cout << " SpeciesHandler: found href in species definition" << endl
           << " name=" << current_name << " href=" << current_href
           << endl;
        sp_.uri_ = current_href;
      }
    }
  }
  else if ( locname == "projector" )
  {
    unsigned int len = attributes.getLength();
    for (unsigned int index = 0; index < len; index++)
    {
      string attrname(XMLString::transcode(attributes.getLocalName(index)));
      if ( attrname == "l")
      {
        current_l = atoi(StrX(attributes.getValue(index)).localForm());
      }
      else if ( attrname == "size" )
      {
        current_size = atoi(StrX(attributes.getValue(index)).localForm());
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void SpeciesHandler::endElement(const XMLCh* const uri,
  const XMLCh* const localname, const XMLCh* const qname, string& content)
{
  string locname(XMLString::transcode(localname));
  istringstream stst(content);
  // cout << " SpeciesHandler::endElement " << StrX(qname)
  //      << " content=" << string(content,0,20) << endl;

  if ( locname == "description")
  {
    // reject ambiguous case where both the href and the definition are given
    if ( current_href != "" )
    {
      cout << " SpeciesHandler: ambiguous definition: uri="
           << StrX(uri) << endl
           << " using local definition (href: " << current_href << " ignored)"
           << endl;
    }
    sp_.description_ = content;
  }
  else if ( locname == "atomic_number" )
  {
    stst >> sp_.atomic_number_;
  }
  else if ( locname == "mass" )
  {
    stst >> sp_.mass_;
  }
  else if ( locname == "symbol" )
  {
    stst >> skipws >> sp_.symbol_;
  }
  else if ( locname == "valence_charge" )
  {
    stst >> sp_.zval_;
  }
  else if ( locname == "lmax" )
  {
    stst >> sp_.lmax_;
  }
  else if ( locname == "llocal" )
  {
    stst >> sp_.llocal_;
  }
  else if ( locname == "nquad" )
  {
    stst >> sp_.nquad_;
  }
  else if ( locname == "rquad" )
  {
    stst >> sp_.rquad_;
  }
  else if ( locname == "mesh_spacing" )
  {
    stst >> sp_.deltar_;
  }
  else if ( locname == "radial_potential" )
  {
    if ( current_l+1  > sp_.vps_.size() )
    {
      sp_.vps_.resize(current_l+1);
      sp_.phi_.resize(current_l+1);
    }
    sp_.vps_[current_l].resize(current_size);
    for ( int i = 0; i < current_size; i++ )
      stst >> sp_.vps_[current_l][i];
  }
  else if ( locname == "radial_function" )
  {
    sp_.phi_[current_l].resize(current_size);
    for ( int i = 0; i < current_size; i++ )
      stst >> sp_.phi_[current_l][i];
  }
}

////////////////////////////////////////////////////////////////////////////////
StructureHandler* SpeciesHandler::startSubHandler(const XMLCh* const uri,
    const XMLCh* const localname, const XMLCh* const qname,
    const Attributes& attributes)
{
  // check if element qname can be processed by another StructureHandler
  // If it can, return a pointer to the StructureHandler, otherwise return 0
  // cout << " SpeciesHandler::startSubHandler " << StrX(qname) << endl;
  return 0;
}

////////////////////////////////////////////////////////////////////////////////
void SpeciesHandler::endSubHandler(const XMLCh* const uri,
    const XMLCh* const localname, const XMLCh* const qname,
    const StructureHandler* const subHandler)
{
  // cout << " SpeciesHandler::endSubHandler " << StrX(qname) << endl;
  // if any StructureHandler was created by startSubHandler, delete it
  // delete subHandler;
}

#endif
