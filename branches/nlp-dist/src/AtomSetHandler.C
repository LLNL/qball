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
// AtomSetHandler.C
//
////////////////////////////////////////////////////////////////////////////////

#if USE_XERCES

#include "AtomSetHandler.h"
#include "AtomSet.h"
#include "Species.h"
#include "SpeciesHandler.h"
#include "SpeciesReader.h"
#include "StrX.h"
#include "SampleReader.h"
using namespace xercesc;
#include <iostream>
#include <cassert>
#include <sstream>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
AtomSetHandler::AtomSetHandler(AtomSet& as) :
  as_(as) {}

////////////////////////////////////////////////////////////////////////////////
AtomSetHandler::~AtomSetHandler() {}

////////////////////////////////////////////////////////////////////////////////
void AtomSetHandler::startElement(const XMLCh* const uri,
  const XMLCh* const localname, const XMLCh* const qname,
  const Attributes& attributes)
{
  // cout << " AtomSetHandler::startElement " << StrX(qname) << endl;
  string locname(XMLString::transcode(localname));

  // consider only elements that are dealt with directly by AtomSetHandler
  // i.e. "atom". The "species" element is delegated to a SpeciesHandler
  if ( locname == "unit_cell")
    {
      D3vector a,b,c;
      unsigned int len = attributes.getLength();
      for (unsigned int index = 0; index < len; index++)
        {
          string attrname(XMLString::transcode(attributes.getLocalName(index)));
          string attrval(XMLString::transcode(attributes.getValue(index)));
          istringstream stst(attrval);
          if ( attrname == "a")
            {
              stst >> a;
            }
          else if ( attrname == "b" )
            {
              stst >> b;
            }
          else if ( attrname == "c" )
            {
              stst >> c;
            }
        }

      as_.set_cell(a,b,c);
    }
  else if ( locname == "atom")
  {
    // set default velocity to zero
    current_atom_velocity = D3vector(0.0,0.0,0.0);

    unsigned int len = attributes.getLength();
    for (unsigned int index = 0; index < len; index++)
    {
      string attrname(XMLString::transcode(attributes.getLocalName(index)));
      if ( attrname == "name")
      {
        current_atom_name = StrX(attributes.getValue(index)).localForm();
      }
      else if ( attrname == "species" )
      {
        current_atom_species = StrX(attributes.getValue(index)).localForm();
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void AtomSetHandler::endElement(const XMLCh* const uri,
  const XMLCh* const localname, const XMLCh* const qname, string& content)
{
  string locname(XMLString::transcode(localname));
  // cout << " AtomSetHandler::endElement " << locname << endl;
  istringstream stst(content);
  if ( locname == "atom")
  {
    // create an instance of Atom using the current values read so far
    // add the Atom to the current AtomSet
    // cout << " AtomSetHandler::endElement: creating Atom(position="
    //      << current_atom_position
    //      << ", velocity=" << current_atom_velocity
    //      << ", name=" << current_atom_name
    //      << ", species=" << current_atom_species
    //      << ")" << endl;

    try
    {
      Atom* a = new Atom(current_atom_name, current_atom_species,
                         current_atom_position, current_atom_velocity);
      as_.addAtom(a);
    }
    catch (...)
    {
      cout << " AtomSetHandler::endElement: and exception occurred"
           << " in AtomSet::addAtom" << endl;
      throw;
    }

    // notify listening nodes and broadcast atom info
    event_type event = atom;
    as_.context().ibcast_send(1,1,(int*)&event,1);
    as_.context().string_bcast(current_atom_name,0);
    as_.context().string_bcast(current_atom_species,0);
    double buf[3];
    buf[0] = current_atom_position.x;
    buf[1] = current_atom_position.y;
    buf[2] = current_atom_position.z;
    as_.context().dbcast_send(3,1,buf,3);
    buf[0] = current_atom_velocity.x;
    buf[1] = current_atom_velocity.y;
    buf[2] = current_atom_velocity.z;
    as_.context().dbcast_send(3,1,buf,3);

  }
  else if ( locname == "position" )
  {
    stst >> current_atom_position;
  }
  else if ( locname == "velocity" )
  {
    stst >> current_atom_velocity;
  }
}

////////////////////////////////////////////////////////////////////////////////
StructureHandler* AtomSetHandler::startSubHandler(const XMLCh* const uri,
    const XMLCh* const localname, const XMLCh* const qname,
    const Attributes& attributes)
{
  // check if element qname can be processed by another StructureHandler
  // If it can, return a pointer to the StructureHandler, otherwise return 0
  // cout << " AtomSetHandler::startSubHandler " << StrX(qname) << endl;

  string locname(XMLString::transcode(localname));
  if ( locname == "species")
  {
    // check for species attributes
    unsigned int len = attributes.getLength();
    for (unsigned int index = 0; index < len; index++)
    {
      string attrname(XMLString::transcode(attributes.getLocalName(index)));
      if ( attrname == "name")
      {
        current_species_name = StrX(attributes.getValue(index)).localForm();
      }
    }

    // delegate to SpeciesHandler
    current_species = new Species(as_.context(),current_species_name);
    return new SpeciesHandler(*current_species);
  }
  else
  {
    return 0;
  }
}

////////////////////////////////////////////////////////////////////////////////
void AtomSetHandler::endSubHandler(const XMLCh* const uri,
    const XMLCh* const localname, const XMLCh* const qname,
    const StructureHandler* const last)
{
  string locname(XMLString::transcode(localname));
  // cout << " AtomSetHandler::endSubHandler " << locname << endl;
  if ( locname == "species" )
  {
    SpeciesReader sp_reader(current_species->context());

    // check if only the uri was provided
    if ( current_species->uri() != "" )
    {
      // href was found in species definition
      // attempt to read the species from that uri

      try
      {
        sp_reader.readSpecies(*current_species,current_species->uri());
      }
      catch ( const SpeciesReaderException& e )
      {
        cout << " SpeciesReaderException caught in AtomSetHandler" << endl;
      }
    }

    // notify listening nodes and broadcast species info
    event_type event = species;
    as_.context().ibcast_send(1,1,(int*)&event,1);
    as_.context().string_bcast(current_species_name,0);
    sp_reader.bcastSpecies(*current_species);

    // cout << "AtomSetHandler::endSubHandler: adding Species:"
    //      << current_species_name << endl;
    try
    {
      as_.addSpecies(current_species,current_species_name);
    }
    catch (...)
    {
      cout << " AtomSetHandler::endSubHandler: and exception occurred"
           << " in AtomSet::addSpecies" << endl;
      throw;
    }
  }
  delete last;
}

#endif
