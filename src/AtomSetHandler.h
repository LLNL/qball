////////////////////////////////////////////////////////////////////////////////
//
// AtomSetHandler.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: AtomSetHandler.h,v 1.3 2009/03/27 00:53:24 draeger1 Exp $

#ifndef ATOMSETHANDLER_H
#define ATOMSETHANDLER_H

#include "StructureHandler.h"
#include "D3vector.h"
#include <string>

class AtomSet;
class Species;

class AtomSetHandler : public StructureHandler
{
  private:

  AtomSet& as_;
  std::string current_atom_name, current_atom_species;
  std::string current_species_name;
  D3vector current_atom_position, current_atom_velocity;
  Species* current_species;

  public:

  // Start of the root element in the structure being handled
  virtual void startElement(const XMLCh* const uri,const XMLCh* const localname,
      const XMLCh* const qname, const Attributes& attributes);

  // End of the root element in the structure being handled
  virtual void endElement(const XMLCh* const uri, const XMLCh* const localname,
      const XMLCh* const qname, std::string& content);

  // start a subhandler
  virtual StructureHandler* startSubHandler(const XMLCh* const uri,
    const XMLCh* const localname, const XMLCh* const qname,
    const Attributes& attributes);

  // end a subhandler
  virtual void endSubHandler(const XMLCh* const uri,
    const XMLCh* const localname, const XMLCh* const qname,
    const StructureHandler* const subHandler);

  AtomSetHandler(AtomSet& as);
  ~AtomSetHandler();
};
#endif
