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
// AtomSet.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef ATOMSET_H
#define ATOMSET_H

#include "Context.h"
#include "Atom.h"
#include "Species.h"
#include "MMSpecies.h"
#include "EmpiricalPotential.h"
#include "SymmetrySet.h"
#include "UnitCell.h"
#include <vector>
#include <string>
#include <list>
#include <map>
#include <string>

class Species;

class AtomSet
{
  private:

  const Context& ctxt_;

  int nel_;
  std::map<std::string,int> na_;  // na_[sp_name]: number of at. of spec. sp_name
  std::map<std::string,int> isp_; // isp_[sp_name]: index of species sp_name
  std::map<std::string,int> ispmm_; // ispmm_[sp_name]: index of MM species sp_name
  std::map<std::string,int> is_; // is_[atom_name]: is index of atom atom_name
  std::map<std::string,int> ia_; // ia_[atom_name]: ia index of atom atom_name
  std::vector<std::string> spname; // spname[is]: name of species is
  vector<vector<vector<int > > > symatomid_;  // symatomid_[is][ia][isym]: index of atom of species
                                             // is from applying symmetry operation isym to atom ia
  UnitCell cell_;
  bool add_fion_ext_;
  vector<vector<double> > fion_ext_;

  public:

  AtomSet(const Context& ctxt) : ctxt_(ctxt), nel_(0) { add_fion_ext_ = false; }
  ~AtomSet(void);

  std::vector<vector<Atom *> > atom_list; // atom_list[is][ia]
  std::vector<vector<Atom *> > mmatom_list;   // mmatom_list[is][ia]
  std::vector<Species *> species_list;    // species_list[is]
  std::vector<MMSpecies *> mmspecies_list;    // mmspecies_list[is]
  std::vector<EmpiricalPotential *> empirical_list; 
  std::vector<vector<int> > usloc_atind;
  std::vector<vector<int> > usloc_atind_t;
  std::vector<int> usloc_nat;
  std::vector<int> usloc_nat_t;
  std::vector<int> naloc_max;
  std::vector<int> naloc_max_t;
  
  const Context& context(void) const { return ctxt_; }
  bool addAtom(Atom *a);
  bool addMMAtom(Atom *a);
  bool delAtom(std::string name);
  bool delMMAtom(std::string name);
  bool addSpecies(Species *sp, std::string name);
  bool addSpecies(Species *sp, std::string name, const double rcpsin);
  bool addMMSpecies(MMSpecies *sp, std::string name);
  bool delSpecies(std::string name);
  bool delMMSpecies(std::string name);
  bool addEmpiricalPotential(EmpiricalPotential *ep);
  bool reset(void); // remove all atoms and species
  Atom *findAtom(std::string name) const;
  Atom *findMMAtom(std::string name) const;
  Species *findSpecies(std::string name) const;
  MMSpecies *findMMSpecies(std::string name) const;
  void listAtoms(void) const;
  void listSpecies(void) const;
  int na(const std::string &spname) const;  // number of atoms of species spname
  int na(int is) const;         // number of atoms of species is
  int na_mm(int is) const;         // number of MM atoms of species is
  int isp(std::string &spname) const; // index of species spname
  int is(const std::string &atom_name) const; // index of atom atom_name
  int ia(const std::string &atom_name) const; // index of atom atom_name
  int isp_mm(std::string spname) const; // index of MM species spname
  int nel(void) const { return nel_; };
  int nsp(void) const { return species_list.size(); }
  int nsp_mm(void) const { return mmspecies_list.size(); }
  double mass(int isp);
  double mass_mm(int isp);
  int atomic_number(int isp) const;

  void get_positions(vector<vector<double> >& tau, bool qmonly) const;
  void get_positions(vector<vector<double> >& tau) const;
  void set_positions(const vector<vector<double> >& tau, bool ignorelock);
  void set_positions(const vector<vector<double> >& tau);
  void get_velocities(vector<vector<double> >& vel, bool qmonly) const;
  void get_velocities(vector<vector<double> >& vel) const;
  void set_velocities(const vector<vector<double> >& vel);
  bool add_fion_ext(void) { return add_fion_ext_; }
  void get_fion_ext(vector<vector<double> >& fion_ext) const;
  void set_fion_ext(vector<vector<double> >& fion_ext);
  D3vector get_fion_ext(int is, int ia);
  void findSymmetricAtoms(const SymmetrySet& symset);
  void set_rcps(const double& ecut) const;
  int symatomid(int is, int ia, int isym) { return symatomid_[is][ia][isym];};
  void printsys(ostream& os) const;
  void print_casino(ostream& os) const;

  const UnitCell& cell(void) const { return cell_; }
  void set_cell(const UnitCell& cell) { cell_ = cell; }
  void set_cell(const D3vector& a, const D3vector& b, const D3vector& c)
  {
    cell_.set(a,b,c);
  }
  void sync(void);
  void reset_velocities(void);
  void randomize_velocities(double temp);
  D3vector vcm(void) const;
  D3vector dipole(void) const;
  void reset_vcm(void);
  void fold_in_ws(void);
  int size(void) const;
 };
std::ostream& operator << ( std::ostream &os, const AtomSet &as );
#endif

// Local Variables:
// mode: c++
// End:
