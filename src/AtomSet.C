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
// AtomSet.C
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#include "AtomSet.h"
#include "NameOf.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
AtomSet::~AtomSet(void)
{
  for ( int is = 0; is < species_list.size(); is++ )
  {
    for ( int ia = 0; ia < atom_list[is].size(); ia++ )
    {
      delete atom_list[is][ia];
    }
    delete species_list[is];
  }
}

////////////////////////////////////////////////////////////////////////////////
bool AtomSet::addSpecies(Species* sp, string name)
{
  if ( findSpecies(name) || findMMSpecies(name) )
  {
    if ( ctxt_.oncoutpe() )
      cout << " AtomSet::addSpecies: species " << name
           << " is already defined" << endl;
    return false;
  }

  const double rcps = 1.5;
  sp->initialize(rcps);

  // create new entry in species list
  species_list.push_back(sp);
  spname.push_back(name);
  isp_[name] = species_list.size()-1;
  na_.insert(map<string,int>::value_type(name,0));
  atom_list.resize(atom_list.size()+1);
  is_[name] = spname.size()-1;

  if ( ctxt_.oncoutpe() )
  {
    cout << endl << " species " << sp->name() << ":" << endl;
    sp->info(cout);
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////
bool AtomSet::addSpecies(Species* sp, string name, const double rcpsin) {
  if ( findSpecies(name) || findMMSpecies(name) ) {
    if ( ctxt_.oncoutpe() )
      cout << " AtomSet::addSpecies: species " << name
           << " is already defined" << endl;
    return false;
  }
  
  sp->initialize(rcpsin);
  sp->fix_rcps = true;

  // create new entry in species list
  species_list.push_back(sp);
  spname.push_back(name);
  isp_[name] = species_list.size()-1;
  na_.insert(map<string,int>::value_type(name,0));
  atom_list.resize(atom_list.size()+1);
  
  return true;
}

////////////////////////////////////////////////////////////////////////////////
bool AtomSet::addMMSpecies(MMSpecies* sp, string name) {
  if ( findSpecies(name) || findMMSpecies(name) ) {
    if ( ctxt_.oncoutpe() )
      cout << " AtomSet::addMMSpecies: species " << name
           << " is already defined" << endl;
    return false;
  }

  // create new entry in species list
  mmspecies_list.push_back(sp);
  ispmm_[name] = mmspecies_list.size()-1;
  mmatom_list.resize(mmatom_list.size()+1);
  
  return true;
}


////////////////////////////////////////////////////////////////////////////////
bool AtomSet::addAtom(Atom *a)
{

  // check atom_list for name
  if ( findAtom(a->name()) )
  {
    // this name is already in the atom list, reject atom definition
    if ( ctxt_.oncoutpe() )
      cout << " AtomSet:addAtom: atom " << a->name()
           << " is already defined" << endl;
    return false;
  }

  // check if species is defined
  string spname = a->species();
  Species *s = findSpecies(spname);
  if ( !s )
  {
    // species not found, cannot define atom
    if ( ctxt_.oncoutpe() )
      cout << " AtomSet:addAtom: species " << spname
           << " is undefined" << endl;
    return false;
  }

  // add an atom to the atom_list
  int is = isp(spname);
  assert ( is >= 0 );
  atom_list[is].push_back(a);
  ia_[a->name()] = atom_list[is].size()-1;
  is_[a->name()] = is;

  // update count of atoms of species spname
  // increment count
  na_[spname]++;

  // update total number of electrons
  nel_ = 0;
  for ( int is = 0; is < species_list.size(); is++ )
  {
    nel_ += species_list[is]->zval() * atom_list[is].size();
  }

  return true;
}

////////////////////////////////////////////////////////////////////////////////
bool AtomSet::addMMAtom(Atom *a) {

  // check atom_list for name
  if ( findMMAtom(a->name()) ) {
    // this name is already in the atom list, reject atom definition
    if ( ctxt_.oncoutpe() )
      cout << " AtomSet:addMMAtom: atom " << a->name()
           << " is already defined" << endl;
    return false;
  }

  // check if species is defined
  string spname = a->species();
  MMSpecies *s = findMMSpecies(spname);
  if ( !s ) {
    // species not found, cannot define atom
    if ( ctxt_.oncoutpe() )
      cout << " AtomSet:addMMAtom: species " << spname
           << " is undefined" << endl;
    return false;
  }
  
  // add an atom to the atom_list
  int is = isp_mm(spname);
  assert ( is >= 0 );
  mmatom_list[is].push_back(a);

  return true;
}

////////////////////////////////////////////////////////////////////////////////
bool AtomSet::delAtom(string name)
{
  vector<vector<Atom*> >::iterator psa = atom_list.begin();
  vector<Species*>::iterator ps = species_list.begin();
  for ( int is = 0; is < species_list.size(); is++ )
  {
    vector<Atom*>::iterator pa = atom_list[is].begin();
    for ( int ia = 0; ia < atom_list[is].size(); ia++ )
    {
      if ( atom_list[is][ia]->name() == name )
      {
        string spname = atom_list[is][ia]->species();
        na_[spname]--;
        delete atom_list[is][ia];

        // remove map entries ia_[name] and is_[name]
        map<string,int>::iterator i = ia_.find(name);
        ia_.erase(i);
        i = is_.find(name);
        is_.erase(i);

        atom_list[is].erase(pa);
        nel_ -= species_list[is]->zval();

        return true;
      }
      pa++;
    }
    psa++;
    ps++;
  }

  // this name was not found in the atom list
  if ( ctxt_.oncoutpe() )
    cout << " AtomSet:delAtom: no such atom: " << name << endl;
  return false;
}

////////////////////////////////////////////////////////////////////////////////
bool AtomSet::addEmpiricalPotential(EmpiricalPotential* ep) {
  string name1 = ep->spname1();
  string name2 = ep->spname2();

  if ( ! (findSpecies(name1) || findMMSpecies(name1)) ) {
    if ( ctxt_.oncoutpe() )
      cout << " AtomSet::addEmpiricalPotential: species " << name1
           << " not defined" << endl;
    return false;
  }
  if ( ! (findSpecies(name2) || findMMSpecies(name2)) ) {
    if ( ctxt_.oncoutpe() )
      cout << " AtomSet::addEmpiricalPotential: species " << name2
           << " not defined" << endl;
    return false;
  }

  // determine is1 and is2, set in ep
  if (!findSpecies(name1))
    ep->is1 = isp_mm(name1) + nsp();
  else 
    ep->is1 = isp(name1);
  if (!findSpecies(name2))
    ep->is2 = isp_mm(name2) + nsp();
  else 
    ep->is2 = isp(name2);

  empirical_list.push_back(ep);

  //ewd DEBUG
  if ( ctxt_.oncoutpe() )
     cout << "<!-- AtomSet::addEmpiricalPotential  is1 = " << ep->is1 << ", is2 = " << ep->is2 << ":  " << empirical_list.size() << " potentials defined.  -->" << endl;


  
  return true;
}


////////////////////////////////////////////////////////////////////////////////
bool AtomSet::delMMAtom(string name) {
  vector<vector<Atom*> >::iterator psa = mmatom_list.begin();
  vector<MMSpecies*>::iterator ps = mmspecies_list.begin();
  for ( int is = 0; is < mmspecies_list.size(); is++ ) {
    vector<Atom*>::iterator pa = mmatom_list[is].begin();
    for ( int ia = 0; ia < mmatom_list[is].size(); ia++ ) {
      if ( mmatom_list[is][ia]->name() == name ) {
        string spname = mmatom_list[is][ia]->species();
        delete mmatom_list[is][ia];
        mmatom_list[is].erase(pa);
        return true;
      }
      pa++;
    }
    psa++;
    ps++;
  }
  
  // this name was not found in the atom list
  if ( ctxt_.oncoutpe() )
    cout << " AtomSet:delMMAtom: no such atom: " << name << endl;
  return false;
}

////////////////////////////////////////////////////////////////////////////////
Atom *AtomSet::findAtom(string name) const
{
  for ( int is = 0; is < species_list.size(); is++ )
  {
    for ( int ia = 0; ia < atom_list[is].size(); ia++ )
    {
      if ( atom_list[is][ia]->name() == name )
        return atom_list[is][ia];
    }
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////////
Atom *AtomSet::findMMAtom(string name) const {
  for ( int is = 0; is < mmspecies_list.size(); is++ ) {
    for ( int ia = 0; ia < mmatom_list[is].size(); ia++ ) {
      if ( mmatom_list[is][ia]->name() == name )
        return mmatom_list[is][ia];
    }
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////////
Species *AtomSet::findSpecies(string name) const
{
  int is = isp(name);
  if ( is >= 0 )
    return species_list[is];
  else
    return 0;
}

////////////////////////////////////////////////////////////////////////////////
MMSpecies *AtomSet::findMMSpecies(string name) const  {
  int is = isp_mm(name);
  if ( is >= 0 )
    return mmspecies_list[is];
  else
    return 0;
}

////////////////////////////////////////////////////////////////////////////////
double AtomSet::mass(int isp) {
  if (isp >= 0)
    return species_list[isp]->mass();
  else
    return 0;
}

////////////////////////////////////////////////////////////////////////////////
double AtomSet::mass_mm(int isp) {
  if (isp >= 0)
    return mmspecies_list[isp]->mass();
  else
    return 0;
}

////////////////////////////////////////////////////////////////////////////////
int AtomSet::atomic_number(int isp) const {
  if (isp >= 0)
    return species_list[isp]->atomic_number();
  else
    return 0;
}

////////////////////////////////////////////////////////////////////////////////
int AtomSet::isp(string& name) const
{
  map<string,int>::const_iterator i = isp_.find(name);
  if ( i != isp_.end() )
    return (*i).second;
  else
    return -1;
}

////////////////////////////////////////////////////////////////////////////////
int AtomSet::isp_mm(string name) const {
  map<string,int>::const_iterator i = ispmm_.find(name);
  if ( i != ispmm_.end() )
    return (*i).second;
  else
    return -1;
}

////////////////////////////////////////////////////////////////////////////////
int AtomSet::is(const string& atom_name) const
{
  map<string,int>::const_iterator i = is_.find(atom_name);
  if ( i != is_.end() )
    return (*i).second;
  else
    return -1;
}

////////////////////////////////////////////////////////////////////////////////
int AtomSet::ia(const string& atom_name) const
{
  map<string,int>::const_iterator i = ia_.find(atom_name);
  if ( i != ia_.end() )
    return (*i).second;
  else
    return -1;
}

////////////////////////////////////////////////////////////////////////////////
void AtomSet::listAtoms(void) const
{
  for ( int is = 0; is < species_list.size(); is++ )
  {
    for ( int ia = 0; ia < atom_list[is].size(); ia++ )
    {
      if ( ctxt_.oncoutpe() )
        cout << *atom_list[is][ia];
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void AtomSet::listSpecies(void) const
{
  for ( int is = 0; is < species_list.size(); is++ )
  {
    if ( ctxt_.oncoutpe() )
    {
      cout << endl << " species " << spname[is] << ":" << endl;
      species_list[is]->info(cout);
      cout << endl;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
int AtomSet::na(const string& spname) const
{
  map<string,int>::const_iterator i = na_.find(spname);
  if ( i != na_.end() )
    return (*i).second;
  else
    return 0;
}

////////////////////////////////////////////////////////////////////////////////
int AtomSet::na(int is) const
{
  assert( is >= 0 && is < atom_list.size() );
  return atom_list[is].size();
}

////////////////////////////////////////////////////////////////////////////////
int AtomSet::na_mm(int is) const {
  assert( is >= 0 && is < mmatom_list.size() );
  return mmatom_list[is].size();
}

////////////////////////////////////////////////////////////////////////////////
int AtomSet::size(void) const {
  int n = 0;
  for ( int is = 0; is < atom_list.size(); is++ )
    n += atom_list[is].size();
  for ( int is = 0; is < mmatom_list.size(); is++ )
    n += mmatom_list[is].size();
  return n;
}

////////////////////////////////////////////////////////////////////////////////
bool AtomSet::reset(void)
{
  // delete all atoms and species

  for ( int is = 0; is < species_list.size(); is++ )
  {
    for ( int ia = 0; ia < atom_list[is].size(); ia++ )
    {
      delete atom_list[is][ia];
    }
    atom_list[is].resize(0);

    delete species_list[is];
  }
  atom_list.resize(0);
  species_list.resize(0);

  for ( int is = 0; is < mmspecies_list.size(); is++ ) {
    for ( int ia = 0; ia < mmatom_list[is].size(); ia++ ) {
      delete mmatom_list[is][ia];
    }
    mmatom_list[is].resize(0);
    
    delete mmspecies_list[is];
  }
  mmatom_list.resize(0);
  mmspecies_list.resize(0);
  
  nel_ = 0;

  return true;
}

////////////////////////////////////////////////////////////////////////////////
void AtomSet::get_positions(vector<vector<double> >& tau, bool qmonly) const {
  if (tau.size() < atom_list.size()) 
    tau.resize(atom_list.size());

  for ( int is = 0; is < atom_list.size(); is++ ) {
    if (tau[is].size() != 3*atom_list[is].size())
      tau[is].resize(3*atom_list[is].size());
    int i = 0;
    for ( int ia = 0; ia < atom_list[is].size(); ia++ ) {
      D3vector t = atom_list[is][ia]->position();
      tau[is][i++] = t.x;
      tau[is][i++] = t.y;
      tau[is][i++] = t.z;
    }
  }

  // get positions of MM atoms, if any
  if (!qmonly) {
    if (tau.size() < atom_list.size()+mmatom_list.size())
      tau.resize(atom_list.size() + mmatom_list.size());

    const int offset = atom_list.size();
    for ( int is = 0; is < mmatom_list.size(); is++ ) {
      if (tau[is+offset].size() != 3*mmatom_list[is].size())
        tau[is+offset].resize(3*atom_list[is].size());
      int i = 0;
      for ( int ia = 0; ia < mmatom_list[is].size(); ia++ ) {
        D3vector t = mmatom_list[is][ia]->position();
        tau[is+offset][i++] = t.x;
        tau[is+offset][i++] = t.y;
        tau[is+offset][i++] = t.z;
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void AtomSet::set_positions(const vector<vector<double> >& tau, bool cellrescale) {

  assert(tau.size() == atom_list.size() || tau.size() == atom_list.size() + mmatom_list.size());
  for ( int is = 0; is < atom_list.size(); is++ ) {
    assert(tau[is].size() == 3*atom_list[is].size());
    int i = 0;
    for ( int ia = 0; ia < atom_list[is].size(); ia++ ) {
      atom_list[is][ia]->set_position(
         D3vector(tau[is][i],tau[is][i+1],tau[is][i+2]),cellrescale);
      i += 3;
    }
  }

  // set positions of MM atoms, if any
  if (tau.size() > atom_list.size()) {
    assert (tau.size() == atom_list.size() + mmatom_list.size());
    const int offset = atom_list.size();
    for ( int is = 0; is < mmatom_list.size(); is++ ) {
      assert(tau[is+offset].size() == 3*mmatom_list[is].size());
      int i = 0;
      for ( int ia = 0; ia < mmatom_list[is].size(); ia++ ) {
        mmatom_list[is][ia]->set_position(
           D3vector(tau[is+offset][i],tau[is+offset][i+1],tau[is+offset][i+2]),cellrescale);
        i += 3;
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void AtomSet::get_positions(vector<vector<double> >& tau) const
{
  // default value of qmonly = false
  get_positions(tau,false);
}

////////////////////////////////////////////////////////////////////////////////
void AtomSet::set_positions(const vector<vector<double> >& tau)
{
  // default value of cellrescale = false
  set_positions(tau,false);
}

////////////////////////////////////////////////////////////////////////////////
void AtomSet::get_velocities(vector<vector<double> >& vel, bool qmonly) const {
  if (vel.size() < atom_list.size())
    vel.resize(atom_list.size());
  for ( int is = 0; is < atom_list.size(); is++ ) {
    if (vel[is].size() != 3*atom_list[is].size())
      vel[is].resize(3*atom_list[is].size());
    int i = 0;
    for ( int ia = 0; ia < atom_list[is].size(); ia++ ) {
      D3vector t = atom_list[is][ia]->velocity();
      vel[is][i++] = t.x;
      vel[is][i++] = t.y;
      vel[is][i++] = t.z;
    }
  }

  // get velocities of MM atoms, if any
  if (!qmonly) {
    if (vel.size() < atom_list.size()+mmatom_list.size())
      vel.resize(atom_list.size() + mmatom_list.size());

    const int offset = atom_list.size();
    for ( int is = 0; is < mmatom_list.size(); is++ ) {
      if (vel[is+offset].size() != 3*mmatom_list[is].size())
        vel[is+offset].resize(3*mmatom_list[is].size());
      int i = 0;
      for ( int ia = 0; ia < mmatom_list[is].size(); ia++ ) {
        D3vector t = mmatom_list[is][ia]->velocity();
        vel[is+offset][i++] = t.x;
        vel[is+offset][i++] = t.y;
        vel[is+offset][i++] = t.z;
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void AtomSet::set_velocities(const vector<vector<double> >& vel) {
  assert(vel.size() == atom_list.size() || vel.size() == atom_list.size() + mmatom_list.size());
  for ( int is = 0; is < atom_list.size(); is++ ) {
    assert(vel[is].size() == 3*atom_list[is].size());
    int i = 0;
    for ( int ia = 0; ia < atom_list[is].size(); ia++ ) {
      atom_list[is][ia]->set_velocity(
        D3vector(vel[is][i],vel[is][i+1],vel[is][i+2]));
      i += 3;
    }
  }

  // set velocities of MM atoms, if any
  if (vel.size() > atom_list.size()) {
    assert (vel.size() == atom_list.size() + mmatom_list.size());
    const int offset = atom_list.size();
    for ( int is = 0; is < mmatom_list.size(); is++ ) {
      assert(vel[is+offset].size() == 3*mmatom_list[is].size());
      int i = 0;
      for ( int ia = 0; ia < mmatom_list[is].size(); ia++ ) {
        mmatom_list[is][ia]->set_velocity(
          D3vector(vel[is+offset][i],vel[is+offset][i+1],vel[is+offset][i+2]));
        i += 3;
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void AtomSet::get_velocities(vector<vector<double> >& vel) const
{
  // default value of qmonly = false
  get_velocities(vel,false);
}

////////////////////////////////////////////////////////////////////////////////
void AtomSet::get_fion_ext(vector<vector<double> >& fion_ext) const {
  if (fion_ext.size() < atom_list.size()) 
    fion_ext.resize(atom_list.size());
  for ( int is = 0; is < atom_list.size(); is++ ) {
    if (fion_ext[is].size() != 3*atom_list[is].size())
      fion_ext[is].resize(3*atom_list[is].size());
    for ( int ia = 0; ia < atom_list[is].size(); ia++ ) {
      fion_ext[is][3*ia] = fion_ext_[is][3*ia];
      fion_ext[is][3*ia+1] = fion_ext_[is][3*ia+1];
      fion_ext[is][3*ia+2] = fion_ext_[is][3*ia+2];
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
D3vector AtomSet::get_fion_ext(int is, int ia) {
  assert(fion_ext_.size() == atom_list.size());
  for ( int ist = 0; ist < atom_list.size(); ist++ )
    assert(fion_ext_[ist].size() == 3*atom_list[ist].size());
  D3vector ftmp(fion_ext_[is][3*ia],fion_ext_[is][3*ia+1],fion_ext_[is][3*ia+2]);
  return ftmp;
}

////////////////////////////////////////////////////////////////////////////////
void AtomSet::set_fion_ext(vector<vector<double> >& fion_ext) {

  assert(fion_ext.size() == atom_list.size());
  if (fion_ext_.size() != fion_ext.size())
    fion_ext_.resize(fion_ext.size());

  for ( int is = 0; is < fion_ext.size(); is++ ) {
    assert (fion_ext[is].size() == 3*atom_list[is].size());
    if (fion_ext_[is].size() != fion_ext[is].size())
      fion_ext_[is].resize(fion_ext[is].size());

    for ( int ia = 0; ia < fion_ext[is].size(); ia++ ) {
      fion_ext_[is][ia] = fion_ext[is][ia];
    }
  }

  add_fion_ext_ = true;
}

////////////////////////////////////////////////////////////////////////////////
void AtomSet::reset_velocities(void) {
  for ( int is = 0; is < atom_list.size(); is++ ) {
    for ( int ia = 0; ia < atom_list[is].size(); ia++ ) {
      atom_list[is][ia]->set_velocity(D3vector(0.0, 0.0, 0.0));
    }
  }

  for ( int is = 0; is < mmatom_list.size(); is++ ) {
    for ( int ia = 0; ia < mmatom_list[is].size(); ia++ ) {
      mmatom_list[is][ia]->set_velocity(D3vector(0.0, 0.0, 0.0));
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void AtomSet::randomize_velocities(double temp) {
   // distribute velocity components w. Maxwell-Boltzmann distribution, i.e.
   // each momentum component should have a normal distribution centered about
   // zero with a variance of m*kT
   const double boltz = 1.0 / ( 11605.0 * 2.0 * 13.6058 ); // convert from Kelvin to Ha
   D3vector vrand;
   for ( int is = 0; is < atom_list.size(); is++ ) {
      double mass = species_list[is]->mass()*1822.89;
      const double norm = sqrt(mass*boltz*temp);
      for ( int ia = 0; ia < atom_list[is].size(); ia++ ) {
         for (int i=0; i<3; i++)
         {
            double u,v;
            double r = 0.0;
            int rcnt = 0;
            while (rcnt < 100 && (r == 0.0 || r > 1.0)) { // force r within 0 < r <= 1
               u = 2.*drand48()-1.;
               v = 2.*drand48()-1.;
               r = u*u + v*v;
               rcnt++;
            }
            if (rcnt >= 100)
               cout << "<ERROR> AtomSet::randomize_velocities failed on task " << ctxt_.mype() << " </ERROR>" << endl;

            double c = sqrt(-2.*log(r)/r);
            vrand[i] = c*u*norm/mass;
         }
         atom_list[is][ia]->set_velocity(vrand);
      }
   }

   for ( int is = 0; is < mmatom_list.size(); is++ ) {
      double mass = mmspecies_list[is]->mass() * 1822.89;
      const double norm = sqrt(mass*boltz*temp);
      for ( int ia = 0; ia < mmatom_list[is].size(); ia++ ) {
         for (int i=0; i<3; i++)
         {
            double u,v;
            double r = 0.0;
            while (r != 0.0 && r <= 1.0) { // force r within 0 < r <= 1
               u = 2.*drand48()-1.;
               v = 2.*drand48()-1.;
               r = u*u + v*v;
            }
            double c = sqrt(-2.*log(r)/r);
            vrand[i] = c*u*norm/mass;
         }
         mmatom_list[is][ia]->set_velocity(vrand);
      }
   }

   // compute kinetic energy, rescale to achieve target temperature
   double ekin_ = 0.0;
   for ( int is = 0; is < atom_list.size(); is++ ) {
      double mass = species_list[is]->mass()*1822.89;
      for ( int ia = 0; ia < atom_list[is].size(); ia++ ) {
         D3vector vv = atom_list[is][ia]->velocity();
         ekin_ += 0.5 * mass * (vv.x*vv.x + vv.y*vv.y + vv.z*vv.z);
      }
   }
   for ( int is = 0; is < mmatom_list.size(); is++ ) {
      double mass = mmspecies_list[is]->mass()*1822.89;
      for ( int ia = 0; ia < mmatom_list[is].size(); ia++ ) {
         D3vector vv = mmatom_list[is][ia]->velocity();
         ekin_ += 0.5 * mass * (vv.x*vv.x + vv.y*vv.y + vv.z*vv.z);
      }
   }

   // compute ndofs = 3*natot ignoring constraints for now
   int natot = 0;
   for ( int is = 0; is < atom_list.size(); is++ )
      natot += atom_list[is].size();

   double scale = sqrt(boltz*temp*3.*natot/(2.*ekin_));
   for ( int is = 0; is < atom_list.size(); is++ ) {
      double mass = species_list[is]->mass()*1822.89;
      for ( int ia = 0; ia < atom_list[is].size(); ia++ ) {
         D3vector vv = atom_list[is][ia]->velocity();
         vv *= scale;
         atom_list[is][ia]->set_velocity(vv);         
      }
   }
   for ( int is = 0; is < mmatom_list.size(); is++ ) {
      double mass = mmspecies_list[is]->mass()*1822.89;
      for ( int ia = 0; ia < mmatom_list[is].size(); ia++ ) {
         D3vector vv = mmatom_list[is][ia]->velocity();
         vv *= scale;
         mmatom_list[is][ia]->set_velocity(vv);         
      }
   }

   //ewd DEBUG:  check kinetic energy from rescaled velocities
   ekin_ = 0.0;
   for ( int is = 0; is < atom_list.size(); is++ ) {
      double mass = species_list[is]->mass()*1822.89;
      for ( int ia = 0; ia < atom_list[is].size(); ia++ ) {
         D3vector vv = atom_list[is][ia]->velocity();
         ekin_ += 0.5 * mass * (vv.x*vv.x + vv.y*vv.y + vv.z*vv.z);
      }
   }
   for ( int is = 0; is < mmatom_list.size(); is++ ) {
      double mass = mmspecies_list[is]->mass()*1822.89;
      for ( int ia = 0; ia < mmatom_list[is].size(); ia++ ) {
         D3vector vv = mmatom_list[is][ia]->velocity();
         ekin_ += 0.5 * mass * (vv.x*vv.x + vv.y*vv.y + vv.z*vv.z);
      }
   }
   if (ctxt_.oncoutpe())
      cout << "<!-- AtomSet::randomize_velocities:  average kinetic energy = " << ekin_/boltz/(1.5*natot) << " K, target temp = " << temp << " K -->" << endl;
   //ewd DEBUG
   



}

////////////////////////////////////////////////////////////////////////////////
D3vector AtomSet::vcm(void) const
{
  D3vector mvsum;
  double msum = 0.0;
  for ( int is = 0; is < atom_list.size(); is++ )
  {
    double mass = species_list[is]->mass();
    for ( int ia = 0; ia < atom_list[is].size(); ia++ )
    {
      D3vector v = atom_list[is][ia]->velocity();
      mvsum += mass * v;
      msum += mass;
    }
  }

  for ( int is = 0; is < mmatom_list.size(); is++ )
  {
    double mass = mmspecies_list[is]->mass();
    for ( int ia = 0; ia < mmatom_list[is].size(); ia++ )
    {
      D3vector v = mmatom_list[is][ia]->velocity();
      mvsum += mass * v;
      msum += mass;
    }
  }

  if ( msum == 0.0 ) return D3vector(0,0,0);
  return mvsum / msum;
}

////////////////////////////////////////////////////////////////////////////////
void AtomSet::reset_vcm(void)
{
  D3vector vc = vcm();
  vector<vector<double> > v;
  get_velocities(v,false);
  // subtract center of mass velocity
  for ( int is = 0; is < atom_list.size(); is++ )
  {
    int i = 0;
    for ( int ia = 0; ia < atom_list[is].size(); ia++ )
    {
      v[is][i++] -= vc.x;
      v[is][i++] -= vc.y;
      v[is][i++] -= vc.z;
    }
  }
  const int offset = atom_list.size();
  for ( int is = 0; is < mmatom_list.size(); is++ )
  {
    int i = 0;
    for ( int ia = 0; ia < mmatom_list[is].size(); ia++ )
    {
      v[is+offset][i++] -= vc.x;
      v[is+offset][i++] -= vc.y;
      v[is+offset][i++] -= vc.z;
    }
  }
  set_velocities(v);
}

////////////////////////////////////////////////////////////////////////////////
D3vector AtomSet::dipole(void) const
{
  D3vector sum;
  for ( int is = 0; is < atom_list.size(); is++ )
  {
    double charge = species_list[is]->zval();
    for ( int ia = 0; ia < atom_list[is].size(); ia++ )
    {
      D3vector p = atom_list[is][ia]->position();
      sum += charge * p;
    }
  }
  // include MM atoms?  For now, no.
  /*
  for ( int is = 0; is < mmatom_list.size(); is++ )
  {
    double charge = mmspecies_list[is]->zval();
    for ( int ia = 0; ia < mmatom_list[is].size(); ia++ )
    {
      D3vector p = mmatom_list[is][ia]->position();
      sum += charge * p;
    }
  }
  */
  return sum;
}

////////////////////////////////////////////////////////////////////////////////
void AtomSet::sync()
{
#if USE_MPI
  // enforce consistency of positions and velocities on all tasks
  // broadcast positions and velocities of task 0 to all tasks
  vector<vector<double> > r,v;
  get_positions(r);
  get_velocities(v);
  for ( int is = 0; is < atom_list.size(); is++ )
  {
    int m = r[is].size();
    double* p = &r[is][0];
    if ( ctxt_.oncoutpe() )
    {
      ctxt_.dbcast_send(m,1,p,m);
    }
    else
    {
      ctxt_.dbcast_recv(m,1,p,m,0,0);
    }
  }
  for ( int is = 0; is < atom_list.size(); is++ )
  {
    int m = v[is].size();
    double* p = &v[is][0];
    if ( ctxt_.oncoutpe() )
    {
      ctxt_.dbcast_send(m,1,p,m);
    }
    else
    {
      ctxt_.dbcast_recv(m,1,p,m,0,0);
    }
  }

  // repeat for MM atoms, if any
  const int offset = atom_list.size();
  for ( int is = 0; is < mmatom_list.size(); is++ )
  {
    int m = r[is+offset].size();
    double* p = &r[is+offset][0];
    if ( ctxt_.oncoutpe() )
    {
      ctxt_.dbcast_send(m,1,p,m);
    }
    else
    {
      ctxt_.dbcast_recv(m,1,p,m,0,0);
    }
  }
  for ( int is = 0; is < mmatom_list.size(); is++ )
  {
    int m = v[is+offset].size();
    double* p = &v[is+offset][0];
    if ( ctxt_.oncoutpe() )
    {
      ctxt_.dbcast_send(m,1,p,m);
    }
    else
    {
      ctxt_.dbcast_recv(m,1,p,m,0,0);
    }
  }
  set_positions(r);
  set_velocities(v);
#endif
}
////////////////////////////////////////////////////////////////////////////////
void AtomSet::fold_in_ws(void)
{
  vector<vector<double> > p;
  get_positions(p);
  for ( int is = 0; is < atom_list.size(); is++ )
  {
    for ( int ia = 0; ia < atom_list[is].size(); ia++ )
    {
      D3vector pos(&p[is][3*ia]);
      cell_.fold_in_ws(pos);
      p[is][3*ia+0] = pos.x;
      p[is][3*ia+1] = pos.y;
      p[is][3*ia+2] = pos.z;
    }
  }

  // MM atoms, if any
  const int offset = atom_list.size();
  for ( int is = 0; is < mmatom_list.size(); is++ )
  {
    for ( int ia = 0; ia < mmatom_list[is].size(); ia++ )
    {
      D3vector pos(&p[is+offset][3*ia]);
      cell_.fold_in_ws(pos);
      p[is+offset][3*ia+0] = pos.x;
      p[is+offset][3*ia+1] = pos.y;
      p[is+offset][3*ia+2] = pos.z;
    }
  }
  set_positions(p);
}

////////////////////////////////////////////////////////////////////////////////
void AtomSet::findSymmetricAtoms(const SymmetrySet& symset) {

  const double threshold = 1.E-6;  // distance threshold used to identify equivalent atoms

  int nsym_ = symset.nsym();
  if (nsym_ <= 0)
    return;

  // resize symatomid_ to current AtomSet, SymmetrySet sizes
  if (symatomid_.size() != atom_list.size())
    symatomid_.resize(atom_list.size());
  for ( int is = 0; is < atom_list.size(); is++ )
    if (symatomid_[is].size() != atom_list[is].size())
      symatomid_[is].resize(atom_list[is].size());
  for ( int is = 0; is < atom_list.size(); is++ )
    for ( int ia = 0; ia < atom_list[is].size(); ia++ )
      if (symatomid_[is][ia].size() != nsym_)
        symatomid_[is][ia].resize(nsym_);

  // initialize array
  for ( int is = 0; is < atom_list.size(); is++ )
    for ( int ia = 0; ia < atom_list[is].size(); ia++ )
      for (int sy = 0; sy < nsym_; sy++)
        symatomid_[is][ia][sy] = -1;

  // fill array by species
  for ( int is = 0; is < atom_list.size(); is++ ) {

    // loop through atoms and calculate their positions in crystal coordinates
    vector<D3vector> r_xtal;
    r_xtal.resize(atom_list[is].size());
    for ( int ia = 0; ia < atom_list[is].size(); ia++ )
      r_xtal[ia] = cell_.cart_to_crystal(atom_list[is][ia]->position());

    // for each atom, calculate the position resulting from each symmetry operation
    // and find the id corresponding to this position

    //ewd:  we're currently just doing the simple N^2 approach -- should be cheap compared
    // with wavefunction costs, but check the timing...

    for ( int ia = 0; ia < atom_list[is].size(); ia++ ) {
      for (int sy = 0; sy < nsym_; sy++) {
        D3vector r_sym = symset.symlist[sy]->applyToVector(r_xtal[ia],true);

        // compare this position against other atoms
        for ( int ja = 0; ja < atom_list[is].size(); ja++ ) {
          // in crystal coordinates, agreement can be off by an integer
          double xdiff = (r_sym.x - r_xtal[ja].x);
          int xdiffnint = (int) ( xdiff > 0.0 ? xdiff+0.5 : xdiff-0.5);
          double ydiff = (r_sym.y - r_xtal[ja].y);
          int ydiffnint = (int) ( ydiff > 0.0 ? ydiff+0.5 : ydiff-0.5);
          double zdiff = (r_sym.z - r_xtal[ja].z);
          int zdiffnint = (int) ( zdiff > 0.0 ? zdiff+0.5 : zdiff-0.5);

          if (abs(xdiff - (double)xdiffnint) < threshold && abs(ydiff - (double)ydiffnint) < threshold && abs(zdiff - (double)zdiffnint) < threshold) {
            if (symatomid_[is][ia][sy] != -1)
              if ( ctxt_.oncoutpe() )
                cout << "<!-- ERROR:  AtomSet.findSymmetricAtoms found multiple symmetry-equivalent atoms!  " << is << ", " << ia << ", " << ja << ", " << sy << "-->" << endl;
            symatomid_[is][ia][sy] = ja;
          }
        }
        if (symatomid_[is][ia][sy] == -1)
          if ( ctxt_.oncoutpe() )
            cout << "<!-- ERROR:  AtomSet.findSymmetricAtoms could not find symmetry-equivalent atom!  " << is << ", " << ia << ", " << sy << "-->" << endl;
      }
    }
  }

}

////////////////////////////////////////////////////////////////////////////////
void AtomSet::set_rcps(const double& ecut) const {

  // choose width of Gaussian charge so that both real-space and reciprocal space
  // terms in Ewald sum converge, and all charge is enclosed in the box

  double ellmin = cell_.min_wsdist();
  const double thresh_recip = 1.E-7;  // fraction of F.T. of charge Gaussian to sqrt(ecut)

  for ( int is = 0; is < species_list.size(); is++ ) {
    if (! species_list[is]->fix_rcps) {

      double rcps = 0.3;  // starting default value
      double rcps_orig = species_list[is]->rcps();
      const double ngrid = 1.0;  // minimum number of 2*pi/L points for 1/rcps to cover in recip. space
      const double zval = (double) species_list[is]->zval();
      const double recipmax = ellmin/(ngrid*2.*M_PI);

      bool recip_cut = false;
      while (!recip_cut) {
        double recip_area = (zval/sqrt(M_PI))*erfc(rcps*sqrt(ecut));
        if (recip_area > thresh_recip*rcps/zval) {
          rcps += 0.05;
        }
        else {
          recip_cut = true;
        }
      }
      if (rcps != rcps_orig)
        if ( ctxt_.oncoutpe() )
          cout << "<!-- AtomSet.set_rcps:  Ewald width for species " << species_list[is]->name() << " is too small for reciprocal sum convergence, increasing rcps from " << rcps_orig << " to " << rcps << " -->" << endl;
      
      
      if (rcps > recipmax) {
        if ( ctxt_.oncoutpe() )
          cout << "<!-- AtomSet.set_rcps:  species " << species_list[is]->name() << " Ewald width larger than reciprocal grid spacing, decreasing rcps from " << rcps << " to " << recipmax << " -->" << endl;
        rcps = recipmax;
      }

      // check that all charge is enclosed 
      const double thresh_charge = 1.E-8;

      bool charge_encl = false;
      double rcps_prev = rcps;
      while (!charge_encl) {
        double charge_spill = zval*erfc(ellmin/(sqrt(2.)*rcps));
        if (charge_spill > thresh_charge) {
          rcps -= 0.05;
        }
        else {
          charge_encl = true;
        }
      }

      if (rcps_prev != rcps) 
        if ( ctxt_.oncoutpe() )
          cout << "<!-- AtomSet.set_rcps:  Ewald width for species " << species_list[is]->name() << " is too large for box size, decreasing radius from " << rcps_prev << " to " << rcps << " so charge outside box is less than " << setprecision(8) << thresh_charge << " -->" << endl;
      
      // update species with new radius
      if (rcps != rcps_orig) 
        species_list[is]->initialize(rcps);
      
    }
    if ( ctxt_.oncoutpe() )
      cout << "<!-- AtomSet.set_rcps:  Ewald width for species " << species_list[is]->name() << " = " << species_list[is]->rcps() << " -->" << endl;
    
  }
  
}
////////////////////////////////////////////////////////////////////////////////
void AtomSet::printsys(ostream& os) const {

  for (int k=0; k < empirical_list.size(); k++) 
    empirical_list[k]->printsys(os);

  string atomcmd = "atom";
  for ( int is = 0; is < species_list.size(); is++ ) {
    for ( int ia = 0; ia < atom_list[is].size(); ia++ ) {
      atom_list[is][ia]->printsys(os,atomcmd);
    }
  }

  atomcmd = "mmatom";
  for ( int is = 0; is < mmspecies_list.size(); is++ ) {
    for ( int ia = 0; ia < mmatom_list[is].size(); ia++ ) {
      mmatom_list[is][ia]->printsys(os,atomcmd);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void AtomSet::print_casino(ostream& os) const {

  os << "GEOMETRY" << endl;
  os << "--------" << endl;
  os << "Number of atoms per primitive cell" << endl;

  int natoms = 0;
  for ( int is = 0; is < species_list.size(); is++ ) 
    natoms += na(is);

  os << natoms << endl;
  os << "Atomic numbers and positions of atoms (au)" << endl;

  for ( int is = 0; is < species_list.size(); is++ ) {
    int atnum = species_list[is]->atomic_number();
    for ( int ia = 0; ia < atom_list[is].size(); ia++ )
      os << " " << atnum << "  " << atom_list[is][ia]->position() << endl;
  }
  os << "Primitive lattice vectors (au)" << endl;
  for (int i=0; i<3; i++)
    os << cell_.a(i) << endl;
  os << endl;
}

////////////////////////////////////////////////////////////////////////////////
ostream& operator << ( ostream &os, const AtomSet &as )
{
  if ( as.context().oncoutpe() )
  {
    os << "<atomset>\n";
    os << as.cell();
    for ( int is = 0; is < as.species_list.size(); is++ )
    {
      os << *as.species_list[is];
    }
    for ( int is = 0; is < as.species_list.size(); is++ )
    {
      for ( int ia = 0; ia < as.atom_list[is].size(); ia++ )
      {
        os << *as.atom_list[is][ia];
      }
    }
  }
  os << "</atomset>\n";

  return os;
}
