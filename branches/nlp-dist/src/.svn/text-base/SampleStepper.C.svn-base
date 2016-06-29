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
// SampleStepper.C
//
////////////////////////////////////////////////////////////////////////////////

#include "SampleStepper.h"
#include "Species.h"

#include <iostream>
#include <iomanip>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
SampleStepper::SampleStepper(Sample& s) : s_(s), atoms_(s_.atoms)
{
  fion.resize(s_.atoms.nsp()+s_.atoms.nsp_mm());
  for ( int is = 0; is < fion.size(); is++ )
    fion[is].resize(3*s_.atoms.na(is));
    
  sigma_eks.resize(6);
  sigma_kin.resize(6);
  sigma_ext.resize(6);
  sigma.resize(6);
  sigma_ext = valarray<double>(s.ctrl.ext_stress,6);
  // ewd:  add more significant figures to conversion
  const double gpa = 29421.0120;
  //const double gpa = 29421.5;
  sigma_ext /= gpa;
}

////////////////////////////////////////////////////////////////////////////////
SampleStepper::~SampleStepper(void)
{
  // print timer map
  for ( TimerMap::iterator i = tmap.begin(); i != tmap.end(); i++ )
  {
    double time = (*i).second.real();
    double tmin = time;
    double tmax = time;
    s_.ctxt_.dmin(1,1,&tmin,1);
    s_.ctxt_.dmax(1,1,&tmax,1);
    if ( s_.ctxt_.myproc()==0 )
    {
      cout << "<!-- timing "
           << setw(15) << (*i).first
           << " : " << setprecision(3) << setw(9) << tmin
           << " "   << setprecision(3) << setw(9) << tmax << " -->" << endl;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void SampleStepper::print_stress(void)
{
  // ewd:  add more significant figures to conversion
  const double gpa = 29421.0120;
  //const double gpa = 29421.5;
  if ( s_.ctxt_.oncoutpe() )
  {
    cout.setf(ios::fixed,ios::floatfield);
    cout.setf(ios::right,ios::adjustfield);
    cout << setprecision(8);
    cout << " <stress_tensor unit=\"GPa\">" << endl;

    cout << "   <sigma_eks_xx> " << setw(12) << sigma_eks[0]*gpa
         << " </sigma_eks_xx>\n"
         << "   <sigma_eks_yy> " << setw(12) << sigma_eks[1]*gpa
         << " </sigma_eks_yy>\n"
         << "   <sigma_eks_zz> " << setw(12) << sigma_eks[2]*gpa
         << " </sigma_eks_zz>\n"
         << "   <sigma_eks_xy> " << setw(12) << sigma_eks[3]*gpa
         << " </sigma_eks_xy>\n"
         << "   <sigma_eks_yz> " << setw(12) << sigma_eks[4]*gpa
         << " </sigma_eks_yz>\n"
         << "   <sigma_eks_xz> " << setw(12) << sigma_eks[5]*gpa
         << " </sigma_eks_xz>\n\n";

    cout << "   <sigma_kin_xx> " << setw(12) << sigma_kin[0]*gpa
         << " </sigma_kin_xx>\n"
         << "   <sigma_kin_yy> " << setw(12) << sigma_kin[1]*gpa
         << " </sigma_kin_yy>\n"
         << "   <sigma_kin_zz> " << setw(12) << sigma_kin[2]*gpa
         << " </sigma_kin_zz>\n"
         << "   <sigma_kin_xy> " << setw(12) << sigma_kin[3]*gpa
         << " </sigma_kin_xy>\n"
         << "   <sigma_kin_yz> " << setw(12) << sigma_kin[4]*gpa
         << " </sigma_kin_yz>\n"
         << "   <sigma_kin_xz> " << setw(12) << sigma_kin[5]*gpa
         << " </sigma_kin_xz>\n\n";

    cout << "   <sigma_ext_xx> " << setw(12) << sigma_ext[0]*gpa
         << " </sigma_ext_xx>\n"
         << "   <sigma_ext_yy> " << setw(12) << sigma_ext[1]*gpa
         << " </sigma_ext_yy>\n"
         << "   <sigma_ext_zz> " << setw(12) << sigma_ext[2]*gpa
         << " </sigma_ext_zz>\n"
         << "   <sigma_ext_xy> " << setw(12) << sigma_ext[3]*gpa
         << " </sigma_ext_xy>\n"
         << "   <sigma_ext_yz> " << setw(12) << sigma_ext[4]*gpa
         << " </sigma_ext_yz>\n"
         << "   <sigma_ext_xz> " << setw(12) << sigma_ext[5]*gpa
         << " </sigma_ext_xz>\n\n";

    cout << "   <sigma_xx> " << setw(12) << sigma[0]*gpa
         << " </sigma_xx>\n"
         << "   <sigma_yy> " << setw(12) << sigma[1]*gpa
         << " </sigma_yy>\n"
         << "   <sigma_zz> " << setw(12) << sigma[2]*gpa
         << " </sigma_zz>\n"
         << "   <sigma_xy> " << setw(12) << sigma[3]*gpa
         << " </sigma_xy>\n"
         << "   <sigma_yz> " << setw(12) << sigma[4]*gpa
         << " </sigma_yz>\n"
         << "   <sigma_xz> " << setw(12) << sigma[5]*gpa
         << " </sigma_xz>\n";

    cout << " </stress_tensor>" << endl;
  }
}

////////////////////////////////////////////////////////////////////////////////
void SampleStepper::compute_sigma(void)
{
  sigma_kin = 0.0;
  // compute kinetic contribution to stress using velocities at time t0
  for ( int is = 0; is < atoms_.atom_list.size(); is++ )
  {
    int i = 0;
    double mass = atoms_.species_list[is]->mass() * 1822.89;
    for ( int ia = 0; ia < atoms_.atom_list[is].size(); ia++ )
    {
      Atom* pa = atoms_.atom_list[is][ia];
      D3vector v = pa->velocity();
      const double vx = v.x;
      const double vy = v.y;
      const double vz = v.z;
 
      sigma_kin[0] += mass * vx * vx;
      sigma_kin[1] += mass * vy * vy;
      sigma_kin[2] += mass * vz * vz;
      sigma_kin[3] += mass * vx * vy;
      sigma_kin[4] += mass * vy * vz;
      sigma_kin[5] += mass * vx * vz;
 
      i += 3;
    }
  }

  // if MM atoms are present, add to kinetic contribution to stress using velocities at time t0
  for ( int is = 0; is < atoms_.mmatom_list.size(); is++ ) {
    int i = 0;
    double mass = atoms_.mmspecies_list[is]->mass() * 1822.89;
    for ( int ia = 0; ia < atoms_.mmatom_list[is].size(); ia++ ) {
      Atom* pa = atoms_.mmatom_list[is][ia];
      D3vector v = pa->velocity();
      const double vx = v.x;
      const double vy = v.y;
      const double vz = v.z;
 
      sigma_kin[0] += mass * vx * vx;
      sigma_kin[1] += mass * vy * vy;
      sigma_kin[2] += mass * vz * vz;
      sigma_kin[3] += mass * vx * vy;
      sigma_kin[4] += mass * vy * vz;
      sigma_kin[5] += mass * vx * vz;
 
      i += 3;
    }
  }

  sigma_kin /= s_.wf.cell().volume();
 
  sigma = sigma_eks + sigma_kin - sigma_ext;
}


