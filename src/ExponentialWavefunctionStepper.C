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
// ExponentialWavefunctionStepper.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ExponentialWavefunctionStepper.C,v 1.8 2011-06-02 15:56:19 schleife Exp $

#include "ExponentialWavefunctionStepper.h"
#include "SelfConsistentPotential.h"
#include "SlaterDet.h"
#include "Sample.h"
#include <iostream>
#include <deque>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
ExponentialWavefunctionStepper::ExponentialWavefunctionStepper(Wavefunction& wf, double tddt, TimerMap& tmap, EnergyFunctional & ef, Sample & s, bool approximated)
  : tddt_(tddt), WavefunctionStepper(wf,tmap), ef_(ef), s_(s), approximated_(approximated)
{
  order_ = 4;
  potential_.resize(3);
  stored_iter_ = 0;
}

////////////////////////////////////////////////////////////////////////////////
void ExponentialWavefunctionStepper::exponential(const double & dt, Wavefunction * dwf){

  // dummy variables to call ef_.energy
  std::vector<std::vector<double> > fion;
  std::valarray<double> sigma;
  bool delete_dwf;
  Wavefunction expwf(wf_);
  
  if (dwf == 0)
  {
    delete_dwf = true;
    dwf = new Wavefunction(wf_);
    ef_.energy(true, *dwf, false, fion, false, sigma);
  }
  else
  {
    delete_dwf = false;
  }

  // exp(A)x ~= x + Ax + A^2x/2! + A^3x/3! + A^4x/4!
  //          = x + Ax + A(Ax)/2 + A(A^2x/2)/3 + A(A^3x/6)/4 
  
  // order 0
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++)
    for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
      expwf.sd(ispin, ikp)->c() = wf_.sd(ispin, ikp)->c();

  complex<double> factor = 1.0;
 
  //order N:
  for(int N = 1; N <= order_; N++)
  {
    factor *= -complex<double>(0.0, 1.0)*dt/double(N);

    if(N != 1)        // for N == 1 this is already done
    {
       //move dwf to wf
       for ( int ispin = 0; ispin < wf_.nspin(); ispin++)
          for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
             wf_.sd(ispin, ikp)->c() = dwf->sd(ispin, ikp)->c();
      
       // apply A
       ef_.energy(true, *dwf, false, fion, false, sigma);
    }
    
    //accumulate the result
    for ( int ispin = 0; ispin < wf_.nspin(); ispin++)
       for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
          expwf.sd(ispin, ikp)->c().axpy(factor, dwf->sd(ispin, ikp)->c());
    
  }
  
  // copy the result back THE wavefunction
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++)
  {
    for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
    {
      wf_.sd(ispin, ikp)->c() = expwf.sd(ispin, ikp)->c();
      s_.hamil_wf->sd(ispin, ikp)->c() = wf_.sd(ispin, ikp)->c();
    }
  }
  
  if(delete_dwf) delete dwf;
  
}


////////////////////////////////////////////////////////////////////////////////
void ExponentialWavefunctionStepper::preupdate()
{
  if (approximated_)
  {
    // save the potential
    potential_[0] = potential_[1];
    potential_[1] = potential_[2];
    potential_[2] = ef_.get_self_consistent_potential();
    stored_iter_++;
  }

  // The propagator is U(t + dt, t) = exp(-i dt/2 H(t + dt)) exp(-i dt/2 H(t))

  // propagate to the half of the step with H(t)
  exponential(0.5*tddt_);

  if( approximated_ && stored_iter_ >= 3 )
  {
     SelfConsistentPotential future_potential;
     future_potential.extrapolate(potential_);
     ef_.set_self_consistent_potential(future_potential);

     // now do the other half of the propagation with H(t + dt)
     exponential(0.5*tddt_);
  }

}

////////////////////////////////////////////////////////////////////////////////
void ExponentialWavefunctionStepper::update(Wavefunction& dwf)
{

   // Now we need to get the Hamiltonian at time t + dt
   //
   // Here we consider that the atoms (and their potential) was
   // propagated to t + dt by the calling routine

   if ( !approximated_  || stored_iter_ < 3 )
   {
      Wavefunction wf_half(wf_);
    
      //save the current status
      for ( int ispin = 0; ispin < wf_.nspin(); ispin++)
         for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
            wf_half.sd(ispin, ikp)->c() = wf_.sd(ispin, ikp)->c();
    
      // do the rest of the propagation step with H(t)
      exponential(0.5*tddt_, &dwf);

      // get the approximated vhxc for the end of the step
      ef_.hamil_cd()->update_density();
      ef_.update_hamiltonian();
      ef_.update_vhxc();

      // restore the wf to the middle of the step
      for ( int ispin = 0; ispin < wf_.nspin(); ispin++)
         for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
            wf_.sd(ispin, ikp)->c() = wf_half.sd(ispin, ikp)->c();

      // now do the other half of the propagation with H(t + dt)
      exponential(0.5*tddt_);
    
   }
   
}
