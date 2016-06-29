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
    : tddt_(tddt), WavefunctionStepper(wf,tmap), ef_(ef), s_(s), approximated_(approximated), expwf_(s.wf), wfhalf_(s.wf), newwf_(s.wf) 
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
  
  if (dwf == 0)
  {
    delete_dwf = true;
    //ewd:  don't want to keep calling Wavefunction constructor, should be able to use wfhalf_ for this call
    //dwf = new Wavefunction(wf_);
    dwf = &wfhalf_;
    tmap_["expowf_ef"].start();
    ef_.energy(true, *dwf, false, fion, false, sigma);
    tmap_["expowf_ef"].stop();
  }
  else
  {
    delete_dwf = false;
  }

  // exp(A)x ~= x + Ax + A^2x/2! + A^3x/3! + A^4x/4!
  //          = x + Ax + A(Ax)/2 + A(A^2x/2)/3 + A(A^3x/6)/4 
  
  // order 0
  tmap_["expowf_copy"].start();
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++)
    for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
      expwf_.sd(ispin, ikp)->c() = wf_.sd(ispin, ikp)->c();
  tmap_["expowf_copy"].stop();

  complex<double> factor = 1.0;
 
  //order N:
  for(int N = 1; N <= order_; N++)
  {
    factor *= -complex<double>(0.0, 1.0)*dt/double(N);

    if(N != 1)        // for N == 1 this is already done
    {
       tmap_["expowf_copy"].start();
       //move dwf to wf
       for ( int ispin = 0; ispin < wf_.nspin(); ispin++)
          for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
             wf_.sd(ispin, ikp)->c() = dwf->sd(ispin, ikp)->c();
       tmap_["expowf_copy"].stop();
      
       // apply A
       tmap_["expowf_ef"].start();
       ef_.energy(true, *dwf, false, fion, false, sigma);
       tmap_["expowf_ef"].stop();
    }
    
    tmap_["expowf_axpy"].start();
    //accumulate the result
    for ( int ispin = 0; ispin < wf_.nspin(); ispin++)
       for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
          expwf_.sd(ispin, ikp)->c().axpy(factor, dwf->sd(ispin, ikp)->c());
    tmap_["expowf_axpy"].stop();
    
  }
  
  tmap_["expowf_copy"].start();
  // copy the result back THE wavefunction
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++)
  {
    for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
    {
      wf_.sd(ispin, ikp)->c() = expwf_.sd(ispin, ikp)->c();
      s_.hamil_wf->sd(ispin, ikp)->c() = wf_.sd(ispin, ikp)->c();
    }
  }
  tmap_["expowf_copy"].stop();
  
  //if(delete_dwf) delete dwf;
  
}


////////////////////////////////////////////////////////////////////////////////
void ExponentialWavefunctionStepper::exponential_condensed(const double & dt1, const double & dt2, Wavefunction * dwf){

  // dummy variables to call ef_.energy
  std::vector<std::vector<double> > fion;
  std::valarray<double> sigma;
  bool delete_dwf;
  
  if (dwf == 0)
  {
    delete_dwf = true;
    //ewd:  don't want to keep calling Wavefunction constructor, should be able to use wfhalf_ for this call
    //dwf = new Wavefunction(wf_);
    dwf = &wfhalf_;
    tmap_["expowf_ef"].start();
    ef_.energy(true, *dwf, false, fion, false, sigma);
    tmap_["expowf_ef"].stop();
  }
  else
  {
    delete_dwf = false;
  }

  // exp(A)x ~= x + Ax + A^2x/2! + A^3x/3! + A^4x/4!
  //          = x + Ax + A(Ax)/2 + A(A^2x/2)/3 + A(A^3x/6)/4 
  
  // order 0
  // For the 0th order term of the Taylor expansion, simply copy the
  // wavefunctions from wf_ to expwf_ for the first exponential calculation
  // and from wf_ to newwf_ for the second exponential calculation
  tmap_["expowf_copy"].start();
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++)
    for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
    {
      expwf_.sd(ispin, ikp)->c() = wf_.sd(ispin, ikp)->c();
      newwf_.sd(ispin, ikp)->c() = wf_.sd(ispin, ikp)->c();
    }
  tmap_["expowf_copy"].stop();

  // initialize factors by which each additional term in the Taylor expansion
  // will be multiplied --- factor 1 for the first exponential and factor 2
  // for the second exponential. 
  complex<double> factor1 = 1.0;
  complex<double> factor2 = 1.0;
 
  //order N:
  for(int N = 1; N <= order_; N++)
  {
    factor1 *= -complex<double>(0.0, 1.0)*dt1/double(N);
    factor2 *= -complex<double>(0.0, 1.0)*dt2/double(N);

    if(N != 1)        // for N == 1 this is already done
    {
       tmap_["expowf_copy"].start();
       //move dwf to wf
       for ( int ispin = 0; ispin < wf_.nspin(); ispin++)
          for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
             wf_.sd(ispin, ikp)->c() = dwf->sd(ispin, ikp)->c();
       tmap_["expowf_copy"].stop();
      
       // apply A
       tmap_["expowf_ef"].start();
       ef_.energy(true, *dwf, false, fion, false, sigma);
       tmap_["expowf_ef"].stop();
    }
    
    tmap_["expowf_axpy"].start();
    //accumulate the result
    for ( int ispin = 0; ispin < wf_.nspin(); ispin++)
       for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
       {
          expwf_.sd(ispin, ikp)->c().axpy(factor1, dwf->sd(ispin, ikp)->c());
          newwf_.sd(ispin, ikp)->c().axpy(factor2, dwf->sd(ispin, ikp)->c());
       } 
    tmap_["expowf_axpy"].stop();
    
  }
  
  tmap_["expowf_copy"].start();
  // copy the result back to THE wavefunction, wf_
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++)
  {
    for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
    {
      wf_.sd(ispin, ikp)->c() = expwf_.sd(ispin, ikp)->c();
      s_.hamil_wf->sd(ispin, ikp)->c() = wf_.sd(ispin, ikp)->c();
    }
  }
  tmap_["expowf_copy"].stop();
  
  //if(delete_dwf) delete dwf;
  
}


////////////////////////////////////////////////////////////////////////////////
void ExponentialWavefunctionStepper::preupdate()
{
  if (approximated_)
  {
    // For AETRS, save the potential
    potential_[0] = potential_[1];
    potential_[1] = potential_[2];
    potential_[2] = ef_.get_self_consistent_potential();
    stored_iter_++;
   
    // Propagate wf_ by half of the time step with H(t)
    exponential(0.5*tddt_);
  } 
  else
  {
    // The propagator is U(t + dt, t) = exp(-i dt/2 H(t + dt)) exp(-i dt/2 H(t))
    // In preupdate(), we propagate the wavefunctions (wf_) to psi(t + dt/2) using
    // only part of this propagator, exp(-i dt/2 H(t)) 

    // For ETRS, propagate the wavefunctions by half a time step and by a full time
    // step in one exponential() call. At the end of the call, wavefunctions at time
    // t + dt will be stored in wf_ and wavefunctions at time t + dt/2 will be
    // stored in newwf_. 
    exponential_condensed(tddt_, 0.5*tddt_);
  } 

  if( approximated_ && stored_iter_ >= 3 )
  {
     // If running AETRS, extrapolate the potential using values of the
     // potential from previous stored iterations 
     tmap_["expowf_copy"].start();
     SelfConsistentPotential future_potential;
     future_potential.extrapolate(potential_);
     ef_.set_self_consistent_potential(future_potential);
     tmap_["expowf_copy"].stop();

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

   if ( approximated_  && stored_iter_ < 3 )
   {
      // The following code block should be executed for AETRS only
      
      tmap_["expowf_copy"].start();
      //save the current status
      for ( int ispin = 0; ispin < wf_.nspin(); ispin++)
         for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
            wfhalf_.sd(ispin, ikp)->c() = wf_.sd(ispin, ikp)->c();
      tmap_["expowf_copy"].stop();
    
      // do the rest of the propagation step with H(t)
      exponential(0.5*tddt_, &dwf);

      tmap_["expowf_ef"].start();
      // get the approximated vhxc for the end of the step
      ef_.hamil_cd()->update_density();
      ef_.update_hamiltonian();
      ef_.update_vhxc();
      tmap_["expowf_ef"].stop();

      tmap_["expowf_copy"].start();
      // restore the wf to the middle of the step
      for ( int ispin = 0; ispin < wf_.nspin(); ispin++)
         for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
            wf_.sd(ispin, ikp)->c() = wfhalf_.sd(ispin, ikp)->c();
      tmap_["expowf_copy"].stop();

      // now do the other half of the propagation with H(t + dt)
      exponential(0.5*tddt_);
    
   }
   else if (!approximated_)
   {
     // For ETRS, update the Hamiltonian from H(t) to H(t + dt) given
     // that the wavefunctions in wf_ contain an approximation to
     // psi(t + dt) (psi was propagated to t + dt using H(t).).
     // After updating the Hamiltonian, use H(t + dt) to propagate
     // psi from t + dt/2 (stored initially in newwf_) to t + dt.
     
     // Get the approximated vhxc for the end of the step and update
     // H(t) to H(t + dt)
     tmap_["expowf_ef"].start();
     ef_.hamil_cd()->update_density();
     ef_.update_hamiltonian();
     ef_.update_vhxc();
     tmap_["expowf_ef"].stop();

     // Copy the wavefunctions at t + dt/2 (currently in newwf_) to wf_
     tmap_["expowf_copy"].start();
     for ( int ispin = 0; ispin < wf_.nspin(); ispin++)
         for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
            wf_.sd(ispin, ikp)->c() = newwf_.sd(ispin, ikp)->c();
     tmap_["expowf_copy"].stop();
    
     // Propagate the wavefunctions in wf_ from t + dt/2 to t + dt
     // using H(t + dt)
     exponential(0.5*tddt_); 
   }
   
}
