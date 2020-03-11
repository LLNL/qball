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

#include <config.h>

#include "ExponentialWavefunctionStepper.h"
#include "SelfConsistentPotential.h"
#include "SlaterDet.h"
#include "Sample.h"
#include <iostream>
#include <deque>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
ExponentialWavefunctionStepper::ExponentialWavefunctionStepper(Wavefunction& wf, CurrentDensity& currd, double tddt, TimerMap& tmap, EnergyFunctional & ef, Sample & s, bool approximated)
    : tddt_(tddt), WavefunctionStepper(wf,tmap), ef_(ef), s_(s), approximated_(approximated), expwf_(s.wf), wfhalf_(s.wf), newwf_(s.wf), currd_(currd), tempvp_(*ef.vp)
{
  order_ = 4;
  potential_.resize(3);
  stored_iter_ = 0;
  merge_exp_ = false;
}

////////////////////////////////////////////////////////////////////////////////
void ExponentialWavefunctionStepper::exponential(int num_exp, double dt1, double dt2, Wavefunction * dwf){

  // The first input argument num_exp determines if the current call to 
  // exponential requires one or two timestep args and exponential calcs 
  
  // dummy variables to call ef_.energy
  std::vector<std::vector<double> > fion;
  std::valarray<double> sigma;

  // if dwf is not explicitly passed, recreate the dwf object with the
  // ef.energy() call
  if (dwf == 0)
  {
    dwf = &wfhalf_;
    tmap_["expowf_ef"].start();
    ef_.energy(wf_, true, *dwf, false, fion, false, sigma);
    tmap_["expowf_ef"].stop();
  }

  // Expand exp(A)x as a 4th order Taylor series
  // exp(A)x ~= x + Ax + A^2x/2! + A^3x/3! + A^4x/4!
  //          = x + Ax + A(Ax)/2 + A(A^2x/2)/3 + A(A^3x/6)/4 
  
  // order 0
  // For the 0th order term of the Taylor expansion, simply copy the
  // wavefunctions from wf_ to expwf_ for the first exponential calculation
  // and conditionally from wf_ to newwf_ for the second exponential calculation
  tmap_["expowf_copy"].start();
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++)
    for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
      expwf_.sd(ispin, ikp)->c() = wf_.sd(ispin, ikp)->c();
  
  if (num_exp == 2)
  {
    for ( int ispin = 0; ispin < wf_.nspin(); ispin++)
      for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
        newwf_.sd(ispin, ikp)->c() = wf_.sd(ispin, ikp)->c();
  } 
  tmap_["expowf_copy"].stop();

  // initialize factors by which each additional term in the Taylor expansion
  // will be multiplied --- factor 1 for the first exponential and factor 2
  // for the second exponential (only used if a second exp() is requested). 
  complex<double> factor1 = 1.0;
  complex<double> factor2 = 1.0;
 
  //order N:
  for(int N = 1; N <= order_; N++)
  {
    factor1 *= -complex<double>(0.0, 1.0)*dt1/double(N);
    if (num_exp == 2)
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
       ef_.energy(wf_, true, *dwf, false, fion, false, sigma);
       tmap_["expowf_ef"].stop();
    }
    
    // accumulate the result of propagating wf_ by the first 
    // time step (dt1) in expwf_. If a second time step is provided, 
    // store the wavefunctions propagated by that time step in newwf_
    tmap_["expowf_axpy"].start();
    for ( int ispin = 0; ispin < wf_.nspin(); ispin++)
       for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
          expwf_.sd(ispin, ikp)->c().axpy(factor1, dwf->sd(ispin, ikp)->c());
    
    if (num_exp == 2)
    {
      for ( int ispin = 0; ispin < wf_.nspin(); ispin++)
         for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
            newwf_.sd(ispin, ikp)->c().axpy(factor2, dwf->sd(ispin, ikp)->c());
    }
    tmap_["expowf_axpy"].stop();
    
  } // Note: this bracket terminates the for loop from 1 to order_
  
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
  
}


////////////////////////////////////////////////////////////////////////////////
void ExponentialWavefunctionStepper::preupdate()
{
  // allocate memory for the input args to exponential()
  int num_exp;
  double dt1;
  double dt2; 
 
  if (approximated_)
  {
    // For AETRS, save the potential
    potential_[0] = potential_[1];
    potential_[1] = potential_[2];
    potential_[2] = ef_.get_self_consistent_potential();
    stored_iter_++;
  }
 
  // The propagator is U(t + dt, t) = exp(-i dt/2 H(t + dt)) exp(-i dt/2 H(t))
  // In preupdate(), we propagate the wavefunctions (wf_) to psi(t + dt/2) using
  // only part of this propagator, exp(-i dt/2 H(t)). Similarly, an approximation
  // to the wavefunctions, psi(t + dt) are obtained by using the propagator,
  // exp(-i dt H(t)). 

  if (merge_exp_) 
  { 
    // Propagate the wavefunctions by half a time step and by a full time
    // step in one exponential() call. At the end of the call, wavefunctions at time
    // t + dt will be stored in wf_ and wavefunctions at time t + dt/2 will be
    // stored in newwf_.
    num_exp = 2;
    dt1 = tddt_;
    dt2 = 0.5*tddt_;
    exponential(num_exp, dt1, dt2);
  }
  else 
  {
    // Propagate the wavefunctions by half a time step and by a full time step
    // in two separate exponential calls.
    num_exp = 1;
    dt1 = 0.5*tddt_;
    dt2 = 0.0;
    // First, propagate to t + dt/2
    exponential(num_exp, dt1, dt2);
    // Then, backup the wavefunctions at t + dt/2 in newwf_ 
    tmap_["expowf_copy"].start();
    for ( int ispin = 0; ispin < wf_.nspin(); ispin++)
       for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
          newwf_.sd(ispin, ikp)->c() = wf_.sd(ispin, ikp)->c(); 
    tmap_["expowf_copy"].stop();
    // Finally, propagate wavefunctions in wf_ from t + dt/2 to t + dt via
    // a second half step
    exponential(num_exp, dt1, dt2);
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
     //
     // At the time the function call to exponential() is made below, psi(t + dt/2)
     // --- what you want to propagate forward by another half step -- is stored
     // in newwf_. Move this version of the wavefunction to wf_ before 
     // propagation.
     tmap_["expowf_copy"].start();
     for ( int ispin = 0; ispin < wf_.nspin(); ispin++)
         for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
            wf_.sd(ispin, ikp)->c() = newwf_.sd(ispin, ikp)->c(); 
     tmap_["expowf_copy"].stop();
     
     // Redefine time step arguments so that only one exponential is calculated
     // within exponential()
     num_exp = 1;
     dt1 = 0.5*tddt_;
     dt2 = 0.0;
     // Propagate psi(t + dt/2) to psi(t + dt) using H(t + dt)
     exponential(num_exp, dt1, dt2);
  } 

}

////////////////////////////////////////////////////////////////////////////////
void ExponentialWavefunctionStepper::update(Wavefunction& dwf,int time)
{
   // Define time step args for both AETRS and ETRS 
   int num_exp = 1;
   double dt1 = 0.5*tddt_;
   double dt2 = 0.0;

   if ( !approximated_ || stored_iter_ < 3 )
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

     // AK: also update current and vp to estimate at t + dt
     if(ef_.vp) {
       tmap_["current"].start();
       currd_.update_current(ef_, false);
       tmap_["current"].stop(); 

       tempvp_ = *ef_.vp;
       ef_.vp->calculate_acceleration(tddt_, currd_.total_current, wf_.cell());
       // first argument is wrong for laser; need to get time from somewhere
       ef_.vp->propagate(time, tddt_);
     }

     // Copy the wavefunctions at t + dt/2 (currently in newwf_) to wf_
     tmap_["expowf_copy"].start();
     for ( int ispin = 0; ispin < wf_.nspin(); ispin++)
         for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
            wf_.sd(ispin, ikp)->c() = newwf_.sd(ispin, ikp)->c();
     tmap_["expowf_copy"].stop();
    
     // Propagate the wavefunctions in wf_ from t + dt/2 to t + dt
     // using H(t + dt)
     exponential(num_exp, dt1, dt2);

     // AK: restore vp to state at t
     if(ef_.vp) *ef_.vp = tempvp_;

   }
   
}
