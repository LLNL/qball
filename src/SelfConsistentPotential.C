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
// SelfConsistentPotential.C
//
////////////////////////////////////////////////////////////////////////////////

#include "SelfConsistentPotential.h"
#include "EnergyFunctional.h"
#include <complex>
#include <vector>
#include <cassert>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
SelfConsistentPotential::SelfConsistentPotential(const EnergyFunctional & ef):
    v_r(ef.v_r), hamil_rhoelg(ef.hamil_rhoelg),rhoelg(ef.rhoelg), eps_(ef.eps_),
    ehart_(ef.ehart_), exc_(ef.exc_),esr_(ef.esr_)
{
}

////////////////////////////////////////////////////////////////////////////////
void SelfConsistentPotential::extrapolate(const std::vector<SelfConsistentPotential> & previous)
{
   assert(previous.size() == 3);

   v_r.resize(previous[0].v_r.size());
   for(int ispin = 0; ispin < v_r.size(); ispin++)
   {
      v_r[ispin].resize(previous[0].v_r[ispin].size());
      for(int ir = 0; ir < v_r[ispin].size(); ir++)
         v_r[ispin][ir] = 3.0*previous[2].v_r[ispin][ir] - 3.0*previous[1].v_r[ispin][ir] + 1.0*previous[0].v_r[ispin][ir];
   }
   
   hamil_rhoelg.resize(previous[0].hamil_rhoelg.size());
   for(int ig = 0; ig <  hamil_rhoelg.size(); ig++)
   {
      hamil_rhoelg[ig] = 3.0*previous[2].hamil_rhoelg[ig]
          - 3.0*previous[1].hamil_rhoelg[ig] + 1.0*previous[0].hamil_rhoelg[ig];
   }
   
   rhoelg.resize(previous[0].rhoelg.size());
   for(int ig = 0; ig <  rhoelg.size(); ig++)
      rhoelg[ig] = 3.0*previous[2].rhoelg[ig] - 3.0*previous[1].rhoelg[ig] + 1.0*previous[0].rhoelg[ig];
   
   eps_ = 3.0*previous[2].eps_ - 3.0*previous[1].eps_ + 1.0*previous[0].eps_;
   ehart_ = 3.0*previous[2].ehart_ - 3.0*previous[1].ehart_ + 1.0*previous[0].ehart_;
   exc_ = 3.0*previous[2].exc_ - 3.0*previous[1].exc_ + 1.0*previous[0].exc_;
   esr_ = 3.0*previous[2].esr_ - 3.0*previous[1].esr_ + 1.0*previous[0].esr_;
}
