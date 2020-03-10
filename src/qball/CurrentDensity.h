////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2014 The Regents of the University of California
//
// This file is part of Qbox
//
// Qbox is distributed under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 2 of
// the License, or (at your option) any later version.
// See the file COPYING in the root directory of this distribution
// or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
// CurrentDensity.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: CurrentDensity.h,v 1.13 2008-09-08 15:56:18 fgygi Exp $

#ifndef CURRENTDENSITY_H
#define CURRENTDENSITY_H

#include "ChargeDensity.h"
#include "FourierTransform.h"
#include "Wavefunction.h"
#include "EnergyFunctional.h"
#include "Sample.h"

class CurrentDensity : private ChargeDensity
{
  private:

  const Wavefunction& wf_;

  public:
  
  std::vector<std::vector<std::vector<double> > > current;
  D3vector total_current;
 
  CurrentDensity (const Sample& s, const Wavefunction & wf);
  ~CurrentDensity (){
  }

  void update_current(EnergyFunctional & energy, bool output=true);

  void plot(const Sample *, const std::string &);
  void plot_vtk(const Sample *, const std::string &);

};

#endif
