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
// SDCellStepper.C
//
////////////////////////////////////////////////////////////////////////////////

#include "SDCellStepper.h"

////////////////////////////////////////////////////////////////////////////////
void SDCellStepper::compute_new_cell(const valarray<double>& sigma)
{
  //cout << " SDCellStepper::compute_new_cell" << endl;
  // multiply stress by A^T to get dE/da_ij
  valarray<double> deda(9);
  
  const UnitCell& cell = s_.wf.cell();
  const double cell_mass = s_.ctrl.cell_mass;
  
  if ( cell_mass <= 0.0 )
  {
    if ( s_.ctxt_.oncoutpe() )
    {
      cout << "<!-- SDCellStepper::compute_new_cell: cell mass is zero\n"
           << "     cannot update cell -->" << endl;
      return;
    }
  }
  
  // deda = - omega * sigma * A^-T
  cell.compute_deda(sigma,deda);
  
  //cout << " SDCellStepper: cell derivatives before constraints" << endl;
  //for ( int i = 0; i < 9; i++ )
  //  cout << " deda[" << i << "] = " << deda[i] << endl;
  
  string cell_lock = s_.ctrl.cell_lock;
  if ( cell_lock != "OFF" )
  {
    // constraints on the cell derivatives
    if ( cell_lock.find("A") != string::npos )
    {
      // vector A is locked
      deda[0] = deda[1] = deda[2] = 0.0;
    }
    if ( cell_lock.find("B") != string::npos )
    {
      // vector B is locked
      deda[3] = deda[4] = deda[5] = 0.0;
    }
    if ( cell_lock.find("C") != string::npos )
    {
      // vector C is locked
      deda[6] = deda[7] = deda[8] = 0.0;
    }
    
    // Check is cell shape should be preserved (if "S" is present in cell_lock)
    // The only changes allowed are renormalizations of a,b,c
    if ( cell_lock.find("S") != string::npos )
    {
      // projection of d in the direction of e
      D3vector d,e;
      
      d = D3vector(deda[0],deda[1],deda[2]);
      e = cell.a(0) / length(cell.a(0));
      d = (d * e) * e;
      deda[0] = d.x; deda[1] = d.y; deda[2] = d.z;
      
      d = D3vector(deda[3],deda[4],deda[5]);
      e = cell.a(1) / length(cell.a(1));
      d = (d * e) * e;
      deda[3] = d.x; deda[4] = d.y; deda[5] = d.z;
      
      d = D3vector(deda[6],deda[7],deda[8]);
      e = cell.a(2) / length(cell.a(2));
      d = (d * e) * e;
      deda[6] = d.x; deda[7] = d.y; deda[8] = d.z;
    }
    
    if ( cell_lock == "R" )
    {
      // preserve aspect ratio
      // deda must be proportional to A
      // All vectors are rescaled by the same constant
      // rescale cell by coefficient alpha, i.e. 
      // deda = alpha * A
      // where alpha * A is the projection of dE/dA in the direction of A
      // alpha = tr (A^T * deda) / || A ||^2  (matrix scalar product)
      const double *a = cell.amat();
      const double num = a[0]*deda[0] + a[1]*deda[1] + a[2]*deda[2] +
                         a[3]*deda[3] + a[4]*deda[4] + a[5]*deda[5] +
                         a[6]*deda[6] + a[7]*deda[7] + a[8]*deda[8];
      const double denom = a[0]*a[0] + a[1]*a[1] + a[2]*a[2] +
                           a[3]*a[3] + a[4]*a[4] + a[5]*a[5] +
                           a[6]*a[6] + a[7]*a[7] + a[8]*a[8];
      const double alpha = num / denom;
      deda = valarray<double>(a,9);
      deda *= alpha;
    }
    //cout << " SDCellStepper: cell derivatives after constraints" << endl;
    //for ( int i = 0; i < 9; i++ )
    //  cout << " deda[" << i << "] = " << deda[i] << endl;
  }
  const double dt = s_.ctrl.dt;
  const double dt2bym = dt*dt/cell_mass;
  
  // cellp = cell - deda * dt^2 / cell_mass
  D3vector a0p = cell.a(0) - dt2bym * D3vector(deda[0],deda[1],deda[2]);
  D3vector a1p = cell.a(1) - dt2bym * D3vector(deda[3],deda[4],deda[5]);
  D3vector a2p = cell.a(2) - dt2bym * D3vector(deda[6],deda[7],deda[8]);

  // maintain constant volume
  if ( cell_lock.find("V") != string::npos ) {

    int ndofs = 3;
    bool scale_axis[3];
    for (int i=0; i<3; i++) 
      scale_axis[i] = true;

    if ( cell_lock.find("A") != string::npos ) {
      ndofs--;
      scale_axis[0] = false;
    }
    if ( cell_lock.find("B") != string::npos ) {
      ndofs--;
      scale_axis[1] = false;
    }
    if ( cell_lock.find("C") != string::npos ) {
      ndofs--;
      scale_axis[2] = false;
    }
    double newvol = a0p * ( a1p ^ a2p );
    double rescale;
    if (ndofs == 3)
      rescale = cbrt(cell.volume()/newvol);
    else if (ndofs == 2)
      rescale = sqrt(cell.volume()/newvol);
    else if (ndofs == 1)
      rescale = cell.volume()/newvol;
    else // ABCV case --> do nothing
      rescale = 1.0;

    if (scale_axis[0])
      a0p *= rescale;
    if (scale_axis[1])
      a1p *= rescale;
    if (scale_axis[2])
      a2p *= rescale;
  }


  cellp = UnitCell(a0p,a1p,a2p);
  
  if ( cell_lock.find("V") != string::npos ) {
    if ( s_.ctxt_.oncoutpe() )
      cout << "<!-- SDCellStepper:  cell volume = " << cellp.volume() << " -->" << endl;
  }

  //cout << " SDCellStepper::compute_new_cell: cellp: " << endl;
  //cout << cellp;
  
}

////////////////////////////////////////////////////////////////////////////////
void SDCellStepper::update_cell(void)
{
  const UnitCell& cell = s_.wf.cell();
  
  // rescale atomic positions in AtomSet
  
  // r_new = A_new A_old^-1 r_old
  vector<vector<double> > r;
  s_.atoms.get_positions(r,false);
  const double* const ainv = cell.amat_inv();
  const double* const ap = cellp.amat();
  
  double tau[3];
  for ( int is = 0; is < r.size(); is++ )
  {
    // transform r to tau: multiply by A^-1
    const int nais = r[is].size()/3;
    for ( int ia = 0; ia < nais; ia++ )
    {
      // multiply r[is][ia] by A_old^-1, result in tau
      cell.vecmult3x3(cell.amat_inv(),&r[is][3*ia],&tau[0]);
      // multiply tau by A_new, result in r[is][3*ia]
      cellp.vecmult3x3(cellp.amat(),&tau[0],&r[is][3*ia]);
    }
  }
  s_.atoms.set_positions(r,true);
  s_.atoms.set_cell(cellp);
  
  // resize wavefunction and basis sets
  
  //cout << " SDCellStepper::update_cell" << endl;
  s_.wf.resize(cellp,s_.wf.refcell(),s_.wf.ecut());
  if ( s_.wfv != 0 )
  {
    s_.wfv->resize(cellp,s_.wf.refcell(),s_.wf.ecut());
    // s_.wfv->clear();
  }
}
