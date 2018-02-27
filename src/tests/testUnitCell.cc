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
// testUnitCell.C
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#include <qball/UnitCell.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
using namespace std;

int main()
{
  UnitCell cell;
  
  cout << " orthorhombic cell: " << endl;
  cell.set(D3vector(4.,0.,0.),D3vector(0.,5.,0.),D3vector(0.,0.,7.));
  cout << cell << endl;
  
  D3vector v;
  
  const int n = 100000;
  const double d = 10.0;
  int count;
  
  count = 0;
  for ( int i = 0; i < n; i++ )
  {
    // generate a random vector in a box of size d x d x d
    double x0 = drand48() - 0.5;
    double x1 = drand48() - 0.5;
    double x2 = drand48() - 0.5;
    v = d * D3vector(x0,x1,x2);
    
    if ( cell.in_ws(v) )
      count++;
    cell.fold_in_ws(v);
    if ( !cell.in_ws(v) )
      cout << v << endl;
  }
  cout << " volume ratio: " << cell.volume() / (d*d*d)
       << "  computed: " << count/((double) n) << endl << endl;;
  
  cout << " FCC cell: " << endl;
  cell.set(D3vector(4,4,0),D3vector(0,4,4),D3vector(4,0,4));  
  cout << cell << endl;
  
  count = 0;
  for ( int i = 0; i < n; i++ )
  {
    // generate a random vector in a box of size d x d x d
    double x0 = drand48() - 0.5;
    double x1 = drand48() - 0.5;
    double x2 = drand48() - 0.5;
    v = d * D3vector(x0,x1,x2);
    
    if ( cell.in_ws(v) )
      count++;
    cell.fold_in_ws(v);
    if ( !cell.in_ws(v) )
      cout << v << endl;
  }
  cout << " volume ratio: " << cell.volume() / (d*d*d)
       << "  computed: " << count/((double) n) << endl << endl;;
  
  cout << " BCC cell: " << endl;
  cell.set(D3vector(4,4,-4),D3vector(-4, 4,4),D3vector(4,-4,4));  
  cout << cell << endl;
  
  count = 0;
  for ( int i = 0; i < n; i++ )
  {
    // generate a random vector in a box of size d x d x d
    double x0 = drand48() - 0.5;
    double x1 = drand48() - 0.5;
    double x2 = drand48() - 0.5;
    v = d * D3vector(x0,x1,x2);
    
    if ( cell.in_ws(v) )
      count++;
    cell.fold_in_ws(v);
    if ( !cell.in_ws(v) )
      cout << v << endl;
  }
  cout << " volume ratio: " << cell.volume() / (d*d*d)
       << "  computed: " << count/((double) n) << endl << endl;;
  
  cout << " Monoclinic cell: " << endl;
  cell.set(D3vector(4,0,0.1),D3vector(0.2, 4,0),D3vector(0.1,0.3,4));  
  cout << cell << endl;
  
  // Check matrix multiplication function
  double ztest[9];
  cell.matmult3x3(cell.amat(),cell.amat_inv(),ztest);
  for ( int i = 0; i < 9; i++ )
    cout << " ztest[" << i << "]=" << ztest[i] << endl;
  
  count = 0;
  for ( int i = 0; i < n; i++ )
  {
    // generate a random vector in a box of size d x d x d
    double x0 = drand48() - 0.5;
    double x1 = drand48() - 0.5;
    double x2 = drand48() - 0.5;
    v = d * D3vector(x0,x1,x2);
    
    if ( cell.in_ws(v) )
      count++;
    cell.fold_in_ws(v);
    if ( !cell.in_ws(v) )
      cout << v << endl;
  }
  cout << " volume ratio: " << cell.volume() / (d*d*d)
       << "  computed: " << count/((double) n) << endl << endl;
       
  // test UnitCell::included_in function
  UnitCell rcell;
  rcell.set(D3vector(4,4,0),D3vector(0,4,4),D3vector(4,0,4));  
  cell.set(D3vector(4,5,0),D3vector(0,4,4),D3vector(4,0,4));
  cout << " rcell: " << rcell << endl;
  cout << " cell: " << cell << endl;
  cout << " cell is";
  if ( !rcell.encloses(cell) ) cout << " not";
  cout << " included in rcell" << endl;
  
  rcell.set(D3vector(4,4,0),D3vector(0,4,4),D3vector(4,0,4));  
  cout << " rcell: " << rcell << endl;
  cell.set(D3vector(3.8,3.8,0),D3vector(0,3.8,3.8),D3vector(3.8,0,3.8));
  cout << " cell: " << cell << endl;
  cout << " cell is";
  if ( !rcell.encloses(cell) ) cout << " not";
  cout << " included in rcell" << endl; 
}
