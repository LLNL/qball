////////////////////////////////////////////////////////////////////////////////  
// Copyright (c) 2016, Lawrence Livermore National Security, LLC. 
// qb@ll:  Qbox at Lawrence Livermore
//
// This file is part of qb@ll.
//
// Produced at the Lawrence Livermore National Laboratory. 
// Written by Xavier Andrade (xavier@llnl.gov).
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
// qbdiff.C:
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <fstream>
#include <vector>
#include <cassert>
#include <sstream>

#include <rapidxml.hpp>

using namespace std;

double energy_tol = 1.0e-6;

template <typename Type>
static Type value(const rapidxml::xml_base<> * node){
  assert(node);
  std::istringstream stst(node->value());
  Type value;
  stst >> value;
  return value;
}

int main(int argc, char** argv)
{

  if(argc != 3){
    cerr << "Usage: " << argv[0] << " file1.xml file2.xml" << endl;
  }

  string filename1 = string(argv[1]);
  string filename2 = string(argv[2]);

  rapidxml::xml_document<> docu1, docu2;
    
  std::ifstream file1(filename1);
  std::vector<char> buffer1((istreambuf_iterator<char>(file1)), istreambuf_iterator<char>());
  docu1.parse<0>(buffer1.data());
  
  std::ifstream file2(filename2);
  std::vector<char> buffer2((istreambuf_iterator<char>(file2)), istreambuf_iterator<char>());
  docu2.parse<0>(buffer2.data());

  rapidxml::xml_node<> * iteration_node1 = docu1.first_node("qbox:simulation")->first_node("run")->first_node("iteration");
  rapidxml::xml_node<> * iteration_node2 = docu2.first_node("qbox:simulation")->first_node("run")->first_node("iteration");

  cout << iteration_node1 << "\t" << iteration_node2 << endl;
  
  if(!iteration_node1 || !iteration_node2){
    cout << "At least one of the files does not contain any iteration." << endl;
    exit(1);
  }

  bool correct = true;
  int iter = 1;
  while(iteration_node1 && iteration_node2){

    cout << "iter " << iter << endl;
    
    bool iter_correct = true;

    assert(iteration_node1->first_node("etotal"));
    assert(iteration_node2->first_node("etotal"));
 
    double etotal1 = value<double>(iteration_node1->first_node("etotal"));
    double etotal2 = value<double>(iteration_node2->first_node("etotal"));

    cout << endl;
    
    cout << "Energy 1   = " << scientific << etotal1 << endl;
    cout << "Energy 2   = " << scientific << etotal2 << endl;
    cout << "Difference = " << scientific << fabs(etotal1 - etotal2) << endl;

    if(fabs(etotal1 - etotal2) > energy_tol){
      iter_correct = false;
    }

    correct = correct && iter_correct;

    cout << endl;

    cout << "Iteration\t" << iter << "\t:\t";
    
    if(iter_correct){
      cout << "[   OK   ]" << endl;
    } else {
      cout << "[  FAIL  ]" << endl;
    }

    iteration_node1 = iteration_node1->next_sibling("iteration");
    iteration_node2 = iteration_node2->next_sibling("iteration");

  }
  
  if(!correct) return 1;

}
