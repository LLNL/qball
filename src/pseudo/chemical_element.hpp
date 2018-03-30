#ifndef PSEUDO_CHEMICALELEMENT_HPP
#define PSEUDO_CHEMICALELEMENT_HPP

/*
 Copyright (C) 2018 Xavier Andrade

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
  
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
  
 You should have received a copy of the GNU Lesser General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <string>
#include <cassert>
#include <cctype>
#include <iostream>

namespace pseudopotential {
  
  class chemical_element {

  public:

    chemical_element(const std::string & symbol = "none"):symbol_(symbol){
      trim(symbol_);
      symbol_[0] = std::toupper(symbol_[0]);
      for(unsigned ii = 1; ii < symbol_.size(); ii++) symbol_[ii] = std::tolower(symbol_[ii]);
    }

    double charge() const {
      return -1.0*atomic_number();
    }

    std::string symbol() const{
      return symbol_;
    }
    
    int atomic_number() const {

      int z = 0;
      
      if(symbol_ == "H")  z =  1;
      if(symbol_ == "D")  z =  1;
      if(symbol_ == "T")  z =  1;
      if(symbol_ == "He") z =  2;
      if(symbol_ == "Li") z =  3;
      if(symbol_ == "Be") z =  4;
      if(symbol_ == "B")  z =  5;
      if(symbol_ == "C")  z =  6;
      if(symbol_ == "N")  z =  7;
      if(symbol_ == "O")  z =  8;
      if(symbol_ == "F")  z =  9;
      if(symbol_ == "Ne") z = 10;
      if(symbol_ == "Na") z = 11;
      if(symbol_ == "Mg") z = 12;
      if(symbol_ == "Al") z = 13;
      if(symbol_ == "Si") z = 14;
      if(symbol_ == "P")  z = 15;
      if(symbol_ == "S")  z = 16;
      if(symbol_ == "Cl") z = 17;
      if(symbol_ == "Ar") z = 18;
      if(symbol_ == "K")  z = 19;
      if(symbol_ == "Ca") z = 20;
      if(symbol_ == "Sc") z = 21;
      if(symbol_ == "Ti") z = 22;
      if(symbol_ == "V")  z = 23;
      if(symbol_ == "Cr") z = 24;
      if(symbol_ == "Mn") z = 25;
      if(symbol_ == "Fe") z = 26;
      if(symbol_ == "Co") z = 27;
      if(symbol_ == "Ni") z = 28;
      if(symbol_ == "Cu") z = 29;
      if(symbol_ == "Zn") z = 30;
      if(symbol_ == "Ga") z = 31;
      if(symbol_ == "Ge") z = 32;
      if(symbol_ == "As") z = 33;
      if(symbol_ == "Se") z = 34;
      if(symbol_ == "Br") z = 35;
      if(symbol_ == "Kr") z = 36;
      if(symbol_ == "Rb") z = 37;
      if(symbol_ == "Sr") z = 38;
      if(symbol_ == "Y")  z = 39;
      if(symbol_ == "Zr") z = 40;
      if(symbol_ == "Nb") z = 41;
      if(symbol_ == "Mo") z = 42;
      if(symbol_ == "Tc") z = 43;
      if(symbol_ == "Ru") z = 44;
      if(symbol_ == "Rh") z = 45;
      if(symbol_ == "Pd") z = 46;
      if(symbol_ == "Ag") z = 47;
      if(symbol_ == "Cd") z = 48;
      if(symbol_ == "In") z = 49;
      if(symbol_ == "Sn") z = 50;
      if(symbol_ == "Sb") z = 51;
      if(symbol_ == "Te") z = 52;
      if(symbol_ == "I")  z = 53;
      if(symbol_ == "Xe") z = 54;
      if(symbol_ == "Cs") z = 55;
      if(symbol_ == "Ba") z = 56;
      if(symbol_ == "La") z = 57;
      if(symbol_ == "Ce") z = 58;
      if(symbol_ == "Pr") z = 59;
      if(symbol_ == "Nd") z = 60;
      if(symbol_ == "Pm") z = 61;
      if(symbol_ == "Sm") z = 62;
      if(symbol_ == "Eu") z = 63;
      if(symbol_ == "Gd") z = 64;
      if(symbol_ == "Tb") z = 65;
      if(symbol_ == "Dy") z = 66;
      if(symbol_ == "Ho") z = 67;
      if(symbol_ == "Er") z = 68;
      if(symbol_ == "Tm") z = 69;
      if(symbol_ == "Yb") z = 70;
      if(symbol_ == "Lu") z = 71;
      if(symbol_ == "Hf") z = 72;
      if(symbol_ == "Ta") z = 73;
      if(symbol_ == "W")  z = 74;
      if(symbol_ == "Re") z = 75;
      if(symbol_ == "Os") z = 76;
      if(symbol_ == "Ir") z = 77;
      if(symbol_ == "Pt") z = 78;
      if(symbol_ == "Au") z = 79;
      if(symbol_ == "Hg") z = 80;
      if(symbol_ == "Tl") z = 81;
      if(symbol_ == "Pb") z = 82;
      if(symbol_ == "Bi") z = 83;
      if(symbol_ == "Po") z = 84;
      if(symbol_ == "At") z = 85;
      if(symbol_ == "Rn") z = 86;
      if(symbol_ == "Fr") z = 87;
      if(symbol_ == "Ra") z = 88;
      if(symbol_ == "Ac") z = 89;
      if(symbol_ == "Th") z = 90;
      if(symbol_ == "Pa") z = 91;
      if(symbol_ == "U")  z = 92;
      if(symbol_ == "Np") z = 93;
      if(symbol_ == "Pu") z = 94;

      assert(z > 0);

      return z;
    }
    
    double mass() const{

      //obtained from http://www.nist.gov/pml/data/comp.cfm

      double element_mass = -1.0;
    
      if(symbol_ == "H")  element_mass = 1.00784      ;
      if(symbol_ == "D")  element_mass = 2.01410177812;
      if(symbol_ == "T")  element_mass = 3.0160492779 ;
      if(symbol_ == "He") element_mass = 4.0026022    ;
      if(symbol_ == "Li") element_mass = 6.938        ;
      if(symbol_ == "Be") element_mass = 9.01218315   ;
      if(symbol_ == "B")  element_mass = 10.806       ;
      if(symbol_ == "C")  element_mass = 12.0096      ;
      if(symbol_ == "N")  element_mass = 14.00643     ;
      if(symbol_ == "O")  element_mass = 15.99903     ;
      if(symbol_ == "F")  element_mass = 18.9984031636;
      if(symbol_ == "Ne") element_mass = 20.17976     ;
      if(symbol_ == "Na") element_mass = 22.989769282 ;
      if(symbol_ == "Mg") element_mass = 24.304       ;
      if(symbol_ == "Al") element_mass = 26.98153857  ;
      if(symbol_ == "Si") element_mass = 28.084       ;
      if(symbol_ == "P")  element_mass = 30.9737619985;
      if(symbol_ == "S")  element_mass = 32.059       ;
      if(symbol_ == "Cl") element_mass = 35.446       ;
      if(symbol_ == "Ar") element_mass = 39.9481      ;
      if(symbol_ == "K")  element_mass = 39.09831     ;
      if(symbol_ == "Ca") element_mass = 40.0784      ;
      if(symbol_ == "Sc") element_mass = 44.9559085   ;
      if(symbol_ == "Ti") element_mass = 47.8671      ;
      if(symbol_ == "V")  element_mass = 50.94151     ;
      if(symbol_ == "Cr") element_mass = 51.99616     ;
      if(symbol_ == "Mn") element_mass = 54.9380443   ;
      if(symbol_ == "Fe") element_mass = 55.8452      ;
      if(symbol_ == "Co") element_mass = 58.9331944   ;
      if(symbol_ == "Ni") element_mass = 58.69344     ;
      if(symbol_ == "Cu") element_mass = 63.5463      ;
      if(symbol_ == "Zn") element_mass = 65.382       ;
      if(symbol_ == "Ga") element_mass = 69.7231      ;
      if(symbol_ == "Ge") element_mass = 72.6308      ;
      if(symbol_ == "As") element_mass = 74.9215956   ;
      if(symbol_ == "Se") element_mass = 78.9718      ;
      if(symbol_ == "Br") element_mass = 79.901       ;
      if(symbol_ == "Kr") element_mass = 83.7982      ;
      if(symbol_ == "Rb") element_mass = 85.46783     ;
      if(symbol_ == "Sr") element_mass = 87.621       ;
      if(symbol_ == "Y")  element_mass = 88.905842    ;
      if(symbol_ == "Zr") element_mass = 91.2242      ;
      if(symbol_ == "Nb") element_mass = 92.906372    ;
      if(symbol_ == "Mo") element_mass = 95.951       ;
      if(symbol_ == "Tc") element_mass = 98           ;
      if(symbol_ == "Ru") element_mass = 101.072      ;
      if(symbol_ == "Rh") element_mass = 102.905502   ;
      if(symbol_ == "Pd") element_mass = 106.421      ;
      if(symbol_ == "Ag") element_mass = 107.86822    ;
      if(symbol_ == "Cd") element_mass = 112.4144     ;
      if(symbol_ == "In") element_mass = 114.8181     ;
      if(symbol_ == "Sn") element_mass = 118.7107     ;
      if(symbol_ == "Sb") element_mass = 121.7601     ;
      if(symbol_ == "Te") element_mass = 127.603      ;
      if(symbol_ == "I")  element_mass = 126.904473   ;
      if(symbol_ == "Xe") element_mass = 131.2936     ;
      if(symbol_ == "Cs") element_mass = 132.905451966;
      if(symbol_ == "Ba") element_mass = 137.3277     ;
      if(symbol_ == "La") element_mass = 138.905477   ;
      if(symbol_ == "Ce") element_mass = 140.1161     ;
      if(symbol_ == "Pr") element_mass = 140.907662   ;
      if(symbol_ == "Nd") element_mass = 144.2423     ;
      if(symbol_ == "Pm") element_mass = 145          ;
      if(symbol_ == "Sm") element_mass = 150.362      ;
      if(symbol_ == "Eu") element_mass = 151.9641     ;
      if(symbol_ == "Gd") element_mass = 157.253      ;
      if(symbol_ == "Tb") element_mass = 158.925352   ;
      if(symbol_ == "Dy") element_mass = 162.5001     ;
      if(symbol_ == "Ho") element_mass = 164.930332   ;
      if(symbol_ == "Er") element_mass = 167.2593     ;
      if(symbol_ == "Tm") element_mass = 168.934222   ;
      if(symbol_ == "Yb") element_mass = 173.0545     ;
      if(symbol_ == "Lu") element_mass = 174.96681    ;
      if(symbol_ == "Hf") element_mass = 178.492      ;
      if(symbol_ == "Ta") element_mass = 180.947882   ;
      if(symbol_ == "W")  element_mass = 183.841      ;
      if(symbol_ == "Re") element_mass = 186.2071     ;
      if(symbol_ == "Os") element_mass = 190.233      ;
      if(symbol_ == "Ir") element_mass = 192.2173     ;
      if(symbol_ == "Pt") element_mass = 195.0849     ;
      if(symbol_ == "Au") element_mass = 196.9665695  ;
      if(symbol_ == "Hg") element_mass = 200.5923     ;
      if(symbol_ == "Tl") element_mass = 204.382      ;
      if(symbol_ == "Pb") element_mass = 207.21       ;
      if(symbol_ == "Bi") element_mass = 208.980401   ;
      if(symbol_ == "Po") element_mass = 209          ;
      if(symbol_ == "At") element_mass = 210          ;
      if(symbol_ == "Rn") element_mass = 222          ;
      if(symbol_ == "Fr") element_mass = 223          ;
      if(symbol_ == "Ra") element_mass = 226          ;
      if(symbol_ == "Ac") element_mass = 227          ;
      if(symbol_ == "Th") element_mass = 232.03774    ;
      if(symbol_ == "Pa") element_mass = 231.035882   ;
      if(symbol_ == "U")  element_mass = 238.028913   ;
      if(symbol_ == "Np") element_mass = 237          ;
      if(symbol_ == "Pu") element_mass = 244          ;

      assert(element_mass > 0.0);
    
      return element_mass;
    }

   
  private:

    static std::string & ltrim(std::string& str, const std::string& chars = "\t\n\v\f\r "){
      str.erase(0, str.find_first_not_of(chars));
      return str;
    }
 
    static std::string & rtrim(std::string& str, const std::string& chars = "\t\n\v\f\r "){
      str.erase(str.find_last_not_of(chars) + 1);
      return str;
    }
 
    static std::string & trim(std::string& str, const std::string& chars = "\t\n\v\f\r "){
      return ltrim(rtrim(str, chars), chars);
    }
    
    std::string symbol_;
  
  };

}

#endif

// Local Variables:
// mode: c++
// coding: utf-8
// End:
