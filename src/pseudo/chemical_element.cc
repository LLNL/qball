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

#include "chemical_element.hpp"

#include <cassert>

namespace pseudopotential {

  void chemical_element::set(const std::string & symbol){

    z = 0;

    if(symbol == "H")  z =  1;
    if(symbol == "He") z =  2;
    if(symbol == "Li") z =  3;
    if(symbol == "Be") z =  4;
    if(symbol == "B")  z =  5;
    if(symbol == "C")  z =  6;
    if(symbol == "N")  z =  7;
    if(symbol == "O")  z =  8;
    if(symbol == "F")  z =  9;
    if(symbol == "Ne") z = 10;
    if(symbol == "Na") z = 11;
    if(symbol == "Mg") z = 12;
    if(symbol == "Al") z = 13;
    if(symbol == "Si") z = 14;
    if(symbol == "P")  z = 15;
    if(symbol == "S")  z = 16;
    if(symbol == "Cl") z = 17;
    if(symbol == "Ar") z = 18;
    if(symbol == "K")  z = 19;
    if(symbol == "Ca") z = 20;
    if(symbol == "Sc") z = 21;
    if(symbol == "Ti") z = 22;
    if(symbol == "V")  z = 23;
    if(symbol == "Cr") z = 24;
    if(symbol == "Mn") z = 25;
    if(symbol == "Fe") z = 26;
    if(symbol == "Co") z = 27;
    if(symbol == "Ni") z = 28;
    if(symbol == "Cu") z = 29;
    if(symbol == "Zn") z = 30;
    if(symbol == "Ga") z = 31;
    if(symbol == "Ge") z = 32;
    if(symbol == "As") z = 33;
    if(symbol == "Se") z = 34;
    if(symbol == "Br") z = 35;
    if(symbol == "Kr") z = 36;
    if(symbol == "Rb") z = 37;
    if(symbol == "Sr") z = 38;
    if(symbol == "Y")  z = 39;
    if(symbol == "Zr") z = 40;
    if(symbol == "Nb") z = 41;
    if(symbol == "Mo") z = 42;
    if(symbol == "Tc") z = 43;
    if(symbol == "Ru") z = 44;
    if(symbol == "Rh") z = 45;
    if(symbol == "Pd") z = 46;
    if(symbol == "Ag") z = 47;
    if(symbol == "Cd") z = 48;
    if(symbol == "In") z = 49;
    if(symbol == "Sn") z = 50;
    if(symbol == "Sb") z = 51;
    if(symbol == "Te") z = 52;
    if(symbol == "I")  z = 53;
    if(symbol == "Xe") z = 54;
    if(symbol == "Cs") z = 55;
    if(symbol == "Ba") z = 56;
    if(symbol == "La") z = 57;
    if(symbol == "Ce") z = 58;
    if(symbol == "Pr") z = 59;
    if(symbol == "Nd") z = 60;
    if(symbol == "Pm") z = 61;
    if(symbol == "Sm") z = 62;
    if(symbol == "Eu") z = 63;
    if(symbol == "Gd") z = 64;
    if(symbol == "Tb") z = 65;
    if(symbol == "Dy") z = 66;
    if(symbol == "Ho") z = 67;
    if(symbol == "Er") z = 68;
    if(symbol == "Tm") z = 69;
    if(symbol == "Yb") z = 70;
    if(symbol == "Lu") z = 71;
    if(symbol == "Hf") z = 72;
    if(symbol == "Ta") z = 73;
    if(symbol == "W")  z = 74;
    if(symbol == "Re") z = 75;
    if(symbol == "Os") z = 76;
    if(symbol == "Ir") z = 77;
    if(symbol == "Pt") z = 78;
    if(symbol == "Au") z = 79;
    if(symbol == "Hg") z = 80;
    if(symbol == "Tl") z = 81;
    if(symbol == "Pb") z = 82;
    if(symbol == "Bi") z = 83;
    if(symbol == "Po") z = 84;
    if(symbol == "At") z = 85;
    if(symbol == "Rn") z = 86;
    if(symbol == "Fr") z = 87;
    if(symbol == "Ra") z = 88;
    if(symbol == "Ac") z = 89;
    if(symbol == "Th") z = 90;
    if(symbol == "Pa") z = 91;
    if(symbol == "U")  z = 92;
    if(symbol == "Np") z = 93;
    if(symbol == "Pu") z = 94;

    assert(z != 0);

  }

  std::string chemical_element::symbol() const{

    if(z == 1)  return "H";
    if(z == 2)  return "He";
    if(z == 3)  return "Li";
    if(z == 4)  return "Be";
    if(z == 5)  return "B";
    if(z == 6)  return "C";
    if(z == 7)  return "N";
    if(z == 8)  return "O";
    if(z == 9)  return "F";
    if(z == 10) return "Ne";
    if(z == 11) return "Na";
    if(z == 12) return "Mg";
    if(z == 13) return "Al";
    if(z == 14) return "Si";
    if(z == 15) return "P";
    if(z == 16) return "S";
    if(z == 17) return "Cl";
    if(z == 18) return "Ar";
    if(z == 19) return "K";
    if(z == 20) return "Ca";
    if(z == 21) return "Sc";
    if(z == 22) return "Ti";
    if(z == 23) return "V";
    if(z == 24) return "Cr";
    if(z == 25) return "Mn";
    if(z == 26) return "Fe";
    if(z == 27) return "Co";
    if(z == 28) return "Ni";
    if(z == 29) return "Cu";
    if(z == 30) return "Zn";
    if(z == 31) return "Ga";
    if(z == 32) return "Ge";
    if(z == 33) return "As";
    if(z == 34) return "Se";
    if(z == 35) return "Br";
    if(z == 36) return "Kr";
    if(z == 37) return "Rb";
    if(z == 38) return "Sr";
    if(z == 39) return "Y ";
    if(z == 40) return "Zr";
    if(z == 41) return "Nb";
    if(z == 42) return "Mo";
    if(z == 43) return "Tc";
    if(z == 44) return "Ru";
    if(z == 45) return "Rh";
    if(z == 46) return "Pd";
    if(z == 47) return "Ag";
    if(z == 48) return "Cd";
    if(z == 49) return "In";
    if(z == 50) return "Sn";
    if(z == 51) return "Sb";
    if(z == 52) return "Te";
    if(z == 53) return "I";
    if(z == 54) return "Xe";
    if(z == 55) return "Cs";
    if(z == 56) return "Ba";
    if(z == 57) return "La";
    if(z == 58) return "Ce";
    if(z == 59) return "Pr";
    if(z == 60) return "Nd";
    if(z == 61) return "Pm";
    if(z == 62) return "Sm";
    if(z == 63) return "Eu";
    if(z == 64) return "Gd";
    if(z == 65) return "Tb";
    if(z == 66) return "Dy";
    if(z == 67) return "Ho";
    if(z == 68) return "Er";
    if(z == 69) return "Tm";
    if(z == 70) return "Yb";
    if(z == 71) return "Lu";
    if(z == 72) return "Hf";
    if(z == 73) return "Ta";
    if(z == 74) return "W";
    if(z == 75) return "Re";
    if(z == 76) return "Os";
    if(z == 77) return "Ir";
    if(z == 78) return "Pt";
    if(z == 79) return "Au";
    if(z == 80) return "Hg";
    if(z == 81) return "Tl";
    if(z == 82) return "Pb";
    if(z == 83) return "Bi";
    if(z == 84) return "Po";
    if(z == 85) return "At";
    if(z == 86) return "Rn";
    if(z == 87) return "Fr";
    if(z == 88) return "Ra";
    if(z == 89) return "Ac";
    if(z == 90) return "Th";
    if(z == 91) return "Pa";
    if(z == 92) return "U";
    if(z == 93) return "Np";
    if(z == 94) return "Pu";

    assert(false);
  }

  double chemical_element::mass() const{

    //obtained from http://www.nist.gov/pml/data/comp.cfm

    double element_mass = -1.0;
    
    if(z == 1 ) element_mass = 1.00784      ;
    if(z == 2 ) element_mass = 4.0026022    ;
    if(z == 3 ) element_mass = 6.938        ;
    if(z == 4 ) element_mass = 9.01218315   ;
    if(z == 5 ) element_mass = 10.806       ;
    if(z == 6 ) element_mass = 12.0096      ;
    if(z == 7 ) element_mass = 14.00643     ;
    if(z == 8 ) element_mass = 15.99903     ;
    if(z == 9 ) element_mass = 18.9984031636;
    if(z == 10) element_mass = 20.17976     ;
    if(z == 11) element_mass = 22.989769282 ;
    if(z == 12) element_mass = 24.304       ;
    if(z == 13) element_mass = 26.98153857  ;
    if(z == 14) element_mass = 28.084       ;
    if(z == 15) element_mass = 30.9737619985;
    if(z == 16) element_mass = 32.059       ;
    if(z == 17) element_mass = 35.446       ;
    if(z == 18) element_mass = 39.9481      ;
    if(z == 19) element_mass = 39.09831     ;
    if(z == 20) element_mass = 40.0784      ;
    if(z == 21) element_mass = 44.9559085   ;
    if(z == 22) element_mass = 47.8671      ;
    if(z == 23) element_mass = 50.94151     ;
    if(z == 24) element_mass = 51.99616     ;
    if(z == 25) element_mass = 54.9380443   ;
    if(z == 26) element_mass = 55.8452      ;
    if(z == 27) element_mass = 58.9331944   ;
    if(z == 28) element_mass = 58.69344     ;
    if(z == 29) element_mass = 63.5463      ;
    if(z == 30) element_mass = 65.382       ;
    if(z == 31) element_mass = 69.7231      ;
    if(z == 32) element_mass = 72.6308      ;
    if(z == 33) element_mass = 74.9215956   ;
    if(z == 34) element_mass = 78.9718      ;
    if(z == 35) element_mass = 79.901       ;
    if(z == 36) element_mass = 83.7982      ;
    if(z == 37) element_mass = 85.46783     ;
    if(z == 38) element_mass = 87.621       ;
    if(z == 39) element_mass = 88.905842    ;
    if(z == 40) element_mass = 91.2242      ;
    if(z == 41) element_mass = 92.906372    ;
    if(z == 42) element_mass = 95.951       ;
    if(z == 43) element_mass = 98           ;
    if(z == 44) element_mass = 101.072      ;
    if(z == 45) element_mass = 102.905502   ;
    if(z == 46) element_mass = 106.421      ;
    if(z == 47) element_mass = 107.86822    ;
    if(z == 48) element_mass = 112.4144     ;
    if(z == 49) element_mass = 114.8181     ;
    if(z == 50) element_mass = 118.7107     ;
    if(z == 51) element_mass = 121.7601     ;
    if(z == 52) element_mass = 127.603      ;
    if(z == 53) element_mass = 126.904473   ;
    if(z == 54) element_mass = 131.2936     ;
    if(z == 55) element_mass = 132.905451966;
    if(z == 56) element_mass = 137.3277     ;
    if(z == 57) element_mass = 138.905477   ;
    if(z == 58) element_mass = 140.1161     ;
    if(z == 59) element_mass = 140.907662   ;
    if(z == 60) element_mass = 144.2423     ;
    if(z == 61) element_mass = 145          ;
    if(z == 62) element_mass = 150.362      ;
    if(z == 63) element_mass = 151.9641     ;
    if(z == 64) element_mass = 157.253      ;
    if(z == 65) element_mass = 158.925352   ;
    if(z == 66) element_mass = 162.5001     ;
    if(z == 67) element_mass = 164.930332   ;
    if(z == 68) element_mass = 167.2593     ;
    if(z == 69) element_mass = 168.934222   ;
    if(z == 70) element_mass = 173.0545     ;
    if(z == 71) element_mass = 174.96681    ;
    if(z == 72) element_mass = 178.492      ;
    if(z == 73) element_mass = 180.947882   ;
    if(z == 74) element_mass = 183.841      ;
    if(z == 75) element_mass = 186.2071     ;
    if(z == 76) element_mass = 190.233      ;
    if(z == 77) element_mass = 192.2173     ;
    if(z == 78) element_mass = 195.0849     ;
    if(z == 79) element_mass = 196.9665695  ;
    if(z == 80) element_mass = 200.5923     ;
    if(z == 81) element_mass = 204.382      ;
    if(z == 82) element_mass = 207.21       ;
    if(z == 83) element_mass = 208.980401   ;
    if(z == 84) element_mass = 209          ;
    if(z == 85) element_mass = 210          ;
    if(z == 86) element_mass = 222          ;
    if(z == 87) element_mass = 223          ;
    if(z == 88) element_mass = 226          ;
    if(z == 89) element_mass = 227          ;
    if(z == 90) element_mass = 232.03774    ;
    if(z == 91) element_mass = 231.035882   ;
    if(z == 92) element_mass = 238.028913   ;
    if(z == 93) element_mass = 237          ;
    if(z == 94) element_mass = 244          ;

    assert(element_mass > 0.0);
    
    return element_mass;
  }

}


// Local Variables:
// mode: c++
// coding: utf-8
// End:
