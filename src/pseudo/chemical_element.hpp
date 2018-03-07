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

namespace pseudopotential {
  
  class chemical_element {

  public:

    chemical_element(const std::string & symbol = "none"){
      this->set(symbol);
    }

    std::string symbol() const;

    double charge() const { return -1.0*z; }
    double mass() const;
    int atomic_number() const { return z; }
    
  private:
  
    char z;
    void set(const std::string & symbol);
  
  };

}

#endif

// Local Variables:
// mode: c++
// coding: utf-8
// End:
