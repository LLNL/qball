#ifndef PSEUDO_SET_HPP
#define PSEUDO_SET_HPP

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
#include <map>
#include <fstream>

#include "element.hpp"

#include <iostream>
#include <dirent.h>
#include "detect_format.hpp"
#include "psml.hpp"
#include "qso.hpp"
#include "upf1.hpp"
#include "upf2.hpp"

namespace pseudopotential {

  class set{

  public:
    
    set(const std::string & dirname){
      DIR * dir = opendir(dirname.c_str());

      struct dirent *ent;
      while ((ent = readdir(dir)) != NULL) {
	const std::string filename(ent->d_name);
	const std::string fullname = dirname + "/" + filename;

	if(filename == "." || filename == "..") continue;

	pseudopotential::format format = detect_format(fullname);
	
	if(format == pseudopotential::format::FILE_NOT_FOUND || format == pseudopotential::format::UNKNOWN) continue;


	// we open the pseudo just to get the species symbol, this could be done in a better way
	pseudopotential::base * pseudo = NULL;
 
	switch(format){
	case pseudopotential::format::QSO:
	  pseudo = new pseudopotential::qso(fullname);
	  break;
	case pseudopotential::format::UPF1:
	  pseudo = new pseudopotential::upf1(fullname, /*uniform_grid = */ true);
	  break;
	case pseudopotential::format::UPF2:
	  pseudo = new pseudopotential::upf2(fullname, /*uniform_grid = */ true);
	  break;
	case pseudopotential::format::PSML:
	  pseudo = new pseudopotential::psml(fullname, /*uniform_grid = */ true);
	  break;
	default:
	  continue;
	}
	
	std::string symbol = pseudo->symbol();
	
	delete pseudo;
	
	element_values vals;

	vals.file_path_ = fullname;
	vals.lmax_ = INVALID_L;
	vals.llocal_ = INVALID_L;
	vals.spacing_ = -1.0;
	vals.radius_ = -1.0;

	map_[symbol] = vals;

      }

    }
    
    bool has(const element & el) const {
      return map_.find(el.symbol()) != map_.end();
    }
    
    const std::string & file_path(const element & el) const {
      return map_.at(el.symbol()).file_path_;
    }
    
    int lmax(const element & el) const {
      return map_.at(el.symbol()).lmax_;
    }
    
    int llocal(const element & el) const {
      return map_.at(el.symbol()).llocal_;
    }
    
    double spacing(const element & el) const {
      return map_.at(el.symbol()).spacing_;
    }
    
    double radius(const element & el) const {
      return map_.at(el.symbol()).radius_;
    }

  private:

    struct element_values {
      std::string file_path_;
      int lmax_;
      int llocal_;
      double spacing_;
      double radius_;
    };
    
    std::map<std::string, element_values> map_;

    
  };

}

#endif
