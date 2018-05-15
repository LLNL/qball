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


namespace pseudopotential {

  class set{

  public:
    
    set(const std::string & filename){
      std::ifstream file(filename);

      if(!file){
	std::cerr << "Internal error: cannot open file '" << filename << "'." << std::endl;
	exit(1);
      }
	
      std::string line, symbol;
      element_values vals;
      int zz;
	
      getline(file, line);

      while(true){
	file >> symbol >> vals.file_path_ >> zz >> vals.lmax_ >> vals.llocal_ >> vals.spacing_ >> vals.radius_;

	if(file.eof()) break;

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
