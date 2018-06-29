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
#include "psp8.hpp"

namespace pseudopotential {

  class set{

  private:

    struct element_values {
      std::string file_path_;
      int lmax_;
      int llocal_;
      double spacing_;
      double radius_;
    };

    typedef std::map<std::string, element_values> element_map;
    
    element_map map_;
    
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

	std::string symbol;
   
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
	case pseudopotential::format::PSP8:
	  pseudo = new pseudopotential::psp8(fullname);
	  break;
	default:
	  //get the symbol from the name
	  for(int ii = 0; ii < 3; ii++){
	    char cc = filename[ii];
	    bool is_letter = (cc >= 'a' && cc <= 'z') || (cc >= 'A' && cc <= 'Z');
	    if(!is_letter) break;
	    symbol.push_back(cc);
	  }
	}

	if(pseudo) symbol = pseudo->symbol();
	
	delete pseudo;
	
	element_values vals;

	vals.file_path_ = fullname;
	vals.lmax_ = INVALID_L;
	vals.llocal_ = INVALID_L;
	vals.spacing_ = -1.0;
	vals.radius_ = -1.0;

	map_[symbol] = vals;

      }

      std::ifstream defaults_file((dirname + "/set_defaults").c_str() );

      if(defaults_file){
	std::string line;

	//first line are comments
	getline(defaults_file, line);

	while(true){
	  std::string symbol;
	  defaults_file >> symbol;
	  if(defaults_file.eof()) break;

	  if(has(symbol)){
	    int z;
	    std::string fname;
	    
	    defaults_file >> fname;
	    defaults_file >> z;
	    defaults_file >> map_[symbol].lmax_;
	    defaults_file >> map_[symbol].llocal_;
	    defaults_file >> map_[symbol].spacing_;
	    defaults_file >> map_[symbol].radius_;
	  }

	  getline(defaults_file, line);
	}
	
	defaults_file.close();
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

    //Iterator interface

    class iterator {

    private:
      element_map::iterator map_it_;

    public:
      iterator(const element_map::iterator & map_it):map_it_(map_it){
      }

      iterator & operator++(){
	++map_it_;
	return *this;
      }

      friend bool operator!=(const iterator & a, const iterator & b){
	return a.map_it_ != b.map_it_;
      }
      
      element operator*(){
	return element(map_it_->first);
      }

    };
    
    iterator begin(){ return iterator(map_.begin()); }
    iterator end(){ return iterator(map_.end()); }
    
  };

}

#endif
