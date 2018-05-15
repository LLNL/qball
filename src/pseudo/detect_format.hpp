#ifndef PSEUDO_DETECT_FORMAT_HPP
#define PSEUDO_DETECT_FORMAT_HPP

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

#include <vector>
#include <string>
#include <rapidxml.hpp>
#include "base.hpp"

namespace pseudopotential {

  pseudopotential::format detect_format(const std::string & filename){

    std::ifstream file(filename.c_str());

    if(!file) return pseudopotential::format::FILE_NOT_FOUND;
    
    std::vector<char> buffer((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    buffer.push_back('\0');

    {
      rapidxml::xml_document<> doc;
      bool is_xml = true;
      try {
	doc.parse<0>(&buffer[0]);
      } catch(rapidxml::parse_error xml_error){
	is_xml = false;
      }
      
      if(is_xml){
	if(doc.first_node("fpmd:species")) return pseudopotential::format::QSO;
	if(doc.first_node("PP_INFO")) return pseudopotential::format::UPF1;
	if(doc.first_node("UPF")) return pseudopotential::format::UPF2;
	if(doc.first_node("psml")) return pseudopotential::format::PSML;
      }
    }

    std::string extension = filename.substr(filename.find_last_of(".") + 1);
    std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);

    if(extension == "psf") return pseudopotential::format::PSF;
    if(extension == "cpi") return pseudopotential::format::CPI;
    if(extension == "fhi") return pseudopotential::format::FHI;
    if(extension == "hgh") return pseudopotential::format::HGH;
    
    return pseudopotential::format::UNKNOWN;
  }
  
}

#endif
