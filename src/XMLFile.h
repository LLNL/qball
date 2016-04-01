////////////////////////////////////////////////////////////////////////////////  
// Copyright (c) 2013-2016, Lawrence Livermore National Security, LLC. 
// qb@ll:  Qbox at Lawrence Livermore
//
// This file is part of qb@ll.
//
// Produced at the Lawrence Livermore National Laboratory. 
// Written by Xavier Andrade (xavier@tddft.org).
//
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
// XMLFile.h
//
////////////////////////////////////////////////////////////////////////////////


#include <cassert>
#include <string>
#include <iostream>

class Tag{

public:
  Tag(const std::string & tag){
    tag_ = tag;
  }

  const std::string & name() const {
    return tag_;
  }
  
  std::string start() const {
    return "<" + tag_ + ">";
  }
  
  std::string end() const {
    return "</" + tag_ + ">";
  }

  std::string text(const std::string & buf) const {

    std::string::size_type pos = 0;

    std::string::size_type start_pos = buf.find(start(), pos);

    assert(start_pos != std::string::npos );
    
    start_pos = buf.find(">", start_pos)+1;

    std::string::size_type end_pos = buf.find(end());

    pos = buf.find(">", end_pos) + 1;

    std::string::size_type len = end_pos - start_pos;
    
    return buf.substr(start_pos, len);
  }

  template <typename Type>
  void get_value(const std::string & buf, Type & value) const {
    istd::stringstream stst(text(buf));
    stst >> value;
  }
  
private:
  std::string tag_;
  
};

// Local Variables:
// mode: c++
// End:
