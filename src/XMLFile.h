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

class XMLFile {

public:

  XMLFile(const std::string & buf)
    :buf_(buf), pos_(0) {
  }
  
  class Tag{

  public:

    const std::string & name() const {
      return tag_;
    }
  
    std::string start() const {
      return "<" + tag_ + ">";
    }
  
    std::string end() const {
      return "</" + tag_ + ">";
    }

    std::string text() const {

      std::string::size_type start_pos = xml_file_->buf_.find(start(), xml_file_->pos_);

      assert(start_pos != std::string::npos );
    
      start_pos = xml_file_->buf_.find(">", start_pos)+1;

      std::string::size_type end_pos = xml_file_->buf_.find(end());

      xml_file_->pos_ = xml_file_->buf_.find(">", end_pos) + 1;

      std::string::size_type len = end_pos - start_pos;
    
      return xml_file_->buf_.substr(start_pos, len);
    }

    template <typename Type>
    void get_value(Type & value) const {
      std::istringstream stst(text());
      stst >> value;
    }

    bool exists() const {
      std::string::size_type start_pos = xml_file_->buf_.find(start(), xml_file_->pos_);
      return (start_pos != std::string::npos);
    }
    
  private:

    friend class XMLFile;

    Tag(XMLFile * xml_file, const std::string & tag)
      : xml_file_(xml_file), tag_(tag){
    }
    
    XMLFile * xml_file_;
    std::string tag_;
  
  };

  Tag next_tag(const std::string & tag) {
    return Tag(this, tag);
  }
  
private:

  friend class Tag;
  
  std::string::size_type pos_;
  std::string buf_;

};

// Local Variables:
// mode: c++
// End:
