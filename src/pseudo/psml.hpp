#ifndef PSEUDO_PSML_HPP
#define PSEUDO_PSML_HPP

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

#include <fstream>
#include <vector>
#include <cassert>
#include <sstream>
#include <iostream>

#include "anygrid.hpp"
#include "base.hpp"
#include "chemical_element.hpp"
#include <rapidxml.hpp>

namespace pseudopotential {

  class psml : public pseudopotential::anygrid {

  public:

    psml(const std::string & filename):
      file_(filename),
      buffer_((std::istreambuf_iterator<char>(file_)), std::istreambuf_iterator<char>()){

      buffer_.push_back('\0');
      doc_.parse<0>(&buffer_[0]);

      root_node_ = doc_.first_node("psml");

      spec_node_ = root_node_->first_node("pseudo-atom-spec");
      //some files do not have "pseudo-atom-spec" but "header"
      if(!spec_node_) spec_node_ = root_node_->first_node("header");
      assert(spec_node_);
      
      //now check the type
      bool has_local_potential = root_node_->first_node("local-potential");
      bool has_nl_projectors = root_node_->first_node("nonlocal-projectors");
      bool has_semilocal_potentials = root_node_->first_node("semilocal-potentials");
      bool has_pseudo_wavefunctions = root_node_->first_node("pseudo-wave-functions");
      if(has_nl_projectors && has_local_potential) {
	type_ = pseudopotential::type::KLEINMAN_BYLANDER;
      } else if(has_semilocal_potentials && has_pseudo_wavefunctions) {
	type_ = pseudopotential::type::SEMILOCAL;
      } else {
	throw status::UNSUPPORTED_TYPE;
      }

      {
	//read lmax
	std::string tag1, tag2;
	if(type_ == pseudopotential::type::KLEINMAN_BYLANDER){
	  tag1 = "nonlocal-projectors";
	  tag2 = "proj";
	} else if(type_ == pseudopotential::type::SEMILOCAL){
	  tag1 = "semilocal-potentials";
	  tag2 = "slps";
	} else {
	  assert(false);
	}
      
	lmax_ = -1;
	rapidxml::xml_node<> * node = root_node_->first_node(tag1.c_str());
	assert(node);
	node = node->first_node(tag2.c_str());
	while(node){
	  int read_l = letter_to_l(node->first_attribute("l")->value());
	  lmax_ = std::max(lmax_, read_l);
	  node = node->next_sibling(tag2.c_str());
	}
	assert(lmax_ >= 0);
	assert(lmax_ < 9);
      }
      
      //read grid
      {
	rapidxml::xml_node<> * node = root_node_->first_node("grid");
	
	assert(node);
	
	int size = value<int>(node->first_attribute("npts"));
	grid_.resize(size);
	std::istringstream stst(node->first_node("grid-data")->value());
	for(int ii = 0; ii < size; ii++) {
	  stst >> grid_[ii];
	  if(ii > 0) assert(grid_[ii] > grid_[ii - 1]);
	}

	assert(fabs(grid_[0]) <= 1e-10);

	mesh_size_ = 0;
	for(double rr = 0.0; rr <= grid_[grid_.size() - 1]; rr += mesh_spacing()) mesh_size_++;
	
      }
      
    }

    pseudopotential::format format() const { return pseudopotential::format::PSML; }
    
    int size() const { return buffer_.size(); };

    std::string description() const {
      return "";
    }
    
    std::string symbol() const {
      return spec_node_->first_attribute("atomic-label")->value();
    }

    int atomic_number() const {
      return value<int>(spec_node_->first_attribute("atomic-number"));
    }

    double mass() const {
      chemical_element el(symbol());
      return el.mass();
    }
    
    int valence_charge() const {
      return value<int>(spec_node_->first_node("valence-configuration")->first_attribute("total-valence-charge"));
    }

    int llocal() const {
      return -1;
    }

    pseudopotential::exchange exchange() const {
      // PSML uses libxc ids, so we just need to read the value
      rapidxml::xml_node<> * node = spec_node_->first_node("exchange-correlation")->first_node("libxc-info")->first_node("functional");
      while(node){
	if(value<std::string>(node->first_attribute("type")) == "exchange") {
	  return pseudopotential::exchange(value<int>(node->first_attribute("id")));
	}
	node = node->next_sibling("functional");
      }
      return pseudopotential::exchange::UNKNOWN;
    }

    pseudopotential::correlation correlation() const {
      // PSML uses libxc ids, so we just need to read the value
      rapidxml::xml_node<> * node = spec_node_->first_node("exchange-correlation")->first_node("libxc-info")->first_node("functional");
      while(node){
	if(value<std::string>(node->first_attribute("type")) == "correlation") {
	  return pseudopotential::correlation(value<int>(node->first_attribute("id")));
	}
	node = node->next_sibling("functional");
      }
      return pseudopotential::correlation::UNKNOWN;
    }
    
    int nchannels() const {
      if(type_ == pseudopotential::type::SEMILOCAL) return 1;
      int nc = 0;
      rapidxml::xml_node<> * node = root_node_->first_node("nonlocal-projectors");
      assert(node);
      node = node->first_node("proj");
      while(node){
	int read_ic = value<int>(node->first_attribute("seq")) - 1;
	nc = std::max(nc, read_ic + 1);
	node = node->next_sibling("proj");
      }
      return nc;
    }
    
    void local_potential(std::vector<double> & val) const {
      read_function(root_node_->first_node("local-potential"), val, true);
    }

    int nprojectors() const {
      rapidxml::xml_node<> * node = root_node_->first_node("nonlocal-projectors")->first_node("proj");
      int count = 0;
      while(node) {
	count++;
	node = node->next_sibling("proj");
      }
      return count;
    }
    
    void projector(int l, int ic, std::vector<double> & val) const {
      rapidxml::xml_node<> * node = root_node_->first_node("nonlocal-projectors")->first_node("proj");
      while(node){
	int read_l = letter_to_l(node->first_attribute("l")->value());
	int read_ic = value<int>(node->first_attribute("seq")) - 1;
	if(l == read_l && ic == read_ic) break;
	node = node->next_sibling("proj");
      }
      read_function(node, val);      
    }

    double d_ij(int l, int ic, int jc) const {
      if(ic != jc) return 0.0;
      
      rapidxml::xml_node<> * node = root_node_->first_node("nonlocal-projectors")->first_node("proj");
      while(node){
	int read_l = letter_to_l(node->first_attribute("l")->value());
	int read_ic = value<int>(node->first_attribute("seq")) - 1;
	if(l == read_l && ic == read_ic) break;
	node = node->next_sibling("proj");
      }
      assert(node);
      return value<double>(node->first_attribute("ekb"));
    }

    bool has_radial_function(int l) const{
      return false;
    }

    void radial_function(int l, std::vector<double> & val) const {
      rapidxml::xml_node<> * node = root_node_->first_node("pseudo-wave-functions")->first_node("pswf");
      while(node){
	int read_l = letter_to_l(node->first_attribute("l")->value());
	if(l == read_l) break;
	node = node->next_sibling("pswf");
      }
      read_function(node, val);
    }

    void radial_potential(int l, std::vector<double> & val) const {
      rapidxml::xml_node<> * node = root_node_->first_node("semilocal-potentials")->first_node("slps");
      while(node){
	int read_l = letter_to_l(node->first_attribute("l")->value());
	if(l == read_l) break;
	node = node->next_sibling("slps");
      }
      read_function(node, val);
    }

    bool has_nlcc() const{
      rapidxml::xml_node<> * node = root_node_->first_node("pseudocore-charge");
      return node;
    }

    void nlcc_density(std::vector<double> & val) const {
      read_function(root_node_->first_node("pseudocore-charge"), val);
      for(unsigned ii = 0; ii < val.size(); ii++) val[ii] /= 4.0*M_PI;
    }
    
    bool has_density(){
      return root_node_->first_node("valence-charge");
    }
    
    void density(std::vector<double> & val) const {
      read_function(root_node_->first_node("valence-charge"), val);
      for(unsigned ii = 0; ii < val.size(); ii++) val[ii] /= 4.0*M_PI;
    }
    
  private:

    void read_function(rapidxml::xml_node<> * base_node, std::vector<double> & val, bool potential_padding = false) const{
      assert(base_node);
      rapidxml::xml_node<> * node = base_node->first_node("radfunc")->first_node("data");
      assert(node);
      int size = grid_.size();
      if(node->first_attribute("npts")) size = value<int>(node->first_attribute("npts"));
      val.resize(grid_.size());
      std::istringstream stst(node->value());
      for(int ii = 0; ii < size; ii++) stst >> val[ii];

      if(potential_padding){
	for(unsigned ii = size; ii < grid_.size(); ii++) val[ii] = -valence_charge()/grid_[ii];
      } else {
	for(unsigned ii = size; ii < grid_.size(); ii++) val[ii] = 0.0;
      }
      
      interpolate(val);
    }
    
    //for some stupid reason psml uses letters instead of numbers for angular momentum
    static int letter_to_l(const std::string & letter){
      if(letter == "s") return 0;
      if(letter == "p") return 1;
      if(letter == "d") return 2;
      if(letter == "f") return 3;
      if(letter == "g") return 4;
      if(letter == "h") return 5;
      return -1;
    }
    
    std::ifstream file_;
    std::vector<char> buffer_;
    rapidxml::xml_document<> doc_;
    rapidxml::xml_node<> * root_node_;
    rapidxml::xml_node<> * spec_node_;
    
  };

}

#endif
