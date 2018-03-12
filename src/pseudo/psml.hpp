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

#include "base.hpp"
#include "chemical_element.hpp"
#include <rapidxml.hpp>

namespace pseudopotential {

  class psml : public pseudopotential::base {

  public:

    psml(const std::string & filename):
      file_(filename),
      buffer_((std::istreambuf_iterator<char>(file_)), std::istreambuf_iterator<char>()){

      buffer_.push_back('\0');
      doc_.parse<0>(&buffer_[0]);

      root_node_ = doc_.first_node("psml");

      type_ = pseudopotential::type::KLEINMAN_BYLANDER;

      //read lmax
      lmax_ = -1;
      for(int l = 0; l <= 10; l++ ){
	if(!has_projectors(l)){
	  lmax_ = l - 1;
	  break;
	}
      }
      assert(lmax_ >= 0);
      assert(lmax_ < 9);

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

    std::string format() const { return "PSML"; }
    
    int size() const { return buffer_.size(); };

    std::string description() const {
      return "";
    }
    
    std::string symbol() const {
      return root_node_->first_node("pseudo-atom-spec")->first_attribute("atomic-label")->value();
    }

    int atomic_number() const {
      return value<int>(root_node_->first_node("pseudo-atom-spec")->first_attribute("atomic-number"));
    }

    double mass() const {
      chemical_element el(symbol());
      return el.mass();
    }
    
    int valence_charge() const {
      return value<int>(root_node_->first_node("pseudo-atom-spec")->first_node("valence-configuration")->first_attribute("total-valence-charge"));
    }

    int llocal() const {
      return -1;
    }

    int nchannels() const {
      int nc = 0;
      rapidxml::xml_node<> * node = root_node_->first_node("nonlocal-projectors")->first_node("proj");
      while(node){
	int read_ic = value<int>(node->first_attribute("seq")) - 1;
	nc = std::max(nc, read_ic + 1);
	node = node->next_sibling("proj");
      }
      return nc;
    }
    
    int nquad() const {
      return 0;
    }

    double rquad() const {
      return 0.0;
    }

    double mesh_spacing() const {
      return 0.01;
    }

    int mesh_size() const {
      return mesh_size_;
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
    
    bool has_projectors(int l) const {
      //note: this function can't use lmax_ or lmax()
      rapidxml::xml_node<> * node = root_node_->first_node("nonlocal-projectors")->first_node("proj");
      while(node){
	int read_l = letter_to_l(node->first_attribute("l")->value());
	if(l == read_l) break;
	node = node->next_sibling("proj");
      }
      return node != NULL;
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

    void radial_function(int l, std::vector<double> & function) const {
    }

    void radial_potential(int l, std::vector<double> & function) const {
    }

    bool has_nlcc() const{
      rapidxml::xml_node<> * node = root_node_->first_node("pseudocore-charge");
      return node;
    }

    void nlcc_density(std::vector<double> & val) const {
      read_function(root_node_->first_node("pseudocore-charge"), val);
      for(unsigned ii = 0; ii < val.size(); ii++) val[ii] /= 4.0*M_PI;
    }
    
    void beta(int index, int & l, std::vector<double> & proj) const {
    }

    void dnm_zero(int nbeta, std::vector<std::vector<double> > & dnm) const {
    }

    bool has_rinner() const {
      return false;
    }
    
    void rinner(std::vector<double> & val) const {
    }

    void qnm(int index, int & l1, int & l2, int & n, int & m, std::vector<double> & val) const {
    }

    void qfcoeff(int index, int ltot, std::vector<double> & val) const {
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
      rapidxml::xml_node<> * node = base_node->first_node("radfunc")->first_node("data");
      assert(node);
      int size = value<int>(node->first_attribute("npts"));
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
    
    void interpolate(std::vector<double> & function) const {
      std::vector<double> function_in_grid = function;
      
      Spline function_spline;
      function_spline.fit(grid_.data(), function_in_grid.data(), function_in_grid.size(), SPLINE_FLAT_BC, SPLINE_NATURAL_BC);
      
      function.clear();
      for(double rr = 0.0; rr <= grid_[grid_.size() - 1]; rr += mesh_spacing()){
	function.push_back(function_spline.value(rr));
      }
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
    std::vector<double> grid_;
    int mesh_size_;
    
  };

}

#endif
