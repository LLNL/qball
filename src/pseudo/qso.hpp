#ifndef PSEUDO_QSO_HPP
#define PSEUDO_QSO_HPP

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
#include <rapidxml.hpp>

namespace pseudopotential {

  template <typename Type>
  static Type value(const rapidxml::xml_base<> * node){
    assert(node);
    std::istringstream stst(node->value());
    Type value;
    stst >> value;
    return value;
  }

  class qso : public pseudopotential::base {

  public:

    qso(const std::string & filename):
      file_(filename),
      buffer_((std::istreambuf_iterator<char>(file_)), std::istreambuf_iterator<char>()){

      buffer_.push_back('\0');
      doc_.parse<0>(&buffer_[0]);

      root_node_ = doc_.first_node("fpmd:species");

      pseudo_node_ = root_node_->first_node("ultrasoft_pseudopotential");
      if(pseudo_node_) type_ = type::ULTRASOFT;

      if(!pseudo_node_){
	pseudo_node_ = root_node_->first_node("norm_conserving_semilocal_pseudopotential");
	if(pseudo_node_) type_ = type::KLEINMAN_BYLANDER;
      }

      if(!pseudo_node_){
	pseudo_node_ = root_node_->first_node("norm_conserving_pseudopotential");
	if(pseudo_node_) type_ = type::NORM_CONSERVING;
      }

      assert(pseudo_node_);

      //read lmax
      lmax_ = -1;
      if(pseudo_node_->first_node("lmax")) {
	lmax_ = value<int>(pseudo_node_->first_node("lmax"));
      } else {
	for(int l = 0; l <= 10; l++ ){
	  if(!has_projectors(l)){
	    lmax_ = l - 1;
	    break;
	  }
	}
      }
      assert(lmax_ >= 0);
      assert(lmax_ < 9);
      
    }

    std::string format() const { return "quantum-simulation.org (XML)"; }
    
    int size() const { return buffer_.size(); };

    std::string description() const {
      return root_node_->first_node("description")->value();
    }
    
    std::string symbol() const {
      return root_node_->first_node("symbol")->value();
    }

    int atomic_number() const {
      return value<int>(root_node_->first_node("atomic_number"));
    }

    double mass() const {
      return value<double>(root_node_->first_node("mass"));
    }
    
    int valence_charge() const {
      return value<int>(pseudo_node_->first_node("valence_charge"));
    }

    int llocal() const {
      if(pseudo_node_->first_node("llocal")){
	return value<int>(pseudo_node_->first_node("llocal"));
      } else {
	return -1;
      }
    }

    int nchannels() const {
      if(type_ == type::ULTRASOFT){
	int np = nbeta();
	int nl = lmax() + 1;
	assert(np%nl == 0);
	return np/nl;
      }
      if(type_ == type::KLEINMAN_BYLANDER) return 2;
      return 1;
    }
    
    int nquad() const {
      return value<int>(pseudo_node_->first_node("nquad"));
    }

    double rquad() const {
      return value<double>(pseudo_node_->first_node("rquad"));
    }

    double mesh_spacing() const {
      return value<double>(pseudo_node_->first_node("mesh_spacing"));
    }

    int mesh_size() const {
      rapidxml::xml_node<> * node = pseudo_node_->first_node("local_potential"); //kleinman bylander
      if(!node) node = pseudo_node_->first_node("vlocal"); //ultrasoft
      if(!node) node = pseudo_node_->first_node("projector"); //norm conserving
      assert(node);
      return value<int>(node->first_attribute("size"));
    }
    
    void local_potential(std::vector<double> & potential) const {
      rapidxml::xml_node<> * node = pseudo_node_->first_node("local_potential");
      if(!node){
	// for ultrasoft, this is called vlocal
	node = pseudo_node_->first_node("vlocal");
      }
      assert(node);
      int size = value<int>(node->first_attribute("size"));
      potential.resize(size);
      std::istringstream stst(node->value());
      for(int ii = 0; ii < size; ii++){
	stst >> potential[ii];
      }
    }

    int nprojectors() const {
      switch(type_){
      case type::ULTRASOFT:
	return nbeta();
      case type::KLEINMAN_BYLANDER: {
	int count = 0;
	rapidxml::xml_node<> * node = pseudo_node_->first_node("projector");
	while(node) {
	  count++;
	  node = node->next_sibling("projector");
	}
	return count;
      }
      case type::NORM_CONSERVING:
	return 0;
      }
      return 0;
    }
    
    bool has_projectors(int l) const {
      rapidxml::xml_node<> * node = pseudo_node_->first_node("projector");
      while(node){
	int read_l = value<int>(node->first_attribute("l"));
	if(l == read_l) break;
	node = node->next_sibling("projector");
      }
      return node != NULL;
    }
    
    void projector(int l, int i, std::vector<double> & proj) const {

      std::string tag = "projector";
      
      rapidxml::xml_node<> * node = pseudo_node_->first_node(tag.c_str());

      while(node){
	int read_l = value<int>(node->first_attribute("l"));
	int read_i = value<int>(node->first_attribute("i")) - 1;
	if(l == read_l && i == read_i) break;
	node = node->next_sibling(tag.c_str());
      }

      assert(node != NULL);

      int size = value<int>(node->first_attribute("size"));
      proj.resize(size);
      std::istringstream stst(node->value());
      for(int ii = 0; ii < size; ii++) stst >> proj[ii];
      
    }
    
    double d_ij(int l, int i, int j) const {
      rapidxml::xml_node<> * node = pseudo_node_->first_node("d_ij");
      
      while(node){
	int read_l = value<int>(node->first_attribute("l"));
	int read_i = value<int>(node->first_attribute("i")) - 1;
	int read_j = value<int>(node->first_attribute("j")) - 1;
	if(l == read_l && i == read_i && j == read_j) break;
	node = node->next_sibling("d_ij");
      }

      assert(node != NULL);

      return value<double>(node);
      
    }

    bool has_radial_function(int l) const{
      rapidxml::xml_node<> * node = pseudo_node_->first_node("projector");
      
      while(node){
	int read_l = value<int>(node->first_attribute("l"));
	if(l == read_l) break;
	node = node->next_sibling("projector");
      }
      
      return node->first_node("radial_function") != NULL;
      
    }

    void radial_function(int l, std::vector<double> & function) const {
      rapidxml::xml_node<> * node = pseudo_node_->first_node("projector");

      while(node){
	int read_l = value<int>(node->first_attribute("l"));
	if(l == read_l) break;
	node = node->next_sibling("projector");
      }

      assert(node != NULL);
      assert(node->first_node("radial_function"));
      
      int size = value<int>(node->first_attribute("size"));
      function.resize(size);
      std::istringstream stst(node->first_node("radial_function")->value());
      for(int ii = 0; ii < size; ii++){
	stst >> function[ii];
      }
      
    }

    void radial_potential(int l, std::vector<double> & function) const {
      rapidxml::xml_node<> * node = pseudo_node_->first_node("projector");

      while(node){
	int read_l = value<int>(node->first_attribute("l"));
	if(l == read_l) break;
	node = node->next_sibling("projector");
      }

      assert(node != NULL);
      assert(node->first_node("radial_potential"));
      
      int size = value<int>(node->first_attribute("size"));
      function.resize(size);
      std::istringstream stst(node->first_node("radial_potential")->value());
      for(int ii = 0; ii < size; ii++){
	stst >> function[ii];
      }
      
    }

    bool has_nlcc() const{
      return pseudo_node_->first_node("rho_nlcc");
    }

    void nlcc_density(std::vector<double> & density) const {
      rapidxml::xml_node<> * node = pseudo_node_->first_node("rho_nlcc");
      assert(node);
      int size = value<int>(node->first_attribute("size"));
      density.resize(size);
      std::istringstream stst(node->value());
      for(int ii = 0; ii < size; ii++){
	stst >> density[ii];
      }
    }
    
    void beta(int index, int & l, std::vector<double> & proj) const {
      rapidxml::xml_node<> * node = pseudo_node_->first_node("beta");

      for(int i = 0; i < index; i++) node = node->next_sibling("beta");

      assert(node != NULL);

      l = value<int>(node->first_attribute("l"));
      int size = value<int>(node->first_attribute("size"));
      proj.resize(size);
      std::istringstream stst(node->value());
      for(int ii = 0; ii < size; ii++){
	stst >> proj[ii];
      }
      
    }

    void dnm_zero(int nbeta, std::vector<std::vector<double> > & dnm) const {
      rapidxml::xml_node<> * node = pseudo_node_->first_node("dnm_zero");
      std::istringstream stst(node->value());

      dnm.resize(nbeta);
      for(int i = 0; i < nbeta; i++){
	dnm[i].resize(nbeta);
	for ( int j = 0; j < nbeta; j++){
	  stst >> dnm[i][j];
	}
      }
    }

    bool has_rinner() const {
      rapidxml::xml_node<> * node = pseudo_node_->first_node("rinner");
      return node;
    }
    
    void rinner(std::vector<double> & val) const {
      rapidxml::xml_node<> * node = pseudo_node_->first_node("rinner");

      assert(node != NULL);

      int size = value<int>(node->first_attribute("size"));
      val.resize(size);
      std::istringstream stst(node->value());
      for(int ii = 0; ii < size; ii++){
	stst >> val[ii];
      }
    }

    void qnm(int index, int & l1, int & l2, int & n, int & m, std::vector<double> & val) const {
      rapidxml::xml_node<> * node = pseudo_node_->first_node("qnm");

      for(int i = 0; i < index; i++) node = node->next_sibling("qnm");

      assert(node != NULL);

      n = value<int>(node->first_attribute("n"));
      m = value<int>(node->first_attribute("m"));
      l1 = value<int>(node->first_attribute("l1"));
      l2 = value<int>(node->first_attribute("l2"));

      int size = value<int>(node->first_attribute("size"));
      val.resize(size);
      std::istringstream stst(node->value());
      for(int ii = 0; ii < size; ii++){
	stst >> val[ii];
      }
      
    }

    void qfcoeff(int index, int ltot, std::vector<double> & val) const {
      rapidxml::xml_node<> * node = pseudo_node_->first_node("qfcoeff");

      while(node){
	int read_index = value<int>(node->first_attribute("i"));
	int read_ltot = value<int>(node->first_attribute("ltot"));
	if(read_index == index && read_ltot == ltot) break;
	node = node->next_sibling("qfcoeff");
      }
      
      assert(node != NULL);

      int size = value<int>(node->first_attribute("size"));
      val.resize(size);
      std::istringstream stst(node->value());
      for(int ii = 0; ii < size; ii++){
	stst >> val[ii];
      }
      
    }
    
  private:

    int nbeta() const {
      return value<int>(pseudo_node_->first_node("nbeta"));
    }


    std::ifstream file_;
    std::vector<char> buffer_;
    rapidxml::xml_document<> doc_;
    rapidxml::xml_node<> * root_node_;
    rapidxml::xml_node<> * pseudo_node_;
    
    
  };

}

#endif
