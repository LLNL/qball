#ifndef PSEUDO_UPF1_HPP
#define PSEUDO_UPF1_HPP

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
#include <cmath>

#include "upf.hpp"
#include "base.hpp"
#include <rapidxml.hpp>

#include "chemical_element.hpp"

namespace pseudopotential {

  class upf1 : public pseudopotential::upf {

  public:
    
    upf1(const std::string & filename):
      file_(filename),
      buffer_((std::istreambuf_iterator<char>(file_)), std::istreambuf_iterator<char>()){

      buffer_.push_back('\0');
      doc_.parse<0>(&buffer_[0]);
      
      std::istringstream header(doc_.first_node("PP_HEADER")->value());

      std::string line;

      int version_number;
      header >> version_number;
      getline(header, line);
      
      header >> symbol_;
      getline(header, line);
      
      std::string pseudo_type;
      header >> pseudo_type;
      getline(header, line);

      //nlcc tag
      getline(header, line);

      getline(header, xc_functional_);

      header >> zval_;
      getline(header, line);

      //total energy
      getline(header, line);

      //cutoff
      getline(header, line);

      //skip lmax
      getline(header, line);

      int size;
      header >> size;
      getline(header, line);

      header >> nwavefunctions_;
      header >> nprojectors_;
      getline(header, line);

      std::transform(pseudo_type.begin(), pseudo_type.end(), pseudo_type.begin(), ::tolower);

      if(pseudo_type == "nc" || pseudo_type == "sl"){
	type_ = pseudopotential::type::KLEINMAN_BYLANDER;
      } else if(pseudo_type == "uspp"){
	throw status::UNSUPPORTED_TYPE_ULTRASOFT;
      } else if(pseudo_type == "paw") {
	throw status::UNSUPPORTED_TYPE_PAW;
      } else {
	throw status::UNSUPPORTED_TYPE;
      }

      // Read the grid
      {
	rapidxml::xml_node<> * node = doc_.first_node("PP_MESH")->first_node("PP_R");
	assert(node);

	std::istringstream stst(node->value());

	//check whether the first point is zero or not
	double xmin;
	stst >> xmin;

	start_point_ = 0;
	if(xmin > 1.0e-10) start_point_ = 1;

	grid_.resize(size + start_point_);

	grid_[0] = 0.0;
	grid_[start_point_] = xmin;
	for(int ii = 0; ii < size - 1; ii++) stst >> grid_[1 + start_point_ + ii];
	
	assert(fabs(grid_[0]) <= 1e-10);

	mesh_size_ = 0;
	for(double rr = 0.0; rr <= grid_[grid_.size() - 1]; rr += mesh_spacing()) mesh_size_++;

      }

          
      //lmax and lloc
      {

	proj_l_.resize(nprojectors());
	proj_c_.resize(nprojectors());

	rapidxml::xml_node<> * node = doc_.first_node("PP_NONLOCAL")->first_node("PP_BETA");

	std::vector<bool> has_l(MAX_L, false);

	lmax_ = 0;
	int iproj = 0;
	while(node){
	  
	  std::string line;
	  std::istringstream stst(node->value());

	  int read_i, read_l;
	  
	  stst >> read_i >> read_l;

	  read_i--;

	  assert(iproj == read_i);
	  
	  lmax_ = std::max(lmax_, read_l);
	  has_l[read_l] = true;
	  proj_l_[iproj] = read_l;
	  proj_c_[iproj] = 0;
	  for(int jproj = 0; jproj < iproj; jproj++) if(read_l == proj_l_[jproj]) proj_c_[iproj]++;
	  
	  node = node->next_sibling("PP_BETA");
	  iproj++;
	}

	assert(lmax_ >= 0);

	llocal_ = -1;
	for(int l = 0; l <= lmax_; l++) if(!has_l[l]) llocal_ = l;
	
      }

      //Read dij
      {
      	rapidxml::xml_node<> * node = doc_.first_node("PP_NONLOCAL")->first_node("PP_DIJ");

	assert(node);
	
	dij_.resize(nchannels()*nchannels()*(lmax_ + 1));
	
	for(unsigned kk = 0; kk < dij_.size(); kk++) dij_[kk] = 0.0;
	
	std::istringstream stst(node->value());

	int nnonzero;

	stst >> nnonzero;
	getline(stst, line);
	
	for(int kk = 0; kk < nnonzero; kk++){
	  int ii, jj;
	  double val;
	  stst >> ii >> jj >> val;
	  val *= 2.0; //convert from 1/Rydberg to 1/Hartree
	  ii--;
	  jj--;

	  assert(proj_l_[ii] == proj_l_[jj]);
	  
	  d_ij(proj_l_[ii], proj_c_[ii], proj_c_[jj]) = val;
	}
      }

    }

    pseudopotential::format format() const { return pseudopotential::format::UPF2; }
    
    int size() const { return buffer_.size(); };

    std::string description() const {
      return doc_.first_node("PP_INFO")->value();
    }
    
    std::string symbol() const {
      return symbol_;
    }

    int atomic_number() const {
      chemical_element el(symbol());
      return el.atomic_number();
    }

    double mass() const {
      chemical_element el(symbol());
      return el.mass();
    }
    
    int valence_charge() const {
      return zval_;
    }

    int llocal() const {
      return llocal_;
    }

    pseudopotential::exchange exchange() const {
      if(xc_functional_ == "PBE") return pseudopotential::exchange::PBE;
      if(xc_functional_ == "PBESOL") return pseudopotential::exchange::PBE_SOL;
      if(xc_functional_ == "SLA  PW   NOGX NOGC") return pseudopotential::exchange::LDA;
      if(xc_functional_ == "BLYP") return pseudopotential::exchange::B88;
      return pseudopotential::exchange::UNKNOWN;
    }

    pseudopotential::correlation correlation() const {
      if(xc_functional_ == "PBE") return pseudopotential::correlation::PBE;
      if(xc_functional_ == "PBESOL") return pseudopotential::correlation::PBE_SOL;
      if(xc_functional_ == "SLA  PW   NOGX NOGC") return pseudopotential::correlation::LDA_PW;
      if(xc_functional_ == "BLYP") return pseudopotential::correlation::LYP;
      return pseudopotential::correlation::UNKNOWN;
    }

    int nchannels() const {
      if(llocal() >= 0){
	return nprojectors()/lmax();
      } else {
	return nprojectors()/(lmax() + 1);
      }
    }
    
    void local_potential(std::vector<double> & potential) const {
      rapidxml::xml_node<> * node = doc_.first_node("PP_LOCAL");

      assert(node);

      potential.resize(grid_.size());
      std::istringstream stst(node->value());
      for(unsigned ii = 0; ii < grid_.size() - start_point_; ii++) {
	stst >> potential[ii + start_point_];
	potential[ii + start_point_] *= 0.5; //Convert from Rydberg to Hartree
      }
      if(start_point_ > 0) extrapolate_first_point(potential);

      interpolate(potential);
      
    }

    int nprojectors() const {
      return nprojectors_;
    }
    
    void projector(int l, int i, std::vector<double> & proj) const {
      rapidxml::xml_node<> * node = doc_.first_node("PP_NONLOCAL")->first_node("PP_BETA");

      assert(node);

      int iproj = 0;
      while(node){
	
	if(l != proj_l_[iproj] || i != proj_c_[iproj]) {
	  iproj++;
	  node = node->next_sibling("PP_BETA");
	  continue;
	}
	
	std::string line;
	std::istringstream stst(node->value());
	
	int read_i, read_l, size;

	stst >> read_i >> read_l;
	getline(stst, line);

	assert(read_l == proj_l_[iproj]);
	
	stst >> size;
	getline(stst, line);

	assert(size >= 0);
	assert(size <= int(grid_.size()));

	proj.resize(grid_.size());

	for(int ii = 0; ii < size; ii++) stst >> proj[ii + start_point_];
	for(unsigned ii = size; ii < grid_.size() - start_point_; ii++) proj[ii + start_point_] = 0.0; 
	    
	break;
      }

      //the projectors come in Rydberg and multiplied by r, so we have to divide and fix the first point
      for(unsigned ii = 1; ii < proj.size(); ii++) proj[ii] /= 2.0*grid_[ii];
      extrapolate_first_point(proj);
      
      interpolate(proj);
    }

    bool has_radial_function(int l) const{
      return false;
    }

    void radial_function(int l, std::vector<double> & function) const {
      function.clear();
    }

    void radial_potential(int l, std::vector<double> & function) const {
      function.clear();
    }

    bool has_nlcc() const{
      return doc_.first_node("PP_NLCC");
    }

    void nlcc_density(std::vector<double> & density) const {
      rapidxml::xml_node<> * node = doc_.first_node("PP_NLCC");
      assert(node);
      std::istringstream stst(node->value());
      
      density.resize(grid_.size());

      for(unsigned ii = 0; ii < grid_.size() - start_point_; ii++) stst >> density[start_point_ + ii];
      extrapolate_first_point(density);
      // this charge does not come multiplied by anything
      
      interpolate(density);
    }
    
    bool has_density() const {
      return doc_.first_node("PP_RHOATOM");
    }
      
    void density(std::vector<double> & val) const {
      rapidxml::xml_node<> * node = doc_.first_node("PP_RHOATOM");
      assert(node);
      
      val.resize(grid_.size());
      
      std::istringstream stst(node->value());
      for(unsigned ii = 0; ii < grid_.size() - start_point_; ii++) stst >> val[start_point_ + ii];

      // the density comes multiplied by 4\pi r
      for(unsigned ii = 1; ii < val.size(); ii++) val[ii] /= 4.0*M_PI*grid_[ii]*grid_[ii];
      extrapolate_first_point(val);
      
      interpolate(val);
    }

    int nwavefunctions() const {
      return nwavefunctions_;
    }
    
    void wavefunction(int index, int & n, int & l, double & occ, std::vector<double> & proj) const {
      rapidxml::xml_node<> * node = doc_.first_node("PP_PSWFC");
      
      assert(node);

      std::istringstream stst(node->value());

      std::string line;
      
      //skip until the correct wavefunction
      for(int ii = 0; ii < index; ii++){
	double tmp;

	stst >> line;
	getline(stst, line);
	for(unsigned ii = 0; ii < grid_.size() - start_point_; ii++) stst >> tmp;
      }

      std::string label;
      stst >> label >> l >> occ;
      getline(stst, line);
      
      n = std::stoi(label.substr(0, 1));
      
      proj.resize(grid_.size());

      for(unsigned ii = 0; ii < grid_.size() - start_point_; ii++) stst >> proj[ii + start_point_];

      //the wavefunctions come multiplied by r, so we have to divide and fix the first point
      for(unsigned ii = 1; ii < grid_.size() - start_point_; ii++) proj[ii] /= grid_[ii];
      extrapolate_first_point(proj);
      
      interpolate(proj);
    }
    
  private:

    std::ifstream file_;
    std::vector<char> buffer_;
    rapidxml::xml_document<> doc_;
    int start_point_;

    std::string symbol_;
    std::string xc_functional_;
    int zval_;
    int nwavefunctions_;
    int nprojectors_;
    std::vector<int> proj_l_;
    std::vector<int> proj_c_;
    
  };

}

#endif
