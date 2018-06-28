#ifndef PSEUDO_PSP8_HPP
#define PSEUDO_PSP8_HPP

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

#include "base.hpp"
#include "element.hpp"

namespace pseudopotential {

  class psp8 : public pseudopotential::base {

  public:

    psp8(const std::string & filename){

      filename_ = filename;
      
      std::ifstream original_file(filename.c_str());
      std::string buffer((std::istreambuf_iterator<char>(original_file)), std::istreambuf_iterator<char>());
      std::replace(buffer.begin(), buffer.end(), 'D', 'E');
      std::replace(buffer.begin(), buffer.end(), 'd', 'e');

      std::istringstream file(buffer);
      
      type_ = pseudopotential::type::KLEINMAN_BYLANDER;
      
      // file size
      file.seekg( 0, std::ios::beg);
      std::streampos file_size_ = file.tellg();
      file.seekg( 0, std::ios::end);
      file_size_ = file.tellg() - file_size_;

      //parse the header
      file.seekg( 0, std::ios::beg);
      std::string line;

      //line 1
      getline(file, description_);

      //line 2
      double val;
      file >> val;
      atomic_number_ = round(val);
      file >> val;
      valence_charge_ = round(val);
      getline(file, line);

      //line 3
      int pspcod;
      file >> pspcod >> ixc_ >> lmax_ >> llocal_ >> mesh_size_;
      if(pspcod != 8) throw status::FORMAT_NOT_SUPPORTED;
      getline(file, line);

      //line 4
      file >> val;
      file >> val;
      nlcc_ = (val > 0.0);
      getline(file, line);
	    
      //line 5
      nprojectors_ = 0;
      nchannels_ = 0;
      for(int l = 0; l <= lmax_; l++){
	int np;
	file >> np;
	nprojl_.push_back(np);
	nprojectors_ += np;
	nchannels_ = std::max(nchannels_, np);
      }
      getline(file, line);

      //line 6
      int extension_switch;
      file >> extension_switch;
      getline(file, line);
      has_density_ = extension_switch == 1;
	
      // there is an extra line for spin orbit stuff
      if(extension_switch == 2) getline(file, line);

      if(extension_switch > 2) throw status::FORMAT_NOT_SUPPORTED;
      
      //the projectors and local potential
      projectors_.resize(lmax_ + 1);
      ekb_.resize(lmax_ + 1);
      for(int l = 0; l <= lmax_; l++){
	projectors_[l].resize(nprojl_[l]);
	ekb_[l].resize(nprojl_[l]);

	
	if(l == llocal_){
	  read_local_potential(file);
	  continue;
	}
	
	if(nprojl_[l] == 0) continue;

	
	int read_l;
	file >> read_l;

	assert(read_l == l);

	for(int iproj = 0; iproj < nprojl_[l]; iproj++){
	  projectors_[l][iproj].resize(mesh_size_);
	  file >> ekb_[l][iproj];
	}
	getline(file, line);

	
	for(int ip = 0; ip < mesh_size_; ip++){
	  int read_ip;
	  double grid_point;
	  file >> read_ip >> grid_point;

    	  assert(read_ip == ip + 1);

	  for(int iproj = 0; iproj < nprojl_[l]; iproj++) file >> projectors_[l][iproj][ip];
	  getline(file, line);
	}

      }

      // the local potential if it was not read before
      if(llocal_ > lmax_) read_local_potential(file);


      //NLCC
      if(nlcc_){
	nlcc_density_.resize(mesh_size_);
	
	for(int ip = 0; ip < mesh_size_; ip++){
	  int read_ip;
	  double grid_point;
	  file >> read_ip >> grid_point >> nlcc_density_[ip];
	  assert(read_ip == ip + 1);
	  getline(file, line);
	}
      }


      if(extension_switch == 1){

	density_.resize(mesh_size_);
	
	for(int ip = 0; ip < mesh_size_; ip++){
	  int read_ip;
	  double grid_point;
	  file >> read_ip >> grid_point >> density_[ip];
	  assert(read_ip == ip + 1);
	  getline(file, line);
	}
      }

      
    }

    pseudopotential::format format() const { return pseudopotential::format::PSP8; }
    
    int size() const {
      return file_size_;
    };

    std::string description() const {
      return description_;
    }
    
    std::string symbol() const {
      pseudopotential::element el(atomic_number_);
      return el.symbol();
    }

    int atomic_number() const {
      return atomic_number_;
    }

    double mass() const {
      pseudopotential::element el(atomic_number_);
      return el.mass();
    }
    
    int valence_charge() const {
      return valence_charge_;
    }

    pseudopotential::exchange exchange() const {
      if(ixc_ > 0){
	if(ixc_ == 1) return pseudopotential::exchange::NONE;
	if(ixc_ >= 2 && ixc_ <= 9) return pseudopotential::exchange::LDA;
	if(ixc_ == 11 || ixc_ == 12) return pseudopotential::exchange::PBE;
      }	else {
	return pseudopotential::exchange(-ixc_);
      }
      
      return pseudopotential::exchange::UNKNOWN;
    
    }

    pseudopotential::correlation correlation() const {
      if(ixc_ > 0){
	if(ixc_ == 1) return pseudopotential::correlation::LDA_XC_TETER93;
	if(ixc_ == 2) return pseudopotential::correlation::LDA_PZ;
	if(ixc_ == 7) return pseudopotential::correlation::LDA_PW;
	if(ixc_ == 7 || ixc_ == 8) return pseudopotential::correlation::NONE;
	if(ixc_ == 11) return pseudopotential::correlation::PBE;
	if(ixc_ == 12) return pseudopotential::correlation::NONE;
      }	else {
	return pseudopotential::correlation(-ixc_/1000);
      }
      
      return pseudopotential::correlation::UNKNOWN;

    }

    int llocal() const {
      if(llocal_ > lmax_) return -1;
      return llocal_;
    }

    int nchannels() const {
      return nchannels_;
    }
    
    double mesh_spacing() const {
      return mesh_spacing_;
    }

    int mesh_size() const {
      return mesh_size_;
    }
    
    void local_potential(std::vector<double> & potential) const {
      potential.resize(mesh_size_);
      assert(mesh_size_ == local_potential_.size());
      for(int ip = 0; ip < mesh_size_; ip++) potential[ip] = local_potential_[ip];
    }

    int nprojectors() const {
      return nprojectors_;
    }
    
    void projector(int l, int i, std::vector<double> & proj) const {
      proj.clear();
      
      if(l > lmax_) return;
      if(i >= nprojl_[l]) return;

      proj.resize(mesh_size_);
      assert(mesh_size_ == projectors_[l][i].size());

      for(int ip = 1; ip < mesh_size_; ip++) proj[ip] = projectors_[l][i][ip]/(mesh_spacing()*ip);

      extrapolate_first_point(proj);
      
    }
    
    double d_ij(int l, int i, int j) const {
      if(i != j) return 0.0;
      if(i >= nprojl_[l]) return 0.0;
      return ekb_[l][i];
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
      return nlcc_;
    }

    void nlcc_density(std::vector<double> & density) const {
      density.resize(mesh_size_);
      assert(mesh_size_ == nlcc_density_.size());
      for(int ip = 0; ip < mesh_size_; ip++) density[ip] = nlcc_density_[ip]/(4.0*M_PI);
    }
    
    bool has_density() const{
      return has_density_;
    }

    void density(std::vector<double> & density) const {
      density.resize(mesh_size_);
      assert(mesh_size_ == density_.size());
      for(int ip = 0; ip < mesh_size_; ip++) density[ip] = density_[ip]/(4.0*M_PI);
    }
    
  private:
    
    void extrapolate_first_point(std::vector<double> & function_) const{

      assert(function_.size() >= 4);

      double x1 = mesh_spacing();
      double x2 = 2*mesh_spacing();
      double x3 = 3*mesh_spacing();
      double f1 = function_[1];
      double f2 = function_[2];
      double f3 = function_[3];


      // obtained from:
      // http://www.wolframalpha.com/input/?i=solve+%7Bb*x1%5E2+%2B+c*x1+%2B+d+%3D%3D+f1,++b*x2%5E2+%2B+c*x2+%2B+d+%3D%3D+f2,+b*x3%5E2+%2B+c*x3+%2B+d+%3D%3D+f3+%7D++for+b,+c,+d
      
      function_[0] = f1*x2*x3*(x2 - x3) + f2*x1*x3*(x3 - x1) + f3*x1*x2*(x1 - x2);
      function_[0] /= (x1 - x2)*(x1 - x3)*(x2 - x3);

    }
    
    void read_local_potential(std::istream & file){
      int read_llocal;
      std::string line;
      
      file >> read_llocal;

      assert(llocal_ == read_llocal);
      getline(file, line);

      local_potential_.resize(mesh_size_);     
      for(int ip = 0; ip < mesh_size_; ip++){
	int read_ip;
	double grid_point;
	file>> read_ip >> grid_point >> local_potential_[ip];
	assert(read_ip == ip + 1);
	getline(file, line);

	if(ip == 1) {
	  mesh_spacing_ = grid_point;
	}
      }
      
    }
    
    size_t file_size_;
    std::string description_;
    int atomic_number_;
    int valence_charge_;
    int ixc_;
    int llocal_;
    int mesh_size_;
    int nchannels_;
    double mesh_spacing_;
    std::vector<int> nprojl_;
    int nprojectors_;
    std::vector<std::vector<std::vector<double> > > projectors_;
    std::vector<std::vector<double> > ekb_;
    std::vector<double> local_potential_;
    bool nlcc_;
    std::vector<double> nlcc_density_;
    bool has_density_;
    std::vector<double> density_;
    
  };

}

#endif
