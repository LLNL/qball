#ifndef PSEUDO_UPF_HPP
#define PSEUDO_UPF_HPP

#include <fstream>
#include <vector>
#include <cassert>

#include <pseudo/base.hpp>
#include <rapidxml.hpp>

#include <pseudo/chemical_element.hpp>
#include <pseudo/spline.h>

namespace pseudopotential {

  class upf : public pseudopotential::base {

  public:
    
    upf(const std::string & filename):
      file_(filename),
      buffer_((istreambuf_iterator<char>(file_)), istreambuf_iterator<char>()){
      
      buffer_.push_back('\0');
      doc_.parse<0>(&buffer_[0]);
      
      root_node_ = doc_.first_node("UPF");
      
      if(!root_node_){
	cerr << "Error: File '" << filename << "' is not a UPF 2 file (version 1 is not supported)." << endl;
	exit(1);
      }
      
      if(root_node_->first_attribute("version")->value()[0] != '2'){
	cerr << "Unsupported UPF pseudopotential, can only read UPF v2." << endl;
	exit(1);
      }
      
      std::string pseudo_type = root_node_->first_node("PP_HEADER")->first_attribute("pseudo_type")->value();
      
      if(pseudo_type == "NC"){
	type_ = type::NORM_CONSERVING_SEMILOCAL;
      } else if(pseudo_type == "USPP"){
	type_ = type::ULTRASOFT;
      } else {
	cerr << "Error: Unsupported UPF pseudopotential." << endl;
	exit(1);
      }
      
      assert(root_node_);

      // Read the grid
      {
	rapidxml::xml_node<> * node = root_node_->first_node("PP_MESH")->first_node("PP_R");
	
	assert(node);
	
	int size = value<int>(node->first_attribute("size"));
	grid_.resize(size);
	std::istringstream stst(node->value());
	for(int ii = 0; ii < size; ii++){
	  stst >> grid_[ii];
	}

	if(fabs(grid_[0]) > 1e-10){
	  cerr << "Unsupported UPF pseudopotential, grid does not start at zero." << endl;
	  exit(1);
	}
	
      }
      
      //Read dij once
      {
      	rapidxml::xml_node<> * node = root_node_->first_node("PP_NONLOCAL")->first_node("PP_DIJ");

	assert(node);
	
	dij_.resize(nprojectors()*nprojectors());
	
	std::istringstream stst(node->value());
	for(unsigned ii = 0; ii < dij_.size(); ii++){
	  stst >> dij_[ii];
	  dij_[ii] *= 0.5; //convert from Rydberg to Hartree
	}
      }
	
    }

    int size() const { return buffer_.size(); };

    std::string description() const {
      return root_node_->first_node("PP_INFO")->value();
    }
    
    std::string symbol() const {
      return root_node_->first_node("PP_HEADER")->first_attribute("element")->value();
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
      return value<int>(root_node_->first_node("PP_HEADER")->first_attribute("z_valence"));
    }

    int lmax() const {
      return value<int>(root_node_->first_node("PP_HEADER")->first_attribute("l_max"));
    }

    int llocal() const {
      return value<int>(root_node_->first_node("PP_HEADER")->first_attribute("l_local"));
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

    int nchannels() const {
      return nprojectors()/(lmax() + 1);
    }
    
    int nbeta() const {
      return 0;
    }

    void local_potential(std::vector<double> & potential) const {
      rapidxml::xml_node<> * node = root_node_->first_node("PP_LOCAL");

      assert(node);

      int size = value<int>(node->first_attribute("size"));

      potential.resize(size);
      std::istringstream stst(node->value());
      for(int ii = 0; ii < size; ii++) {
	stst >> potential[ii];
	potential[ii] *= 0.5; //Convert from Rydberg to Hartree
      }
      
      interpolate(potential);
      
    }

    int nprojectors() const {
      return value<int>(root_node_->first_node("PP_HEADER")->first_attribute("number_of_proj"));
    }
    
    bool has_projectors(int l) const {
      return l >=0 && l <= lmax();
    }
    
    void projector(int l, int i, std::vector<double> & proj) const {
      rapidxml::xml_node<> * node = NULL;

      for(int iproj = 1; iproj <= nprojectors(); iproj++){
	std::string tag = "PP_BETA." + std::to_string(iproj);
	node = root_node_->first_node("PP_NONLOCAL")->first_node(tag.c_str());

	assert(node);
	
	int read_l = value<int>(node->first_attribute("angular_momentum"));
	int read_i = (value<int>(node->first_attribute("index")) - 1)%nchannels();
	if(l == read_l && i == read_i) break;
      }

      assert(node);

      int size = value<int>(node->first_attribute("size"));
      proj.resize(size);
      std::istringstream stst(node->value());
      for(int ii = 0; ii < size; ii++) stst >> proj[ii];

      //the projectors come multiplied by r, so we have to divide and fix the first point
      for(int ii = 1; ii < size; ii++) proj[ii] /= grid_[ii];
      extrapolate_first_point(proj);
      
      interpolate(proj);
    }
    
    double d_ij(int l, int i, int j) const {
      int n = l*nchannels() + i;
      int m = l*nchannels() + j;

      return dij_[n*nprojectors() + m];
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
      return root_node_->first_node("PP_NLCC");
    }

    void nlcc_density(std::vector<double> & density) const {
      rapidxml::xml_node<> * node = root_node_->first_node("PP_NLCC");
      assert(node);
      
      int size = value<int>(node->first_attribute("size"));
      density.resize(size);

      std::istringstream stst(node->value());
      for(int ii = 0; ii < size; ii++) stst >> density[ii];

      interpolate(density);
    }
    
    void beta(int index, int & l, std::vector<double> & proj) const {
      proj.clear();
    }

    void dnm_zero(int nbeta, std::vector<std::vector<double> > & dnm) const {
      dnm.clear();
    }

    bool has_rinner() const {
      return false;
    }
    
    void rinner(std::vector<double> & val) const {
      val.clear();
    }

    void qnm(int index, int & l1, int & l2, int & n, int & m, std::vector<double> & val) const {
      val.clear();
    }

    void qfcoeff(int index, int ltot, std::vector<double> & val) const {
      val.clear();
    }
    
  private:

    void interpolate(std::vector<double> & function) const {
      std::vector<double> function_in_grid = function;
      
      Spline function_spline;
      function_spline.fit(grid_.data(), function_in_grid.data(), function_in_grid.size(), SPLINE_FLAT_BC, SPLINE_NATURAL_BC);

      function.clear();
      for(double rr = 0.0; rr <= grid_[grid_.size() - 1]; rr += mesh_spacing()){
	function.push_back(function_spline.value(rr));
      }
    }

    void extrapolate_first_point(std::vector<double> & function_) const{
      double x1 = grid_[1];
      double x2 = grid_[2];
      double x3 = grid_[3];
      double f1 = function_[1];
      double f2 = function_[2];
      double f3 = function_[3];

      function_[0] = f1*x2*x3*(x2 - x3) + f2*x1*x3*(x3 - x1) + f3*x1*x2*(x1 - x2);
      function_[0] /= (x1 - x2)*(x1 - x3)*(x2 - x3);

    }
    
    ifstream file_;
    vector<char> buffer_;
    rapidxml::xml_document<> doc_;
    rapidxml::xml_node<> * root_node_;
    std::vector<double> grid_;
    std::vector<double> dij_;
    
  };

}

#endif
