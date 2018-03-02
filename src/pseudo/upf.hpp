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

      //TODO: check for version

      if(root_node_->first_node("PP_HEADER")->first_attribute("pseudo_type")->value() == string("NC")){
	type_ = type::NORM_CONSERVING_SEMILOCAL;
      } else {
	cerr << "Unsupported UPF pseudopotential" << endl;
	exit(1);
      }
      
      assert(root_node_);

      // Read the grid
      
      rapidxml::xml_node<> * node = root_node_->first_node("PP_MESH")->first_node("PP_R");

      assert(node != NULL);

      int size = value<int>(node->first_attribute("size"));
      grid_.resize(size);
      std::istringstream stst(node->value());
      for(int ii = 0; ii < size; ii++){
	stst >> grid_[ii];
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
      return number_of_projectors()/(lmax() + 1);
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
      for(int ii = 0; ii < size; ii++) stst >> potential[ii];

      interpolate(potential);
      
    }

    int number_of_projectors() const {
      return value<int>(root_node_->first_node("PP_HEADER")->first_attribute("number_of_proj"));
    }
    
    bool has_projectors(int l) const {
      return l >=0 && l <= lmax();
    }
    
    void projector(int l, int i, std::vector<double> & proj) const {
      rapidxml::xml_node<> * node;

      for(int iproj = 1; iproj <= number_of_projectors(); iproj++){
	std::string tag = "PP_BETA." + std::to_string(iproj);
	node = root_node_->first_node("PP_NONLOCAL")->first_node(tag.c_str());

	assert(node);
	
	int read_l = value<int>(node->first_attribute("angular_momentum"));
	int read_i = (value<int>(node->first_attribute("index")) - 1)%nchannels();
	if(l == read_l && i == read_i) break;
      }

      assert(node != NULL);

      int size = value<int>(node->first_attribute("size"));
      proj.resize(size);
      std::istringstream stst(node->value());
      for(int ii = 0; ii < size; ii++) stst >> proj[ii];

      interpolate(proj);

    }
    
    double d_ij(int l, int i, int j) const {
#if 0
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
#endif 
    }

    bool has_radial_function(int l) const{
#if 0
      rapidxml::xml_node<> * node = pseudo_node_->first_node("projector");
      
      while(node){
	int read_l = value<int>(node->first_attribute("l"));
	if(l == read_l) break;
	node = node->next_sibling("projector");
      }
      
      return node->first_node("radial_function") != NULL;
#endif 
    }

    void radial_function(int l, std::vector<double> & function) const {
#if 0
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
#endif 
    }

    void radial_potential(int l, std::vector<double> & function) const {
#if 0
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
#endif 
    }

    bool has_nlcc() const{
#if 0
      return pseudo_node_->first_node("rho_nlcc");
#endif
    }

    void nlcc_density(std::vector<double> & density) const {
#if 0
      rapidxml::xml_node<> * node = pseudo_node_->first_node("rho_nlcc");
      assert(node);
      int size = value<int>(node->first_attribute("size"));
      density.resize(size);
      std::istringstream stst(node->value());
      for(int ii = 0; ii < size; ii++){
	stst >> density[ii];
      }
#endif
    }
    
    void beta(int index, int & l, std::vector<double> & proj) const {
#if 0
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
#endif 
    }

    void dnm_zero(int nbeta, std::vector<std::vector<double> > & dnm) const {
#if 0
      rapidxml::xml_node<> * node = pseudo_node_->first_node("dnm_zero");
      std::istringstream stst(node->value());

      dnm.resize(nbeta);
      for(int i = 0; i < nbeta; i++){
	dnm[i].resize(nbeta);
	for ( int j = 0; j < nbeta; j++){
	  stst >> dnm[i][j];
	}
      }
#endif
    }

    bool has_rinner() const {
#if 0
      rapidxml::xml_node<> * node = pseudo_node_->first_node("rinner");
      return node;
#endif
    }
    
    void rinner(std::vector<double> & val) const {
#if 0
      rapidxml::xml_node<> * node = pseudo_node_->first_node("rinner");

      assert(node != NULL);

      int size = value<int>(node->first_attribute("size"));
      val.resize(size);
      std::istringstream stst(node->value());
      for(int ii = 0; ii < size; ii++){
	stst >> val[ii];
      }
#endif
    }

    void qnm(int index, int & l1, int & l2, int & n, int & m, std::vector<double> & val) const {
#if 0
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
#endif 
    }

    void qfcoeff(int index, int ltot, std::vector<double> & val) const {
#if 0
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
#endif 
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
    
    ifstream file_;
    vector<char> buffer_;
    rapidxml::xml_document<> doc_;
    rapidxml::xml_node<> * root_node_;
    std::vector<double> grid_;
    
    
  };

}

#endif
