#include <fstream>
#include <vector>
#include <cassert>

#include <pseudo/base.hpp>
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
      buffer_((istreambuf_iterator<char>(file_)), istreambuf_iterator<char>()){

      buffer_.push_back('\0');
      doc_.parse<0>(&buffer_[0]);

      root_node_ = doc_.first_node("fpmd:species");

      pseudo_node_ = root_node_->first_node("ultrasoft_pseudopotential");
      if(pseudo_node_) type_ = type::ULTRASOFT;

      if(!pseudo_node_){
	pseudo_node_ = root_node_->first_node("norm_conserving_semilocal_pseudopotential");
	if(pseudo_node_) type_ = type::NORM_CONSERVING_SEMILOCAL;
      }

      if(!pseudo_node_){
	pseudo_node_ = root_node_->first_node("norm_conserving_pseudopotential");
	if(pseudo_node_) type_ = type::NORM_CONSERVING;
      }

      assert(pseudo_node_);
      
    }

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

    int lmax() const {
      return value<int>(pseudo_node_->first_node("lmax"));
    }

    int llocal() const {
      return value<int>(pseudo_node_->first_node("llocal"));
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

    int nbeta() const {
      return value<int>(pseudo_node_->first_node("nbeta"));
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

    bool has_projectors(int l){
      rapidxml::xml_node<> * node = pseudo_node_->first_node("projector");
      while(node){
	int read_l = value<int>(node->first_attribute("l"));
	if(l == read_l) break;
	node = node->next_sibling();
      }
      return node != NULL;
    }
    
    void projector(int l, int i, std::vector<double> & proj) const {
      rapidxml::xml_node<> * node = pseudo_node_->first_node("projector");

      while(node){
	int read_l = value<int>(node->first_attribute("l"));
	int read_i = value<int>(node->first_attribute("i")) - 1;
	if(l == read_l && i == read_i) break;
	node = node->next_sibling();
      }

      assert(node != NULL);

      int size = value<int>(node->first_attribute("size"));
      proj.resize(size);
      std::istringstream stst(node->value());
      for(int ii = 0; ii < size; ii++){
	stst >> proj[ii];
      }
      
    }
    
    double d_ij(int l, int i, int j) const {
      rapidxml::xml_node<> * node = pseudo_node_->first_node("d_ij");
      
      while(node){
	int read_l = value<int>(node->first_attribute("l"));
	int read_i = value<int>(node->first_attribute("i")) - 1;
	int read_j = value<int>(node->first_attribute("j")) - 1;
	if(l == read_l && i == read_i && j == read_j) break;
	node = node->next_sibling();
      }

      assert(node != NULL);

      return value<double>(node);
      
    }

    bool has_radial_function(int l){
      rapidxml::xml_node<> * node = pseudo_node_->first_node("projector");
      
      while(node){
	int read_l = value<int>(node->first_attribute("l"));
	if(l == read_l) break;
	node = node->next_sibling();
      }
      
      return node->first_node("radial_function") != NULL;
      
    }

    void radial_function(int l, std::vector<double> & function) const {
      rapidxml::xml_node<> * node = pseudo_node_->first_node("projector");

      while(node){
	int read_l = value<int>(node->first_attribute("l"));
	if(l == read_l) break;
	node = node->next_sibling();
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
	node = node->next_sibling();
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

    bool has_nlcc(){
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
    

  private:

    ifstream file_;
    vector<char> buffer_;
    rapidxml::xml_document<> doc_;
    rapidxml::xml_node<> * root_node_;
    rapidxml::xml_node<> * pseudo_node_;
    
    
  };

}

