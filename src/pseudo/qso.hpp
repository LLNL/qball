#include <fstream>
#include <vector>
#include <cassert>

#include <pseudo/base.hpp>
#include <rapidxml.hpp>

namespace pseudopotential {

  template <typename Type>
  static Type value(const rapidxml::xml_node<> * node){
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

      pseudo_node_ = root_node_->first_node("norm_conserving_semilocal_pseudopotential");
      if(pseudo_node_) type_ = type::NORM_CONSERVING_SEMILOCAL;
      
      pseudo_node_ = root_node_->first_node("norm_conserving_pseudopotential");
      if(pseudo_node_) type_ = type::NORM_CONSERVING;
      
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
    
  private:

    ifstream file_;
    vector<char> buffer_;
    rapidxml::xml_document<> doc_;
    rapidxml::xml_node<> * root_node_;
    rapidxml::xml_node<> * pseudo_node_;
    
    
  };

}

