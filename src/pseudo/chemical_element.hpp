#ifndef PSEUDO_CHEMICALELEMENT_HPP
#define PSEUDO_CHEMICALELEMENT_HPP

#include <string>

namespace pseudopotential {
  
  class chemical_element {

  public:

    chemical_element(const std::string & symbol = "none"){
      this->set(symbol);
    }

    std::string symbol() const;

    double charge() const { return -1.0*z; }
    double mass() const;
    int atomic_number() const { return z; }
    
  private:
  
    char z;
    void set(const std::string & symbol);
  
  };

}

#endif

// Local Variables:
// mode: c++
// coding: utf-8
// End:
