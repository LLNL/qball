#ifndef PSEUDO_BASE_HPP
#define PSEUDO_BASE_HPP

#include <vector>

namespace pseudopotential {

  enum class type { ULTRASOFT, NORM_CONSERVING, NORM_CONSERVING_SEMILOCAL };
  
  class base {

  public:

    pseudopotential::type type() const { return type_; }
    
  protected:

    pseudopotential::type type_;
    
  };

}

#endif
