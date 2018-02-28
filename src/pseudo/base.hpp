
namespace pseudopotential {

  enum class type { ULTRASOFT, NORM_CONSERVING, NORM_CONSERVING_SEMILOCAL };
  
  class base {
    
  protected:

    pseudopotential::type type_;
    
  };

}
