#ifndef PSEUDO_BASE_HPP
#define PSEUDO_BASE_HPP

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

#include <vector>
#include <string>
#include <rapidxml.hpp>

namespace pseudopotential {

  enum class type {
    ULTRASOFT         = 30,
    SEMILOCAL         = 31,
    KLEINMAN_BYLANDER = 32,
    PAW               = 33
  };
  
  enum class status {
    SUCCESS                    = 0,
    FILE_NOT_FOUND             = 455,
    FORMAT_NOT_SUPPORTED       = 456,
    UNKNOWN_FORMAT             = 457,
    UNSUPPORTED_TYPE_ULTRASOFT = 458,
    UNSUPPORTED_TYPE_PAW       = 459,
    UNSUPPORTED_TYPE           = 460
  };

  enum class format {
    UPF1                       = 775,
    UPF2                       = 776,
    QSO                        = 777,
    PSML                       = 778,
    PSF                        = 779,
    CPI                        = 780,
    FHI                        = 781,
    HGH                        = 782
  };

  // these values match libxc convention
  enum class exchange {
    UNKNOWN                    =  -2,
    ANY                        =  -1,
    NONE                       =   0,    
    LDA                        =   1,
    PBE                        = 101,
    PBE_SOL                    = 116,
    B88                        = 106
  };

  enum class correlation {
    UNKNOWN                    =  -2,
    ANY                        =  -1,
    NONE                       =   0,
    LDA_PW                     =  12,
    PBE                        = 130,
    PBE_SOL                    = 133,
    LYP                        = 131
  };
  
  class base {

  public:

    virtual ~base(){}
    virtual pseudopotential::type type() const { return type_; }
    virtual int lmax() const { return lmax_; }

    //Pure virtual functions
    virtual pseudopotential::format format() const = 0;
    virtual int size() const = 0;
    virtual std::string description() const = 0;
    virtual std::string symbol() const = 0;
    virtual int atomic_number() const = 0;
    virtual double mass() const = 0;
    virtual int valence_charge() const = 0;
    virtual int llocal() const = 0;
    virtual int nchannels() const = 0;
    virtual int nquad() const = 0;
    virtual double rquad() const = 0;
    virtual double mesh_spacing() const = 0;
    virtual int mesh_size() const = 0;
    virtual void local_potential(std::vector<double> & potential) const = 0;
    virtual int nprojectors() const = 0;
    virtual void projector(int l, int i, std::vector<double> & proj) const = 0;
    virtual double d_ij(int l, int i, int j) const = 0;
    virtual bool has_radial_function(int l) const= 0;
    virtual void radial_function(int l, std::vector<double> & function) const = 0;
    virtual void radial_potential(int l, std::vector<double> & function) const = 0;
    virtual bool has_nlcc() const= 0;
    virtual void nlcc_density(std::vector<double> & density) const = 0;
    virtual void beta(int index, int & l, std::vector<double> & proj) const = 0;
    virtual void dnm_zero(int nbeta, std::vector<std::vector<double> > & dnm) const = 0;
    virtual bool has_rinner() const = 0;
    virtual void rinner(std::vector<double> & val) const = 0;
    virtual void qnm(int index, int & l1, int & l2, int & n, int & m, std::vector<double> & val) const = 0;
    virtual void qfcoeff(int index, int ltot, std::vector<double> & val) const = 0;

    //Functions for things that might not be provided
    virtual bool has_density() const { return false; }
    virtual void density(std::vector<double> & val) const { val.clear(); }
    virtual int nwavefunctions() const { return 0; }
    virtual void wavefunction(int index, int & n, int & l, double & occ, std::vector<double> & val) const { val.clear(); }
    virtual pseudopotential::exchange exchange() const { return pseudopotential::exchange::UNKNOWN; }
    virtual pseudopotential::correlation correlation() const { return pseudopotential::correlation::UNKNOWN; }
    
  protected:

    template <typename Type>
    static Type value(const rapidxml::xml_base<> * node){
      assert(node);
      std::istringstream stst(node->value());
      Type value;
      stst >> value;
      return value;
    }
    
    pseudopotential::type type_;
    int lmax_;
    
  };

}

#endif
