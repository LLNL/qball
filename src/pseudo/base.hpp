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
#include <sstream>

#define MAX_L 10
#define INVALID_L 333

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
    FILE_NOT_FOUND             = 773,
    UNKNOWN                    = 774,
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
    virtual double mesh_spacing() const = 0;
    virtual int mesh_size() const = 0;
    virtual void local_potential(std::vector<double> & potential) const = 0;
    virtual int nprojectors() const = 0;
    virtual void projector(int l, int i, std::vector<double> & proj) const = 0;
    virtual double d_ij(int l, int i, int j) const = 0;
    virtual bool has_radial_function(int l) const= 0;
    virtual void radial_function(int l, std::vector<double> & function) const = 0;
    virtual void radial_potential(int l, std::vector<double> & function) const = 0;

    virtual void grid(std::vector<double> & val) const {
      val.resize(mesh_size());
      for(unsigned ii = 0; ii < val.size(); ii++) val[ii] = ii*mesh_spacing();
    }

    virtual void grid_weights(std::vector<double> & val) const {
      val.resize(mesh_size());
      for(unsigned ii = 1; ii < val.size() - 1; ii++) val[ii] = mesh_spacing();
      val[0] = 0.5*mesh_spacing();
      val[val.size() - 1] = 0.5*mesh_spacing();
    }

    //Functions for things that might not be provided
    virtual int nquad() const { return 0; }
    virtual double rquad() const { return 0.0; }
    virtual bool has_nlcc() const { return false; }
    virtual void nlcc_density(std::vector<double> & density) const { density.clear(); }
    virtual void beta(int index, int & l, std::vector<double> & proj) const { l = 0; proj.clear(); }
    virtual void dnm_zero(int nbeta, std::vector<std::vector<double> > & dnm) const { dnm.clear(); }
    virtual bool has_rinner() const { return false; }
    virtual void rinner(std::vector<double> & val) const { val.clear(); }
    virtual void qnm(int index, int & l1, int & l2, int & n, int & m, std::vector<double> & val) const { val.clear(); }
    virtual void qfcoeff(int index, int ltot, std::vector<double> & val) const { val.clear(); }
    virtual bool has_density() const { return false; }
    virtual void density(std::vector<double> & val) const { val.clear(); }
    virtual int nwavefunctions() const { return 0; }
    virtual void wavefunction(int index, int & n, int & l, double & occ, std::vector<double> & val) const { val.clear(); }
    virtual pseudopotential::exchange exchange() const { return pseudopotential::exchange::UNKNOWN; }
    virtual pseudopotential::correlation correlation() const { return pseudopotential::correlation::UNKNOWN; }

    virtual bool has_total_angular_momentum() const { return false; }
    virtual int projector_2j(int l, int ic) const { return 0; } // returns j multiplied by 2
    virtual int wavefunction_2j(int ii) const { return 0; } // returns j multiplied by 2
    
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
