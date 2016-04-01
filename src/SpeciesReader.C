////////////////////////////////////////////////////////////////////////////////  
// Copyright (c) 2013-2016, Lawrence Livermore National Security, LLC. 
// qb@ll:  Qbox at Lawrence Livermore
//
// This file is part of qb@ll.
//
// Produced at the Lawrence Livermore National Laboratory. 
// Written by Erik Draeger (draeger1@llnl.gov) and Francois Gygi (fgygi@ucdavis.edu)
// and Xavier Andrade (xavier@tddft.org).
//
// Based on the Qbox code by Francois Gygi Copyright (c) 2008 
// LLNL-CODE-635376. All rights reserved. 
//
// qb@ll is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details, in the file COPYING in the
// root directory of this distribution or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
// SpeciesReader.C:
//
////////////////////////////////////////////////////////////////////////////////

#include "Species.h"
#include "SpeciesReader.h"
#include <cassert>
#include <string>
#include <iostream>
#include <vector>
using namespace std;

#if 0
//#if USE_XERCES
#include "StructuredDocumentHandler.h"
#include "SpeciesHandler.h"
#include <xercesc/util/XMLUniDefs.hpp>
#include <xercesc/sax2/Attributes.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>
using namespace xercesc;
#else
#include <sstream>
#include <cstdio>
#include <sys/stat.h>
#endif

class Tag{

public:
  Tag(const string & tag){
    tag_ = tag;
  }

  const string & name() const {
    return tag_;
  }
  
  string start() const {
    return "<" + tag_ + ">";
  }
  
  string end() const {
    return "</" + tag_ + ">";
  }

  string text(const string & buf) const {

    string::size_type pos = 0;

    string::size_type start_pos = buf.find(start(), pos);

    assert(start_pos != string::npos );
    
    start_pos = buf.find(">", start_pos)+1;

    string::size_type end_pos = buf.find(end());

    pos = buf.find(">", end_pos) + 1;

    string::size_type len = end_pos - start_pos;
    
    return buf.substr(start_pos, len);
  }

  template <typename Type>
  void get_value(const string & buf, Type & value) const {
    istringstream stst(text(buf));
    stst >> value;
  }
  
private:
  string tag_;
  
};

////////////////////////////////////////////////////////////////////////////////
SpeciesReader::SpeciesReader(const Context& ctxt) : ctxt_(ctxt) {}

////////////////////////////////////////////////////////////////////////////////
void SpeciesReader::readSpecies (Species& sp, const string uri)
{
  // Simplified parsing not using Xerces
  // parse file uri
  
  if ( ctxt_.oncoutpe() )
  {
    struct stat statbuf;
    bool found_file = !stat(uri.c_str(),&statbuf);
    assert(found_file);
      cout << "  <!-- SpeciesReader opening file "
           << uri << " size: "
           << statbuf.st_size << " -->" << endl;

    FILE* infile;
    infile = fopen(uri.c_str(),"r");
    if ( !infile )
    {
      cout << "  <!-- SpeciesReader::readSpecies could not open file "
           << uri << " for reading" << " -->" << endl;
      return;
    }
    off_t sz = statbuf.st_size;
    string buf;
    buf.resize(sz);
    fread(&buf[0],sizeof(char),sz,infile);

    string::size_type pos = 0;

    string tag, start_tag, end_tag;
    string::size_type start, end, len;

    //ewd:  determine type of pseudopotential
    bool ultrasoft = false;

    {
      Tag tag("ultrasoft_pseudopotential");
      start = buf.find(tag.start(), pos);
      if (start != string::npos) 
	ultrasoft = true;
      sp.usoft_ = ultrasoft;
      if (ultrasoft) 
	cout << "  <!-- SpeciesReader::readSpecies: potential type:  ultrasoft -->" << endl;
      else 
	cout << "  <!-- SpeciesReader::readSpecies: potential type:  norm-conserving -->" << endl;
    }

    {
      Tag tag("description");
      sp.description_ = tag.text(buf);
      cout << "  <!-- SpeciesReader::readSpecies: read " << tag.name() << " "
	   << sp.description_
	   << " -->" << endl;
    }

    {
      Tag tag("symbol");
      tag.get_value(buf, sp.symbol_);
      cout << "  <!-- SpeciesReader::readSpecies: read " << tag.name() << " "
	   << sp.symbol_
	   << " -->" << endl;
    }

    {
      Tag tag("atomic_number");
      tag.get_value(buf, sp.atomic_number_);
      cout << "  <!-- SpeciesReader::readSpecies: read " << tag.name() << " "
	   << sp.atomic_number_
	   << " -->" << endl;
    }

    {
      Tag tag("mass");
      tag.get_value(buf, sp.mass_);
      cout << "  <!-- SpeciesReader::readSpecies: read " << tag.name() << " "
	   << sp.mass_
	   << " -->" << endl;
    }

    {
      Tag tag("valence_charge");
      tag.get_value(buf, sp.zval_);
      cout << "  <!-- SpeciesReader::readSpecies: read " << tag.name() << " "
	   << sp.zval_
	   << " -->" << endl;
    }

    {
      Tag tag("lmax");
      tag.get_value(buf, sp.lmax_);
      cout << "  <!-- SpeciesReader::readSpecies: read " << tag.name() << " "
	   << sp.lmax_
	   << " -->" << endl;
    }
    
    if (!ultrasoft) { 
      tag = "llocal";
      start_tag = string("<") + tag + string(">");
      end_tag = string("</") + tag + string(">");
      start = buf.find(start_tag,pos);
      assert(start != string::npos );
      start = buf.find(">",start)+1;
      end = buf.find(end_tag);
      pos = buf.find(">",end)+1;
      len = end - start;
      {
        istringstream stst(buf.substr(start,len));
        stst >> sp.llocal_;
      }
      cout << "  <!-- SpeciesReader::readSpecies: read " << tag << " "
           << sp.llocal_
           << " -->" << endl;
 
      tag = "nquad";
      start_tag = string("<") + tag + string(">");
      end_tag = string("</") + tag + string(">");
      start = buf.find(start_tag,pos);
      assert(start != string::npos );
      start = buf.find(">",start)+1;
      end = buf.find(end_tag);
      pos = buf.find(">",end)+1;
      len = end - start;
      {
        istringstream stst(buf.substr(start,len));
        stst >> sp.nquad_;
      }
      cout << "  <!-- SpeciesReader::readSpecies: read " << tag << " "
           << sp.nquad_
           << " -->" << endl;
 
      tag = "rquad";
      start_tag = string("<") + tag + string(">");
      end_tag = string("</") + tag + string(">");
      start = buf.find(start_tag,pos);
      assert(start != string::npos );
      start = buf.find(">",start)+1;
      end = buf.find(end_tag);
      pos = buf.find(">",end)+1;
      len = end - start;
      {
        istringstream stst(buf.substr(start,len));
        stst >> sp.rquad_;
      }
      cout << "  <!-- SpeciesReader::readSpecies: read " << tag << " "
           << sp.rquad_
           << " -->" << endl;
    }
    
    tag = "mesh_spacing";
    start_tag = string("<") + tag + string(">");
    end_tag = string("</") + tag + string(">");
    start = buf.find(start_tag,pos);
    assert(start != string::npos );
    start = buf.find(">",start)+1;
    end = buf.find(end_tag);
    pos = buf.find(">",end)+1;
    len = end - start;
    {
      istringstream stst(buf.substr(start,len));
      stst >> sp.deltar_;
    }
    cout << "  <!-- SpeciesReader::readSpecies: read " << tag << " "
         << sp.deltar_
         << " -->" << endl;

    if (ultrasoft) { 
      tag = "nbeta";
      start_tag = string("<") + tag + string(">");
      end_tag = string("</") + tag + string(">");
      start = buf.find(start_tag,pos);
      assert(start != string::npos );
      start = buf.find(">",start)+1;
      end = buf.find(end_tag);
      pos = buf.find(">",end)+1;
      len = end - start;
      {
        istringstream stst(buf.substr(start,len));
        stst >> sp.nbeta_;
      }
      cout << "  <!-- SpeciesReader::readSpecies: read " << tag << " "
           << sp.nbeta_
           << " -->" << endl;
    }

    if (!ultrasoft) { 
      for ( int l = 0; l < sp.lmax_ + 1; l++ )
      {
        // read projector
        int size;
        tag = "projector";
        start_tag = string("<") + tag;
        start = buf.find(start_tag,pos);
        assert(start != string::npos );
 
        pos = buf.find("l=",pos)+3;
        start = pos;
        end = buf.find("\"",start);
        len = end - start;
      
        {
          istringstream stst(buf.substr(start,len));
          int lread;
          stst >> lread;
          //cout << " lread=" << lread << endl;
          //cout << " l=" << l << endl;
          assert(l==lread);
        }
 
        pos = buf.find("size=",pos)+6;
        start = pos;
        end = buf.find("\"",start);
        len = end - start;
        {
          istringstream stst(buf.substr(start,len));
          stst >> size;
          // cout << " size=" << size << endl;
        }
        // read radial potential
        sp.vps_.resize(sp.vps_.size()+1);
        sp.vps_[l].resize(size);
 
        tag = "radial_potential";
        start_tag = string("<") + tag + string(">");
        end_tag = string("</") + tag + string(">");
        start = buf.find(start_tag,pos);
        assert(start != string::npos );
        start = buf.find(">",start)+1;
        end = buf.find(end_tag,start);
        pos = buf.find(">",end)+1;
        len = end - start;
        {
          istringstream stst(buf.substr(start,len));
          for ( int i = 0; i < size; i++ )
          {
            stst >> sp.vps_[l][i];
            //cout << sp.vps_[l][i] << endl;
          }
        }
        cout << "  <!-- SpeciesReader::readSpecies: read " << tag << " l="
             << l << " size=" << size << " -->" << endl;
 
        sp.phi_.resize(sp.phi_.size()+1);
        sp.phi_[l].resize(size);
 
        // read radial function only if the radial_function tag was found, for nonlocal potentials
 
        if ( l != sp.llocal_ )
        {
           tag = "radial_function";
           start_tag = string("<") + tag + string(">");
           end_tag = string("</") + tag + string(">");
           start = buf.find(start_tag,pos);

        /*
        if ( l != sp.llocal_ )
        {
          // if l is not the local potential, there must be a radial function
          assert(start != string::npos );
        }
        */
        
           if ( start != string::npos )
           {
              start = buf.find(">",start)+1;
              end = buf.find(end_tag,start);
              pos = buf.find(">",end)+1;
              len = end - start;
              {
                 istringstream stst(buf.substr(start,len));
                 for ( int i = 0; i < size; i++ )
                 {
                    stst >> sp.phi_[l][i];
                    //cout << sp.phi_[l][i] << endl;
                 }
              }
              cout << "  <!-- SpeciesReader::readSpecies: read " << tag << " l="
                   << l << " size=" << size << " -->" << endl;
           }
        }
        
      }
      
      // read rho_nlcc
      tag = "rho_nlcc";
      start_tag = string("<") + tag;
      start = buf.find(start_tag,pos);
      if (start != string::npos ) {
         pos = buf.find("size=",pos)+6;
         start = pos;
         end = buf.find("\"",start);
         len = end - start;
         int size;
         {
            istringstream stst(buf.substr(start,len));
            stst >> size;
         }

         sp.rhor_nlcc_.resize(size);
         sp.nlcc_ = true;
         cout << "  <!-- SpeciesReader::readSpecies: nlcc found. -->" << endl;
         
         ostringstream oss;
         oss << size;
         start_tag = string("<") + tag + string(" size=\"") + oss.str() + string("\"") + string(">");
         end_tag = string("</") + tag + string(">");
         pos -= 128;  // back position index up a line or two
         start = buf.find(start_tag,pos);
         if ( start != string::npos )
         {
            start = buf.find(">",start)+1;
            end = buf.find(end_tag,start);
            pos = buf.find(">",end)+1;
            len = end - start;
            {
               istringstream stst(buf.substr(start,len));
               for ( int i = 0; i < size; i++ )
               {
                  stst >> sp.rhor_nlcc_[i];
               }
            }
            cout << "  <!-- SpeciesReader::readSpecies: read " << tag << " size=" << size << " -->" << endl;
         }
      }
    }
    else {   // read ultrasoft functions

      // read vlocal
      int size;
      tag = "vlocal";
      start_tag = string("<") + tag;
      start = buf.find(start_tag,pos);
      assert(start != string::npos );
 
      pos = buf.find("size=",pos)+6;
      start = pos;
      end = buf.find("\"",start);
      len = end - start;
      {
        istringstream stst(buf.substr(start,len));
        stst >> size;
        // cout << " size=" << size << endl;
      }

      sp.llocal_ = 0;
      sp.vps_.resize(1);
      sp.vps_[0].resize(size);

      ostringstream oss;
      oss << size;
      start_tag = string("<") + tag + string(" size=\"") + oss.str() + string("\"") + string(">");
      end_tag = string("</") + tag + string(">");
      pos -= 128;  // back position index up a line or two
      start = buf.find(start_tag,pos);
      if ( start != string::npos )
      {
        start = buf.find(">",start)+1;
        end = buf.find(end_tag,start);
        pos = buf.find(">",end)+1;
        len = end - start;
        {
          istringstream stst(buf.substr(start,len));
          for ( int i = 0; i < size; i++ )
          {
            stst >> sp.vps_[0][i];
            //cout << sp.vps_[0][i] << endl;
          }
        }
        cout << "  <!-- SpeciesReader::readSpecies: read " << tag << " size=" << size << " -->" << endl;
      }

      sp.betar_.resize(sp.nbeta_);
      sp.betal_.resize(sp.nbeta_);
      for ( int b = 0; b < sp.nbeta_; b++ )
      {
        // read beta function
        int size;
        tag = "beta";
        start_tag = string("<") + tag;
        start = buf.find(start_tag,pos);
        assert(start != string::npos );
 
        pos = buf.find("l=",pos)+3;
        start = pos;
        end = buf.find("\"",start);
        len = end - start;
      
        {
          istringstream stst(buf.substr(start,len));
          int lread;
          stst >> lread;
          sp.betal_[b] = lread;
        }
 
        pos = buf.find("size=",pos)+6;
        start = pos;
        end = buf.find("\"",start);
        len = end - start;
        {
          istringstream stst(buf.substr(start,len));
          stst >> size;
          // cout << " size=" << size << endl;
        }
        sp.betar_[b].resize(size);

        start_tag = string(">");
        start = buf.find(start_tag,pos) + 1;
        end_tag = string("</") + tag + string(">");
        end = buf.find(end_tag,start);
        pos = buf.find(">",end)+1;
        len = end - start;
        {
          istringstream stst(buf.substr(start,len));
          for ( int i = 0; i < size; i++ )
          {
            stst >> sp.betar_[b][i];
          }
        }
        cout << "  <!-- SpeciesReader::readSpecies: read " << tag << " b="
             << b << " size=" << size << " -->" << endl;
      }

      // read dnm_zero (D_nm^0)
      sp.dzero_.resize(sp.nbeta_);
      for ( int b = 0; b < sp.nbeta_; b++ )
        sp.dzero_[b].resize(sp.nbeta_);

      tag = "dnm_zero";
      start_tag = string("<") + tag + string(">");
      end_tag = string("</") + tag + string(">");
      start = buf.find(start_tag,pos);
      int nbetasq = sp.nbeta_*sp.nbeta_;
      if ( start != string::npos )
      {
        start = buf.find(">",start)+1;
        end = buf.find(end_tag,start);
        pos = buf.find(">",end)+1;
        len = end - start;
        {
          istringstream stst(buf.substr(start,len));
          for ( int i = 0; i < sp.nbeta_; i++ )
            for ( int j = 0; j < sp.nbeta_; j++ )
              stst >> sp.dzero_[i][j];
          cout << "  <!-- SpeciesReader::readSpecies: read " << tag << " -->" << endl;
        }
      }

      // read rinner:  radii inside which Qnm^L are replaced by polynomial expansions
      //   (see Laasonen, Phys. Rev. B 47, 10142 (1993).)
      tag = "rinner";
      start_tag = string("<") + tag;
      start = buf.find(start_tag,pos);
      if (start != string::npos ) {
        pos = buf.find("size=",pos)+6;
        start = pos;
        end = buf.find("\"",start);
        len = end - start;
        {
          istringstream stst(buf.substr(start,len));
          stst >> size;
        }
        assert(size <= 2*sp.lmax_+1);
        sp.rinner_.resize(2*sp.lmax_+1);
        for (int i=0; i<2*sp.lmax_+1; i++)
          sp.rinner_[i] = 0.0;
        start_tag = string(">");
        start = buf.find(start_tag,pos) + 1;
        end_tag = string("</") + tag + string(">");
        end = buf.find(end_tag,start);
        pos = buf.find(">",end)+1;
        len = end - start;
        {
          istringstream stst(buf.substr(start,len));
          for ( int i = 0; i < size; i++ )
          {
            stst >> sp.rinner_[i];
          }
        }
        cout << "  <!-- SpeciesReader::readSpecies: read " << tag
             << " size=" << size << " -->" << endl;
      }
      
      // read Q_nm
      int nqnm = 0;
      for (int i=1; i<=sp.nbeta_; i++)
        nqnm += i;
      sp.nqfun_ = nqnm;
      sp.qfunr_.resize(nqnm);
      sp.qfunl1_.resize(nqnm);
      sp.qfunl2_.resize(nqnm);
      sp.qfunb1_.resize(nqnm);
      sp.qfunb2_.resize(nqnm);
      sp.qfcoeff_.resize(nqnm);
      
      for ( int q = 0; q < nqnm; q++ )
      {
        // read qnm function
        int size;
        int iread, nread, mread, lread1, lread2;
        tag = "qnm";
        start_tag = string("<") + tag;
        start = buf.find(start_tag,pos);
        assert(start != string::npos );
        end_tag = string("</") + tag + string(">");
 
        pos = buf.find("i=",pos)+3;
        start = pos;
        end = buf.find("\"",start);
        len = end - start;
        {
          istringstream stst(buf.substr(start,len));
          stst >> iread;
        }
        pos = buf.find("n=",pos)+3;
        start = pos;
        end = buf.find("\"",start);
        len = end - start;
        {
          istringstream stst(buf.substr(start,len));
          stst >> nread;
        }
        pos = buf.find("m=",pos)+3;
        start = pos;
        end = buf.find("\"",start);
        len = end - start;
        {
          istringstream stst(buf.substr(start,len));
          stst >> mread;
        }
        pos = buf.find("l1=",pos)+4;
        start = pos;
        end = buf.find("\"",start);
        len = end - start;
        {
          istringstream stst(buf.substr(start,len));
          stst >> lread1;
        }
        pos = buf.find("l2=",pos)+4;
        start = pos;
        end = buf.find("\"",start);
        len = end - start;
        {
          istringstream stst(buf.substr(start,len));
          stst >> lread2;
        }
 
        pos = buf.find("size=",pos)+6;
        start = pos;
        end = buf.find("\"",start);
        len = end - start;
        {
          istringstream stst(buf.substr(start,len));
          stst >> size;
        }
        sp.qfunr_[q].resize(size);

        start_tag = string(">");
        start = buf.find(start_tag,pos) + 1;
        end = buf.find(end_tag,start);
        //pos = buf.find(">",end)+1;
        len = end - start;

        {
          istringstream stst(buf.substr(start,len));
          for ( int i = 0; i < size; i++ )
          {
            stst >> sp.qfunr_[q][i];
          }
        }
        sp.qfunl1_[q] = lread1;
        sp.qfunl2_[q] = lread2;
        sp.qfunb1_[q] = nread;
        sp.qfunb2_[q] = mread;
        cout << "  <!-- SpeciesReader::readSpecies: read " << tag << " q="
             << q << " l1=" << lread1 << " l2=" << lread2 << " size=" << size << " -->" << endl;
        pos = buf.find(end_tag,pos);

        // read qcoeffs
        if (sp.rinner_.size() > 0) {
          sp.qfcoeff_[q].resize(2*sp.lmax_+1);
          for (int ltot=0; ltot<2*sp.lmax_+1; ltot++) {
            tag = "qfcoeff";
            start_tag = string("<") + tag;
            start = buf.find(start_tag,pos);
            assert(start != string::npos );
            end_tag = string("</") + tag + string(">");
            
            pos = buf.find("i=",pos)+3;
            start = pos;
            end = buf.find("\"",start);
            len = end - start;
            {
              istringstream stst(buf.substr(start,len));
              stst >> iread;
            }
            pos = buf.find("n=",pos)+3;
            start = pos;
            end = buf.find("\"",start);
            len = end - start;
            {
              istringstream stst(buf.substr(start,len));
              stst >> nread;
            }
            pos = buf.find("m=",pos)+3;
            start = pos;
            end = buf.find("\"",start);
            len = end - start;
            {
              istringstream stst(buf.substr(start,len));
              stst >> mread;
            }
            pos = buf.find("ltot=",pos)+4;
            start = pos;
            end = buf.find("\"",start);
            len = end - start;
            {
              istringstream stst(buf.substr(start,len));
              stst >> lread1;
            }
            pos = buf.find("size=",pos)+6;
            start = pos;
            end = buf.find("\"",start);
            len = end - start;
            {
              istringstream stst(buf.substr(start,len));
              stst >> size;
            }
            sp.qfcoeff_[q][ltot].resize(size);

            start_tag = string(">");
            start = buf.find(start_tag,pos) + 1;
            end = buf.find(end_tag,start);
            //pos = buf.find(">",end)+1;
            len = end - start;
          
            {
              istringstream stst(buf.substr(start,len));
              for ( int i = 0; i < size; i++ )
              {
                stst >> sp.qfcoeff_[q][ltot][i];
              }
            }
            cout << "  <!-- SpeciesReader::readSpecies: read " << tag << " q="
                 << q << " ltot=" << lread1 << " size=" << size << " -->" << endl;
            pos = buf.find(end_tag,pos);
          }
        }
      }

      // read rho_nlcc
      tag = "rho_nlcc";
      start_tag = string("<") + tag;
      start = buf.find(start_tag,pos);
      if (start != string::npos ) {
         pos = buf.find("size=",pos)+6;
         start = pos;
         end = buf.find("\"",start);
         len = end - start;
         {
            istringstream stst(buf.substr(start,len));
            stst >> size;
         }
         sp.rhor_nlcc_.resize(size);
         sp.nlcc_ = true;
         cout << "  <!-- SpeciesReader::readSpecies: nlcc found. -->" << endl;
         
         start_tag = string("<") + tag + string(" size=\"") + oss.str() + string("\"") + string(">");
         end_tag = string("</") + tag + string(">");
         pos -= 128;  // back position index up a line or two
         start = buf.find(start_tag,pos);
         if ( start != string::npos )
         {
            start = buf.find(">",start)+1;
            end = buf.find(end_tag,start);
            pos = buf.find(">",end)+1;
            len = end - start;
            {
               istringstream stst(buf.substr(start,len));
               for ( int i = 0; i < size; i++ )
               {
                  stst >> sp.rhor_nlcc_[i];
               }
            }
            cout << "  <!-- SpeciesReader::readSpecies: read " << tag << " size=" << size << " -->" << endl;
         }
      }
      
    }
    
    sp.uri_ = uri;
  }
}

////////////////////////////////////////////////////////////////////////////////
void SpeciesReader::bcastSpecies(Species& sp)
{
  //cout << ctxt_.mype() << ": starting bcastSpecies" << endl;
  if ( ctxt_.oncoutpe() )
  {
    ctxt_.ibcast_send(1,1,&sp.atomic_number_,1);
    ctxt_.dbcast_send(1,1,&sp.mass_,1);
    ctxt_.ibcast_send(1,1,&sp.zval_,1);
    ctxt_.ibcast_send(1,1,&sp.lmax_,1);
    ctxt_.ibcast_send(1,1,&sp.llocal_,1);
    ctxt_.dbcast_send(1,1,&sp.deltar_,1);
    int iusoft = (sp.ultrasoft() ? 1 : 0 );
    ctxt_.ibcast_send(1,1,&iusoft,1);
    if (!sp.ultrasoft()) {
      ctxt_.ibcast_send(1,1,&sp.nquad_,1);
      ctxt_.dbcast_send(1,1,&sp.rquad_,1);
    }
    else {
      ctxt_.ibcast_send(1,1,&sp.nbeta_,1);
      ctxt_.ibcast_send(1,1,&sp.nqfun_,1);
    }
    int inlcc = (sp.nlcc() ? 1 : 0 );
    ctxt_.ibcast_send(1,1,&inlcc,1);
  }
  else
  {
    // calculate row and col indices of process oncoutpe
    int irow = ctxt_.coutpe();
    while (irow >= ctxt_.nprow())
      irow -= ctxt_.nprow();
    int icol = int (ctxt_.coutpe()/ctxt_.nprow());
    assert(ctxt_.pmap(irow,icol) == ctxt_.coutpe());
            
    ctxt_.ibcast_recv(1,1,&sp.atomic_number_,1,irow,icol);
    ctxt_.dbcast_recv(1,1,&sp.mass_,1,irow,icol);
    ctxt_.ibcast_recv(1,1,&sp.zval_,1,irow,icol);
    ctxt_.ibcast_recv(1,1,&sp.lmax_,1,irow,icol);
    ctxt_.ibcast_recv(1,1,&sp.llocal_,1,irow,icol);
    ctxt_.dbcast_recv(1,1,&sp.deltar_,1,irow,icol);
    int iusoft;
    ctxt_.ibcast_recv(1,1,&iusoft,1,irow,icol);
    sp.usoft_ = ( iusoft == 1 ? true : false );
    if (!sp.ultrasoft()) {
      ctxt_.ibcast_recv(1,1,&sp.nquad_,1,irow,icol);
      ctxt_.dbcast_recv(1,1,&sp.rquad_,1,irow,icol);
      sp.vps_.resize(sp.lmax_+1);
      sp.phi_.resize(sp.lmax_+1);
    }
    else {
      ctxt_.ibcast_recv(1,1,&sp.nbeta_,1,irow,icol);
      ctxt_.ibcast_recv(1,1,&sp.nqfun_,1,irow,icol);
      sp.vps_.resize(1);
      sp.betar_.resize(sp.nbeta_);
      sp.betal_.resize(sp.nbeta_);
      sp.dzero_.resize(sp.nbeta_);
      sp.qfunr_.resize(sp.nqfun_);
      sp.qfunl1_.resize(sp.nqfun_);
      sp.qfunl2_.resize(sp.nqfun_);
      sp.qfunb1_.resize(sp.nqfun_);
      sp.qfunb2_.resize(sp.nqfun_);
      sp.qfcoeff_.resize(sp.nqfun_);
    }
    int inlcc;
    ctxt_.ibcast_recv(1,1,&inlcc,1,irow,icol);
    sp.nlcc_ = ( inlcc == 1 ? true : false );
  }

  ctxt_.string_bcast(sp.symbol_,ctxt_.coutpe());
  ctxt_.string_bcast(sp.description_,ctxt_.coutpe());
  ctxt_.string_bcast(sp.uri_,ctxt_.coutpe());
  
  if (!sp.ultrasoft()) {
  
    for ( int l = 0; l <= sp.lmax_; l++ )
    {
      int np_vps;
      if ( ctxt_.oncoutpe() )
      {
        np_vps = sp.vps_[l].size();
        ctxt_.ibcast_send(1,1,&np_vps,1);
        ctxt_.dbcast_send(np_vps,1,&sp.vps_[l][0],np_vps);
      }
      else
      {
        // calculate row and col indices of process oncoutpe
        int irow = ctxt_.coutpe();
        while (irow >= ctxt_.nprow())
          irow -= ctxt_.nprow();
        int icol = int (ctxt_.coutpe()/ctxt_.nprow());
        assert(ctxt_.pmap(irow,icol) == ctxt_.coutpe());
            
        ctxt_.ibcast_recv(1,1,&np_vps,1,irow,icol);
        sp.vps_[l].resize(np_vps);
        ctxt_.dbcast_recv(np_vps,1,&sp.vps_[l][0],np_vps,irow,icol);
      }
    }

    // broadcast atomic orbitals
    for ( int l = 0; l <= sp.lmax_; l++ )
    {
      int np_phi;
      if ( ctxt_.oncoutpe() )
      {
        np_phi = sp.phi_[l].size();
        ctxt_.ibcast_send(1,1,&np_phi,1);
        ctxt_.dbcast_send(np_phi,1,&sp.phi_[l][0],np_phi);
      }
      else
      {
        // calculate row and col indices of process oncoutpe
        int irow = ctxt_.coutpe();
        while (irow >= ctxt_.nprow())
          irow -= ctxt_.nprow();
        int icol = int (ctxt_.coutpe()/ctxt_.nprow());
        assert(ctxt_.pmap(irow,icol) == ctxt_.coutpe());
      
        ctxt_.ibcast_recv(1,1,&np_phi,1,irow,icol);
        sp.phi_[l].resize(np_phi);
        ctxt_.dbcast_recv(np_phi,1,&sp.phi_[l][0],np_phi,irow,icol);
      }
    }

    // broadcast nlcc density, if any
    int np_nlcc;
    if ( ctxt_.oncoutpe() )
    {
      np_nlcc = sp.rhor_nlcc_.size();
      ctxt_.ibcast_send(1,1,&np_nlcc,1);
      if (np_nlcc > 0)
         ctxt_.dbcast_send(np_nlcc,1,&sp.rhor_nlcc_[0],np_nlcc);
    }
    else
    {
      // calculate row and col indices of process oncoutpe
      int irow = ctxt_.coutpe();
      while (irow >= ctxt_.nprow())
        irow -= ctxt_.nprow();
      int icol = int (ctxt_.coutpe()/ctxt_.nprow());
      assert(ctxt_.pmap(irow,icol) == ctxt_.coutpe());
      
      ctxt_.ibcast_recv(1,1,&np_nlcc,1,irow,icol);
         
      if (np_nlcc > 0) {
         sp.rhor_nlcc_.resize(np_nlcc);
         ctxt_.dbcast_recv(np_nlcc,1,&sp.rhor_nlcc_[0],np_nlcc,irow,icol);         
      }
    }

  }
  else {  // ultrasoft

    // broadcast local potential, stored on coutpe in vps_[0]
    int np_vps;
    if ( ctxt_.oncoutpe() )
    {
      np_vps = sp.vps_[0].size();
      ctxt_.ibcast_send(1,1,&np_vps,1);
      ctxt_.dbcast_send(np_vps,1,&sp.vps_[0][0],np_vps);
    }
    else
    {
      // calculate row and col indices of process oncoutpe
      int irow = ctxt_.coutpe();
      while (irow >= ctxt_.nprow())
        irow -= ctxt_.nprow();
      int icol = int (ctxt_.coutpe()/ctxt_.nprow());
      assert(ctxt_.pmap(irow,icol) == ctxt_.coutpe());
      
      ctxt_.ibcast_recv(1,1,&np_vps,1,irow,icol);
      sp.vps_[0].resize(np_vps);
      ctxt_.dbcast_recv(np_vps,1,&sp.vps_[0][0],np_vps,irow,icol);
    }

      // broadcast beta fns
    for ( int b = 0; b < sp.nbeta_; b++ )
    {
      int np_beta;
      if ( ctxt_.oncoutpe() )
      {
        np_beta = sp.betar_[b].size();
        ctxt_.ibcast_send(1,1,&np_beta,1);
        ctxt_.dbcast_send(np_beta,1,&sp.betar_[b][0],np_beta);
        ctxt_.ibcast_send(1,1,&sp.betal_[b],1);
      }
      else
      {
        // calculate row and col indices of process oncoutpe
        int irow = ctxt_.coutpe();
        while (irow >= ctxt_.nprow())
          irow -= ctxt_.nprow();
        int icol = int (ctxt_.coutpe()/ctxt_.nprow());
        assert(ctxt_.pmap(irow,icol) == ctxt_.coutpe());            
        ctxt_.ibcast_recv(1,1,&np_beta,1,irow,icol);
        sp.betar_[b].resize(np_beta);
        ctxt_.dbcast_recv(np_beta,1,&sp.betar_[b][0],np_beta,irow,icol);
        ctxt_.ibcast_recv(1,1,&sp.betal_[b],1,irow,icol);
      }
    }

    // broadcast dzero
    for ( int b = 0; b < sp.nbeta_; b++ )
    {
      int np_dzero;
      if ( ctxt_.oncoutpe() )
      {
        np_dzero = sp.dzero_[b].size();
        ctxt_.ibcast_send(1,1,&np_dzero,1);
        ctxt_.dbcast_send(np_dzero,1,&sp.dzero_[b][0],np_dzero);
      }
      else
      {
        // calculate row and col indices of process oncoutpe
        int irow = ctxt_.coutpe();
        while (irow >= ctxt_.nprow())
          irow -= ctxt_.nprow();
        int icol = int (ctxt_.coutpe()/ctxt_.nprow());
        assert(ctxt_.pmap(irow,icol) == ctxt_.coutpe());            
        ctxt_.ibcast_recv(1,1,&np_dzero,1,irow,icol);
        sp.dzero_[b].resize(np_dzero);
        ctxt_.dbcast_recv(np_dzero,1,&sp.dzero_[b][0],np_dzero,irow,icol);
      }
    }

    // broadcast rinner
    int np_rinner;
    if ( ctxt_.oncoutpe() )
    {
      np_rinner = sp.rinner_.size();
      ctxt_.ibcast_send(1,1,&np_rinner,1);
      if (np_rinner > 0) 
        ctxt_.dbcast_send(np_rinner,1,&sp.rinner_[0],np_rinner);
    }
    else
    {
      // calculate row and col indices of process oncoutpe
      int irow = ctxt_.coutpe();
      while (irow >= ctxt_.nprow())
        irow -= ctxt_.nprow();
      int icol = int (ctxt_.coutpe()/ctxt_.nprow());
      assert(ctxt_.pmap(irow,icol) == ctxt_.coutpe());            
      ctxt_.ibcast_recv(1,1,&np_rinner,1,irow,icol);
      if (np_rinner > 0) {
        sp.rinner_.resize(np_rinner);
        ctxt_.dbcast_recv(np_rinner,1,&sp.rinner_[0],np_rinner,irow,icol);
      }
    }

    // broadcast qfun, qfcoeff
    for ( int q = 0; q < sp.nqfun_; q++ )
    {
      int np_qfun, np_qfcoeff;
      if ( ctxt_.oncoutpe() )
      {
        np_qfun = sp.qfunr_[q].size();
        ctxt_.ibcast_send(1,1,&np_qfun,1);
        ctxt_.dbcast_send(np_qfun,1,&sp.qfunr_[q][0],np_qfun);
        ctxt_.ibcast_send(1,1,&sp.qfunl1_[q],1);
        ctxt_.ibcast_send(1,1,&sp.qfunl2_[q],1);
        ctxt_.ibcast_send(1,1,&sp.qfunb1_[q],1);
        ctxt_.ibcast_send(1,1,&sp.qfunb2_[q],1);
        int qfcosize = sp.qfcoeff_[q].size();
        ctxt_.ibcast_send(1,1,&qfcosize,1);
        if (qfcosize > 0) {
          for (int ltot=0; ltot<2*sp.lmax_+1; ltot++) {
            np_qfcoeff = sp.qfcoeff_[q][ltot].size();
            ctxt_.ibcast_send(1,1,&np_qfcoeff,1);
            ctxt_.dbcast_send(np_qfcoeff,1,&sp.qfcoeff_[q][ltot][0],np_qfcoeff);
          }
        }
      }
      else
      {
        // calculate row and col indices of process oncoutpe
        int irow = ctxt_.coutpe();
        while (irow >= ctxt_.nprow())
          irow -= ctxt_.nprow();
        int icol = int (ctxt_.coutpe()/ctxt_.nprow());
        assert(ctxt_.pmap(irow,icol) == ctxt_.coutpe());            
        ctxt_.ibcast_recv(1,1,&np_qfun,1,irow,icol);
        sp.qfunr_[q].resize(np_qfun);
        ctxt_.dbcast_recv(np_qfun,1,&sp.qfunr_[q][0],np_qfun,irow,icol);
        ctxt_.ibcast_recv(1,1,&sp.qfunl1_[q],1,irow,icol);
        ctxt_.ibcast_recv(1,1,&sp.qfunl2_[q],1,irow,icol);
        ctxt_.ibcast_recv(1,1,&sp.qfunb1_[q],1,irow,icol);
        ctxt_.ibcast_recv(1,1,&sp.qfunb2_[q],1,irow,icol);
        int qfcosize;
        ctxt_.ibcast_recv(1,1,&qfcosize,1,irow,icol);
        if (qfcosize > 0) { 
          assert(qfcosize == 2*sp.lmax_+1);
          sp.qfcoeff_[q].resize(2*sp.lmax_+1);
          for (int ltot=0; ltot<2*sp.lmax_+1; ltot++) {
            ctxt_.ibcast_recv(1,1,&np_qfcoeff,1,irow,icol);
            sp.qfcoeff_[q][ltot].resize(np_qfcoeff);
            ctxt_.dbcast_recv(np_qfcoeff,1,&sp.qfcoeff_[q][ltot][0],np_qfcoeff,irow,icol);
          }
        }
      }
    }

    int np_nlcc;
    if ( ctxt_.oncoutpe() )
    {
      np_nlcc = sp.rhor_nlcc_.size();
      ctxt_.ibcast_send(1,1,&np_nlcc,1);
      if (np_nlcc > 0)
         ctxt_.dbcast_send(np_nlcc,1,&sp.rhor_nlcc_[0],np_nlcc);
    }
    else
    {
      // calculate row and col indices of process oncoutpe
      int irow = ctxt_.coutpe();
      while (irow >= ctxt_.nprow())
        irow -= ctxt_.nprow();
      int icol = int (ctxt_.coutpe()/ctxt_.nprow());
      assert(ctxt_.pmap(irow,icol) == ctxt_.coutpe());
      
      ctxt_.ibcast_recv(1,1,&np_nlcc,1,irow,icol);
         
      if (np_nlcc > 0) {
         sp.rhor_nlcc_.resize(np_nlcc);
         ctxt_.dbcast_recv(np_nlcc,1,&sp.rhor_nlcc_[0],np_nlcc,irow,icol);         
      }
    }
  }
}
