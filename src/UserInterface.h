////////////////////////////////////////////////////////////////////////////////  
// Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
// qb@ll:  Qbox at Lawrence Livermore
//
// This file is part of qb@ll.
//
// Produced at the Lawrence Livermore National Laboratory. 
// Written by Erik Draeger (draeger1@llnl.gov) and Francois Gygi (fgygi@ucdavis.edu).
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
// UserInterface.h:
//
////////////////////////////////////////////////////////////////////////////////
// $ Id: $

#include <config.h>

#ifndef USER_INTERFACE_H
#define USER_INTERFACE_H

#include <iostream>
#include <string>
#include <cstring>
#include <iomanip>
#include <list>
#include <algorithm>
#include "Context.h"
#include "Dimensions.h"

using namespace std;

class UserInterface;

class Cmd
{
  public:
  UserInterface *ui;
  virtual char const*name(void) const = 0;
  virtual char const*help_msg(void) const = 0;
  virtual int action(int argc, char **argv) = 0;
};

class Var {
public:
  UserInterface *ui;
  virtual char const*name ( void ) const = 0;
  virtual int set ( int argc, char **argv ) = 0;
  virtual string print ( void ) const = 0;
  virtual const Dimensions & dimensions() const {return Dimensions::one;}
};

class UserInterface
{
  private:
  
  const Context& ctxt_;
  int coutpe_;
  bool oncoutpe_;

  char *readCmd(char *s, int max, istream &fp, bool echo);
  bool terminate_;

  public: 

  list<Cmd*> cmdlist;
  list<Var*> varlist;

  // Used in SaveSysCmd to allow onfirstpe processor to output .sys file 
  // corresponding to a given qbLink input file.  
  void set_coutpe(int num) { 
    coutpe_ = num;
    oncoutpe_ = (ctxt_.mype() == coutpe_);
  }
  int get_coutpe(void) {
   return coutpe_;
  } 

  void addCmd(Cmd *newcmd)
  {
    newcmd->ui = this;
    cmdlist.push_back( newcmd );
  };

  Cmd *findCmd(char const*cmdname)
  {
    list<Cmd*>::iterator cmd;
    for ( cmd = cmdlist.begin();
          (cmd != cmdlist.end() && (strcmp((*cmd)->name(),cmdname)));
          cmd++ );

    if ( cmd != cmdlist.end() )
    {
      return (*cmd);
    }
    else
    {
      return 0;
    }
  };

  void addVar(Var *newvar)
  {
    newvar->ui = this;
    varlist.push_back( newvar );
  };

  Var *findVar(char const*varname)
  {
    list<Var*>::iterator var;
    for ( var = varlist.begin();
          (var != varlist.end() && (strcmp((*var)->name(),varname)));
          var++ );

    if ( var != varlist.end() )
    {
      return (*var);
    }
    else
    {
      return 0;
    }
  };

  void processCmds(istream &cmdstream, char const*prompt, bool echo, bool interactive = true);
  
  void terminate(void) { terminate_ = true; };

  bool oncoutpe(void) const { return oncoutpe_; };
  
  UserInterface(const Context& ctxt);

  void error(const string & message) const {
    if(oncoutpe()) cout << endl << "<ERROR>" << endl << " " << message << endl << "</ERROR>"<< endl << endl;
  }
  
  void warning(const string & message) const {
    if(oncoutpe()) cout << endl << "<WARNING>" << endl << " " << message << endl << "</WARNING>"<< endl << endl;
  }
  
};

// Local Variables:
// mode: c++
// End:

#endif

// Local Variables:
// mode: c++
// End:
