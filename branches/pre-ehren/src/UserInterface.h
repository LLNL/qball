////////////////////////////////////////////////////////////////////////////////
//
// UserInterface.h:
//
////////////////////////////////////////////////////////////////////////////////
// $ Id: $

#ifndef USER_INTERFACE_H
#define USER_INTERFACE_H

#include <iostream>
#include <string>
#include <cstring>
#include <iomanip>
#include <list>
#include <algorithm>
#include "Context.h"
using namespace std;

class UserInterface;

class Cmd
{
  public:
  UserInterface *ui;
  virtual char *name(void) const = 0;
  virtual char *help_msg(void) const = 0;
  virtual int action(int argc, char **argv) = 0;
};

class Var
{
  public:
  UserInterface *ui;
  virtual char *name ( void ) const = 0;
  virtual int set ( int argc, char **argv ) = 0;
  virtual string print ( void ) const = 0;
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

  Cmd *findCmd(char *cmdname)
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

  Var *findVar(char *varname)
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

  void processCmds(istream &cmdstream, char *prompt, bool echo);
  
  void terminate(void) { terminate_ = true; };

  bool oncoutpe(void) { return oncoutpe_; };
  
  UserInterface(const Context& ctxt);
};
#endif
