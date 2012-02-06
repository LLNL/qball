////////////////////////////////////////////////////////////////////////////////
//
// UserInterface.C: definition of readCmd and processCmds
//
////////////////////////////////////////////////////////////////////////////////
// $Id: UserInterface.C,v 1.7 2009/03/25 22:30:34 draeger1 Exp $

#include "UserInterface.h"
#include "Context.h"
#include <string>
#include <list>
#include <unistd.h> // isatty
#include <fstream>
#include <cassert>
using namespace std;

#if USE_MPI
#include <mpi.h>
#else
typedef int MPI_Comm;
#endif

////////////////////////////////////////////////////////////////////////////////
UserInterface::UserInterface(const Context& ctxt) : ctxt_(ctxt),terminate_(false),coutpe_(ctxt_.coutpe()),oncoutpe_(ctxt_.oncoutpe())
{
}
////////////////////////////////////////////////////////////////////////////////
char *UserInterface::readCmd(char *s, int max, istream &fp, bool echo)
{
  int ch, i = 0;
  while ( (ch = fp.get()) != EOF && !( ch == '\n' || ch ==';' || ch == '#') )
  {
    if ( ch == '\\' ) // line continuation character
    {
      // check if backslash is followed by a newline
      ch = fp.get();
      if ( ch == '\n' )
      {
        // backslash followed by newline, do nothing
      }
      else
      {
        // backslash not followed by newline
        if ( i < max - 1 )
          s[i++] = '\\';
        if ( i < max - 1 )
          s[i++] = ch;
      }
    }
    else
    {
      if (i < max - 1)
        s[i++] = ch;
    }
  }
  if (max > 0) s[i] = '\0';  /* add terminating NULL */

  if ( !(ch == '\n' || ch == ';' || ch == '#') )
    return NULL;             /* return NULL for end of file */
    
  // output command line if reading from a script
  if ( echo ) cout << s;
  
  if ( ch == '#' )
  {
    if ( echo ) cout << '#';
    while ( (ch = fp.get()) != EOF && !( ch == '\n' ) )
    {
      if ( echo ) cout << (char) ch;
    }
    if ( !(ch == '\n') )
      return NULL;             /* return NULL for end of file */
  }
  
  return s;
}

///////////////////////////////////////////////////////////////////////////////
void UserInterface::processCmds ( istream &cmdstream, char *prompt, bool echo )
{
  // read and process commands from cmdstream until end of file is reached

  char cmdline[256];
  list<Cmd*>::iterator cmd;
  char *tok;
  const char *separators = " ;\t";
  int i,done,status;
  if ( ctxt_.oncoutpe() )
    cout << "<!-- " << prompt << " ";
    
  // read a command terminated by '\n' or ';'
  if ( ctxt_.oncoutpe() )
  {
    done = terminate_ || !readCmd(cmdline, 256, cmdstream, echo );
    // readCmd returns cmdline if a command is read, NULL at EOF
  }
  string cmdlinestr(cmdline);
  ctxt_.string_bcast(cmdlinestr,ctxt_.coutpe());
  cmdlinestr.copy(cmdline,cmdlinestr.length(),0);
  cmdline[cmdlinestr.length()] = 0;

  if ( ctxt_.oncoutpe() ) {
    ctxt_.ibcast_send(1,1,&done,1);
  }
  else {
    // calculate row and col indices of process oncoutpe
    int irow = ctxt_.coutpe();
    while (irow >= ctxt_.nprow())
      irow -= ctxt_.nprow();
    int icol = int (ctxt_.coutpe()/ctxt_.nprow());
    assert(ctxt_.pmap(irow,icol) == ctxt_.coutpe());
            
    ctxt_.ibcast_recv(1,1,&done,1,irow,icol);
  }

  while ( !done )
  {
    if ( ctxt_.oncoutpe() )
      cout << " -->" << endl;
      
    // cmdline contains a string of tokens terminated by '\0'
    // cout << " command line is: " << cmdline << endl;

    // comment lines: start with '#'
    if ( cmdline[i=strspn(cmdline," ")] == '#' )
    {
      // cout << " comment line" << endl;
      // do nothing, not even write prompt
    }
    else if ( cmdline[i=strspn(cmdline," ")] == '!' )
    {
      // shell escape commands start with '!'
      // cout << " shell escape" << endl;
      if ( ctxt_.oncoutpe() )
      {
        cout << "<!-- ";
        system ( &cmdline[i+1] );
        cout << prompt << " ";
      }
    }
    else
    {
      // cout << " command split in the following tokens:" << endl;

      // scan tokens and build argument list
      list<char*> arglist;
      int ntok = 0;
      tok = strtok(cmdline, separators);
      while ( tok != 0 )
      {
        arglist.push_back(tok);
        ntok++;
        // cout << "\"" << tok << "\"" << endl;
        tok = strtok(0,separators);
      }
      // cout << " total of " << ntok << " tokens" << endl;
      // arglist.dump();

      // build ac and av
      int ac = ntok;
      char **av = new char *[ntok+1];
      i = 0;
      list<char*>::iterator iarg = arglist.begin();
      while ( iarg != arglist.end() )
      {
        av[i++] = *iarg++;
      }
      av[ntok] = 0;
  
      // write arguments
      //cout << "mype = " << ctxt_.mype() << ", cmdline = " << cmdline << endl;
      // for ( i = 0; i < ac; i++ )
      // {
      //   cout << "mype = " << ctxt_.mype() << ", av[" << i << "] = " << av[i] << endl;
      // }

      // search cmdlist for command

      tok = av[0]; 

      // check for empty command line
      if ( tok != 0 )
      {
        Cmd *cmdptr = findCmd(tok);

        if ( cmdptr )
        {
          ctxt_.barrier();
          cmdptr->action(ac,av);
          ctxt_.barrier();
        }
        else
        {
          // command is not in the command list, check for script files
          ifstream cmdstr;
          if ( ctxt_.oncoutpe() )
          {
            cmdstr.open(av[0],ios::in);
            status = !cmdstr;
          }

          if ( ctxt_.oncoutpe() ) {
            ctxt_.ibcast_send(1,1,&status,1);
          }
          else {
            // calculate row and col indices of process oncoutpe
            int irow = ctxt_.coutpe();
            while (irow >= ctxt_.nprow())
              irow -= ctxt_.nprow();
            int icol = int (ctxt_.coutpe()/ctxt_.nprow());
            assert(ctxt_.pmap(irow,icol) == ctxt_.coutpe());

            ctxt_.ibcast_recv(1,1,&status,1,irow,icol);
          }
          if ( !status )
          {
            // create new prompt in the form: prompt<filename>
            char *newprompt=0;
            if ( ctxt_.oncoutpe() )
            {
              newprompt = new char[strlen(prompt)+strlen(av[0])+4];
              newprompt = strcpy(newprompt,prompt);
              newprompt = strcat(newprompt,"[");
              newprompt = strcat(newprompt,av[0]);
              newprompt = strcat(newprompt,"]");
            }
              // MPI: process commands on all processes.
              // Note: newprompt is 0 on all processes > 0
              // Note: echo == true for scripts
            processCmds (cmdstr, newprompt, true);
            if ( ctxt_.oncoutpe() )
              delete newprompt;
          }
          else
          {
            if ( ctxt_.oncoutpe() )
              cout << "<WARNING> no such command or file name: " << tok << " </WARNING>"
                   << endl;
          }
        }
      }
      delete [] av;
      if ( ctxt_.oncoutpe() )      
        cout << "<!-- " << prompt << " ";
    }
    
    // read a command terminated by '\n' or ';'
    if ( ctxt_.oncoutpe() )
    {
      done = terminate_ || !readCmd(cmdline, 256, cmdstream, echo );
    }

    string cmdlinestr(cmdline);
    ctxt_.string_bcast(cmdlinestr,ctxt_.coutpe());
    cmdlinestr.copy(cmdline,cmdlinestr.length(),0);
    cmdline[cmdlinestr.length()] = 0;

    if ( ctxt_.oncoutpe() ) {
      ctxt_.ibcast_send(1,1,&done,1);
    }
    else {
      // calculate row and col indices of process oncoutpe
      int irow = ctxt_.coutpe();
      while (irow >= ctxt_.nprow())
        irow -= ctxt_.nprow();
      int icol = int (ctxt_.coutpe()/ctxt_.nprow());
      assert(ctxt_.pmap(irow,icol) == ctxt_.coutpe());

      ctxt_.ibcast_recv(1,1,&done,1,irow,icol);
    }

  }
  if ( ctxt_.oncoutpe() )          
    cout << " -->" << endl << "<!-- end of command stream -->" << endl;
}
