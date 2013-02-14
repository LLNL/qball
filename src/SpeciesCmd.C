////////////////////////////////////////////////////////////////////////////////
//
// SpeciesCmd.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SpeciesCmd.C,v 1.4 2009/03/25 22:30:34 draeger1 Exp $

#include "SpeciesCmd.h"
#include "SpeciesReader.h"
#include "Species.h"
using namespace std;

////////////////////////////////////////////////////////////////////////////////
int SpeciesCmd::action(int argc, char **argv) {
  if (! (argc == 3 || argc == 4)) {
    if ( ui->oncoutpe() )
      cout << "  <!-- use: species name uri [ewald_width] -->" << endl;
    return 1;
  }
  
  if ( ui->oncoutpe() )
    cout << "  <!-- SpeciesCmd: defining species " << argv[1]
         << " as " << argv[2] << " -->" << endl;

  SpeciesReader sp_reader(s->ctxt_);
  
  Species* sp = new Species(s->ctxt_,argv[1]);
  
  try {
    sp_reader.readSpecies(*sp,argv[2]);
    sp_reader.bcastSpecies(*sp);
    if (argc == 4) {
      const double rcpsin = atof(argv[3]);
      s->atoms.addSpecies(sp,argv[1],rcpsin);
    }
    else
      s->atoms.addSpecies(sp,argv[1]);

    if (sp->ultrasoft()) {
      s->ctrl.ultrasoft = true;
      s->wf.set_ultrasoft(true);
    }
    
    if (sp->nlcc()) {
       s->ctrl.nlcc = true;
    }
  }
  catch ( const SpeciesReaderException& e ) {
    cout << " SpeciesReaderException caught in SpeciesCmd" << endl;
    cout << " SpeciesReaderException: cannot define Species" << endl;
  }
  catch (...) {
    cout << " SpeciesCmd: cannot define Species" << endl;
  }
  
  return 0;
}
