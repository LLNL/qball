////////////////////////////////////////////////////////////////////////////////
//
// BOSampleStepper.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: BOSampleStepper.C,v 1.76 2010/08/30 20:55:43 draeger1 Exp $

#include "BOSampleStepper.h"
#include "EnergyFunctional.h"
#include "SlaterDet.h"
#include "Basis.h"
#include "WavefunctionStepper.h"
#include "SDWavefunctionStepper.h"
#include "JDWavefunctionStepper.h"
#include "PSDWavefunctionStepper.h"
#include "PSDAWavefunctionStepper.h"
#include "SOTDWavefunctionStepper.h"
#include "SORKTDWavefunctionStepper.h"
#include "FORKTDWavefunctionStepper.h"
#include "TDEULERWavefunctionStepper.h"
#include "SDIonicStepper.h"
#include "SDAIonicStepper.h"
#include "CGIonicStepper.h"
#include "MDIonicStepper.h"
#include "BMDIonicStepper.h"
#include "SDCellStepper.h"
#include "Preconditioner.h"
#include "AndersonMixer.h"
#include "MLWFTransform.h"
#include "SimpleConvergenceDetector.h"
#include "Hugoniostat.h"
#include "PrintMem.h"
#include "profile.h"
#ifdef USE_APC
#include "apc.h"
#endif
#include <iostream>
#include <iomanip>
#include <deque>
#ifdef HPM
#include <bgpm/include/bgpm.h>
extern "C" void HPM_Start(char *);
extern "C" void HPM_Stop(char *);
#endif
using namespace std;

////////////////////////////////////////////////////////////////////////////////
BOSampleStepper::BOSampleStepper(Sample& s, int nitscf, int nite) :
  SampleStepper(s), cd_(s), ef_(s,s.wf,cd_),dwf(s.wf), wfv(s.wfv), nitscf_(nitscf),
  nite_(nite), initial_atomic_density(false)
{
   const string wf_dyn = s_.ctrl.wf_dyn;
   tddft_involved_ = s_.ctrl.tddft_involved;
   if (tddft_involved_)
   {
      Wavefunction* wf0 = new Wavefunction(dwf);
      Wavefunction* wf1 = new Wavefunction(s.wf);
      Wavefunction* wf2 = new Wavefunction(s.wf);
      Wavefunction* wf3 = new Wavefunction(s.wf);
      Wavefunction* wf4 = new Wavefunction(s.wf);

      wfdeque.push_back ( wf0 );
      wfdeque.push_back ( wf1 );
      wfdeque.push_back ( wf2 );
      wfdeque.push_back ( wf3 );
      wfdeque.push_back ( wf4 );
   }
}
////////////////////////////////////////////////////////////////////////////////
BOSampleStepper::~BOSampleStepper()
{
   // delete any created Wavefunctions
   if (tddft_involved_)
      for (deque<Wavefunction*>::const_iterator it = wfdeque.begin(); it != wfdeque.end(); ++it)
         delete *it;
}
////////////////////////////////////////////////////////////////////////////////
void BOSampleStepper::initialize_density(void)
{
  // initialize cd_ with a sum of atomic densities

  double atom_radius[] =
      {
          0.00,  0.80,  0.59,  3.16,  2.12,  // null H  He  Li  Be
          1.64,  1.27,  1.06,  0.91,  0.79,  //  B   C   N   O   F
          0.72,  3.59,  2.74,  2.23,  2.10,  // Ne  Na  Mg  Al  Si
          1.85,  1.66,  1.49,  1.34,  4.59,  //  P   S  Cl  Ar   K
          3.67,  3.48,  3.33,  3.23,  3.14,  // Ca  Sc  Ti   V  Cr
          3.04,  2.95,  2.87,  2.82,  2.74,  // Mn  Fe  Co  Ni  Cu
          2.68,  2.57,  2.36,  2.15,  1.95,  // Zn  Ga  Ge  As  Se
          1.78,  1.66,  5.01,  4.14,  4.01,  // Br  Kr  Rb  Sr   Y
          3.89,  3.74,  3.59,  3.46,  3.36,  // Zr  Nb  Mo  Tc  Ru
          3.27,  3.19,  3.12,  3.04,  2.95,  // Rh  Pd  Ag  Cd  In
          2.74,  2.51,  2.32,  2.17,  2.04,  // Sn  Sb  Te   I  Xe
          5.63,  4.78,  0.00,  0.00,  4.67,  // Cs  Ba  La  Ce  Pr
          3.89,  3.87,  4.50,  4.37,  4.40,  // Nd  Pm  Sm  Eu  Gd
          4.25,  4.31,  0.00,  4.27,  4.20,  // Tb  Dy  Ho  Er  Tm
          4.20,  4.10,  3.93,  3.78,  3.65,  // Yb  Lu  Hf  Ta   W
          4.25,  4.31,  0.00,  4.27,  4.20,  // Tb  Dy  Ho  Er  Tm
          4.20,  4.10,  3.93,  3.78,  3.65,  // Yb  Lu  Hf  Ta   W
          3.55,  3.50,  3.40,  3.34,  3.29,  // Re  Os  Ir  Pt  Au
          3.23,  2.95,  2.91,  2.70,  2.55,  // Hg  Tl  Pb  Bi  Po
          4.00,  2.27,  4.00,  4.00,  4.00,  // At  Rn  Fr  Ra  Ac
          4.00,  4.00,  4.00,  4.00,  4.00,  // Th  Pa   U  Np  Pu
          4.00,  4.00,  4.00,  4.00,  4.00,  // Am  Cm  Bk  Cf  Es
          4.00,  4.00,  4.00,  4.00          // Fm  Md  No  Lr
      };
  

  if (!s_.wf.hasdata())
    s_.wf.set_hasdata(true);

  const AtomSet& atoms = s_.atoms;
  const Basis* const vbasis = cd_.vbasis();
  const int ngloc = vbasis->localsize();
  vector<vector<complex<double> > > rhops;
  const int nsp = atoms.nsp();
  rhops.resize(nsp);
  vector<complex<double> > rhopst(ngloc);
  const double * const g2 = vbasis->g2_ptr();
  
  for ( int is = 0; is < nsp; is++ )
  {
    rhops[is].resize(ngloc);
    Species *s = atoms.species_list[is];
    const int zval = s->zval();
    const int atomic_number = s->atomic_number();
    assert(atomic_number < sizeof(atom_radius)/sizeof(double));
    // scaling factor 2.0 in next line is empirically adjusted
    double rc = 2.0 * atom_radius[atomic_number];
    for ( int ig = 0; ig < ngloc; ig++ )
    {
      double arg = 0.25 * rc * rc * g2[ig];
      rhops[is][ig] = zval * exp( -arg );
    }
  }
  
  vector<vector<double> > tau0;
  tau0.resize(nsp);
  for ( int is = 0; is < nsp; is++ )
    tau0.resize(3*atoms.na(is));
  atoms.get_positions(tau0);
  StructureFactor sf;
  sf.init(tau0,*vbasis);
  sf.update(tau0,*vbasis);
  
  memset( (void*)&rhopst[0], 0, 2*ngloc*sizeof(double) );
  for ( int is = 0; is < nsp; is++ )
  {
    complex<double> *s = &sf.sfac[is][0];
    for ( int ig = 0; ig < ngloc; ig++ )
    {
      const complex<double> sg = s[ig];
      rhopst[ig] += sg * rhops[is][ig];
    }
  }

  for (int ispin = 0; ispin < s_.wf.nspin(); ispin++) {
    //if (s_.wf.spinactive(ispin)) {
      cd_.rhog[ispin].resize(ngloc);
      for (int i = 0; i < ngloc; i++) {
        cd_.rhog[ispin][i] = rhopst[i];
      }
      //}
  }
  initial_atomic_density = true;
}

////////////////////////////////////////////////////////////////////////////////
void BOSampleStepper::step(int niter)
{
  const bool oncoutpe = s_.ctxt_.oncoutpe();

  const bool anderson_charge_mixing = ( s_.ctrl.charge_mix_ndim > 0 );
  
  // determine whether eigenvectors must be computed
  // eigenvectors are computed if explicitly requested with wf_diag==T
  // or if the SlaterDet has fractionally occupied states
  //const bool fractional_occ = (s_.wf.nspin() * s_.wf.nel() != 2 * s_.wf.nst());
  int occtest = (2 * s_.wf.nst()) - s_.wf.nspin() * s_.wf.nel();
  const bool fractional_occ = (occtest != 0 && occtest != 1);
  const bool compute_eigvec = fractional_occ || s_.ctrl.wf_diag == "T";
  const bool compute_mlwf = s_.ctrl.wf_diag == "MLWF";
  const bool compute_mlwfc = s_.ctrl.wf_diag == "MLWFC";
  enum ortho_type { GRAM, LOWDIN, ORTHO_ALIGN, RICCATI };

  if (fractional_occ && oncoutpe)
       cout << "<!-- BOSampleStepper:  fractional occupation detected. -->" << endl;
  else if (oncoutpe)
       cout << "<!-- BOSampleStepper:  fractional occupation not detected. -->" << endl;
   
  
  if (s_.ctrl.reshape_context)
    s_.wf.set_reshape_context(s_.ctrl.reshape_context);
  
  AtomSet& atoms = s_.atoms;
  Wavefunction& wf = s_.wf;
  const int nspin = wf.nspin();

  const UnitCell& cell = wf.cell();

  const double dt = s_.ctrl.dt;

  const string wf_dyn = s_.ctrl.wf_dyn;
  const string atoms_dyn = s_.ctrl.atoms_dyn;
  const string cell_dyn = s_.ctrl.cell_dyn;

  // AS: flag used below to control initial step for SOTD
  bool wfv_is_new = false;
  const bool extrapolate_wf = ( atoms_dyn == "MD" ) && ( nite_ == 1 ) && (!tddft_involved_);


  //ewd DEBUG
  //if (tddft_involved_ && oncoutpe)
  //   cout << "DEBUG:  TDDFT_INVOLVED SET TO TRUE!" << endl;
  //ewd DEBUG

  
  const bool ntc_extrapolation =
    s_.ctrl.debug.find("NTC_EXTRAPOLATION") != string::npos;
  const bool asp_extrapolation =
    s_.ctrl.debug.find("ASP_EXTRAPOLATION") != string::npos;

  //ewd check for case where PSDA used (incorrectly) with nite = 1 and empty states
  if (wf_dyn == "PSDA" && nite_ == 1 && compute_eigvec) {
    if ( oncoutpe ) {
      cout << "<ERROR> BOSampleStepper:  PSDA unstable with empty states and nite = 1. </ERROR>" << endl;
      cout << "<ERROR> BOSampleStepper:  Increase nite or use wf_dyn = PSD. </ERROR>" << endl;
    }
    return;
  }

  
  Wavefunction* wfmm;
  if ( extrapolate_wf && ( ntc_extrapolation || asp_extrapolation ) )
    wfmm = new Wavefunction(wf);

  // Next lines: special value of niter = 0: GS calculation only
  const bool atoms_move = ( niter > 0 && atoms_dyn != "LOCKED" );
  const bool compute_stress = ( s_.ctrl.stress == "ON" );
  const bool cell_moves = ( niter > 0 && compute_stress &&
                            cell_dyn != "LOCKED" );
  // GS-only calculation:
  const bool gs_only = !atoms_move && !cell_moves;
  const bool use_confinement = ( s_.ctrl.ecuts > 0.0 );

  const string charge_mixing = s_.ctrl.charge_mixing;

  const bool ultrasoft = s_.ctrl.ultrasoft;
  const bool usdiag = (ultrasoft && atoms_move);
  const bool nlcc = s_.ctrl.nlcc;
  cd_.set_nlcc(nlcc);
  
  //ewd check that MLWF not being used with ultrasoft (not yet implemented)
  if (ultrasoft && (compute_mlwf || compute_mlwfc)) {
    if ( oncoutpe ) 
      cout << "<ERROR> BOSampleStepper:  Maximally-localized Wannier Functions not yet implemented with ultrasoft. </ERROR>" << endl;
    return;
  }
  if (oncoutpe && ultrasoft && compute_stress)
      cout << "<WARNING> BOSampleStepper:  stress not yet implemented with ultrasoft!  Results WILL be wrong. </WARNING>" << endl;
  if (oncoutpe && nlcc && compute_stress)
      cout << "<WARNING> BOSampleStepper:  stress not yet implemented with non-linear core corrections!  Results WILL be wrong. </WARNING>" << endl;
  
  // use extra memory for SlaterDets if memory variable = normal, large or huge
  if (s_.ctrl.extra_memory >= 3) 
    wf.set_highmem();
  
  if (!gs_only && s_.ctrl.dft_plus_u && oncoutpe)
    cout << "<WARNING> Forces and stress not currently implemented with DFT+U! </WARNING>" << endl;
  
  Timer tm_iter;

  const bool use_preconditioner = wf_dyn == "PSD" || wf_dyn == "PSDA" || wf_dyn == "JD";
  Preconditioner *preconditioner = 0;
  if ( use_preconditioner )
  {
    // create a preconditioner using the information about wf in s_.wf
    // and the information about the hessian in df
    preconditioner = new Preconditioner(s_,s_.wf,ef_);
  }

  // initialize occupation
  wf.update_occ(s_.ctrl.smearing_width,s_.ctrl.smearing_ngauss);

  WavefunctionStepper* wf_stepper = 0;
  if ( wf_dyn == "SD" )
  {
    const double emass = s_.ctrl.emass;
    double dt2bye = (emass == 0.0) ? 0.5 / wf.ecut() : dt*dt/emass;

    // divide dt2bye by facs coefficient if stress == ON
    const double facs = 2.0;
    if ( s_.ctrl.stress == "ON" )
    {
      dt2bye /= facs;
    }
    wf_stepper = new SDWavefunctionStepper(wf,dt2bye,tmap);
  }
  else if ( wf_dyn == "TDEULER" )
     wf_stepper = new TDEULERWavefunctionStepper(wf,s_.ctrl.tddt,tmap);
  else if ( wf_dyn == "SOTD" )
  {
     // AS: prepare the array to store the |psi(t-tddt)> for the second-order propagation
     if ( s_.wfv == 0 )
     {
        s_.wfv = new Wavefunction(wf);
        s_.wfv->clear();
        wfv_is_new = true;
     }

     // AS: initially prepare the deque with the respective wave functions
     wfdeque.push_back ( s_.wfv );
     *wfdeque[0]=*(s_.wfv);

     wfdeque.push_back ( s_.wfv );
     *wfdeque[1]=wf;

     wf_stepper = new SOTDWavefunctionStepper(wf,s_.ctrl.tddt,tmap,&wfdeque);
  }
  else if ( wf_dyn == "SORKTD" )
     wf_stepper = new SORKTDWavefunctionStepper(wf,s_.ctrl.tddt,tmap,&wfdeque);
  else if ( wf_dyn == "FORKTD" )
     wf_stepper = new FORKTDWavefunctionStepper(wf,s_.ctrl.tddt,tmap,&wfdeque);
  else if ( wf_dyn == "PSD" )
    wf_stepper = new PSDWavefunctionStepper(wf,*preconditioner,tmap);
  else if ( wf_dyn == "PSDA" )
    wf_stepper = new PSDAWavefunctionStepper(wf,*preconditioner,tmap);  
  else if ( wf_dyn == "JD" )
    wf_stepper = new JDWavefunctionStepper(wf,*preconditioner,ef_,tmap);  
  // wf_stepper == 0 indicates that wf_dyn == LOCKED

  IonicStepper* ionic_stepper = 0;
  if ( atoms_dyn == "SD" )
    ionic_stepper = new SDIonicStepper(s_);
  else if ( atoms_dyn == "SDA" )
    ionic_stepper = new SDAIonicStepper(s_);
  else if ( atoms_dyn == "CG" )
    ionic_stepper = new CGIonicStepper(s_);
  // AS: IMPULSIVE does MD with constant velocities
  else if ( atoms_dyn == "MD" || atoms_dyn == "IMPULSIVE" )
    ionic_stepper = new MDIonicStepper(s_);
  else if ( atoms_dyn == "BMD" )
    ionic_stepper = new BMDIonicStepper(s_);

  if ( ionic_stepper )
    ionic_stepper->setup_constraints();

  CellStepper* cell_stepper = 0;
  if ( cell_dyn == "SD" )
    cell_stepper = new SDCellStepper(s_);

  // Allocate wavefunction velocity if not available
  if ( atoms_move && extrapolate_wf )
  {
    if ( s_.wfv == 0 )
    {
      s_.wfv = new Wavefunction(wf);
      s_.wfv->clear();
    }
  }

  MLWFTransform* mlwft=0;

  if ( compute_mlwf || compute_mlwfc )
  {
    // MLWF can be computed at the gamma point only
    // There must be a single k-point, and it must be gamma
    if ( wf.nkp() > 1 || ( wf.nkp()==1 && wf.kpoint(0) != D3vector(0,0,0) ) )
    {
      if ( oncoutpe )
      {
        cout << "<ERROR> BOSampleStepper::step: MLWF can be computed at k=0 only </ERROR>"
             << endl;
        cout << "<ERROR> BOSampleStepper::step: cannot run </ERROR>" << endl;
      }
      return;
    }
    if (wf.nspin() > 1 && s_.ctxt_.oncoutpe()) 
      cout << "<ERROR> nspin > 1!  MLWF doesn't currently work with spin-polarized systems </ERROR>" << endl;
    assert(wf.nspin()==1);
    mlwft = new MLWFTransform(*wf.sd(0,0));
  }

  //ewd:  experimental Hugoniostat for Kyle and Sebastien
  Hugoniostat* hugstat = 0;
  if (s_.ctrl.hugoniostat == "ON") {
    hugstat = new Hugoniostat(s_.ctrl.hug_etot,s_.ctrl.hug_volume,s_.ctrl.hug_pressure,ionic_stepper->temp(),s_.ctxt_.oncoutpe());
    hugstat->set_deltatemp(s_.ctrl.hug_deltatemp);
    hugstat->set_updatefreq(s_.ctrl.hug_freq);
  }

  // Charge mixing variables
  vector<vector<complex<double> > > rhog_in;
  vector<vector<complex<double> > > drhog;
  vector<vector<complex<double> > > rhobar;
  vector<vector<complex<double> > > drhobar;

  rhog_in.resize(nspin);
  drhog.resize(nspin);
  rhobar.resize(nspin);
  drhobar.resize(nspin);
  for (int ispin = 0; ispin < nspin; ispin++) {
    rhog_in[ispin].resize(cd_.rhog[ispin].size());
    drhog[ispin].resize(cd_.rhog[ispin].size());
    rhobar[ispin].resize(cd_.rhog[ispin].size());
    drhobar[ispin].resize(cd_.rhog[ispin].size());
  }

  vector<double> wkerker(rhog_in[0].size());
  vector<double> wls(rhog_in[0].size());

  const int anderson_ndim = s_.ctrl.charge_mix_ndim;

  vector<AndersonMixer*> mixer;
  mixer.resize(nspin);
  for (int ispin = 0; ispin < nspin; ispin++)
    //if (s_.wf.spinactive(ispin))
      mixer[ispin] = new AndersonMixer(2*rhog_in[ispin].size(),anderson_ndim,&cd_.vcontext());
      
  // compute Kerker preconditioning
  // real space Kerker cutoff in a.u.
  const double rc_Kerker = s_.ctrl.charge_mix_rcut;
  const double *const g2 = cd_.vbasis()->g2_ptr();

  // define q1 cutoff for row weighting of LS charge mixing
  // Use rc1 = 3 a.u. default cutoff
  double rc1 = 3.0;
  // check if override from the debug variable
  // use: set debug RC1 <value>
  if ( s_.ctrl.debug.find("RC1") != string::npos )
  {
    istringstream is(s_.ctrl.debug);
    string s;
    is >> s >> rc1;
    cout << " override rc1 value: rc1 = " << rc1 << endl;
    assert(rc1 >= 0.0);
  }

  if ( rc1 != 0.0 )
  {
    const double q1 = 2.0 * M_PI / rc1;
    for ( int i=0; i < wls.size(); i++ )
    {
      if ( g2[i] != 0.0 )
        wls[i] = sqrt(g2[i] / ( g2[i] + q1*q1 ));
      else
        wls[i] = 1.0;
    }
  }
  else
  {
    for ( int i=0; i < wls.size(); i++ )
      wls[i] = 1.0;
  }

  if ( rc_Kerker > 0.0 )
  {
    const double q0_kerker = 2 * M_PI / rc_Kerker;
    const double q0_kerker2 = q0_kerker * q0_kerker;
    for ( int i=0; i < wkerker.size(); i++ )
      wkerker[i] = g2[i] / ( g2[i] + q0_kerker2 );
  }
  else
  {
    for ( int i=0; i < wkerker.size(); i++ )
      wkerker[i] = 1.0;
  }

  // if ultrasoft, calculate position-dependent functions
  if (ultrasoft) {
     tmap["usfns"].start();
     cd_.update_usfns();
     wf.update_usfns();
     tmap["usfns"].stop();
  }
  // if non-linear core correction defined, calculate position-dependent density
  if (nlcc)
     cd_.update_nlcc();
  
  double tmpthresh = (atoms_move ? s_.ctrl.threshold_force : 0.0 );
  SimpleConvergenceDetector conv_force(s_.ctrl.threshold_force_nsteps, tmpthresh);
  tmpthresh = (compute_stress ? s_.ctrl.threshold_stress : 0.0 );
  SimpleConvergenceDetector conv_stress(s_.ctrl.threshold_stress_nsteps, tmpthresh);

  // calculate time available to avoid exceeding run_timer
  double tbase, tleft;
  bool testtimer = true;
  if (s_.ctrl.run_timer > 0.0) {
    tbase = MPI_Wtime();
    tleft = s_.ctrl.run_timer - (tbase - s_.ctrl.time_init);
  }
  
  // Next line: special case of niter=0: compute GS only        
  for ( int iter = 0; iter < max(niter,1); iter++ )
  {

    // check timing
    if (s_.ctrl.run_timer > 0.0 && niter > 1 && iter > 1 && testtimer) {
      double tnow = MPI_Wtime();
      double sofar = tnow - tbase;
      double tleft = s_.ctrl.run_timer - ( tnow - s_.ctrl.time_init); 
      double titer = sofar/iter; // avg. time per iteration
      double tmaxiter;
      MPI_Allreduce(&titer, &tmaxiter, 1, MPI_DOUBLE, MPI_MAX, s_.ctxt_.comm());
      double maxtime = tmaxiter*(niter-iter);
      if (maxtime > tleft) { // we'll exceed timer, lower niter
        s_.ctrl.timer_hit = true;
        int newiter = (int)(tleft/tmaxiter);
        int tmpniter = iter + newiter - 1;
        if (tmpniter != niter) {
          niter = tmpniter;
          if ( s_.ctxt_.oncoutpe() )
            cout << "<!-- estimated ionic iteration time = " << setprecision(3) << tmaxiter << " sec will exceed run_timer, lowering niter to " << niter << " -->" << endl;
        }
        else {
          testtimer = false;
          if ( s_.ctxt_.oncoutpe() )
            cout << "<!-- estimated ionic iteration time = " << setprecision(3) << tmaxiter << " sec is stable, using niter = " << niter << " -->" << endl;
        }
      }
    }
    
    // ionic iteration
    // test convergence of forces and/or stress
    bool ionic_converge = true;
    if (atoms_move && !conv_force.isConverged())
      ionic_converge = false;
    if (compute_stress && !conv_stress.isConverged())
      ionic_converge = false;
    if (gs_only)
      ionic_converge = false;

    if (ionic_converge) {
      if ( s_.ctxt_.oncoutpe() ) {
        cout.setf(ios::scientific,ios::floatfield);
        cout << "  <!-- BOSampleStepper: ionic convergence reached:  -->" << endl;
        if (atoms_move)
          cout << "  <!-- BOSampleStepper:  ionic convergence, maximum forces varied by less than " << conv_force.threshold() << " a.u. over last " << conv_force.nsteps() << " steps -->" << endl;
        if (compute_stress)
          cout << "  <!-- BOSampleStepper:  ionic convergence, maximum stress varied by less than " << conv_stress.threshold() << " GPa over last " << conv_stress.nsteps() << " steps -->" << endl;
      }
      iter = niter;
    }
    else {

      // check for change in basis size which can cause errors in charge mixing between iterations
      for (int ispin = 0; ispin < nspin; ispin++) {
        //if (s_.wf.spinactive(ispin)) {
          if (rhog_in[ispin].size() != cd_.rhog[ispin].size()) {
            if ( s_.ctxt_.oncoutpe() ) 
              cout << "<ERROR> ChargeDensity basis mismatch!  A ref_cell may be needed. </ERROR>" << endl;
            s_.ctxt_.abort(1);
          }
          //}
      }

      if (iter > 0) {
        // calculate max force and stress values, add to convergence detectors
        if (atoms_move) {
          double maxforce = 0.0;
          for ( int is = 0; is < atoms.atom_list.size(); is++ ) 
            for ( int ia = 0; ia < atoms.atom_list[is].size(); ia++ ) 
              for (int q=0; q<3; q++) 
                if (fabs(fion[is][3*ia+q]) > maxforce) maxforce = fabs(fion[is][3*ia+q]);

          conv_force.addValue(maxforce);
        }
        if (compute_stress) {
          // ewd:  add more significant figures to conversion
          const double gpa = 29421.0120;
          //const double gpa = 29421.5;
          double maxstress = 0.0;
          for (int q=0; q<6; q++) 
            if (fabs(sigma_eks[q]) > maxstress) maxstress = fabs(sigma_eks[q]);

          conv_stress.addValue(gpa*maxstress);
        }
      }
      

      tm_iter.start();
#ifdef USE_APC
      ApcStart(1);
#endif
      
      if ( oncoutpe )
        cout << "<iteration count=\"" << iter+1 << "\">\n";

      if ( ionic_stepper )
        atoms.sync();

      // compute energy and ionic forces using existing wavefunction
      
      if ( !gs_only )
      {
        tmap["charge"].start();
        if ( initial_atomic_density ) {
          cd_.update_rhor();
          if (ultrasoft)
            //ewd DEBUG:  need to renormalize ultrasoft density so tot_charge is correct?
            assert(false);
        }
        else 
          cd_.update_density();
        if (tddft_involved_)
           ( ef_.hamil_cd() )->update_density();
        tmap["charge"].stop();

        if (tddft_involved_)
           ef_.update_hamiltonian();
        ef_.update_vhxc();
        const bool compute_forces = true;

        // if symmetry is used, need to calculate set of symmetry-equivalent atoms for 
        // force averaging
        if ( s_.symmetries.nsym() > 0) {
          atoms.findSymmetricAtoms(s_.symmetries);
        }

        // need eigenvalues to compute forces w. ultrasoft
        if (ultrasoft) { 
          tmap["gram"].start();
          s_.wf.gram();
          tmap["gram"].stop();
          ef_.energy(true,dwf,false,fion,false,sigma_eks);
          tmap["diag"].start();
          s_.wf.diag(dwf,true);
          tmap["diag"].stop();
          tmap["usfns"].start();
          s_.wf.update_usfns();
          tmap["usfns"].stop();
          s_.wf.printeig();
        }
        
        double energy =
            ef_.energy(false,dwf,compute_forces,fion,compute_stress,sigma_eks);

        // average forces over symmetric atoms
        if ( compute_forces && s_.symmetries.nsym() > 0) {
          const int nsym_ = s_.symmetries.nsym();
          for ( int is = 0; is < atoms.atom_list.size(); is++ ) {
            for ( int ia = 0; ia < atoms.atom_list[is].size(); ia++ ) {
              // start with identity symmetry operation
              D3vector fsym(fion[is][3*ia],fion[is][3*ia+1],fion[is][3*ia+2]);
              for ( int isym = 0; isym < nsym_; isym++) {
                int ja = atoms.symatomid(is,ia,isym);
                D3vector ftmp(fion[is][3*ja],fion[is][3*ja+1],fion[is][3*ja+2]);
                //fsym = fsym + s_.symmetries.symlist[isym]->applyToVector(ftmp,false);
                //ewd need to convert force to crystal coordinates, as symmetry matrix
                // is in crystal coords
                D3vector ftmp_xtal = cell.cart_to_crystal(ftmp);
                D3vector fsym_xtal = s_.symmetries.symlist[isym]->applyToVector(ftmp_xtal,false);
                fsym = fsym + cell.crystal_to_cart(fsym_xtal);
              }
              fion[is][3*ia] = fsym.x/(double)(nsym_+1);
              fion[is][3*ia+1] = fsym.y/(double)(nsym_+1);
              fion[is][3*ia+2] = fsym.z/(double)(nsym_+1);
            }
          }       
        }

        if (compute_forces && atoms.add_fion_ext()) {
          for ( int is = 0; is < atoms.atom_list.size(); is++ ) {
            for ( int ia = 0; ia < atoms.atom_list[is].size(); ia++ ) {
              D3vector ftmp = atoms.get_fion_ext(is,ia);
              fion[is][3*ia] += ftmp.x;
              fion[is][3*ia+1] += ftmp.y;
              fion[is][3*ia+2] += ftmp.z;
            }
          }
        }

        if ( oncoutpe )
        {
          cout.setf(ios::fixed,ios::floatfield);
          cout.setf(ios::right,ios::adjustfield);
          cout << "  <ekin>   " << setprecision(8)
               << setw(15) << ef_.ekin() << " </ekin>\n";
          if ( use_confinement )
            cout << "  <econf>  " << setw(15) << ef_.econf() << " </econf>\n";
          cout << "  <eps>    " << setw(15) << ef_.eps() << " </eps>\n"
               << "  <enl>    " << setw(15) << ef_.enl() << " </enl>\n"
               << "  <ecoul>  " << setw(15) << ef_.ecoul() << " </ecoul>\n"
               << "  <exc>    " << setw(15) << ef_.exc() << " </exc>\n"
               << "  <esr>    " << setw(15) << ef_.esr() << " </esr>\n"
               << "  <eself>  " << setw(15) << ef_.eself() << " </eself>\n"
               << "  <ets>    " << setw(15) << ef_.ets() << " </ets>\n"
               << "  <etotal> " << setw(15) << ef_.etotal() << " </etotal>\n";
          if ( compute_stress )
          {
            const double pext = (sigma_ext[0]+sigma_ext[1]+sigma_ext[2])/3.0;
            const double enthalpy = ef_.etotal() + pext * cell.volume();
            cout << "  <pv>     " << setw(15) << pext * cell.volume()
                 << " </pv>" << endl;
            cout << "  <enthalpy> " << setw(15) << enthalpy << " </enthalpy>\n"
                 << flush;
          }
        }
        
        // Add empirical forces (if any) to ions
        double epot_mmion = 0.0;
        if (atoms.empirical_list.size() > 0) {
          const int nspqm_ = atoms.atom_list.size();
          for (int k=0; k<atoms.empirical_list.size(); k++) {
          
            int is1 = atoms.empirical_list[k]->is1;
            int is2 = atoms.empirical_list[k]->is2;

            double prefac = 1.0;
            if (is1 == is2) prefac = 0.5;
            
            vector<Atom*>::iterator pa1; 
            vector<Atom*>::iterator pa2; 
            if (is1 >= nspqm_) 
              pa1 = atoms.mmatom_list[is1-nspqm_].begin();
            else 
              pa1 = atoms.atom_list[is1].begin();
            for (int ia=0; ia<atoms.na(is1); ia++) {
              if (is2 >= nspqm_)
                pa2 = atoms.mmatom_list[is2-nspqm_].begin();
              else 
                pa2 = atoms.atom_list[is2].begin();
              for (int ja=0; ja<atoms.na(is2); ja++) {
                D3vector r12 = (*pa1)->position() - (*pa2)->position();
                cell.fold_in_ws(r12);
                
                double r12val = length(r12);
                if (r12val > 0.0) 
                  epot_mmion += prefac*atoms.empirical_list[k]->pot(r12val);
              
                if ( compute_forces) {
                  D3vector emp_force = prefac*atoms.empirical_list[k]->force(r12);
                  fion[is1][3*ia] += emp_force.x;
                  fion[is1][3*ia+1] += emp_force.y;
                  fion[is1][3*ia+2] += emp_force.z;
                  fion[is2][3*ja] -= emp_force.x;
                  fion[is2][3*ja+1] -= emp_force.y;
                  fion[is2][3*ja+2] -= emp_force.z;
                }
                pa2++;
              }
              pa1++;
            }
          }
        }
      
        if ( iter > 0 && ionic_stepper )
        {
          ionic_stepper->compute_v(energy,fion);
        }
        // at this point, positions r0, velocities v0 and forces fion are
        // consistent
        double ekin_ion = 0.0, temp_ion = 0.0;
        if ( ionic_stepper )
        {
          ekin_ion = ionic_stepper->ekin();
          temp_ion = ionic_stepper->temp();
        }
        
        // if hugoniostat is turned on, add current total energy, pressure and volume
        if (s_.ctrl.hugoniostat == "ON") { 
          double pext = 0.0;
          if ( compute_stress )
            pext = (sigma_ext[0]+sigma_ext[1]+sigma_ext[2])/3.0;
          hugstat->addValues(ef_.etotal(),cell.volume(),pext,ionic_stepper->temp());
          if (hugstat->updatenow) {
            double currthtemp = s_.ctrl.th_temp;
            hugstat->updateTemp(currthtemp);
            s_.ctrl.th_temp = currthtemp;
          }
        }

        // print positions, velocities and forces at time t0
        if ( oncoutpe )
        {
          cout << "<atomset>" << endl;
          cout << atoms.cell();
          for ( int is = 0; is < atoms.atom_list.size(); is++ )
          {
            int i = 0;
            for ( int ia = 0; ia < atoms.atom_list[is].size(); ia++ )
            {
              Atom* pa = atoms.atom_list[is][ia];
              cout << "  <atom name=\"" << pa->name() << "\""
                   << " species=\"" << pa->species()
                   << "\">\n"
                   << "    <position> " << pa->position() << " </position>\n"
                   << "    <velocity> " << pa->velocity() << " </velocity>\n"
                   << "    <force> "
                   << fion[is][i] << " "
                   << fion[is][i+1] << " "
                   << fion[is][i+2]
                   << " </force>\n";
              cout << "  </atom>" << endl;
              i += 3;
            }
          }
          // MM atoms
          for ( int is = 0; is < atoms.mmatom_list.size(); is++ ) {
            int i = 0;
            const int offset = atoms.atom_list.size();
            for ( int ia = 0; ia < atoms.mmatom_list[is].size(); ia++ ) {
              Atom* pa = atoms.mmatom_list[is][ia];
              cout << "  <mmatom name=\"" << pa->name() << "\""
                   << " species=\"" << pa->species()
                   << "\">\n"
                   << "    <position> " 
                   << ionic_stepper->r0(is+offset,i) << " "
                   << ionic_stepper->r0(is+offset,i+1) << " " 
                   << ionic_stepper->r0(is+offset,i+2) << " </position>\n"
                   << "    <velocity> " 
                   << ionic_stepper->v0(is+offset,i) << " "
                   << ionic_stepper->v0(is+offset,i+1) << " " 
                   << ionic_stepper->v0(is+offset,i+2) << " </velocity>\n"
                   << "    <force> " 
                   << fion[is+offset][i] << " "
                   << fion[is+offset][i+1] << " " 
                   << fion[is+offset][i+2]
                   << " </force>\n  </mmatom>" << endl;
              
              i += 3;
            }
          }
          cout << "</atomset>" << endl;
          cout << "  <econst> " << energy+ekin_ion << " </econst>\n";
          if (atoms.mmatom_list.size() > 0) 
            cout << "  <epot_mmion> " << epot_mmion << " </epot_mmion>\n";
          cout << "  <ekin_ion> " << ekin_ion << " </ekin_ion>\n";
          cout << "  <temp_ion> " << temp_ion << " </temp_ion>\n";
        }
        
        if ( atoms_move )
        {
           if ( s_.constraints.size() > 0 )
          {
            s_.constraints.compute_forces(ionic_stepper->r0(), fion);
            if ( oncoutpe )
            {
              s_.constraints.list_constraints(cout);
            }
          }
          // move atoms to new position: r0 <- r0 + v0*dt + dt2/m * fion
          ionic_stepper->compute_r(energy,fion);
          ef_.atoms_moved();
          if (ultrasoft) {
             tmap["usfns"].start();
             cd_.update_usfns();
             wf.update_usfns();
             tmap["usfns"].stop();
          }
          tmap["charge"].start();
          if (nlcc)
             cd_.update_nlcc();
          tmap["charge"].stop();
        }

        if ( compute_stress )
        {
          compute_sigma();
          print_stress();

          if ( cell_moves )
          {
            cell_stepper->compute_new_cell(sigma);

            // Update cell
            cell_stepper->update_cell();

            ef_.cell_moved();
            ef_.atoms_moved(); // modifications of the cell also move ions
            if (ultrasoft) {
              tmap["usfns"].start();
              cd_.update_usfns();
              wf.update_usfns();
              tmap["usfns"].stop();
            }
            tmap["charge"].start();
            if (nlcc)
               cd_.update_nlcc();
            tmap["charge"].stop();

            if ( use_preconditioner )
              preconditioner->update();
          }
        }
      } // if !gs_only

      // Recalculate ground state wavefunctions
      
      // wavefunction extrapolation
      if ( atoms_move && extrapolate_wf )
      {
        for ( int ispin = 0; ispin < nspin; ispin++ )
        {
          if (s_.wf.spinactive(ispin))
          {
            for ( int ikp = 0; ikp < s_.wf.nkp(); ikp++ )
            {
              if (s_.wf.kptactive(ikp))
              {
                assert(s_.wf.sd(ispin,ikp) != 0);
                if ( ntc_extrapolation )
                {
                  double* c = (double*) s_.wf.sd(ispin,ikp)->c().cvalptr();
                  double* cv = (double*) s_.wfv->sd(ispin,ikp)->c().cvalptr();
                  double* cmm = (double*) wfmm->sd(ispin,ikp)->c().cvalptr();
                  const int mloc = s_.wf.sd(ispin,ikp)->c().mloc();
                  const int nloc = s_.wf.sd(ispin,ikp)->c().nloc();
                  const int len = 2*mloc*nloc;
                  if ( iter == 0 )
                  {
                    // copy c on cv
                    for ( int i = 0; i < len; i++ )
                    {
                      const double x = c[i];
                      const double v = cv[i];
                      // extrapolation using velocity in cv
                      c[i] = x + dt * v;
                      cv[i] = x;
                    }
                    if (ultrasoft) {
                      tmap["usfns"].start();
                      s_.wf.sd(ispin,ikp)->update_usfns();
                      tmap["usfns"].stop();
                    }
                    tmap["gram"].start();
                    s_.wf.sd(ispin,ikp)->gram();
                    tmap["gram"].stop();
                  }
                  else if ( iter == 1 )
                  {
                     s_.wfv->align(s_.wf);
                    for ( int i = 0; i < len; i++ )
                    {
                      const double x = c[i];
                      const double xm = cv[i];
                      c[i] = 2.0 * x - xm;
                      cv[i] = x;
                      cmm[i] = xm;
                    }
                    if (ultrasoft) {
                      tmap["usfns"].start();
                      s_.wf.sd(ispin,ikp)->update_usfns();
                      tmap["usfns"].stop();
                    }
                    tmap["gram"].start();
                    s_.wf.sd(ispin,ikp)->gram();
                    tmap["gram"].stop();
                  }
                  else
                  {
                    // align wf with wfmm before extrapolation
                    // s_.wf.align(*wfmm);
                    wfmm->align(s_.wf);
                  
                    // extrapolate
                    for ( int i = 0; i < len; i++ )
                    {
                      const double x = c[i];   // current wf (scf converged) at t
                      const double xm = cv[i]; // extrapolated wf at t
                      const double xmm = cmm[i]; // extrapolated wf at t-dt
                      c[i] = 2.0 * x - xmm;
                      // save extrapolated value at t in cmm
                      cmm[i] = xm;
                    }

                    // orthogonalize the extrapolated value
                    if (ultrasoft) {
                      tmap["usfns"].start();
                      s_.wf.sd(ispin,ikp)->update_usfns();
                      tmap["usfns"].stop();
                    }
                    tmap["gram"].start();
                    s_.wf.sd(ispin,ikp)->gram();
                    tmap["gram"].stop();
                    //tmap["lowdin"].start();
                    //s_.wf.sd(ispin,ikp)->lowdin();
                    //tmap["lowdin"].stop();
                  
                    // c[i] now contains the extrapolated value
                    // save a copy in cv[i]
                    for ( int i = 0; i < len; i++ )
                    {
                      cv[i] = c[i];
                    }
                  }
                  // c[i] is now ready for electronic iterations
                }
                else if ( asp_extrapolation )
                {
                  double* c = (double*) s_.wf.sd(ispin,ikp)->c().cvalptr();
                  double* cv = (double*) s_.wfv->sd(ispin,ikp)->c().cvalptr();
                  double* cmm = (double*) wfmm->sd(ispin,ikp)->c().cvalptr();
                  const int mloc = s_.wf.sd(ispin,ikp)->c().mloc();
                  const int nloc = s_.wf.sd(ispin,ikp)->c().nloc();
                  const int len = 2*mloc*nloc;
                  if ( iter == 0 )
                  {
                    for ( int i = 0; i < len; i++ )
                    {
                      const double x = c[i];
                      const double v = cv[i];
                      // extrapolation using velocity in cv
                      c[i] = x + dt * v;
                      cv[i] = x;
                    }
                    if (ultrasoft) {
                      tmap["usfns"].start();
                      s_.wf.sd(ispin,ikp)->update_usfns();
                      tmap["usfns"].stop();
                    }
                    tmap["gram"].start();
                    s_.wf.sd(ispin,ikp)->gram();
                    tmap["gram"].stop();
                  }
                  else if ( iter == 1 )
                  {
                    //s_.wfv->align(s_.wf);
                    for ( int i = 0; i < len; i++ )
                    {
                      const double x = c[i];
                      const double xm = cv[i];
                      c[i] = 2.0 * x - xm;
                      cv[i] = x;
                      cmm[i] = xm;
                    }
                    if (ultrasoft) {
                      tmap["usfns"].start();
                      s_.wf.sd(ispin,ikp)->update_usfns();
                      tmap["usfns"].stop();
                    }
                    tmap["gram"].start();
                    s_.wf.sd(ispin,ikp)->gram();
                    tmap["gram"].stop();
                  }
                  else
                  {
                    // align wf with wfmm before extrapolation
                    // s_.wf.align(*wfmm);
                    // wfmm->align(s_.wf);
                  
                    // extrapolate
                    for ( int i = 0; i < len; i++ )
                    {
                      const double x = c[i];   // current wf (scf converged) at t
                      const double xm = cv[i]; // extrapolated wf at t
                      const double xmm = cmm[i]; // extrapolated wf at t-dt
                      const double asp_a1 = 0.5;
                      c[i] = 2.0 * x - xm +
                          asp_a1 * ( x - 2.0 * xm + xmm );
                      //c[i] = 2.5 * x - 2.0 * xm + 0.5 * xmm;
                      cmm[i] = xm;
                      cv[i] = x;
                    }
                    
                    // orthogonalize the extrapolated value
                    if (ultrasoft) {
                      tmap["usfns"].start();
                      s_.wf.sd(ispin,ikp)->update_usfns();
                      tmap["usfns"].stop();
                    }
                    tmap["gram"].start();
                    s_.wf.sd(ispin,ikp)->gram();
                    tmap["gram"].stop();
                    //tmap["lowdin"].start();
                    //s_.wf.sd(ispin,ikp)->lowdin();
                    //tmap["lowdin"].stop();
                  
                    // c[i] now contains the extrapolated value
                  }
                  // c[i] is now ready for electronic iterations
                }
                else // normal extrapolation
                {
                  double* c = (double*) s_.wf.sd(ispin,ikp)->c().cvalptr();
                  double* cv = (double*) s_.wfv->sd(ispin,ikp)->c().cvalptr();
                  const int mloc = s_.wf.sd(ispin,ikp)->c().mloc();
                  const int nloc = s_.wf.sd(ispin,ikp)->c().nloc();
                  const int len = 2*mloc*nloc;
                  if ( iter == 0 )
                  {
                    // copy c to cv
                    for ( int i = 0; i < len; i++ )
                    {
                      const double x = c[i];
                      const double v = cv[i];
                      c[i] = x + dt * v;
                      cv[i] = x;
                    }
                    //tmap["lowdin"].start();
                    //s_.wf.sd(ispin,ikp)->lowdin();
                    //tmap["lowdin"].stop();
                    if (ultrasoft) {
                      tmap["usfns"].start();
                      s_.wf.sd(ispin,ikp)->update_usfns();
                      tmap["usfns"].stop();
                    }
                    tmap["gram"].start();
                    s_.wf.sd(ispin,ikp)->gram();
                    tmap["gram"].stop();
                  }
                  else
                  {
                    tmap["align"].start();
                    s_.wfv->align(s_.wf);
                    tmap["align"].stop();
                  
                    // linear extrapolation
                    for ( int i = 0; i < len; i++ )
                    {
                      const double x = c[i];
                      const double xm = cv[i];
                      c[i] = 2.0 * x - xm;
                      cv[i] = x;
                    }
                    //tmap["ortho_align"].start();
                    //s_.wf.sd(ispin,ikp)->ortho_align(*s_.wfv->sd(ispin,ikp));
                    //tmap["ortho_align"].stop();
                    
                    //tmap["riccati"].start();
                    //s_.wf.sd(ispin,ikp)->riccati(*s_.wfv->sd(ispin,ikp));
                    //tmap["riccati"].stop();
                    
                    if (ultrasoft) {
                      tmap["usfns"].start();
                      s_.wf.sd(ispin,ikp)->update_usfns();
                      tmap["usfns"].stop();
                    }
                    //ewd: lowdin doesn't yet work correctly with ultrasoft
                    //tmap["lowdin"].start();
                    //s_.wf.sd(ispin,ikp)->lowdin();
                    //tmap["lowdin"].stop();
                    tmap["gram"].start();
                    s_.wf.sd(ispin,ikp)->gram();
                    tmap["gram"].stop();
                  }
                }
              }
            }
          }
        }
      } // atoms_move && extrapolate_wf

      // do nitscf self-consistent iterations, each with nite electronic steps
      if ( wf_stepper != 0 )
      {
        wf_stepper->preprocess();
        if ( anderson_charge_mixing )
          for (int ispin = 0; ispin < nspin; ispin++)
            //if (s_.wf.spinactive(ispin))
              mixer[ispin]->restart();

        SimpleConvergenceDetector conv_scf(s_.ctrl.threshold_scf_nsteps, s_.ctrl.threshold_scf);
        bool convflag = false;

#ifdef HPM  
        HPM_Start("scfloop");
#endif
#ifdef TAU
        QB_Pstart(14,scfloop);
#endif
        // SCF LOOP
        for ( int itscf = 0; itscf < nitscf_; itscf++ )
        {

          // check timing
          if (niter <= 1 && s_.ctrl.run_timer > 0.0 && nitscf_ > 1 && itscf > 1 && testtimer) {
            double tnow = MPI_Wtime();
            double sofar = tnow - tbase;
            double tleft = s_.ctrl.run_timer - ( tnow - s_.ctrl.time_init); 
            double tscf = sofar/itscf; // avg. time per iteration
            double tmaxiter;
            MPI_Allreduce(&tscf, &tmaxiter, 1, MPI_DOUBLE, MPI_MAX, s_.ctxt_.comm());
            double maxtime = tmaxiter*(nitscf_-itscf);
            if (maxtime > tleft) { // we'll exceed timer, lower nitscf_
              s_.ctrl.timer_hit = true;
              int newiter = (int)(tleft/tmaxiter);
              int tmpnitscf = itscf + newiter - 1;
              if (tmpnitscf != nitscf_) {
                nitscf_ = tmpnitscf;
                if ( s_.ctxt_.oncoutpe() )
                  cout << "<!-- estimated scf iteration time = " << setprecision(3) << tmaxiter << " sec will exceed run_timer, changing nscf to " << nitscf_ << " -->" << endl;
              }
              else {
                testtimer = false;
                if ( s_.ctxt_.oncoutpe() )
                  cout << "<!-- estimated scf iteration time = " << setprecision(3) << tmaxiter << " sec is stable, using nscf = " << nitscf_ << " -->" << endl;
              }
            }
          }
          
          // check convergence
          if (conv_scf.isConverged()) {
            if ( s_.ctxt_.oncoutpe() ) {
              cout.setf(ios::scientific,ios::floatfield);
              cout << "  <!-- BOSampleStepper: scf convergence at itscf = " << itscf << ", energy varied by less than " << setprecision(2) 
                   << conv_scf.threshold() << " a.u. over " << conv_scf.nsteps() 
                   << " scf steps. -->" << endl;
            }
            itscf = nitscf_;
            convflag = true;
          }          
          // continue itscf loop
          else {
            if (itscf > 0)
              conv_scf.addValue(ef_.etotal());
            
            if ( nite_ > 1 && oncoutpe )
              cout << "  <!-- BOSampleStepper: start scf iteration -->" << endl;

            // compute new density in cd_.rhog
            tmap["charge"].start();
            if ( itscf==0 && initial_atomic_density ) {
              cd_.update_rhor();
              // AS: TODO: double check what to use here
              if (tddft_involved_)
                 ( ef_.hamil_cd() )->update_rhor();
              if (ultrasoft)
                assert(false);  // ewd DEBUG: need to implement this for ultrasoft
            }
            else
              cd_.update_density();
            tmap["charge"].stop();

            if ( charge_mixing != "off" && nite_ > 1) {
               if ( itscf == 0) {
                  //ewd:  read rhog_in from checkpoint if possible
                  for ( int ispin = 0; ispin < nspin; ispin++ ) {
                     int rhogflag = 1;

                     if (cell_dyn == "LOCKED" && (atoms_dyn == "LOCKED" || dt == 0.0) && s_.rhog_last.size() == rhog_in[ispin].size()) {
                        if (s_.rhog_last[ispin].size() == rhog_in[ispin].size()) {
                           rhogflag = 0;
                           for ( int i=0; i < rhog_in[ispin].size(); i++ ) 
                              rhog_in[ispin][i] = s_.rhog_last[ispin][i];
                        }
                     }
                  
                     if (rhogflag) {
                        for ( int i=0; i < rhog_in[ispin].size(); i++ )
                        {
                           if (tddft_involved_)
                              rhog_in[ispin][i] = ( ef_.hamil_cd() )->rhog[ispin][i];
                           else
                              rhog_in[ispin][i] = cd_.rhog[ispin][i];
                        }
                     }
                  }
               }
               // itscf > 0
               else {

                  for ( int ispin = 0; ispin < nspin; ispin++ ) {
                     // compute correction drhog
                     for ( int i=0; i < rhog_in[ispin].size(); i++ )
                     {
                        drhog[ispin][i] = (cd_.rhog[ispin][i] - rhog_in[ispin][i]);
                     }
                     
                     const double alpha = s_.ctrl.charge_mix_coeff;
                     // Anderson acceleration
                     if ( anderson_charge_mixing )
                     {
                        // row weighting of LS calculation
                        for ( int i=0; i < drhog[ispin].size(); i++ )
                           drhog[ispin][i] /= wls[i];
                    
                        //ewd:  try running this on every task
                        mixer[ispin]->update((double*)&rhog_in[ispin][0],(double*)&drhog[ispin][0],
                                             (double*)&rhobar[ispin][0],(double*)&drhobar[ispin][0]);

                        for ( int i=0; i < drhobar[ispin].size(); i++ )
                           drhobar[ispin][i] *= wls[i];
                    
                        for ( int i=0; i < rhog_in[ispin].size(); i++ )
                           rhog_in[ispin][i] = rhobar[ispin][i] + alpha * drhobar[ispin][i] * wkerker[i];
                     }
                     else {
                        for ( int i=0; i < rhog_in[ispin].size(); i++ )
                           rhog_in[ispin][i] += alpha * drhog[ispin][i] * wkerker[i];
                     }

                     if (tddft_involved_)
                     {
                        for ( int i=0; i < rhog_in[ispin].size(); i++ )
                           ( ef_.hamil_cd() )->rhog[ispin][i] = rhog_in[ispin][i];
                        ( ef_.hamil_cd() )->update_rhor();                     
                     }
                     
                     // Apply correction
                     for ( int i=0; i < rhog_in[ispin].size(); i++ )
                        cd_.rhog[ispin][i] = rhog_in[ispin][i];
                  }
                  cd_.update_rhor();
               }              
            }

            // AS: update the Hamiltonian, the potential, and the energy before propagation
            // starts and/or after the mixing
            if (tddft_involved_)
               ef_.update_hamiltonian();
              
            ef_.update_vhxc();
            //ef_.update_exc_ehart_eps();

            // AS: the first step is an EULER one in the case of the second-order propagation
            if ( ( wf_dyn == "SOTD" ) && wfv_is_new ) {
               // AS: |psi(t-tddt)> is the wave function in the very beginning
               *wfdeque[0]=wf;

               // AS: |psi(t)> is calculated by doing one TDEULER step using the EULER stepper
               WavefunctionStepper* wf_init_stepper = new TDEULERWavefunctionStepper(wf,s_.ctrl.tddt,tmap);
               
               cout << wf_dyn << " initialization energy: " << ef_.energy(true,dwf,false,fion,false,sigma_eks) << endl;
               cout << wf_dyn << " initialization expectation value: " << s_.wf.dot(dwf) << endl;

               wf_init_stepper->update(dwf);
               // AS: now we have |psi(t)>
               *wfdeque[1]=wf;
               (*s_.hamil_wf)=s_.wf;

               // AS: the wave functions used in the Hamiltonian are NOT updated here
               cd_.update_density();
               ( ef_.hamil_cd() )->update_density();

               // AS: update the Hamiltonian after the first EULER step
               ef_.update_hamiltonian();
               ef_.update_vhxc();
               // AS: modified version of EnergyFunctional::update_vhxc(void) which leaves the potential
               // untouched and only recalculates the energy terms
               //ef_.update_exc_ehart_eps();   
               // AS: we are ready to go with the SOTD
               if ( wf_dyn == "SOTD" ) wfv_is_new=false;
            }
          
            if (wf_dyn=="FORKTD") {
               // AS: initialize wfdeque[0] with 0, it will be k_1
               *wfdeque[0]=dwf;
               // AS: initialize wfdeque[1] with 0, it will be k_2
               *wfdeque[1]=s_.wf;
               // AS: initialize wfdeque[2] with 0, it will be k_3
               *wfdeque[2]=s_.wf;
               // AS: initialize wfdeque[3] with 0, it will be k_4
               *wfdeque[3]=s_.wf;
               // AS: initialize wfdeque[4] with wf, to keep a copy
               *wfdeque[4]=s_.wf;
            }
            
            // reset stepper only if multiple non-selfconsistent steps
            if ( nite_ > 1 ) wf_stepper->preprocess();
            double lastnonscfetot;
            const double nonscfthresh = s_.ctrl.threshold_nonscf;
            for ( int ite = 0; ite < nite_; ite++ )
            {
              double energy = ef_.energy(true,dwf,false,fion,false,sigma_eks);
              double delta_etotnonscf = abs(energy-lastnonscfetot);
              if (ite > 0 && delta_etotnonscf < nonscfthresh) {
                if ( s_.ctxt_.oncoutpe() ) {
                  cout.setf(ios::fixed,ios::floatfield);
                  cout.setf(ios::right,ios::adjustfield);
                  cout << "  <etotal_int> " << setw(15) << setprecision(8)
                       << energy << " </etotal_int>\n";
                  if ( compute_stress ) {
                    const double pext = (sigma_ext[0]+sigma_ext[1]+sigma_ext[2])/3.0;
                    const double enthalpy = energy + pext * cell.volume();
                    cout << "  <enthalpy_int> " << setw(15) 
                         << enthalpy << " </enthalpy_int>\n"
                         << flush;
                  }
                  cout.setf(ios::scientific,ios::floatfield);
                  cout << "  <!-- BOSampleStepper: non-scf threshold " << setprecision(2) << nonscfthresh << " reached. -->" << endl;
                }
                ite = nite_;
              }
              else {
                lastnonscfetot = energy;
                // compute the sum of eigenvalues (with fixed weight)
                // to measure convergence of the subspace update
                // compute trace of the Hamiltonian matrix Y^T H Y
                // scalar product of Y and (HY): tr Y^T (HY) = sum_ij Y_ij (HY)_ij
                // Note: since the hamiltonian is hermitian and dwf=H*wf
                // the dot product in the following line is real

                if (ultrasoft) { 
                  const double us_eigenvalue_sum = s_.wf.sdot(dwf);
                  if ( oncoutpe )
                    cout  << "  <eigenvalue_sum> "
                          << us_eigenvalue_sum << " </eigenvalue_sum>" << endl;
                }
                else {
                  const double eigenvalue_sum = s_.wf.dot(dwf);
                  if ( oncoutpe )
                    cout << "  <eigenvalue_sum> "
                         << eigenvalue_sum << " </eigenvalue_sum>" << endl;
                }



                // AS: calculation and output of < psi(t) | H(t_0) | psi(t) > , i.e., the expectation values
                // AS: in the case the wave functions are propagated in time
                // AS: code taken from Wavefunction::diag(Wavefunction& dwf, bool eigvec)
                if ( tddft_involved_ && s_.ctrl.wf_diag == "EIGVAL" )
                {
                   if ( oncoutpe ) cout << "<" << wf_dyn << " expectation set>" << endl;

                   // AS: this is the new version, not copying code, but copying the wave function instead
                   Wavefunction* to_diag_wf1 = new Wavefunction(s_.wf);
                   (*to_diag_wf1) = (s_.wf);
                   (*to_diag_wf1).diag(dwf,false);

                   if ( oncoutpe )
                   {
                      for ( int ispin = 0; ispin < (*to_diag_wf1).nspin(); ispin++ )
                      {
                         for ( int ikp = 0; ikp < (*to_diag_wf1).nkp(); ikp++ )
                         {
                            const int nst = (*to_diag_wf1).sd(ispin,ikp)->nst();
                            const double eVolt = 2.0 * 13.6058;
                            cout <<    "  <" << wf_dyn << " expectation values spin=\"" << ispin
                                 << "\" kpoint=\"" << s_.wf.sd(ispin,ikp)->kpoint()
                                 << "\" n=\"" << nst << "\">" << endl;
                            cout << "  ";
                            for ( int i = 0; i < nst; i++ )
                            {
                               cout << setw(15) << setprecision(8)
                                    << (*to_diag_wf1).sd(ispin,ikp)->eig(i)*eVolt;
                               if ( i%5 == 4 ) cout << endl;
                            }
                            if ( nst%5 != 0 ) cout << endl;
                            cout << "  </" << wf_dyn << " expectation values>" << endl;
                         }
                      }
                   }

                   if ( oncoutpe ) cout << "</" << wf_dyn << " expectation set>" << endl;
                   delete(to_diag_wf1);
                }

                // AS: calculation and output of < psi(t) | psi(t) > , i.e., the orthonormalization
                // AS: code adopted from SlaterDet::gram(void)
                if ( tddft_involved_ && s_.ctrl.wf_diag == "EIGVAL" )
                {
                   for ( int ispin = 0; ispin < (wf).nspin(); ispin++ )
                   {
                      for ( int ikp = 0; ikp < (wf).nkp(); ikp++ )
                      {
                         ComplexMatrix ortho(wf.sd(ispin,ikp)->context(),(wf.sd(ispin,ikp)->c()).n(),(wf.sd(ispin,ikp)->c()).n(),(wf.sd(ispin,ikp)->c()).nb(),(wf.sd(ispin,ikp)->c()).nb());
                                                                                                      
                         ortho.gemm('c','n',1.0,(wf).sd(ispin,ikp)->c(),(wf).sd(ispin,ikp)->c(),0.0);
                         if ( oncoutpe )
                         {
                            cout << "ortho: " << endl;
                         }
                         //ewd:  this is not going to print correctly, as tasks will not print data in order
                         cout << ortho;   
                      }
                   }
                }

                if (wf_dyn=="SORKTD") {
                   wfdeque.push_back( s_.wfv );
                   // AS: initialize wfdeque[0] with 0, it will be k_1
                   (*wfdeque[0]).clear();

                   wfdeque.push_back( s_.wfv );
                   // AS: initialize wfdeque[1] with 0, it will be k_2
                   (*wfdeque[1]).clear();

                   wfdeque.push_back( s_.wfv );
                   // AS: initialize wfdeque[2] with wf, to keep a copy
                   (*wfdeque[2])=s_.wf;

                   for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
                   {
                      for ( int ikp = 0; ikp < wf.nkp(); ikp++ )
                      {
                         // AS: here we turn wfdeque[0] into k_1
                         ((*wfdeque[0]).sd(ispin,ikp)->c()).axpy(-1.0*s_.ctrl.tddt*(complex<double>(0,1)), dwf.sd(ispin,ikp)->c() );
                         // AS: here we turn wf into |psi(t)+0.5*k_1>
                         ((s_.wf).sd(ispin,ikp)->c()).axpy(0.5, (*wfdeque[0]).sd(ispin,ikp)->c() );
                      }
                   }

                   // AS: for the correct output of the nonscf energy
                   if (s_.ctrl.non_selfc_energy)
                   {
                      // AS: the wave functions used in the Hamiltonian are NOT updated here
                      cd_.update_density();
                      // AS: modified version of EnergyFunctional::update_vhxc(void) which leaves the
                      // potential untouched and only recalculates the energy terms
                      //ef_.update_exc_ehart_eps();
                   }
                   
                   // AS: update the Hamiltonian with psi(t)+0.5*k_1
                   (*s_.hamil_wf)=s_.wf;
                   ( ef_.hamil_cd() )->update_density();
                   ef_.update_hamiltonian();
                   ef_.update_vhxc();

                   // AS: apply the Hamiltonian to |psi(t)+0.5*k_1>
                   ef_.energy(true,dwf,false,fion,false,sigma_eks);

                   for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
                   {
                      for ( int ikp = 0; ikp < wf.nkp(); ikp++ )
                      {
                         // AS: here we turn wfdeque[1] into k_2
                         ((*wfdeque[1]).sd(ispin,ikp)->c()).axpy(-1.0*s_.ctrl.tddt*(complex<double>(0,1)), dwf.sd(ispin,ikp)->c() );
                      }
                   }
                }

                if (wf_dyn=="FORKTD") {
                   //ewd DEBUG
                   tmap["forktd_tot"].start();
                   tmap["forktd_push"].start();

                   // AS: initialize wfdeque[0] with 0, it will be k_1
                   (*wfdeque[0]).clear();
                   
                   // AS: initialize wfdeque[1] with 0, it will be k_2
                   (*wfdeque[1]).clear();

                   // AS: initialize wfdeque[2] with 0, it will be k_3
                   (*wfdeque[2]).clear();

                   // AS: initialize wfdeque[3] with 0, it will be k_4
                   (*wfdeque[3]).clear();

                   // AS: initialize wfdeque[4] with wf, to keep a copy
                   (*wfdeque[4])=s_.wf;

                   //ewd DEBUG
                   tmap["forktd_push"].stop();
                   tmap["forktd_axpy"].start();
                   //ewd DEBUG

                   for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
                   {
                      for ( int ikp = 0; ikp < wf.nkp(); ikp++ )
                      {
                         // AS: here we turn wfdeque[0] into k_1
                         ((*wfdeque[0]).sd(ispin,ikp)->c()).axpy(-1.0*s_.ctrl.tddt*(complex<double>(0,1)), dwf.sd(ispin,ikp)->c() );
                         // AS: here we turn wf into |psi(t)+0.5*k_1>
                         ((s_.wf).sd(ispin,ikp)->c()).axpy(0.5, (*wfdeque[0]).sd(ispin,ikp)->c() );
                      }
                   }

                   //ewd DEBUG
                   tmap["forktd_axpy"].stop();
                   //ewd DEBUG

                   // AS: for the correct output of the nonscf energy
                   if (s_.ctrl.non_selfc_energy)
                   {
                      // AS: the wave functions used in the Hamiltonian are NOT updated here
                      cd_.update_density();

                      // AS: modified version of EnergyFunctional::update_vhxc(void) which leaves the
                      // potential untouched and only recalculates the energy terms
                      //ef_.update_exc_ehart_eps();
                   }

                   //ewd DEBUG
                   tmap["forktd_wfcp"].start();
                   //ewd DEBUG

                   // AS: set the Hamiltonian accordingly
                   (*s_.hamil_wf)=s_.wf;

                   //ewd DEBUG
                   tmap["forktd_wfcp"].stop();
                   tmap["forktd_cd"].start();
                   //ewd DEBUG

                   ( ef_.hamil_cd() )->update_density();

                   //ewd DEBUG
                   tmap["forktd_cd"].stop();
                   tmap["forktd_efup"].start();
                   //ewd DEBUG

                   ef_.update_hamiltonian();
                   ef_.update_vhxc();

                   //ewd DEBUG
                   tmap["forktd_efup"].stop();
                   tmap["forktd_efen"].start();
                   //ewd DEBUG

                   // AS: apply the Hamiltonian to |psi(t)+0.5*k_1>
                   ef_.energy(true,dwf,false,fion,false,sigma_eks);

                   //ewd DEBUG
                   tmap["forktd_efen"].stop();
                   tmap["forktd_wfcp"].start();
                   //ewd DEBUG

                   // AS: setting wf back to |psi(t)>
                   s_.wf=(*wfdeque[4]);

                   //ewd DEBUG
                   tmap["forktd_wfcp"].stop();
                   tmap["forktd_axpy"].start();
                   //ewd DEBUG                                            
                   for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
                   {
                      for ( int ikp = 0; ikp < wf.nkp(); ikp++ )
                      {
                         // AS: here we turn wfdeque[1] into k_2
                         ((*wfdeque[1]).sd(ispin,ikp)->c()).axpy(-1.0*s_.ctrl.tddt*(complex<double>(0,1)), dwf.sd(ispin,ikp)->c() );
                         // AS: here we turn wf into |psi(t)+0.5*k_2>
                         ((s_.wf).sd(ispin,ikp)->c()).axpy(0.5, (*wfdeque[1]).sd(ispin,ikp)->c() );
                      }
                   }


                   //ewd DEBUG
                   tmap["forktd_axpy"].stop();
                   //ewd DEBUG

                   // AS: for the correct output of the energy
                   if (s_.ctrl.non_selfc_energy)
                   {
                      // AS: the wave functions used in the Hamiltonian are NOT updated here
                      cd_.update_density();

                      // AS: modified version of EnergyFunctional::update_vhxc(void) which leaves the
                      // potential untouched and only recalculates the energy terms
                      //ef_.update_exc_ehart_eps();
                   }
                   
                   // AS: set the Hamiltonian accordingly

                   //ewd DEBUG
                   tmap["forktd_wfcp"].start();
                   //ewd DEBUG

                   (*s_.hamil_wf)=s_.wf;

                   //ewd DEBUG
                   tmap["forktd_wfcp"].stop();
                   tmap["forktd_cd"].start();
                   //ewd DEBUG

                   ( ef_.hamil_cd() )->update_density();

                   //ewd DEBUG
                   tmap["forktd_cd"].stop();
                   tmap["forktd_efup"].start();
                   //ewd DEBUG

                   ef_.update_hamiltonian();
                   ef_.update_vhxc();

                   //ewd DEBUG
                   tmap["forktd_efup"].stop();
                   tmap["forktd_efen"].start();
                   //ewd DEBUG

                   // AS: apply the Hamiltonian to |psi(t)+0.5*k_2>
                   ef_.energy(true,dwf,false,fion,false,sigma_eks);

                   //ewd DEBUG
                   tmap["forktd_efen"].stop();
                   tmap["forktd_wfcp"].start();
                   //ewd DEBUG

                   // AS: setting wf back to |psi(t)>
                   s_.wf=(*wfdeque[4]);

                   tmap["forktd_wfcp"].stop();
                   tmap["forktd_axpy"].start();
                   //ewd DEBUG

                   for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
                   {
                      for ( int ikp = 0; ikp < wf.nkp(); ikp++ )
                      {
                         // AS: here we turn wfdeque[2] into k_3
                         ((*wfdeque[2]).sd(ispin,ikp)->c()).axpy(-1.0*s_.ctrl.tddt*(complex<double>(0,1)), dwf.sd(ispin,ikp)->c() );
                         // AS: here we turn wf into |psi(t)+k_3>
                         ((s_.wf).sd(ispin,ikp)->c()).axpy(1.0, (*wfdeque[2]).sd(ispin,ikp)->c() );
                      }
                   }
                   //ewd DEBUG
                   tmap["forktd_axpy"].stop();
                   //ewd DEBUG

                   // AS: for the correct output of the energy
                   if (s_.ctrl.non_selfc_energy)
                   {
                      // AS: the wave functions used in the Hamiltonian are NOT updated here
                      cd_.update_density();

                      // AS: modified version of EnergyFunctional::update_vhxc(void) which leaves the
                      // potential untouched and only recalculates the energy terms
                      //ef_.update_exc_ehart_eps();
                   }

                   // AS: set the Hamiltonian accordingly

                   //ewd DEBUG
                   tmap["forktd_wfcp"].start();
                   //ewd DEBUG

                   (*s_.hamil_wf)=s_.wf;

                   //ewd DEBUG
                   tmap["forktd_wfcp"].stop();
                   tmap["forktd_cd"].start();
                   //ewd DEBUG

                   ( ef_.hamil_cd() )->update_density();

                   //ewd DEBUG
                   tmap["forktd_cd"].stop();
                   tmap["forktd_efup"].start();
                   //ewd DEBUG

                   ef_.update_hamiltonian();
                   ef_.update_vhxc();

                   //ewd DEBUG
                   tmap["forktd_efup"].stop();
                   tmap["forktd_efen"].start();
                   //ewd DEBUG

                   // AS: apply the Hamiltonian to |psi(t)+k_3>
                   ef_.energy(true,dwf,false,fion,false,sigma_eks);

                   //ewd DEBUG
                   tmap["forktd_efen"].stop();
                   tmap["forktd_axpy"].start();
                   //ewd DEBUG

                   for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
                   {
                      for ( int ikp = 0; ikp < wf.nkp(); ikp++ )
                      {
                         // AS: here we turn wfdeque[3] into k_4
                         ((*wfdeque[3]).sd(ispin,ikp)->c()).axpy(-1.0*s_.ctrl.tddt*(complex<double>(0,1)), dwf.sd(ispin,ikp)->c() );
                      }
                   }

                   //ewd DEBUG
                   tmap["forktd_axpy"].stop();
                   tmap["forktd_tot"].stop();
                   //ewd DEBUG
                }
                // AS: energy renormalization is currently disabled
                // output of the renormalized energy
                // double wf_dyn_energy;
                // double wf_dyn_eigenvalue_sum;
                // if ( tddft_involved ) {
                // wf_dyn_energy = ef_.energy(true,dwf,false,fion,false,sigma_eks,true);
                // wf_dyn_eigenvalue_sum = real(s_.wf.dot(dwf));
                // if ( oncoutpe )
                // {
                //   cout << wf_dyn << " energy: " << wf_dyn_energy << endl;
                //   cout << wf_dyn << " expectation value: " << wf_dyn_eigenvalue_sum << endl;
                // }
                // }    
                
                wf_stepper->update(dwf);
                if (ultrasoft)
                   wf.update_usfns();
                
                // update ultrasoft functions if needed, call gram
                for ( int ispin = 0; ispin < s_.wf.nspin(); ispin++ ) {
                  if (s_.wf.spinactive(ispin)) {
                    for ( int ikp = 0; ikp < s_.wf.nkp(); ikp++ ) {
                      if (s_.wf.kptactive(ikp)) {
                        assert(s_.wf.sd(ispin,ikp) != 0);
                        //if (ultrasoft) { 
                        //  tmap["usfns"].start();
                        //  s_.wf.sd(ispin,ikp)->update_usfns(); // calculate betapsi, spsi
                        //  tmap["usfns"].stop();
                        //}
                        tmap["gram"].start();
                        s_.wf.sd(ispin,ikp)->gram();
                        tmap["gram"].stop();

                        //if (ultrasoft) { 
                        //  tmap["usfns"].start();
                        //  s_.wf.sd(ispin,ikp)->calc_betapsi(); // calculate betapsi
                        //  tmap["usfns"].stop();
                        //}
                      }
                    }
                  }
                }

                // AS: change the phase of the wave function if the respective variable is set
                if ( s_.wf.phase_real_set() ) {
                   s_.wf.phase_wf_real();

                   //ewd DEBUG
                   //if (oncoutpe)
                   //cout << "DEBUG:  PHASE_WF_REAL SET TO TRUE!" << endl;
                   //ewd DEBUG
  
                }
                
                
                if ( oncoutpe )
                {
                  cout.setf(ios::fixed,ios::floatfield);
                  cout.setf(ios::right,ios::adjustfield);
                  cout << "  <etotal_int> " << setw(15) << setprecision(8)
                       << energy << " </etotal_int>\n";
                  if ( compute_stress )
                  {
                    const double pext = (sigma_ext[0]+sigma_ext[1]+sigma_ext[2])/3.0;
                    const double enthalpy = energy + pext * cell.volume();
                    cout << "  <enthalpy_int> " << setw(15)
                         << enthalpy << " </enthalpy_int>\n"
                         << flush;
                  }                  
                }
                
                if ( ( nite_ > 1 ) && ( (wf_dyn=="SORKTD") || (wf_dyn=="FORKTD") ) ) {
                   (*s_.hamil_wf)=s_.wf;
                   ( ef_.hamil_cd() )->update_density();
                   cd_.update_density();
                   ef_.update_hamiltonian();
                   ef_.update_vhxc();
                   // AS: this should be only needed if cd_ has changed and no update_vhxc call is
                   //     done and the energy is written somewhere after
                   // ef_.update_exc_ehart_eps();
                }

                // AS: for the correct output of the energy
                if ( s_.ctrl.non_selfc_energy && (wf_dyn!="SORKTD") && (wf_dyn!="FORKTD") )
                {

                   //ewd DEBUG
                   //if (oncoutpe)
                   //   cout << "DEBUG:  NON_SELFC_ENERGY SET TO TRUE!" << endl;
                   //ewd DEBUG

                   // AS: the wave functions used in the Hamiltonian are NOT updated here
                   cd_.update_density();

                   // AS: modified version of EnergyFunctional::update_vhxc(void) which leaves the
                   // potential untouched and only
                   // AS: recalculates the energy terms
                   ef_.update_exc_ehart_eps();
                }
              }
            } // for ite
            // subspace diagonalization
            if ( compute_eigvec || s_.ctrl.wf_diag == "EIGVAL" || usdiag && (!tddft_involved_))
            {
              ef_.energy(true,dwf,false,fion,false,sigma_eks);
              tmap["diag"].start();
              s_.wf.diag(dwf,compute_eigvec);
              tmap["diag"].stop();

              // update ultrasoft functions w. new eigenvectors
              if (compute_eigvec && ultrasoft)
                 s_.wf.update_usfns();

              if (itscf%s_.ctrl.iprint == 0)
                s_.wf.printeig();
            }

            // update occupation numbers if fractionally occupied states
            if ( fractional_occ )
            {
              wf.update_occ(s_.ctrl.smearing_width,s_.ctrl.smearing_ngauss);
              if (itscf%s_.ctrl.iprint == 0)
                s_.wf.printocc();
              const double wf_entropy = wf.entropy();
              if ( oncoutpe )
              {
                cout << "  <!-- Wavefunction entropy: " << wf_entropy << " -->" << endl;
                //const double boltz = 1.0 / ( 11605.0 * 2.0 * 13.6058 );
                cout << "  <!-- Entropy contribution to free energy: "
                     << - wf_entropy * s_.ctrl.smearing_width * 2.0 << " -->" << endl;
              }
            }

            // AS: after the first set of non-selfconsistent steps the Hamiltonian is no longer fixed
            // AS: this makes the Hamiltonian time dependent
            if (tddft_involved_)
            {
               (*s_.hamil_wf)=s_.wf;
               ( ef_.hamil_cd() )->update_density();
            }
            
            if ( nite_ > 1 && oncoutpe )
              cout << "  <!-- BOSampleStepper: end scf iteration -->" << endl;
          }
        } // for itscf

#ifdef TAU  
        QB_Pstop(scfloop);
#endif
#ifdef HPM
        HPM_Stop("scfloop");
#endif

        // AS: keep the previous wave function
        if ( wf_dyn == "SOTD" ) *(s_.wfv)=*wfdeque[0];

        // AS: non-adiabatic overlap for succeeding steps of the Born-Oppenheimer trajectory is calculated
        if ( ( atoms_move ) && ( s_.ctrl.na_overlap_min >= 0 ) )
        {
           // AS: after the first step previous_wf is simply set to the current wf
           if ( s_.previous_wf == 0 )
           {
              s_.previous_wf = new Wavefunction(wf);
              (*s_.previous_wf) = s_.wf;
           }
           
           // AS: ATTENTION: currently the min/max indices are not used yet! Overlap is calculated for all the states!
           for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
           {
              for ( int ikp = 0; ikp < wf.nkp(); ikp++ )
              {
                 if ( wf.force_complex_set() ) {
                    ComplexMatrix overlap1(wf.sd(ispin,ikp)->context(),(wf.sd(ispin,ikp)->c()).n(),(wf.sd(ispin,ikp)->c()).n(),(wf.sd(ispin,ikp)->c()).nb(),(wf.sd(ispin,ikp)->c()).nb());
                    ComplexMatrix overlap2(wf.sd(ispin,ikp)->context(),(wf.sd(ispin,ikp)->c()).n(),(wf.sd(ispin,ikp)->c()).n(),(wf.sd(ispin,ikp)->c()).nb(),(wf.sd(ispin,ikp)->c()).nb());

                    // AS: The following two lines should give something like the d_jk
                    // overlap1.gemm('c','n',-1.0*(complex<double>(0,1))/2.0/s_.ctrl.dt,(*s_.previous_wf).sd(ispin,ikp)->c(),wf.sd(ispin,ikp)->c(),0.0);
                    // overlap1.gemm('c','n',1.0*(complex<double>(0,1))/2.0/s_.ctrl.dt,wf.sd(ispin,ikp)->c(),(*s_.previous_wf).sd(ispin,ikp)->c(),1.0);

                    // AS: The following two lines are the same stuff without the complex i
                    overlap1.gemm('c','n',-1.0/2.0/s_.ctrl.dt,(*s_.previous_wf).sd(ispin,ikp)->c(),wf.sd(ispin,ikp)->c(),0.0);
                    overlap2.gemm('c','n',1.0/2.0/s_.ctrl.dt,wf.sd(ispin,ikp)->c(),(*s_.previous_wf).sd(ispin,ikp)->c(),0.0);

                    // AS: this is the normalization
                    // overlap1.gemm('c','n',1.0,wf.sd(ispin,ikp)->c(),wf.sd(ispin,ikp)->c(),0.0);
                    
                    cout << "OVERLAP complex: " << endl;
                    cout << "overlap1" << overlap1 << endl;
                    cout << "overlap2" << overlap2 << endl;
                    overlap1 += overlap2;
                    cout << "overlap1+overlap2" << overlap1 << endl;
                 }
                 else
                 {
                    DoubleMatrix overlap1(wf.sd(ispin,ikp)->context(),(wf.sd(ispin,ikp)->c()).n(),(wf.sd(ispin,ikp)->c()).n(),(wf.sd(ispin,ikp)->c()).nb(),(wf.sd(ispin,ikp)->c()).nb());
                    DoubleMatrix overlap2(wf.sd(ispin,ikp)->context(),(wf.sd(ispin,ikp)->c()).n(),(wf.sd(ispin,ikp)->c()).n(),(wf.sd(ispin,ikp)->c()).nb(),(wf.sd(ispin,ikp)->c()).nb());

                    // DoubleMatrix proxies
                    DoubleMatrix previous_wf_proxy( (*s_.previous_wf).sd(ispin,ikp)->c() );
                    DoubleMatrix wf_proxy( wf.sd(ispin,ikp)->c() );

                    // AS: The following two lines are something like the d_jk without the complex i
                    // AS: 2*(-1.0/2.0/s_.ctrl.dt) is the factor
                    overlap1.gemm('c','n',-1.0/s_.ctrl.dt,previous_wf_proxy,wf_proxy,0.0);
                    // AS: -1.0*(-1.0/2.0/s_.ctrl.dt) is the factor
                    overlap1.ger(1.0/2.0/s_.ctrl.dt,previous_wf_proxy,0,wf_proxy,0);
                    overlap2.gemm('c','n',1.0/s_.ctrl.dt,wf_proxy,previous_wf_proxy,0.0);
                    overlap2.ger(-1.0/2.0/s_.ctrl.dt,wf_proxy,0,previous_wf_proxy,0);

                    // AS: the following three lines reproduce the correct normalization
                    // overlap1.gemm('c','n',2.0,wf_proxy,wf_proxy,0.0);
                    // rank-1 update using first row of sdc_proxy() and c_proxy
                    // overlap1.ger(-1.0,wf_proxy,0,wf_proxy,0);

                    cout << "OVERLAP real: " << endl;
                    cout << "overlap1" << overlap1 << endl;
                    cout << "overlap2" << overlap2 << endl;
                    overlap1 += overlap2;
                    cout << "overlap1+overlap2" << overlap1 << endl;
                 }
              }
           }

           // AS: keep a copy of the wave function to calculate the overlap during the following iteration
           delete s_.previous_wf;
           s_.previous_wf = new Wavefunction(s_.wf);
           (*s_.previous_wf) = s_.wf;
        }

        // AS: this implements a function to print the density every N number of MD steps
        if ( (s_.ctrl.print_density_every > 0) && ((iter % s_.ctrl.print_density_every) == 0) ) {
           PlotCmd plot_density_from_stepper(&s_);
           std::ostringstream oss; oss << std::setfill('0') << std::setw(6) << iter;
           string arg1 = "plot", arg2 = "-density", arg3 = "plots/density-" + oss.str() + ".cub";
           char* argv[] = {const_cast<char*>(arg1.c_str()), const_cast<char*>(arg2.c_str()), const_cast<char*>(arg3.c_str())};
           plot_density_from_stepper.action(3,argv);
        }

            
     

  
        if (!convflag && s_.ctrl.threshold_scf > 0.0) 
          if ( s_.ctxt_.oncoutpe() ) 
            cout << "<WARNING> Ionic iteration finished without triggering scf convergence threshold, consider increasing nscf = " << nitscf_ << " </WARNING>" << endl;
                      
        if ( compute_mlwf || compute_mlwfc )
        {
          SlaterDet& sd = *(wf.sd(0,0));
          mlwft->compute_transform();

          if ( compute_mlwf )
            mlwft->apply_transform(sd);
          
          if ( oncoutpe )
          {
            cout << " <mlwf_set size=\"" << sd.nst() << "\">" << endl;
            for ( int i = 0; i < sd.nst(); i++ )
            {
              D3vector ctr = mlwft->center(i);
              double sp = mlwft->spread(i);
              cout.setf(ios::fixed, ios::floatfield);
              cout.setf(ios::right, ios::adjustfield);
              cout << "   <mlwf center=\"" << setprecision(6)
                   << setw(12) << ctr.x
                   << setw(12) << ctr.y
                   << setw(12) << ctr.z
                   << " \" spread=\" " << sp << " \"/>"
                   << endl;
            }
            cout << " </mlwf_set>" << endl;
            D3vector edipole = mlwft->dipole();
            cout << " <electronic_dipole> " << edipole
                 << " </electronic_dipole>" << endl;
            D3vector idipole = atoms.dipole();
            cout << " <ionic_dipole> " << idipole
                 << " </ionic_dipole>" << endl;
            cout << " <total_dipole> " << idipole + edipole
                 << " </total_dipole>" << endl;
            cout << " <total_dipole_length> " << length(idipole + edipole)
                 << " </total_dipole_length>" << endl;
          }
        }


        // If GS calculation only, print energy and atomset at end of iterations
        bool fastend = false;
#ifdef BGQ        
        fastend = true;
#endif        

        if ( gs_only && !fastend)
        {
           // need eigenvalues to compute forces w. ultrasoft
          if (ultrasoft) { 
            ef_.energy(true,dwf,false,fion,false,sigma_eks);
            tmap["diag"].start();
            s_.wf.diag(dwf,compute_eigvec);
            //s_.wf.diag(dwf,true);  // ewd:  why true?  Why did I do this??
            tmap["diag"].stop();

            // update ultrasoft functions w. new eigenvectors
            if (compute_eigvec && ultrasoft) { 
               tmap["usfns"].start();
               s_.wf.update_usfns();
               tmap["usfns"].stop();
            }
            s_.wf.printeig();
          }

           tmap["charge"].start();
           cd_.update_density();
           tmap["charge"].stop();

           if (tddft_involved_)
           {
              ( ef_.hamil_cd() )->update_density();
              ef_.update_hamiltonian();
           }
           
          ef_.update_vhxc();
          const bool compute_forces = true;
          ef_.energy(false,dwf,compute_forces,fion,compute_stress,sigma_eks);

          if (s_.atoms.add_fion_ext()) {
            for ( int is = 0; is < atoms.atom_list.size(); is++ ) {
              for ( int ia = 0; ia < atoms.atom_list[is].size(); ia++ ) {
                D3vector ftmp = s_.atoms.get_fion_ext(is,ia);
                fion[is][3*ia] += ftmp.x;
                fion[is][3*ia+1] += ftmp.y;
                fion[is][3*ia+2] += ftmp.z;
              }
            }
          }

          if ( oncoutpe )
          {
            cout << ef_;
            cout << "<atomset>" << endl;
            cout << atoms.cell();
            for ( int is = 0; is < atoms.atom_list.size(); is++ )
            {
              int i = 0;
              for ( int ia = 0; ia < atoms.atom_list[is].size(); ia++ )
              {
                Atom* pa = atoms.atom_list[is][ia];
                cout << "  <atom name=\"" << pa->name() << "\""
                     << " species=\"" << pa->species()
                     << "\">\n"
                     << "    <position> " << pa->position() << " </position>\n"
                     << "    <velocity> " << pa->velocity() << " </velocity>\n"
                     << "    <force> "
                     << fion[is][i] << " "
                     << fion[is][i+1] << " "
                     << fion[is][i+2]
                     << " </force>\n";
                cout << "  </atom>" << endl;
                i += 3;
              }
            }
            cout << "</atomset>" << endl;
            if ( compute_stress )
            {
              compute_sigma();
              print_stress();
            }
          }
        }
        else if (gs_only && fastend) {
           if ( oncoutpe )
              cout << ef_;
        }
        wf_stepper->postprocess();
      }
      else
      {
        // wf_stepper == 0, wf_dyn == LOCKED
        // evaluate and print energy
        tmap["charge"].start();
        cd_.update_density();
        tmap["charge"].stop();
        ef_.update_vhxc();
        ef_.energy(true,dwf,false,fion,false,sigma_eks);
        if ( oncoutpe )
        {
          cout << ef_;
        }
      }

#ifdef USE_APC
      ApcStop(1);
#endif
      // print iteration time
      double time = tm_iter.real();
      double tmin = time;
      double tmax = time;
      s_.ctxt_.dmin(1,1,&tmin,1);
      s_.ctxt_.dmax(1,1,&tmax,1);
      if ( oncoutpe )
      {
       cout << left << setw(34) << "<timing where=\"run\""
           << setw(24) << " name=\" iteration\""
             << " min=\"" << setprecision(3) << setw(9) << tmin << "\""
             << " max=\"" << setprecision(3) << setw(9) << tmax << "\"/>"
             << endl;
        cout << "</iteration>" << endl;
      }
      if ( atoms_move )
        s_.constraints.update_constraints(dt);
    } // else (not converged)
  } // for iter
  // print memory usage of main data objects
  if (s_.ctxt_.oncoutpe()) {
    double memtot = 0.0;
    double loctot = 0.0;

    int kmult = s_.wf.nspin()*s_.wf.nkp();
    int kmultloc = s_.wf.nkptloc();
    ef_.print_memory(cout,memtot,loctot);
    s_.wf.sd(0,0)->print_memory(cout,kmult,kmultloc,memtot,loctot);

    cd_.print_memory(cout,memtot,loctot);
    
    PrintMem pm;
    string totunit = pm.memunit(memtot);
    string locunit = pm.memunit(loctot);
    cout << "<!-- memory total       :  " << setw(7) << memtot << totunit << "  (" << loctot << locunit << " local) -->" << endl;
  }
  

#ifdef PRINTALL
  //ewd print timing
  cd_.print_timing();
  ef_.print_timing();
  // print timer map
  for ( TimerMap::iterator i = tmap.begin(); i != tmap.end(); i++ )
  {
     double time = (*i).second.real();
     double tmin = time;
     double tmax = time;
     s_.ctxt_.dmin(1,1,&tmin,1);
     s_.ctxt_.dmax(1,1,&tmax,1);
     if ( s_.ctxt_.mype()==0 )
     {
       cout << left << setw(34) << "<timing where=\"run\""
           << setw(8) << " name=\""
           << setw(15) << (*i).first << "\""
           << " min=\"" << setprecision(3) << setw(9) << tmin << "\""
           << " max=\"" << setprecision(3) << setw(9) << tmax << "\"/>"
           << endl;
     }
  }
#endif  
  
  
  //ewd save last charge density 
  if (cell_dyn == "LOCKED") {
    s_.rhog_last.resize(s_.wf.nspin());
    for ( int ispin = 0; ispin < s_.wf.nspin(); ispin++ )
    {
      if (s_.wf.spinactive(ispin))
      {
        s_.rhog_last[ispin].resize(rhog_in[ispin].size());
        for ( int i=0; i < rhog_in[ispin].size(); i++ )
          s_.rhog_last[ispin][i] = rhog_in[ispin][i];
      }
    }
  }
  
  if ( atoms_move )
  {

    // need eigenvalues to compute forces w. ultrasoft
    if (ultrasoft) { 
      ef_.energy(true,dwf,false,fion,false,sigma_eks);
      tmap["diag"].start();
      //s_.wf.diag(dwf,compute_eigvec);
      s_.wf.diag(dwf,true);
      tmap["diag"].stop();
      s_.wf.printeig();
    }

    // compute ionic forces at last position to update velocities
    // consistently with last position
    tmap["charge"].start();
    cd_.update_density();
    tmap["charge"].stop();

    ef_.update_vhxc();
    const bool compute_forces = true;
    double energy =
      ef_.energy(false,dwf,compute_forces,fion,compute_stress,sigma_eks);

    // average forces over symmetric atoms
    if ( compute_forces && s_.symmetries.nsym() > 0) {
      int nsym_ = s_.symmetries.nsym();
      for ( int is = 0; is < atoms.atom_list.size(); is++ ) {
        for ( int ia = 0; ia < atoms.atom_list[is].size(); ia++ ) {
          // start with identity symmetry operation
          D3vector fsym(fion[is][3*ia],fion[is][3*ia+1],fion[is][3*ia+2]);
          for ( int isym = 0; isym < nsym_; isym++) {
            int ja = s_.atoms.symatomid(is,ia,isym);
            D3vector ftmp(fion[is][3*ja],fion[is][3*ja+1],fion[is][3*ja+2]);
            fsym = fsym + s_.symmetries.symlist[isym]->applyToVector(ftmp,false);
          }
          fion[is][3*ia] = fsym.x/(double)(nsym_+1);
          fion[is][3*ia+1] = fsym.y/(double)(nsym_+1);
          fion[is][3*ia+2] = fsym.z/(double)(nsym_+1);
        }
      }       
    }

    if (compute_forces && s_.atoms.add_fion_ext()) {
      for ( int is = 0; is < atoms.atom_list.size(); is++ ) {
        for ( int ia = 0; ia < atoms.atom_list[is].size(); ia++ ) {
          D3vector ftmp = s_.atoms.get_fion_ext(is,ia);
          fion[is][3*ia] += ftmp.x;
          fion[is][3*ia+1] += ftmp.y;
          fion[is][3*ia+2] += ftmp.z;
        }
      }
    }
      
    ionic_stepper->compute_v(energy,fion);
    // positions r0 and velocities v0 are consistent
  }

  if ( atoms_move && extrapolate_wf )
  {
    // compute wavefunction velocity after last iteration
    // s_.wfv contains the previous wavefunction

    s_.wfv->align(s_.wf);

    if (ultrasoft) {
      for ( int ispin = 0; ispin < s_.wf.nspin(); ispin++ ) {
        if (s_.wf.spinactive(ispin)) {
          for ( int ikp = 0; ikp < s_.wf.nkp(); ikp++ ) {
            if (s_.wf.kptactive(ikp)) {
              assert(s_.wf.sd(ispin,ikp) != 0);
              tmap["usfns"].start();
              s_.wf.sd(ispin,ikp)->update_usfns();
              tmap["usfns"].stop();
            }
          }
        }
      }
    }

    for ( int ispin = 0; ispin < s_.wf.nspin(); ispin++ )
    {
      if (s_.wf.spinactive(ispin))
      {
        for ( int ikp = 0; ikp < s_.wf.nkp(); ikp++ )
        {
          if (s_.wf.kptactive(ikp))
          {
            assert(s_.wf.sd(ispin,ikp) != 0);
            double* c = (double*) s_.wf.sd(ispin,ikp)->c().cvalptr();
            double* cm = (double*) s_.wfv->sd(ispin,ikp)->c().cvalptr();
            const int mloc = s_.wf.sd(ispin,ikp)->c().mloc();
            const int nloc = s_.wf.sd(ispin,ikp)->c().nloc();
            const int len = 2*mloc*nloc;
            double dt_inv;
            if (dt == 0.0) dt_inv = 0.0;
            else dt_inv = 1.0 / dt;
            if ( ntc_extrapolation )
            {
              double* cmm = (double*) wfmm->sd(ispin,ikp)->c().cvalptr();
              for ( int i = 0; i < len; i++ )
              {
                const double x = c[i];
                const double xmm = cmm[i];
                cm[i] = dt_inv * ( x - xmm );
              }
              tmap["gram"].start();
              s_.wf.sd(ispin,ikp)->gram();
              tmap["gram"].stop();
            }
            else // normal extrapolation or asp_extrapolation
            {
              for ( int i = 0; i < len; i++ )
              {
                const double x = c[i];
                const double xm = cm[i];
                cm[i] = dt_inv * ( x - xm );
              }
              tmap["gram"].start();
              s_.wf.sd(ispin,ikp)->gram();
              tmap["gram"].stop();
            }
          }
        }
      }
    }

    // compute ionic forces at last position to update velocities
    // consistently with last position
    tmap["charge"].start();
    cd_.update_density();
    tmap["charge"].stop();

    ef_.update_vhxc();
    const bool compute_forces = true;
    double energy =
      ef_.energy(false,dwf,compute_forces,fion,compute_stress,sigma_eks);

    // average forces over symmetric atoms
    if ( compute_forces && s_.symmetries.nsym() > 0) {
      int nsym_ = s_.symmetries.nsym();
      for ( int is = 0; is < atoms.atom_list.size(); is++ ) {
        for ( int ia = 0; ia < atoms.atom_list[is].size(); ia++ ) {
          // start with identity symmetry operation
          D3vector fsym(fion[is][3*ia],fion[is][3*ia+1],fion[is][3*ia+2]);
          for ( int isym = 0; isym < nsym_; isym++) {
            int ja = s_.atoms.symatomid(is,ia,isym);
            D3vector ftmp(fion[is][3*ja],fion[is][3*ja+1],fion[is][3*ja+2]);
            fsym = fsym + s_.symmetries.symlist[isym]->applyToVector(ftmp,false);
          }
          fion[is][3*ia] = fsym.x/(double)(nsym_+1);
          fion[is][3*ia+1] = fsym.y/(double)(nsym_+1);
          fion[is][3*ia+2] = fsym.z/(double)(nsym_+1);
        }
      }       
    }

    if (compute_forces && s_.atoms.add_fion_ext()) {
      for ( int is = 0; is < atoms.atom_list.size(); is++ ) {
        for ( int ia = 0; ia < atoms.atom_list[is].size(); ia++ ) {
          D3vector ftmp = s_.atoms.get_fion_ext(is,ia);
          fion[is][3*ia] += ftmp.x;
          fion[is][3*ia+1] += ftmp.y;
          fion[is][3*ia+2] += ftmp.z;
        }
      }
    }
      
    ionic_stepper->compute_v(energy,fion);
    // positions r0 and velocities v0 are consistent
  }
  else if (wf_dyn != "SOTD")
  {
    // delete wavefunction velocities
     // AS: only if we don't need them anymore; for SOTD we do need them later on
     if ( s_.wfv != 0 )
      delete s_.wfv;
    s_.wfv = 0;
  }

  delete mlwft;

  // delete steppers
  delete wf_stepper;
  delete ionic_stepper;
  delete cell_stepper;

  // delete preconditioner
  if ( use_preconditioner ) delete preconditioner;
  if ( ntc_extrapolation || asp_extrapolation ) delete wfmm;

  // delete Hugoniostat
  if (hugstat != 0) delete hugstat;

  for (int ispin = 0; ispin < nspin; ispin++) 
    if (s_.wf.spinactive(ispin))
      delete mixer[ispin];

  initial_atomic_density = false;
}

////////////////////////////////////////////////////////////////////////////////
void BOSampleStepper::get_forces(vector<vector<double> >& f) const {
  AtomSet& atoms = s_.atoms;
  if (f.size() < atoms.atom_list.size()) 
    f.resize(atoms.atom_list.size());

  for ( int is = 0; is < atoms.atom_list.size(); is++ ) {
    if (f[is].size() != 3*atoms.atom_list[is].size())
      f[is].resize(3*atoms.atom_list[is].size());
    int i = 0;
    for ( int ia = 0; ia < atoms.atom_list[is].size(); ia++ ) {
      f[is][i] = fion[is][i];
      f[is][i+1] = fion[is][i+1];
      f[is][i+2] = fion[is][i+2];
      i+=3;
    }
  }

  // get forces on MM atoms, if any
  if (f.size() < atoms.atom_list.size()+atoms.mmatom_list.size())
    f.resize(atoms.atom_list.size() + atoms.mmatom_list.size());
  const int offset = atoms.atom_list.size();
  for ( int is = 0; is < atoms.mmatom_list.size(); is++ ) {
    if (f[is+offset].size() != 3*atoms.mmatom_list[is].size())
      f[is+offset].resize(3*atoms.atom_list[is].size());
    int i = 0;
    for ( int ia = 0; ia < atoms.mmatom_list[is].size(); ia++ ) {
      D3vector t = atoms.mmatom_list[is][ia]->position();
      f[is+offset][i] = fion[is+offset][i];
      f[is+offset][i+1] = fion[is+offset][i+1];
      f[is+offset][i+2] = fion[is+offset][i+2];
      i+=3;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
double BOSampleStepper::get_energy(string ename) {
  if (ename == "ekin") { return ef_.ekin(); }
  else if (ename == "econf") { return ef_.econf(); }
  else if (ename == "eps") { return ef_.eps(); }
  else if (ename == "enl") { return ef_.enl(); }
  else if (ename == "ehart") { return ef_.ehart(); }
  else if (ename == "ecoul") { return ef_.ecoul(); }
  else if (ename == "exc") { return ef_.exc(); }
  else if (ename == "esr") { return ef_.esr(); }
  else if (ename == "eself") { return ef_.eself(); }
  else if (ename == "ets") { return ef_.ets(); }
  else if (ename == "etotal") { return ef_.etotal(); }
  else { return ef_.etotal(); }
}

////////////////////////////////////////////////////////////////////////////////
valarray<double> BOSampleStepper::get_stress(string sname) {
  if (sname == "total") { return sigma; }
  else if (sname == "ext") { return sigma_ext; }
  else if (sname == "eks") { return sigma_eks; }
  else if (sname == "kin") { return sigma_kin; }
  else { 
    valarray<double> nullv(1);
    nullv[0] = 0.0;
    return nullv; 
  }
}
