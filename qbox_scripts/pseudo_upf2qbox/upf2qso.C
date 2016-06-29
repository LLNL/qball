//
// upf2qso.C: transform a UPF pseudopotential to QSO format
//
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <cassert>
#include "spline.h"
#include "PeriodicTable.h"
using namespace std;

int seek_tag(string tag)
{
  // skip to line following tag
  bool done = false;
  string buf, s;
  istringstream is;

  while ( !done )
  {
    getline(cin,buf);
    is.clear();
    is.str(buf);
    is >> s;
    done = ( s == tag );
    if ( cin.eof() )
    {
      cerr << " found EOF before " << tag << endl;
      return 1;
    }
  }
  return 0;
}

int main()
{

  int nplin = 1401;
  const double mesh_spacing = 0.005;
  bool ignore_rinner = false;  // set true to skip rinner, qfcoef sections

  PeriodicTable pt;
  string buf,s;
  istringstream is;
  int status;
  bool ultrasoft;
  bool nlcc;
  
  string upf_pp_info;
  status = seek_tag("<PP_INFO>");
  bool done = false;
  while (!done)
  {
    getline(cin,buf);
    is.clear();
    is.str(buf);
    is >> s;
    done = ( s == "</PP_INFO>" );
    if ( !done )
    {
      upf_pp_info += buf + '\n';
    }
  }

  // remove all '<' and '>' characters from the PP_INFO field
  // for XML compatibility
  string::size_type p = upf_pp_info.find_first_of("<>");
  while ( p != string::npos )
  {
    upf_pp_info[p] = ' '; 
    p = upf_pp_info.find_first_of("<>");
  }

  status = seek_tag("<PP_HEADER>");

  // version number (ignore)
  getline(cin,buf);

  // element symbol
  string upf_symbol;
  getline(cin,buf);
  is.clear();
  is.str(buf);
  is >> upf_symbol;

  // get atomic number and mass
  const int atomic_number = pt.z(upf_symbol);
  const double mass = pt.mass(upf_symbol);

  // NC vs. US flag
  string upf_ncflag;
  getline(cin,buf);
  is.clear();
  is.str(buf);
  is >> upf_ncflag;
  if ( upf_ncflag == "US" )
    ultrasoft = true;
  else if ( upf_ncflag == "NC" )
    ultrasoft = false;
  else {
    cerr << " not a supported potential type" << endl;
    cerr << " NC/US flag: " << upf_ncflag << endl;
    return 1;
  }
  // NLCC flag
  string upf_nlcc_flag;
  getline(cin,buf);
  is.clear();
  is.str(buf);
  is >> upf_nlcc_flag;
  nlcc = false;
  if ( upf_nlcc_flag != "F" )
  {
     nlcc = true;
  }

  // XC functional (add in description)
  string upf_xcf[4];
  getline(cin,buf);
  is.clear();
  is.str(buf);
  is >> upf_xcf[0] >> upf_xcf[1] >> upf_xcf[2] >> upf_xcf[3];

  // add XC functional information to description
  upf_pp_info += upf_xcf[0] + ' ' + upf_xcf[1] + ' ' + 
                 upf_xcf[2] + ' ' + upf_xcf[3] + '\n';

  // Z valence
  double upf_zval;
  getline(cin,buf);
  is.clear();
  is.str(buf);
  is >> upf_zval;
  
  // Total energy (ignore)
  getline(cin,buf);
  
  // suggested cutoff (ignore)
  getline(cin,buf);

  // max angular momentum
  int upf_lmax;
  getline(cin,buf);
  is.clear();
  is.str(buf);
  is >> upf_lmax;

  // number of points in mesh
  int upf_mesh_size;
  getline(cin,buf);
  is.clear();
  is.str(buf);
  is >> upf_mesh_size;

  // number of wavefunctions, number of projectors
  int upf_nwf, upf_nproj;
  getline(cin,buf);
  is.clear();
  is.str(buf);
  is >> upf_nwf >> upf_nproj;

  // Wavefunctions
  vector<string> upf_shell(upf_nwf);
  vector<int> upf_l(upf_nwf);
  vector<double> upf_occ(upf_nwf);
  // skip header
  getline(cin,buf);
  for ( int ip = 0; ip < upf_nwf; ip++ )
  {
    getline(cin,buf);
    is.clear();
    is.str(buf);
    is >> upf_shell[ip] >> upf_l[ip] >> upf_occ[ip];
  }
  status = seek_tag("</PP_HEADER>");
  if ( status != 0 ) return status;

  // read mesh
  status = seek_tag("<PP_MESH>");
  if ( status != 0 ) return status;
  status = seek_tag("<PP_R>");
  if ( status != 0 ) return status;
  vector<double> upf_r(upf_mesh_size);
  for ( int i = 0; i < upf_mesh_size; i++ )
   cin >> upf_r[i];
  status = seek_tag("</PP_R>");
  if ( status != 0 ) return status;
  status = seek_tag("<PP_RAB>");
  if ( status != 0 ) return status;
  vector<double> upf_rab(upf_mesh_size);
  for ( int i = 0; i < upf_mesh_size; i++ )
   cin >> upf_rab[i];
  status = seek_tag("</PP_RAB>");
  if ( status != 0 ) return status;
  status = seek_tag("</PP_MESH>");
  if ( status != 0 ) return status;

  vector<double> upf_nlcc;
  if (nlcc)
  {
     status = seek_tag("<PP_NLCC>");
     if ( status != 0 ) return status;
     upf_nlcc.resize(upf_mesh_size);
     for ( int i = 0; i < upf_mesh_size; i++ )
        cin >> upf_nlcc[i];
     status = seek_tag("</PP_NLCC>");
     if ( status != 0 ) return status;
  }
  
  status = seek_tag("<PP_LOCAL>");
  if ( status != 0 ) return status;
  vector<double> upf_vloc(upf_mesh_size);
  for ( int i = 0; i < upf_mesh_size; i++ )
    cin >> upf_vloc[i];
  status = seek_tag("</PP_LOCAL>");
  if ( status != 0 ) return status;

  status = seek_tag("<PP_NONLOCAL>");
  if ( status != 0 ) return status;
  vector<vector<double> > upf_vnl;
  upf_vnl.resize(upf_nproj);
  vector<int> upf_proj_l(upf_nproj);
  for ( int j = 0; j < upf_nproj; j++ )
  {
    status = seek_tag("<PP_BETA>");
    if ( status != 0 ) return status;
    int ip, l, np;
    cin >> ip >> l;
    while ( cin.get() != '\n' );
    assert(ip-1 < upf_nproj);
    assert(l <= upf_lmax);
    upf_proj_l[ip-1] = l;
    cin >> np;
    upf_vnl[j].resize(np);
    for ( int i = 0; i < np; i++ )
      cin >> upf_vnl[j][i];
    status = seek_tag("</PP_BETA>");
    if ( status != 0 ) return status;
  }
  status = seek_tag("<PP_DIJ>");
  if ( status != 0 ) return status;
  int upf_ndij;
  cin >> upf_ndij;
  while ( cin.get() != '\n' );
  if ( !ultrasoft && upf_ndij != upf_nproj )
  {
    cerr << " Norm-conserving potential:  number of non-zero Dij differs from number of projectors" << endl;
    return 1;
  }
  vector<vector<double> > upf_d;
  upf_d.resize(upf_nproj);
  for (int i=0; i<upf_nproj; i++)
    upf_d[i].resize(upf_nproj);
  for (int i=0; i<upf_nproj; i++)
    for (int j=0; j<upf_nproj; j++)
      upf_d[i][j] = 0.0;
  
  for ( int i = 0; i < upf_ndij; i++ )
  {
    int m,n;
    double dread;
    cin >> m >> n >> dread;
    upf_d[m-1][n-1] = dread;
    if ( !ultrasoft && m != n )
    {
      cerr << " Norm-conserving potential has off-diagonal Dij elements" << endl;
      cerr << " m=" << m << " n=" << n << endl;
      return 1;
    }
  }
  //status = seek_tag("</PP_DIJ>");
  //if ( status != 0 ) return status;

  vector<vector<double> > upf_qij;
  vector<double> upf_rinner;
  vector<vector<double> > upf_qfcoef;
  int upf_nqf;
  int upf_nqfcoef;
  vector<double> upf_qint;
  vector<int> upf_qind1;
  vector<int> upf_qind2;
  vector<int> upf_qij_l;
  int upf_nqij = 0; // total number of QIJ fns

  if (ultrasoft) {
    status = seek_tag("<PP_QIJ>");
    if ( status != 0 ) return status;
    cin >> upf_nqf;
    while ( cin.get() != '\n' );

    int rsize = 2*upf_lmax+1;

    if (!ignore_rinner) { // use to disable rinner for hydrogen

      status = seek_tag("<PP_RINNER>");
      if ( status != 0 ) return status;
      upf_rinner.resize(rsize);
      for ( int i = 0; i<rsize; i++ ) {
        int n;
        getline(cin,buf);
        is.clear();
        is.str(buf);
        is >> n >> upf_rinner[i];
      }
      status = seek_tag("</PP_RINNER>");
      if ( status != 0 ) return status;
      
    }

  
    for (int i=1; i<=upf_nproj; i++)
      upf_nqij += i;

    upf_qij.resize(upf_nqij);
    upf_qfcoef.resize(upf_nqij);
    upf_qint.resize(upf_nqij);
    upf_qind1.resize(upf_nqij);
    upf_qind2.resize(upf_nqij);
    upf_qij_l.resize(upf_nqij);
    for ( int i = 0; i < upf_nqij; i++ )
    {
      int m,n,l;
      cin >> m >> n >> l;
      while ( cin.get() != '\n' );
      upf_qind1[i] = m;
      upf_qind2[i] = n;
      upf_qij_l[i] = l;

      // read Q_int
      cin >> upf_qint[i];
      while ( cin.get() != '\n' );
      
      // Read Q_ij
      upf_qij[i].resize(upf_mesh_size);
      for (int j = 0; j < upf_mesh_size; j++)
        cin >> upf_qij[i][j];

      // Read qfcoef
      if (!ignore_rinner) { // use to disable rinner for hydrogen
        status = seek_tag("<PP_QFCOEF>");
        if ( status != 0 ) return status;
        upf_nqfcoef = rsize*upf_nqf;
        upf_qfcoef[i].resize(upf_nqfcoef);
        vector<double> qftmp(upf_nqfcoef);
        for (int j = 0; j < upf_nqfcoef; j++)
          cin >> qftmp[j];
        for (int j = 0; j < upf_nqfcoef; j++)
          upf_qfcoef[i][j] = qftmp[j];
        status = seek_tag("</PP_QFCOEF>");
        if ( status != 0 ) return status;
      }
    }
    status = seek_tag("</PP_QIJ>");
    if ( status != 0 ) return status;
  }

  
  status = seek_tag("</PP_NONLOCAL>");
  if ( status != 0 ) return status;


  vector<int> iproj;
  int upf_llocal = 0;
  if (ultrasoft) {
    // make table iproj[l] mapping l to iproj
    // vnl(l) is in vnl[iproj[l]] if iproj[l] > -1
    // vlocal if iproj[llocal] = -1
    iproj.resize(upf_lmax+1);
    for ( int l = 0; l <= upf_lmax; l++ )
      iproj[l] = -1;
    for ( int l = 0; l <= upf_lmax; l++ )
      for ( int j = 0; j < upf_nproj; j++ )
        iproj[upf_proj_l[j]] = j;

    // determine angular momentum of local potential in UPF file
    for ( int l = 0; l <= upf_lmax; l++ )
      if ( iproj[l] == -1 )
        upf_llocal = l;
  }

  status = seek_tag("<PP_PSWFC>");
  vector<vector<double> > upf_wf;
  vector<int> upf_wf_l(upf_nwf);
  vector<double> upf_wf_occ(upf_nwf);
  upf_wf.resize(upf_nwf);
  for ( int j = 0; j < upf_nwf; j++ )
  {
    upf_wf[j].resize(upf_mesh_size);
    string label;
    cin >> label >> upf_wf_l[j] >> upf_wf_occ[j];
    while ( cin.get() != '\n' );
    for ( int i = 0; i < upf_mesh_size; i++ )
      cin >> upf_wf[j][i];
  }
  status = seek_tag("</PP_PSWFC>");
  if ( status != 0 ) return status;

  // print summary
  cerr << "PP_INFO:" << endl << upf_pp_info << endl;
  cerr << "Element: " << upf_symbol << endl;
  cerr << "NC or US: " << upf_ncflag << endl;
  cerr << "NLCC: " << upf_nlcc_flag << endl;
  cerr << "XC: " << upf_xcf[0] << " " << upf_xcf[1] << " "
       << upf_xcf[2] << " " << upf_xcf[3] << endl;
  cerr << "Zv: " << upf_zval << endl;
  cerr << "lmax: " << upf_lmax << endl;
  cerr << "llocal: " << upf_llocal << endl;
  cerr << "nwf: " << upf_nwf << endl;
  cerr << "mesh_size: " << upf_mesh_size << endl;

  // interpolate functions on linear mesh

  // compute delta_vnl[l][i] on the upf log mesh

  // divide the projector function by the wavefunction, except if
  // the wavefunction amplitude is smaller than tol, outside of rcut_divide.
  const double tol = 1.e-3;
  const double rcut_divide = 1.0;
  vector<vector<double> > delta_vnl;
  if (!ultrasoft) {
    delta_vnl.resize(upf_nproj);
    for ( int j = 0; j < upf_nproj; j++ )
    {
      delta_vnl[j].resize(upf_wf[j].size());
      for ( int i = 0; i < delta_vnl[j].size(); i++ )
      {
        double den = upf_wf[j][i];
        if ( upf_r[i] < rcut_divide || fabs(den) > tol )
          delta_vnl[j][i] = upf_vnl[j][i] / upf_wf[j][i];
        else
          delta_vnl[j][i] = 0.0;
      }
    }
  }

  vector<vector<double> > vps;
  if (!ultrasoft) {
    vps.resize(upf_nproj+1);
    for ( int j = 0; j < upf_nproj; j++ )
    {
      vps[j].resize(upf_mesh_size);
      for ( int i = 0; i < delta_vnl[j].size(); i++ )
        vps[j][i] = upf_vloc[i] + delta_vnl[j][i];
    }
  }
  
  // interpolate vloc
  vector<double> f(upf_mesh_size), fspl(upf_mesh_size);

  // factor 0.5: convert from Ry in UPF to Hartree in QSO
  for ( int i = 0; i < upf_vloc.size(); i++ )
    f[i] = 0.5 * upf_vloc[i];

  int n = upf_vloc.size();
  int bcnat_left = 0;
  double yp_left = 0.0;
  int bcnat_right = 1;
  double yp_right = 0.0;
  spline(n,&upf_r[0],&f[0],yp_left,yp_right,
         bcnat_left,bcnat_right,&fspl[0]);
  
  vector<double> vloc_lin(nplin);
  for ( int i = 1; i < nplin; i++ )
  {
    double r = i * mesh_spacing;
    splint(n,&upf_r[0],&f[0],&fspl[0],r,&vloc_lin[i]);
  }
  // use value closest to origin for r=0
  vloc_lin[0] = 0.5 * upf_vloc[0];

  // interpolate vps[j], j=0, nproj-1 
  vector<vector<double> > vps_lin;
  vector<vector<double> > qij_lin;
  vector<double> rhonlcc_lin;
  if (!ultrasoft) { 
    vps_lin.resize(vps.size());
    for ( int j = 0; j < vps.size(); j++ )
    {
      vps_lin[j].resize(nplin);
    }

    for ( int j = 0; j < upf_nproj; j++ )
    {
      // factor 0.5: convert from Ry in UPF to Hartree in QSO
      for ( int i = 0; i < upf_vloc.size(); i++ )
        f[i] = 0.5 * vps[j][i];

      int n = upf_vloc.size();
      int bcnat_left = 0;
      double yp_left = 0.0;
      int bcnat_right = 1;
      double yp_right = 0.0;
      spline(n,&upf_r[0],&f[0],yp_left,yp_right,
             bcnat_left,bcnat_right,&fspl[0]);
    
      for ( int i = 1; i < nplin; i++ )
      {
        double r = i * mesh_spacing;
        splint(n,&upf_r[0],&f[0],&fspl[0],r,&vps_lin[j][i]);
      }
      vps_lin[j][0] = 0.5 * vps[j][0];
    }
  }
  else { // ultrasoft


    //ewd DEBUG
    ofstream brlog("betar.log.dat");
    for ( int j = 0; j < upf_nproj; j++ ) {
      brlog << "# betar " << j << endl << endl;
      for ( int i = 0; i < upf_vnl[j].size(); i++ )
        brlog << upf_r[i] << "  " << upf_vnl[j][i] << endl;
    }
    ofstream blog("beta.log.dat");
    for ( int j = 0; j < upf_nproj; j++ ) {
      blog << "# betar " << j << endl << endl;
      for ( int i = 1; i < upf_vnl[j].size(); i++ )
        blog << upf_r[i] << "  " << upf_vnl[j][i]/upf_r[i] << endl;
    }
    //ewd DEBUG

      
    // interpolate beta fns
    vps_lin.resize(upf_nproj+1);
    for ( int j = 0; j < upf_nproj; j++ )
      vps_lin[j].resize(nplin);
    
    for ( int j = 0; j < upf_nproj; j++ )
    {
      // factor 0.5: convert from Ry in UPF to Hartree in QSO <-- EWD, do we need this??
      for ( int i = 0; i < upf_vnl[j].size(); i++ )
        //f[i] = 0.5 * upf_vnl[j][i];
        f[i] = upf_vnl[j][i];

      int n = upf_vnl[j].size();
      int bcnat_left = 0;
      double yp_left = 0.0;
      int bcnat_right = 0;  //ewd:  use flat boundary conditions for beta fns
      double yp_right = 0.0;
      spline(n,&upf_r[0],&f[0],yp_left,yp_right,
             bcnat_left,bcnat_right,&fspl[0]);
    
      double rmax = upf_r[n-1];
      for ( int i = 1; i < nplin; i++ )
      {
        double r = i * mesh_spacing;
        if (r < rmax) {
          double vspl;
          splint(n,&upf_r[0],&f[0],&fspl[0],r,&vspl);
          vps_lin[j][i] = vspl/r;
        }
        else
          vps_lin[j][i] = 0.0;
      }
      //vps_lin[j][0] = 0.5 * upf_vnl[j][0];
      //vps_lin[j][0] = upf_vnl[j][0];
      vps_lin[j][0] = 0.0;
    }

    // interpolate upf_qij = Q_ij(r)*r^2 
    qij_lin.resize(upf_nqij);
    for ( int j = 0; j < upf_nqij; j++ )
      qij_lin[j].resize(nplin);

    for ( int j = 0; j < upf_nqij; j++ )
    {
      for ( int i = 0; i < upf_mesh_size; i++ )
        f[i] = upf_qij[j][i];
      
      int n = upf_mesh_size;
      int bcnat_left = 0;
      double yp_left = 0.0;
      int bcnat_right = 1;
      double yp_right = 0.0;
      spline(n,&upf_r[0],&f[0],yp_left,yp_right,
             bcnat_left,bcnat_right,&fspl[0]);
      
      double rmax = upf_r[n-1];
      for ( int i = 1; i < nplin; i++ )
      {
        double r = i * mesh_spacing;
        if (r < rmax) {
          double qspl;
          splint(n,&upf_r[0],&f[0],&fspl[0],r,&qspl);
          qij_lin[j][i] = qspl/(r*r);
        }
        else
          qij_lin[j][i] = 0.0;
      }
      qij_lin[j][0] = 0.0;
    }

    //interpolate non-linear core correction
    if (nlcc)
    {
       rhonlcc_lin.resize(nplin);
       for ( int i = 0; i < upf_mesh_size; i++ )
          f[i] = upf_nlcc[i];
      
       int n = upf_mesh_size;
       int bcnat_left = 0;
       double yp_left = 0.0;
       int bcnat_right = 1;
       double yp_right = 0.0;
       spline(n,&upf_r[0],&f[0],yp_left,yp_right,
              bcnat_left,bcnat_right,&fspl[0]);
      
       double rmax = upf_r[n-1];
       for ( int i = 1; i < nplin; i++ )
       {
          double r = i * mesh_spacing;
          if (r < rmax) {
             double qspl;
             splint(n,&upf_r[0],&f[0],&fspl[0],r,&qspl);
             rhonlcc_lin[i] = qspl;
          }
          else
             rhonlcc_lin[i] = 0.0;
       }
       rhonlcc_lin[0] = 0.0;
    }

  }
  
  // write potentials in gnuplot format on file vlin.dat
  ofstream vlin("vlin.dat");
  if (ultrasoft) {
    vlin << "# vloc" << endl;
    for ( int i = 0; i < nplin; i++ )
      vlin << i*mesh_spacing << " " << vloc_lin[i] << endl;
    vlin << endl << endl;

    // write out beta fns
    ofstream brlin("betar.lin.dat");
    for ( int j = 0; j < upf_nproj ; j++) {
      brlin << "# beta " << j << endl;
      for ( int i = 0; i < nplin; i++ )
        brlin << i*mesh_spacing << " " << setprecision(12) << (double)i*mesh_spacing*vps_lin[j][i] << endl;
      brlin << endl << endl;
    }
    brlin.close();
    ofstream blin("beta.lin.dat");
    for ( int j = 0; j < upf_nproj ; j++) {
      blin << "# beta " << j << endl;
      for ( int i = 0; i < nplin; i++ )
        blin << i*mesh_spacing << " " << setprecision(12) << vps_lin[j][i] << endl;
      blin << endl << endl;
    }
    blin.close();
    
    //write out qij fns
    ofstream qijrsqlin("qijrsqlin.dat");
    ofstream qijrsqlog("qijrsqlog.dat");
    for ( int j = 0; j < upf_nqij; j++ ) {
      qijrsqlin << "# qij " << j << endl << endl;
      for ( int i = 1; i < nplin; i++ )
      {
        double r = i * mesh_spacing;
        qijrsqlin << r << "   " << r*r*qij_lin[j][i] << endl;
      }
      qijrsqlog << "# qij " << j << endl << endl;
      for ( int i = 0; i < upf_mesh_size; i++ )
      {
        qijrsqlog << upf_r[i] << "   " << upf_qij[j][i] << endl;
      }
    }
    ofstream qijlin("qijlin.dat");
    ofstream qijlog("qijlog.dat");
    for ( int j = 0; j < upf_nqij; j++ ) {
      qijlin << "# qij " << j << endl << endl;
      for ( int i = 1; i < nplin; i++ )
      {
        double r = i * mesh_spacing;
        qijlin << r << "   " << qij_lin[j][i] << endl;
      }
      qijlog << "# qij " << j << endl << endl;
      for ( int i = 1; i < upf_mesh_size; i++ )
      {
        qijlog << upf_r[i] << "   " << upf_qij[j][i]/(upf_r[i]*upf_r[i]) << endl;
      }
    }
  }
  else { // norm-conserving
    for ( int l = 0; l <= upf_lmax; l++ )
    {
      vlin << "# v, l=" << l << endl;
      if ( l == upf_llocal )
      {
        // l == llocal
        for ( int i = 0; i < nplin; i++ )
          vlin << i*mesh_spacing << " " << vloc_lin[i] << endl;
        vlin << endl << endl;
      }
      else
      {
        for ( int i = 0; i < nplin; i++ )
          vlin << i*mesh_spacing << " " << vps_lin[l][i] << endl;
        vlin << endl << endl;
      }
    }
  }
  
  // interpolate wavefunctions on the linear mesh
  

  vector<vector<double> > wf_lin;
  if (!ultrasoft)
  {
    wf_lin.resize(upf_nwf);
    for ( int j = 0; j < upf_nwf; j++ )
    {
      wf_lin[j].resize(nplin);
      for ( int i = 0; i < upf_wf[j].size(); i++ )
        f[i] = upf_wf[j][i] / upf_r[i];

      int n = upf_wf[j].size();
      int bcnat_left = 1;
      double yp_left = 0.0;
      int bcnat_right = 1;
      double yp_right = 0.0;
      spline(n,&upf_r[0],&f[0],yp_left,yp_right,
             bcnat_left,bcnat_right,&fspl[0]);
    
      for ( int i = 1; i < nplin; i++ )
      {
        double r = i * mesh_spacing;
        splint(n,&upf_r[0],&f[0],&fspl[0],r,&wf_lin[j][i]);
      }
      // compute value at origin, depending on angular momentum
      if ( upf_wf_l[j] == 0 )
        // take value closest to r=0
        wf_lin[j][0] = upf_wf[j][0]/upf_r[0];
      else
        wf_lin[j][0] = 0.0;

      vlin << "# phi, l=" << upf_l[j] << endl;
      for ( int i = 0; i < nplin; i++ )
        vlin << i*mesh_spacing << " " << wf_lin[j][i] << endl;
      vlin << endl << endl;
    }
  }
  
  cerr << " interpolation done" << endl;

  if (ultrasoft)
  {
    // output potential on log mesh
    ofstream vout("vlog.dat");
    for ( int i = 0; i < upf_vloc.size(); i++ )
      vout << upf_r[i] << " " << 0.5*upf_vloc[i] << endl;
  }

#if 0
  // output potential on log mesh
  ofstream vout("v.dat");
  for ( int i = 0; i < upf_vloc.size(); i++ )
    vout << upf_r[i] << " " << upf_vloc[i] << endl;
  vout << endl << endl;
  for ( int j = 0; j < upf_nproj; j++ )
  {
    for ( int i = 0; i < upf_vnl[j].size(); i++ )
      vout << upf_r[i] << " " << upf_vnl[j][i] << endl;
    vout << endl << endl;
  }
#endif
  
  // Generate QSO file

  // output potential in QSO format
  cout << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
  cout << "<fpmd:species xmlns:fpmd=\"http://www.quantum-simulation.org/ns/fpmd/fpmd-1.0\"" << endl;
  cout << "  xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"" << endl;
  cout << "  xsi:schemaLocation=\"http://www.quantum-simulation.org/ns/fpmd/fpmd-1.0"  << endl;
  cout << "  species.xsd\">" << endl;
  cout << "<description>" << endl;
  cout << "Translated from UPF format by upf2qso" << endl;
  cout << upf_pp_info;
  cout << "</description>" << endl;
  cout << "<symbol>" << upf_symbol << "</symbol>" << endl;
  cout << "<atomic_number>" << atomic_number << "</atomic_number>" << endl;
  cout << "<mass>" << mass << "</mass>" << endl;


  if (!ultrasoft) { 
    cout << "<norm_conserving_pseudopotential>" << endl;
    cout << "<valence_charge>" << upf_zval << "</valence_charge>" << endl;
    cout << "<lmax>" << upf_lmax << "</lmax>" << endl;
    cout << "<llocal>" << upf_llocal << "</llocal>" << endl;
    cout << "<nquad>0</nquad>" << endl;
    cout << "<rquad>0.0</rquad>" << endl;
    cout << "<mesh_spacing>" << mesh_spacing << "</mesh_spacing>" << endl;

    for ( int l = 0; l <= upf_lmax; l++ )
    {
      cout << "<projector l=\"" << l << "\" size=\"" << nplin << "\">" 
           << endl;
      cout << "<radial_potential>" << endl;
      if ( l == upf_llocal ) 
      {
        // l == llocal
        cout.setf(ios::fixed,ios::floatfield);
        for ( int i = 0; i < nplin; i++ )
          cout << setprecision(12) << vloc_lin[i] << endl;
      }
      else
      {
        int lind = -1;
        for (int i=0; i<upf_nproj; i++)
          if (l == upf_proj_l[i])
            lind = i;

        //ewd DEBUG
        cerr << "l = " << l << ", lind = " << lind << endl;

        cout.setf(ios::fixed,ios::floatfield);
        for ( int i = 0; i < nplin; i++ )
          cout << setprecision(12) << vps_lin[lind][i] << endl;
      }
      cout << "</radial_potential>" << endl;
      //if ( l != upf_llocal )
      //{
        cout << "<radial_function>" << endl;
        for ( int i = 0; i < nplin; i++ )
          cout << setprecision(12) << wf_lin[l][i] << endl;
        cout << "</radial_function>" << endl;
        //}
      cout << "</projector>" << endl;
    }
    cout << "</norm_conserving_pseudopotential>" << endl;
  }
  else {
    cout << "<ultrasoft_pseudopotential>" << endl;
    cout << "<valence_charge>" << upf_zval << "</valence_charge>" << endl;
    cout << "<lmax>" << upf_lmax << "</lmax>" << endl;
    cout << "<mesh_spacing>" << mesh_spacing << "</mesh_spacing>" << endl;
    cout << "<nbeta>" << upf_nproj << "</nbeta>" << endl;
    cout << "<vlocal size=\"" << nplin << "\">" << endl;
    cout.setf(ios::fixed,ios::floatfield);
    for ( int i = 0; i < nplin; i++ )
      cout << setprecision(12) << vloc_lin[i] << endl;
    cout << "</vlocal>" << endl;
    for ( int j = 0; j < upf_nproj ; j++) {
      cout << "<beta i=\"" << j << "\" l=\"" << upf_proj_l[j] << "\" size=\"" << nplin << "\">" << endl;
      for ( int i = 0; i < nplin; i++ )
        cout << setprecision(12) << vps_lin[j][i] << endl;
      cout << "</beta>" << endl;
    }    
    // print out D_ij^0 terms
    cout.setf(ios::fixed,ios::floatfield);
    cout << "<dnm_zero>" << endl;
    for (int i=0; i<upf_nproj; i++) {
      for (int j=0; j<upf_nproj; j++)
        cout << setprecision(10) << 0.5*upf_d[i][j] << "  ";
      cout << endl;
    }
    cout << "</dnm_zero>" << endl;
    // print out rinner values
    if (!ignore_rinner) {
      int rsize = 2*upf_lmax+1;
      cout << "<rinner size=\"" << rsize << "\">" << endl;
      for ( int i = 0; i<rsize; i++ ) 
        cout << upf_rinner[i] << endl;
      cout << "</rinner>" << endl;
    }
      
    for (int j=0; j<upf_nqij; j++) {
      cout << "<qnm i=\"" << j << "\" n=\"" << upf_qind1[j]-1 << "\" m=\"" << upf_qind2[j]-1 << "\" l1=\"" << upf_proj_l[upf_qind1[j]-1] << "\" l2=\"" << upf_proj_l[upf_qind2[j]-1] << "\" size=\"" << nplin << "\">" << endl;
      for ( int i = 0; i < nplin; i++ )
        cout << setprecision(12) << qij_lin[j][i] << endl;
      cout << "</qnm>" << endl;
      if (!ignore_rinner) {        
        int cnt = 0;
        for (int ltot=0; ltot<2*upf_lmax+1; ltot++) {
          cout << "<qfcoeff i=\"" << j << "\" n=\"" << upf_qind1[j]-1 << "\" m=\"" << upf_qind2[j]-1 << "\" ltot=\"" << ltot << "\" size=\"" << upf_nqf << "\">" << endl;
          for (int iq = 0; iq < upf_nqf; iq++) 
            cout << setprecision(12) << upf_qfcoef[j][iq+cnt] << endl;
          cnt+=upf_nqf;
          cout << "</qfcoeff>" << endl;
        }
      }
    }

    if (nlcc)
    {
       cout << "<rho_nlcc size=\"" << nplin << "\">" << endl;
       cout.setf(ios::fixed,ios::floatfield);
       for ( int i = 0; i < nplin; i++ )
          cout << setprecision(12) << rhonlcc_lin[i] << endl;
       cout << "</rho_nlcc>" << endl;
    }

    cout << "</ultrasoft_pseudopotential>" << endl;
  }
    
  cout << "</fpmd:species>" << endl;
  return 0;
}
