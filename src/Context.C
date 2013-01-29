////////////////////////////////////////////////////////////////////////////////
//
// Context.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Context.C,v 1.6 2008/07/02 17:06:51 draeger1 Exp $

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <string>
using namespace std;

#ifdef SCALAPACK
#include "blacs.h"
#endif

#include "Context.h"

#ifndef SCALAPACK
void Cblacs_pinfo(int *mypnum, int *nprocs)
{
    *mypnum = 0;
    *nprocs = 1;
    return;
}

void Cblacs_get(int icontxt, int what, int *val)
{
    *val=0;
    return;
}

void Cblacs_gridinit(int *icontxt, char order[], int nprow, int npcol)
{
    *icontxt=0;
    return;
}

void Cblacs_gridmap(int *icontxt, int *pmap, int ldpmap, int nprow, int npcol)
{
    pmap[0]=0;
    return;
}

void Cblacs_barrier(int icontxt, char* scope)
{
    return;
}

void Cblacs_abort(int icontxt, int errornum)
{
    return;
}

void Cblacs_gridexit(int icontxt)
{
    return;
}

void Cblacs_gridinfo(int icontxt, int *nprow, int *npcol, 
                                  int *myprow, int *mypcol)
{
    *nprow=1; *npcol=1; *myprow=0; *mypcol=0; return;
}

int Cblacs_pnum(int icontxt, int prow, int pcol)
{ return 0;}

int Csys2blacs_handle(int comm) {
  return 0;
}
void Cdgesd2d(int icontxt,int m,int n,double *A,int lda,int rdest,int cdest)
{ return; }

void Cdgerv2d(int icontxt,int m,int n,double *A,int lda,int rdest,int cdest)
{ return; }

void Cdgsum2d(int icontxt, char* scope, char* top, 
              int m, int n, double *a, int lda, int rdest, int cdest)
{ return; }

void Cdgamx2d(int icontxt, char scope[], char top[],int m,int n,double *A,
              int lda, int *ra, int *ca, int rcflag, int rdest, int cdest)
{ return; }

void Cdgamn2d(int icontxt, char scope[], char top[],int m,int n,double *A,
              int lda, int *ra, int *ca, int rcflag, int rdest, int cdest)
{ return; }

void Cdgebs2d(int icontxt, char scope[], char top[],int m,int n,double *A,
              int lda)
{ return; }

void Cdgebr2d(int icontxt, char scope[], char top[],int m,int n,double *A,
              int lda, int rsrc, int csrc)
{ return; }

void Cigesd2d(int icontxt,int m,int n,int *A,int lda,int rdest,int cdest)
{ return; }

void Cigerv2d(int icontxt,int m,int n,int *A,int lda,int rdest,int cdest)
{ return; }

void Cigsum2d(int icontxt, char* scope, char* top, 
              int m, int n, int *a, int lda, int rdest, int cdest)
{ return; }

void Cigamx2d(int icontxt, char scope[], char top[],int m,int n,int *A,int lda, 
              int *ra, int *ca, int rcflag, int rdest, int cdest)
{ return; }

void Cigamn2d(int icontxt, char scope[], char top[],int m,int n,int *A,int lda, 
              int *ra, int *ca, int rcflag, int rdest, int cdest)
{ return; }

void Cigebs2d(int icontxt, char scope[], char top[],int m,int n,int *A,
              int lda)
{ return; }

void Cigebr2d(int icontxt, char scope[], char top[],int m,int n,int *A,
              int lda, int rsrc, int csrc)
{ return; }

#endif

#if USE_MPI
#include <mpi.h>
#else
typedef int MPI_Comm;
const int MPI_COMM_NULL = 0;
#endif

struct ContextRep
{ 
  private:
  
  int ictxt_;
  int myrow_;
  int mycol_;
  int nprow_;
  int npcol_;
  int size_;
  int myproc_;
  int mype_;
  bool onpe0_;
  int coutpe_;
  bool oncoutpe_;
  bool active_;
 
  vector<int> pmap_;
  MPI_Comm comm_;
 
  // keep assignment and copy constructors private
  ContextRep& operator=(const Context& c);
  ContextRep(const Context& c);

  public:
  
  int ictxt() const { return ictxt_; }
  int myrow() const { return myrow_; }
  int mycol() const { return mycol_; }
  int nprow() const { return nprow_; }
  int npcol() const { return npcol_; }
 
  // number of processes in the context
  // returns -1 if current process is not part of this context
  int size() const { return size_; }
  // position of current process in row-major order
  // returns -1 if current process is not part of this context
  int myproc() const { return myproc_; }
  int mype() const { return mype_; }
  int pmap(int irow, int icol) const { return pmap_[irow+nprow_*icol]; }
 
  bool onpe0(void) const { return onpe0_; }
  bool oncoutpe(void) const { return oncoutpe_; }
  int coutpe(void) const { return coutpe_; }
  void set_coutpe(int num) { 
    coutpe_ = num;  
    oncoutpe_ = ( mype_ == coutpe_); 
  }
  bool active(void) const { return active_; }
  void abort(int ierr) const { Cblacs_abort(ictxt_,ierr); }
  void barrier(void) const { Cblacs_barrier(ictxt_,"A"); }
  void barrier(char scope) const { Cblacs_barrier(ictxt_,&scope); }
  
  void dsend(int m, int n, double* a, int lda, int rdest, int cdest) const
  {
     Cdgesd2d(ictxt_,m,n,a,lda,rdest,cdest);
  }
 
  void drecv(int m, int n, double* a, int lda, int rsrc, int csrc) const
  {
     Cdgerv2d(ictxt_,m,n,a,lda,rsrc,csrc);
  }
 
  void dsum(char scope, char topology, int m, int n, double* a, int lda,
            int rdest, int cdest) const
  {
    Cdgsum2d(ictxt_,&scope,&topology,m,n,a,lda,rdest,cdest);
  }
 
  void dmax(char scope, char topology, int m, int n, double* a, int lda,
            int rdest, int cdest) const
  {
    Cdgamx2d(ictxt_,&scope,&topology,m,n,a,lda,(int*)0,(int*)0,-1,rdest,cdest);
  }
 
  void dmax(char scope, char topology, int m, int n, double* a, int lda,
            int* ra, int* ca, int rcflag, int rdest, int cdest) const
  {
    Cdgamx2d(ictxt_,&scope,&topology,m,n,a,lda,ra,ca,rcflag,rdest,cdest);
  }
 
  void dmin(char scope, char topology, int m, int n, double* a, int lda,
            int* ra, int* ca, int rcflag, int rdest, int cdest) const
  {
    Cdgamn2d(ictxt_,&scope,&topology,m,n,a,lda,ra,ca,rcflag,rdest,cdest);
  }
 
  void dmin(char scope, char topology, int m, int n, double* a, int lda,
            int rdest, int cdest) const
  {
    Cdgamn2d(ictxt_,&scope,&topology,m,n,a,lda,(int*)0,(int*)0,-1,rdest,cdest);
  }
 
  void dbcast_send(char scope, char topology, 
                   int m, int n, double* a,int lda) const
  {
    Cdgebs2d(ictxt_,&scope,&topology,m,n,a,lda);
  }
 
  void dbcast_recv(char scope, char topology, int m, int n, double* a, int lda,
               int rsrc, int csrc) const
  {
    Cdgebr2d(ictxt_,&scope,&topology,m,n,a,lda,rsrc,csrc);
  }
 
  void isend(int m, int n, int* a, int lda, int rdest, int cdest) const
  {
    Cigesd2d(ictxt_,m,n,a,lda,rdest,cdest);
  }
 
  void irecv(int m, int n, int* a, int lda, int rsrc, int csrc) const
  {
    Cigerv2d(ictxt_,m,n,a,lda,rsrc,csrc);
  }
 
  void isum(char scope, char topology, int m, int n, int* a, int lda,
            int rdest, int cdest) const
  {
    Cigsum2d(ictxt_,&scope,&topology,m,n,a,lda,rdest,cdest);
  }
 
  void imax(char scope, char topology, int m, int n, int* a, int lda,
            int* ra, int* ca, int rcflag, int rdest, int cdest) const
  {
    Cigamx2d(ictxt_,&scope,&topology,m,n,a,lda,ra,ca,rcflag,rdest,cdest);
  }
 
  void imax(char scope, char topology, int m, int n, int* a, int lda,
            int rdest, int cdest) const
  {
    Cigamx2d(ictxt_,&scope,&topology,m,n,a,lda,(int*)0,(int*)0,-1,rdest,cdest);
  }
 
  void imin(char scope, char topology, int m, int n, int* a, int lda,
            int* ra, int* ca, int rcflag, int rdest, int cdest) const
  {
    Cigamn2d(ictxt_,&scope,&topology,m,n,a,lda,ra,ca,rcflag,rdest,cdest);
  }
 
  void imin(char scope, char topology, int m, int n, int* a, int lda,
            int rdest, int cdest) const
  {
    Cigamn2d(ictxt_,&scope,&topology,m,n,a,lda,(int*)0,(int*)0,-1,rdest,cdest);
  }
 
  void ibcast_send(char scope, char topology, 
                   int m, int n, int* a,int lda) const
  {
    Cigebs2d(ictxt_,&scope,&topology,m,n,a,lda);
  }
 
  void ibcast_recv(char scope, char topology, int m, int n, int* a, int lda,
               int rsrc, int csrc) const
  {
    Cigebr2d(ictxt_,&scope,&topology,m,n,a,lda,rsrc,csrc);
  }
  
  void string_send(string& s, int rdest, int cdest) const
  {
#if USE_MPI
    int len = s.size();
    isend(1,1,&len,1,rdest,cdest);
    int ilen = len/(sizeof(int)/sizeof(char));
    if ( len%(sizeof(int)/sizeof(char)) != 0 ) ilen++;
    int* ibuf = new int[ilen];
    s.copy((char*)ibuf,string::npos);
    //isend(ilen,1,ibuf,1,rdest,cdest);
    isend(ilen,1,ibuf,ilen,rdest,cdest);
    delete [] ibuf;
#endif
  }
  
  void string_recv(string& s, int rsrc, int csrc) const
  {
#if USE_MPI
    int len = -1;
    irecv(1,1,&len,1,rsrc,csrc);
    int ilen = len/(sizeof(int)/sizeof(char));
    if ( len%(sizeof(int)/sizeof(char)) != 0 ) ilen++;
    int* ibuf = new int[ilen];
    //irecv(ilen,1,ibuf,1,rsrc,csrc);
    irecv(ilen,1,ibuf,ilen,rsrc,csrc);
    s.resize(len);
    s.assign((char*)ibuf,len);
    delete [] ibuf;
#endif
  }
  
  void string_bcast(string& s, int isrc) const
  {
#if USE_MPI
    int len;
    if ( mype() == isrc )
    {
      len = s.length();
    }
    MPI_Bcast(&len,1,MPI_INT,isrc,comm());
    char* buf = new char[len+1];
    // s is initialized only on task isrc
    if ( mype() == isrc )
    {
      s.copy(buf,string::npos);
      buf[len]=0;
      assert(buf[len]=='\0');
    }
    MPI_Bcast(buf,len+1,MPI_CHAR,isrc,comm());
    s = buf;
    delete [] buf;
#endif
  }
 
  bool operator==(const ContextRep& c) const
  { return (ictxt_==c.ictxt());}
 
  // MPI communicator for this context. Returns MPI_COMM_NULL if
  // this process is not part of the context
  MPI_Comm comm(void) const { return comm_; }
 
  // Constructors

  // default global context: construct a single-row global ContextRep
  explicit ContextRep();
 
  // global ContextRep of size nprow * npcol with column major order
  explicit ContextRep(int nprow, int npcol);
 
  // construct a ContextRep of size nprow*npcol using the processes
  // in context c starting at process (irstart,icstart) 
  explicit ContextRep(const ContextRep &c, int nprow, int npcol, 
    int irstart, int icstart);
  // construct a ContextRep of same size as c, different shape
  explicit ContextRep(const ContextRep &c, int nprow, int npcol);
  // construct a Context of shape nprow x npcol from subcommunicator comm
  //  --> this avoids all calls to MPI_COMM_WORLD, allowing one to run on subset of pes
#if USE_MPI
  explicit ContextRep(const MPI_Comm &comm, int nprow, int npcol);
#endif
 
  ~ContextRep();
 
  void print(ostream& os) const;
};

////////////////////////////////////////////////////////////////////////////////
ContextRep::ContextRep() : ictxt_(-1), myrow_(-1), mycol_(-1)
{
  // construct a single row global context
  int nprocs;
  char order='R';
  Cblacs_pinfo( &mype_, &nprocs );
  nprow_ = 1;
  npcol_ = nprocs;

  Cblacs_get(0, 0, &ictxt_ );
  Cblacs_gridinit( &ictxt_, &order, nprow_, npcol_ );

  // get values of nprow_, npcol_, myrow_ and mycol_ in the new context
  if ( ictxt_ >= 0 )
    Cblacs_gridinfo(ictxt_, &nprow_, &npcol_, &myrow_, &mycol_);
  
  size_ = nprow_ * npcol_;
  myproc_ = myrow_ < 0 ? -1 : mycol_ + npcol_ * myrow_;
  onpe0_ = ( mype_ == 0 );
  coutpe_ = 0;
  oncoutpe_ = onpe0_;
  active_ = ( ictxt_ >= 0 );

  pmap_.resize(size_);
  for ( int i = 0; i < size_; i++ )
    pmap_[i] = i;

#if USE_MPI
  MPI_Group group_world, subgroup;
  MPI_Comm_group(MPI_COMM_WORLD,&group_world);
  MPI_Group_incl(group_world,size_,&pmap_[0],&subgroup);
  MPI_Comm_create(MPI_COMM_WORLD,subgroup,&comm_);
  MPI_Group_free(&group_world);
  MPI_Group_free(&subgroup);
#else
  comm_ = 0;
#endif
}

////////////////////////////////////////////////////////////////////////////////
ContextRep::ContextRep(int nprow, int npcol) :
  ictxt_(-1), myrow_(-1), mycol_(-1), nprow_(nprow), npcol_(npcol)
{
  int nprocs;
  char order = 'C';
  Cblacs_pinfo( &mype_, &nprocs );
  if ( nprocs < nprow * npcol )
  {
    cout << " nprocs=" << nprocs << endl;
    cout << " Context nprow*npcol > nprocs" << endl;
    Cblacs_abort(ictxt_, 1);
  }
  Cblacs_get(0, 0, &ictxt_ );
  Cblacs_gridinit( &ictxt_, &order, nprow, npcol );

  // get values of nprow_, npcol_, myrow_ and mycol_ in the new context
  if ( ictxt_ >= 0 )
    Cblacs_gridinfo(ictxt_, &nprow_, &npcol_, &myrow_, &mycol_);
  
  size_ = nprow_ * npcol_;
  myproc_ = Cblacs_pnum(ictxt_,myrow_,mycol_);
  onpe0_ = ( mype_ == 0 );
  coutpe_ = 0;
  oncoutpe_ = onpe0_;
  active_ = ( ictxt_ >= 0 );
  
  pmap_.resize(size_);
  // column-major order
  int i = 0;
  for ( int ic = 0; ic < npcol; ic++ )
    for ( int ir = 0; ir < nprow; ir++ )
    {
      pmap_[ir+nprow*ic] = i;
      i++;
    }

#if USE_MPI
  MPI_Group group_world, subgroup;
  MPI_Comm_group(MPI_COMM_WORLD,&group_world);
  MPI_Group_incl(group_world,size_,&pmap_[0],&subgroup);
  MPI_Comm_create(MPI_COMM_WORLD,subgroup,&comm_);
  MPI_Group_free(&group_world);
  MPI_Group_free(&subgroup);
#else
  comm_ = 0;
#endif
}
    
////////////////////////////////////////////////////////////////////////////////
ContextRep::ContextRep(const ContextRep& c, int nprow, int npcol, 
  int irstart, int icstart) :
ictxt_(-1), myrow_(-1), mycol_(-1), nprow_(nprow), npcol_(npcol)
{
  assert(c.active());
  vector<int> gmap_;
  int nprocs;

  //Cblacs_pinfo( &mype_, &nprocs );
  mype_ = c.mype();
  nprocs = c.size();

  // construct a (nprow*npcol) context using the processes in c
  // located at (irstart,icstart)

  if ( irstart < 0 || icstart < 0 || 
       irstart+nprow > c.nprow() || icstart+npcol > c.npcol() ) {
    cout << " Context::Context: cut rectangle: invalid parameters" << endl;
    Cblacs_abort(ictxt_, 1);
  }
  pmap_.resize(nprow*npcol);
  gmap_.resize(nprow*npcol);
  // build pmap
  int i = 0;
  for ( int ic = icstart; ic < icstart+npcol; ic++ )
    for ( int ir = irstart; ir < irstart+nprow; ir++ )
    {
      pmap_[i] = c.pmap(ir,ic)-c.pmap(0,0);
      gmap_[i] = ic + c.npcol()*ir;
      i++;
    }

  Cblacs_get(c.ictxt(), 10, &ictxt_ );
  Cblacs_gridmap(&ictxt_,&gmap_[0],nprow,nprow,npcol);
  // get values of nprow_, npcol_, myrow_ and mycol_ in the new context
  if ( ictxt_ >= 0 )
    Cblacs_gridinfo(ictxt_, &nprow_, &npcol_, &myrow_, &mycol_);
  size_ = nprow_ * npcol_;
  myproc_ = myrow_ < 0 ? -1 : mycol_ + npcol_ * myrow_;
  onpe0_ = ( mype_ == 0 );
  coutpe_ = c.coutpe();
  oncoutpe_ = ( mype_ == coutpe_);

  active_ = ( ictxt_ >= 0 );
#if USE_MPI
  MPI_Group c_group, subgroup;
  MPI_Comm_group(c.comm(),&c_group);
  MPI_Group_incl(c_group,size_,&pmap_[0],&subgroup);
  MPI_Comm_create(c.comm(),subgroup,&comm_);
  MPI_Group_free(&c_group);
  MPI_Group_free(&subgroup);
#else
  comm_ = 0;
#endif
}

////////////////////////////////////////////////////////////////////////////////
ContextRep::ContextRep(const ContextRep& c, int nprow, int npcol) :
ictxt_(-1), myrow_(-1), mycol_(-1), nprow_(nprow), npcol_(npcol) {
  // create context with same number of processors, different shape as c
  assert(c.active());
  vector<int> gmap_;
  int nprocs;
  Cblacs_pinfo( &mype_, &nprocs );
  // construct a (nprow*npcol) context using the processes in c
  // located at (irstart,icstart)

  assert(nprow*npcol == c.size());

  pmap_.resize(nprow*npcol);
  gmap_.resize(nprow*npcol);
  // build pmap
  int i = 0;
  for ( int ic = 0; ic < npcol; ic++ )
    for ( int ir = 0; ir < nprow; ir++ ) {
      pmap_[ir+nprow*ic] = i;
      gmap_[i] = ic + c.npcol()*ir;
      i++;
    }
  // build gmap for blacs_gridmap
  i=0;
  for ( int ic = 0; ic < c.npcol(); ic++ )
    for ( int ir = 0; ir < c.nprow(); ir++ ) {
      gmap_[i] = ic + c.npcol()*ir;
      i++;
    }

  Cblacs_get(c.ictxt(), 10, &ictxt_ );
  Cblacs_gridmap(&ictxt_,&gmap_[0],nprow,nprow,npcol);
  
  // get values of nprow_, npcol_, myrow_ and mycol_ in the new context
  if ( ictxt_ >= 0 )
    Cblacs_gridinfo(ictxt_, &nprow_, &npcol_, &myrow_, &mycol_);
  
  size_ = nprow_ * npcol_;
  myproc_ = myrow_ < 0 ? -1 : mycol_ + npcol_ * myrow_;
  onpe0_ = ( mype_ == 0 );
  coutpe_ = c.coutpe();
  oncoutpe_ = ( mype_ == coutpe_);
  active_ = ( ictxt_ >= 0 );
  
#if USE_MPI
  MPI_Group c_group, subgroup;
  MPI_Comm_group(c.comm(),&c_group);
  MPI_Group_incl(c_group,size_,&pmap_[0],&subgroup);
  MPI_Comm_create(c.comm(),subgroup,&comm_);
  MPI_Group_free(&c_group);
  MPI_Group_free(&subgroup);
#else
  comm_ = 0;
#endif
}

////////////////////////////////////////////////////////////////////////////////
#if USE_MPI
ContextRep::ContextRep(const MPI_Comm& tcomm, int nprow, int npcol) :
ictxt_(-1), myrow_(-1), mycol_(-1), nprow_(nprow), npcol_(npcol) {
  // create context of shape nprow x npcol from subcommunicator comm, avoiding
  // all Blacs calls which use MPI_COMM_WORLD (e.g. CBlacs_get)

  int ierr = MPI_Comm_dup(tcomm,&comm_);
  assert(ierr == 0);

  int nprocs;
  char order = 'C';
  MPI_Comm_size(comm_,&nprocs);
  MPI_Comm_rank(comm_,&mype_);
  int bhandle = Csys2blacs_handle(comm_);
  ictxt_ = bhandle;
  Cblacs_gridinit( &ictxt_, &order, nprow, npcol );

  // get values of nprow_, npcol_, myrow_ and mycol_ in the new context
  if ( ictxt_ >= 0 )
    Cblacs_gridinfo(ictxt_, &nprow_, &npcol_, &myrow_, &mycol_);
  
  size_ = nprow_ * npcol_;
  myproc_ = Cblacs_pnum(ictxt_,myrow_,mycol_);
  onpe0_ = ( mype_ == 0 );
  coutpe_ = 0;
  oncoutpe_ = onpe0_;
  active_ = ( ictxt_ >= 0 );

  pmap_.resize(size_);
  // column-major order
  int i = 0;
  for ( int ic = 0; ic < npcol; ic++ )
    for ( int ir = 0; ir < nprow; ir++ )
    {
      pmap_[ir+nprow*ic] = i;
      i++;
    }
}

#endif

////////////////////////////////////////////////////////////////////////////////
ContextRep::~ContextRep() {
  if ( myrow_ != -1 ) {
    Cblacs_gridexit( ictxt_ );
#if USE_MPI
    if (comm_ != MPI_COMM_NULL)
      MPI_Comm_free(&comm_);
#endif
  }
}

////////////////////////////////////////////////////////////////////////////////
void ContextRep::print(ostream& os) const
{
  if ( active_ )
  {
    os << " " << nprow_ << "x" << npcol_ << " ";
    os << "{ ";
    for ( int ir = 0; ir < nprow_; ir++ )
    {
      os << "{ ";
      for ( int ic = 0; ic < npcol_; ic++ )
      {
        os << pmap_[ir+nprow_*ic] << " ";
      }
      os << "} ";
    }
    //os << "} (ictxt=" << ictxt_ << ")" << endl;
    os << "} (ictxt=" << ictxt_ << ")";
    int handle;
    Cblacs_get(ictxt_,10,&handle);
    os << "(handle=" << handle << ")";
    os << "(comm=" << comm_<< ")" << endl;
  }
  else
  {
    os << " (inactive)" << endl;
  }
}

////////////////////////////////////////////////////////////////////////////////
class ContextImpl
{
  private:
  
  ContextRep* rep;
  int* pcount;
  
  public:
  
  ContextRep* operator->() { return rep; }
  ContextRep* get_rep(void) { return rep; }
  ContextImpl(ContextRep* pp) : rep(pp), pcount(new int(1)) {}
  ContextImpl(const ContextImpl& r) : rep(r.rep), pcount(r.pcount)
  { (*pcount)++; }
  
  ContextImpl& operator=(const ContextImpl& r)
  {
    if ( rep == r.rep ) return *this;
    if ( --(*pcount) == 0 )
    {
      delete rep;
      delete pcount;
    }
    rep = r.rep;
    pcount = r.pcount;
    (*pcount)++;
    return *this;
  }
  
  ~ContextImpl(void)
  {
    if ( --(*pcount) == 0 )
    {
      delete rep;
      delete pcount;
    }
  }
};

////////////////////////////////////////////////////////////////////////////////
Context::Context(void) : pimpl_(new ContextImpl(new ContextRep())) {}

////////////////////////////////////////////////////////////////////////////////
Context::Context(int nprow, int npcol):  
pimpl_(new ContextImpl(new ContextRep(nprow,npcol))) {}
 
////////////////////////////////////////////////////////////////////////////////
Context::Context(const Context &c, int nprow, int npcol, 
        int irstart, int icstart) : 
pimpl_(new ContextImpl(new ContextRep(*(c.pimpl_)->get_rep(),nprow,npcol,
irstart,icstart))) { 
 (*pimpl_)->set_coutpe(c.coutpe()); 
}
 
////////////////////////////////////////////////////////////////////////////////
Context::Context(const Context &c, int nprow, int npcol) : 
pimpl_(new ContextImpl(new ContextRep(*(c.pimpl_)->get_rep(),nprow,npcol))) {
 (*pimpl_)->set_coutpe(c.coutpe()); 
}
 
////////////////////////////////////////////////////////////////////////////////
#if USE_MPI
Context::Context(const MPI_Comm &comm, int nprow, int npcol) : 
pimpl_(new ContextImpl(new ContextRep(comm,nprow,npcol))) {
}
#endif 
////////////////////////////////////////////////////////////////////////////////
Context::~Context() { delete pimpl_; }

////////////////////////////////////////////////////////////////////////////////
Context::Context(const Context& c) : pimpl_(new ContextImpl(*(c.pimpl_))) {
 (*pimpl_)->set_coutpe(c.coutpe()); 
}

////////////////////////////////////////////////////////////////////////////////
Context& Context::operator=(const Context& rhs)
{
  if ( &rhs == this ) return *this;
  *pimpl_ = *(rhs.pimpl_);
  return *this;
}

////////////////////////////////////////////////////////////////////////////////
int Context::ictxt() const { return (*pimpl_)->ictxt(); }

////////////////////////////////////////////////////////////////////////////////
int Context::myrow() const { return (*pimpl_)->myrow(); }

////////////////////////////////////////////////////////////////////////////////
int Context::mycol() const { return (*pimpl_)->mycol(); }

////////////////////////////////////////////////////////////////////////////////
int Context::nprow() const { return (*pimpl_)->nprow(); }

////////////////////////////////////////////////////////////////////////////////
int Context::npcol() const { return (*pimpl_)->npcol(); }

////////////////////////////////////////////////////////////////////////////////
int Context::size() const { return (*pimpl_)->size(); }

////////////////////////////////////////////////////////////////////////////////
int Context::myproc() const { return (*pimpl_)->myproc(); }

////////////////////////////////////////////////////////////////////////////////
int Context::mype() const { return (*pimpl_)->mype(); }

////////////////////////////////////////////////////////////////////////////////
int Context::pmap(int irow, int icol) const
{ return (*pimpl_)->pmap(irow,icol); }

////////////////////////////////////////////////////////////////////////////////
bool Context::oncoutpe(void) const { return (*pimpl_)->oncoutpe(); }

////////////////////////////////////////////////////////////////////////////////
int Context::coutpe(void) const { return (*pimpl_)->coutpe(); }

////////////////////////////////////////////////////////////////////////////////
void Context::set_coutpe(int num) { (*pimpl_)->set_coutpe(num); }

////////////////////////////////////////////////////////////////////////////////
bool Context::onpe0(void) const { return (*pimpl_)->onpe0(); }

////////////////////////////////////////////////////////////////////////////////
bool Context::active(void) const { return (*pimpl_)->active(); }

////////////////////////////////////////////////////////////////////////////////
void Context::abort(int ierr) const { (*pimpl_)->abort(ierr); }

////////////////////////////////////////////////////////////////////////////////
void Context::barrier(void) const { (*pimpl_)->barrier(); }

////////////////////////////////////////////////////////////////////////////////
void Context::barrier(char scope) const { (*pimpl_)->barrier(scope); }

////////////////////////////////////////////////////////////////////////////////
void Context::dsend(int m, int n, double* a, 
  int lda, int rdest, int cdest) const
{ (*pimpl_)->dsend(m,n,a,lda,rdest,cdest); }

void Context::drecv(int m, int n, double* a, 
  int lda, int rsrc, int csrc) const
{ (*pimpl_)->drecv(m,n,a,lda,rsrc,csrc); }

void Context::dsum(char scope, char topology, 
  int m, int n, double* a, int lda, int rdest, int cdest) const
{ (*pimpl_)->dsum(scope,topology,m,n,a,lda,rdest,cdest); }

void Context::dsum(char scope, int m, int n, double* a, int lda) const
{ (*pimpl_)->dsum(scope,' ',m,n,a,lda,-1,-1); }

void Context::dsum(int m, int n, double* a, int lda) const
{ (*pimpl_)->dsum('A',' ',m,n,a,lda,-1,-1); }

void Context::dmax(char scope, char topology, 
  int m, int n, double* a, int lda, int* ra, int* ca, int rcflag, 
  int rdest, int cdest) const
{ (*pimpl_)->dmax(scope,topology,m,n,a,lda,ra,ca,rcflag,rdest,cdest); }

void Context::dmax(char scope, char topology, 
  int m, int n, double* a, int lda, int rdest, int cdest) const
{ (*pimpl_)->dmax(scope,topology,m,n,a,lda,rdest,cdest); }

void Context::dmax(char scope, int m, int n, double* a, int lda) const
{ (*pimpl_)->dmax(scope,' ',m,n,a,lda,-1,-1); }

void Context::dmax(int m, int n, double* a, int lda) const
{ (*pimpl_)->dmax('A',' ',m,n,a,lda,-1,-1); }

void Context::dmin(char scope, char topology, 
  int m, int n, double* a, int lda, int* ra, int* ca, int rcflag, 
  int rdest, int cdest) const
{ (*pimpl_)->dmin(scope,topology,m,n,a,lda,ra,ca,rcflag,rdest,cdest); }

void Context::dmin(char scope, char topology, 
  int m, int n, double* a, int lda, int rdest, int cdest) const
{ (*pimpl_)->dmin(scope,topology,m,n,a,lda,rdest,cdest); }

void Context::dmin(char scope, int m, int n, double* a, int lda) const
{ (*pimpl_)->dmin(scope,' ',m,n,a,lda,-1,-1); }

void Context::dmin(int m, int n, double* a, int lda) const
{ (*pimpl_)->dmin('A',' ',m,n,a,lda,-1,-1); }

void Context::dbcast_send(char scope, char topology, 
  int m, int n, double* a, int lda) const
{ (*pimpl_)->dbcast_send(scope,topology,m,n,a,lda); }

void Context::dbcast_send(char scope, int m, int n, double* a, int lda) const
{ (*pimpl_)->dbcast_send(scope,' ',m,n,a,lda); }

void Context::dbcast_send(int m, int n, double* a, int lda) const
{ (*pimpl_)->dbcast_send('A',' ',m,n,a,lda); }

void Context::dbcast_recv(char scope, char topology, 
  int m, int n, double* a, int lda, int rsrc,int csrc) const
{ (*pimpl_)->dbcast_recv(scope,topology,m,n,a,lda,rsrc,csrc); }

void Context::dbcast_recv(char scope, 
  int m, int n, double* a, int lda, int rsrc, int csrc) const
{ (*pimpl_)->dbcast_recv(scope,' ',m,n,a,lda,rsrc,csrc); }

void Context::dbcast_recv(int m, int n, double* a, int lda, 
  int rsrc,int csrc) const
{ (*pimpl_)->dbcast_recv('A',' ',m,n,a,lda,rsrc,csrc); }

////////////////////////////////////////////////////////////////////////////////
void Context::isend(int m, int n, int* a, 
  int lda, int rdest, int cdest) const
{ (*pimpl_)->isend(m,n,a,lda,rdest,cdest); }

void Context::irecv(int m, int n, int* a, 
  int lda, int rsrc, int csrc) const
{ (*pimpl_)->irecv(m,n,a,lda,rsrc,csrc); }

void Context::isum(char scope, char topology, 
  int m, int n, int* a, int lda, int rdest, int cdest) const
{ (*pimpl_)->isum(scope,topology,m,n,a,lda,rdest,cdest); }

void Context::isum(char scope, int m, int n, int* a, int lda) const
{ (*pimpl_)->isum(scope,' ',m,n,a,lda,-1,-1); }

void Context::isum(int m, int n, int* a, int lda) const
{ (*pimpl_)->isum('A',' ',m,n,a,lda,-1,-1); }

void Context::imax(char scope, char topology, 
  int m, int n, int* a, int lda, int* ra, int* ca, int rcflag, 
  int rdest, int cdest) const
{ (*pimpl_)->imax(scope,topology,m,n,a,lda,ra,ca,rcflag,rdest,cdest); }

void Context::imax(char scope, char topology, 
  int m, int n, int* a, int lda, int rdest, int cdest) const
{ (*pimpl_)->imax(scope,topology,m,n,a,lda,rdest,cdest); }

void Context::imax(char scope, int m, int n, int* a, int lda) const
{ (*pimpl_)->imax(scope,' ',m,n,a,lda,-1,-1); }

void Context::imax(int m, int n, int* a, int lda) const
{ (*pimpl_)->imax('A',' ',m,n,a,lda,-1,-1); }

void Context::imin(char scope, char topology, 
  int m, int n, int* a, int lda, int* ra, int* ca, int rcflag, 
  int rdest, int cdest) const
{ (*pimpl_)->imin(scope,topology,m,n,a,lda,ra,ca,rcflag,rdest,cdest); }

void Context::imin(char scope, char topology, 
  int m, int n, int* a, int lda, int rdest, int cdest) const
{ (*pimpl_)->imin(scope,topology,m,n,a,lda,rdest,cdest); }

void Context::imin(char scope, int m, int n, int* a, int lda) const
{ (*pimpl_)->imin(scope,' ',m,n,a,lda,-1,-1); }

void Context::imin(int m, int n, int* a, int lda) const
{ (*pimpl_)->imin('A',' ',m,n,a,lda,-1,-1); }

void Context::ibcast_send(char scope, char topology, 
  int m, int n, int* a, int lda) const
{ (*pimpl_)->ibcast_send(scope,topology,m,n,a,lda); }

void Context::ibcast_send(char scope, int m, int n, int* a, int lda) const
{ (*pimpl_)->ibcast_send(scope,' ',m,n,a,lda); }

void Context::ibcast_send(int m, int n, int* a, int lda) const
{ (*pimpl_)->ibcast_send('A',' ',m,n,a,lda); }

void Context::ibcast_recv(char scope, char topology, 
  int m, int n, int* a, int lda, int rsrc,int csrc) const
{ (*pimpl_)->ibcast_recv(scope,topology,m,n,a,lda,rsrc,csrc); }

void Context::ibcast_recv(char scope, 
  int m, int n, int* a, int lda, int rsrc, int csrc) const
{ (*pimpl_)->ibcast_recv(scope,' ',m,n,a,lda,rsrc,csrc); }

void Context::ibcast_recv(int m, int n, int* a, int lda, 
  int rsrc,int csrc) const
{ (*pimpl_)->ibcast_recv('A',' ',m,n,a,lda,rsrc,csrc); }

void Context::string_send(string& s, int rdest, int cdest) const
{ (*pimpl_)->string_send(s,rdest,cdest); }

void Context::string_recv(string& s, int rsrc, int csrc) const
{ (*pimpl_)->string_recv(s,rsrc,csrc); }

void Context::string_bcast(string& s, int isrc) const
{ (*pimpl_)->string_bcast(s,isrc); }

////////////////////////////////////////////////////////////////////////////////
bool Context::operator==(const Context& ctxt) const
{ return ( (*pimpl_)->ictxt() == ctxt.ictxt() ); }

////////////////////////////////////////////////////////////////////////////////
MPI_Comm Context::comm(void) const { return (*pimpl_)->comm(); }
 
////////////////////////////////////////////////////////////////////////////////
void Context::print(ostream& os) const { (*pimpl_)->print(os);}

////////////////////////////////////////////////////////////////////////////////
ostream& operator<<(ostream& os, const Context& c)
{
  c.print(os);
  return os;
}

