////////////////////////////////////////////////////////////////////////////////
//
// blacs.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: blacs.h,v 1.2 2008/07/02 17:06:51 draeger1 Exp $

#ifndef BLACS_H
#define BLACS_H

extern "C"{
void igesd2d(int*,int*,int*, int*, int*,int*,int*);
void sgesd2d(int*,int*,int*, double*, int*,int*,int*);
void igerv2d(int*,int*,int*, int*, int*,int*,int*);
void sgerv2d(int*,int*,int*, double*, int*,int*,int*);
void sgsum2d(int*,char*,char*,int*,int*,double*,int*,int*,int*);
void igamn2d(int*,char*,char*, int*, 
     int*,int*,int*, int*, int*,int*,int*, int*);
void blacs_pinfo(int*, int*);
void blacs_get(int*, int*, int*);
void blacs_barrier(int*, char*);
void blacs_gridinfo(int*, int *, int *, int *, int *);
void blacs_gridinit(int *, char*, int*, int*);
void blacs_gridmap(int*, int *, int*, int*, int*);
void blacs_abort(int*, int*);
void blacs_gridexit(int*);
int  blacs_pnum(int*, int*, int*);
int sys2blacs_handle(int);
}


#ifdef SCALAPACK
extern "C"{
#endif
// C interface to the BLACS
void Cdgesd2d(int,int,int, double*, int,int,int);
void Cdgerv2d(int,int,int, double*, int,int,int);
void Cdgsum2d(int,char*,char*,int,int,double*,int,int,int);
void Cdgamx2d(int,char*,char*,int,int,double*,int,int*,int*,int,int,int);
void Cdgamn2d(int,char*,char*,int,int,double*,int,int*,int*,int,int,int);
void Cdgebs2d(int,char*,char*,int,int,double*,int);
void Cdgebr2d(int,char*,char*,int,int,double*,int,int,int);

void Cigesd2d(int,int,int, int*, int,int,int);
void Cigerv2d(int,int,int, int*, int,int,int);
void Cigsum2d(int,char*,char*,int,int,int*,int,int,int);
void Cigamx2d(int,char*,char*,int,int,int*,int,int*,int*,int,int,int);
void Cigamn2d(int,char*,char*,int,int,int*,int,int*,int*,int,int,int);
void Cigebs2d(int,char*,char*,int,int,int*,int);
void Cigebr2d(int,char*,char*,int,int,int*,int,int,int);

void Cblacs_pinfo(int*, int*);
void Cblacs_get(int, int, int*);
void Cblacs_barrier(int, char*);
void Cblacs_gridinfo(int, int*, int*, int*, int*);
void Cblacs_gridinit(int*, char [], int, int);
void Cblacs_gridmap(int*, int*, int, int, int);
void Cblacs_abort(int, int);
void Cblacs_gridexit(int);
int Cblacs_pnum(int, int, int);
int Csys2blacs_handle(int);

#ifdef SCALAPACK
}
#endif

#endif
