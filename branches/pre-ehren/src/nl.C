////////////////////////////////////////////////////////////////////////////////
double NonLocalPotential::energy(bool compute_hpsi, SlaterDet& dsd, 
    bool compute_forces, vector<vector<double> >& fion, 
    bool compute_stress, valarray<double>& sigma_enl)
{
  const bool compute_anl = false;

  const vector<double>& occ = sd_.occ();
  const int ngwl = basis_.localsize();
  // define atom block size
  const int na_block_size = 32;
  valarray<double> gr(na_block_size*ngwl); // gr[ig+ia*ngwl]
  valarray<double> cgr(na_block_size*ngwl); // cgr[ig+ia*ngwl]
  valarray<double> sgr(na_block_size*ngwl); // sgr[ig+ia*ngwl]
  vector<vector<double> > tau;
  atoms_.get_positions(tau);

  double enl = 0.0;
  double tsum[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  if ( nspnl == 0 ) return 0.0;
  const double omega = basis_.cell().volume();
  assert(omega != 0.0);
  const double omega_inv = 1.0 / omega;
  
  for ( int is = 0; is < nsp; is++ )   
  {  
    if ( npr[is] > 0 ) // species is is non-local
    {
      if ( compute_anl )
      {
        // define number of atom blocks
        const int na_blocks = na[is] / na_block_size + 
                              ( na[is] % na_block_size == 0 ? 0 : 1 );
                              
        valarray<double> anl_loc(npr[is]*na_block_size*2*ngwl);
        const int nstloc = sd_.nstloc();
        // fnl_loc[ipra][n]
        valarray<double> fnl_loc(npr[is]*na_block_size*nstloc);
        valarray<double> fnl_buf(npr[is]*na_block_size*nstloc);
        for ( int ia_block = 0; ia_block < na_blocks; ia_block++ )
        {
          // process projectors of atoms in block ia_block
          
          const int iastart = ia_block * na_block_size;
          const int iaend = (ia_block+1) * na_block_size < na[is] ?
                          (ia_block+1) * na_block_size :
                          na[is];
          const int ia_block_size = iaend - iastart; 
          
          // compute cgr[is][ia][ig], sgr[is][ia][ig]
          int k = 3;
          double mone = -1.0, zero = 0.0;
          char cn='n';
        
          // next line: const cast is ok since dgemm_ does not modify argument
          double* gx = const_cast<double*>(basis_.gx_ptr(0));
          dgemm(&cn,&cn,(int*)&ngwl,(int*)&ia_block_size,&k,&mone,
                gx,(int*)&ngwl, &tau[is][3*iastart],&k,
                &zero,&gr[0],(int*)&ngwl);

          int len = ia_block_size * ngwl;
#if AIX || BGL
          vsincos(&sgr[is][0],&cgr[is][0],&gr[0],&len);
#else
          for ( int i = 0; i < len; i++ )
          {
            const double arg = gr[i];
            sgr[i] = sin(arg);
            cgr[i] = cos(arg);
          }
#endif

          // compute anl_loc
          for ( int ipr = 0; ipr < npr[is]; ipr++ )
          {
            // twnl[is][ig+ngwl*ipr]
            const double * t = &twnl[is][ngwl*ipr];
            const int l = lproj[is][ipr];

            // anl_loc[ig+ipra*ngwl]
            double * a = &anl_loc[ipr*ia_block_size*ngwl];

            if ( l == 0 )
            {
              for ( int ia = 0; ia < ia_block_size; ia++ )
              {
                for ( int ig = 0; ig < ngwl; ig++ )
                {
                  a[ig+ia*ngwl]   = t[ig] * cgr[ig+ia*ngwl];
                  a[ig+1+ia*ngwl] = t[ig] * sgr[ig+ia*ngwl];
                }
              }
            }
            else if ( l == 1 )
            {
              for ( int ia = 0; ia < ia_block_size; ia++ )
              {
                for ( int ig = 0; ig < ngwl; ig++ )
                {
                  /* Next line: -i * eigr */
                  /* -i * (a+i*b) = b - i*a */
                  a[ig+ia*ngwl]   =  t[ig] * sgr[ig+ia*ngwl];
                  a[ig+1+ia*ngwl] = -t[ig] * cgr[ig+ia*ngwl];
                }
              }
            }
            else if ( l == 2 )
            {
              for ( int ia = 0; ia < ia_block_size; ia++ )
              {
                for ( int ig = 0; ig < ngwl; ig++ )
                {
                  // Next line: (-) sign for -eigr
                  a[ig+ia*ngwl]   = -t[ig] * cgr[ig+ia*ngwl];
                  a[ig+1+ia*ngwl] = -t[ig] * sgr[ig+ia*ngwl];
                }
              }
            }
          } // ipr
          
          // array anl_loc is complete
          
          // compute fnl[npra][nstloc] = anl^T * c
          double one=1.0;
          char ct='t';
          int twongwl = 2 * ngwl;
          int nprnaloc = ia_block_size * npr[is];
          const complex<double>* c = sd_.c().cvalptr();
          dgemm(&ct,&cn,&nprnaloc,(int*)&nstloc,&twongwl,&one,
                &anl_loc[0],&twongwl, (double*)c, &twongwl, 
                &zero,&fnl_loc[0],&nprnaloc);
                
          // correct for double counting if ctxt_.myrow() == 0
          if ( ctxt_.myrow() == 0 )
          {
            // rank-one update
            // dger(m,n,alpha,x,incx,y,incy,a,lda);
            // a += alpha * x * transpose(y)
            // x = first row of anl_loc
            // y^T = first row of c
            double alpha = -0.5;
            dger(&nprnaloc,(int*)&nstloc,&alpha,&anl_loc[0],&twongwl,
                 (double*)c,&twongwl,&fnl_loc[0],&nprnaloc);
          }
          
          // Allreduce fnl partial sum
          MPI_Comm basis_comm = basis_.context().comm();
          double fnl_size = nprnaloc*nstloc;
          MPI_Allreduce(&fnl_loc[0],&fnl_buf[0],fnl_size,
                        MPI_DOUBLE,MPI_SUM,basis_comm);
                        
          // factor 2.0 in next line is: counting G, -G
          fnl_loc = 2.0 * fnl_buf;
          
          // accumulate Enl contribution
          const int nbase = ctxt_.mycol() * sd_.c().nb();
          for ( int ipr = 0; ipr < npr[is]; ipr++ )
          {
            const double fac = wt[is][ipr] * omega_inv;
            for ( int n = 0; n < nstloc; n++ )
            {
              const double facn = fac * occ[n + nbase];
              for ( int ia = 0; ia < ia_block_size; ia++ )
              {
                const int i = ia + ipr*ia_block_size + n * nprnaloc;
                cout << "fnl_loc[ipr=" << ipr << ",ia=" << ia
                     << ",n=" << n << "]: " << fnl_loc[i] << endl;
                const double tmp = fnl_loc[i];
                enl += facn * tmp * tmp;
                fnl_loc[i] = fac * tmp;
              }
            }
          }
          
          if ( compute_hpsi )
          {
            // compute cp += anl * fnl
            complex<double>* cp = dsd.c().valptr();
            dgemm(&cn,&cn,&twongwl,(int*)&nstloc,&nprnaloc,&one,
                  &anl_loc[0],&twongwl, &fnl_loc[0],&nprnaloc,
                  &one,(double*)cp, &twongwl);
          }
          
          assert(compute_forces==false);
          assert(compute_stress==false);
          
        } // ia_block
      }
      else
      {
      // compute fnl
      // block distribution for fnl: same as SlaterDet for nst
      DoubleMatrix fnl(ctxt_,anl[is]->n(),sd_.c().n(),
                       anl[is]->nb(),sd_.c().nb());
      
      const DoubleMatrix c_proxy(sd_.c());
      
      tmap["fnl_gemm"].start();
      fnl.gemm('t','n',2.0,*anl[is],c_proxy,0.0);
      tmap["fnl_gemm"].stop();
      
      // correct for double counting of G=0 components
      // rank-1 update using first row of *anl[is] and c_proxy
      fnl.ger(-1.0,*anl[is],0,c_proxy,0);
      
      cout << fnl << endl;
      
      // compute the non-local energy
      // multiply fnl[ipra+nprna*n] by fac = wt[is][ipr] * omega_inv;
      // block sizes: npr*nalocmax x c().nb()
      // loop over local array
      double*f = fnl.valptr(0);
      const int mb = fnl.mb();
      const int nb = fnl.nb();
      const int mloc = fnl.mloc();
      for ( int li=0; li < fnl.mblocks(); li++)
      {
        const int mbs = fnl.mbs(li);
        for ( int lj=0; lj < fnl.nblocks(); lj++)
        {
          const int nbs = fnl.nbs(lj);
          for ( int ii=0; ii < mbs; ii++)
          {
            assert(mbs%npr[is]==0);
            // mbs/npr[is] is the number of atoms in the block li
            const int ipr = ii / (mbs/npr[is]);
            const double fac = wt[is][ipr] * omega_inv;
            for ( int jj=0; jj < nbs; jj++)
            {
              // global index: i(li,ii), j(lj,jj)
              const int nglobal = fnl.j(lj,jj);
              const double facn = fac * occ[nglobal];
              const int iii = ii+li*mb;
              const int jjj = jj+lj*nb;
              const double tmp = f[iii+mloc*jjj];
              enl += facn * tmp * tmp;
              f[iii+mloc*jjj] = fac * tmp;
            }
          }
        }
      }

      if ( compute_hpsi )
      {
        tmap["enl_hpsi"].start();
        // Apply operator to electronic states and accumulate in dsd
        DoubleMatrix cp_proxy(dsd.c());
        cp_proxy.gemm('n','n',1.0,*anl[is],fnl,1.0);
        tmap["enl_hpsi"].stop();
      }

      // ionic forces
      if ( compute_forces )
      {
        tmap["enl_fion"].start();
        double *tmpfion = new double[3*na[is]];
        for ( int i = 0; i < 3*na[is]; i++ )
          tmpfion[i] = 0.0;
          
        DoubleMatrix danl(ctxt_,anl[is]->m(),anl[is]->n(),
          anl[is]->mb(),anl[is]->nb());
        DoubleMatrix dfnl(ctxt_,fnl.m(),fnl.n(),fnl.mb(),fnl.nb());
        const int ngwl = basis_.localsize();
          
        for ( int j = 0; j < 3; j++ )
        {
          const double *const gxj = basis_.gx_ptr(j);
          for ( int ipr = 0; ipr < npr[is]; ipr++ )
          {
            const int l = lproj[is][ipr];
 
            // twnl[is][ig+ngwl*ipr]
            const double *t = &twnl[is][ngwl*ipr];
            for ( int ia = 0; ia < naloc[is]; ia++ )
            {
              // danl[ig+ipra*ngwl]
              // index = ig+cmloc_anl*(ia+nais*ipr), ig=0
              const int ipra = ia+naloc[is]*ipr;
              double *da = danl.valptr(2*(sd_.c().mloc()*ipra));
              const double *c = &cosgr[is][ia*ngwl];
              const double *s = &singr[is][ia*ngwl];
 
              if ( l == 0 )
              {
                for ( int ig = 0; ig < ngwl; ig++ )
                {
                  const double tt = gxj[ig] * t[ig];
                  // Next lines: -i * ( a + ib ) = b - ia
                  *da++ =  tt * *s++;
                  *da++ = -tt * *c++;
                }
              }
              else if ( l == 1 )
              {
                for ( int ig = 0; ig < ngwl; ig++ )
                {
                  // Next lines: (-i)**2 * ( a + ib ) = - a - ib
                  const double tt = - gxj[ig] * t[ig];
                  *da++ = tt * *c++;
                  *da++ = tt * *s++;
                }
              }
              else if ( l == 2 )
              {
                for ( int ig = 0; ig < ngwl; ig++ )
                {
                  // Next lines: (-i) * - ( a + ib ) = i*(a+ib) = - b + ia
                  const double tt = gxj[ig] * t[ig];
                  *da++ = -tt * *s++;
                  *da++ =  tt * *c++;
                }
              }
            } // ia
          } // ipr

          // compute dfnl
          const DoubleMatrix c_proxy(sd_.c());
 
          dfnl.gemm('t','n',2.0,danl,c_proxy,0.0);
 
          // Note: no need to correct for double counting of the
          // G=0 component which is always zero

          // non-local forces
 
          // loop over local array
          // block sizes: npr*nalocmax x c().nb()
          const double*f = fnl.valptr(0);
          const double*df = dfnl.valptr(0);
          const int mloc = fnl.mloc();
          const int mb = fnl.mb();
          const int nb = fnl.nb();
          for ( int li=0; li < fnl.mblocks(); li++)
          {
            // find index of first atom in block li
            const int ia_first = nalocmax[is] * 
              ( li * fnl.context().nprow() + fnl.context().myrow() );
            const int mbs = fnl.mbs(li);
            for ( int lj=0; lj < fnl.nblocks(); lj++)
            {
              const int nbs = fnl.nbs(lj);
              for ( int ii=0; ii < mbs; ii++)
              {
                // ia_local: index of atom within block li
                const int ia_local = ii % ( mbs / npr[is] );
                const int ia_global = ia_local + ia_first;
                assert(3*ia_global+j < 3*na[is]);
                for ( int jj=0; jj < nbs; jj++)
                {
                  const int nglobal = fnl.j(lj,jj);
                  // Factor 2.0 in next line from derivative of |Fnl|^2
                  const double facn = 2.0 * occ[nglobal];
                  const int iii = ii+li*mb;
                  const int jjj = jj+lj*nb;
                  tmpfion[3*ia_global+j] -= facn *
                    f[iii+mloc*jjj] * df[iii+mloc*jjj];
                }
              }
            }
          }
        } // j
        
        ctxt_.dsum(3*na[is],1,tmpfion,3*na[is]);
        for ( int ia = 0; ia < na[is]; ia++ )
        {
          fion[is][3*ia+0] += tmpfion[3*ia];
          fion[is][3*ia+1] += tmpfion[3*ia+1];
          fion[is][3*ia+2] += tmpfion[3*ia+2];
        }
        delete [] tmpfion;
        tmap["enl_fion"].stop();
      } // compute_forces
        
      if ( compute_stress )
      {
        const int ngwl = basis_.localsize();
        DoubleMatrix danl(ctxt_,anl[is]->m(),anl[is]->n(),
          anl[is]->mb(),anl[is]->nb());
        DoubleMatrix dfnl(ctxt_,fnl.m(),fnl.n(),fnl.mb(),fnl.nb());
        
        for ( int ij = 0; ij < 6; ij++ )
        {
          int ipr = 0;
          while ( ipr < npr[is] )
          {
            const int l = lproj[is][ipr];
            if ( l == 0 )
            {
              // dtwnl[is][ipr][ij][ngwl]
              // index = ig + ngwl * ( ij + 6 * ipr))
              // ipr = iquad + nquad[is] * ilm, where ilm = 0
              const double *const dt0 = &dtwnl[is][ngwl*(ij+6*ipr)];
              for ( int ia = 0; ia < naloc[is]; ia++ )
              {
                const int ipra0 = ia+naloc[is]*ipr;
                double *da0 = danl.valptr(2*(sd_.c().mloc()*ipra0));
                const double *c = &cosgr[is][ia*ngwl];
                const double *s = &singr[is][ia*ngwl];
                for ( int ig = 0; ig < ngwl; ig++ )
                {
                  const double d0 = dt0[ig];
                  // danl[is][ipr][iquad][ia][ig].re =
                  //   dtwnl[is][ipr][iquad][j][ig] * cosgr[is][ia][ig]
                  *da0++ = *c++ * d0;
                  // danl[is][ipr][iquad][ia][ig].im =
                  //   dtwnl[is][ipr][iquad][j][ig] * singr[is][ia][ig]
                  *da0++ = *s++ * d0;
                }
              }
            }
            else if ( l == 1 )
            {
              const int ipr1 = ipr;
              const int ipr2 = ipr + 1;
              const int ipr3 = ipr + 2;
              // dtwnl[is][ipr][ij][ngwl]
              // index = ig + ngwl * ( ij + 6 * iprx ))
              const double *dt1 = &dtwnl[is][ngwl*(ij+6*ipr1)];
              const double *dt2 = &dtwnl[is][ngwl*(ij+6*ipr2)];
              const double *dt3 = &dtwnl[is][ngwl*(ij+6*ipr3)];
              for ( int ia = 0; ia < naloc[is]; ia++ )
              {
                const int ipra1 = ia+naloc[is]*ipr1;
                const int ipra2 = ia+naloc[is]*ipr2;
                const int ipra3 = ia+naloc[is]*ipr3;
                double *da1 = danl.valptr(2*(sd_.c().mloc()*ipra1));
                double *da2 = danl.valptr(2*(sd_.c().mloc()*ipra2));
                double *da3 = danl.valptr(2*(sd_.c().mloc()*ipra3));
 
                const double *c = &cosgr[is][ia*ngwl];
                const double *s = &singr[is][ia*ngwl];
                for ( int ig = 0; ig < ngwl; ig++ )
                {
                  const double d1 = dt1[ig];
                  const double d2 = dt2[ig];
                  const double d3 = dt3[ig];
                  // Next line: (-i)^l factor is -i
                  // Next line: -i * eigr
                  // -i * (a+i*b) = b - i*a
                  const double tc = -*c++; //  -cosgr[is][ia][ig]
                  const double ts =  *s++; //   singr[is][ia][ig]
                  *da1++ = d1 * ts;
                  *da1++ = d1 * tc;
                  *da2++ = d2 * ts;
                  *da2++ = d2 * tc;
                  *da3++ = d3 * ts;
                  *da3++ = d3 * tc;
                }
              }
            }
            else if ( l == 2 )
            {
              const int ipr4 = ipr;                                            
              const int ipr5 = ipr + 1;                                        
              const int ipr6 = ipr + 2;                                        
              const int ipr7 = ipr + 3;                                        
              const int ipr8 = ipr + 4;                                        
              // dtwnl[is][ipr][iquad][ij][ngwl]
              // index = ig + ngwl * ( ij + 6 * ( iquad + nquad[is] * ipr ))
              const double *dt4 = &dtwnl[is][ngwl*(ij+6*ipr4)];
              const double *dt5 = &dtwnl[is][ngwl*(ij+6*ipr5)];
              const double *dt6 = &dtwnl[is][ngwl*(ij+6*ipr6)];
              const double *dt7 = &dtwnl[is][ngwl*(ij+6*ipr7)];
              const double *dt8 = &dtwnl[is][ngwl*(ij+6*ipr8)];
              for ( int ia = 0; ia < naloc[is]; ia++ )
              {
                const int ipra4 = ia+naloc[is]*ipr4;
                const int ipra5 = ia+naloc[is]*ipr5;
                const int ipra6 = ia+naloc[is]*ipr6;
                const int ipra7 = ia+naloc[is]*ipr7;
                const int ipra8 = ia+naloc[is]*ipr8;
                double *da4 = danl.valptr(2*(sd_.c().mloc()*ipra4));
                double *da5 = danl.valptr(2*(sd_.c().mloc()*ipra5));
                double *da6 = danl.valptr(2*(sd_.c().mloc()*ipra6));
                double *da7 = danl.valptr(2*(sd_.c().mloc()*ipra7));
                double *da8 = danl.valptr(2*(sd_.c().mloc()*ipra8));
 
                const double *c = &cosgr[is][ia*ngwl];
                const double *s = &singr[is][ia*ngwl];
                for ( int ig = 0; ig < ngwl; ig++ )
                {
                  const double d4 = dt4[ig];
                  const double d5 = dt5[ig];
                  const double d6 = dt6[ig];
                  const double d7 = dt7[ig];
                  const double d8 = dt8[ig];
                  // Next lines: (-i)^2 * ( a + ib ) =  - ( a + ib )
                  const double tc = -*c++;
                  const double ts = -*s++;
                  *da4++ = d4 * tc;
                  *da4++ = d4 * ts;
                  *da5++ = d5 * tc;
                  *da5++ = d5 * ts;
                  *da6++ = d6 * tc;
                  *da6++ = d6 * ts;
                  *da7++ = d7 * tc;
                  *da7++ = d7 * ts;
                  *da8++ = d8 * tc;
                  *da8++ = d8 * ts;
                }
              }
            }                                                                
            else                                                               
            {                                                                  
              assert(false);                                                   
            } // l
            ipr += 2*l+1;
          } // while ipr
      
          // compute dfnl
          const DoubleMatrix c_proxy(sd_.c());
 
          dfnl.gemm('t','n',2.0,danl,c_proxy,0.0);
 
          // Note: no need to correct for double counting of the
          // G=0 component which is always zero
          
           // partial contributions to the stress sigma_ij
          // Note: fnl was already premultiplied by the factor
          // fac = wt[is][ipr][iquad] * omega_inv;
          const double *const f = fnl.cvalptr(0);
          const double *const df = dfnl.cvalptr(0);
          const int mb = fnl.mb();
          const int nb = fnl.nb();
          const int mloc = fnl.mloc();
          for ( int li=0; li < fnl.mblocks(); li++)
          {
            const int mbs = fnl.mbs(li);
            for ( int lj=0; lj < fnl.nblocks(); lj++)
            {
              const int nbs = fnl.nbs(lj);
              for ( int ii=0; ii < mbs; ii++)
              {
                for ( int jj=0; jj < nbs; jj++)
                {
                  // global index: i(li,ii), j(lj,jj)
                  const int nglobal = fnl.j(lj,jj);
                  const double facn = 2.0 * occ[nglobal];
                  const int iii = ii+li*mb;
                  const int jjj = jj+lj*nb;
                  const double tmp = f[iii+mloc*jjj];
                  const double dtmp = df[iii+mloc*jjj];
                  tsum[ij] += facn * tmp * dtmp;
                }
              }
            }
          }
        } // ij
      } // compute_stress
      } // compute_anl
    } // npr[is]>0
  } // is

  ctxt_.dsum(1,1,&enl,1);

  sigma_enl = 0.0;
  if ( compute_stress )
  {
    ctxt_.dsum(6,1,&tsum[0],6);
    sigma_enl[0] = ( enl + tsum[0] ) * omega_inv;
    sigma_enl[1] = ( enl + tsum[1] ) * omega_inv;
    sigma_enl[2] = ( enl + tsum[2] ) * omega_inv;
    sigma_enl[3] = + tsum[3] * omega_inv;
    sigma_enl[4] = + tsum[4] * omega_inv;
    sigma_enl[5] = + tsum[5] * omega_inv;
  }

  return enl;
}

