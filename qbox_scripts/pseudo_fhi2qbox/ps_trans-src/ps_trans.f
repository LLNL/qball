MODULE nrtype
  INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
  INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)
  INTEGER, PARAMETER :: SP = KIND(1.0)
  INTEGER, PARAMETER :: DP = KIND(1.0D0)
  INTEGER, PARAMETER :: SPC = KIND((1.0,1.0))
  INTEGER, PARAMETER :: DPC = KIND((1.0D0,1.0D0))
  INTEGER, PARAMETER :: LGT = KIND(.true.)
  REAL(SP), PARAMETER :: PI=3.141592653589793238462643383279502884197_sp
  REAL(SP), PARAMETER :: PIO2=1.57079632679489661923132169163975144209858_sp
  REAL(SP), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_sp
  REAL(SP), PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_sp
  REAL(SP), PARAMETER :: EULER=0.5772156649015328606065120900824024310422_sp
  REAL(DP), PARAMETER :: PI_D=3.141592653589793238462643383279502884197_dp
  REAL(DP), PARAMETER :: PIO2_D=1.57079632679489661923132169163975144209858_dp
  REAL(DP), PARAMETER :: TWOPI_D=6.283185307179586476925286766559005768394_dp
  TYPE sprs2_sp
     INTEGER(I4B) :: n,len
     REAL(SP), DIMENSION(:), POINTER :: val
     INTEGER(I4B), DIMENSION(:), POINTER :: irow
     INTEGER(I4B), DIMENSION(:), POINTER :: jcol
  END TYPE sprs2_sp
  TYPE sprs2_dp
     INTEGER(I4B) :: n,len
     REAL(DP), DIMENSION(:), POINTER :: val
     INTEGER(I4B), DIMENSION(:), POINTER :: irow
     INTEGER(I4B), DIMENSION(:), POINTER :: jcol
  END TYPE sprs2_dp
END MODULE nrtype
MODULE nrutil
  USE nrtype
  IMPLICIT NONE
  INTEGER(I4B), PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8
  INTEGER(I4B), PARAMETER :: NPAR_GEOP=4,NPAR2_GEOP=2
  INTEGER(I4B), PARAMETER :: NPAR_CUMSUM=16
  INTEGER(I4B), PARAMETER :: NPAR_CUMPROD=8
  INTEGER(I4B), PARAMETER :: NPAR_POLY=8
  INTEGER(I4B), PARAMETER :: NPAR_POLYTERM=8
  INTERFACE array_copy
     MODULE PROCEDURE array_copy_r, array_copy_d, array_copy_i
  END INTERFACE
  INTERFACE swap
     MODULE PROCEDURE swap_i,swap_r,swap_rv,swap_c, &
          swap_cv,swap_cm,swap_z,swap_zv,swap_zm, &
          masked_swap_rs,masked_swap_rv,masked_swap_rm
  END INTERFACE
  INTERFACE reallocate
     MODULE PROCEDURE reallocate_rv,reallocate_rm,&
          reallocate_iv,reallocate_im,reallocate_hv
  END INTERFACE
  INTERFACE imaxloc
     MODULE PROCEDURE imaxloc_r,imaxloc_i
  END INTERFACE
  INTERFACE assert
     MODULE PROCEDURE assert1,assert2,assert3,assert4,assert_v
  END INTERFACE
  INTERFACE assert_eq
     MODULE PROCEDURE assert_eq2,assert_eq3,assert_eq4,assert_eqn
  END INTERFACE
  INTERFACE arth
     MODULE PROCEDURE arth_r, arth_d, arth_i
  END INTERFACE
  INTERFACE geop
     MODULE PROCEDURE geop_r, geop_d, geop_i, geop_c, geop_dv
  END INTERFACE
  INTERFACE cumsum
     MODULE PROCEDURE cumsum_r,cumsum_i
  END INTERFACE
  INTERFACE poly
     MODULE PROCEDURE poly_rr,poly_rrv,poly_dd,poly_ddv,&
          poly_rc,poly_cc,poly_msk_rrv,poly_msk_ddv
  END INTERFACE
  INTERFACE poly_term
     MODULE PROCEDURE poly_term_rr,poly_term_cc
  END INTERFACE
  INTERFACE outerprod
     MODULE PROCEDURE outerprod_r,outerprod_d
  END INTERFACE
  INTERFACE outerdiff
     MODULE PROCEDURE outerdiff_r,outerdiff_d,outerdiff_i
  END INTERFACE
  INTERFACE scatter_add
     MODULE PROCEDURE scatter_add_r,scatter_add_d
  END INTERFACE
  INTERFACE scatter_max
     MODULE PROCEDURE scatter_max_r,scatter_max_d
  END INTERFACE
  INTERFACE diagadd
     MODULE PROCEDURE diagadd_rv,diagadd_r
  END INTERFACE
  INTERFACE diagmult
     MODULE PROCEDURE diagmult_rv,diagmult_r
  END INTERFACE
  INTERFACE get_diag
     MODULE PROCEDURE get_diag_rv, get_diag_dv
  END INTERFACE
  INTERFACE put_diag
     MODULE PROCEDURE put_diag_rv, put_diag_r
  END INTERFACE
CONTAINS
  !BL
  SUBROUTINE array_copy_r(src,dest,n_copied,n_not_copied)
    REAL(SP), DIMENSION(:), INTENT(IN) :: src
    REAL(SP), DIMENSION(:), INTENT(OUT) :: dest
    INTEGER(I4B), INTENT(OUT) :: n_copied, n_not_copied
    n_copied=min(size(src),size(dest))
    n_not_copied=size(src)-n_copied
    dest(1:n_copied)=src(1:n_copied)
  END SUBROUTINE array_copy_r
  !BL
  SUBROUTINE array_copy_d(src,dest,n_copied,n_not_copied)
    REAL(DP), DIMENSION(:), INTENT(IN) :: src
    REAL(DP), DIMENSION(:), INTENT(OUT) :: dest
    INTEGER(I4B), INTENT(OUT) :: n_copied, n_not_copied
    n_copied=min(size(src),size(dest))
    n_not_copied=size(src)-n_copied
    dest(1:n_copied)=src(1:n_copied)
  END SUBROUTINE array_copy_d
  !BL
  SUBROUTINE array_copy_i(src,dest,n_copied,n_not_copied)
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: src
    INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: dest
    INTEGER(I4B), INTENT(OUT) :: n_copied, n_not_copied
    n_copied=min(size(src),size(dest))
    n_not_copied=size(src)-n_copied
    dest(1:n_copied)=src(1:n_copied)
  END SUBROUTINE array_copy_i
  !BL
  !BL
  SUBROUTINE swap_i(a,b)
    INTEGER(I4B), INTENT(INOUT) :: a,b
    INTEGER(I4B) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_i
  !BL
  SUBROUTINE swap_r(a,b)
    REAL(SP), INTENT(INOUT) :: a,b
    REAL(SP) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_r
  !BL
  SUBROUTINE swap_rv(a,b)
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
    REAL(SP), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_rv
  !BL
  SUBROUTINE swap_c(a,b)
    COMPLEX(SPC), INTENT(INOUT) :: a,b
    COMPLEX(SPC) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_c
  !BL
  SUBROUTINE swap_cv(a,b)
    COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: a,b
    COMPLEX(SPC), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_cv
  !BL
  SUBROUTINE swap_cm(a,b)
    COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
    COMPLEX(SPC), DIMENSION(size(a,1),size(a,2)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_cm
  !BL
  SUBROUTINE swap_z(a,b)
    COMPLEX(DPC), INTENT(INOUT) :: a,b
    COMPLEX(DPC) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_z
  !BL
  SUBROUTINE swap_zv(a,b)
    COMPLEX(DPC), DIMENSION(:), INTENT(INOUT) :: a,b
    COMPLEX(DPC), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_zv
  !BL
  SUBROUTINE swap_zm(a,b)
    COMPLEX(DPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
    COMPLEX(DPC), DIMENSION(size(a,1),size(a,2)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_zm
  !BL
  SUBROUTINE masked_swap_rs(a,b,mask)
    REAL(SP), INTENT(INOUT) :: a,b
    LOGICAL(LGT), INTENT(IN) :: mask
    REAL(SP) :: swp
    if (mask) then
       swp=a
       a=b
       b=swp
    end if
  END SUBROUTINE masked_swap_rs
  !BL
  SUBROUTINE masked_swap_rv(a,b,mask)
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
    LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
    REAL(SP), DIMENSION(size(a)) :: swp
    where (mask)
       swp=a
       a=b
       b=swp
    end where
  END SUBROUTINE masked_swap_rv
  !BL
  SUBROUTINE masked_swap_rm(a,b,mask)
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a,b
    LOGICAL(LGT), DIMENSION(:,:), INTENT(IN) :: mask
    REAL(SP), DIMENSION(size(a,1),size(a,2)) :: swp
    where (mask)
       swp=a
       a=b
       b=swp
    end where
  END SUBROUTINE masked_swap_rm
  !BL
  !BL
  FUNCTION reallocate_rv(p,n)
    REAL(SP), DIMENSION(:), POINTER :: p, reallocate_rv
    INTEGER(I4B), INTENT(IN) :: n
    INTEGER(I4B) :: nold,ierr
    allocate(reallocate_rv(n),stat=ierr)
    if (ierr /= 0) call &
         nrerror('reallocate_rv: problem in attempt to allocate memory')
    if (.not. associated(p)) RETURN
    nold=size(p)
    reallocate_rv(1:min(nold,n))=p(1:min(nold,n))
    deallocate(p)
  END FUNCTION reallocate_rv
  !BL
  FUNCTION reallocate_iv(p,n)
    INTEGER(I4B), DIMENSION(:), POINTER :: p, reallocate_iv
    INTEGER(I4B), INTENT(IN) :: n
    INTEGER(I4B) :: nold,ierr
    allocate(reallocate_iv(n),stat=ierr)
    if (ierr /= 0) call &
         nrerror('reallocate_iv: problem in attempt to allocate memory')
    if (.not. associated(p)) RETURN
    nold=size(p)
    reallocate_iv(1:min(nold,n))=p(1:min(nold,n))
    deallocate(p)
  END FUNCTION reallocate_iv
  !BL
  FUNCTION reallocate_hv(p,n)
    CHARACTER(1), DIMENSION(:), POINTER :: p, reallocate_hv
    INTEGER(I4B), INTENT(IN) :: n
    INTEGER(I4B) :: nold,ierr
    allocate(reallocate_hv(n),stat=ierr)
    if (ierr /= 0) call &
         nrerror('reallocate_hv: problem in attempt to allocate memory')
    if (.not. associated(p)) RETURN
    nold=size(p)
    reallocate_hv(1:min(nold,n))=p(1:min(nold,n))
    deallocate(p)
  END FUNCTION reallocate_hv
  !BL
  FUNCTION reallocate_rm(p,n,m)
    REAL(SP), DIMENSION(:,:), POINTER :: p, reallocate_rm
    INTEGER(I4B), INTENT(IN) :: n,m
    INTEGER(I4B) :: nold,mold,ierr
    allocate(reallocate_rm(n,m),stat=ierr)
    if (ierr /= 0) call &
         nrerror('reallocate_rm: problem in attempt to allocate memory')
    if (.not. associated(p)) RETURN
    nold=size(p,1)
    mold=size(p,2)
    reallocate_rm(1:min(nold,n),1:min(mold,m))=&
         p(1:min(nold,n),1:min(mold,m))
    deallocate(p)
  END FUNCTION reallocate_rm
  !BL
  FUNCTION reallocate_im(p,n,m)
    INTEGER(I4B), DIMENSION(:,:), POINTER :: p, reallocate_im
    INTEGER(I4B), INTENT(IN) :: n,m
    INTEGER(I4B) :: nold,mold,ierr
    allocate(reallocate_im(n,m),stat=ierr)
    if (ierr /= 0) call &
         nrerror('reallocate_im: problem in attempt to allocate memory')
    if (.not. associated(p)) RETURN
    nold=size(p,1)
    mold=size(p,2)
    reallocate_im(1:min(nold,n),1:min(mold,m))=&
         p(1:min(nold,n),1:min(mold,m))
    deallocate(p)
  END FUNCTION reallocate_im
  !BL
  FUNCTION ifirstloc(mask)
    LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
    INTEGER(I4B) :: ifirstloc
    INTEGER(I4B), DIMENSION(1) :: loc
    loc=maxloc(merge(1,0,mask))
    ifirstloc=loc(1)
    if (.not. mask(ifirstloc)) ifirstloc=size(mask)+1
  END FUNCTION ifirstloc
  !BL
  FUNCTION imaxloc_r(arr)
    REAL(SP), DIMENSION(:), INTENT(IN) :: arr
    INTEGER(I4B) :: imaxloc_r
    INTEGER(I4B), DIMENSION(1) :: imax
    imax=maxloc(arr(:))
    imaxloc_r=imax(1)
  END FUNCTION imaxloc_r
  !BL
  FUNCTION imaxloc_i(iarr)
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: iarr
    INTEGER(I4B), DIMENSION(1) :: imax
    INTEGER(I4B) :: imaxloc_i
    imax=maxloc(iarr(:))
    imaxloc_i=imax(1)
  END FUNCTION imaxloc_i
  !BL
  FUNCTION iminloc(arr)
    REAL(SP), DIMENSION(:), INTENT(IN) :: arr
    INTEGER(I4B), DIMENSION(1) :: imin
    INTEGER(I4B) :: iminloc
    imin=minloc(arr(:))
    iminloc=imin(1)
  END FUNCTION iminloc
  !BL
  SUBROUTINE assert1(n1,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL, INTENT(IN) :: n1
    if (.not. n1) then
       write (*,*) 'nrerror: an assertion failed with this tag:', &
            string
       STOP 'program terminated by assert1'
    end if
  END SUBROUTINE assert1
  !BL
  SUBROUTINE assert2(n1,n2,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL, INTENT(IN) :: n1,n2
    if (.not. (n1 .and. n2)) then
       write (*,*) 'nrerror: an assertion failed with this tag:', &
            string
       STOP 'program terminated by assert2'
    end if
  END SUBROUTINE assert2
  !BL
  SUBROUTINE assert3(n1,n2,n3,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL, INTENT(IN) :: n1,n2,n3
    if (.not. (n1 .and. n2 .and. n3)) then
       write (*,*) 'nrerror: an assertion failed with this tag:', &
            string
       STOP 'program terminated by assert3'
    end if
  END SUBROUTINE assert3
  !BL
  SUBROUTINE assert4(n1,n2,n3,n4,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL, INTENT(IN) :: n1,n2,n3,n4
    if (.not. (n1 .and. n2 .and. n3 .and. n4)) then
       write (*,*) 'nrerror: an assertion failed with this tag:', &
            string
       STOP 'program terminated by assert4'
    end if
  END SUBROUTINE assert4
  !BL
  SUBROUTINE assert_v(n,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL, DIMENSION(:), INTENT(IN) :: n
    if (.not. all(n)) then
       write (*,*) 'nrerror: an assertion failed with this tag:', &
            string
       STOP 'program terminated by assert_v'
    end if
  END SUBROUTINE assert_v
  !BL
  FUNCTION assert_eq2(n1,n2,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2
    INTEGER :: assert_eq2
    if (n1 == n2) then
       assert_eq2=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq2'
    end if
  END FUNCTION assert_eq2
  !BL
  FUNCTION assert_eq3(n1,n2,n3,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2,n3
    INTEGER :: assert_eq3
    if (n1 == n2 .and. n2 == n3) then
       assert_eq3=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq3'
    end if
  END FUNCTION assert_eq3
  !BL
  FUNCTION assert_eq4(n1,n2,n3,n4,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2,n3,n4
    INTEGER :: assert_eq4
    if (n1 == n2 .and. n2 == n3 .and. n3 == n4) then
       assert_eq4=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq4'
    end if
  END FUNCTION assert_eq4
  !BL
  FUNCTION assert_eqn(nn,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, DIMENSION(:), INTENT(IN) :: nn
    INTEGER :: assert_eqn
    if (all(nn(2:) == nn(1))) then
       assert_eqn=nn(1)
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eqn'
    end if
  END FUNCTION assert_eqn
  !BL
  SUBROUTINE nrerror(string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    write (*,*) 'nrerror: ',string
    STOP 'program terminated by nrerror'
  END SUBROUTINE nrerror
  !BL
  FUNCTION arth_r(first,increment,n)
    REAL(SP), INTENT(IN) :: first,increment
    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), DIMENSION(n) :: arth_r
    INTEGER(I4B) :: k,k2
    REAL(SP) :: temp
    if (n > 0) arth_r(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth_r(k)=arth_r(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth_r(k)=arth_r(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth_r(k+1:min(k2,n))=temp+arth_r(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
  END FUNCTION arth_r
  !BL
  FUNCTION arth_d(first,increment,n)
    REAL(DP), INTENT(IN) :: first,increment
    INTEGER(I4B), INTENT(IN) :: n
    REAL(DP), DIMENSION(n) :: arth_d
    INTEGER(I4B) :: k,k2
    REAL(DP) :: temp
    if (n > 0) arth_d(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth_d(k)=arth_d(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth_d(k)=arth_d(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth_d(k+1:min(k2,n))=temp+arth_d(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
  END FUNCTION arth_d
  !BL
  FUNCTION arth_i(first,increment,n)
    INTEGER(I4B), INTENT(IN) :: first,increment,n
    INTEGER(I4B), DIMENSION(n) :: arth_i
    INTEGER(I4B) :: k,k2,temp
    if (n > 0) arth_i(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth_i(k)=arth_i(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth_i(k)=arth_i(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
  END FUNCTION arth_i
  !BL
  !BL
  FUNCTION geop_r(first,factor,n)
    REAL(SP), INTENT(IN) :: first,factor
    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), DIMENSION(n) :: geop_r
    INTEGER(I4B) :: k,k2
    REAL(SP) :: temp
    if (n > 0) geop_r(1)=first
    if (n <= NPAR_GEOP) then
       do k=2,n
          geop_r(k)=geop_r(k-1)*factor
       end do
    else
       do k=2,NPAR2_GEOP
          geop_r(k)=geop_r(k-1)*factor
       end do
       temp=factor**NPAR2_GEOP
       k=NPAR2_GEOP
       do
          if (k >= n) exit
          k2=k+k
          geop_r(k+1:min(k2,n))=temp*geop_r(1:min(k,n-k))
          temp=temp*temp
          k=k2
       end do
    end if
  END FUNCTION geop_r
  !BL
  FUNCTION geop_d(first,factor,n)
    REAL(DP), INTENT(IN) :: first,factor
    INTEGER(I4B), INTENT(IN) :: n
    REAL(DP), DIMENSION(n) :: geop_d
    INTEGER(I4B) :: k,k2
    REAL(DP) :: temp
    if (n > 0) geop_d(1)=first
    if (n <= NPAR_GEOP) then
       do k=2,n
          geop_d(k)=geop_d(k-1)*factor
       end do
    else
       do k=2,NPAR2_GEOP
          geop_d(k)=geop_d(k-1)*factor
       end do
       temp=factor**NPAR2_GEOP
       k=NPAR2_GEOP
       do
          if (k >= n) exit
          k2=k+k
          geop_d(k+1:min(k2,n))=temp*geop_d(1:min(k,n-k))
          temp=temp*temp
          k=k2
       end do
    end if
  END FUNCTION geop_d
  !BL
  FUNCTION geop_i(first,factor,n)
    INTEGER(I4B), INTENT(IN) :: first,factor,n
    INTEGER(I4B), DIMENSION(n) :: geop_i
    INTEGER(I4B) :: k,k2,temp
    if (n > 0) geop_i(1)=first
    if (n <= NPAR_GEOP) then
       do k=2,n
          geop_i(k)=geop_i(k-1)*factor
       end do
    else
       do k=2,NPAR2_GEOP
          geop_i(k)=geop_i(k-1)*factor
       end do
       temp=factor**NPAR2_GEOP
       k=NPAR2_GEOP
       do
          if (k >= n) exit
          k2=k+k
          geop_i(k+1:min(k2,n))=temp*geop_i(1:min(k,n-k))
          temp=temp*temp
          k=k2
       end do
    end if
  END FUNCTION geop_i
  !BL
  FUNCTION geop_c(first,factor,n)
    COMPLEX(SP), INTENT(IN) :: first,factor
    INTEGER(I4B), INTENT(IN) :: n
    COMPLEX(SP), DIMENSION(n) :: geop_c
    INTEGER(I4B) :: k,k2
    COMPLEX(SP) :: temp
    if (n > 0) geop_c(1)=first
    if (n <= NPAR_GEOP) then
       do k=2,n
          geop_c(k)=geop_c(k-1)*factor
       end do
    else
       do k=2,NPAR2_GEOP
          geop_c(k)=geop_c(k-1)*factor
       end do
       temp=factor**NPAR2_GEOP
       k=NPAR2_GEOP
       do
          if (k >= n) exit
          k2=k+k
          geop_c(k+1:min(k2,n))=temp*geop_c(1:min(k,n-k))
          temp=temp*temp
          k=k2
       end do
    end if
  END FUNCTION geop_c
  !BL
  FUNCTION geop_dv(first,factor,n)
    REAL(DP), DIMENSION(:), INTENT(IN) :: first,factor
    INTEGER(I4B), INTENT(IN) :: n
    REAL(DP), DIMENSION(size(first),n) :: geop_dv
    INTEGER(I4B) :: k,k2
    REAL(DP), DIMENSION(size(first)) :: temp
    if (n > 0) geop_dv(:,1)=first(:)
    if (n <= NPAR_GEOP) then
       do k=2,n
          geop_dv(:,k)=geop_dv(:,k-1)*factor(:)
       end do
    else
       do k=2,NPAR2_GEOP
          geop_dv(:,k)=geop_dv(:,k-1)*factor(:)
       end do
       temp=factor**NPAR2_GEOP
       k=NPAR2_GEOP
       do
          if (k >= n) exit
          k2=k+k
          geop_dv(:,k+1:min(k2,n))=geop_dv(:,1:min(k,n-k))*&
               spread(temp,2,size(geop_dv(:,1:min(k,n-k)),2))
          temp=temp*temp
          k=k2
       end do
    end if
  END FUNCTION geop_dv
  !BL
  !BL
  RECURSIVE FUNCTION cumsum_r(arr,seed) RESULT(ans)
    REAL(SP), DIMENSION(:), INTENT(IN) :: arr
    REAL(SP), OPTIONAL, INTENT(IN) :: seed
    REAL(SP), DIMENSION(size(arr)) :: ans
    INTEGER(I4B) :: n,j
    REAL(SP) :: sd
    n=size(arr)
    if (n == 0_i4b) RETURN
    sd=0.0_sp
    if (present(seed)) sd=seed
    ans(1)=arr(1)+sd
    if (n < NPAR_CUMSUM) then
       do j=2,n
          ans(j)=ans(j-1)+arr(j)
       end do
    else
       ans(2:n:2)=cumsum_r(arr(2:n:2)+arr(1:n-1:2),sd)
       ans(3:n:2)=ans(2:n-1:2)+arr(3:n:2)
    end if
  END FUNCTION cumsum_r
  !BL
  RECURSIVE FUNCTION cumsum_i(arr,seed) RESULT(ans)
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: arr
    INTEGER(I4B), OPTIONAL, INTENT(IN) :: seed
    INTEGER(I4B), DIMENSION(size(arr)) :: ans
    INTEGER(I4B) :: n,j,sd
    n=size(arr)
    if (n == 0_i4b) RETURN
    sd=0_i4b
    if (present(seed)) sd=seed
    ans(1)=arr(1)+sd
    if (n < NPAR_CUMSUM) then
       do j=2,n
          ans(j)=ans(j-1)+arr(j)
       end do
    else
       ans(2:n:2)=cumsum_i(arr(2:n:2)+arr(1:n-1:2),sd)
       ans(3:n:2)=ans(2:n-1:2)+arr(3:n:2)
    end if
  END FUNCTION cumsum_i
  !BL
  !BL
  RECURSIVE FUNCTION cumprod(arr,seed) RESULT(ans)
    REAL(SP), DIMENSION(:), INTENT(IN) :: arr
    REAL(SP), OPTIONAL, INTENT(IN) :: seed
    REAL(SP), DIMENSION(size(arr)) :: ans
    INTEGER(I4B) :: n,j
    REAL(SP) :: sd
    n=size(arr)
    if (n == 0_i4b) RETURN
    sd=1.0_sp
    if (present(seed)) sd=seed
    ans(1)=arr(1)*sd
    if (n < NPAR_CUMPROD) then
       do j=2,n
          ans(j)=ans(j-1)*arr(j)
       end do
    else
       ans(2:n:2)=cumprod(arr(2:n:2)*arr(1:n-1:2),sd)
       ans(3:n:2)=ans(2:n-1:2)*arr(3:n:2)
    end if
  END FUNCTION cumprod
  !BL
  !BL
  FUNCTION poly_rr(x,coeffs)
    REAL(SP), INTENT(IN) :: x
    REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs
    REAL(SP) :: poly_rr
    REAL(SP) :: pow
    REAL(SP), DIMENSION(:), ALLOCATABLE :: vec
    INTEGER(I4B) :: i,n,nn
    n=size(coeffs)
    if (n <= 0) then
       poly_rr=0.0_sp
    else if (n < NPAR_POLY) then
       poly_rr=coeffs(n)
       do i=n-1,1,-1
          poly_rr=x*poly_rr+coeffs(i)
       end do
    else
       allocate(vec(n+1))
       pow=x
       vec(1:n)=coeffs
       do
          vec(n+1)=0.0_sp
          nn=ishft(n+1,-1)
          vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
          if (nn == 1) exit
          pow=pow*pow
          n=nn
       end do
       poly_rr=vec(1)
       deallocate(vec)
    end if
  END FUNCTION poly_rr
  !BL
  FUNCTION poly_dd(x,coeffs)
    REAL(DP), INTENT(IN) :: x
    REAL(DP), DIMENSION(:), INTENT(IN) :: coeffs
    REAL(DP) :: poly_dd
    REAL(DP) :: pow
    REAL(DP), DIMENSION(:), ALLOCATABLE :: vec
    INTEGER(I4B) :: i,n,nn
    n=size(coeffs)
    if (n <= 0) then
       poly_dd=0.0_dp
    else if (n < NPAR_POLY) then
       poly_dd=coeffs(n)
       do i=n-1,1,-1
          poly_dd=x*poly_dd+coeffs(i)
       end do
    else
       allocate(vec(n+1))
       pow=x
       vec(1:n)=coeffs
       do
          vec(n+1)=0.0_dp
          nn=ishft(n+1,-1)
          vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
          if (nn == 1) exit
          pow=pow*pow
          n=nn
       end do
       poly_dd=vec(1)
       deallocate(vec)
    end if
  END FUNCTION poly_dd
  !BL
  FUNCTION poly_rc(x,coeffs)
    COMPLEX(SPC), INTENT(IN) :: x
    REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs
    COMPLEX(SPC) :: poly_rc
    COMPLEX(SPC) :: pow
    COMPLEX(SPC), DIMENSION(:), ALLOCATABLE :: vec
    INTEGER(I4B) :: i,n,nn
    n=size(coeffs)
    if (n <= 0) then
       poly_rc=0.0_sp
    else if (n < NPAR_POLY) then
       poly_rc=coeffs(n)
       do i=n-1,1,-1
          poly_rc=x*poly_rc+coeffs(i)
       end do
    else
       allocate(vec(n+1))
       pow=x
       vec(1:n)=coeffs
       do
          vec(n+1)=0.0_sp
          nn=ishft(n+1,-1)
          vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
          if (nn == 1) exit
          pow=pow*pow
          n=nn
       end do
       poly_rc=vec(1)
       deallocate(vec)
    end if
  END FUNCTION poly_rc
  !BL
  FUNCTION poly_cc(x,coeffs)
    COMPLEX(SPC), INTENT(IN) :: x
    COMPLEX(SPC), DIMENSION(:), INTENT(IN) :: coeffs
    COMPLEX(SPC) :: poly_cc
    COMPLEX(SPC) :: pow
    COMPLEX(SPC), DIMENSION(:), ALLOCATABLE :: vec
    INTEGER(I4B) :: i,n,nn
    n=size(coeffs)
    if (n <= 0) then
       poly_cc=0.0_sp
    else if (n < NPAR_POLY) then
       poly_cc=coeffs(n)
       do i=n-1,1,-1
          poly_cc=x*poly_cc+coeffs(i)
       end do
    else
       allocate(vec(n+1))
       pow=x
       vec(1:n)=coeffs
       do
          vec(n+1)=0.0_sp
          nn=ishft(n+1,-1)
          vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
          if (nn == 1) exit
          pow=pow*pow
          n=nn
       end do
       poly_cc=vec(1)
       deallocate(vec)
    end if
  END FUNCTION poly_cc
  !BL
  FUNCTION poly_rrv(x,coeffs)
    REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs,x
    REAL(SP), DIMENSION(size(x)) :: poly_rrv
    INTEGER(I4B) :: i,n,m
    m=size(coeffs)
    n=size(x)
    if (m <= 0) then
       poly_rrv=0.0_sp
    else if (m < n .or. m < NPAR_POLY) then
       poly_rrv=coeffs(m)
       do i=m-1,1,-1
          poly_rrv=x*poly_rrv+coeffs(i)
       end do
    else
       do i=1,n
          poly_rrv(i)=poly_rr(x(i),coeffs)
       end do
    end if
  END FUNCTION poly_rrv
  !BL
  FUNCTION poly_ddv(x,coeffs)
    REAL(DP), DIMENSION(:), INTENT(IN) :: coeffs,x
    REAL(DP), DIMENSION(size(x)) :: poly_ddv
    INTEGER(I4B) :: i,n,m
    m=size(coeffs)
    n=size(x)
    if (m <= 0) then
       poly_ddv=0.0_dp
    else if (m < n .or. m < NPAR_POLY) then
       poly_ddv=coeffs(m)
       do i=m-1,1,-1
          poly_ddv=x*poly_ddv+coeffs(i)
       end do
    else
       do i=1,n
          poly_ddv(i)=poly_dd(x(i),coeffs)
       end do
    end if
  END FUNCTION poly_ddv
  !BL
  FUNCTION poly_msk_rrv(x,coeffs,mask)
    REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs,x
    LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
    REAL(SP), DIMENSION(size(x)) :: poly_msk_rrv
    poly_msk_rrv=unpack(poly_rrv(pack(x,mask),coeffs),mask,0.0_sp)
  END FUNCTION poly_msk_rrv
  !BL
  FUNCTION poly_msk_ddv(x,coeffs,mask)
    REAL(DP), DIMENSION(:), INTENT(IN) :: coeffs,x
    LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
    REAL(DP), DIMENSION(size(x)) :: poly_msk_ddv
    poly_msk_ddv=unpack(poly_ddv(pack(x,mask),coeffs),mask,0.0_dp)
  END FUNCTION poly_msk_ddv
  !BL
  !BL
  RECURSIVE FUNCTION poly_term_rr(a,b) RESULT(u)
    REAL(SP), DIMENSION(:), INTENT(IN) :: a
    REAL(SP), INTENT(IN) :: b
    REAL(SP), DIMENSION(size(a)) :: u
    INTEGER(I4B) :: n,j
    n=size(a)
    if (n <= 0) RETURN
    u(1)=a(1)
    if (n < NPAR_POLYTERM) then
       do j=2,n
          u(j)=a(j)+b*u(j-1)
       end do
    else
       u(2:n:2)=poly_term_rr(a(2:n:2)+a(1:n-1:2)*b,b*b)
       u(3:n:2)=a(3:n:2)+b*u(2:n-1:2)
    end if
  END FUNCTION poly_term_rr
  !BL
  RECURSIVE FUNCTION poly_term_cc(a,b) RESULT(u)
    COMPLEX(SPC), DIMENSION(:), INTENT(IN) :: a
    COMPLEX(SPC), INTENT(IN) :: b
    COMPLEX(SPC), DIMENSION(size(a)) :: u
    INTEGER(I4B) :: n,j
    n=size(a)
    if (n <= 0) RETURN
    u(1)=a(1)
    if (n < NPAR_POLYTERM) then
       do j=2,n
          u(j)=a(j)+b*u(j-1)
       end do
    else
       u(2:n:2)=poly_term_cc(a(2:n:2)+a(1:n-1:2)*b,b*b)
       u(3:n:2)=a(3:n:2)+b*u(2:n-1:2)
    end if
  END FUNCTION poly_term_cc
  !BL
  !BL
  FUNCTION zroots_unity(n,nn)
    INTEGER(I4B), INTENT(IN) :: n,nn
    COMPLEX(SPC), DIMENSION(nn) :: zroots_unity
    INTEGER(I4B) :: k
    REAL(SP) :: theta
    zroots_unity(1)=1.0
    theta=TWOPI/n
    k=1
    do
       if (k >= nn) exit
       zroots_unity(k+1)=cmplx(cos(k*theta),sin(k*theta),SPC)
       zroots_unity(k+2:min(2*k,nn))=zroots_unity(k+1)*&
            zroots_unity(2:min(k,nn-k))
       k=2*k
    end do
  END FUNCTION zroots_unity
  !BL
  FUNCTION outerprod_r(a,b)
    REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
    REAL(SP), DIMENSION(size(a),size(b)) :: outerprod_r
    outerprod_r = spread(a,dim=2,ncopies=size(b)) * &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerprod_r
  !BL
  FUNCTION outerprod_d(a,b)
    REAL(DP), DIMENSION(:), INTENT(IN) :: a,b
    REAL(DP), DIMENSION(size(a),size(b)) :: outerprod_d
    outerprod_d = spread(a,dim=2,ncopies=size(b)) * &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerprod_d
  !BL
  FUNCTION outerdiv(a,b)
    REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
    REAL(SP), DIMENSION(size(a),size(b)) :: outerdiv
    outerdiv = spread(a,dim=2,ncopies=size(b)) / &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerdiv
  !BL
  FUNCTION outersum(a,b)
    REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
    REAL(SP), DIMENSION(size(a),size(b)) :: outersum
    outersum = spread(a,dim=2,ncopies=size(b)) + &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outersum
  !BL
  FUNCTION outerdiff_r(a,b)
    REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
    REAL(SP), DIMENSION(size(a),size(b)) :: outerdiff_r
    outerdiff_r = spread(a,dim=2,ncopies=size(b)) - &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerdiff_r
  !BL
  FUNCTION outerdiff_d(a,b)
    REAL(DP), DIMENSION(:), INTENT(IN) :: a,b
    REAL(DP), DIMENSION(size(a),size(b)) :: outerdiff_d
    outerdiff_d = spread(a,dim=2,ncopies=size(b)) - &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerdiff_d
  !BL
  FUNCTION outerdiff_i(a,b)
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: a,b
    INTEGER(I4B), DIMENSION(size(a),size(b)) :: outerdiff_i
    outerdiff_i = spread(a,dim=2,ncopies=size(b)) - &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerdiff_i
  !BL
  FUNCTION outerand(a,b)
    LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: a,b
    LOGICAL(LGT), DIMENSION(size(a),size(b)) :: outerand
    outerand = spread(a,dim=2,ncopies=size(b)) .and. &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerand
  !BL
  SUBROUTINE scatter_add_r(dest,source,dest_index)
    REAL(SP), DIMENSION(:), INTENT(OUT) :: dest
    REAL(SP), DIMENSION(:), INTENT(IN) :: source
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: dest_index
    INTEGER(I4B) :: m,n,j,i
    n=assert_eq2(size(source),size(dest_index),'scatter_add_r')
    m=size(dest)
    do j=1,n
       i=dest_index(j)
       if (i > 0 .and. i <= m) dest(i)=dest(i)+source(j)
    end do
  END SUBROUTINE scatter_add_r
  SUBROUTINE scatter_add_d(dest,source,dest_index)
    REAL(DP), DIMENSION(:), INTENT(OUT) :: dest
    REAL(DP), DIMENSION(:), INTENT(IN) :: source
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: dest_index
    INTEGER(I4B) :: m,n,j,i
    n=assert_eq2(size(source),size(dest_index),'scatter_add_d')
    m=size(dest)
    do j=1,n
       i=dest_index(j)
       if (i > 0 .and. i <= m) dest(i)=dest(i)+source(j)
    end do
  END SUBROUTINE scatter_add_d
  SUBROUTINE scatter_max_r(dest,source,dest_index)
    REAL(SP), DIMENSION(:), INTENT(OUT) :: dest
    REAL(SP), DIMENSION(:), INTENT(IN) :: source
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: dest_index
    INTEGER(I4B) :: m,n,j,i
    n=assert_eq2(size(source),size(dest_index),'scatter_max_r')
    m=size(dest)
    do j=1,n
       i=dest_index(j)
       if (i > 0 .and. i <= m) dest(i)=max(dest(i),source(j))
    end do
  END SUBROUTINE scatter_max_r
  SUBROUTINE scatter_max_d(dest,source,dest_index)
    REAL(DP), DIMENSION(:), INTENT(OUT) :: dest
    REAL(DP), DIMENSION(:), INTENT(IN) :: source
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: dest_index
    INTEGER(I4B) :: m,n,j,i
    n=assert_eq2(size(source),size(dest_index),'scatter_max_d')
    m=size(dest)
    do j=1,n
       i=dest_index(j)
       if (i > 0 .and. i <= m) dest(i)=max(dest(i),source(j))
    end do
  END SUBROUTINE scatter_max_d
  !BL
  SUBROUTINE diagadd_rv(mat,diag)
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
    REAL(SP), DIMENSION(:), INTENT(IN) :: diag
    INTEGER(I4B) :: j,n
    n = assert_eq2(size(diag),min(size(mat,1),size(mat,2)),'diagadd_rv')
    do j=1,n
       mat(j,j)=mat(j,j)+diag(j)
    end do
  END SUBROUTINE diagadd_rv
  !BL
  SUBROUTINE diagadd_r(mat,diag)
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
    REAL(SP), INTENT(IN) :: diag
    INTEGER(I4B) :: j,n
    n = min(size(mat,1),size(mat,2))
    do j=1,n
       mat(j,j)=mat(j,j)+diag
    end do
  END SUBROUTINE diagadd_r
  !BL
  SUBROUTINE diagmult_rv(mat,diag)
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
    REAL(SP), DIMENSION(:), INTENT(IN) :: diag
    INTEGER(I4B) :: j,n
    n = assert_eq2(size(diag),min(size(mat,1),size(mat,2)),'diagmult_rv')
    do j=1,n
       mat(j,j)=mat(j,j)*diag(j)
    end do
  END SUBROUTINE diagmult_rv
  !BL
  SUBROUTINE diagmult_r(mat,diag)
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
    REAL(SP), INTENT(IN) :: diag
    INTEGER(I4B) :: j,n
    n = min(size(mat,1),size(mat,2))
    do j=1,n
       mat(j,j)=mat(j,j)*diag
    end do
  END SUBROUTINE diagmult_r
  !BL
  FUNCTION get_diag_rv(mat)
    REAL(SP), DIMENSION(:,:), INTENT(IN) :: mat
    REAL(SP), DIMENSION(size(mat,1)) :: get_diag_rv
    INTEGER(I4B) :: j
    j=assert_eq2(size(mat,1),size(mat,2),'get_diag_rv')
    do j=1,size(mat,1)
       get_diag_rv(j)=mat(j,j)
    end do
  END FUNCTION get_diag_rv
  !BL
  FUNCTION get_diag_dv(mat)
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: mat
    REAL(DP), DIMENSION(size(mat,1)) :: get_diag_dv
    INTEGER(I4B) :: j
    j=assert_eq2(size(mat,1),size(mat,2),'get_diag_dv')
    do j=1,size(mat,1)
       get_diag_dv(j)=mat(j,j)
    end do
  END FUNCTION get_diag_dv
  !BL
  SUBROUTINE put_diag_rv(diagv,mat)
    REAL(SP), DIMENSION(:), INTENT(IN) :: diagv
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
    INTEGER(I4B) :: j,n
    n=assert_eq2(size(diagv),min(size(mat,1),size(mat,2)),'put_diag_rv')
    do j=1,n
       mat(j,j)=diagv(j)
    end do
  END SUBROUTINE put_diag_rv
  !BL
  SUBROUTINE put_diag_r(scal,mat)
    REAL(SP), INTENT(IN) :: scal
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
    INTEGER(I4B) :: j,n
    n = min(size(mat,1),size(mat,2))
    do j=1,n
       mat(j,j)=scal
    end do
  END SUBROUTINE put_diag_r
  !BL
  SUBROUTINE unit_matrix(mat)
    REAL(SP), DIMENSION(:,:), INTENT(OUT) :: mat
    INTEGER(I4B) :: i,n
    n=min(size(mat,1),size(mat,2))
    mat(:,:)=0.0_sp
    do i=1,n
       mat(i,i)=1.0_sp
    end do
  END SUBROUTINE unit_matrix
  !BL
  FUNCTION upper_triangle(j,k,extra)
    INTEGER(I4B), INTENT(IN) :: j,k
    INTEGER(I4B), OPTIONAL, INTENT(IN) :: extra
    LOGICAL(LGT), DIMENSION(j,k) :: upper_triangle
    INTEGER(I4B) :: n
    n=0
    if (present(extra)) n=extra
    upper_triangle=(outerdiff(arth_i(1,1,j),arth_i(1,1,k)) < n)
  END FUNCTION upper_triangle
  !BL
  FUNCTION lower_triangle(j,k,extra)
    INTEGER(I4B), INTENT(IN) :: j,k
    INTEGER(I4B), OPTIONAL, INTENT(IN) :: extra
    LOGICAL(LGT), DIMENSION(j,k) :: lower_triangle
    INTEGER(I4B) :: n
    n=0
    if (present(extra)) n=extra
    lower_triangle=(outerdiff(arth_i(1,1,j),arth_i(1,1,k)) > -n)
  END FUNCTION lower_triangle
  !BL
  FUNCTION vabs(v)
    REAL(SP), DIMENSION(:), INTENT(IN) :: v
    REAL(SP) :: vabs
    vabs=sqrt(dot_product(v,v))
  END FUNCTION vabs
  !BL
END MODULE nrutil
MODULE nr
  USE nrtype

  INTERFACE tridag
     module procedure tridag_par
  END INTERFACE

contains
  FUNCTION locate(xx,x)
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: xx
    REAL(DP), INTENT(IN) :: x
    INTEGER(I4B) :: locate
    INTEGER(I4B) :: n,jl,jm,ju
    LOGICAL :: ascnd
    n=size(xx)
    ascnd = (xx(n) >= xx(1))
    jl=0
    ju=n+1
    do
       if (ju-jl <= 1) exit
       jm=(ju+jl)/2
       if (ascnd .eqv. (x >= xx(jm))) then
          jl=jm
       else
          ju=jm
       end if
    end do
    if (x == xx(1)) then
       locate=1
    else if (x == xx(n)) then
       locate=n-1
    else
       locate=jl
    end if
  END FUNCTION locate

  FUNCTION splin2(x1a,x2a,ya,y2a,x1,x2)
    USE nrutil, ONLY : assert_eq
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: x1a,x2a
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: ya,y2a
    REAL(DP), INTENT(IN) :: x1,x2
    REAL(DP) :: splin2
    INTEGER(I4B) :: j,m,ndum
    REAL(DP), DIMENSION(size(x1a)) :: yytmp,y2tmp2
    m=assert_eq(size(x1a),size(ya,1),size(y2a,1),'splin2: m')
    ndum=assert_eq(size(x2a),size(ya,2),size(y2a,2),'splin2: ndum')
    do j=1,m
       yytmp(j)=splint(x2a,ya(j,:),y2a(j,:),x2)
    end do
    call spline(x1a,yytmp,1.0e30_dp,1.0e30_dp,y2tmp2)
    splin2=splint(x1a,yytmp,y2tmp2,x1)
  END FUNCTION splin2
  SUBROUTINE splie2(x1a,x2a,ya,y2a)
    USE nrutil, ONLY : assert_eq
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: x1a,x2a
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: ya
    REAL(DP), DIMENSION(:,:), INTENT(OUT) :: y2a
    INTEGER(I4B) :: j,m,ndum
    m=assert_eq(size(x1a),size(ya,1),size(y2a,1),'splie2: m')
    ndum=assert_eq(size(x2a),size(ya,2),size(y2a,2),'splie2: ndum')
    do j=1,m
       call spline(x2a,ya(j,:),1.0e30_dp,1.0e30_dp,y2a(j,:))
    end do
  END SUBROUTINE splie2
  SUBROUTINE tridag_ser(a,b,c,r,u)
    USE nrutil, ONLY : assert_eq,nrerror
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: a,b,c,r
    REAL(DP), DIMENSION(:), INTENT(OUT) :: u
    REAL(DP), DIMENSION(size(b)) :: gam
    INTEGER(I4B) :: n,j
    REAL(DP) :: bet
    n=assert_eq((/size(a)+1,size(b),size(c)+1,size(r),size(u)/),'tridag_ser')
    bet=b(1)
    if (bet == 0.0) call nrerror('tridag_ser: Error at code stage 1')
    u(1)=r(1)/bet
    do j=2,n
       gam(j)=c(j-1)/bet
       bet=b(j)-a(j-1)*gam(j)
       if (bet == 0.0) &
            call nrerror('tridag_ser: Error at code stage 2')
       u(j)=(r(j)-a(j-1)*u(j-1))/bet
    end do
    do j=n-1,1,-1
       u(j)=u(j)-gam(j+1)*u(j+1)
    end do
  END SUBROUTINE tridag_ser

  RECURSIVE SUBROUTINE tridag_par(a,b,c,r,u)
    USE nrutil, ONLY : assert_eq,nrerror
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: a,b,c,r
    REAL(DP), DIMENSION(:), INTENT(OUT) :: u
    INTEGER(I4B), PARAMETER :: NPAR_TRIDAG=4
    INTEGER(I4B) :: n,n2,nm,nx
    REAL(DP), DIMENSION(size(b)/2) :: y,q,piva
    REAL(DP), DIMENSION(size(b)/2-1) :: x,z
    REAL(DP), DIMENSION(size(a)/2) :: pivc
    n=assert_eq((/size(a)+1,size(b),size(c)+1,size(r),size(u)/),'tridag_par')
    if (n < NPAR_TRIDAG) then
       call tridag_ser(a,b,c,r,u)
    else
       if (maxval(abs(b(1:n))) == 0.0) &
            call nrerror('tridag_par: possible singular matrix')
       n2=size(y)
       nm=size(pivc)
       nx=size(x)
       piva = a(1:n-1:2)/b(1:n-1:2)
       pivc = c(2:n-1:2)/b(3:n:2)
       y(1:nm) = b(2:n-1:2)-piva(1:nm)*c(1:n-2:2)-pivc*a(2:n-1:2)
       q(1:nm) = r(2:n-1:2)-piva(1:nm)*r(1:n-2:2)-pivc*r(3:n:2)
       if (nm < n2) then
          y(n2) = b(n)-piva(n2)*c(n-1)
          q(n2) = r(n)-piva(n2)*r(n-1)
       end if
       x = -piva(2:n2)*a(2:n-2:2)
       z = -pivc(1:nx)*c(3:n-1:2)
       call tridag_par(x,y,z,q,u(2:n:2))
       u(1) = (r(1)-c(1)*u(2))/b(1)
       u(3:n-1:2) = (r(3:n-1:2)-a(2:n-2:2)*u(2:n-2:2) &
            -c(3:n-1:2)*u(4:n:2))/b(3:n-1:2)
       if (nm == n2) u(n)=(r(n)-a(n-1)*u(n-1))/b(n)
    end if
  END SUBROUTINE tridag_par
  SUBROUTINE spline(x,y,yp1,ypn,y2)
    USE nrutil, ONLY : assert_eq
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: x,y
    REAL(DP), INTENT(IN) :: yp1,ypn
    REAL(DP), DIMENSION(:), INTENT(OUT) :: y2
    INTEGER(I4B) :: n
    REAL(DP), DIMENSION(size(x)) :: a,b,c,r
    n=assert_eq(size(x),size(y),size(y2),'spline')
    c(1:n-1)=x(2:n)-x(1:n-1)
    r(1:n-1)=6.0_dp*((y(2:n)-y(1:n-1))/c(1:n-1))
    r(2:n-1)=r(2:n-1)-r(1:n-2)
    a(2:n-1)=c(1:n-2)
    b(2:n-1)=2.0_dp*(c(2:n-1)+a(2:n-1))
    b(1)=1.0
    b(n)=1.0
    if (yp1 > 0.99e30_dp) then
       r(1)=0.0
       c(1)=0.0
    else
       r(1)=(3.0_dp/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
       c(1)=0.5
    end if
    if (ypn > 0.99e30_dp) then
       r(n)=0.0
       a(n)=0.0
    else
       r(n)=(-3.0_dp/(x(n)-x(n-1)))*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)
       a(n)=0.5
    end if
    call tridag(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n))
  END SUBROUTINE spline
  FUNCTION splint(xa,ya,y2a,x)
    USE nrutil, ONLY : assert_eq,nrerror
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: xa,ya,y2a
    REAL(DP), INTENT(IN) :: x
    REAL(DP) :: splint
    INTEGER(I4B) :: khi,klo,n
    REAL(DP) :: a,b,h
    n=assert_eq(size(xa),size(ya),size(y2a),'splint')
    klo=max(min(locate(xa,x),n-1),1)
    khi=klo+1
    h=xa(khi)-xa(klo)
    if (h == 0.0) call nrerror('bad xa input in splint')
    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
    splint=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0_dp
  END FUNCTION splint
  SUBROUTINE gauleg(x1,x2,x,w)
    USE nrutil, ONLY : arth,assert_eq,nrerror
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: x1,x2
    REAL(DP), DIMENSION(:), INTENT(OUT) :: x,w
    REAL(DP), PARAMETER :: EPS=3.0e-14_dp
    INTEGER(I4B) :: its,j,m,n
    INTEGER(I4B), PARAMETER :: MAXIT=10
    REAL(DP) :: xl,xm
    REAL(DP), DIMENSION((size(x)+1)/2) :: p1,p2,p3,pp,z,z1
    LOGICAL(LGT), DIMENSION((size(x)+1)/2) :: unfinished
    n=assert_eq(size(x),size(w),'gauleg')
    m=(n+1)/2
    xm=0.5_dp*(x2+x1)
    xl=0.5_dp*(x2-x1)
    z=cos(PI_D*(arth(1,1,m)-0.25_dp)/(n+0.5_dp))
    unfinished=.true.
    do its=1,MAXIT
       where (unfinished)
          p1=1.0
          p2=0.0
       end where
       do j=1,n
          where (unfinished)
             p3=p2
             p2=p1
             p1=((2.0_dp*j-1.0_dp)*z*p2-(j-1.0_dp)*p3)/j
          end where
       end do
       where (unfinished)
          pp=n*(z*p1-p2)/(z*z-1.0_dp)
          z1=z
          z=z1-p1/pp
          unfinished=(abs(z-z1) > EPS)
       end where
       if (.not. any(unfinished)) exit
    end do
    if (its == MAXIT+1) call nrerror('too many iterations in gauleg')
    x(1:m)=xm-xl*z
    x(n:n-m+1:-1)=xm+xl*z
    w(1:m)=2.0_dp*xl/((1.0_dp-z**2)*pp**2)
    w(n:n-m+1:-1)=w(1:m)
  END SUBROUTINE gauleg
END MODULE nr
module m_error_func
  use nrtype
  real(kind=DP), parameter :: o_sqpi = -1.12837916709551261672d0
				  ! -2.d0/sqrt(pi)

  interface erfc
     module procedure erfc_d
  end interface

  interface erfc_f
     module procedure erfc_fd
  end interface

  interface d_erfc
     module procedure d_erfc_d
  end interface

  interface erf
     module procedure erf_f
  end interface

contains
  function erfc_d(z)
    implicit none
    real(kind=DP), intent(in) :: z
    real(kind=DP) :: erfc_d
    erfc_d = 1d0/(1.d0 + &
         0.0705230784d0*z   + 0.0422820123d0*z**2 + &
         0.0092705272d0*z**3+0.0001520143d0*z**4 + & 
         0.0002765672d0*z**5+0.0000430638d0*z**6)**16
  end function erfc_d
  function erfc_fd(z)
    implicit none
    real(kind=DP), intent(in) :: z
    real(kind=DP) :: erfc_fd, z2, z3, z4, z5, z6, a, b, c, d, e, f, s, t, &
         zz, sgn,sft
    zz = z
    sgn = 1.d0
    sft = 0.d0
    if(z < 0.d0) then
       zz = -z
       sgn = -1.d0
       sft = 2.d0
    end if
    z2 = zz*zz
    z3 = z2*zz
    z4 = z3*zz
    z5 = z4*zz
    z6 = z5*zz
    a = 0.0705230784d0*zz
    s = 1.d0 + a
    b = 0.0422820123d0*z2
    s = s + b
    c = 0.0092705272d0*z3
    s = s + c
    d = 0.0001520143d0*z4
    s = s + d
    e = 0.0002765672d0*z5
    s = s + e
    f = 0.0000430638d0*z6
    s = s + f
    t = s**16
    erfc_fd = sgn/t+sft
  end function erfc_fd
  function d_erfc_d(z)
    real(kind=DP), intent(in) :: z
    real(kind=DP) :: d_erfc_d
    d_erfc_d = exp(-z**2) * o_sqpi
  end function d_erfc_d
  function erf_f(z)
    implicit none
    real(kind=DP), intent(in) :: z
    real(kind=DP) :: erf_f, z2, z3, z4, z5, z6, a, b, c, d, e, f, s, t, &
         zz
    zz = z
    z2 = zz*zz
    z3 = z2*zz
    z4 = z3*zz
    z5 = z4*zz
    z6 = z5*zz
    a = 0.0705230784d0*zz
    s = 1.d0 + a
    b = 0.0422820123d0*z2
    s = s + b
    c = 0.0092705272d0*z3
    s = s + c
    d = 0.0001520143d0*z4
    s = s + d
    e = 0.0002765672d0*z5
    s = s + e
    f = 0.0000430638d0*z6
    s = s + f
    t = s**16
    erf_f = 1.d0/t
  end function erf_f
end module m_error_func
module m_vld
  use nrtype
  use m_error_func

  type coefpot
     real(kind=DP) :: a1, a2, c1
  end type coefpot

contains  
 function vl1(g,cc)
   implicit none
   real (kind=DP) :: vl1
   intent(in) :: g, cc
   real (kind=DP) :: g, g2, gg, x1, x2, exp1, exp2
   type (coefpot) :: cc

   g2 = -g**2
   gg = g2*0.25d0

   x1 = gg/cc%a1
   x2 = gg/cc%a2
   exp1 = exp(x1)
   exp2 = exp(x2)

   vl1 = (cc%c1*(exp1-exp2) + exp2)/g2

 end function vl1
 function vlr(r,cc,zval)
   implicit none
   real (kind=DP) :: vlr
   intent(in) :: r, cc, zval
   real (kind=DP) :: r, g2, rr, x1, x2, erf1, erf2, zval, rtmp
   type (coefpot) :: cc


!   print *,r,erf1,erf2
!   IF (dabs(R) < 1.0D-15) THEN 
!      vlr = -ZVAL*(CC%C1*(SQRT(CC%A1)-SQRT(CC%A2))+SQRT(CC%A2))
   IF (dabs(R) < 1.0D-10) THEN 
      rtmp = 1d-10
      x1 = rtmp*sqrt(cc%a1)
      x2 = rtmp*sqrt(cc%a2)

      erf1 = 1.d0-erfc(x1)
      erf2 = 1.d0-erfc(x2)
      vlr = -zval*(cc%c1*(erf1-erf2) + erf2)/rtmp
   else
      x1 = r*sqrt(cc%a1)
      x2 = r*sqrt(cc%a2)

      erf1 = 1.d0-erfc(x1)
      erf2 = 1.d0-erfc(x2)
      vlr = -zval*(cc%c1*(erf1-erf2) + erf2)/r
   end IF

 end function vlr
 SUBROUTINE VLD(atom_kind,CC)
   implicit none
   CHARACTER (len=*), intent(in) :: atom_kind
   type (coefpot), intent(out) :: CC
   character(len=10) :: atom
   integer(kind=I4B) :: iln
   ! -----< data for PP local part >-----
   atom = atom_kind
   iln = len(atom)

   select case(atom(1:iln))
   case('H')
      CC = coefpot(16.2200d0, 5.5500d0, 1.1924d0)
   case('Mu') ! must be same with H
      CC = coefpot(16.2200d0, 5.5500d0, 1.1924d0)
   case('D') ! must be same with H
      CC = coefpot(16.2200d0, 5.5500d0, 1.1924d0)
   case('He')
      CC = coefpot(56.2300D0, 19.240D0, 1.1998D0)
   case('Li')
      CC = coefpot(1.8400d0, 0.7300d0, 2.9081d0)
   case('Be')
      CC = coefpot(2.6100d0, 1.0000d0, 1.5280d0)
   case('B')
      CC = coefpot(6.2100d0, 2.4700d0, 1.6546d0)
   case('C')
      CC = coefpot(9.2800d0, 3.6900d0, 1.5222d0)
   case('N')
      CC = coefpot(12.8700d0, 5.1200d0, 1.4504d0)
   case('O')
      CC = coefpot(18.0900d0, 7.1900d0, 1.4224d0)
   case('F')
      CC = coefpot(23.7800d0, 9.4500d0, 1.3974d0)
   case('Na')
      CC = coefpot(1.7100d0, 0.5000d0, 5.1815d0)
   case('Mg')
      CC = coefpot(2.0400d0, 0.8100d0, 3.5602d0)
   case('Al')
      CC = coefpot(1.7700d0, 0.7000d0, 1.7905d0)
   case('Si')
      CC = coefpot(2.1600d0, 0.8600d0, 1.6054d0)
   case('P')
      CC = coefpot(2.5900d0, 1.0300d0, 1.4995d0)
   case('S')
      CC = coefpot(2.9900d0, 1.1900d0, 1.4261d0)
   case('Cl')
      CC = coefpot(3.4800d0, 1.3800d0, 1.3860d0)
   case('Ar')
      CC = coefpot(3.9900D0, 1.5900D0, 1.3622D0)
   case('K')
      CC = coefpot(3.9900D0, 1.5900D0, 1.3622D0)
      !          CC = coefpot(1.4200d0, 0.2600d0, 6.3140d0)
   case('Ca')
      CC = coefpot(1.6100d0, 0.4500d0, 4.8360d0)
   case('Sc')
      CC = coefpot(3.9600d0, 0.6900d0, 3.7703d0)
   case('Ti')
      CC = coefpot(4.6800d0, 0.9400d0, 3.3889d0)
   case('V')
      CC = coefpot(5.1400d0, 1.1100d0, 2.9680d0)
   case('Cr')
      CC = coefpot(5.1900d0, 1.3700d0, 2.8897d0)
   case('Mn')
      CC = coefpot(6.0300d0, 1.6300d0, 2.7024d0)
   case('Fe')
      CC = coefpot(6.5100d0, 1.9100d0, 2.6179d0)
   case('Co')
      CC = coefpot(6.9500d0, 2.3800d0, 2.7407d0)
   case('Ni')
      CC = coefpot(7.6000d0, 2.7400d0, 2.6949d0)
   case('Cu')
      CC = coefpot(7.5900d0, 3.0200d0, 2.6959d0)
   case('Zu')
      CC = coefpot(8.7800d0, 3.4900d0, 2.6313d0)
   case('Ga')
      CC = coefpot(2.0100d0, 0.8000d0, 4.0433d0)
   case('Ge')
      CC = coefpot(2.2800d0, 0.9100d0, 3.1110d0)
   case('As')
      CC = coefpot(2.6000d0, 1.0300d0, 2.6218d0)
   case('Se')
      CC = coefpot(2.8800D0, 1.1400D0, 2.2934D0)
   case('Br')
      CC = coefpot(3.2000D0, 1.2700D0, 2.1007D0)
   case('Kr')
      CC = coefpot(3.4900D0, 1.3900D0, 1.9478D0)
   case('Rb')
      CC = coefpot(3.4900D0, 1.3900D0, 1.9478D0)
      !          CC = coefpot(1.3700D0, 0.2100D0, 6.8301D0)
   case('Sr')
      CC = coefpot(1.5200D0, 0.3300D0, 4.8514D0)
   case('Y')
      CC = coefpot(2.0600D0, 0.4900D0, 4.1719D0)
   case('Zr')
      CC = coefpot(2.2800D0, 0.6600D0, 3.9162D0)
   case('Nb')
      CC = coefpot(2.4100D0, 0.8200D0, 3.7419D0)
   case('Mo')
      CC = coefpot(2.5700D0, 1.0200D0, 3.8044D0)
   case('Tc')
      CC = coefpot(2.8200D0, 1.1200D0, 3.3669D0)
   case('Ru')
      CC = coefpot(3.0000d0, 1.1900d0, 3.0213d0)
   case('Rh')
      CC = coefpot(3.2100d0, 1.2800d0, 2.7857d0)
   case('Pd')
      CC = coefpot(3.3100d0, 1.3200d0, 2.5256d0)
   case('Ag')
      CC = coefpot(3.5300d0, 1.4100d0, 2.3857d0)
   case('Cd')
      CC = coefpot(3.9100D0, 1.5600D0, 2.3128D0)
   case('In')
      CC = coefpot(1.7900D0, 0.7100D0, 6.7251D0)
   case('Sn')
      CC = coefpot(1.9700D0, 0.7800D0, 5.0086D0)
   case('Sb')
      CC = coefpot(2.1200D0, 0.8500D0, 4.0534D0)
   case('Te')
      CC = coefpot(2.3700D0, 0.9500D0, 3.5696D0)
   case('I')
      CC = coefpot(2.5200D0, 1.0100D0, 3.0856D0)
   case('Xe')
      CC = coefpot(2.6300D0, 1.0500D0, 2.6837D0)
   case('Cs')
      CC = coefpot(2.6300D0, 1.0500D0, 2.6837D0)
      !          CC = coefpot(1.2900d0, 0.1700d0, 7.8924d0)
   case('Pt')
      CC = coefpot(2.7100d0, 1.0800d0, 2.5166d0)
   case('Au')
      CC = coefpot(2.8500d0, 1.1400d0, 2.3778d0)
   case default
      print *, ' FOR ATOMIC NUMBER = ',atom
      stop
   END select

   !        X1   = DSQRT(CCA1)*RAD(I)
   !        X2   = DSQRT(CCA2)*RAD(I)
   !        ERF1 = DERF(X1)
   !        ERF2 = DERF(X2)
   !        IF (RAD(I).EQ.0.0d0) THEN
   !            VCORE(I) = -ZVAL*(CCC1*(DSQRT(CCA1)
   !     &                             -DSQRT(CCA2))+DSQRT(CCA2))
   !          ELSE
   !            VCORE(I) = -ZVAL*(CCC1*(ERF1-ERF2)+ERF2)/RAD(I)
   !        END IF

 END SUBROUTINE VLD
end module m_vld
module m_read_input
  use nrtype

  character(len=10), private, dimension(11) :: exc_type = &
       (/'Wigner    ','HL        ','CA        ','PW91      ','BP86      ', &
       'PBEGGA    ','CAPW91_REL','CAPW91    ','BLYP      ','PWLYP     ', &
       'EXC_ONLY  '/)
!       (/'Wigner','HL','CA','PW91','BP86','PBEGGA','CAPW91_REL', &
!       'CAPW91','BLYP','PWLYP','EXC_ONLY'/)

  !  iexc   exchange(X)-correlation(C) type, LDA's use Dirac exchange
  !          1  LDA   Wigner
  !          2  LDA   Hedin/Lundqvist
  !          3  LDA   Ceperley/Alder Perdew/Zunger (1980)           
  !          4  GGA   Perdew/Wang (1991)
  !          5  GGA   Becke (1988) X, Perdew (1986) C
  !          6  GGA   Perdew/Burke/Ernzerhof (1996)
  !          7  LDA   like 8 + relativistic correction 
  !          8  LDA   Ceperley/Alder Perdew/Wang (1991)
  !          9  GGA   Becke (1988) X, Lee/Yang/Parr (1988) C
  !         10  GGA   Perdew/Wang (1991) X, Lee/Yang/Parr (1988) C
  !         11  LDA   exchange only
  !         default is 8

  real(kind=DP), private :: z, zv, rnlc
  integer(kind=I4B), private :: nc, nv, iexc, l_max, lines
 
  logical, private :: i_am_not_called = .true.

  character(len=1), private :: ps_type ! t = TM, h = Hamann

  character(len=79), private, allocatable :: input(:)

contains
  subroutine read_input(atom)
    implicit none
    character(len=*), intent(in) :: atom 
    integer(kind=I4B) nn, ll, i, ln
    real(kind=DP) ff
    character(len=10) dum
    integer(kind=I4B), parameter :: unit = 10, max_line = 1000
    integer :: values(8)
    character (len=8) :: date
    character (len=10) :: time
    character (len=5) :: zone

    i_am_not_called = .false.

    open(unit,file=trim(atom)//'.ini',form='formatted')
    read(unit,*) z,nc,nv,iexc,rnlc
    zv = z
    do i = 1, nc
       read(unit,*) nn, ll, ff
       zv = zv - ff
    end do
    do i = 1, nv
       read(unit,*) nn, ll, ff
    end do
    read(unit,*)l_max,ps_type
!    print *,'Atom, Zv = ', atom,zv

    rewind unit

    ln = 0
    do while(ln < max_line)
       read(unit,*,end=99,err=99) dum
       ln = ln + 1
    end do
    stop 'input too long'
99  allocate(input(ln+1))
    rewind unit

    do i = 1, ln
       read(unit,'(a78)')input(i)
!       print '(a80)','#'//input(i)
    end do

    call date_and_time(date,time,zone,values)

    input(ln+1) = 'Generated on '//date//', at '//time//', time zone '//zone

    close(unit)

    lines = ln + 1

  end subroutine read_input
  function get_exc_type()
    implicit none
    character(len=10) get_exc_type
    if(i_am_not_called) stop 'i am not called'
    get_exc_type = exc_type(iexc)
  end function get_exc_type
  function get_l_max()
    implicit none
    integer(kind=I4B) get_l_max
    if(i_am_not_called) stop 'i am not called'
    get_l_max = l_max + 1
  end function get_l_max
  function get_zv()
    implicit none
    real(kind=DP) get_zv
    if(i_am_not_called) stop 'i am not called'
    get_zv = zv
  end function get_zv
  function get_ps_type()
    character(len=1) get_ps_type
    get_ps_type = ps_type
  end function get_ps_type
  function get_pcc()
    implicit none
    character(len=3) get_pcc
    if(rnlc > 0.d0) then
       get_pcc = 'yes'
    else
       get_pcc = 'no'
    end if
  end function get_pcc
  function get_lines_header()
    implicit none
    integer get_lines_header
    get_lines_header = lines
  end function get_lines_header
  function get_headers(i)
    implicit none
    integer(kind=I4B), intent(in) :: i
    character(len=80) get_headers
    if(i > size(input)) stop 'Too many lines for the header'
    get_headers = '#'//input(i)
  end function get_headers
end module m_read_input
module m_fitit
  use nrtype
  use nr, only : spline, splint
  
  real(kind=DP), private, allocatable :: ro(:),wo(:,:),po(:,:), y2_p(:,:), &
       y2_w(:,:)

  integer(kind=I4B), private :: l_max, ngrid
  integer(kind=I4B), private, parameter :: unit = 40

contains
  subroutine init_fitit(ll_max,atom,rc_smooth)
    implicit none
    integer(kind=I4B), intent(in) :: ll_max !real one + 1
    character(len=*), intent(in) :: atom
    real(kind=DP), intent(in), optional :: rc_smooth

    real(kind=DP)  x1, x2, x3, rcsmooth, amesh, al, pp
    real(kind=DP) alpha, beta
    integer n, i, j, ism, ll, l, lmax
    real(kind=DP) dd, sum, rdum
    character(len=80) dum

    l_max = ll_max
    
    open(unit,file=trim(atom)//'.cpi',form='formatted')

    do i = 1, 11
       read(unit,*) dum
    end do

    read(unit,*)ngrid, rdum
    print *,ngrid

    backspace unit

    allocate(ro(ngrid), wo(ngrid,l_max), po(ngrid,l_max), &
         y2_p(ngrid,l_max), y2_w(ngrid,l_max))

    if(present(rc_smooth)) then
       rcsmooth = rc_smooth
    else
       rcsmooth = 0.0D0
    end if

    do l = 0, l_max - 1
       !
       ! Read in the data
       !
       ll = l + 1

       
       read(unit,*) n, rdum
       
       if(n /= ngrid) stop 'n != ngrid'
       
       do i=1,n
          read(unit,*) j,x1,x2,x3
          ro(i)=x1
          wo(i,ll)=x2/x1
          po(i,ll)=x3
       end do

!       dd = log(ro(2)/ro(1))
       
!       sum = 0.d0
!       do i = 1, n
!          sum = sum + (ro(i)*wo(i,ll))**2*dd*ro(i)
!       end do

!       print *,'Normalization = ',sum, l

       !
       ! Smooth the pseudopotential?
       !
       if (rcsmooth.gt.0.0D0) then
          amesh=ro(2)/ro(1)
          al=log(amesh)
          do i=1,n
             ism=i
             if(ro(i).gt.rcsmooth) goto 1000
          enddo
          stop 'rcsmooth > maximum r on grid!'
1000      continue
          if (ism.le.2) stop 'rcsmooth < r(2)'
          pp = (2.0D0*po(ism-2,ll)-16.0D0*po(ism-1,ll)+ &
               16.0D0*po(ism+1,ll)-2.0D0*po(ism+2,ll))/(24.0D0*al*ro(ism))
          alpha = pp / (2.0D0 * ro(ism))
          beta  = po(ism,ll)-alpha*ro(ism)*ro(ism)
          do i=1,ism-1
             po(i,ll)=alpha*ro(i)*ro(i)+beta
          enddo
       endif
       !
       ! Pseudopotential
       !

       call spline(ro,po(:,ll),0.0D0,1.0D30,y2_p(:,ll))

       if(l.eq.0) then
          call spline(ro,wo(:,ll),1.0D030,1.0D30,y2_w(:,ll))
       else
          call spline(ro,wo(:,ll),0.0D0,1.0D30,y2_w(:,ll))
       endif
    end do

    close(unit)

  end subroutine init_fitit
  function Phi(r,l)
    implicit none
    real(kind=DP), intent(in) :: r
    integer(kind=I4B), intent(in) :: l
    real(kind=DP) Phi
    integer ll
    ll = l + 1

    if(ro(1) <= r .and. r <= ro(ubound(ro,1))) then
       Phi = splint(ro,wo(:,ll),y2_w(:,ll),r)
    elseif(r < ro(1)) then
       Phi = wo(1,ll)
    else
       stop 'wrong input'
    end if

  end function Phi
  function V_pseudo(r,l)
    implicit none
    real(kind=DP), intent(in) :: r
    integer(kind=I4B), intent(in) :: l
    real(kind=DP) V_pseudo
    integer ll
    ll = l + 1

    if(ro(1) <= r .and. r <= ro(ubound(ro,1))) then
       V_pseudo = splint(ro,po(:,ll),y2_p(:,ll),r)
    elseif(r < ro(1)) then
       V_pseudo = po(1,ll)
    else
       stop 'wrong input'
    end if

  end function V_pseudo
end module m_fitit
module m_potential
  use nrtype
  use nr, only : gauleg
  use m_vld
  use m_read_input
  use m_fitit
  save

  real(kind=DP), private, allocatable :: rr(:), V_ps(:,:), v_nloc(:,:), &
       psi(:,:), v_loc(:), v_lc(:)

  !for spline tables

  real(kind=DP), private, allocatable :: y2a(:)

  ! for gaussian integration
  !  real(kind=DP), private, allocatable :: xx(:), w_x(:)

  integer, private, parameter :: mesh_i = 2000

  real(kind=DP), private ::xx(mesh_i), w_x(mesh_i)

  integer, private :: l_loc

contains
  subroutine init_potential(atom, loc, r_0, r_max)
    implicit none
    character(len=*), intent(in) :: atom
    integer(kind=I4B), intent(in) :: loc
    real(kind=DP), intent(in) :: r_0, r_max

    integer i, j, l, mesh_tmp, mesh
    real(kind=DP) :: gauss, mass, zv
    !for gauleg (gaussian integration)

    character(len=1) dum
    type(coefpot) cc

    real(kind=DP) rrr

!    integer(kind=I4B) :: nc, nv, iexc, nn, ll
!    real(kind=DP) :: z, ff

    l_loc = loc

    call gauleg(r_0,r_max,xx,w_x)

    allocate(psi(mesh_i,get_l_max()),v_lc(mesh_i),v_loc(mesh_i), &
         v_nloc(mesh_i,get_l_max()))

    do i = 1, mesh_i
       v_lc(i) = V_pseudo(xx(i),l_loc)
    end do

    call vld(atom,cc)

    do i = 1, mesh_i
       v_loc(i) = v_lc(i) - vlr(xx(i),cc,get_zv()) !subtract core potential
    end do

    do l = 1, get_l_max()
!       print *,l,get_nonlocal(l-1)
       if(get_nonlocal(l-1)) then
          do i = 1, mesh_i
             v_nloc(i,l) = V_pseudo(xx(i),l-1) - v_lc(i)
!             if(l == 2) then
!                print *, i, xx(i), V_pseudo(xx(i),l-1)-v_lc(i)
!             end if
          end do
       end if
       do i = 1, mesh_i
          psi(i,l) = Phi(xx(i),l-1)
       end do
    end do

  end subroutine init_potential
  function v_ps_moment(l,g,nonlocal)
    implicit none
    integer(kind=I4B), intent(in) :: l
    real(kind=DP), intent(in) :: g
    logical, intent(in) :: nonlocal

    real(kind=DP) v_ps_moment
    integer(kind=I4B) :: i
    real(kind=DP) sum, gr

    if(g == 0.d0) stop 'this routine does not cover g = 0'

    sum = 0.d0
    if(nonlocal) then
       do i = 1, mesh_i
          gr = g*xx(i)
          sum = sum +  xx(i)**2*psi(i,l)*v_nloc(i,l)*bj(gr,l)*w_x(i)
       end do
       sum = sum/g**(l-1)
    else
       do i = 1, mesh_i
          gr = g*xx(i)
          sum = sum +  xx(i)**2*v_loc(i)*bj(gr,1)*w_x(i)
       end do
    end if
    v_ps_moment = sum
  contains
    function bj(gr,l)
      implicit none
      real(kind=DP), intent(in) :: gr
      integer(kind=I4B), intent(in) :: l
      real(kind=DP) bj

      if(l == 1) then
         bj = sin(gr)/gr
      else if(l == 2) then
         bj = sin(gr)/gr**2 - cos(gr)/gr
      else if(l == 3) then
         bj = ((3.d0-gr**2)*sin(gr) - 3.d0*gr*cos(gr))/gr**3
      end if

    end function bj
  end function v_ps_moment
  subroutine print_v_ps_moment(l,nonlocal)
    implicit none
    integer(kind=I4B), intent(in) :: l
    logical, intent(in) :: nonlocal

    integer(kind=I4B) :: i
    character(len=1) nm

    write(nm,'(i1.1)')l
    open(12,file='ps.'//nm,form='formatted')
    if(nonlocal) then
       do i = 1, mesh_i
          write(12,'(2f15.7)') xx(i),psi(i,l)*v_nloc(i,l)
       end do
    else
       do i = 1, mesh_i
          write(12,'(2e15.7)') xx(i),v_loc(i)
       end do
    end if
    close(12)

  end subroutine print_v_ps_moment
  function v_ps_moment_0(l,nonlocal)
    implicit none
    integer(kind=I4B), intent(in) :: l
    logical, intent(in) :: nonlocal

    real(kind=DP) v_ps_moment_0
    integer(kind=I4B) :: i
    real(kind=DP) sum

    sum = 0.d0
    if(nonlocal) then
       do i = 1, mesh_i
          sum = sum +  xx(i)**(l+1)*psi(i,l)*v_nloc(i,l)*bj_0(l)*w_x(i)
       end do
    else
       do i = 1, mesh_i
          sum = sum +  xx(i)**2*v_loc(i)*w_x(i)
       end do
    end if
    v_ps_moment_0 = sum
  contains
    function bj_0(l)
      implicit none
      integer(kind=I4B), intent(in) :: l
      real(kind=DP) bj_0

      if(l == 1) then
         bj_0 = 1.d0
      else if(l == 2) then
         bj_0 = 1.d0/3.d0
      else if(l == 3) then
         bj_0 = 1.d0/15.d0
      end if

    end function bj_0
  end function v_ps_moment_0
  function v_ps_psi2(l)
    implicit none
    integer(kind=I4B), intent(in) :: l

    real(kind=DP) v_ps_psi2
    integer(kind=I4B) :: i
    real(kind=DP) th, ph
    real(kind=DP) sum

!    print *, ubound(v_nloc,1),ubound(v_nloc,2),ubound(psi,1),ubound(psi,2)
!    print *, ubound(xx,1),ubound(w_x,1)
!    print *, v_nloc(82,2),psi(82,2),xx(82),w_x(82)

!    call exit()
    sum = 0.d0
    do i = 1, mesh_i
       sum = sum +  v_nloc(i,l)*(psi(i,l)*xx(i))**2*w_x(i)
    end do

    v_ps_psi2 = sum
  end function v_ps_psi2
  function norm_test(l)
    implicit none
    integer(kind=I4B), intent(in) :: l

    real(kind=DP) norm_test
    integer(kind=I4B) :: i
    real(kind=DP) sum


    sum = 0.d0
    do i = 1, mesh_i
       sum = sum +  xx(i)**2*psi(i,l)**2*w_x(i)
    end do

    norm_test = sum
  end function norm_test
  function get_nonlocal(ll)
    implicit none
    integer(kind=I4B), intent(in) :: ll
    logical get_nonlocal

    get_nonlocal = ll /= l_loc

  end function get_nonlocal
end module m_potential
module m_atomic_data
  use nrtype
  save
  private
  integer(kind=I4B), parameter :: max_atom = 200

  type atomic_data
     character(len=2) :: name
     real(kind=DP) :: mass
  end type atomic_data

  type(atomic_data) ::  atom_data(max_atom)



  public get_atom_mass, init_mass_data
contains
  function get_atom_mass(atm)
    implicit none
    character(len=*), intent(in) :: atm
    integer(kind=I4B) i
    real(kind=DP) get_atom_mass
    do i = 1, ubound(atom_data,1)
       if(atom_data(i)%name == atm) then
!          print *,'atom ',atm,' ',atom_data(i)%mass
          get_atom_mass = atom_data(i)%mass
       end if
    end do
  end function get_atom_mass
  subroutine init_mass_data()
    implicit none
    atom_data(1) = atomic_data('H', 1.0079400d0)
    atom_data(2) = atomic_data('He', 4.0002602d0)
    atom_data(3) = atomic_data('Li', 6.9410000d0)
    atom_data(4) = atomic_data('Be', 9.0121800d0)
    atom_data(5) = atomic_data('B',  10.8220000d0)
    atom_data(6) = atomic_data('C', 12.0110000d0)
    atom_data(7) = atomic_data('N', 14.0076000d0)
    atom_data(8) = atomic_data('O', 15.9994000d0)
    atom_data(9) = atomic_data('F', 18.9984030d0)
    atom_data(10) = atomic_data('Ne', 20.1790000d0)
    atom_data(11) = atomic_data('Na', 22.9897700d0)
    atom_data(12) = atomic_data('Mg', 24.3050000d0)
    atom_data(13) = atomic_data('Al', 26.9815400d0)
    atom_data(14) = atomic_data('Si', 28.0860000d0)
    atom_data(15) = atomic_data('P', 30.9737600d0)
    atom_data(16) = atomic_data('S', 32.0660000d0)
    atom_data(17) = atomic_data('Cl', 35.4530000d0)
    atom_data(18) = atomic_data('Ar', 39.9480000d0)
    atom_data(19) = atomic_data('K', 39.0983000d0)
    atom_data(20) = atomic_data('Ca', 40.0780000d0)
    atom_data(21) = atomic_data('Sc', 44.9559100d0)
    atom_data(22) = atomic_data('Ti', 47.8800000d0)
    atom_data(23) = atomic_data('V', 47.8800000d0)
    atom_data(24) = atomic_data('Cr', 51.9961000d0)
    atom_data(25) = atomic_data('Mn', 54.9380000d0)
    atom_data(26) = atomic_data('Fe', 55.8470000d0)
    atom_data(27) = atomic_data('Co', 58.9332000d0)
    atom_data(28) = atomic_data('Ni', 58.6900000d0)
    atom_data(29) = atomic_data('Cu', 63.5460000d0)
    atom_data(30) = atomic_data('Zn', 65.3900000d0)
    atom_data(31) = atomic_data('Ga', 69.7230000d0)
    atom_data(32) = atomic_data('Ge', 72.5900000d0)
    atom_data(33) = atomic_data('As', 74.9216000d0)
    atom_data(34) = atomic_data('Se', 78.9600000d0)
    atom_data(35) = atomic_data('Br', 79.9040000d0)
    atom_data(36) = atomic_data('Kr', 83.8000000d0)
    atom_data(37) = atomic_data('Rb', 85.4678000d0)
    atom_data(38) = atomic_data('Sr', 87.6200000d0)
    atom_data(39) = atomic_data('Y', 88.9059000d0)
    atom_data(40) = atomic_data('Zr', 91.2240000d0)
    atom_data(41) = atomic_data('Nb', 92.9064000d0)
    atom_data(42) = atomic_data('Mo', 95.9400000d0)
    !  atom_data(43) = atomic_data('Tc', 0.0000000d0)
    atom_data(44) = atomic_data('Ru',101.0700000d0)
    atom_data(45) = atomic_data('Rh', 102.9055000d0)
    atom_data(46) = atomic_data('Pd', 106.4200000d0)
    atom_data(47) = atomic_data('Ag', 107.8682000d0)
    atom_data(48) = atomic_data('Cd', 112.4100000d0)
    atom_data(49) = atomic_data('In', 114.8200000d0)
    atom_data(50) = atomic_data('Sn', 118.7100000d0)
    atom_data(51) = atomic_data('Sb', 121.7500000d0)
    atom_data(52) = atomic_data('Te', 127.6000000d0)
    atom_data(53) = atomic_data('I', 126.9045000d0)
    atom_data(54) = atomic_data('Xe', 131.2900000d0)
    atom_data(55) = atomic_data('Mu', 0.1134940d0)
    atom_data(56) = atomic_data('BM', 28.0860000d0)
    atom_data(101) = atomic_data('D', 2.01400d0) !from BUTSURIGAKUJITEN p892
    atom_data(102) = atomic_data('T', 3.01605d0) !from BUTSURIGAKUJITEN p892
  end subroutine init_mass_data
end module m_atomic_data
module m_ps_conf
  use nrtype

contains
  function get_ps_conf(loc,lmax)
    implicit none
    integer(kind=I4B), intent(in) :: loc, lmax
    character(len=4) :: get_ps_conf, ps_conf = 'eeee'
    integer i

    do i = 0, lmax
       ps_conf(i+1:i+1) = 'n'
    end do

    ps_conf(loc+1:loc+1) = 'l'

    get_ps_conf = ps_conf

  end function get_ps_conf
end module m_ps_conf
program main
  use nrtype
  use m_potential
  use m_read_input
  use m_fitit
  use m_atomic_data
  use m_ps_conf
  implicit none

  real(kind=DP), parameter :: r_0 = 0.d0, r_max = 7.d0
  real(kind=DP) :: rc_smooth = 0.d0

  integer, parameter :: n_g_max = 5001, ngrid = 701
  real(kind=DP), parameter :: g_max = 1d3

  integer :: l_loc
  character(len=2) atom
  character(len=10) exctype, out
  character(len=4) psconf
  character(len=1) name!, ang_mom

!  read(*,*)atom, ang_mom, out
  read(*,*)atom, l_loc, out

!  if(ang_mom == 'S') then
!     l_loc = 0
!  elseif(ang_mom == 'P') then
!     l_loc = 1
!  elseif(ang_mom == 'D') then
!     l_loc = 2
!  else
!     stop 'INVALID ANGULAR MOMENTUM STATE IS SPECIFIED!'
!  end if

  call read_input(atom)

  if(l_loc + 1 > get_l_max()) stop 'INVALID CHOICE OF LOCAL POT'

  psconf = get_ps_conf(l_loc,get_l_max()-1)
  print *,'Exchange type: ',get_exc_type()
  print *,'PS CONF: ',psconf
  print *,'NON LINEAR CORE: ',get_pcc()

  open(10,file='fit.prm',form='formatted')
  read(10,*)rc_smooth
  close(10)
  call init_fitit(get_l_max(),atom,rc_smooth) !read ps data and set up function

  call init_mass_data()

  if(out == 'MP_DFT') then
     !MP_DFT

     call generate_MP_DFT(atom,r_0,r_max,g_max,n_g_max,l_loc,get_exc_type(), &
          psconf,get_pcc())

  elseif(out == 'JEEP') then
     !JEEP fomat

     call generate_JEEP(atom,r_0,r_max,ngrid,l_loc,get_exc_type())

  elseif(out == 'BOTH') then

     call generate_MP_DFT(atom,r_0,r_max,g_max,n_g_max,l_loc,get_exc_type(), &
          psconf,get_pcc())
     call generate_JEEP(atom,r_0,r_max,ngrid,l_loc,get_exc_type())

  else

     stop 'INVALID OUTPUT OPTION'
  end if

contains
  subroutine generate_MP_DFT(atom,r_0,r_max,g_max,n_g_max,l_loc,exctype,&
       psconf,pcc)
    implicit none
    character(len=*), intent(in) :: atom, exctype, psconf, pcc
    real(kind=DP), intent(in) :: r_0, r_max, g_max
    integer, intent(in) :: n_g_max, l_loc

    integer l, ig, i
    real(kind=DP) g
    real(kind=DP), allocatable :: v_ps(:,:), v_denom(:), g_lst(:)

    call init_potential(atom,l_loc,r_0,r_max)

    ! data on gaussian int mesh is made.
    !l_loc should be a real angular momentum number

    allocate(v_ps(n_g_max,get_l_max()),v_denom(get_l_max()),g_lst(n_g_max))

    v_denom = 0.d0
    do l = 1, get_l_max()
       v_ps(1,l) = v_ps_moment_0(l,get_nonlocal(l-1))
       g_lst(1) = 0.d0
       do ig = 1, n_g_max - 1
          g = real(ig,DP)*sqrt(g_max)/real(n_g_max-1,DP)
          g_lst(ig+1) = g
          v_ps(ig+1,l) = v_ps_moment(l,g,get_nonlocal(l-1))
       end do
       if(get_nonlocal(l-1)) v_denom(l) = v_ps_psi2(l)
       call print_v_ps_moment(l,get_nonlocal(l-1))
!       print *,'Norm = ', norm_test(l)
    end do
     
    open(10,file=trim(exctype)//'.'//trim(atom),form='formatted')

    do i = 1, get_lines_header()
       write(10,'(a80)') get_headers(i)
    end do

    write(10,'(i3,1x,a2,1x,f5.1,1x,f13.7,1x,a4,1x,a10,1x,a3)') 0, atom, &
         get_zv(), get_atom_mass(atom),psconf, exctype, pcc
    write(10,'(D24.12,2i8)')g_max, n_g_max, get_l_max()
    write(10,'(3D24.12)')((v_ps(ig,l),ig=1,n_g_max),l = 1, get_l_max())
    write(10,'(3D24.12)')(v_denom(l),l=1,get_l_max())
    close(10)
     
!    open(11,file=trim(atom)//'_ps.dat',form='formatted')
!    do ig = 1, n_g_max
!       write(11,'(4f15.7)'),g_lst(ig),(v_ps(ig,l),l=1,get_l_max())
!    end do
!    close(11)

  end subroutine generate_MP_DFT
  subroutine generate_JEEP(atom,r_0,r_max,ngrid,l_loc,exctype)
    implicit none
    character(len=*), intent(in) :: atom, exctype
    real(kind=DP), intent(in) :: r_0, r_max
    integer, intent(in) :: ngrid, l_loc
    
    
    integer l, i
    real(kind=DP) :: rr, d_rr, ps, wv
    character(len=20) ps_name

    if(get_ps_type() == 't') then
       ps_name = 'TM_'//trim(exctype)//'_'//trim(atom)
    elseif(get_ps_type() == 'h') then
       ps_name = 'Hamann_'//trim(exctype)//'_'//trim(atom)
    else
       stop 'Invalid ps type'
    end if

    open(10,file=trim(ps_name),form='formatted')
    do i = 1, get_lines_header()
       write(10,'(a80)') get_headers(i)
    end do

    write(10,'(i5,1x,i3,1x,f4.1,1x,f13.7,1x,i2,1x,a20,1x,a13,i2,1x,i2)') &
         ngrid, int(get_zv(),I2B),1.0,get_atom_mass(atom),get_l_max(), &
         ps_name,'Green 1.0 1.0',l_loc,1 !the last number for KB?
    do l = 1, get_l_max()
       write(name,'(i1.1)')l-1
       open(11,file=trim(atom)//'_ps.'//name,form='formatted')
       open(12,file=trim(atom)//'_wv.'//name,form='formatted')
       if(l /= 1) then
          write(10,*) ngrid
       end if
       d_rr = (r_max - r_0)/real(ngrid-1,DP)
       rr = 0.d0
       do i = 1, ngrid
          ps = V_pseudo(rr,l-1)
          write(10,'(2e23.12)')rr, ps
          write(11,'(2e23.12)')rr, ps
          rr = rr + d_rr
       end do
       write(10,*) ngrid
       rr = 0.d0
       do i = 1, ngrid
          wv = Phi(rr,l-1)
          write(10,'(2e23.12)')rr, wv
          write(12,'(2e23.12)')rr, wv
          rr = rr + d_rr
       end do
       close(11)
       close(12)
    end do
    close(10)
  end subroutine generate_JEEP
end program main
