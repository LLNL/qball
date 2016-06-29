module m_fitit
  use m_dft_op_types_def
  use nr, only : spline, splint
  
  real(kind=DP), private, allocatable :: ro(:),wo(:,:),po(:,:), y2_p(:,:), &
       y2_w(:,:)

  integer(kind=INTG), private :: l_max, ngrid
  integer(kind=INTG), private, parameter :: baseunit = 40

contains
  subroutine init_fitit(ll_max)
    implicit none
    integer(kind=INTG), intent(in) :: ll_max !real one + 1
    character(len=1) dum
    integer(kind=INTG), intent(in) ::  ! this is real angular momentum number
    real(kind=DP  x1, x2, x3, rcsmooth, amesh, al, pp
    real(kind=DP) alpha, beta
    integer n, i, j, ism, ll, l, lmax, iunit
    real(kind=DP) dd, sum

    l_max = ll_max
    
    read(baseunit,*) dum, ngrid

    backspace baseunit

    print *,ngrid

    allocate(ro(ngrid), wo(ngrid,l_max), po(ngrid,l_max), &
         y2_p(ngrid,l_max), y2_w(ngrid,l_max))


    rcsmooth = 0.02D0

    do l = 0, l_max - 1
       !
       ! Read in the data
       !
       ll = l + 1
       iunit = baseunit + l
       
       read(iunit,*) dum, n
       
       if(n /= ngrid) stop 'n != ngrid'
       
       do i=1,n
          read(iunit,*) j,x1,x2,x3
          ro(i)=x1
          wo(i,ll)=x2/x1
          po(i,ll)=x3
       end do

       dd = log(r(2)/r(1))
       
       sum = 0.d0
       do i = 1, n
          sum = sum + (ro(i)*wo(i,ll))**2*dd*ro(i)
       end do

       print *,'Normalization = ',sum, l

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
             po(i,ll)=alpha*ro(i,ll)*ro(i,ll)+beta
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
  end subroutine init_fitit
  function Phi(r,l)
    implicit none
    real(kind=DP), intent(in) :: r
    integer(kind=INTG), intent(in) :: l
    real(kind=DP) Phi
    integer ll
    ll = l + 1

    if(r >= ro(1)) then
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
    integer(kind=INTG), intent(in) :: l
    real(kind=DP) V_pseudo
    integer ll
    ll = l + 1

    if(r >= ro(1)) then
       V_pseudo = splint(ro,po(:,ll),y2_p(:,ll),r)
    elseif(r < ro(1)) then
       V_pseudo = po(1,ll)
    else
       stop 'wrong input'
    end if

  end function V_pseudo
end module m_fitit
