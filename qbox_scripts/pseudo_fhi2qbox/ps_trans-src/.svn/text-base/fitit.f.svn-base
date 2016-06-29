      program fitit
      implicit none
      integer maxn,ngrid,maxl
      parameter (maxn=1000)
      parameter (ngrid=700)
      parameter (maxl=3)
      real*8 x1,x2,x3,delta,rcsmooth,amesh,al,pp
      real*8 alpha,beta
      real*8 r(maxn),w(maxn),p(maxn),y2(maxn)
      real*8 ro(ngrid),wo(ngrid,maxl),po(ngrid,maxl)
      integer n,i,j,ism,l,ll,lmax,iunit
      character dum*1, atom*2

      real*8 dd, sum

      delta = 0.01D0
      rcsmooth = 0.02D0

c
c Read in the data
c
      write(*,*)'l max?'
      read(*,*)lmax
      write(*,*)'atom?'
      read(*,*)atom

      open(70,file='pseudo.'//trim(atom),form='formatted')
      l = 0
      do ll = 1, lmax
      iunit = 39 + ll
      read(iunit,*) dum, n
      if (n.gt.maxn) stop 'maxn too small'
      do i=1,n
        read(iunit,*) j,x1,x2,x3
        r(i)=x1
        w(i)=x2/x1
        p(i)=x3
      enddo


      dd = log(r(2)/r(1))

      sum = 0.d0
      do i = 1, n
         sum = sum + (r(i)*w(i))**2*dd*r(i)
      end do

      print *,'Norm = ',sum,l
      pause

c
c Smooth the pseudopotential?
c
      if (rcsmooth.gt.0.0D0) then
        amesh=r(2)/r(1)
        al=log(amesh)
        do i=1,n
          ism=i
          if(r(i).gt.rcsmooth) goto 1000
        enddo
        stop 'rcsmooth > maximum r on grid!'
 1000   continue
        if (ism.le.2) stop 'rcsmooth < r(2)'
        pp = (2.0D0*p(ism-2)-16.0D0*p(ism-1)+
     >        16.0D0*p(ism+1)-2.0D0*p(ism+2))/(24.0D0*al*r(ism))
        alpha = pp / (2.0D0 * r(ism))
        beta  = p(ism)-alpha*r(ism)*r(ism)
        do i=1,ism-1
          p(i)=alpha*r(i)*r(i)+beta
        enddo
      endif
c
c Pseudopotential
c
      call spline(r,p,n,0.0D0,1.0D30,y2)
      ro(1)=0.0D0
      po(1,ll)=p(1)
c      x1=delta
c      do i=2,ngrid
      x1=0
      do i=1,ngrid
        call splint(r,p,y2,n,x1,x2)
        ro(i)=x1
        po(i,ll)=x2
        x1=x1+delta
      enddo
c
c Wavefunction
c
      if(l.eq.0) then
        call spline(r,w,n,1.0D030,1.0D30,y2)
      else
        call spline(r,w,n,0.0D0,1.0D30,y2)
      endif
      x1=0.0D0
      do i=1,ngrid
        call splint(r,w,y2,n,x1,x2)
        wo(i,ll)=x2!/x1
        x1=x1+delta
      enddo
      if (l.eq.0) then
        wo(1,ll)=w(1)!/r(1)
      else
        wo(1,ll)=0.0D0
      endif

      sum = 0.d0
      do i = 1, ngrid
         sum = sum + (ro(i)*wo(i,ll))**2*delta
      end do
!      print *,'Norm on grid = ', sum

c
c Write to disk
c
      if(l == 0) then
         write(70,*) ngrid,lmax
      else
         write(70,*) ngrid
      endif

      do i=1,ngrid
        write(70,9000) ro(i),po(i,ll)
      enddo
      write(70,*) ngrid
      do i=1,ngrid
        write(70,9000) ro(i),wo(i,ll)
      enddo


      l = l + 1
      end do ! l loop
      close(70)

      do i=1,ngrid
         write(71,9100) ro(i),(po(i,ll),ll=1,lmax)
      enddo

      do i=1,ngrid
         write(72,9100) ro(i),(wo(i,ll)*ro(i),ll=1,lmax)
      enddo

      stop
 9000 format(f10.6,3x,e23.15)
 9100 format(f10.6,3x,3e23.15)
      end 
