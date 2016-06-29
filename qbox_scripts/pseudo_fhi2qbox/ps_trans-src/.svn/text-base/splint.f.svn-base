      subroutine splint(xa,ya,y2a,n,x,y)
      implicit none
      integer n
      real*8 x,y,xa(n),y2a(n),ya(n)
      integer k,khi,klo
      real*8 a,b,h

      if(x .le. xa(1)) then
         y = ya(1)
         return
      end if

      klo=1
      khi=n
 1000 if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x) then
          khi=k
        else 
          klo=k
        endif
        goto 1000
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0) then
        write(*,*) "Bad xa in splint"
        stop
      endif
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+
     >  ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0D0
      return
      end


      subroutine splint2(xa,ya,y2a,n,nx,x,y)
c
      double precision a, b, h, x, xa, y, y2a, ya
      integer khi, klo, n, nx
c
      dimension xa(n),ya(n),y2a(n)
      khi=nx+1
      klo=nx
      if(xa(klo) .gt. x .or. xa(khi) .lt. x) stop 'splint'
      h=xa(khi)-xa(klo)
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+
     &      ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.d0
      return
      end

