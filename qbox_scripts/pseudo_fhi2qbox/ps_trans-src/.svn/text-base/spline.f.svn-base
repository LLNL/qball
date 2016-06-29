      subroutine spline2(x,y,n,yp1,ypn,y2)
      integer n,nmax
      real*8 yp1,ypn,x(n),y(n),y2(n)
      parameter (nmax=500)
      integer i,k
      real*8 p,qn,sig,un,u(nmax)
      if(yp1.gt.0.99D30) then
        y2(1)=0.0D0
        u(1)=0.0D0
      else
        y2(1)=-0.5D0
        u(1)=(3.0D0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.0D0
        y2(i)=(sig-1.0D0)/p
        u(i)=(6.0D0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     >       /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      enddo
      if (ypn.gt.0.99D30) then 
        qn=0.0D0
        un=0.0D0
      else
        qn=0.5D0
        un=(3.0D0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.0D0)
      do k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
      enddo
      return
      end


      subroutine spline(x,y,n,yp1,ypn,y2)
      implicit none
c
      double precision p, qn, sig, u, un, x, y, y2, yp1, ypn
      integer i, k, n, nmax
c
c
      parameter (nmax=600)
      dimension x(n),y(n),y2(n),u(nmax)
      if (yp1 .gt. 0.99d30) then
        y2(1)=0.d0
        u(1)=0.d0
      else
        y2(1)=-0.5d0
        u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.d0
        y2(i)=(sig-1.d0)/p
        u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     *      /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
11    continue
      if (ypn.gt..99d30) then
        qn=0.d0
        un=0.d0
      else
        qn=0.5d0
        un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      end

