C Finite difference fractional step code 
C Staggered grid
C
c Lid velocity at y=H is 1.0
c Lid height H = 1.0
c Square cavity of dimension 1 by 1
c Reynolds number is an input parameter
C Verzicco & Orlandi's method used for pressure BC, etc.
C Time integration by 3rd-order Runge-Kutta
C
C References:
C  (1) Verzicco R & Orlandi P, A finite-difference scheme for 
C      three-dimensional incompressible flows in cylindrical coordinates.
C      J. Comp. Phys. 123, 402-414 (1996).
C
C  (2) M. M. Rai, and P. Moin, Direct simulations of turbulent flow using
C      finite-difference schemes. J. Comput. Phys. 96, 15 (1991).
C
C  (3) Kim J and Moin P, Application of a fractional-step method to incompressibe
C     Navier-Stokes Equations. J. Comp. Phys. 59, 308-323 (1985).
C ------------------------------------------------------

       program incomp2d
* the following is the contents of the include file:      
*      parameter (n1=32,n1p2=n1+2,n1p1=n1+1,n1m1=n1-1,
*     >           n2=32,n2p2=n2+2,n2p1=n2+1,n2m1=n2-1)
      include 'param.inc'
        PARAMETER ( n1h=n1/2 , n2h=n2/2 )
      common/one/qx(n1p1,n2p2),qy(n1p2,n2p1),dx,dy,dx2,dy2,
     >           dt,re,aaa(3),ggg(3),rrr(3)
      common/two/hx(2:n1,2:n2p1),hy(2:n1p1,2:n2),
     >           hxp(2:n1,2:n2p1),hyp(2:n1p1,2:n2)
      common/three/ax(2:n1,2:n2p1),ay(2:n1p1,2:n2)
      common/six/pr(n1,n2),phi(n1,n2),phip(n1,n2),prp(n1,n2)
      common/seven/x(n1p1),y(n2p1),xc(n1),yc(n2)
      dimension qxold(n1p1,n2p2)
        INTEGER IFAXX(13)
        REAL TRIGSX(3*n1h+1)
        COMMON /TRIG/  TRIGSX
        COMMON /FAX/   IFAXX
      real len
c     NAMELIST/INPUT/niter,iread,re,len,sar,dt,rerror
      NAMELIST/INPUT/niter,nprof,iread,re,len,dt,rerror
* FFT setup
      CALL FFTFAX(n1,IFAXX,TRIGSX)
      open(10,file='velprof1.dat',form='formatted',
     >     status='unknown')
      open(11,file='velprof2.dat',form='formatted',
     >     status='unknown')
      open(25,file='uc.dat',form='formatted',status='unknown')
      open(100,file='eps.dat',form='formatted',status='unknown')
      open(5,file='data.inp',form='formatted',status='old')
      open(7,file='data.out',form='formatted',status='unknown')
      open(201,file='restart.flo',form='unformatted',
     >          status='unknown')
      open(202,file='endrun.flo',form='unformatted',
     >          status='unknown')
* DX(),DY(),arrays are not needed for dx,dy uniform
* here, x,y are locations of cell boundarys (this will change 
* for dx,dy not uniform
*
      zero = 1.0e-30
c set up Runge-Kutta coefficients
      ggg(1) = 8./15.
      ggg(2) = 5./12.
      ggg(3) = 3./4.
      rrr(1) = 0.
      rrr(2) = -17./60.
      rrr(3) = -5./12.      
      aaa(1) = 8./15.
      aaa(2) = 2./15.
      aaa(3) = 1./3.
* Len = length,width of cavity
* input Reynolds number
      read(5,INPUT)
      write(7,INPUT)
c if IREAD = 0, new start (else RESTART)
      IF (iread .eq. 0) THEN
      itbeg = 0
      tau = 0.0
      dx = len / float(n1)
      dy = len / float(n2)
      dx2 = dx ** 2.0
      dy2 = dy ** 2.0
* fill x,y vectors (grid line coords)
      do i=1,n1p1
        x(i) = float(i)*dx - dx
      enddo
      do j=1,n2p1
        y(j) = float(j)*dy - dy
      enddo

      do i=1,n1p1
      write(21,*)i,x(i),y(i)
      enddo
* fill the cell center x,y vectors 
      xc(1) = dx/2.
      do i=2,n1
        xc(i) = xc(i-1) + dx
      enddo
      yc(1) = dy/2.
      do j=2,n2
        yc(j) = yc(j-1) + dy
      enddo
*
* stability criterion (Peyret & Taylor, 1983)...
* there are two: eqs. 6.2.10 & 6.2.11
* for eq. 6.2.10, assume |u0| & |v0| are equal to 1.0
*
      c1 = dt * re
      write(*,*) 'Stability criterion 1 = ',c1
      if(c1 .gt. 1.0) then
        write(*,*) 'Stability criterion #1 greater than 1.0'
c        write(*,*) 'Program stopped'
c        stop
      endif
      c2 = 4*dt/(re*dx2)
      write(*,*) 'Stability criterion 2 = ',c2              
      if(c2 .gt. 1.0) then
        write(*,*) 'Stability criterion #2 greater than 1.0'
c        write(*,*) 'Program stopped'
c        stop
      endif
*
* initialize all arrays at t=zero...
* (can include boundary for velocity, but make general 
* for boundary velocity not equal to zero)
* for this case, hxp,hyp are initially zero
*
      do j=2,n2p1
        do i=2,n1
          qx(i,j) = zero
          hxp(i,j) = zero
        enddo
      enddo
      do i=2,n1p1
        do j=2,n2
          qy(i,j) = zero
          hyp(i,j) = zero
        enddo
      enddo
      do j=1,n2
        do i=1,n1
          pr(i,j) = zero
          phi(i,j) = zero
        enddo
      enddo
*
* fill bc for velocities defined along the solid boundary within
* the domain
      do j=2,n2p1
        qx(1,j) = zero
        qx(n1p1,j) = zero
      enddo
      do i=2,n1p1
        qy(i,1) = zero
        qy(i,n2p1) = zero
      enddo
*
* fill bc for velocities defined outside the domain...
      do i=1,n1p1
        qx(i,1) = -qx(i,2)
        qx(i,n2p2) = 2.0 - qx(i,n2p1)
      enddo
      do j=1,n2p1
        qy(1,j) = -qy(2,j)
        qy(n1p2,j) = -qy(n1p1,j)
      enddo
*
c fill in the corner points
c              don't need for 2-D
c
c
      ELSE
      read(201) itbeg
      read(201) tau,dx,dy,dx2,dy2
      read(201) x,y,xc,yc
      read(201) qx,qy,hxp,hyp
      read(201) pr,phi
      ENDIF
*
* here is the time-stepping loop...
*
      itbeg = itbeg+1
      niter = niter+itbeg-1
      do itime=itbeg,niter
*
      if(mod(itime-1,nprof).eq.0)then
       write(*,*)'itime-1=', itime-1, 'writing'
       call output
      end if
*      
      write(25,*)tau,0.5*(qx(49,49)+qx(49,50)),
     >    0.5*(qy(49,49)+qy(50,49))

      tau = tau + dt
*
* save old qx,qy values to determine if steady-state 
* has been reached.
*
      do j=1,n2p2
        do i=1,n1p1
          qxold(i,j) = qx(i,j)
        enddo
      enddo
c
c start of Runge-Kutta steps (3)...
c
      do irk=1,3      
*
* calculate non-linear (hx,hy) and viscous (ax,ay) terms
* (sub-step 1)
*
        call hxycalc
        call axycalc
*
* calculate intermediate velocities..
* (sub-step 2)
        call qhat(irk)
*
* calculate the phi from the poisson equation...
* (sub-step 3)
        call trigon(irk)
* 
* calculate the final velocities (at current t-step)...
* (sub-step 4)
        call ufinal(irk)
*
* update the pressure (assume wanted)...
* (sub-step 5)
        call press(irk)
*
c end of Runge-Kutta steps...
c
      enddo
* determine if steady-state has been reached...
*
      epsmax = 0.
      sum = 0.
      sumtot = 0.
      do j=1,n2p2
        do i=1,n1p1
          sum = sum + (qx(i,j) - qxold(i,j))**2.
          sumtot = sumtot + (qx(i,j))**2.
        enddo
      enddo
      epsmax = sqrt(sum/sumtot)
* write steady state parameter to file      
      if(mod(itime,50).eq.0.) write(100,*) tau,epsmax
      if(mod(itime,10).eq.0.) write(*,*) itime,tau,epsmax
      if(itime .gt. 5 .and. epsmax .le. rerror) then
        write(*,*) 'Steady state reached at tau = ',tau
        write(*,*) 'Iteration Number: ',itime
        goto 999
      endif
*
* here is the end of the time-stepping loop...
      enddo        
*
  999 continue
c
c write restart file:
c
      write(202) niter
      write(202) tau,dx,dy,dx2,dy2
      write(202) x,y,xc,yc
      write(202) qx,qy,hxp,hyp
      write(202) pr,phi
*
      end
*
*******************************************************
*
      subroutine hxycalc
      include 'param.inc'
* calculate the non-linear terms (x-dir)
      common/one/qx(n1p1,n2p2),qy(n1p2,n2p1),dx,dy,dx2,dy2,
     >           dt,re,aaa(3),ggg(3),rrr(3)
      common/two/hx(2:n1,2:n2p1),hy(2:n1p1,2:n2),
     >           hxp(2:n1,2:n2p1),hyp(2:n1p1,2:n2)
*
      do j=2,n2p1
        do i=2,n1
          temp1 = qx(i+1,j)**2.0 - qx(i-1,j)**2.0
          temp2 = (qx(i,j+1)+qx(i,j))*(qy(i+1,j)+qy(i,j))
     >            -(qx(i,j)+qx(i,j-1))*(qy(i+1,j-1)+qy(i,j-1))
          hx(i,j) = -temp1/(2.*dx) - temp2/(4.*dy)
        enddo
      enddo
*
      do j=2,n2
        do i=2,n1p1
          temp1 = (qx(i,j+1)+qx(i,j))*(qy(i+1,j)+qy(i,j))
     >            -(qx(i-1,j+1)+qx(i-1,j))*(qy(i-1,j)+qy(i,j))
          temp2 = qy(i,j+1)**2.0 - qy(i,j-1)**2.0
          hy(i,j) = -temp1/(4.*dx) -temp2/(2.*dy)
        enddo
      enddo
*
      return
      end
*
*********************************************************            
*
      subroutine axycalc
      include 'param.inc'
c calculate viscous terms      
      common/one/qx(n1p1,n2p2),qy(n1p2,n2p1),dx,dy,dx2,dy2,
     >           dt,re,aaa(3),ggg(3),rrr(3)
      common/three/ax(2:n1,2:n2p1),ay(2:n1p1,2:n2)
c
      cc = 1.0/(2.0*re)
      ccx = cc/dx2
      ccy = cc/dy2
c      
      do j=2,n2p1
        do i=2,n1
          ax1 = qx(i+1,j) + qx(i-1,j) - 2.*qx(i,j)
          ax2 = qx(i,j+1) + qx(i,j-1) - 2.*qx(i,j)
          ax(i,j) = ax1*ccx + ax2*ccy
        enddo
      enddo
      do j=2,n2
        do i=2,n1p1
          ay1 = qy(i+1,j) + qy(i-1,j) - 2.*qy(i,j)
          ay2 = qy(i,j+1) + qy(i,j-1) - 2.*qy(i,j)
          ay(i,j) = ay1*ccx + ay2*ccy
        enddo
      enddo
*
      return
      end
*
*********************************************************            
*
      subroutine qhat(irk)
      include 'param.inc'
* calculate the intermediate-time velocities
      common/one/qx(n1p1,n2p2),qy(n1p2,n2p1),dx,dy,dx2,dy2,
     >           dt,re,aaa(3),ggg(3),rrr(3)
      common/two/hx(2:n1,2:n2p1),hy(2:n1p1,2:n2),
     >           hxp(2:n1,2:n2p1),hyp(2:n1p1,2:n2)
      common/three/ax(2:n1,2:n2p1),ay(2:n1p1,2:n2)
      common/six/pr(n1,n2),phi(n1,n2),phip(n1,n2),prp(n1,n2)
*note: max dim. on a,b,c,d is max(n1p1,n2p1)...     
      dimension a(n1p1),b(n1p1),c(n1p1),d(n1p1)
*
* the qx's
c          
      do j=2,n2p1
        do i=2,n1
          ax(i,j) = dt*(ggg(irk)*hx(i,j)+rrr(irk)*hxp(i,j)
     >              - aaa(irk)*(pr(i,j-1)-pr(i-1,j-1))/dx)
     >                 + 2.*aaa(irk)*dt*ax(i,j)
        enddo
      enddo
*
* set up vectors for tri-diag calculation
*
      aa = -dt*aaa(irk)/(2.*re*dx2)
      bb = 1. + dt*aaa(irk)/(re*dx2)
      do j=2,n2p1
        do i=2,n1
* note: a,b,c are constant here - but leave for generality        
          a(i) = aa
          b(i) = bb
          c(i) = a(i)
          d(i) = ax(i,j)
        enddo
        call sy(2,n1,a,b,c,d)
        do i=2,n1
          ax(i,j) = d(i)
        enddo
      enddo
      aa = -dt*aaa(irk)/(2.*re*dy2)
      bb = 1. + dt*aaa(irk)/(re*dy2)
      do i=2,n1
        do j=2,n2p1
          a(j) = aa
          b(j) = bb
          c(j) = a(j)
          d(j) = ax(i,j)
        enddo
        call sy(2,n2p1,a,b,c,d)
        do j=2,n2p1
          ax(i,j) = d(j)
          qx(i,j) = ax(i,j) + qx(i,j)
        enddo
      enddo
*   
* the qy's
*
        do j=2,n2
          do i=2,n1p1
            ay(i,j) = dt*(ggg(irk)*hy(i,j)+rrr(irk)*hyp(i,j)
     >              - aaa(irk)*(pr(i-1,j)-pr(i-1,j-1))/dy)
     >                 + 2.*aaa(irk)*dt*ay(i,j)
          enddo
        enddo
* 
* set up the vectors for tri-diag calculation
      aa = -dt*aaa(irk)/(2.*re*dx2)
      bb = 1. + dt*aaa(irk)/(re*dx2)
      do j=2,n2
        do i=2,n1p1
          a(i) = aa
          b(i) = bb
          c(i) = a(i)
          d(i) = ay(i,j)
        enddo
        call sy(2,n1p1,a,b,c,d)
        do i=2,n1p1
          ay(i,j) = d(i)
        enddo
      enddo
      aa = -dt*aaa(irk)/(2.*re*dy2)
      bb = 1. + dt*aaa(irk)/(re*dy2)
      do i=2,n1p1
        do j=2,n2
          a(j) = aa
          b(j) = bb
          c(j) = a(j)
          d(j) = ay(i,j)
        enddo
        call sy(2,n2,a,b,c,d)
        do j=2,n2
          ay(i,j) = d(j)
          qy(i,j) = ay(i,j) + qy(i,j)
        enddo
      enddo
*
* save previous hx,hy values in hxp,hyp
*
      do j=2,n2p1
        do i=2,n1
          hxp(i,j) = hx(i,j)
        enddo
      enddo
      do j=2,n2
        do i=2,n1p1
          hyp(i,j) = hy(i,j)
        enddo
      enddo
*      
      return
      end
*                
*********************************************************            
*
      subroutine ufinal(irk)
      include 'param.inc'
      common/one/qx(n1p1,n2p2),qy(n1p2,n2p1),dx,dy,dx2,dy2,
     >           dt,re,aaa(3),ggg(3),rrr(3)
      common/six/pr(n1,n2),phi(n1,n2),phip(n1,n2),prp(n1,n2)
*
* calculate the vel. at t=n+1...
      do j=2,n2p1
        do i=2,n1
          qx(i,j)=qx(i,j)-(ggg(irk)*(phi(i,j-1)-phi(i-1,j-1))
     >            + rrr(irk)*(phip(i,j-1)-phip(i-1,j-1)))*dt/dx
        enddo
      enddo
*
      do j=2,n2
        do i=2,n1p1
          qy(i,j)=qy(i,j)-(ggg(irk)*(phi(i-1,j)-phi(i-1,j-1))
     >            + rrr(irk)*(phip(i-1,j)-phip(i-1,j-1)))*dt/dy
        enddo
      enddo
*
* fill bc for velocities defined outside the domain...
      do i=2,n1
        qx(i,1) = -qx(i,2)
        qx(i,n2p2) = 2.0 - qx(i,n2p1)
      enddo
      do j=2,n2
        qy(1,j) = -qy(2,j)
        qy(n1p2,j) = -qy(n1p1,j)
      enddo
*
c fill in the corner points
c              don't need for 2-D
c
*
      return
      end
*                
*********************************************************            
*
      subroutine press(irk)
      include 'param.inc'
      common/one/qx(n1p1,n2p2),qy(n1p2,n2p1),dx,dy,dx2,dy2,
     >           dt,re,aaa(3),ggg(3),rrr(3)
      common/two/hx(2:n1,2:n2p1),hy(2:n1p1,2:n2),
     >           hxp(2:n1,2:n2p1),hyp(2:n1p1,2:n2)
      common/three/ax(2:n1,2:n2p1),ay(2:n1p1,2:n2)
      common/six/pr(n1,n2),phi(n1,n2),phip(n1,n2),prp(n1,n2)
*
      dt2re = dt/(2.*re)
*     
* internal... 
      do j=2,n2m1
        do i=2,n1m1
          pr(i,j)=ggg(irk)*(phi(i,j)-dt2re*
     >            ((phi(i-1,j)-2.*phi(i,j)+phi(i+1,j))/dx2 +
     >             (phi(i,j-1)-2.*phi(i,j)+phi(i,j+1))/dy2)) +
     >            rrr(irk)*(phip(i,j)-dt2re*pr(i,j))
          pr(i,j)=prp(i,j) + pr(i,j)/aaa(irk)
        enddo
      enddo
* sides...      
      do j=2,n2m1
        pr(1,j)=ggg(irk)*(phi(1,j)-dt2re*
     >          ((phi(2,j)-phi(1,j))/dx2 +
     >           (phi(1,j-1)-2.*phi(1,j)+phi(1,j+1))/dy2)) +
     >          rrr(irk)*(phip(1,j)-dt2re*pr(1,j))
        pr(1,j)=prp(1,j) + pr(1,j)/aaa(irk)
        pr(n1,j)=ggg(irk)*(phi(n1,j)-dt2re*
     >           ((phi(n1m1,j)-phi(n1,j))/dx2 +
     >            (phi(n1,j-1)-2.*phi(n1,j)+phi(n1,j+1))/dy2)) +
     >           rrr(irk)*(phip(n1,j)-dt2re*pr(n1,j))
        pr(n1,j)=prp(n1,j) + pr(n1,j)/aaa(irk)
      enddo
      do i=2,n1m1
        pr(i,1)=ggg(irk)*(phi(i,1)-dt2re*
     >          ((phi(i-1,1)-2.*phi(i,1)+phi(i+1,1))/dx2 +
     >           (phi(i,2)-phi(i,1))/dy2)) +
     >          rrr(irk)*(phip(i,1)-dt2re*pr(i,1))
        pr(i,1)=prp(i,1) + pr(i,1)/aaa(irk)
        pr(i,n2)=ggg(irk)*(phi(i,n2)-dt2re*
     >           ((phi(i-1,n2)-2.*phi(i,n2)+phi(i+1,n2))/dx2 +
     >            (phi(i,n2m1)-phi(i,n2))/dy2)) +
     >           rrr(irk)*(phip(i,n2)-dt2re*pr(i,n2))
        pr(i,n2)=prp(i,n2) + pr(i,n2)/aaa(irk)
      enddo
* corners
        pr(1,1)=ggg(irk)*(phi(1,1)-dt2re*
     >          ((phi(2,1)-phi(1,1))/dx2 +
     >           (phi(1,2)-phi(1,1))/dy2)) +
     >          rrr(irk)*(phip(1,1)-dt2re*pr(1,1))
        pr(1,1)=prp(1,1) + pr(1,1)/aaa(irk)
        pr(n1,1)=ggg(irk)*(phi(n1,1)-dt2re*
     >           ((phi(n1m1,1)-phi(n1,1))/dx2 +
     >            (phi(n1,2)-phi(n1,1))/dy2)) +
     >           rrr(irk)*(phip(n1,1)-dt2re*pr(n1,1))
        pr(n1,1)=prp(n1,1) + pr(n1,1)/aaa(irk)
        pr(1,n2)=ggg(irk)*(phi(1,n2)-dt2re*
     >           ((phi(2,n2)-phi(1,n2))/dx2 +
     >            (phi(1,n2m1)-phi(1,n2))/dy2)) +
     >           rrr(irk)*(phip(1,n2)-dt2re*pr(1,n2))
        pr(1,n2)=prp(1,n2) + pr(1,n2)/aaa(irk)
        pr(n1,n2)=ggg(irk)*(phi(n1,n2)-dt2re*
     >            ((phi(n1m1,n2)-phi(n1,n2))/dx2 +
     >             (phi(n1,n2m1)-phi(n1,n2))/dy2)) +
     >            rrr(irk)*(phip(n1,n2)-dt2re*pr(n1,n2))
        pr(n1,n2)=prp(n1,n2) + pr(n1,n2)/aaa(irk)
*
      return
      end     
*      
*                
*********************************************************            
*
      subroutine output
      include 'param.inc'
      common/one/qx(n1p1,n2p2),qy(n1p2,n2p1),dx,dy,dx2,dy2,
     >           dt,re,aaa(3),ggg(3),rrr(3)
      common/two/hx(2:n1,2:n2p1),hy(2:n1p1,2:n2),
     >           hxp(2:n1,2:n2p1),hyp(2:n1p1,2:n2)
      common/three/ax(2:n1,2:n2p1),ay(2:n1p1,2:n2)
      common/six/pr(n1,n2),phi(n1,n2),phip(n1,n2),prp(n1,n2)
      common/seven/x(n1p1),y(n2p1),xc(n1),yc(n2)
      dimension uvel(n1p1,n2p1),vvel(n1p1,n2p1),strmfn(n1p1,n2p1),
     >          vortfn(n1p1,n2p1),temp(n1p1)
*

      do j=1,n2p1
          u9 = 0.5*(qx(n1/2+1,j)+qx(n1/2+1,j+1))
          v9 = 0.5*(qy(n1/2+1,j)+qy(n1/2+2,j))
C velocity at the vertical line
        write(10,102) y(j),u9,v9
      enddo

      do i=1,n1p1
          u9 = 0.5*(qx(i,n2/2+1)+qx(i,n2/2+2))
          v9 = 0.5*(qy(i,n2/2+1)+qy(i+1,n2/2+1))
C velocity at the vertical line
        write(11,102) x(i),u9,v9
      enddo

      write(10,101)
      write(11,101)
101   format(1x)
102   format(2x,3(1pe15.4))
*
* save values of velocity at mesh vertices...
*
      do i=1,n1p1
        uvel(i,1) = 0.5*(qx(i,1)+qx(i,2))
        uvel(i,n2p1) = 0.5*(qx(i,n2p2)+qx(i,n2p1))
        do j=2,n2
          uvel(i,j) = 0.5*(qx(i,j+1)+qx(i,j))
        enddo
      enddo
      do j=1,n2p1
        vvel(1,j) = 0.5*(qy(1,j)+qy(2,j))
        vvel(n1p1,j) = 0.5*(qy(n1p2,j)+qy(n1p1,j))
        do i=2,n1
          vvel(i,j) = 0.5*(qy(i+1,j)+qy(i,j))
        enddo
      enddo
*
* calculate streamfunction and vorticity
*
      do i=1,n1p1
        temp1 = 0.
        strmfn(i,1) = 0.
        do j=2,n2p1
          temp1 = temp1 + qx(i,j)*dy
          strmfn(i,j) = temp1
        enddo
      enddo
*
      do j=1,n2p1
	do i=1,n1p1
	  vortfn(i,j)=(qy(i+1,j)-qy(i,j))/dx -
     >                (qx(i,j+1)-qx(i,j))/dy
        enddo
      enddo
*
c      open(11,file='output.dat',form='formatted',
c     >      status='unknown')
c* this is Tecplot format...
c      write(11,'(''TITLE = "2D Cavity"'')')
c      write(11,'(''VARIABLES = X, Y, U, V, VORT, STRM'')')
c*23456789012345678901234567890123456789012345678901234567890123456789012
c      write(11,'(''ZONE T="Zone-one", I='',i4,''  J='',i4,''  F=BLOCK'')
c     >    ') n1p1,n2p1
c      do j=1,n2p1
c        write(11,905) (x(i),i=1,n1p1)
c      enddo
c      do j=1,n2p1
c        do i=1,n1p1
c          temp(i) = y(j)
c        enddo
c        write(11,905) (temp(i),i=1,n1p1)
c      enddo
c      write(11,905) ((uvel(i,j),i=1,n1p1),j=1,n2p1)
c      write(11,905) ((vvel(i,j),i=1,n1p1),j=1,n2p1)
c      write(11,905) ((vortfn(i,j),i=1,n1p1),j=1,n2p1)
c      write(11,905) ((strmfn(i,j),i=1,n1p1),j=1,n2p1)
      write(19,908) ((vortfn(i,j),i=1,n1p1),j=1,n2p1)
      write(18,907) ((strmfn(i,j),i=1,n1p1),j=1,n2p1)
*
* write in format for Gnuplot
      open(12,file='vortfn.dat',form='formatted',status='unknown')
      do j=1,n2p1
	do i=1,n1p1
	  write(12,*) x(i),y(j),vortfn(i,j)
        enddo
        write(12,906) 
      enddo
*
      open(13,file='strmfn.dat',form='formatted',status='unknown')
      do j=1,n2p1
	do i=1,n1p1
c  write(13,*) x(i),y(j),strmfn(i,j)
C Rescale to match the LBM code
        write(13,*) strmfn(i,j)*96.0/10.0
        enddo
        write(13,906) 
      enddo
*
      open(14,file='uvel.dat',form='formatted',status='unknown')
      do j=1,n2p1
	do i=1,n1p1
	  write(14,*) x(i),y(j),uvel(i,j)
        enddo
        write(14,906) 
      enddo
*
      open(15,file='vvel.dat',form='formatted',status='unknown')
      do j=1,n2p1
	do i=1,n1p1
	  write(15,*) x(i),y(j),vvel(i,j)
        enddo
        write(15,906) 
      enddo
*
* added 10/8/96
      open(16,file='pres.dat',form='formatted',status='unknown')
      do j=1,n2
	do i=1,n1
	  write(16,*) xc(i),yc(j),pr(i,j)
        enddo
        write(16,906) 
      enddo
*
  905 format(6(e12.4))
  906 format(x)
  907 format(2x,8f8.4)
  908 format(2x,8f8.2)
*
      return
      end      
*                
**********************************************************************
* Subroutine SY, for solving a tridiagonal system of equations.  Based
* on the Thomas algorithm, from "Computational Fluid Mechanics and
* Heat Transfer," by Anderson, Tannehill, and Pletcher.
* 27 May 1992   Jim DeSpirito
*
      subroutine sy(il,iu,bb,dd,aa,cc)
*
*      double precision aa(iu),bb(iu),cc(iu),dd(iu)
      real aa(iu),bb(iu),cc(iu),dd(iu)
*
*.....il = subscript of first equation
*.....iu = subscript of last equation
*.....bb = coefficient behind diagonal
*.....dd = coefficient on diagonal
*.....aa = coefficient ahead of diagonal
*.....cc = element of constant vector
*
*.....establish upper triangular matrix
*
      lp = il + 1
      do 10 i = lp,iu
	 r = bb(i)/dd(i-1)
	 dd(i) = dd(i) - r * aa(i-1)
	 cc(i) = cc(i) - r * cc(i-1)
   10 continue
*
*.....back substitution
*
      cc(iu) = cc(iu) / dd(iu)
      do 20 i = lp,iu
	 j = iu - i + il
	 cc(j) = (cc(j) - aa(j) * cc(j+1)) / dd(j)
   20 continue
*
*.....solution stored in cc
*
      return
      end
*
*********************************************************            
*
      subroutine trigon(irk)
      include 'param.inc'
        PARAMETER ( n1h=n1/2 , n2h=n2/2 )
* use trigonometric series method to calculate the poisson equation 
      common/one/qx(n1p1,n2p2),qy(n1p2,n2p1),dx,dy,dx2,dy2,
     >           dt,re,aaa(3),ggg(3),rrr(3)
      common/six/pr(n1,n2),phi(n1,n2),phip(n1,n2),prp(n1,n2)
      dimension a(n2),b(n2),c(n2),d(n2),qrhs(n1,n2)
        INTEGER IFAXX(13)
        REAL TRIGSX(3*n1h+1)
        COMMON /TRIG/  TRIGSX
        COMMON /FAX/   IFAXX
*
* set RHS array containing known values
c first calc Laplacian of old phi & store in pr() for use in press sub...
c save old pressure value first (for pressure calc only)
c
        do j=1,n2
          do i=1,n1
            prp(i,j) = pr(i,j)
          enddo
        enddo
c
      do j=2,n2m1
        do i=2,n1m1
          pr(i,j)=((phi(i-1,j)-2.*phi(i,j)+phi(i+1,j))/dx2 +
     >             (phi(i,j-1)-2.*phi(i,j)+phi(i,j+1))/dy2)
        enddo
      enddo
      do j=2,n2m1
        pr(1,j)=((phi(2,j)-phi(1,j))/dx2 +
     >           (phi(1,j-1)-2.*phi(1,j)+phi(1,j+1))/dy2)
        pr(n1,j)=((phi(n1m1,j)-phi(n1,j))/dx2 +
     >            (phi(n1,j-1)-2.*phi(n1,j)+phi(n1,j+1))/dy2)
      enddo
      do i=2,n1m1
        pr(i,1)=((phi(i-1,1)-2.*phi(i,1)+phi(i+1,1))/dx2 +
     >           (phi(i,2)-phi(i,1))/dy2)
        pr(i,n2)=((phi(i-1,n2)-2.*phi(i,n2)+phi(i+1,n2))/dx2 +
     >            (phi(i,n2m1)-phi(i,n2))/dy2)
      enddo
      pr(1,1)=(phi(2,1)-phi(1,1))/dx2 + (phi(1,2)-phi(1,1))/dy2
      pr(n1,1)=(phi(n1m1,1)-phi(n1,1))/dx2+(phi(n1,2)-phi(n1,1))/dy2
      pr(1,n2)=(phi(2,n2)-phi(1,n2))/dx2+(phi(1,n2m1)-phi(1,n2))/dy2
      pr(n1,n2)=(phi(n1m1,n2)-phi(n1,n2))/dx2 + 
     >          (phi(n1,n2m1)-phi(n1,n2))/dy2
c
      do j=1,n2
        do i=1,n1
          phip(i,j) = phi(i,j)
          qrhs(i,j) = ((qx(i+1,j+1)-qx(i,j+1))/dx +
     >                (qy(i+1,j+1)-qy(i+1,j))/dy)/(ggg(irk)*dt) -
     >                 rrr(irk)*pr(i,j)/ggg(irk)
          qrhs(i,j) = dx2*qrhs(i,j)
        enddo
      enddo
*
* transform the rhs in x-direction
      pi=6.*asin(0.5)
      CALL R2FC1(qrhs)

      phi(1,1)=0.0
      phi(1,2)=qrhs(1,1)/(dx2/dy2)
      do j=3,n2
      phi(1,j)=-phi(1,j-2)+2.0*phi(1,j-1)
     1  +qrhs(1,j-1)/(dx2/dy2)
      enddo

      call sy(1,n2,a,b,c,d)

      do i=2,n1

      waven=2.*(1.0-cos(pi*(i-1)/n1))
      do j=1,n2
      a(j)=dx2/dy2
      b(j)=-2.*dx2/dy2-waven
      c(j)=a(j)
      d(j)=qrhs(i,j)
      enddo
      b(1)=b(1)+c(1)
      b(n2)=b(n2)+c(1)
      call sy(1,n2,a,b,c,d)
      do j=1,n2
      phi(i,j)=d(j)
      enddo

      enddo
*
* transform the solution back to physical space
      CALL F2RC1(phi)
*
      return
      end
*
**************************************************************
