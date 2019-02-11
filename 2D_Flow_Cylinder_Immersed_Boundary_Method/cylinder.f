       program incomp2d
C
C
C      Immersed boundary method using the forcing method of Goldstein et al. (1993)
C
c the following is the contents of the include file:      
c      parameter (n1=32,n1p2=n1+2,n1p1=n1+1,n1m1=n1-1,
c     >           n2=32,n2p2=n2+2,n2p1=n2+1,n2m1=n2-1)
      include 'param.inc'
        PARAMETER ( n1h=n1/2 , n2h=n2/2 )
      common/one/qx(n1p1,n2p2),qy(n1p2,n2p1),dx,dy,dx2,dy2,
     >           dt,re,aaa(3),ggg(3),rrr(3)
      common/two/hx(2:n1,2:n2p1),hy(2:n1p1,2:n2),
     >           hxp(2:n1,2:n2p1),hyp(2:n1p1,2:n2)
      common/three/ax(2:n1,2:n2p1),ay(2:n1p1,2:n2)
      common/six/pr(n1,n2),phi(n1,n2),phip(n1,n2),prp(n1,n2)
      common/seven/x(n1p1),y(n2p1),xc(n1),yc(n2)
      common/force/qxf(n1p1,n2p2),qyf(n1p2,n2p1)
      dimension qxold(n1p1,n2p2)
      dimension xbub(1500),ybub(1500),qxs(1500),qys(1500),
     >          ffx(1500),ffy(1500),fsumx(1500),fsumy(1500),
     >          spr(1500)
      dimension ijqx(n1p1,n2p1),ijqy(n1p1,n2p1)
      dimension ddx(2,2),ddy(2,2),tmp(n1,n2)
      dimension uvel(n1p1,n2p1),vvel(n1p1,n2p1)
      dimension qxf0(n1p1,n2p1),qyf0(n1p1,n2p1)
        INTEGER IFAXX(13)
        REAL TRIGSX(3*n1h+1)
        COMMON /TRIG/  TRIGSX
        COMMON /FAX/   IFAXX
      real lenx,leny
      logical ipress
      NAMELIST/INPUT/niter,nshort,nlong,iread,re,lenx,leny,dt,
     >              rerror,ipress,iflow,alphaf,betaf,xcent,ycent,
     >               dbub,ddd,pout,uc
* FFT setup
      CALL FFTFAX(n1,IFAXX,TRIGSX)
      open(5,file='data.inp',form='formatted',status='old')
      open(7,file='data.out',form='formatted',status='unknown')
      open(501,file='restart.flo',form='unformatted',
     >          status='unknown')
      open(502,file='endrun.flo',form='unformatted',
     >          status='unknown')
* DX(),DY(),arrays are not needed for dx,dy uniform
* here, x,y are locations of cell boundarys (this will change 
* for dx,dy not uniform
*
      zero = 1.0e-30
      pi = acos(-1.0)
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
      ndump=100
* Lenx,Leny = length,width of cavity
* input Reynolds number
      read(5,INPUT)
      write(7,INPUT)
c if IREAD = 0, new start (else RESTART)
      IF (iread .eq. 0) THEN
      itbeg = -1
      tau = 0.0
      dx = lenx / float(n1)
      dy = leny / float(n2)
      dx2 = dx ** 2.0
      dy2 = dy ** 2.0
* fill x,y vectors (grid line coords)
      do i=1,n1p1
        x(i) = float(i)*dx - dx
      enddo
      do j=1,n2p1
        y(j) = float(j)*dy - dy
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
c
c===> Forcing set-up
c
ccc this was cylinder
c read in xcent,ycent from input
      rad = 0.5
c      rad = 0.375
      delth = 2*pi/1500
      theta = 0.
      do i=1,1500
cc        xbub(i) = xcent + rad * cos(theta)
cc        ybub(i) = ycent + rad * sin(theta)
        xbub(i) = xcent - rad * sin(theta)
        ybub(i) = ycent - rad * cos(theta)
        theta = theta + delth
        write(11,*) i,xbub(i),ybub(i)
        fsumx(i) = 0.
        fsumy(i) = 0.
      enddo
      ij = 1500
c  The trailing edge
      npr = ij/2 + 1
c
c this is the bubble 
c
      goto 1111
c read in xcent,ycent,dbub,ddd from input
c dia, D, = 1.0
      coef = 15.0
      xbub(1) = xcent
      ybub(1) = ycent
      ij=1
12    y1 = ybub(ij) + ddd
       do ijk=1,5
        temp = 1. - 1./sqrt(1.+coef*(y1-ycent))
        temp = temp*dbub/2.0
        x1 = xcent - temp
        dis = sqrt((x1-xbub(ij))**2 + (y1-ybub(ij))**2)
        y1=ddd/dis*(y1-ybub(ij))+ybub(ij)
       enddo
       ybub(ij+1)=y1
       temp = 1. - 1./sqrt(1.+coef*(y1-ycent))
       temp = temp*dbub/2.0
        xbub(ij+1) = xcent - temp
        ij=ij+1
        if( (dbub-(ybub(ij)-ycent)).lt. ddd)goto 13
        goto 12
13      continue
        xbub(ij+1)=xcent-0.75*dbub/2.
        ybub(ij+1)=ycent+dbub
        kkk=ij
        ij =ij +1
        mmm=nint(0.75/ddd)
        dxdx=0.75/mmm
        do i=1,mmm
        ij=ij+1
        ybub(ij)=ycent+dbub
        xbub(ij)=xbub(ij-1)+dxdx
        enddo

        do i=1,kkk
        ij=ij+1
        xbub(ij)=1.0 - xbub(kkk+1-i)
        ybub(ij)=ybub(kkk+1-i)
        enddo
C
c save number of points along side for pressure calcs
        npr = kkk+1
c set integral variables to zero                 
        do i=1,ij
          fsumx(i) = 0.
          fsumy(i) = 0.
          spr(i) = 0.
        enddo
C
 1111 continue      
       write(*,*)'the geometry is finished',ij
      do i=1,ij
        write(99,909) xbub(i),ybub(i)
      enddo
  909 format(2(1pe12.4))
c
c===> end of forcing set-up      
c      
*
* stability criterion (Peyret & Taylor, 1983)...
* there are two: eqs. 6.2.10 & 6.2.11
* for eq. 6.2.10, assume |u0| & |v0| are equal to 1.0
*
c      c1 = dt * re
c      write(*,*) 'Stability criterion 1 = ',c1
c      if(c1 .gt. 1.0) then
c        write(*,*) 'Stability criterion #1 greater than 1.0'
c        write(*,*) 'Program stopped'
c        stop
c      endif
c      c2 = 4*dt/(re*dx2)
c      write(*,*) 'Stability criterion 2 = ',c2              
c      if(c2 .gt. 1.0) then
c        write(*,*) 'Stability criterion #2 greater than 1.0'
c        write(*,*) 'Program stopped'
c        stop
c      endif
c
c calculate CFL number based on mean flow of 1.0
c
      CFL = dt/min(dx,dy)
      write(*,*) 'CFL number = ',CFL
c
*
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
          pr(i,j) = pout
          phi(i,j) = zero
        enddo
      enddo
c
c free stream velocity at inlet 
c zero gradient at exit
      do i=2,n1p1
        qy(i,1) = uc
        qy(i,n2p1) = qy(i,n2)
      enddo
*
* fill bc for velocities defined outside the domain...
c zero vel at inlet, zero gradient at exit
      do i=1,n1p1
        qx(i,1) = -qx(i,2)
        qx(i,n2p2) = qx(i,n2p1)
      enddo
c far-field at i=1 & i=imax
      do j=1,n2p1
        qy(1,j) = qy(2,j)
        qy(n1p2,j) = qy(n1p1,j)
      enddo
c
c fill bc for velocities defined along the boundary within
c the domain
c zero vel at far-field (i=1) & (i=imax)
      do j=2,n2p1
        qx(1,j) = zero
        qx(n1p1,j) = zero
      enddo
c
c Call poisson equation to solve for initial pressure field
c
      call hxycalc
      call initpr
c
c
      ELSE
      read(501) itbeg,ndump,epsmax
      read(501) tau,dx,dy,dx2,dy2
      read(501) x,y,xc,yc
      read(501) qx,qy,hxp,hyp
      read(501) pr,phi
      read(501) xbub,ybub,fsumx,fsumy,npr,ij,spr
      ENDIF
c
c Here open files to append if restart
      if(tau .gt. 0.)then
        open(14,file='v139_154.dat',form='formatted',
     >           status='unknown',access='append')
        open(40,file='eps.dat',form='formatted',
     >           status='unknown',access='append')
        open(12,file='norm.dat',form='formatted',
     >           status='unknown',access='append')
        open(15,file='bound.dat',form='formatted',
     >           status='unknown',access='append')
      else
      open(14,file='v139_154.dat',form='formatted',status='unknown')
        open(40,file='eps.dat',form='formatted',status='unknown')
        open(12,file='norm.dat',form='formatted',status='unknown')
        open(15,file='bound.dat',form='formatted',status='unknown')
      endif
c
* here is the time-stepping loop...
*
      itbeg = itbeg+1
      niter = niter+itbeg-1
      do itime=itbeg,niter
*      
      write(14,*)tau,qy(139,154)

        if(mod(itime,nshort).eq.0.) then
c
c     compute drag coefficient
c 
      amax= ( qy(1,n2p1)+qy(n1p2,n2p1) ) /2.
      cdcd=0.0
      do i=1,n1p2
      cdcd=cdcd + 2.*qy(i,n2p1)/amax *(1.0-qy(i,n2p1)/amax)*dx
      enddo
c
c ---> calculate the pressure at the boundary
        DO K=1,npr
          ix = int(xbub(k)/dx + 0.5)
          iy = int(ybub(k)/dy + 0.5)
          ixp = ix + 1
          iyp = iy + 1
          xs = xbub(k) 
          ys = ybub(k)
          aa = (xs-x(ix)-0.5*dx)/dx
          bb = (ys-y(iy)-0.5*dy)/dy
          spr(k) = pr(ix,iy)*(1.-aa)*(1.-bb) + pr(ixp,iy)*aa*(1.-bb)
     1          + pr(ix,iyp)*(1.-aa)*bb + pr(ixp,iyp)*aa*bb
cccccc for check         
         if(mod(itime,nlong).eq.0) then
          if(k.eq.1)then
          write(*,*)ix,iy,pr(ix,iy),pr(ix,iyp),pr(ixp,iy),pr(ixp,iyp)
          endif
         endif
cccccc for check         
        ENDDO
         delp = spr(1) - spr(npr)
         p1=0.
         p2=0.
         p3=0.
         do i=1,32
         p1=p1+(pr(i,32)+pr(i,33))/2.
         p3=p3+(pr(i,64)+pr(i,65))/2.
         enddo
         p1=p1/32.
         p3=p3/32.
         do i=1,4
         p2=p2+(pr(i,64)+pr(i,65)+pr(33-i,64)+pr(33-i,65))/4.
         enddo
         p2=p2/4.0
         p12=p1-p2
         p13=p1-p3
cccccc for check         
         pp11=0.25*(pr(16,32)+pr(16,33)+pr(17,32)+pr(17,33))
        if(mod(itime,nlong).eq.0)then
         write(*,*)'pp11=',pp11,' pp1, pp2,pp3,pp4 = ',
     >      pr(16,32),pr(16,33),pr(17,32),pr(17,33)    
         endif
cccccc for check         
          write(40,777) itime,tau,epsmax,cdcd
          write(*,777) itime,tau,epsmax,cdcd
777       format(2x,i5,3(e14.4))
          if(mod(itime,nlong).eq.0)then
            call output(ndump)
c calculate the divergence field
            divmax=0
            do j=1,n2
              do i=1,n1
                div=(qx(i+1,j)-qx(i,j))/dx
     >              + (qy(i,j+1)-qy(i,j))/dy
                if(abs(div).gt.divmax)divmax=abs(div)
                tmp(i,j) = div
              enddo
            enddo
            write(*,*)'divmax= ',divmax
          endif
        endif
c
c write out divergence field for checking only
      if(itime .eq. niter)then
        write(600,900) ((tmp(i,j),i=1,n1),j=1,n2)
      endif
  900 format(8(1pe12.4))
c
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
c
        do j=1,n2p1
          do i=1,n1p1
            qxf0(i,j) = 0.
            qyf0(i,j) = 0.
            ijqx(i,j) = 0.
            ijqy(i,j) = 0.
          enddo
        enddo
c
c interpolate the field velocity (assumes uniform mesh)
c
c save values of velocity at mesh vertices...
C
C Treat forcing at the immersed boundary
      do j=1,n2p1
        do i=1,n1p1
         uvel(i,j) = 0.5*(qx(i,j+1)+qx(i,j))
         vvel(i,j) = 0.5*(qy(i+1,j)+qy(i,j))
        enddo
      enddo

        ub2 = 0.
        DO K=1,ij
cc for check
          ix = int(xbub(k)/dx)+1
          iy = int(ybub(k)/dy)+1
          ixp = ix + 1
          iyp = iy + 1
          xs = xbub(k)
          ys = ybub(k)
          aa = (xs-x(ix))/dx
          bb = (ys-y(iy))/dy
          qxs(k) = uvel(ix,iy)*(1.-aa)*(1.-bb) + uvel(ixp,iy)*aa*(1.-bb)
     1          + uvel(ix,iyp)*(1.-aa)*bb + uvel(ixp,iyp)*aa*bb
          qys(k) = vvel(ix,iy)*(1.-aa)*(1.-bb) + vvel(ixp,iy)*aa*(1.-bb)
     1          + vvel(ix,iyp)*(1.-aa)*bb + vvel(ixp,iyp)*aa*bb
C compute the forcing magnitude
          fsumx(k) = fsumx(k) + qxs(k)*aaa(irk)*dt
          fsumy(k) = fsumy(k) + qys(k)*aaa(irk)*dt
          ffx(k) = alphaf*fsumx(k) + betaf*qxs(k)
          ffy(k) = alphaf*fsumy(k) + betaf*qys(k)
c distribute back
          qxf0(ix,iy) = qxf0(ix,iy) + (1.-aa)*(1.-bb)*ffx(k)
          qxf0(ix,iyp) = qxf0(ix,iyp) + (1.-aa)*bb*ffx(k)
          qxf0(ixp,iy) = qxf0(ixp,iy) + aa*(1.-bb)*ffx(k)
          qxf0(ixp,iyp) = qxf0(ixp,iyp) + aa*bb*ffx(k)
          qyf0(ix,iy) = qyf0(ix,iy) + (1.-aa)*(1.-bb)*ffy(k)
          qyf0(ix,iyp) = qyf0(ix,iyp) + (1.-aa)*bb*ffy(k)
          qyf0(ixp,iy) = qyf0(ixp,iy) + aa*(1.-bb)*ffy(k)
          qyf0(ixp,iyp) = qyf0(ixp,iyp) + aa*bb*ffy(k)
          ijqx(ix,iy) = ijqx(ix,iy)+1
          ijqx(ixp,iy) = ijqx(ixp,iy)+1
          ijqx(ix,iyp) = ijqx(ix,iyp)+1
          ijqx(ixp,iyp) = ijqx(ixp,iyp)+1
          ijqy(ix,iy) = ijqy(ix,iy)+1
          ijqy(ixp,iy) = ijqy(ixp,iy)+1
          ijqy(ix,iyp) = ijqy(ix,iyp)+1
          ijqy(ixp,iyp) = ijqy(ixp,iyp)+1
c
c calculate the L2 norm for parameter check
          temp = qxs(k)**2 + qys(k)**2
          ub2 = ub2 + temp
        ENDDO
c
c calculate the L2 norm for parameter check
        anorm = sqrt(ub2/float(ij))
        if(mod(itime,5).eq.0)then
          write(12,901) itime,tau,anorm
        endif
  901 format(x,i5,2(1pe12.4))
c
c --> finish calculating the force
        do j=1,n2p1
          do i=1,n1p1
            if(ijqx(i,j).gt.0) then
              qxf0(i,j)=qxf0(i,j)/ijqx(i,j)
            endif
          enddo
        enddo
        do j=1,n2p1
          do i=1,n1p1
            if(ijqy(i,j).gt.0) then
            qyf0(i,j)=qyf0(i,j)/ijqy(i,j)
            endif
          enddo
        enddo
c
        do j=2,n2p1
          do i=1,n1p1
              qxf(i,j)=(qxf0(i,j)+qxf0(i,j-1))/2.0
          enddo
        enddo
c
        do i=1,n1p1
         qxf(i,1)=0.0
         qxf(i,n2p2)=0.0
        enddo
c
        do j=1,n2p1
          do i=2,n1p1
              qyf(i,j)=(qyf0(i,j)+qyf0(i-1,j))/2.0
          enddo
        enddo
c
        do j=1,n2p1
         qyf(1,j)=0.0
         qyf(n1p2,j)=0.0
        enddo
c        
c
c ===> end of forcing setup
c
*
        call hxycalc
        call axycalc
c
* calculate intermediate velocities..
* (sub-step 2)
        call qhat(irk,uc)
*
* calculate the phi from the poisson equation...
* (sub-step 3)
        call trigon(irk)
* 
* calculate the final velocities (at current t-step)...
* (sub-step 4)
        call ufinal(irk,uc)
*
* update the pressure (assume wanted)...
* (sub-step 5)
        call press(irk)
*
c end of Runge-Kutta steps...
c
      enddo
c
        if(mod(itime,nlong).eq.0) then
          do k=1,npr
            write(15,907) k,xbub(k),ybub(k),spr(k)
          enddo
          write(15,906)
        endif
c
  902 format(3(1pe12.4))
  906 format(1x)
  907 format(i6,3(1pe12.4))
* determine if steady-state has been reached...
*
      if(mod(itime,nshort).eq.0.)then
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
      endif
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
c write output for last time step
cccc            call output(ndump)
c
c write restart file:
c
      write(502) niter,ndump,epsmax
      write(502) tau,dx,dy,dx2,dy2
      write(502) x,y,xc,yc
      write(502) qx,qy,hxp,hyp
      write(502) pr,phi
      write(502) xbub,ybub,fsumx,fsumy,npr,ij,spr
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
      subroutine qhat(irk,uc)
      include 'param.inc'
* calculate the intermediate-time velocities
      common/one/qx(n1p1,n2p2),qy(n1p2,n2p1),dx,dy,dx2,dy2,
     >           dt,re,aaa(3),ggg(3),rrr(3)
      common/two/hx(2:n1,2:n2p1),hy(2:n1p1,2:n2),
     >           hxp(2:n1,2:n2p1),hyp(2:n1p1,2:n2)
      common/three/ax(2:n1,2:n2p1),ay(2:n1p1,2:n2)
      common/six/pr(n1,n2),phi(n1,n2),phip(n1,n2),prp(n1,n2)
      common/force/qxf(n1p1,n2p2),qyf(n1p2,n2p1)
*note: max dim. on a,b,c,d is max(n1p1,n2p1)...    
      dimension a(max(n1p1,n2p1)),b(max(n1p1,n2p1)),
     >          c(max(n1p1,n2p1)),d(max(n1p1,n2p1))
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
c
c===> forcing
c
        do j=2,n2p1
          do i=2,n1
            ax(i,j) = ax(i,j) + qxf(i,j)
          enddo
        enddo
c
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
        b(2)=b(2)-a(2)
c zero-gradient
        b(n2p1)=b(n2p1)+c(n2p1)
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
c 
c===> forcing
c
        do j=2,n2
          do i=2,n1p1
            ay(i,j) = ay(i,j) + qyf(i,j)
          enddo
        enddo
c        
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
c for far-field and symmetry
        b(2)=b(2)+c(2)
        b(n1p1)=b(n1p1)+c(n1p1)
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
c for zero-gradient        
        b(n2)=b(n2)+c(n2)
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
c
c fill bc for velocities defined outside the domain...
c first exit BC
      do i=2,n1p1
        qy(i,n2p1) = qy(i,n2)
      enddo
c zero vel at inlet, zero gradient at exit
      do i=1,n1p1
        qx(i,1) = -qx(i,2)
        qx(i,n2p2) =  qx(i,n2p1)
      enddo
c symmetry at i=imax & i=1
      do j=1,n2p1
        qy(1,j) = qy(2,j)
        qy(n1p2,j) = qy(n1p1,j)
      enddo
c
      return
      end
*                
*********************************************************            
*
      subroutine ufinal(irk,uc)
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
c
c fill bc for velocities defined outside the domain...
c first exit BC
      do i=2,n1p1
        qy(i,n2p1) = qy(i,n2)
      enddo
c zero vel at inlet, zero gradient at exit
      do i=1,n1p1
        qx(i,1) = -qx(i,2)
        qx(i,n2p2) =  qx(i,n2p1)
      enddo
c far-field at i=1 & i=imax
      do j=1,n2p1
        qy(1,j) = qy(2,j)
        qy(n1p2,j) = qy(n1p1,j)
      enddo
c
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
      subroutine output(ndump)
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
c this is v-vel profile
      nh = n1p2/2
      open(10,file='velprof.dat',form='formatted',
     >     status='unknown')
      do j=1,n2p1
        write(10,909) j,y(j),vvel(nh,j)
      enddo
      write(10,906)
  909 format(x,i3,2e12.4)
*
* calculate streamfunction and vorticity
*
c calc in j-dir
c      do i=1,n1p1
c        temp1 = 0.
c        strmfn(i,1) = 0.
c        do j=2,n2p1
c          temp1 = temp1 + qx(i,j)*dy
c          strmfn(i,j) = temp1
c        enddo
c      enddo
*
c calc in i-dir
      do j=1,n2p1
        temp1 = 0.
        strmfn(1,j) = 0.
        do i=2,n1p1
          temp1 = temp1 + qy(i,j)*dx
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
      write(ndump,901)((vortfn(i,k),i=1,n1p1),k=1,n2p1)
      write(ndump+100,901)((strmfn(i,k),i=1,n1p1),k=1,n2p1)
      write(ndump+200,901) ((uvel(i,j),i=1,n1p1),j=1,n2p1)
      write(ndump+200,901) ((vvel(i,j),i=1,n1p1),j=1,n2p1)
      write(ndump+300,901) ((pr(i,j),i=1,n1),j=1,n2)
 901   format(8(1pe12.4))
      close(ndump-10)
      close(ndump)
      close(ndump+10)
      close(ndump+20)
906    format(1x)
c
      ndump = ndump+1      
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
      dimension a(max(n1,n2)),b(max(n1,n2)),c(max(n1,n2)),
     >          d(max(n1,n2)),qrhs(n1,n2)
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

c      call sy(1,n2,a,b,c,d)

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
*********************************************************            
*
      subroutine initpr
      include 'param.inc'
        PARAMETER ( n1h=n1/2 , n2h=n2/2 )
* use trigonometric series method to calculate the poisson equation 
      common/one/qx(n1p1,n2p2),qy(n1p2,n2p1),dx,dy,dx2,dy2,
     >           dt,re,aaa(3),ggg(3),rrr(3)
      common/two/hx(2:n1,2:n2p1),hy(2:n1p1,2:n2),
     >           hxp(2:n1,2:n2p1),hyp(2:n1p1,2:n2)
      common/six/pr(n1,n2),phi(n1,n2),phip(n1,n2),prp(n1,n2)
      dimension a(max(n1,n2)),b(max(n1,n2)),c(max(n1,n2)),
     >          d(max(n1,n2)),qrhs(n1,n2)
        INTEGER IFAXX(13)
        REAL TRIGSX(3*n1h+1)
        COMMON /TRIG/  TRIGSX
        COMMON /FAX/   IFAXX
*
* set RHS array containing known values
c
      do j=2,n2-1
        do i=2,n1-1
          phi(i,j) = (hx(i+1,j)-hx(i,j))/dx +
     >               (hy(i,j+1)-hy(i,j))/dy
        enddo
      enddo
      do j=2,n2-1
        phi(1,j) = hx(2,j+1)/dx + (hy(2,j+1)-hy(2,j))/dy
        phi(n1,j) = -hx(n1,j+1)/dx + (hy(n1p1,j+1)-hy(n1p1,j))/dy
      enddo
      do i=2,n1-1
        phi(i,1) = (hx(i+1,2)-hx(i,2))/dx + hy(i+1,2)/dy
        phi(i,n2) = (hx(i+1,n2p1)-hx(i,n2p1))/dx - hy(i+1,n2)/dy
      enddo
      phi(1,1) = hx(2,2)/dx + hy(2,2)/dy
      phi(n1,n2) = -hx(n1,n2p1)/dx - hy(n1p1,n2)/dy
      phi(1,n2) = hx(2,n2p1)/dx - hy(2,n2)/dy
      phi(n1,1) = -hx(n1,2)/dx + hy(n1p1,2)
c
      do j=1,n2
        do i=1,n1
          phi(i,j) = dx2*phi(i,j)
        enddo
      enddo
*
* transform the rhs in x-direction
      pi=6.*asin(0.5)
      CALL R2FC1(phi)

      pr(1,1)=0.0
      pr(1,2)=phi(1,1)/(dx2/dy2)
      do j=3,n2
      pr(1,j)=-pr(1,j-2)+2.0*pr(1,j-1)
     1  + phi(1,j-1)/(dx2/dy2)
      enddo

ccc      call sy(1,n2,a,b,c,d)

      do i=2,n1

      waven=2.*(1.0-cos(pi*(i-1)/n1))
      do j=1,n2
      a(j)=dx2/dy2
      b(j)=-2.*dx2/dy2-waven
      c(j)=a(j)
      d(j)=phi(i,j)
      enddo
      b(1)=b(1)+c(1)
      b(n2)=b(n2)+c(1)
      call sy(1,n2,a,b,c,d)
      do j=1,n2
      pr(i,j)=d(j)
      enddo

      enddo
*
* transform the solution back to physical space
      CALL F2RC1(pr)
*
      return
      end
*
**************************************************************
