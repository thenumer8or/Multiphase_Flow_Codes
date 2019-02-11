!RANDOM
        Function ranpp(idum,iy1,iv1)
!
!       Minimal random number generator of Park and Miller with
!       Bays-Durham shuffle and added safegards, see Numerical Recipes
!
        Integer idum, IA,IM,IQ,IR,NTAB,NDIV
        Real ranpp,AM,EPS,RNMX
        Parameter (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,	&
             NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
        INTEGER j,k,iv1(NTAB),iy1
!       Save iv1,iy1
!       Data iv1/NTAB*0/,iy1 /0/

        if (idum.le.0 .or. iy1.eq.0) then
                idum=max(-idum,1)
                do j=NTAB+8,1,-1
                   k=idum/IQ
                   idum=IA*(idum-k*IQ)-IR*k
                   if (idum.lt.0) idum=idum+IM
                   if (j.le.NTAB) iv1(j)=idum
                enddo
                iy1=iv1(1)
        endif
        k=idum/IQ
        idum=IA*(idum-k*IQ)-IR*k
        if(idum.lt.0) idum=idum+IM
        j=1+iy1/NDIV
        iy1=iv1(j)
        iv1(j)=idum
        ranpp=min(AM*iy1,RNMX)

        return
        END



        Function ranff(idum,iy,iv)
!
!       Minimal random number generator of Park and Miller with
!       Bays-Durham shuffle and added safegards, see Numerical Recipes
!
        Integer idum, IA,IM,IQ,IR,NTAB,NDIV
        Real ranff,AM,EPS,RNMX
        Parameter (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,	&
             NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
        INTEGER j,k,iv(NTAB),iy
!       Save iv,iy
!       Data iv/NTAB*0/,iy /0/

        if (idum.le.0 .or. iy.eq.0) then
                idum=max(-idum,1)
                do j=NTAB+8,1,-1
                   k=idum/IQ
                   idum=IA*(idum-k*IQ)-IR*k
                   if (idum.lt.0) idum=idum+IM
                   if (j.le.NTAB) iv(j)=idum
                enddo
                iy=iv(1)
        endif
        k=idum/IQ
        idum=IA*(idum-k*IQ)-IR*k
        if(idum.lt.0) idum=idum+IM
        j=1+iy/NDIV
        iy=iv(j)
        iv(j)=idum
        ranff=min(AM*iy,RNMX)

        return
        END

!     subroutine random
!RANDOM
      subroutine random(iseedf,iyf,ivf)
!
!   Use this routine to time step the random forcing and
!   generate the random, incompressible, real-valued body
!   force (a1r,a2r,a3r) for the next time level.
!   Calculations are in spectral space.
!
      implicit real(a-h,o-z)
!
      include 'params.inc'
!RANDOM
      parameter (NTAB=32)
      parameter ( n1pp=n1+2)
!
      dimension fac1(n1pp),fac2(n2),fac3(n3)
      dimension fsq1(n1pp),fsq2(n2),fsq3(n3)
!
      dimension a1r(6,5,5), a2r(6,5,5), a3r(6,5,5),	&
               b1r(6,5,5), b2r(6,5,5), b3r(6,5,5)
      dimension fr(6,5,5,3)
!RANDOM
      integer iseedf,iyf,ivf(NTAB)
!     double precision rand
!     external function rand
!
      common /noise/a1r,a2r,a3r,b1r,b2r,b3r
      common /waveno/fac1,fac2,fac3,fsq1,fsq2,fsq3
      common /parms/xkf,var,tf,dt
!
      pi2=12.0*asin(0.5)
      smallx=1.0e-18
!
!   Generate white noise
      do mc=1,3
        do  k=1,5
          do  j=1,5
            do  i=1,6,2
              ip=i+1
!             cont1=ranf()
!             cont1=rand()
!RANDOM
              cont1 = ranff(iseedf,iyf,ivf)
              if(cont1.lt.1.e-6)cont1=1.e-6
              rmag=sqrt(-4.0*var*dt*alog(cont1)/tf)
!             rph=ranf()
!             rph=rand()
!RANDOM
              rph =  ranff(iseedf,iyf,ivf)
              fr(i ,j,k,mc)=rmag*cos(pi2*rph)
              fr(ip,j,k,mc)=rmag*sin(pi2*rph)
!           write(75,*)cont1,rph
            end do
          end do
        end do
       end do
!
!   Overwrite components for k1=0 to ensure real forcing
!   in physical space with no spatial mean.
      do mc=1,3
        do  k=1,5
          kk=7-k
          if (k.eq.1) kk=1
          do  j=2,3
            jj=7-j
            fr(1,jj,kk,mc)=  fr(1,j,k,mc)
            fr(2,jj,kk,mc)= -fr(2,j,k,mc)
          end do
        end do
        do  k=2,3
          kk=7-k
          fr(1,1,kk,mc)=  fr(1,1,k,mc)
          fr(2,1,kk,mc)= -fr(2,1,k,mc)
        end do
        fr(1,1,1,mc)=0.0
        fr(2,1,1,mc)=0.0
      end do
!
!  Euler step for Langevin type process with the
!  noise forcing. Then obtain incompressible part
!  and save as (a1r,a2r,a3r). Possible wavenos. are
!  k1=0,1,2 and k2,k3 = 0,1,2,-2,-1 .
      do  k=1,5
        kk=k
        if (k.gt.3) kk=k+n3-5
        do  j=1,5
        jj=j
        if (j.gt.3) jj=j+n2-5
          do  i=1,6
            rad=sqrt( fsq1(i)+fsq2(jj)+fsq3(kk))
            if ( rad.lt.xkf ) then
              b1r(i,j,k)=b1r(i,j,k)*(1.0-dt/tf)+fr(i,j,k,1)
              b2r(i,j,k)=b2r(i,j,k)*(1.0-dt/tf)+fr(i,j,k,2)
              b3r(i,j,k)=b3r(i,j,k)*(1.0-dt/tf)+fr(i,j,k,3)
              temp=1.0/((rad+smallx)**2)
              divb=temp*(fac1(i )*b1r(i,j,k)+		&
                        fac2(jj)*b2r(i,j,k)+		&
                        fac3(kk)*b3r(i,j,k))
              a1r(i,j,k)=b1r(i,j,k)-fac1(i )*divb
              a2r(i,j,k)=b2r(i,j,k)-fac2(jj)*divb
              a3r(i,j,k)=b3r(i,j,k)-fac3(kk)*divb
            else
              a1r(i,j,k)=0.0
              a2r(i,j,k)=0.0
              a3r(i,j,k)=0.0
            end if
          end do
        end do
      end do
!
!   Eliminate forcing of a mean flow at zero
!   waveno.
      a1r(1,1,1)=0.0
      a1r(2,1,1)=0.0
      a2r(1,1,1)=0.0
      a2r(2,1,1)=0.0
      a3r(1,1,1)=0.0
      a3r(2,1,1)=0.0
!
!   Ensure that the forcing is real-valued. This step
!   should only clean up rounding errors if fr is ok.
      do  k=1,5
        kk=7-k
        if ( k.eq.1) kk=1
        do  j=2,3
          jj=7-j
          a1r(1,jj,kk)=  a1r(1,j,k)
          a1r(2,jj,kk)= -a1r(2,j,k)
          a2r(1,jj,kk)=  a2r(1,j,k)
          a2r(2,jj,kk)= -a2r(2,j,k)
          a3r(1,jj,kk)=  a3r(1,j,k)
          a3r(2,jj,kk)= -a3r(2,j,k)
        end do
      end do
!
      do  k=2,3
        kk=7-k
        a1r(1,1,kk)=  a1r(1,1,k)
        a1r(2,1,kk)= -a1r(2,1,k)
        a2r(1,1,kk)=  a2r(1,1,k)
        a2r(2,1,kk)= -a2r(2,1,k)
        a3r(1,1,kk)=  a3r(1,1,k)
        a3r(2,1,kk)= -a3r(2,1,k)
      end do
      return
      end
!

!
      subroutine clean (cur)
!
!   For k_1=0 ensure that rounding and other errors are corrected
!   to ensure the physical data is real valued. For k1>0 this is
!   automatic by the FFT's but must be explicitly checked for k1=0.
!
      implicit real(a-h,o-z)
!
      include 'params.inc'
      parameter ( n1p=n1+1, n1pp=n1+2)
      parameter ( n2h=n2/2, n2hp=n2h+1, n3h=n3/2, n3hp=n3h+1)
!
      dimension cur(n1pp,n2,n3)
!     real*4, dimension(:,:,:), allocatable :: cur
!     allocate ( cur(n1pp,n2,n3) )
!
!   Real-value check for k2,k3 plane but k2 not zero
!
      do  k=1,n3
        kk=n3+2-k
        if ( k .eq. 1 ) kk=1
        do  j=2,n2h
          jj=n2+2-j
          cur(1,j,k)=0.5*(cur(1,j,k)+cur(1,jj,kk))
          cur(2,j,k)=0.5*(cur(2,j,k)-cur(2,jj,kk))
          cur(1,jj,kk)=  cur(1,j,k)
          cur(2,jj,kk)= -cur(2,j,k)
        end do
      end do
!
!   Real-value check for k1=k2=0 and k3 not zero
!
      do  k=2,n3h
        kk=n3+2-k
        cur(1,1,k)=0.5*(cur(1,1,k)+cur(1,1,kk))
        cur(2,1,k)=0.5*(cur(2,1,k)-cur(2,1,kk))
        cur(1,1,kk)=  cur(1,1,k)
        cur(2,1,kk)= -cur(2,1,k)
      end do
!
!   Truncate the n/2 wave number as this is noise
!   usually. This is separate from possible filtering
!   later on.     
!
      do  k=1,n3
        do  i=1,n1pp
          cur(i,n2hp,k)=0.0
        end do
      end do
!
      do  j=1,n2
        do  i=1,n1pp
          cur(i,j,n3hp)=0.0
        end do
      end do
!
      do  k=1,n3
        do  j=1,n2
          cur(n1p,j,k)=0.0
          cur(n1pp,j,k)=0.0
        end do
      end do
!
      return
      end
!
!
      subroutine value( u, up, icpt, lhnode,bg )
!
!   This routine does the 6-pt Lagrange interpolation of
!   the data field u using the coefficients set up in the
!   main program. One component field is done at a time.
!
      implicit real(a-h,o-z)
!
      include 'params.inc'
      parameter ( n1pp=n1+2)
!
      dimension u(n1pp,n2,n3)
!     real*4, dimension(:,:,:), allocatable :: u
      dimension up(npart,3), bg(npart,6,3), lhnode(npart,3)
      dimension ix(6),iy(6),iz(6)
!
!     common /interp/lhnode,bg
!
!     allocate ( cur(n1pp,n2,n3) )
!
!   Nodes are labelled 1 to 6 with the particle centered
!   between nodes 3 and 4. lhnode is an index to locate
!   node 3 on the grid, and has a value between 0 and n-1.
!   The code below keeps jx in the range 0 to n1-1, and then
!   increments by 1 to give ix, which matches the standard
!   indexing range of 1 to n1.
!
      do  ip=1,npart
        up(ip,icpt)=0.0
      end do
!
      mn1=n1
      mn2=n2
      mn3=n3

      do  ip=1,npart

              jx=lhnode(ip,1)-2
              ix(1)=( mod(jx,mn1) + mn1*(1-isign(1,jx))/2) + 1
!             ix(1) = jx
!             if(jx.lt.0) ix(1) = jx + mn1
!             if(jx.ge.mn1) ix(1) = jx - mn1
!             ix(1) = ix(1) + 1

              if(ix(1).le.(mn1-5) )then
              ix(2) = ix(1) + 1
              ix(3) = ix(2) + 1
              ix(4) = ix(3) + 1
              ix(5) = ix(4) + 1
              ix(6) = ix(5) + 1
              else
              ix(2) = mod(ix(1),mn1) + 1 
              ix(3) = mod(ix(2),mn1) + 1 
              ix(4) = mod(ix(3),mn1) + 1 
              ix(5) = mod(ix(4),mn1) + 1 
              ix(6) = mod(ix(5),mn1) + 1 
              endif

              jy=lhnode(ip,2)-2
              iy(1)=( mod(jy,mn2) + mn2*(1-isign(1,jy))/2) + 1
!             iy(1) = jy
!             if(jy.lt.0) iy(1) = jy + mn2
!             if(jy.ge.mn2) iy(1) = jy - mn2
!             iy(1) = iy(1) + 1

              if(iy(1).le.(mn2-5) )then
              iy(2) = iy(1) + 1
              iy(3) = iy(2) + 1
              iy(4) = iy(3) + 1
              iy(5) = iy(4) + 1
              iy(6) = iy(5) + 1
              else
              iy(2) = mod(iy(1),mn2) + 1
              iy(3) = mod(iy(2),mn2) + 1 
              iy(4) = mod(iy(3),mn2) + 1
              iy(5) = mod(iy(4),mn2) + 1
              iy(6) = mod(iy(5),mn2) + 1
              endif


              jz=lhnode(ip,3)-2
              iz(1)=( mod(jz,mn3) + mn3*(1-isign(1,jz))/2) + 1
!             iz(1) = jz
!             if(jz.lt.0) iz(1) = jz + mn3
!             if(jz.ge.mn3) iz(1) = jz - mn3
!             iz(1) = iz(1) + 1

              if(iz(1).le.(mn3-5) ) then
              iz(2) = iz(1) + 1
              iz(3) = iz(2) + 1
              iz(4) = iz(3) + 1
              iz(5) = iz(4) + 1
              iz(6) = iz(5) + 1
              else
              iz(2) = mod(iz(1),mn3) + 1
              iz(3) = mod(iz(2),mn3) + 1 
              iz(4) = mod(iz(3),mn3) + 1
              iz(5) = mod(iz(4),mn3) + 1
              iz(6) = mod(iz(5),mn3) + 1
              endif

       do  ndx=1,6
       ix1 = ix(ndx)
         do  ndy=1,6
         iy1 = iy(ndy)
          do  ndz=1,6
          iz1 = iz(ndz)
      up(ip,icpt)=up(ip,icpt)+u(ix1,iy1,iz1)*bg(ip,ndx,1)*	&
                         bg(ip,ndy,2)*bg(ip,ndz,3)
            end do
          end do
        end do
      end do
      return
      end
!
!
      subroutine ptdat(numb,npset,ttime,yp,vp,up,pertvel,pssp,vortp,yp0)
!
!   This routine gives short data analysis of the instantaneous
!   particle statistics such as mean rise velocity, etc. Data can be
!   reread later to form time averages.
!
      include 'params.inc'
      dimension yp(npart,3),vp(npart,3),up(npart,3),pssp(npart,3)
      dimension pertvel(npart,3)
      dimension vortp(npart,3),yp0(npart,3)
      dimension xbar1(3),x2bar1(3),disp(3)
      dimension xbar2(3),x2bar2(3)
      dimension xbar3(3),x2bar3(3)
      dimension xbar4(3),x2bar4(3)
!     common /part1/yp,vp,up,pssp,vortp
      common /boxsize/hx,hy,hz
!     common /location/yp0
!
!-Change 6/31/95
! 1 for x, 2 for r, 3 for total cross-axis dispersion
!
!   For each group of particles write out data averages to a
!   separate file. Read these later to obtain time averages.
!
!
      do  iset=1,npset
        iunit1=20+iset
        iunit2=30+iset
        jset=(iset-1)*numb
!
!  Estimate mean square particle displacement
!
          do ic=1,3
          sum1=0.0

          do  ip=1,numb
          jp=ip+jset
          sum1 =sum1+(yp(jp,ic)-yp0(jp,ic))**2
          end do

          disp(ic)=sum1/float(numb)
          enddo
!       write(iunit,101)disp
!
!   Estimate mean and mean square particle velocities
!
        do  ic=1,3
          sum1=0.0
          sum2=0.0

          do  ip=1,numb
            jp=ip+jset
            sum1=sum1+vp(jp,ic)
            sum2=sum2+vp(jp,ic)**2
          end do

          xbar1(ic)=sum1/float(numb)
          x2bar1(ic)=sum2/float(numb)
        end do

!       write(iunit,101) xbar1,x2bar1
!
!   Estimate mean and mean square flow velocities
!
        do  ic=1,3
          sum1=0.0
          sum2=0.0

          do  ip=1,numb
            jp=ip+jset
            sum1=sum1+up(jp,ic)
            sum2=sum2+up(jp,ic)**2
          end do

          xbar2(ic)=sum1/float(numb)
          x2bar2(ic)=sum2/float(numb)
        end do
!       write(iunit,101) xbar2,x2bar2
!
!   Estimate mean square relative velocities
!
        do  ic=1,3
          sum2=0.0

          do  ip=1,numb
            jp=ip+jset
            sum2=sum2+(vp(jp,ic)-up(jp,ic)-pertvel(jp,ic))**2
          end do

          x2bar3(ic)=sum2/float(numb)
        end do
!       write(iunit,101) x2bar3
!
!   Estimate mean and mean square values of pssp cpts.
!
!       do  ic=1,3
!         xbar(ic)=0.0
!         x2bar(ic)=0.0
!         do  ip=1,numb
!           jp=ip+jset
!           xbar(ic)=xbar(ic)+pssp(jp,ic)
!           x2bar(ic)=x2bar(ic)+pssp(jp,ic)**2
!         end do
!         xbar(ic)=xbar(ic)/float(numb)
!         x2bar(ic)=x2bar(ic)/float(numb)
!       end do
!       write(iunit,101) xbar,x2bar
!
!   Estimate mean and mean square vorticity
!
        do  ic=1,3
          sum1=0.0
          sum2=0.0

          do  ip=1,numb
            jp=ip+jset
            sum1=sum1+vortp(jp,ic)
            sum2=sum2+vortp(jp,ic)**2
          end do

          xbar4(ic)=sum1/float(numb)
          x2bar4(ic)=sum2/float(numb)
        end do

!-Change 6/31/95
        write(iunit1,101) ttime,disp,xbar1,x2bar1,x2bar3
        write(iunit2,101) ttime,xbar2,x2bar2,xbar4,x2bar4
!
      end do
!
 101  format(2x, 13(1pe12.4) )
!
      return
      end
!
!
!
      subroutine PERTURBVELOC(vp,yp,up,pertvel,nnodes)


      implicit real(a-h,o-z)

      include 'params.inc'

      integer head(npvcd,npvcd,npvcd),list(npart)
      integer headd(npvcd,npvcd,npvcd,0:nnodes-1)
      integer ncheck(npart),my_thread
      dimension vp(npart,3),up(npart,3),yp(npart,3)
      dimension pertvel(npart,3),pertvelold(npart,3),r21(3)
      dimension pvint(npart,3)
      common /geomrad/rad1,rad2,rcoll11,rcoll12,rcoll22,  &
                        numb,wpvcd
      common /particle/w0diff
      common /boxsize/hx,hy,hz
      common /tiempo/ttime,istep,nhalt
      common /machine/nproc

        write(*,*) 'I am entering PERTURBVELOC'

	xmaxerror=errorval

	write(*,*) 'errorval =',errorval

      do ic0=1,3
        do ip0=1,npart
        pertvelold(ip0,ic0)=pertvel(ip0,ic0)
        pvint(ip0,ic0)=pertvel(ip0,ic0)
        pertvel(ip0,ic0)=0.0
        end do
      end do

        T1=TIMEF( )/1000.00

        my_thread = 1
        do i=1,npvcd
        do j=1,npvcd
        do k=1,npvcd
        headd(i,j,k,my_thread)=0
        enddo
        enddo
        enddo
        T2=TIMEF( )/1000.00
        Thead=T2-T1

        T1=TIMEF( )/1000.00

      my_thread=1
      do ip=1,npart
      ix=1+int(yp(ip,1)/wpvcd)
      iy=1+int(yp(ip,2)/wpvcd)
      iz=1+int(yp(ip,3)/wpvcd)
!
      if(ix.gt.npvcd ) then
      ix = npvcd
      endif

      if(iy.gt.npvcd)then
      iy = npvcd
      endif

      if(iz.gt.npvcd)then
      iz = npvcd
      endif
!
        list(ip)=headd(ix,iy,iz,my_thread)
        headd(ix,iy,iz,my_thread)=ip
      end do

        do l1=1,npvcd
        do l2=1,npvcd
        do l3=1,ncd
        head(l1,l2,l3)=headd(l1,l2,l3,1)
        enddo
        enddo
        enddo

        T2=TIMEF( )/1000.00
        Theadlist=T2-T1

        write(*,*) 'Thead,Theadlist'
        write(*,999) Thead,Theadlist
999     format(2f10.6)

      nnsu=0

100    continue

      do ic=1,npart
        ncheck(ic)=0
      enddo


      DO ip1=1,NPART
         ix=1+int(yp(ip1,1)/wpvcd)
         iy=1+int(yp(ip1,2)/wpvcd)
         iz=1+int(yp(ip1,3)/wpvcd)
!
      if(ix.gt.npvcd ) then
!      write(103,*)ix,iy,iz,ip
!      write(103,*)ttime,yp(ip,1),yp(ip,2),		&
!        yp(ip,3),wpvcd,hx,hy,hz,'PV subroutine'
      ix = npvcd
      endif

      if(iy.gt.npvcd)then
!      write(103,*)ix,iy,iz,ip
!      write(103,*)ttime,yp(ip,1),yp(ip,2),		&
!        yp(ip,3),wpvcd,hx,hy,hz,'PV subroutine'
      iy = npvcd
      endif

      if(iz.gt.npvcd)then
!      write(103,*)ix,iy,iz,ip
!      write(103,*)ttime,yp(ip,1),yp(ip,2),		&
!        yp(ip,3),wpvcd,hx,hy,hz,'PV subroutine'
      iz = npvcd
      endif
!
        do i=-1,1
        do j=-1,1
        do k=-1,1

        ix2=ix+i
        sx=0.0
        if(ix2.lt.1) then
        ix2=npvcd
        sx=-1.
        endif
        if(ix2.gt.npvcd)then
        ix2=1
        sx=1.
        endif
!
        iy2=iy+j
        sy=0.0
        if(iy2.lt.1) then
        iy2=npvcd
        sy=-1.
        endif
        if(iy2.gt.npvcd)then
        iy2=1
        sy=1.
        endif
!
        iz2=iz+k
        sz=0.0
        if(iz2.lt.1) then
        iz2=npvcd
        sz=-1.
        endif
        if(iz2.gt.npvcd)then
        iz2=1
        sz=1.
        endif

        ip2=head(ix2,iy2,iz2)
        if(ip2.eq.0)goto 201

199     continue

       IF (ip1.ne.ip2) THEN

       ypb1=yp(ip2,1)+sx*hx
       ypb2=yp(ip2,2)+sy*hy
       ypb3=yp(ip2,3)+sz*hz

       r21(1)= yp(ip1,1) - ypb1
       r21(2)= yp(ip1,2) - ypb2
       r21(3)= yp(ip1,3) - ypb3
       dnij=sqrt(r21(1)**2+r21(2)**2+r21(3)**2) 

       if (ip2.le.numb) then
         rad=rad1
       else
         rad=rad2
       endif

        if(ip1.le.numb) then
                if(ip2.le.numb) then
                     rcoll=rcoll11
                else
                     rcoll=rcoll12
                endif
        endif
        if(ip1.gt.numb) then
                if(ip2.gt.numb) then
                     rcoll=rcoll22
                else
                     rcoll=rcoll12
                endif
        endif
       if(dnij.lt.rcoll) then 
       write(*,*) 'WARNING dnij<rcoll ',dnij,rcoll,ip1,ip2
       endif

       if((dnij/rad).lt.reldist) then

       fact1=rad/dnij
       fact2=fact1**3
       fact1=0.75*fact1  

       al1=fact1-0.75*fact2
       al2=fact1+0.25*fact2

       relv1=vp(ip2,1)-up(ip2,1)-pvint(ip2,1)
       relv2=vp(ip2,2)-up(ip2,2)-pvint(ip2,2)
       relv3=vp(ip2,3)-up(ip2,3)-pvint(ip2,3)  

        aloop=relv1*r21(1)+relv2*r21(2)+relv3*r21(3)
        aloop=al1*aloop/dnij**2

       pertvel(ip1,1)=pertvel(ip1,1)+aloop*r21(1)+al2*relv1
       pertvel(ip1,2)=pertvel(ip1,2)+aloop*r21(2)+al2*relv2
       pertvel(ip1,3)=pertvel(ip1,3)+aloop*r21(3)+al2*relv3

       
       endif

       ENDIF

200     ip2=list(ip2)
        if(ip2.eq.0)goto 201
        goto 199

201     continue

      enddo
      enddo
      enddo

      pvint(ip1,1)=pertvel(ip1,1)
      pvint(ip1,2)=pertvel(ip1,2) 
      pvint(ip1,3)=pertvel(ip1,3)

      ENDDO

! CONVERGENCY CRITERIA
      do ie=1,npart
       error1=(pertvel(ie,1)-pertvelold(ie,1))/w0diff
       error2=(pertvel(ie,2)-pertvelold(ie,2))/w0diff
       error3=(pertvel(ie,3)-pertvelold(ie,3))/w0diff

       IF((istep.eq.10).or.(istep.eq.(nhalt/2))		&
                      .or.(istep.eq.(nhalt-100))) THEN

       if(istep.eq.10)nfile=325
       if(istep.eq.(nhalt/2))nfile=328
       if(istep.eq.(nhalt-100))nfile=331

       if(ie.eq.npart) then
       write(nfile+1,*) istep,nnsu,error1,error2,error3 
       endif

       if(ie.eq.(npart/2)) then
       write(nfile+2,*) istep,nnsu,error1,error2,error3 
       endif

       if(ie.eq.1) then
       write(nfile+3,*) istep,nnsu,error1,error2,error3 
       endif

       ENDIF

        if((abs(error1).le.xmaxerror).and.	&
          (abs(error2).le.xmaxerror).and.	&
          (abs(error3).le.xmaxerror)) then
         ncheck(ie)=0
        else
         ncheck(ie)=1
        endif
      enddo

      nsumcheck=0
      do ipo=1,npart
        nsumcheck=nsumcheck+ncheck(ipo)
      enddo

      if( (nnsu.ge.( nmaxiter/2   )) .and.	&
         (nnsu.lt.((nmaxiter/2)+1)) ) then
          xmaxerror=xmaxerror*10.
      endif

      if(nnsu.ge.(nmaxiter-2)) then
          xmaxerror=xmaxerror*10.
      endif


!	write(*,*) xmaxerror

      if((nsumcheck.eq.0).or.(nnsu.ge.nmaxiter))goto 203

!        write(*,*) 'nsumcheck=', nsumcheck

           nnsu=nnsu+1   
	do ic=1,3
           do ip=1,npart
             pertvelold(ip,ic)=pertvel(ip,ic)
             pertvel(ip,ic)=0.0
           enddo
	enddo

      goto 100

203   continue

      IF (nsumcheck.ne.0) then
 
       write(*,*) 'WARNING ***** nnsu= ',nnsu,' and nsumcheck= '	&
           ,nsumcheck, ' xmaxerror= ',xmaxerror

!      do ipo=1,npart
!        if (ncheck(ipo).eq.1)then
!        write(*,*) 'Particle with problem= ',ipo
!        write(*,*) 'Position= ',yp(ipo,1),yp(ipo,2),yp(ipo,3)
!        write(*,*) 'PertVelo= ',pertvel(ipo,1),pertvel(ipo,2)		&
!                              ,pertvel(ipo,3)
!        write(*,*) 'PertVeloOLD= ',pertvelold(ipo,1),pertvelold(ipo,2)	&
!                              ,pertvelold(ipo,3)
!        endif
!      enddo
  
      ENDIF
        write(*,*) 'I am leaving PERTURBVELOC',nnsu

      return
      end

!==================================================================
      SUBROUTINE gaussian(u,pi)

      INCLUDE 'params.inc'

      integer                                      :: n1pp
      integer                                      :: idum
      REAL,DIMENSION   (n1+2,n2,n3)  :: u

      idum = -8765497213.
      n1pp = n1+2

      do  k=1,n3
        do  i=1,n1pp,2
          ipp=i+1
          do  j=1,n2
            t1 = randfs(idum)
            t2 = randfs(idum)
            if(t1.le.1.e-10) t1 = 1.e-10
            if(t2.le.1.e-10) t2 = 1.e-10
            t2 = 2.*pi*t2
            u(i,j,k) = sqrt(-2.0*log(t1))*cos(t2)
            u(ipp,j,k)= sqrt(-2.0*log(t1))*sin(t2)
          end do
        end do
      end do

      RETURN
      END SUBROUTINE gaussian


!==================================================================

! random generator

      Function randfs(idum)

      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL randfs,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,IA1=40014,IA2=40692,IQ1=53668,IQ2=52774 &
        ,IR1=12211, &
      IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)

! Long period (> 2.1018) random number generator of L'Ecuyer with Bays-Durham shu.e
! and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive
! of the endpoint values). Call with idum a negative integer to initialize; thereafter, do not
! alter idum between successive deviates in a sequence. RNMX should approximate the largest
! floating value that is less than 1.

      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/

      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum

        do j=NTAB+8,1,-1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
        enddo

        iy=iv(1)
      endif

      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1

      if (idum.lt.0) idum=idum+IM1

      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2

      if (idum2.lt.0) idum2=idum2+IM2

      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1

      randfs=min(AM*iy,RNMX)

      return
      END

!==================================================================

!==================================================================

      subroutine supfordet (cur,cvr,cwr,fac1,fac2,fac3,force)

      include 'params.inc'

      real*4, dimension(n1+2,n2,n3)      :: cur,cvr,cwr
      real*4, dimension(n1+2)            :: fac1
      real*4, dimension(n2)              :: fac2
      real*4, dimension(n3)              :: fac3
      real*4                             :: k2,force(2),tmp(2),tots1,tots2,tot,ek
      integer                            :: ik2

      n1pp=n1+2

        tots1 = 0.0
        tots2 = 0.0
        do  k=1,n3
          do  i=1,5,2
            do  j=1,n2

              k2 = real( fac1(i)**2 + fac2(j)**2 + fac3(k)**2 )
              ik2 = int( sqrt(k2)+0.5 )

              if(ik2.eq.1)then
              tot = cur(i,j,k)**2 + cur(i+1,j,k)**2 +         &
                    cvr(i,j,k)**2 + cvr(i+1,j,k)**2 +         &
                    cwr(i,j,k)**2 + cwr(i+1,j,k)**2
              if (i.eq.1) tot = tot / 2.
              tots1 = tots1 + tot
              end if

              if(ik2.eq.2)then
              tot = cur(i,j,k)**2 + cur(i+1,j,k)**2 +         &
                    cvr(i,j,k)**2 + cvr(i+1,j,k)**2 +         &
                    cwr(i,j,k)**2 + cwr(i+1,j,k)**2
              if (i.eq.1) tot = tot / 2.
              tots2 = tots2 + tot
              end if

            end do
          end do
        end do
        tmp(1) = tots1
        tmp(2) = tots2

        do  kk=1,5
          do  i=1,5,2
            do jj=1,5

              j=jj
              k=kk
              if (kk.eq.4) k=n3-1
              if (jj.eq.4) j=n2-1
              if (kk.eq.5) k=n3
              if (jj.eq.5) j=n2

              k2 = real( fac1(i)**2 + fac2(j)**2 + fac3(k)**2 )
              ik2 = int( sqrt(k2)+0.5 )

              do is=1,2
                if (ik2.eq.is) then

                  ek = sqrt(force(is)/tmp(is))
                  cur(i,j,k)  = cur(i,j,k)  *ek
                  cur(i+1,j,k)= cur(i+1,j,k)*ek
                  cvr(i,j,k)  = cvr(i,j,k)  *ek
                  cvr(i+1,j,k)= cvr(i+1,j,k)*ek
                  cwr(i,j,k)  = cwr(i,j,k)  *ek
                  cwr(i+1,j,k)= cwr(i+1,j,k)*ek

                endif
              enddo

            end do
          end do
        end do

       return
       end subroutine supfordet

!==================================================================

      SUBROUTINE supforstoch(cur,cvr,cwr,iseedf,iyf,ivf,dt)

      implicit real(a-h,o-z)

      include 'params.inc'
      REAL*4,DIMENSION (n1+2,n2,n3) :: cur, cvr, cwr
      dimension a1r(6,5,5), a2r(6,5,5), a3r(6,5,5),     &
               b1r(6,5,5), b2r(6,5,5), b3r(6,5,5)
      real    dt
      integer iseedf,iyf,ivf(32)
      common /noise/a1r,a2r,a3r,b1r,b2r,b3r

      call random(iseedf,iyf,ivf)
      do  k=1,5
        kk=k
        if (k.gt.3) kk=k+n3-5
        do  j=1,5
          jj=j
          if (j.gt.3) jj=j+n2-5
          do  i=1,6
             cur(i,jj,kk)=cur(i,jj,kk)+dt*a1r(i,j,k)
             cvr(i,jj,kk)=cvr(i,jj,kk)+dt*a2r(i,j,k)
             cwr(i,jj,kk)=cwr(i,jj,kk)+dt*a3r(i,j,k)
          end do
        end do
      end do

      RETURN
      END SUBROUTINE supforstoch

!==================================================================

      SUBROUTINE projection(cur,cvr,cwr,fac1,fac2,fac3)

      include 'params.inc'

      REAL*4,DIMENSION (n1+2,n2,n3)   :: cur, cvr, cwr
      REAL*4,DIMENSION (n1+2)         :: fac1
      REAL*4,DIMENSION (n2)           :: fac2
      REAL*4,DIMENSION (n3)           :: fac3
      REAL                            :: k2 ,tmpr, tmpi


      do  k=1,n3
       do  i=1,n1+2,2
        do  j=1,n2

         k2 = real ( fac1(i)**2 + fac2(j)**2 + fac3(k)**2 )

         tmpr = fac1(i) * cur(i,j,k)   + fac2(j) * cvr(i,j,k)   + fac3(k) * cwr(i,j,k)
         tmpi = fac1(i) * cur(i+1,j,k) + fac2(j) * cvr(i+1,j,k) + fac3(k) * cwr(i+1,j,k)
         if (k2.gt.0.5) then
          tmpr = tmpr / k2
          tmpi = tmpi / k2
         endif

         cur(i,j,k)   = cur(i,j,k)   - fac1(i) * tmpr
         cur(i+1,j,k) = cur(i+1,j,k) - fac1(i) * tmpi
         cvr(i,j,k)   = cvr(i,j,k)   - fac2(j) * tmpr
         cvr(i+1,j,k) = cvr(i+1,j,k) - fac2(j) * tmpi
         cwr(i,j,k)   = cwr(i,j,k)   - fac3(k) * tmpr
         cwr(i+1,j,k) = cwr(i+1,j,k) - fac3(k) * tmpi

        enddo
       enddo
      enddo

      RETURN
      END SUBROUTINE projection

!==================================================================
