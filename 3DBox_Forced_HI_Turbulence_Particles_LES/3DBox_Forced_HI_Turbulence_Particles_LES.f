c           
c            Implement fftw 
c
c    Periodic box code to study the geometric collision rates between two groups
c    of particles in an isotropic turbulence or other periodic flow. 
c    The turbulence is (locally) homogeneous and isotropic with forcing at
c    at low wavenumbers to sustain the flow. An AB2/CN2 scheme is used to advance the flow.
c    Calculations are in Fourier space and pseudospectral.
c
c       The particles move in this version as 
c          dV/dt = ( u(Y,t) - V(t) + W )/ tau .
c    The parameters tau and W characterise the motion. Values of
c    u  are found by Lagrange interpolation. Axes are fixed
c    so that  W = ( wp , 0, 0 ) and wp >0 is the still fluid settling
c    velocity.
c
c    Collision detections are used to count individual collision events.
c
c    OpenMP implemented for SGI/Origin at NCAR.
c
c     VERSION March, 2003
c
c	IMCLUDE HYDRODYNAMIC INTERACTIONS  (09/16/03)
c
      USE OMP_LIB  
      implicit real(a-h,o-z)
c
      include 'params.inc'
!      include 'plans.h'
c
      parameter(n1p=n1+1)
      parameter(n1pp=n1+2,n2d=n2*2,n3d=n3*2)
      parameter(nwds=n1pp*n2*n3,npts=n1*n2*n3)
      parameter(n1h=n1/2,n2h=n2/2,n3h=n3/2)
      parameter(n1hp=n1h+1,n2hp=n2h+1,n3hp=n3h+1)
      parameter(nparthalf=npart/2)
CRANDOM
      parameter (NTAB=32)
c

      real*4, dimension(:,:,:), allocatable :: cur,cvr,cwr
      real*4, dimension(:,:,:), allocatable :: tur,tvr,twr

      real*4, dimension(:,:,:), allocatable :: tur_filtered
     1        ,tvr_filtered,twr_filtered
      real*4, dimension(:,:,:), allocatable :: sxp,syp,szp
      real*4, dimension(:,:,:), allocatable :: sxr,syr,szr
      real*4, dimension(:,:,:), allocatable :: sxr_filtered
     1        ,syr_filtered,szr_filtered
      real*4, dimension(:,:,:), allocatable :: filter,tempa
c
      dimension fac1(n1pp), fac2(n2), fac3(n3)
      dimension fsq1(n1pp), fsq2(n2), fsq3(n3)
      dimension gr11(0:179),gr12(0:179),gr22(0:179)
      real      gr11w,gr12w,gr22w
c
      dimension probe(8),efx(6),phi1(n1pp),phi2(n1pp),
     1   phi3(n1pp),phi4(n1pp),phi5(n1pp),vdist(200)
      dimension  eddy_time(n1pp)

c
c   Particle arrays.
c   ----------------
      real*4, dimension(:), allocatable ::tau,wp,dmoved
      real*4, dimension(:,:), allocatable :: yp,ypp,vp,up,dmove
      real*4, dimension(:,:), allocatable :: pertvel
      real*4, dimension(:,:), allocatable :: vortp,pssp,yp0
      real*4, dimension(:,:,:), allocatable :: fvp,fdp,bg
      integer, dimension(:,:), allocatable :: lhnode
      integer, dimension(:), allocatable :: removpart
      integer, dimension(:), allocatable :: list
      integer*8, dimension(:), allocatable :: plan1DF,plan1DB
      integer*8, dimension(:), allocatable :: plan2DRC,plan2DCR
      integer, dimension(:,:), allocatable :: nir11n,nir12n,nir22n
      real*4, dimension(:,:), allocatable :: wrbin11n,wrbin12n,wrbin22n
      real*4, dimension(:), allocatable :: sgs_visck ,sgs_visck_old
      integer, dimension(:,:,:), allocatable :: con1,con2  

      dimension dxyz(3)
c
      integer head(ncd,ncd,ncd)
      integer nir11(0:179),nir12(0:179),nir22(0:179)
      real wrbin11(0:179),wrbin12(0:179),wrbin22(0:179)
!      integer con1(128,128,128),con2(128,128,128)
      real pc1(0:30),pc2(0:30),pcr(0:30) !pcr is probability distribution of random
      real dc1,dc2, Factorial_Ic
      real average_con,common_factor
      integer Ng
      real Ck,rm  !Kolomogorov constant obtained from DNS
      real sumk,sumlogei,sumlogki,sumlogki2,sumlogeilogki
      real C_kolm, totaldiss, eta_les,f_eta_k,x,y,z,thetakpq,btemp,temp2
      integer  ip2,ip, iq, iq2,model_sgs
      real Tbegin,Tend,Tduration
c
c   Flow forcing arrays.
c   --------------------
      dimension a1r(6,5,5), a2r(6,5,5), a3r(6,5,5),
     &          b1r(6,5,5), b2r(6,5,5), b3r(6,5,5)
c
c
c     double precision rand
c     external function rand,srand
CRANDOM
       integer iseedp,iyp,ivp(NTAB)
       integer iseedf,iyf,ivf(NTAB)
       real*8 T1,T2
       character*22 fconc,fvort,fparticle
       character*50  readtemp     

!       character*24   fdate
!       external       fdate

      logical newflo, newpart, forced, pflow, HYDINT, NOTURB,
     1         OVERL
      logical fdns,les
      real cutoff_num, sgs_visc
      real sumabswr11,sumabswr12,sumabswr22
c
      common /waveno/fac1,fac2,fac3,fsq1,fsq2,fsq3
      common /noise/a1r,a2r,a3r,b1r,b2r,b3r
      common /parms/xkf,var,tf,dt
c     common /interp/lhnode,bg
c     common /part1/yp,vp,up,pssp,vortp
c     common /part2/tau,wp
c     common /rannum/iseedf
      common /geomrad/rad1,rad2,rcoll11,rcoll12,rcoll22,numb,wpvcd
      common /particle/w0diff
      common /boxsize/hx,hy,hz
      common /tiempo/ttime,istep,nhalt
c     common /location/yp0
c
c   Initialize arrays for fft
      call OMP_SET_NUM_THREADS(8)

c$omp parallel shared(nnodes)
       nnodes =  omp_get_num_threads()
c$omp end parallel

       print*,"nnodes=",nnodes
       allocate ( plan1DF(0:(nnodes-1)) )
       allocate ( plan1DB(0:(nnodes-1)) )
       allocate ( plan2DRC(0:(nnodes-1)) )
       allocate ( plan2DCR(0:(nnodes-1)) )
c
c$omp parallel private(my_thread)
c$omp+ shared(plan2DRC,plan2DCR)
c$omp+ shared(plan1DF,plan1DB)

      my_thread = omp_get_thread_num()
       call set_plan2DRC(n1,n2,plan2DRC(my_thread))
       call set_plan1DF (n3,plan1DF(my_thread))
       call set_plan2DCR (n1,n2,plan2DCR(my_thread))
       call set_plan1DB (n3,plan1DB(my_thread))

c$omp end parallel
c
      pi=4.*atan(1.)
      pi2=2.*pi
c
      allocate ( tau(npart) )
      allocate ( wp(npart) )
      allocate ( dmoved(npart) )
      allocate ( yp(npart,3) )
      allocate ( ypp(npart,3) )
      allocate ( dmove(npart,3) )
      allocate ( vp(npart,3) )
      allocate ( up(npart,3) )
      allocate ( pertvel(npart,3) )
      allocate ( removpart(nparthalf) )
      allocate ( vortp(npart,3) )
      allocate ( pssp(npart,3) )
      allocate ( yp0(npart,3) )
      allocate ( lhnode(npart,3) )
      allocate ( fvp(npart,3,0:2) )
      allocate ( fdp(npart,3,0:3) )
      allocate ( bg(npart,6,3) )
      allocate ( list(npart) )
c      allocate ( ipa(npart,3) )
c      allocate ( ipb(npart,3) )
c
      allocate ( cur(n1pp,n2,n3) )
      allocate ( cvr(n1pp,n2,n3) )
      allocate ( cwr(n1pp,n2,n3) )

      allocate ( tur(n1pp,n2,n3) )
      allocate ( tvr(n1pp,n2,n3) )
      allocate ( twr(n1pp,n2,n3) )




      allocate ( sxp(n1pp,n2,n3) )
      allocate ( syp(n1pp,n2,n3) )
      allocate ( szp(n1pp,n2,n3) )
      allocate ( sxr(n1pp,n2,n3) )
      allocate ( syr(n1pp,n2,n3) )
      allocate ( szr(n1pp,n2,n3) )




      allocate ( filter(n1pp,n2,n3) )
      allocate ( tempa(n1pp,n2,n3) )
      allocate ( nir11n(0:179,0:nnodes-1) )
      allocate ( nir12n(0:179,0:nnodes-1) )
      allocate ( nir22n(0:179,0:nnodes-1) )
      allocate ( wrbin11n(0:179,0:nnodes-1) )
      allocate ( wrbin12n(0:179,0:nnodes-1) )
      allocate ( wrbin22n(0:179,0:nnodes-1) )
      allocate ( sgs_visck(n2h) )
      allocate ( sgs_visck_old(n2h) )
      allocate ( con1(128,128,128) )
      allocate ( con2(128,128,128) )

c   Set up Fourier grid
c   -------------------
c   hx, hy, hz = box size
c   fx, fy, fz = fundamental wave numbers
c   If n1 > n2 this corresponds to an elongated box in 
c   the vertical.
      hx = pi2
      hy = pi2
      hz = pi2
      DXYZ(1)=HX/real(N1)
      DXYZ(2)=HY/real(N2)
      DXYZ(3)=HZ/real(N3)
      fx = pi2/hx
      fy = pi2/hy
      fz = pi2/hz
c
c   fac1, fac2, fac3 = arrays of wave numbers
c   -----------------------------------------
c -  x : ( 0 --> n1h ) x fx in real/imag pairs
c -  y : ( 0 --> n2h, -(n2h-1) --> -1 ) x fy
c -  z : ( 0 --> n3h, -(n3h-1) --> -1 ) x fz
c
c   fsq1, fsq2, fsq3 = arrays of squares of wave numbers
c
c$omp parallel do
      do  i=1,n1pp,2
        fac1(i)=float((i-1)/2)*fx
        fsq1(i)=fac1(i)**2
        fac1(i+1)=fac1(i)
        fsq1(i+1)=fsq1(i)
      end do
c$omp end parallel do
c
c$omp parallel do private(jj)
      do  j=1,n2
        jj=j-1
        if ( j .gt. n2hp ) jj=jj-n2
        fac2(j)=float(jj)*fy
        fsq2(j)=fac2(j)**2
      end do
c$omp end parallel do
c
c$omp parallel do private(kk)
      do  k=1,n3
        kk=k-1
        if ( k .gt. n3hp ) kk=kk-n3
        fac3(k)=float(kk)*fz
        fsq3(k)=fac3(k)**2
      end do
c$omp end parallel do
c
c
c   Set up spectral filter
c   ----------------------
c   specify filter(i,j,k) here
c
!      fkstar=fac1(n1pp)-1.5
       fkstar=n2/3.
       itrunc=int(fkstar)
c$omp parallel do private(j,i,coeff)
      do  k=1,n3
        do  j=1,n2
          do  i=1,n1pp
            coeff=sqrt(fsq1(i)+fsq2(j)+fsq3(k))
            if ( coeff.lt.fkstar ) then
            filter(i,j,k)=1.0
            else
            filter(i,j,k)=0.0
            endif
          end do
        end do
      end do
c$omp end parallel do
c   Set up forcing parameters
c   -------------------------
c   This scheme assumes a cubic box!
c
c   xkf= forcing radius
c   var= variance of white noise
c   tf = time constant
c   xnf= no. of modes forced
      xkf=sqrt(8.0)
      var=447.31
      tf=0.038
      xnf=0.0
c
c   Calculate no. of modes forced
c
      index=0
      do  k=1,5
        kk=k
        if ( k.gt.3 ) kk=k+n3-5
        do  j=1,5
          jj=j
          if ( j.gt.3 ) jj=j+n2-5
          do  i=1,6,2
            rad=sqrt( fsq1(i)+fsq2(jj)+fsq3(kk))
            if ( rad.lt.xkf ) then
              if ( i.le.2) then
                index=index+1
              else
                index=index+2
              endif
            endif
          end do
        end do
      end do
      xnf=float(index-1)
c
c   Open data files for regular input/output
c   ----------------------------------------
c  input.dat contains basic information to control the job,
c  output.dat is a monitoring log file for the job.
c  probe.dat contains the short set of flow statistics.
c    
      open (unit=107, file='pc1.dat')
      open (unit=108, file='pc2.dat')

      open (unit=109, file='ncgrwr11.dat')
      open (unit=110, file='ncgrwr12.dat')
      open (unit=115, file='ncgrwr22.dat')


      open ( unit=103, file= 'index-check.dat')
      open ( unit=106, file= 'vpmean-var.dat')
      open ( unit=15, file='partin.dat')
      open ( unit=16, file='partout.dat')

      open ( unit=8, file='collision.dat')
      open ( unit=80,file='detecpairs.dat',status='unknown')
      open ( unit=71,file='Wr11.dat', status='unknown')
      open ( unit=81,file='Vp11.dat', status='unknown')
      open ( unit=72,file='Wr12.dat', status='unknown')
      open ( unit=82,file='Vp12.dat', status='unknown')
      open ( unit=73,file='Wr22.dat', status='unknown')
      open ( unit=83,file='Vp22.dat', status='unknown')

      open ( unit=74,file='Wr11c.dat', status='unknown')
      open ( unit=84,file='Vp11c.dat', status='unknown')
      open ( unit=75,file='Wr12c.dat', status='unknown')
      open ( unit=85,file='Vp12c.dat', status='unknown')
      open ( unit=76,file='Wr22c.dat', status='unknown')
      open ( unit=86,file='Vp22c.dat', status='unknown')

      open ( unit=111,file='nir11a.dat',status='unknown')
      open ( unit=112,file='nir12a.dat',status='unknown')
      open ( unit=113,file='nir22a.dat',status='unknown')

      open ( unit=211,file='nir11b.dat',status='unknown')
      open ( unit=212,file='nir12b.dat',status='unknown')
      open ( unit=213,file='nir22b.dat',status='unknown')

      open ( unit=311,file='nir11c.dat',status='unknown')
      open ( unit=312,file='nir12c.dat',status='unknown')
      open ( unit=313,file='nir22c.dat',status='unknown')

      open ( unit=121,file='wrbin11a.dat',status='unknown')
      open ( unit=122,file='wrbin12a.dat',status='unknown')
      open ( unit=123,file='wrbin22a.dat',status='unknown')

      open ( unit=221,file='wrbin11b.dat',status='unknown')
      open ( unit=222,file='wrbin12b.dat',status='unknown')
      open ( unit=223,file='wrbin22b.dat',status='unknown')

      open ( unit=321,file='wrbin11c.dat',status='unknown')
      open ( unit=322,file='wrbin12c.dat',status='unknown')
      open ( unit=323,file='wrbin22c.dat',status='unknown')
c

      Tbegin=SECNDS(0.0)
c   Define physical parameters
c   --------------------------
c
c   ndrag=1 for nonlinear drag force, 0 for Stokes drag
c
c   arad1=particle radius of group 1/dx
c   arad2=particle radius of group 2/dx
c
c   rnu = viscosity
c
c   dt = time step size 
c   dtr, dtp = weighted time steps for adams-bashforth
c   rnuf used for viscosity term used in time-stepping
c
c   newflo=.true. for start of new flow, .false. for restart
c   newpart=.true. for start of new particle system
c   forced=.true. if flow is driven by random forcing
c
c   istep = time step counter
c   nhalt = number of time step to be run
c   nlong= number of steps between long data output
c   nshort= number of steps between short data output
c
      dx=pi2/float(n2)
c
      read(15,*) readtemp, ndrag
      read(15,*) readtemp,rnu
      read(15,*) readtemp,dt
c
c
      dtr=1.5*dt
      dtp=0.5*dt
      rnuf=0.5*dt*rnu
      ttime=0.0
      smallx=1.0e-18
c
      read (15,*) readtemp, newflo
      read (15,*) readtemp, forced
c
      read (15,*) readtemp,newpart
      read (15,*) readtemp,pflow
      read (15,*) readtemp,HYDINT
      read (15,*) readtemp,NOTURB
      read (15,*) readtemp,nfactor
      read (15,*) readtemp,OVERL

	if(NOTURB) then
	  hx=hx/float(nfactor)
	  hy=hy/float(nfactor)
          hz=hz/float(nfactor)
	endif

 	if(HYDINT) OVERL=.false.

      read (15,*) readtemp,nhalt,nlong,nshort
      read (15,*) readtemp,itvort
      read (15,*) readtemp,ttime
c
      if(pflow) then
      open ( unit=18, file='probe.dat')
      open ( unit=19, file='spectrum.dat')
      open ( unit=17, file='flowstat.dat')
      open ( unit=13, file='transfer.dat')
      open ( unit=14, file='vortpdf.dat')
      endif
c
      ifield = nhalt/itvort
      do iout = 0,ifield
      idump = 400 + iout
      i0d=idump/100
      i1d=(idump-i0d*100)/10
      i2d=idump-i0d*100-i1d*10
      fconc = './HI2/fort'//char(i0d+48)//char(i1d+48)
     &    //char(i2d+48)//'.dat'
      open( unit=idump,file=fconc,status='unknown')
      write(idump,*) 'variables="x","y","z","u","v","w","ip"'

c
      idump = 600 + iout

      i0d=idump/100
      i1d=(idump-i0d*100)/10
      i2d=idump-i0d*100-i1d*10      
      fvort = './HI2/fort'//char(i0d+48)//char(i1d+48)
     & //char(i2d+48)//'.dat'

      open( unit=idump,file=fvort,status='unknown')
      write(idump,*) 'variables="x","y","z","<greek>w</greek>"'
      write(idump,*) 'zone i=',n1+1,'j=',n2+1,'k=',1
      enddo

      icdump=400
      itdump=600
c     ifdump=80
c
c    Set iseeds: iseedp for seeding particles in a new particle
c   system and iseedf for the random forcing.
c
      read (15,*) readtemp, iseedp,iseedf

CRANDOM
      if(newpart)then
      iyp = 0
      iseedp = - iseedp
      endif

      if(newflo)then
      iyf = 0
      iseedf = - iseedf
      endif

c
c    Start the output file for job
c    -----------------------------
c
      write (16,*) 'Nproc= ',nnodes
      write (16,1001)
      write (16,1010)  n1,n2,n3
      write (16,1015)  rnu,dt,var,tf
      write (16,1020)  ttime
      write (16,1025)  newflo,newpart,forced
      write (16,1030)  nhalt,nlong,nshort
      write (16,1035)  iseedp,iseedf

c
c    Set up particle parameters
c    ------------------------
c   Particles are grouped into npset equal lots of numb
c   particles. The total should be npart. Each set will have
c   a different tau and wp value.
c   Assume for now that npset is .le. 8 .
c

      read (15,*) readtemp, etk0,tauk0,vk0
      write (16,*)'etk0,tauk0,vk0=',etk0,tauk0,vk0
      read (15,*) readtemp, les
      write(16,*) "les=",les
      read (15,*) readtemp, model_sgs
      write (16,*) readtemp, model_sgs
      read (15,*) readtemp, Ck
      write(16,*) "Kolmpgorov constant=",Ck   
      read (15,*) readtemp, fdns
      write(*,*) "fdns=",fdns
      write(16,*) "fdns=",fdns
      read (15,*) readtemp, cutoff_num
      write(16,*) "cutoff_num=",cutoff_num



c radius normalized by kolmogorov length scale
      read(15,*)readtemp,arad1,arad2
      write(16,*)'arad1,arad2= ',arad1,arad2
      rad1 = arad1*etk0
      rad2 = arad2*etk0

      read (15,*) readtemp, npset
      numb=npart/npset
      write (17,1080) npart,npset,numb
      if ( (numb*npset).ne.npart)  then
        write (16,1100)
        stop
      end if
c
c    Open data files for short particle statistics
c
      if (npset.ge.1) open( unit=21,file='part1a.dat',status='unknown')
      if (npset.ge.1) open( unit=31,file='part1b.dat',status='unknown')
      if (npset.ge.2) open( unit=22,file='part2a.dat',status='unknown')
      if (npset.ge.2) open( unit=32,file='part2b.dat',status='unknown')
c
      write (16,1110)
c
      do  iset=1,npset
        read (15,*) readtemp, tauset,wpset
        write (16,*) tauset,wpset

c$omp parallel do private(ip,k)
        do  jp=1,numb
          ip=jp+(iset-1)*numb
          tau(ip)=tauset*tauk0
          wp(ip)=wpset*vk0
        end do
c$omp end parallel do
      end do

      if(fdns) then
        allocate ( tur_filtered(n1pp,n2,n3) )
        allocate ( tvr_filtered(n1pp,n2,n3) )
        allocate ( twr_filtered(n1pp,n2,n3) )
        allocate ( sxr_filtered(n1pp,n2,n3) )
        allocate ( syr_filtered(n1pp,n2,n3) )
        allocate ( szr_filtered(n1pp,n2,n3) )
      endif


        if(wp(1).eq.wp(npart)) then
         w0diff=wp(1)
        else
	 w0diff=abs(wp(1)-wp(npart))
         if( (w0diff.gt.wp(1)).or.(w0diff.gt.wp(npart)) ) then
         w0diff=wp(1)
         if(w0diff.gt.wp(npart))w0diff=wp(npart)
         endif
        endif
c
c
c    If this is a restart, load data cur,..,szp
c    ------------------------------------------
c
      if ( newflo ) then
      continue
      else
        open ( unit=11, file='./HI2/restart.flo',
     &  form='unformatted',status='old' )
        read (11) cur, cvr, cwr
        read (11) sxp, syp, szp
CRANDOM
        read (11) iseedf,iyf,ivf
        if ( forced ) read (11) a1r,a2r,a3r
        if ( forced ) read (11) b1r,b2r,b3r
        close(11)
        if(model_sgs.eq.2) then
          open(unit=11, file='./HI2/les_sgs_visc_old.dat',
     &  form='unformatted',status='old' )
          do i=1,itrunc
          read(11,*) sgs_visck_old(i)
          enddo
          close(11)
         endif
       endif

c
c parameters for collision detection
c
       wcd=hx/float(ncd)
       wpvcd=hx/float(npvcd)
       rcoll11= rad1 + rad1
       rcoll12= rad1 + rad2
       rcoll22= rad2 + rad2
       write(16,*)'wcd, rcoll11, rcoll12, rcoll22 = '
     1 , wcd, rcoll11, rcoll12, rcoll22

        do ir=0,179
        nir11(ir)=0
        nir12(ir)=0
        nir22(ir)=0
        wrbin11(ir)=0.0
        wrbin12(ir)=0.0
        wrbin22(ir)=0.0
        enddo

       rcoll11a = 0.975*rcoll11
       rcoll11b = 1.025*rcoll11
       rcoll12a = 0.975*rcoll12
       rcoll12b = 1.025*rcoll12
       rcoll22a = 0.975*rcoll22
       rcoll22b = 1.025*rcoll22

      istep=0

      if ( newpart ) then

      do  iset=1,npset
cc$omp parallel do private(ip,ar1,ar2,ar3)
cc$omp+ shared(iseedp,iyp,ivp)
        do  jp=1,numb
        ip = jp + (iset-1)*numb
CRANDOM
cc$omp critical
          ar1=ranpp(iseedp,iyp,ivp)
          ar2=ranpp(iseedp,iyp,ivp)
          ar3=ranpp(iseedp,iyp,ivp)
cc$omp end critical

          yp(ip,1)= ar1*hx
          yp(ip,2)= ar2*hy
          yp(ip,3)= ar3*hz
        end do
cc$omp end parallel do
        end do
c
c   Overlapping Detection and Correction
c   ***********************************
c
c   Indentifying the particle location relative to
c   the ovarelapping detection grid
c
1981 	continue
c
c  First ensure that particles remain in box, use periodicity
c
c$omp parallel do
      do  ip=1,npart
        if (  yp(ip,1).ge.hx ) yp(ip,1)=yp(ip,1)-hx
        if (  yp(ip,1).lt.0.0) yp(ip,1)=yp(ip,1)+hx
        if (  yp(ip,2).ge.hy ) yp(ip,2)=yp(ip,2)-hy
        if (  yp(ip,2).lt.0.0) yp(ip,2)=yp(ip,2)+hy
        if (  yp(ip,3).ge.hz ) yp(ip,3)=yp(ip,3)-hz
        if (  yp(ip,3).lt.0.0) yp(ip,3)=yp(ip,3)+hz
      end do
c$omp end parallel do
c
	if(OVERL)goto 809
c	write(*,*) 'Before Overl Det'

        do i=1,ncd
        do j=1,ncd
        do k=1,ncd
        head(i,j,k)=0
        enddo
        enddo
        enddo

      do ip=1,npart
      ix=1+int(yp(ip,1)/wcd)
      iy=1+int(yp(ip,2)/wcd)
      iz=1+int(yp(ip,3)/wcd)
c
      if(ix.gt.ncd ) then
      ix = ncd
      endif

      if(iy.gt.ncd)then
      iy = ncd
      endif

      if(iz.gt.ncd)then
      iz = ncd
      endif
c
        list(ip)=head(ix,iy,iz)
        head(ix,iy,iz)=ip
      end do

	noverl=0
c
c$omp parallel do private(i,j,k)
c$omp+ reduction(+:noverl)
c$omp+ private(ix,iy,iz)
c$omp+ private(ix2,sx,iy2,sy,iz2,sz)
c$omp+ private(ip2,ypb1,ypb2,ypb3)
c$omp+ private(dnij,ar1,ar2,ar3)
c$omp+ private(rcoll)
c$omp+ shared(iseedp,iyp,ivp)
        DO IP1=1, NPART
         ix=1+int(yp(ip1,1)/wcd)
         iy=1+int(yp(ip1,2)/wcd)
         iz=1+int(yp(ip1,3)/wcd)
c
      if(ix.gt.ncd ) then
      write(103,*)ix,iy,iz,ip1
      write(103,*)ttime,yp(ip1,1),yp(ip1,2),
     1   yp(ip1,3),wcd,hx,hy,hz
      ix = ncd
      endif

      if(iy.gt.ncd)then
      write(103,*)ix,iy,iz,ip1
      write(103,*)ttime,yp(ip1,1),yp(ip1,2),
     1   yp(ip1,3),wcd,hx,hy,hz
      iy = ncd
      endif

      if(iz.gt.ncd)then
      write(103,*)ix,iy,iz,ip1
      write(103,*)ttime,yp(ip1,1),yp(ip1,2),
     1   yp(ip1,3),wcd,hx,hy,hz
      iz = ncd
      endif
c
        do i=-1,1
        do j=-1,1
        do k=-1,1
        ix2=ix+i
        sx=0.0
        if(ix2.lt.1) then
        ix2=ncd
        sx=-1.
        endif
        if(ix2.gt.ncd)then
        ix2=1
        sx=1.
        endif
c
        iy2=iy+j
        sy=0.0
        if(iy2.lt.1) then
        iy2=ncd
        sy=-1.
        endif
        if(iy2.gt.ncd)then
        iy2=1
        sy=1.
        endif
c
        iz2=iz+k
        sz=0.0
        if(iz2.lt.1) then
        iz2=ncd
        sz=-1.
        endif
        if(iz2.gt.ncd)then
        iz2=1
        sz=1.
        endif
        ip2=head(ix2,iy2,iz2)
c
c       ip2 has to be always greater than ip1 to avoid double counting
c
        if(ip2.eq.0)goto 2001
        if(ip2.le.ip1)goto 911
c
1991     if(ip1.le.numb) then
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

       ypb1=yp(ip2,1)+sx*hx
       ypb2=yp(ip2,2)+sy*hy
       ypb3=yp(ip2,3)+sz*hz

       dnij=sqrt((yp(ip1,1)-ypb1)**2
     1 +(yp(ip1,2)-ypb2)**2+(yp(ip1,3)-ypb3)**2)

        if(dnij.le.rcoll) then
           noverl=noverl+1
           ar1=ranpp(iseedp,iyp,ivp)
           ar2=ranpp(iseedp,iyp,ivp)
           ar3=ranpp(iseedp,iyp,ivp)
           yp(ip2,1)= ar1*hx
           yp(ip2,2)= ar2*hy
           yp(ip2,3)= ar3*hz
        endif

911    ip2=list(ip2)
        if(ip2.eq.0)goto 2001
        if(ip2.le.ip1)goto 911
        goto 1991
2001     continue
      enddo
      enddo
      enddo
      ENDDO

c	write(*,*) 'noverl= ',noverl

	if(noverl.ne.0) goto 1981
c
c  END OF DETECTION OF OVERLAPPING
c

c	write(*,*) 'After Overl Det'
809	continue
c
      else
        open ( unit=12, file='./HI2/restart.part', 
     &  form='unformatted', status='old')
        read (12) yp,vp,fvp,fdp
CRANDOM
        read(12)iseedp,iyp,ivp
        close(12)

      endif

      do  iset=1,npset
c$omp parallel do private(ip)
        do  jp=1,numb
        ip = jp + (iset-1)*numb
          yp0(ip,1)=yp(ip,1)
          yp0(ip,2)=yp(ip,2)
          yp0(ip,3)=yp(ip,3)
          ypp(ip,1)=yp(ip,1)
          ypp(ip,2)=yp(ip,2)
          ypp(ip,3)=yp(ip,3)
        end do
c$omp end parallel do
        end do
c$omp parallel do private(j,i)
        do k=1, 128
         do j=1, 128
          do i=1,128
          con1(i,j,k)=0
          con2(i,j,k)=0
          enddo
         enddo
        enddo
c$omp end parallel do
       do i=0, 30
        pc1(i)=0.0
        pc2(i)=0.0
        pcr(i)=0.0
       enddo

c
c
c   Start of time step loop
c   -----------------------
c
c     Note at this point the flow data at time level n is
c   available and can be used to advance the particles to time
c   level n+1. We need to interpolate the flow data, so
c   some initial setup for this is done and can be applied
c   as needed below. 
c     Once the particles are moved the flow is advanced.
c
c     Note the particles are reintroduced into the box by
c   periodicity if they move out of it. So we can assume that
c   0<= y1 <= hx, etc.
c
c
      itocol11=0
      itocol12=0
      itocol22=0
      icc1=0
      icc2=0
      icc3=0
      nR11=0
      nR12=0
      nR22=0
C
      open(8, file='sgs_vis.dat')
 100  continue
      write(*,*)'begining istep=',istep  !, fdate()

        do ir=0,179
        nir11(ir)=0
        nir12(ir)=0
        nir22(ir)=0
        wrbin11(ir)=0.0
        wrbin12(ir)=0.0
        wrbin22(ir)=0.0
        enddo


c       write(*,*)'istep=',istep
c
c-- write particle concentration field
c     if ( mod(istep,itvort).eq.0) then

cc$omp parallel do private(j,i)
c     do  k=1,n3
c       do j=1,n2
c         do i=1,n1
c           tur( i,j,k) = 0.0
c         end do
c       end do
c     end do
cc$omp end parallel do

cc$omp parallel do private(jx,ix1,jy,ix2,jz,ix3)
c       do  ip=1,numb
c       jx=nint(yp(Ip,1)/DXYZ(1))
c       iX1=( mod(jx,n1) + n1*(1-isign(1,jx))/2) + 1
c       jy=nint(yp(Ip,2)/DXYZ(2))
c       IX2=( mod(jy,n2) + n2*(1-isign(1,jy))/2) + 1
c       jz=nint(yp(Ip,3)/DXYZ(3))
c       IX3=( mod(jz,n3) + n3*(1-isign(1,jz))/2) + 1
c       TUR(IX1,IX2,IX3)=TUR(IX1,IX2,IX3)+1.
c       end do
cc$omp end parallel do
C
c     write(icdump,707)(((tur(i,j,k),i=1,n1),j=1,n2),k=1,n3)
c     close(icdump)
c     icdump=icdump+1
c     endif

CTiming

       T1 = SECNDS(0.0)
       TTT1 = T1

	if(NOTURB) goto 800 

c   Set up interpolation factors bg for the particles and
c   locate grid points. lhnode is an index for the nearest
c   grid point coordinate to the left of the particle. The
c   value of lhnode is between 0 and n-1.
c
      do  ix=1,3
c$omp parallel do private(node,z,z2,z3,z4,z5)
        do ip=1,npart
          node=int( yp(ip,ix)/dxyz(ix) )
          if(yp(ip,ix).lt.0.0)node=node-1
          lhnode(ip,ix)=node
          z= ( yp(ip,ix)-float(node)*dxyz(ix) )/dxyz(ix)
          z2=z*z
          z3=z2*z
          z4=z3*z
          z5=z4*z
          bg(ip,1,ix)=( 6.0*z - 5.0*z2 - 5.0*z3 + 5.0*z4 - z5 )/120.0
          bg(ip,2,ix)=( - 12.0*z + 16.0*z2 - z3 - 4.0*z4 + z5 )/24.0
          bg(ip,3,ix)=( 12.0 - 4.0*z - 15.0*z2 + 5.0*z3 + 3.0*z4 
     &                 - z5 )/12.0
          bg(ip,4,ix)=( 12.0*z + 8.0*z2 - 7.0*z3 - 2.0*z4 + z5 )/12.0
          bg(ip,5,ix)=( - 6.0*z - z2 + 7.0*z3 + z4 - z5 )/24.0
          bg(ip,6,ix)=( 4.0*z - 5.0*z3 + z5 )/120.0
        end do
c$omp end parallel do
      end do
CTiming

800	continue

       T2 = SECNDS(0.0)
       if(T2.ge.T1) then
       Tpart1 = T2-T1
       else
       Tpart1=(24*3600-T1)+T2
       endif

CTiming
       T1 = SECNDS(0.0)
        if(NOTURB) goto 801
c
c   Copy fourier coefficients and obtain velocity
c   in physical space. Store in tur, tvr, twr.
c
c$omp parallel do private(j,i)
      do  k=1,n3
        do j=1,n2
          do i=1,n1pp
            tur( i,j,k) = cur( i,j,k)
            tvr( i,j,k) = cvr( i,j,k)
            twr( i,j,k) = cwr( i,j,k)
          end do
        end do
      end do
c$omp end parallel do
c      
       if(fdns) then
c$omp parallel do private(j,i,coeff)
      do  k=1,n3
        do j=1,n2
          do i=1,n1pp
            coeff=sqrt(fsq1(i)+fsq2(j)+fsq3(k))
            tur_filtered( i,j,k) = cur( i,j,k)
            tvr_filtered( i,j,k) = cvr( i,j,k)
            twr_filtered( i,j,k) = cwr( i,j,k)
            if(coeff.ge.fkstar*cutoff_num) then            
              tur_filtered( i,j,k) = 0.0
              tvr_filtered( i,j,k) = 0.0
              twr_filtered( i,j,k) = 0.0
            endif
          end do
        end do
      end do
c$omp end parallel do
        print*,"before transfer complex to real space1"
        call fft3DCR(tur_filtered,plan1DB,plan2DCR,nnodes)
        call fft3DCR(tvr_filtered,plan1DB,plan2DCR,nnodes)
        call fft3DCR(twr_filtered,plan1DB,plan2DCR,nnodes)
        print*,"after transfer complex to real space1"
       endif

        print*,"before transfer complex to real space1"       
        call fft3DCR(tur,plan1DB,plan2DCR,nnodes)
        call fft3DCR(tvr,plan1DB,plan2DCR,nnodes)
        call fft3DCR(twr,plan1DB,plan2DCR,nnodes)

        print*,"after transfer complex to real space1"

c     if ( mod(istep,itvort).eq.0 .and. istep.gt.0) then
c     write(ifdump)cur,cvr,cwr
c     close(ifdump)
c     ifdump=ifdump+1
c     endif

c
c   Calculate vorticity in fourier space and store in sxr,syr,szr
c   -------------------------------------------------------------
c   wz  = dv/dx - du/dy
c  
c   wy  = du/dz - dw/dx
c  
c   wx  = dw/dy - dv/dz
c  
c
c$omp parallel do private(j,i,ip,zur,zvr,zwr,zui,zvi,zwi)
      do  k=1,n3
        do  i=1,n1pp,2
          ip=i+1
          do  j=1,n2 
            zur=cur(i,j,k)
            zvr=cvr(i,j,k)
            zwr=cwr(i,j,k)
            zui=cur(ip,j,k)
            zvi=cvr(ip,j,k)
            zwi=cwr(ip,j,k)
c
            szr(i,j,k)  = -fac1(i)*zvi+fac2(j)*zui
            szr(ip,j,k) = +fac1(i)*zvr-fac2(j)*zur
c
            syr(i,j,k)  = -fac3(k)*zui+fac1(i)*zwi
            syr(ip,j,k) = +fac3(k)*zur-fac1(i)*zwr
c
            sxr(i,j,k)  = -fac2(j)*zwi+fac3(k)*zvi
            sxr(ip,j,k) = +fac2(j)*zwr-fac3(k)*zvr
          end do
        end do
      end do
c$omp end parallel do


       if(fdns.and.mod(istep,itvort).eq.0) then
c$omp parallel do private(j,i,coeff)
      do  k=1,n3
        do j=1,n2
          do i=1,n1pp
            coeff=sqrt(fsq1(i)+fsq2(j)+fsq3(k))
            sxr_filtered( i,j,k) = sxr( i,j,k)
            syr_filtered( i,j,k) = syr( i,j,k)
            szr_filtered( i,j,k) = szr( i,j,k)
            if(coeff.ge.fkstar*cutoff_num) then
              sxr_filtered( i,j,k) = 0.0
              syr_filtered( i,j,k) = 0.0
              szr_filtered( i,j,k) = 0.0
            endif
          end do
        end do
      end do
c$omp end parallel do

      print*,"before Transform vorticity to physical space"
      call fft3DCR(sxr_filtered,plan1DB,plan2DCR,nnodes)
      call fft3DCR(syr_filtered,plan1DB,plan2DCR,nnodes)
      call fft3DCR(szr_filtered,plan1DB,plan2DCR,nnodes)
      print*,"after Transform vorticity to physical space"

       endif
c
c
c   Transform vorticity to physical space
c   -------------------------------------
c
      print*,"before Transform vorticity to physical space"
      call fft3DCR(sxr,plan1DB,plan2DCR,nnodes)
      call fft3DCR(syr,plan1DB,plan2DCR,nnodes)
      call fft3DCR(szr,plan1DB,plan2DCR,nnodes)
      print*,"after Transform vorticity to physical space"

       if(les) then
        do i=1,n1pp
           phi1(i)=0.0
        enddo
cc reduction .... $omp parallel do private(j,i,temp,itemp,xmult,usqr)

        do k=1,n3
           do j=1,n2
              temp=fsq3(k)+fsq2(j)
              do i=1,n1pp
                 itemp=int(sqrt(temp+fsq1(i))+0.5)
                 if ( itemp .gt. 0 ) then
                 if ( i .le. 2 ) then
                 xmult=0.5
                 else
                 xmult=1.
                 endif
                 usqr=xmult*(cur(i,j,k)**2+cvr(i,j,k)**2
     &                      +cwr(i,j,k)**2)
                 phi1(itemp)=phi1(itemp)+usqr
                 endif
              enddo
           enddo
        enddo

cc$omp end parallel do

       itrunc=int(fkstar)

       do i=1,n2h
         sgs_visck(i)=0.0
       enddo

       if(model_sgs.eq.1) then
         sgs_visc=sqrt(abs(phi1(itrunc))/float(itrunc))
         do i=1, itrunc
            visc_temp2=1+34.467*exp(-3.03*float(itrunc)/float(i))
            sgs_visck(i)=0.441*Ck**(-3./2.)*sgs_visc* visc_temp2
         enddo
       else
         totaldiss=0.0
         do i=1, itrunc
           totaldiss=totaldiss+2*(rnu+sgs_visck_old(i))*i**2*phi1(i)
         enddo                      
         eta_les=(rnu**3/totaldiss)**0.25
	 sumk=(itrunc-3)-itrunc/2+1

	 sumlogei=0.0
	 sumlogki=0.0
	 sumlogki2=0.0
	 sumlogeilogki=0.0

	 do i=itrunc/2, (itrunc-3)
           sumlogei=sumlogei+log(phi1(i))
	   sumlogki=sumlogki+log(float(i))
	   sumlogki2=sumlogki2+(log(float(i)))**2
           sumlogeilogki=sumlogeilogki+log(phi1(i))*log(float(i))
	 enddo

         temp1=sumlogei*sumlogki/sumk-sumlogeilogki
	 temp2=sumlogki2-sumlogki**2/sumk
	 rm=temp1/temp2

         C_kolm=0.0

         do i=itrunc/2, (itrunc-3)
           C_kolm=C_kolm+phi1(i)*totaldiss**(-2./3.)*i**rm
         enddo
         C_kolm=C_kolm/sumk

!         write(8,*)" C_kolm=",C_kolm

         do i=itrunc+1,n1pp
           f_eta_k=exp(-5.2*(((i*eta_les)**4+0.3**4)**0.25-0.3))
           phi1(i)=C_kolm*totaldiss**(2./3.)*i**(-rm)*f_eta_k
         enddo

        if(istep.eq.int(istep/10)*10) then
          do i=1, n1pp
            f_eta_k=exp(-5.2*(((i*eta_les)**4+0.3**4)**0.25-0.3))
            write(8,101) i, phi1(i),totaldiss,C_kolm,eta_les,f_eta_k
          enddo
        endif  

        sumk=0.0
	do i=1, n1pp
         sumk=sumk+i**2*phi1(i)
         eddy_time(i)=0.218*C_kolm**1.5*sumk**0.5
	enddo  

!	open(8,file="result.dat")
	      
	do k=1, itrunc
 	  k2=k**2
        do ip=1, n1pp
          ip2=ip**2
	    do iq=1, n1pp
            iq2=iq**2   
            if(ip.gt.itrunc.or.iq.gt.itrunc) then
	        costheta=float(ip2+iq2-k2)/float(2*ip*iq)
	        if(abs(costheta).le.1.0) then
                  x=float(ip2+iq2-k2)/float(2*ip*iq)
  		  y=float(k2+iq2-ip2)/float(2*k*iq)
	          z=float(k2+ip2-iq2)/float(2*k*ip)
                  thetakpq=1./(eddy_time(k)+eddy_time(ip)+eddy_time(iq))
	          btemp=(x*y+z**3)/float(iq)
	          temp2=phi1(iq)*(k**2*phi1(ip)-ip**2*phi1(k))
	          sgs_visck(k)=sgs_visck(k)+thetakpq*btemp*temp2
		endif
	      endif
	    enddo
	  enddo
        sgs_visck(k)=-sgs_visck(k)/(2*k2*phi1(k))
        sgs_visck_old(k)=sgs_visck(k)
!         write(8,*) "i=",k,"sgs_visck(",k,")=",sgs_visck(k)
      enddo
      endif

!       if(rm.gt.3) rm=3.0

!        do i=1, itrunc
!          temp_k=float(i)/float(itrunc) 
!          visc_temp1=0.3055*((5-rm)/(rm+1))*
!     1            sqrt(3-rm)*Ck**(-3./2.)

!          visc_temp2=(0.23656-0.18146*temp_k+0.89461*temp_k**4)*
!     1         exp(-1.88539*temp_k+1.24564*temp_k**2)/temp_CK 

!          C_plus=(0.24+2.28805*temp_k-14.862*temp_k**2+69.0*temp_k**3)*
!     1           exp(-5.7726*temp_k+1.13543*temp_k**2)

!          C_ratio_backscatter=C_plus*temp_k**2*phi1(itrunc)/phi1(i)
                              
!         visc_temp2=0.01+34.467*exp(-3.03*float(itrunc)/float(i))
!          eddy_visc_forward=visc_temp1*visc_temp2*sgs_visc
          
!          sgs_visck(i)=(1.-C_ratio_backscatter)*eddy_visc_forward
!        enddo
!             print*,"dynamic"
!          else
!             print*,"old"
!              do i=1, itrunc
!               visc_temp2=0.05+34.467*exp(-3.03*float(itrunc)/float(i))
!               sgs_visck(i)=0.441*Ck**(-3./2.)*sgs_visc* visc_temp2
!          temp_k=float(i)/float(itrunc)
!          visc_temp2=(0.23656-0.18146*temp_k+0.89461*temp_k**4)*
!     1         exp(-1.88539*temp_k+1.24564*temp_k**2)
!          sgs_visck(i)=visc_temp1*visc_temp2*sgs_visc
!              enddo
!          endif   
         endif
 
c-- save vorticity field for visualization
      if(mod(istep,itvort).eq.0) then

      vrms=0.0

c$omp parallel do 
c$omp+ private(j,i)
        do k=1,n3
           do j=1,n2
              do i=1,n1
         if(fdns) then
           tempa(i,j,k)=sqrt(sxr_filtered(i,j,k)**2+syr_filtered(i,j,k)**2+
     &                      szr_filtered(i,j,k)**2)
           vrms=vrms+tempa(i,j,k)**2         
         else
           tempa(i,j,k)=sqrt(sxr(i,j,k)**2+syr(i,j,k)**2+szr(i,j,k)**2)
           vrms=vrms+tempa(i,j,k)**2
         endif      
             enddo
            enddo
        enddo
c$omp end parallel do

      vrms=sqrt(vrms/3./real(n1*n2*n3))

c$omp parallel do private(j,i)
      do k=1,n3
      do j=1,n2
      do i=1,n1               
        tempa(i,j,k)=tempa(i,j,k)/vrms
      if(tempa(i,j,k).gt.999.0)tempa(i,j,k)=999.0
      enddo
      enddo
      enddo
c$omp end parallel do

c      write(itdump,707)(((tempa(i,j,k),i=1,n1),j=1,n2),k=1,n3)
c      close(itdump)
c      itdump=itdump+1
      k=n2/2+1
      do i=1,n1+1
        do j=1,n2+1
          i_temp=i
          j_temp=j
          if(i.eq.n1+1) i_temp=1
          if(j.eq.n2+1) j_temp=1
        write(itdump,707) (i-1)*dx,(j-1)*dx,0.0,tempa(i_temp,j_temp,k)
        enddo
      enddo
      close(itdump)
      itdump=itdump+1
      endif
c
c   Short flow statistics if needed.
c   -------------------------------
c    Do this in physical space for mean square velocity
c    and vorticity. KE dissipation rate probe(7) is found
c    from the vorticity using the periodic bc's/homogeneity.
c    The rate of working by the random forcing is probe(8).
c

      if(pflow)then
      if ( mod(istep,nshort).eq.0 .and. istep.gt.0) then
c     compute energy and dissipation spectrum
c$omp parallel do
        do i=1,n1pp
           phi1(i)=0.0
           phi2(i)=0.0
           phi3(i)=0.0
           phi4(i)=0.0
           phi5(i)=0.0
        enddo
c$omp end parallel do
C
C
cc reduction .... $omp parallel do private(j,i,temp,itemp,xmult,usqr)
        do k=1,n3
           do j=1,n2
              temp=fsq3(k)+fsq2(j)
              do i=1,n1pp
                 itemp=int(sqrt(temp+fsq1(i))+0.5)
                 if ( itemp .gt. 0 ) then
                 if ( i .le. 2 ) then
                 xmult=0.5
                 else
                 xmult=1.
                 endif
                 usqr=xmult*(cur(i,j,k)**2+cvr(i,j,k)**2
     &                                       +cwr(i,j,k)**2)
                 phi1(itemp)=phi1(itemp)+usqr
           phi4(itemp)=phi4(itemp)+usqr*(temp+fsq1(i))
C
           phi3(itemp)=phi3(itemp)+xmult*(cur(i,j,k)**2*fsq1(i)
     &              +cvr(i,j,k)**2*fsq2(j)+cwr(i,j,k)**2*fsq3(k))
           phi2(itemp)=phi2(itemp)+2.*rnu*usqr*(temp+fsq1(i))
                 endif
              enddo
           enddo
        enddo
cc$omp end parallel do
C
C
      itrunc=int(fkstar)
      if(les) then       
        do i=1,itrunc          
          write(19,102) i,phi1(i),phi2(i),
     1                  (phi2(i)/rnu)*(rnu+sgs_visck(i)),sgs_visck(i)
     1                  ,rm    
        enddo       
      else
        do i=1,itrunc
          write(19,102) i,phi1(i),phi2(i)
        enddo
      endif
C   
        do  ic=1,8
          probe(ic)=0.0
        end do

cc reduction .. . $omp parallel do private(j,i,ip)
        do  k=1,n3
          do  j=1,n2
            do  i=1,n1
              probe(1)=probe(1)+tur(i,j,k)*tur(i,j,k)
              probe(2)=probe(2)+tvr(i,j,k)*tvr(i,j,k)
              probe(3)=probe(3)+twr(i,j,k)*twr(i,j,k)
              probe(4)=probe(4)+sxr(i,j,k)*sxr(i,j,k)
              probe(5)=probe(5)+syr(i,j,k)*syr(i,j,k)
              probe(6)=probe(6)+szr(i,j,k)*szr(i,j,k)
            end do
          end do
        end do
cc$omp end parallel do

        do  ic=1,6
        probe(ic)=probe(ic)/float(n1*n2*n3)
        end do
        probe(7)=0.0
        tmse=0.0
        do i=1,itrunc
        probe(7)=probe(7)+phi2(i)
        tmse=tmse+phi3(i)
        enddo
        tmse=tmse*2.0/3.0
c
        do  ic=1,6
          efx(ic)=0.0
        end do
        do  k=1,5
          kk=k
          if(k.gt.3) kk=k+n3-5
          do  j=1,5
            jj=j
            if(j.gt.3) jj=j+n2-5
            do  i=1,6
              efx(i)=a1r(i,j,k)*cur(i,jj,kk)+
     &               a2r(i,j,k)*cvr(i,jj,kk)+
     &               a3r(i,j,k)*cwr(i,jj,kk)+efx(i)
            end do
          end do
        end do
        efsum= (efx(1)+efx(2))+2.0*(efx(3)+efx(4)+efx(5)+
     &          efx(6))
        probe(8)=efsum
c
        write(18,101) istep,ttime,probe

        if ( mod(istep,nlong).eq.0) then
        ediss=probe(7)
        eta=(rnu**3/ediss)**(1./(6.-2.0))
        tk = (rnu/ediss)**(1./(3.-1.0))
        vk = eta/tk
c       write(20,*)'ediss=  ',ediss
c       write(20,*)'eta=    ',eta
c       write(20,*)'vk=     ',vk
c       write(20,*)'tk=     ',tk
        uprime=sqrt( (probe(1)+probe(2)+probe(3))/3.0)
        tmse=sqrt(uprime*uprime/tmse)
        ree=uprime*tmse**(2.-1.0)/rnu
        xl=uprime**3/ediss
        et=xl/uprime
        resol=fkstar*eta

        temp=0.0
c$omp parallel do reduction(MAX:temp)
c$omp+ private(j,i,outa)
        do k=1,n3
           do j=1,n2
              do i=1,n1
               outa=abs(tur(i,j,k))+abs(tvr(i,j,k))
     &              +abs(twr(i,j,k))
               temp = max(temp,outa)
          enddo
         enddo
        enddo
c$omp end parallel do

C----  Computing pdf distribution for vorticity magnitude
        vmax=0.0

c$omp parallel do reduction(MAX:vmax)
c$omp+ private(j,i,sss)
        do k=1,n3
           do j=1,n2
              do i=1,n1
           sss=sqrt(sxr(i,j,k)**2+syr(i,j,k)**2+szr(i,j,k)**2)
         vmax = max(vmax,sss)
          enddo
         enddo
        enddo
c$omp end parallel do

        delta=vmax*1.00001/200.
        do i=1,200
        vdist(i)=0.0
        enddo
        do k=1,n3
           do j=1,n2
              do i=1,n1
         sss=sqrt(sxr(i,j,k)**2+syr(i,j,k)**2+szr(i,j,k)**2)
         isss=int(sss/delta+1.0)
         vdist(isss)=vdist(isss)+1.0
         enddo
         enddo
        enddo
        do i=1,200
        xxx=delta*(real(i)-0.5)
        vdist(i)=vdist(i)/real(n1*n2*n3)
        write(14,103)xxx,vdist(i)
        enddo
C
        CFL=temp*dt/dx
        cfl2=rnu*dt/dx**2
c       write(20,*)'  tmse=  ',tmse
c       write(20,*)'  ree=  ',ree
c       write(20,*)'  L_f=  ',xl
c       write(20,*)'  T_e=  ',et
c       write(20,*)'uprime=  ',uprime
c       write(20,*)'  resol=  ',resol
c       write(20,*)'  CFL=  ',CFL
c       write(20,*)'  CFL2=  ',cfl2


        write(17,105)ttime,uprime,xl,ediss,et,tmse,
     .       eta,vk,tk,ree,resol,CFL,cfl2,sgs_visc
        write(*,*)istep,resol,CFL
        endif
c
      end if
      end if

CTiming

801	continue

       T2 = SECNDS(0.0)
       if(T2.ge.T1) then
       Tflow1 = T2-T1
       else
       Tflow1=(24*3600-T1)+T2
       endif
CTiming
        T1 = SECNDS(0.0)
c
c
c   Interpolate flow velocity ( and vorticity if needed)
c   ---------------------------------------------------
c
	if (NOTURB) then
c$omp parallel do
        do  ip=1,npart
          up(ip,1)=0.0
          up(ip,2)=0.0
          up(ip,3)=0.0
        end do
c$omp end parallel do
	else
         if(.not.fdns) then
           icpt=1
           call value(tur,up,icpt,lhnode,bg)
           if ( mod(istep,nshort).eq.0)
     1     call value(sxr,vortp,icpt,lhnode,bg)
           jcpt=2
           call value(tvr,up,jcpt,lhnode,bg)
           if ( mod(istep,nshort).eq.0)
     1     call value(syr,vortp,jcpt,lhnode,bg)
           kcpt=3
           call value(twr,up,kcpt,lhnode,bg)
           if ( mod(istep,nshort).eq.0)
     1     call value(szr,vortp,kcpt,lhnode,bg)
         else
           icpt=1
           call value(tur_filtered,up,icpt,lhnode,bg)
           if ( mod(istep,nshort).eq.0)
     1     call value(sxr_filtered,vortp,icpt,lhnode,bg)
           jcpt=2
           call value(tvr_filtered,up,jcpt,lhnode,bg)
           if ( mod(istep,nshort).eq.0)
     1     call value(syr_filtered,vortp,jcpt,lhnode,bg)
           kcpt=3
           call value(twr_filtered,up,kcpt,lhnode,bg)
           if ( mod(istep,nshort).eq.0)
     1     call value(szr_filtered,vortp,kcpt,lhnode,bg)        
         endif
	endif
c
c   For a new run set the initial particle velocity
c
      if ( newpart ) then
c$omp parallel do
        do  ip=1,npart
!          vp(ip,1)=wp(ip)+up(ip,1)
!          vp(ip,2)=up(ip,2)
!          vp(ip,3)=up(ip,3)

          vp(ip,1)=wp(ip)
          vp(ip,2)=0.0
          vp(ip,3)=0.0

        end do
c$omp end parallel do
      end if
CTiming
       T2 = SECNDS(0.0)
       if(T2.ge.T1) then
       Tpart2 = T2-T1
       else
       Tpart2=(24*3600-T1)+T2
       endif

CTiming

       T1 = SECNDS(0.0)
c
c INCLUDE HERE THE CHANGE ON vp for ALREADY COLLIDED PARTICLES!!!!

	if(OVERL) goto 807
       if((istep.gt.0).and.(newcoltot.gt.0)) then
           do ipos=1,ncollpart
           ip=removpart(ipos)
           vp(ip,1)=wp(ip)+up(ip,1)
           vp(ip,2)=up(ip,2)
           vp(ip,3)=up(ip,3)
           enddo
       endif
807	continue
c
c to get the relative velocity including hydrodynamic interactions

       if (HYDINT) then
         if ( newpart ) then
c$omp parallel do
          do  ip=1,npart
          pertvel(ip,1)=0.0
          pertvel(ip,2)=0.0
          pertvel(ip,3)=0.0
          enddo
c$omp end parallel do
         endif
         call PERTURBVELOC(vp,yp,up,pertvel)
       else
c$omp parallel do
         do  ip=1,npart
          pertvel(ip,1)=0.0
          pertvel(ip,2)=0.0
          pertvel(ip,3)=0.0
         end do
c$omp end parallel do
       endif

CTiming
       T2 = SECNDS(0.0)

       if(T2.ge.T1) then
       Tpertvel = T2-T1
       else
       Tpertvel=(24*3600-T1)+T2
       endif
CTiming
       T1 = SECNDS(0.0)
	if (NOTURB) goto 802
c
c   Calculate cross product: u x w in physical space
c   ------------------------------------------------
c   x-cpt  =  tvr*szr - twr*syr 
c   y-cpt  =  twr*sxr - tur*szr 
c   z-cpt  =  tur*syr - tvr*sxr 
c   then overwrite, saving these as sxr, syr, szr
c
c$omp parallel do private(j,i,zu,zv,zw,zwx,zwy,zwz)
      do  k=1,n3
        do  j=1,n2
          do  i=1,n1
            zu=tur(i,j,k)
            zv=tvr(i,j,k)
            zw=twr(i,j,k)
            zwx=sxr(i,j,k)
            zwy=syr(i,j,k)
            zwz=szr(i,j,k)
            sxr(i,j,k)=zv*zwz-zw*zwy
            syr(i,j,k)=zw*zwx-zu*zwz
            szr(i,j,k)=zu*zwy-zv*zwx
          end do
        end do
      end do
c$omp end parallel do
c
      call fft3DRC(sxr,plan1DF,plan2DRC,nnodes)
      call fft3DRC(syr,plan1DF,plan2DRC,nnodes)
      call fft3DRC(szr,plan1DF,plan2DRC,nnodes)
c
c   filter the nonlinear term
c   -----------------------------
c
c$omp parallel do private(j,i)
      do  k=1,n3
        do  j=1,n2
          do  i=1,n1pp
            sxr(i,j,k)=sxr(i,j,k)*filter(i,j,k)
            syr(i,j,k)=syr(i,j,k)*filter(i,j,k)
            szr(i,j,k)=szr(i,j,k)*filter(i,j,k)
          end do
        end do
      end do
c$omp end parallel do
c
c     compute the energy transfer spectrum
c
      if ( mod(istep,nlong).eq.0 .and. istep.gt.0) then
        do k=1,n3
           do j=1,n2
              temp=fsq3(k)+fsq2(j)
              do i=1,n1pp
                 itemp=int(sqrt(temp+fsq1(i))+0.5)
                 if ( itemp .gt. 0 ) then
                 if ( i .le. 2 ) then
                 xmult=1.0
                 else
                 xmult=2.
                 endif
                 usqr=xmult*(cur(i,j,k)*sxr(i,j,k)+cvr(i,j,k)*syr(i,j,k)
     &                      +cwr(i,j,k)*szr(i,j,k) )
                 phi5(itemp)=phi5(itemp)+usqr
                 endif
              enddo
           enddo
        enddo
      itrunc=int(fkstar)
      do i=1,itrunc
      write(13,102) i,phi4(i),phi5(i)
      enddo
      endif
c
c
c   More optional data for the pressure.
c   -----------------------------------
c
      if(pflow)then
      if ( mod(istep,nshort).eq.0 ) then
c
c   Compute now q=(u*u)/2 from tur,tvr,twr still saved
c   in physical space, and save.
c
      do  k=1,n3
        do  j=1,n2
          do  i=1,n1
            tur(i,j,k)=0.5*(tur(i,j,k)**2 + tvr(i,j,k)**2 +
     &                      twr(i,j,k)**2)
          end do
        end do
      end do
c
c
c   Compute the modified pressure P from  div ( s ) in spectral form.
c
      do  k=1,n3
        do  j=1,n2
          tempk=fsq3(k)+fsq2(j)
          do  i=1,n1pp
            tvr(i,j,k)=(fac1(i)*sxr(i,j,k)+fac2(j)*syr(i,j,k)+fac3(k)*
     &            szr(i,j,k))/(tempk+fsq1(i)+smallx)
          end do
        end do
      end do
c
c   Now transform to physical space and combine with q to get the
c   pressure field in physical space.Save P, and p in twr.
c
      call fft3DCR(tvr,plan1DB,plan2DCR,nnodes)
      do  k=1,n3
        do  j=1,n2
          do  i=1,n1
            twr(i,j,k)=tvr(i,j,k)-tur(i,j,k)
          end do
        end do
      end do
c
      icpt=1
      call value(tur,pssp,icpt, lhnode,bg)
      jcpt=2
      call value(tvr,pssp,jcpt, lhnode,bg)
      kcpt=3
      call value(twr,pssp,kcpt, lhnode,bg)
c
      end if
      end if
CTiming
802	continue

       T2 = SECNDS(0.0)
       if(T2.ge.T1) then
       Tflow2  = T2-T1
       else
       Tflow2 =(24*3600-T1)+T2
       endif
c
CTiming

       T1 = SECNDS(0.0)
c
c   Short particle statistics
c
      if ( mod(istep,nshort).eq.0)
     1   call ptdat(numb,npset,ttime,yp,vp,up,pertvel,pssp,vortp,yp0)

      if ( mod(istep,nshort).eq.0) then
c$omp parallel do private(j,i)
        do k=1, 128
         do j=1, 128
          do i=1,128
          con1(i,j,k)=0
          con2(i,j,k)=0
          enddo
         enddo
        enddo
c$omp end parallel do

       do i=0, 30
        pc1(i)=0.0
        pc2(i)=0.0
        pcr(i)=0.0
       enddo



        Ng=128*128*128

        average_con=float(numb)/float(ng)
        common_factor=exp(-average_con)

        do Ic=0, 30
          if(Ic.eq.0) then
            Factorial_Ic=1.0
          elseif(Ic.eq.1) then
            Factorial_Ic=1.0
          else
            Factorial_Ic=Factorial_Ic*float(Ic)
          endif
          pcr(Ic)=average_con**(Ic)*common_factor/Factorial_Ic
        enddo

        do  ip=1,npart

          ii=int(yp(ip,1)/(pi2/128))+1
          if(ii.gt.128) ii=ii-128 
          if(ii.le.0) ii=ii+128
           
          jj=int(yp(ip,2)/(pi2/128))+1
          if(jj.gt.128) jj=jj-128
          if(jj.le.0) jj=jj+128
        

          kk=int(yp(ip,3)/(pi2/128))+1
          if(kk.gt.128) kk=kk-128
          if(kk.le.0) kk=kk+128

          if(ip.le.numb) con1(ii,jj,kk)=con1(ii,jj,kk)+1
          if(ip.gt.numb) con2(ii,jj,kk)=con2(ii,jj,kk)+1
          if(con1(ii,jj,kk).gt.30)  con1(ii,jj,kk)=30
          if(con2(ii,jj,kk).gt.30)  con2(ii,jj,kk)=30
!          print*, "ip=",ip,"ii=",ii,"jj=",jj,"kk=",kk
        end do

        do k=1,128
        do j=1,128
        do i=1,128          
           pc1(con1(i,j,k))= pc1(con1(i,j,k))+1.
           pc2(con2(i,j,k))= pc2(con2(i,j,k))+1.
        enddo
        enddo
        enddo

        do i=0, 30
          pc1(i)=pc1(i)/float(Ng)
          pc2(i)=pc2(i)/float(Ng)
        enddo

        dc1=0.0
        dc2=0.0

        do i=0,30
          dc1=dc1+(pc1(i)-pcr(i))**2
          dc2=dc2+(pc2(i)-pcr(i))**2
        enddo

        write(107, 10001) ttime,dc1, pc1(0),pc1(1),pc1(2),pc1(3)
     1        ,pc1(4),pc1(5),pc1(6),pc1(7),pc1(8),pc1(9),pc1(10)
        write(108, 10001) ttime,dc2, pc2(0),pc2(1),pc2(2),pc2(3)
     1        ,pc2(4),pc2(5),pc2(6),pc2(7),pc2(8),pc2(9),pc2(10)
10001   format(1x,15(1x,f9.6))
       endif
        
c Advance the particle velocity
c   ---------------------------
c    Equations of motion set up here are
c
c        dV/dt= (  u(Y,t) )/tau +  ( W - V )/tau
c
c        dY/dt= V(t)
c
c    The velocity is time-stepped first, the first group of terms is
c   advanced with an AB4 scheme while an AM4 scheme is used for the
c   second group. As W is fixed its contribution is simpler.
c    The position is advanced with an AM4 scheme afterwards.
c    The AB4 scheme has the format for  dz/dt=g  of:
c      z(n+1) = z(n) + (dt/24) * ( 55*g(n) - 59*g(n-1) + 37*g(n-2)
c                                  - 9*g(n-3) )
c    The AM4 scheme has the format 
c      z(n+1) = z(n) + (dt/24) * ( 9*g(n+1) + 19*g(n) - 5*g(n-1)
c                                    + g(n-2) )
c
      do  ic=1,3
c$omp parallel do
        do ip=1,npart
          fvp(ip,ic,0)=vp(ip,ic)
          fdp(ip,ic,0)=up(ip,ic)+pertvel(ip,ic)
        end do
c$omp end parallel do

	if(OVERL)goto 808
       if((istep.gt.0).and.(newcoltot.gt.0)) then
           do ipos=1,ncollpart
           ip=removpart(ipos)
           fvp(ip,ic,1)=fvp(ip,ic,0)
           fvp(ip,ic,2)=fvp(ip,ic,0)
           fdp(ip,ic,1)=fdp(ip,ic,0)
           fdp(ip,ic,2)=fdp(ip,ic,0)
           fdp(ip,ic,3)=fdp(ip,ic,0)
           enddo
       endif
808  	continue

        if ( newpart ) then
c$omp parallel do
          do  ip=1,npart
            fvp(ip,ic,1)=fvp(ip,ic,0)
            fvp(ip,ic,2)=fvp(ip,ic,0)
            fdp(ip,ic,1)=fdp(ip,ic,0)
            fdp(ip,ic,2)=fdp(ip,ic,0)
            fdp(ip,ic,3)=fdp(ip,ic,0)
          end do
c$omp end parallel do
        end if
      end do
c
c   Set newpart=.false. now as all initialisation is done.
c
      newpart=.false.
c
C
c - write particle motion for later detailed analysis
c
c    Advance the particle velocity
c
c$omp parallel do private(dtau,tempv,tempd,ic,Vrel,Rep,fSN,unitv)
        do  ip=1,npart
        Vrel=sqrt((vp(ip,1)-up(ip,1)-pertvel(ip,1))**2+
     :            (vp(ip,2)-up(ip,2)-pertvel(ip,2))**2+
     :            (vp(ip,3)-up(ip,3)-pertvel(ip,3))**2)

        fSN=1.0

        if(ndrag.eq.1) then
        if(ip.le.numb) Rep=Vrel*2.*rad1/rnu
        if(ip.gt.numb) Rep=Vrel*2.*rad2/rnu
        fSN=1.+0.15*Rep**0.687
        endif

      do  ic=1,3
        unitv=0.0
        if ( ic.eq.1 ) unitv=1.0
          dtau=dt/(24.0*tau(ip))*fSN
          tempv= 19.0*fvp(ip,ic,0) - 5.0*fvp(ip,ic,1)
     &           + fvp(ip,ic,2)
          tempd= 55.0*fdp(ip,ic,0) - 59.0*fdp(ip,ic,1)
     &               +37.0*fdp(ip,ic,2) - 9.0*fdp(ip,ic,3)
!          vp(ip,ic)=( vp(ip,ic) + dt*unitv*wp(ip)/tau(ip)*fSN
          vp(ip,ic)=( vp(ip,ic) + dt*unitv*wp(ip)/tau(ip)
     :   - dtau*tempv + dtau*tempd)/( 1.0 + 9.0*dtau )
        end do
      end do
c$omp end parallel do
c
c   Advance the particle position
c   ---------------------------
c
      dtc=dt/24.0
c
      do  ic=1,3
c$omp parallel do
        do  ip=1,npart
        dmove(ip,ic)=dtc*( 9.0*vp(ip,ic) +19.0*fvp(ip,ic,0)
     &               -5.0*fvp(ip,ic,1) + fvp(ip,ic,2))
          yp(ip,ic)=yp(ip,ic)+dmove(ip,ic)
        end do
c$omp end parallel do
      end do
c
c   Ensure that particles remain in box, use periodicity
c
c$omp parallel do
      do  ip=1,npart
        if (  yp(ip,1).ge.hx ) yp(ip,1)=yp(ip,1)-hx
        if (  yp(ip,1).lt.0.0) yp(ip,1)=yp(ip,1)+hx
        if (  yp(ip,2).ge.hy ) yp(ip,2)=yp(ip,2)-hy
        if (  yp(ip,2).lt.0.0) yp(ip,2)=yp(ip,2)+hy
        if (  yp(ip,3).ge.hz ) yp(ip,3)=yp(ip,3)-hz
        if (  yp(ip,3).lt.0.0) yp(ip,3)=yp(ip,3)+hz
!        print*,ip, yp(ip,1),yp(ip,2),yp(ip,3)
      end do
c$omp end parallel do
c
         
c-- save particle in a slice dx for visualization

      if(mod(istep,itvort).eq.0) then
        do  ip=1,npart
            if(yp(ip,3).ge.3.1415926-0.5*dx*0.25.and.
     &         yp(ip,3).le.3.1415926+0.5*dx*0.25) then     
              write(icdump,707) yp(ip,1),yp(ip,2),0.0,
     &                 vp(ip,1),vp(ip,2),0.0,float(ip)
            endif
        enddo
         close(icdump)
         icdump=icdump+1
      endif

C --    computing particle mean square velocity
      if(mod(istep,nshort).eq.0) then

      vpmean1x=0.0
      vpmean1y=0.0
      vpmean1z=0.0
      vpmean2x=0.0
      vpmean2y=0.0
      vpmean2z=0.0
      vpvar1x=0.0
      vpvar1y=0.0
      vpvar1z=0.0
      vpvar2x=0.0
      vpvar2y=0.0
      vpvar2z=0.0
c$omp parallel do reduction(+:vpvar1x,vpvar1y,vpvar1z)
c$omp+ reduction(+:vpvar2x,vpvar2y,vpvar2z)
c$omp+ reduction(+:vpmean1x,vpmean1y,vpmean1z)
c$omp+ reduction(+:vpmean2x,vpmean2y,vpmean2z)
c$omp+ private(ip2)
      do ip=1,numb
      vpmean1x=vpmean1x+vp(ip,1)
      vpmean1y=vpmean1y+vp(ip,2)
      vpmean1z=vpmean1z+vp(ip,3)
      vpvar1x=vpvar1x+vp(ip,1)**2
      vpvar1y=vpvar1y+vp(ip,2)**2
      vpvar1z=vpvar1z+vp(ip,3)**2
      if(npset.gt.1) then
        ip2=numb+ip      
        vpmean2x=vpmean2x+vp(ip2,1)
        vpmean2y=vpmean2y+vp(ip2,2)
        vpmean2z=vpmean2z+vp(ip2,3)
        vpvar2x=vpvar2x+vp(ip2,1)**2
        vpvar2y=vpvar2y+vp(ip2,2)**2
        vpvar2z=vpvar2z+vp(ip2,3)**2
      endif
      enddo
c$omp end parallel do
      vpmean1x=vpmean1x/float(numb)
      vpmean1y=vpmean1y/float(numb)
      vpmean1z=vpmean1z/float(numb)
      vpmean2x=vpmean2x/float(numb)
      vpmean2y=vpmean2y/float(numb)
      vpmean2z=vpmean2z/float(numb)

c	write(106,*) 'vpvar1x,vpvar2x= ',vpvar1x,vpvar2x
c	write(106,*) 'vpvar1x/numb,vpvar2x/numb= ',
c     1 vpvar1x/float(numb),vpvar2x/float(numb)
c	write(106,*) 'vpmean1x,vpmean2x= ',vpmean1x,vpmean2x
c        write(106,*) 'vpmean1x**2,vpmean2x**2= ',vpmean1x**2,
c     1  vpmean2x**2

      vpvar1x=vpvar1x/float(numb)
      vpvar1y=vpvar1y/float(numb)
      vpvar1z=vpvar1z/float(numb)
      vpvar2x=vpvar2x/float(numb)
      vpvar2y=vpvar2y/float(numb)
      vpvar2z=vpvar2z/float(numb)

c	write(106,*) 'vpvar1x,vpvar2x= ',vpvar1x,vpvar2x

      endif

        if(wp(1).eq.wp(npart)) then
         w0diff=wp(1)
        else
         w0diff=abs(vpmean1x-vpmean2x)
         if( (w0diff.gt.wp(1)).or.(w0diff.gt.wp(npart)) ) then
         w0diff=wp(1)
         if(w0diff.gt.wp(npart))w0diff=wp(npart)
         endif
        endif

c
CTiming
       T2 = SECNDS(0.0)
       if(T2.ge.T1) then
       Tpart3 = T2-T1
       else
       Tpart3 =(24*3600-T1)+T2
       endif
CTiming
       T1 = SECNDS(0.0)
c
c   Collision Detection
c   ********************
c
c   Indentifying the particle location relative to
c   the collision detection grid
c
      icc1e=0
      icc2e=0
      icc3e=0
      newcol11=0
      newcol12=0
      newcol22=0
      newcoltot=0
      ncollpart=0
       nR11e=0
       nR12e=0
       nR22e=0

        do i=1,ncd
        do j=1,ncd
        do k=1,ncd
        head(i,j,k)=0
        enddo
        enddo
        enddo

      do ip=1,npart
      ix=1+int(ypp(ip,1)/wcd)
      iy=1+int(ypp(ip,2)/wcd)
      iz=1+int(ypp(ip,3)/wcd)
!      print*,"ip=",ip,"ix=",ix,"iy=",iy,"iz=",iz
c
      if(ix.gt.ncd ) then
      write(103,*)ix,iy,iz,ip
      write(103,*)ttime,ypp(ip,1),ypp(ip,2),
     1   ypp(ip,3),wcd,hx,hy,hz,'COLLDET'
      ix = ncd
      endif

      if(iy.gt.ncd)then
      write(103,*)ix,iy,iz,ip
      write(103,*)ttime,ypp(ip,1),ypp(ip,2),
     1   ypp(ip,3),wcd,hx,hy,hz,'COLLDET'
      iy = ncd
      endif

      if(iz.gt.ncd)then
      write(103,*)ix,iy,iz,ip
      write(103,*)ttime,ypp(ip,1),ypp(ip,2),
     1   ypp(ip,3),wcd,hx,hy,hz,'COLLDET'
      iz = ncd
      endif
c
        list(ip)=head(ix,iy,iz)
        head(ix,iy,iz)=ip
      end do

c$omp parallel do
      do ip=1,npart
      dmoved(ip)=sqrt(dmove(ip,1)**2+dmove(ip,2)**2
     &     +dmove(ip,3)**2)
      end do
c$omp end parallel do
CC
c$omp parallel private(i,j,k,my_thread)
        my_thread = omp_get_thread_num()
        do i=0,179
        nir11n(i,my_thread)=0
        nir12n(i,my_thread)=0
        nir22n(i,my_thread)=0

        wrbin11n(i,my_thread)=0.0
        wrbin12n(i,my_thread)=0.0
        wrbin22n(i,my_thread)=0.0
        enddo
c$omp end parallel
cc     
       sumabswr11=0.0
       sumabswr12=0.0
       sumabswr22=0.0

       print*, "before the hard part!! not parallel now!!!!" 
!c$omp parallel do private(i,j,k,icheck)
!c$omp+ reduction(+:icc1e,icc2e,icc3e)
!c$omp+ reduction(+:newcol11,newcol12,newcol22)
!c$omp+ reduction(+:nR11e,nR12e,nR22e)
!c$omp+ private(dnij1WRGR,xa1,xa2,xa3,xb1,xb2,xb3)
!c$omp+ private(va1,va2,va3,vb1,vb2,vb3)
!c$omp+ private(xr,yr,zr,Vr1,Vr2,Wr)
!c$omp+ private(rcolla,rcollb)
!c$omp+ private(ix,iy,iz,ypa1,ypa2,ypa3)
!c$omp+ private(ix2,sx,iy2,sy,iz2,sz)
!c$omp+ private(ip2,yppb1,yppb2,yppb3,ypb1,ypb2,ypb3)
!c$omp+ private(dnij0,dnij1,tmove,x01,x02,x03,x11,x12,x13)
!c$omp+ private(v01,v02,v03,v11,v12,v13,bx01,bx02,bx03)
!c$omp+ private(bx11,bx12,bx13,bv01,bv02,bv03,bv11,bv12,bv13)
!c$omp+ private(tfac,tfac1,tfac2,tfac3,xa,ya,za,xb,yb,zb,dnij)
!c$omp+ private(rcoll,ncollip1)
!c$omp+ private(ir,my_thread)
        DO IP1=1, NPART

	ncollip1=0

         ix=1+int(ypp(ip1,1)/wcd)
         iy=1+int(ypp(ip1,2)/wcd)
         iz=1+int(ypp(ip1,3)/wcd)
c
      if(ix.gt.ncd ) then
      write(103,*)ix,iy,iz,ip1
      write(103,*)ttime,ypp(ip1,1),ypp(ip1,2),
     1   ypp(ip1,3),wcd,hx,hy,hz
      ix = ncd
      endif

      if(iy.gt.ncd)then
      write(103,*)ix,iy,iz,ip1
      write(103,*)ttime,ypp(ip1,1),ypp(ip1,2),
     1   ypp(ip1,3),wcd,hx,hy,hz
      iy = ncd
      endif

      if(iz.gt.ncd)then
      write(103,*)ix,iy,iz,ip1
      write(103,*)ttime,ypp(ip1,1),ypp(ip1,2),
     1   ypp(ip1,3),wcd,hx,hy,hz
      iz = ncd
      endif
c
       ypa1=ypp(ip1,1)+dmove(ip1,1)
       ypa2=ypp(ip1,2)+dmove(ip1,2)
       ypa3=ypp(ip1,3)+dmove(ip1,3)

        xa1=ypp(ip1,1)
        xa2=ypp(ip1,2)
        xa3=ypp(ip1,3)
       va1=fvp(ip1,1,0)
       va2=fvp(ip1,2,0)
       va3=fvp(ip1,3,0)

        do i=-1,1
        do j=-1,1
        do k=-1,1
        ix2=ix+i
        sx=0.0
        if(ix2.lt.1) then
        ix2=ncd
        sx=-1.
        endif
        if(ix2.gt.ncd)then
        ix2=1
        sx=1.
        endif
c
        iy2=iy+j
        sy=0.0
        if(iy2.lt.1) then
        iy2=ncd
        sy=-1.
        endif
        if(iy2.gt.ncd)then
        iy2=1
        sy=1.
        endif
c
        iz2=iz+k
        sz=0.0
        if(iz2.lt.1) then
        iz2=ncd
        sz=-1.
        endif
        if(iz2.gt.ncd)then
        iz2=1
        sz=1.
        endif
        ip2=head(ix2,iy2,iz2)
c
c	ip2 has to be always greater than ip1 to avoid double counting
c
        if(ip2.eq.0)goto 200
	if(ip2.le.ip1)goto 91

c       if(ip1.eq.5720) then
c        ijole=ip2
c528     write(525,*) ttime,istep,ijole,ix,iy,iz,ip2,ix2,iy2,iz2
c        ijole=list(ijole)
c        if(ijole.eq.0) goto 199
c        goto 528
c       endif

c
199   	if(ip1.le.numb) then
		if(ip2.le.numb) then
		  rcoll=rcoll11
                  rcolla=rcoll11a
                  rcollb=rcoll11b
		else
		  rcoll=rcoll12
                  rcolla=rcoll12a
                  rcollb=rcoll12b
		endif
	endif
	if(ip1.gt.numb) then
                if(ip2.gt.numb) then
                  rcoll=rcoll22
                  rcolla=rcoll22a
                  rcollb=rcoll22b
                else
                  rcoll=rcoll12
                  rcolla=rcoll12a
                  rcollb=rcoll12b
		endif
        endif

       xb1=ypp(ip2,1)+sx*hx
       xb2=ypp(ip2,2)+sy*hy
       xb3=ypp(ip2,3)+sz*hz
       vb1=fvp(ip2,1,0)
       vb2=fvp(ip2,2,0)
       vb3=fvp(ip2,3,0)

       dnij1WRGR=sqrt((xa1-xb1)**2+(xa2-xb2)**2+(xa3-xb3)**2)

c-add bin statistics
C      if(dnij1WRGR.lt.(10.*rcoll))then
       if(dnij1WRGR.ge.rcoll .and. dnij1WRGR.le.(10.*rcoll))then
C       if(dnij1WRGR.ge.(10.*rcoll) .and. dnij1WRGR.le.(19.*rcoll))then
        my_thread = omp_get_thread_num()
c modified on Oct. 9, 2003
        xr=(xa1-xb1)/dnij1WRGR
        yr=(xa2-xb2)/dnij1WRGR
        zr=(xa3-xb3)/dnij1WRGR
        Vr1=va1*xr+va2*yr+va3*zr
        Vr2=vb1*xr+vb2*yr+vb3*zr
        Wr=Vr1-Vr2

        ir=INT( (dnij1WRGR-rcoll) /(0.05*rcoll) )
c        ir=INT( (dnij1WRGR-10.0*rcoll) /(0.05*rcoll) )

        if(ir.gt.179)ir=179
c        if(ir.lt.0) ir=0
        if(ir.ge.0) then
c
          if(ip1.le.numb) then
             if(ip2.le.numb) then
                nir11n(ir,my_thread)=nir11n(ir,my_thread)+1
                wrbin11n(ir,my_thread)=wrbin11n(ir,my_thread)+abs(Wr)
             else
CC I made changes here and below !!!
               nir12n(ir,my_thread)=nir12n(ir,my_thread)+1
               wrbin12n(ir,my_thread)=wrbin12n(ir,my_thread)+abs(Wr)
             endif
          endif

          if(ip1.gt.numb) then
             if(ip2.gt.numb) then
                nir22n(ir,my_thread)=nir22n(ir,my_thread)+1
                wrbin22n(ir,my_thread)=wrbin22n(ir,my_thread)+abs(Wr)
             else
                nir12n(ir,my_thread)=nir12n(ir,my_thread)+1
                wrbin12n(ir,my_thread)=wrbin12n(ir,my_thread)+abs(Wr)
             endif
          endif
        endif
        endif
        
c       if(dnij1WRGR.ge.rcoll .and. dnij1WRGR.le.rcollb )then
c       write(110,*) ttime, 'ip1= ',ip1,' ip2= ',ip2 
c       endif

       if(dnij1WRGR.ge.rcolla .and. dnij1WRGR.le.rcollb )then
c modified on Oct. 22, 1997
        xr=(xa1-xb1)/dnij1WRGR
        yr=(xa2-xb2)/dnij1WRGR
        zr=(xa3-xb3)/dnij1WRGR
        Vr1=va1*xr+va2*yr+va3*zr
        Vr2=vb1*xr+vb2*yr+vb3*zr
        Wr=Vr1-Vr2
c
        if(ip1.le.numb) then
                if(ip2.le.numb) then
                  nR11e=nR11e+1
                  sumabswr11=sumabswr11+abs(Wr)
                  write(81,191)ttime,va1,va2,va3,vb1,vb2,vb3
                  write(71,191)ttime,xr,yr,zr,Wr
                else
                  nR12e=nR12e+1
                  sumabswr12=sumabswr12+abs(Wr)
                  write(82,191)ttime,va1,va2,va3,vb1,vb2,vb3
                  write(72,191)ttime,xr,yr,zr,Wr
                endif
        endif
        if(ip1.gt.numb) then
                if(ip2.gt.numb) then
                  nR22e=nR22e+1
                  sumabswr22=sumabswr22+abs(Wr)
                  write(83,191)ttime,va1,va2,va3,vb1,vb2,vb3
                  write(73,191)ttime,xr,yr,zr,Wr
                else
                  nR12e=nR12e+1
                  sumabswr12=sumabswr12+abs(Wr)
                  write(82,191)ttime,va1,va2,va3,vb1,vb2,vb3
                  write(72,191)ttime,xr,yr,zr,Wr
                endif
        endif

191   format(f10.6,6(1x,f11.5))
       endif

 
       yppb1=ypp(ip2,1)+sx*hx
       yppb2=ypp(ip2,2)+sy*hy
       yppb3=ypp(ip2,3)+sz*hz
c
       ypb1=yppb1+dmove(ip2,1)
       ypb2=yppb2+dmove(ip2,2)
       ypb3=yppb3+dmove(ip2,3)
c---   coarse check
       dnij0=sqrt((ypp(ip1,1)-yppb1)**2
     1 +(ypp(ip1,2)-yppb2)**2+(ypp(ip1,3)-yppb3)**2)
       tmove=rcoll+dmoved(ip1)+dmoved(ip2)
       if(dnij0.gt.tmove)goto 91
       dnij1=sqrt((ypa1-ypb1)**2+(ypa2-ypb2)**2+(ypa3-ypb3)**2)

c        np1=5720 
c        np2=25473 
c       if((ip1.eq.np1).and.(ip2.eq.np2)) then
c      sep=sqrt((yp(np1,1)-yp(np2,1))**2+(yp(np1,2)-yp(np2,2))**2+
c     1         (yp(np1,3)-yp(np2,3))**2)
c      sep=sep/rcoll
c      write(425,*) ttime,yp(np1,1),yp(np1,2),yp(np1,3),yp(np2,1),
c     1                   yp(np2,2),yp(np2,3),istep,sep
c      write(425,*) ttime,ypp(np1,1),ypp(np1,2),ypp(np1,3),ypp(np2,1),
c     1                   ypp(np2,2),ypp(np2,3)
c      write(425,*) yppb1,yppb2,yppb3,'      ',ypb1,ypb2,ypb3
c      write(425,*) ypa1,ypa2,ypa3
c      write(425,*) tmove,'     ',dnij0,dnij1,rcoll
c      write(425,*)
c       endif

C
C Searching for type III collisions
       if(dnij0.le.rcoll )then
       if(dnij1.le.rcoll )then
       x01=ypp(ip1,1)
       x02=ypp(ip1,2)
       x03=ypp(ip1,3)
       x11=ypa1-x01
       x12=ypa2-x02
       x13=ypa3-x03
       v01=fvp(ip1,1,0)
       v02=fvp(ip1,2,0)
       v03=fvp(ip1,3,0)
       v11=vp(ip1,1)
       v12=vp(ip1,2)
       v13=vp(ip1,3)
       bx01=yppb1
       bx02=yppb2
       bx03=yppb3
       bx11=ypb1-bx01
       bx12=ypb2-bx02
       bx13=ypb3-bx03
       bv01=fvp(ip2,1,0)
       bv02=fvp(ip2,2,0)
       bv03=fvp(ip2,3,0)
       bv11=vp(ip2,1)
       bv12=vp(ip2,2)
       bv13=vp(ip2,3)

        do icheck=1,19
       tfac=0.05*float(icheck)
       tfac1=dt*tfac
       tfac2=tfac**2
       tfac3=tfac2*tfac

       xa=x01+v01*tfac1+(3.*x11-(2.*v01+v11)*dt)*tfac2
     1 +((v11+v01)*dt-2.*x11)*tfac3
       ya=x02+v02*tfac1+(3.*x12-(2.*v02+v12)*dt)*tfac2
     1 +((v12+v02)*dt-2.*x12)*tfac3
       za=x03+v03*tfac1+(3.*x13-(2.*v03+v13)*dt)*tfac2
     1 +((v13+v03)*dt-2.*x13)*tfac3

       xb=bx01+bv01*tfac1+(3.*bx11-(2.*bv01+bv11)*dt)*tfac2
     1 +((bv11+bv01)*dt-2.*bx11)*tfac3
       yb=bx02+bv02*tfac1+(3.*bx12-(2.*bv02+bv12)*dt)*tfac2
     1 +((bv12+bv02)*dt-2.*bx12)*tfac3
       zb=bx03+bv03*tfac1+(3.*bx13-(2.*bv03+bv13)*dt)*tfac2
     1 +((bv13+bv03)*dt-2.*bx13)*tfac3

       dnij=sqrt( (xa-xb)**2+(ya-yb)**2+(za-zb)**2 )

       if(dnij.gt.rcoll )then
       icc3e=icc3e+1
       ncollip1=ncollip1+1
c
       goto 92
       endif

        enddo

       endif

        if(OVERL)then
        continue
        else
        write(*,*) 'WARNING --> TYPE III COLLISION ',ip1,ip2,istep
c      write(325,*) ttime,yp(ip1,1),yp(ip1,2),yp(ip1,3),yp(ip2,1),
c     1                   yp(ip2,2),yp(ip2,3),istep
c      write(325,*) dnij0,dnij1,dnij,rcoll

        goto 92
        endif

       goto 91

       endif
C
C searching for type I collisions
c       if(dnij0.gt.rcoll) then
       if(dnij1.le.rcoll) then
	 icc1e=icc1e+1
	 ncollip1=ncollip1+1
         goto 92
       endif
c       endif

C searching for type  II collisions

       x01=ypp(ip1,1)
       x02=ypp(ip1,2)
       x03=ypp(ip1,3)
       x11=ypa1-x01
       x12=ypa2-x02
       x13=ypa3-x03
       v01=fvp(ip1,1,0)
       v02=fvp(ip1,2,0)
       v03=fvp(ip1,3,0)
       v11=vp(ip1,1)
       v12=vp(ip1,2)
       v13=vp(ip1,3)
       bx01=yppb1
       bx02=yppb2
       bx03=yppb3
       bx11=ypb1-bx01
       bx12=ypb2-bx02
       bx13=ypb3-bx03
       bv01=fvp(ip2,1,0)
       bv02=fvp(ip2,2,0)
       bv03=fvp(ip2,3,0)
       bv11=vp(ip2,1)
       bv12=vp(ip2,2)
       bv13=vp(ip2,3)

        do icheck=1,20
       tfac=0.05*float(icheck)
       tfac1=dt*tfac
       tfac2=tfac**2
       tfac3=tfac2*tfac

       xa=x01+v01*tfac1+(3.*x11-(2.*v01+v11)*dt)*tfac2
     1 +((v11+v01)*dt-2.*x11)*tfac3
       ya=x02+v02*tfac1+(3.*x12-(2.*v02+v12)*dt)*tfac2
     1 +((v12+v02)*dt-2.*x12)*tfac3
       za=x03+v03*tfac1+(3.*x13-(2.*v03+v13)*dt)*tfac2
     1 +((v13+v03)*dt-2.*x13)*tfac3

       xb=bx01+bv01*tfac1+(3.*bx11-(2.*bv01+bv11)*dt)*tfac2
     1 +((bv11+bv01)*dt-2.*bx11)*tfac3
       yb=bx02+bv02*tfac1+(3.*bx12-(2.*bv02+bv12)*dt)*tfac2
     1 +((bv12+bv02)*dt-2.*bx12)*tfac3
       zb=bx03+bv03*tfac1+(3.*bx13-(2.*bv03+bv13)*dt)*tfac2
     1 +((bv13+bv03)*dt-2.*bx13)*tfac3

       dnij=sqrt( (xa-xb)**2+(ya-yb)**2+(za-zb)**2 )
       if(dnij.lt.rcoll )then
       icc2e=icc2e+1
       ncollip1=ncollip1+1
c
       goto 92
       endif
        enddo

       goto 91
92     continue
        if(ip1.le.numb) then
                if(ip2.le.numb) then
                  newcol11=newcol11+1
                else
                  newcol12=newcol12+1
		endif
        endif
        if(ip1.gt.numb) then
                if(ip2.gt.numb) then
                  newcol22=newcol22+1
                else
                  newcol12=newcol12+1
		endif
        endif

c       write(109,*) ttime, 'ip1= ',ip1,' ip2= ',ip2

c This section is NEW. We want to get the approaching angle when
c particles are colliding

        v01=fvp(ip1,1,0)
        v02=fvp(ip1,2,0)
        v03=fvp(ip1,3,0)
        bv01=fvp(ip2,1,0)
        bv02=fvp(ip2,2,0)
        bv03=fvp(ip2,3,0)
        xr=(ypp(ip1,1)-yppb1)/dnij0
        yr=(ypp(ip1,2)-yppb2)/dnij0
        zr=(ypp(ip1,3)-yppb3)/dnij0
        Vr1=v01*xr+v02*yr+v03*zr
        Vr2=bv01*xr+bv02*yr+bv03*zr
        Wr=Vr1-Vr2

        if(ip1.le.numb) then
                if(ip2.le.numb) then
                  write(84,191)ttime,v01,v02,v03,bv01,bv02,bv03
                  write(74,191)ttime,xr,yr,zr,Wr
                else
                  write(85,191)ttime,v01,v02,v03,bv01,bv02,bv03
                  write(75,191)ttime,xr,yr,zr,Wr
                endif
        endif
        if(ip1.gt.numb) then
                if(ip2.gt.numb) then
                  write(86,191)ttime,v01,v02,v03,bv01,bv02,bv03
                  write(76,191)ttime,xr,yr,zr,Wr
                else
                  write(85,191)ttime,v01,v02,v03,bv01,bv02,bv03
                  write(75,191)ttime,xr,yr,zr,Wr
                endif
        endif


!c$omp critical
        ncollpart=ncollpart+1
        removpart(ncollpart)=ip2
!c$omp end critical

c  We assume that only one collision is possible for a given particle
C     goto 99   WE WANT TO CHECK EVEN IF MORE THAN ONE COLLISION 
c
91    ip2=list(ip2)
        if(ip2.eq.0)goto 200
	if(ip2.le.ip1)goto 91
        goto 199
200     continue
      enddo
      enddo
      enddo
99    continue

       if(ncollip1.ne.0) then
!c$omp critical
        ncollpart=ncollpart+1
        removpart(ncollpart)=ip1
!c$omp end critical
       endif

      ENDDO
!c$omp end parallel do
c     
      print*,"after the hard part!"
      do i=0,179
       do j=0,nnodes-1
      nir11(i) = nir11(i) + nir11n(i,j)
      nir12(i) = nir12(i) + nir12n(i,j)
      nir22(i) = nir22(i) + nir22n(i,j)

      wrbin11(i) = wrbin11(i) + wrbin11n(i,j)
      wrbin12(i) = wrbin12(i) + wrbin12n(i,j)
      wrbin22(i) = wrbin22(i) + wrbin22n(i,j)
       enddo
      enddo

      if(mod(istep,nshort).eq.0) then
          id4=mod(istep,100000)/10000
          id3=mod(istep,10000)/1000
          id2=mod(istep,1000)/100
          id1=mod(istep,100)/10
          id0=mod(istep,10)
          open (200,file='./gr/gr'//char(48+id4)//char(48+id3)
     1//char(48+id2)//char(48+id1)//char(48+id0)//'.dat')

        if(npset.eq.1) then

          collrate11=newcol11/(2.*pi)**3/dt
          averagecon11=0.5*npart*(npart-1)/(2*pi)**3
          sumabswr11=sumabswr11/nR11e
          vshell=(4./3.)*pi*rcoll11**3*(1.025**3-0.975**3)
          gr11w=nR11e/vshell/averagecon11
          simucollrate11=2*pi*rcoll11**2*sumabswr11*gr11w*
     1                      (npart/(2*pi)**3)**2/2
          write(109,707)ttime,gr11w,sumabswr11,collrate11,simucollrate11

          do i=0,179
            wrbin11(i) = wrbin11(i)/nir11(i)
            vs11=(4./3.)*pi*rcoll11**3*((1+(i+1)*0.05)**3-(1+i*0.05)**3)
            gr11(i)=nir11(i)/vs11/averagecon11
            write(200,707) (i+0.5)*(0.05)+1, wrbin11(i), gr11(i)
          enddo

       else
          collrate11=newcol11/(2.*pi)**3/dt
          collrate12=newcol12/(2.*pi)**3/dt
          collrate22=newcol22/(2.*pi)**3/dt
          halfpart=float(npart)/2
          averagecon11=0.5*(halfpart)*(halfpart-1)/(2*pi)**3
          averagecon22=0.5*(halfpart)*(halfpart)/(2*pi)**3
          averagecon12=(halfpart)*(halfpart)/(2*pi)**3

          sumabswr11=sumabswr11/nR11e
          vshell=(4./3.)*pi*rcoll11**3*(1.025**3-0.975**3)
          gr11w=nR11e/vshell/averagecon11

          sumabswr22=sumabswr22/nR22e
          vshell=(4./3.)*pi*rcoll22**3*(1.025**3-0.975**3)
          gr22w=nR22e/vshell/averagecon22

          sumabswr12=sumabswr12/nR12e
          vshell=(4./3.)*pi*rcoll12**3*(1.025**3-0.975**3)
          gr12w=nR12e/vshell/averagecon12

          simucollrate11=2*pi*rcoll11**2*sumabswr11*gr11w*
     &                   (npart/2/(2*pi)**3)**2/2
          simucollrate22=2*pi*rcoll22**2*sumabswr22*gr22w*
     &                   (npart/2/(2*pi)**3)**2/2
          simucollrate12=2*pi*rcoll12**2*sumabswr12*gr12w*
     &                   (npart/2/(2*pi)**3)**2

          write(109,707)ttime,gr11w,sumabswr11,collrate11,simucollrate11
          write(110,707)ttime,gr12w,sumabswr12,collrate12,simucollrate12
          write(115,707)ttime,gr22w,sumabswr22,collrate22,simucollrate22


         do i=0,179
            wrbin11(i) = wrbin11(i)/nir11(i)
            wrbin12(i) = wrbin12(i)/nir12(i)
            wrbin22(i) = wrbin22(i)/nir22(i)
            vs11=(4./3.)*pi*rcoll11**3*((1+(i+1)*0.05)**3-(1+i*0.05)**3)
            vs12=(4./3.)*pi*rcoll12**3*((1+(i+1)*0.05)**3-(1+i*0.05)**3)
            vs22=(4./3.)*pi*rcoll22**3*((1+(i+1)*0.05)**3-(1+i*0.05)**3)

            gr11(i)=nir11(i)/vs11/averagecon11
            gr22(i)=nir22(i)/vs22/averagecon22
            gr12(i)=nir12(i)/vs12/averagecon12

            write(200,707) (i+0.5)*(0.05)+1, wrbin11(i), gr11(i)
     1      ,wrbin12(i), gr12(i),wrbin22(i), gr22(i)
        enddo
       endif 
       close(200)  
      endif
      
c
      itocol11=itocol11+newcol11
      itocol12=itocol12+newcol12
      itocol22=itocol22+newcol22
      icc1 = icc1 + icc1e
      icc2 = icc2 + icc2e
      icc3 = icc3 + icc3e
      newcoltot=newcol11+newcol12+newcol22
      nR11 = nR11 + nR11e
      nR12 = nR12 + nR12e
      nR22 = nR22 + nR22e


c Relocate particles that just collided

	if(OVERL) goto 805
        if(ncollpart.gt.0) then
         do ippc=1,ncollpart
          ipID=removpart(ippc)
           ar1=ranpp(iseedp,iyp,ivp)
           ar2=ranpp(iseedp,iyp,ivp)
           ar3=ranpp(iseedp,iyp,ivp)
           yp(ipID,1)= ar1*hx
           yp(ipID,2)= ar2*hy
           yp(ipID,3)= ar3*hz
         enddo
        endif
805	continue

c Check for Overlapping of the newly relocated particles

	if(OVERL) goto 806
       IF(newcoltot.ne.0) THEN
c
c   Overlapping Detection and Correction
c   ***********************************
c
c   Indentifying the particle location relative to
c   the ovarelapping detection grid
c
1982    continue
c
        do i=1,ncd
        do j=1,ncd
        do k=1,ncd
        head(i,j,k)=0
        enddo
        enddo
        enddo

      do ip=1,npart
      ix=1+int(yp(ip,1)/wcd)
      iy=1+int(yp(ip,2)/wcd)
      iz=1+int(yp(ip,3)/wcd)
c
      if(ix.gt.ncd ) then
      ix = ncd
      endif

      if(iy.gt.ncd)then
      iy = ncd
      endif

      if(iz.gt.ncd)then
      iz = ncd
      endif
c
        list(ip)=head(ix,iy,iz)
        head(ix,iy,iz)=ip
      end do

        noverl=0
c
c$omp parallel do private(i,j,k)
c$omp+ reduction(+:noverl)
c$omp+ private(ip1,ix,iy,iz)
c$omp+ private(ix2,sx,iy2,sy,iz2,sz)
c$omp+ private(ip2,ypb1,ypb2,ypb3)
c$omp+ private(dnij,ar1,ar2,ar3)
c$omp+ private(rcoll)
c$omp+ shared(iseedp,iyp,ivp)
        DO IPOS=1, NCOLLPART

          ip1=removpart(ipos)

         ix=1+int(yp(ip1,1)/wcd)
         iy=1+int(yp(ip1,2)/wcd)
         iz=1+int(yp(ip1,3)/wcd)
c
      if(ix.gt.ncd ) then
      write(103,*)ix,iy,iz,ip1
      write(103,*)ttime,yp(ip1,1),yp(ip1,2),
     1   yp(ip1,3),wcd,hx,hy,hz
      ix = ncd
      endif

      if(iy.gt.ncd)then
      write(103,*)ix,iy,iz,ip1
      write(103,*)ttime,yp(ip1,1),yp(ip1,2),
     1   yp(ip1,3),wcd,hx,hy,hz
      iy = ncd
      endif

      if(iz.gt.ncd)then
      write(103,*)ix,iy,iz,ip1
      write(103,*)ttime,yp(ip1,1),yp(ip1,2),
     1   yp(ip1,3),wcd,hx,hy,hz
      iz = ncd
      endif
c
        do i=-1,1
        do j=-1,1
        do k=-1,1
        ix2=ix+i
        sx=0.0
        if(ix2.lt.1) then
        ix2=ncd
        sx=-1.
        endif
        if(ix2.gt.ncd)then
        ix2=1
        sx=1.
        endif
c
        iy2=iy+j
        sy=0.0
        if(iy2.lt.1) then
        iy2=ncd
        sy=-1.
        endif
        if(iy2.gt.ncd)then
        iy2=1
        sy=1.
        endif
c
        iz2=iz+k
        sz=0.0
        if(iz2.lt.1) then
        iz2=ncd
        sz=-1.
        endif
        if(iz2.gt.ncd)then
        iz2=1
        sz=1.
        endif
        ip2=head(ix2,iy2,iz2)
c
c
        if(ip2.eq.0)goto 2002
        if(ip2.eq.ip1)goto 912
c
1992     if(ip1.le.numb) then
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

       ypb1=yp(ip2,1)+sx*hx
       ypb2=yp(ip2,2)+sy*hy
       ypb3=yp(ip2,3)+sz*hz

       dnij=sqrt((yp(ip1,1)-ypb1)**2
     1 +(yp(ip1,2)-ypb2)**2+(yp(ip1,3)-ypb3)**2)

        if(dnij.le.rcoll) then
           noverl=noverl+1
           ar1=ranpp(iseedp,iyp,ivp)
           ar2=ranpp(iseedp,iyp,ivp)
           ar3=ranpp(iseedp,iyp,ivp)
           yp(ip2,1)= ar1*hx
           yp(ip2,2)= ar2*hy
           yp(ip2,3)= ar3*hz
        endif

912    ip2=list(ip2)
        if(ip2.eq.0)goto 2002
        if(ip2.eq.ip1)goto 912
        goto 1992
2002     continue
      enddo
      enddo
      enddo
      ENDDO

c        write(*,*) 'noverlAFTERCOLL= ',noverl

        if(noverl.ne.0) goto 1982
c
c  END OF DETECTION OF OVERLAPPING
c

       ENDIF
806	continue 
c
CTiming
       T2 = SECNDS(0.0)
       if(T2.ge.T1) then
       Tcollision = T2-T1
       else
       Tcollision=(24*3600-T1)+T2
       endif

CTiming
       T1 = SECNDS(0.0 )

	if (NOTURB) goto 803
C-Change 6/31/95
c
c   ADVANCE THE FLOW FIELD
c
c   If this is a new run replace AB2 by an Euler step
c   -------------------------------------------------
c
      if ( newflo ) then
c
c$omp parallel do private(j,i)
      do  k=1,n3
        do j=1,n2
          do i=1,n1pp
            sxp(i,j,k)=sxr(i,j,k)
            syp(i,j,k)=syr(i,j,k)
            szp(i,j,k)=szr(i,j,k)
          end do
        end do
      end do
c$omp end parallel do
c
      newflo=.false.
c
      endif
c
c   Advance the flow field - 1st part non-linear terms advanced
c   -----------------------------------------------------------
c                  2
c              nu k  dt  ^      3 dt     ^           dt      ^
c    bu = (1 - --------) u    + ---- [   s   ]   - ---- [    s   ]
c                 2       (n)     2           1(n)    2           1(n-1)
c
c                  2
c              nu k  dt  ^      3 dt     ^           dt      ^
c    bv = (1 - --------) v    + ---- [   s   ]   - ---- [    s   ]
c                 2       (n)     2           2(n)    2           2(n-1)
c
c
c                  2
c              nu k  dt  ^      3 dt     ^           dt      ^
c    bw = (1 - --------) w    + ---- [   s   ]   - ---- [   s   ]
c                 2       (n)     2           3(n)    2           3(n-1)
c

       if(les) then
c$omp parallel do private(j,i,facsq,temp,itemp,sgs_visc1)
        do k=1,n3
           do j=1,n2
              facsq=fsq3(k)+fsq2(j)
              do i=1,n1pp
                 itemp=int(sqrt(facsq+fsq1(i))+0.5)
                 if(itemp.le.itrunc) then
                   sgs_visc1=(rnuf/rnu)*(rnu+sgs_visck(itemp))
                 else
                   sgs_visc1=rnuf
                 endif
            temp=1.0-sgs_visc1*(fsq1(i)+facsq)
            cur(i,j,k)=temp*cur(i,j,k)+dtr*sxr(i,j,k)-dtp*sxp(i,j,k)
            cvr(i,j,k)=temp*cvr(i,j,k)+dtr*syr(i,j,k)-dtp*syp(i,j,k)
            cwr(i,j,k)=temp*cwr(i,j,k)+dtr*szr(i,j,k)-dtp*szp(i,j,k)
         end do
        end do
       end do
c$omp end parallel do
        else
c$omp parallel do private(j,i,facsq,temp)
      do  k=1,n3
        do  j=1,n2
          facsq=fsq3(k)+fsq2(j)
          do  i=1,n1pp
            temp=1.0-rnuf*(fsq1(i)+facsq)
            cur(i,j,k)=temp*cur(i,j,k)+dtr*sxr(i,j,k)-dtp*sxp(i,j,k)
            cvr(i,j,k)=temp*cvr(i,j,k)+dtr*syr(i,j,k)-dtp*syp(i,j,k)
            cwr(i,j,k)=temp*cwr(i,j,k)+dtr*szr(i,j,k)-dtp*szp(i,j,k)
          end do
        end do
      end do
c$omp end parallel do
      endif
c
c
c   Add the random forcing
c   ----------------------
c   Note that we could have combined the forcing a1r,a2,a3r with
c   the terms sxr,syr,szr at each step ( and so sxp,etc. ), but
c   as this is random forcing this Euler step is as effective.
c
      if ( forced ) then
CRANDOM
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
      end if
c             
c   Advance the flow field   -  linear terms
c   ----------------------------------------
c   The presure term is eliminated by continuity eqn. Note
c   the pressure here is a pseudo-pressure not the real one.
c   To get the pressure explicitly one must solve a Poisson eq.
c   as done in the earlier data step.
c                                                               
c
c   Finally the velocity field in spectral form at (n+1) step is
c                                        2
c       ^              bu - k1 ( k.bu )/k
c       u       = ----------------------
c        (n+1)                  2
c                  (1 + 0.5 nu k  dt)
c
c   etc.
c
c

       if(les) then
c$omp parallel do private(j,i,temp,itemp,temp1,temp2,pres,sgs_visc1)
        do k=1,n3
           do j=1,n2
              temp=fsq3(k)+fsq2(j)
              do i=1,n1pp
                 itemp=int(sqrt(temp+fsq1(i))+0.5)
                 if(itemp.le.itrunc) then
                   sgs_visc1=(rnuf/rnu)*(rnu+sgs_visck(itemp)) 
                 else
                   sgs_visc1=rnuf
                 endif
                   temp1=1.0/(1.0+sgs_visc1*(fsq1(i)+temp))
                   temp2=1.0/(fsq1(i)+temp+smallx)
            pres=temp2*(fac1(i)*cur(i,j,k)+fac2(j)*cvr(i,j,k)
     &           +fac3(k)*cwr(i,j,k))
            cur(i,j,k)=temp1*(cur(i,j,k)-fac1(i)*pres)
            cvr(i,j,k)=temp1*(cvr(i,j,k)-fac2(j)*pres)
            cwr(i,j,k)=temp1*(cwr(i,j,k)-fac3(k)*pres)
          end do
        end do
       end do
Cc$omp end parallel do
        else
c$omp parallel do private(j,i,temp,temp1,temp2,pres)
       do  k=1,n3
        do  j=1,n2
          temp=fsq3(k)+fsq2(j)+smallx
          do  i=1,n1pp
            temp1=1.0/(1.0+rnuf*(fsq1(i)+temp))
            temp2=1.0/(fsq1(i)+temp)
            pres=temp2*(fac1(i)*cur(i,j,k)+fac2(j)*cvr(i,j,k)
     &           +fac3(k)*cwr(i,j,k))
            cur(i,j,k)=temp1*(cur(i,j,k)-fac1(i)*pres)
            cvr(i,j,k)=temp1*(cvr(i,j,k)-fac2(j)*pres)
            cwr(i,j,k)=temp1*(cwr(i,j,k)-fac3(k)*pres)
          end do
        end do
      end do
c$omp end parallel do
      endif 
c
c   Clean and filter data
c   ---------------------
c   Enforce symmetry of fourier coefficients in k_1=0 plane
c   and truncate the n/2 wavenos.
c
      call clean(cur)
      call clean(cvr)
      call clean(cwr)
c
c
c   Force a zero mean flow
c   ----------------------
c
      cur(1,1,1)=0.0
      cur(2,1,1)=0.0
      cvr(1,1,1)=0.0
      cvr(2,1,1)=0.0
      cwr(1,1,1)=0.0
      cwr(2,1,1)=0.0
c
c   Use the filter here if needed
c   -----------------------------
c
c$omp parallel do private(j,i)
      do  k=1,n3
        do  j=1,n2
          do  i=1,n1pp
            cur(i,j,k)=cur(i,j,k)*filter(i,j,k)
            cvr(i,j,k)=cvr(i,j,k)*filter(i,j,k)
            cwr(i,j,k)=cwr(i,j,k)*filter(i,j,k)
          end do
        end do
      end do
c$omp end parallel do
c
c   Update (u x w + f )
c   -------------------
c
c$omp parallel do private(j,i)
      do  k=1,n3
        do   j=1,n2
          do  i =1,n1pp
            sxp(i,j,k)=sxr(i,j,k)
            syp(i,j,k)=syr(i,j,k)
            szp(i,j,k)=szr(i,j,k)
          end do
        end do
      end do
c$omp end parallel do
CTiming
803	continue
       T2 = SECNDS(0.0 )

       if(T2.ge.T1) then
       Tflow3 = T2-T1
       else
       Tflow3 =(24*3600-T1)+T2
       endif
CTiming
       T1 = SECNDS(0.0 )

c
c   Update particle history
c   ---------------------
c SECNDS(0.0 )
      do  ic=1,3
c$omp parallel do
        do  ip=1,npart
          fvp(ip,ic,2)=fvp(ip,ic,1)
          fvp(ip,ic,1)=fvp(ip,ic,0)
          fdp(ip,ic,3)=fdp(ip,ic,2)
          fdp(ip,ic,2)=fdp(ip,ic,1)
          fdp(ip,ic,1)=fdp(ip,ic,0)
          ypp(ip,ic)=yp(ip,ic)
        end do
c$omp end parallel do
      end do
c
CTiming
       T2 = SECNDS(0.0 )
       if(T2.ge.T1) then
       Tpart4 = T2-T1
       else
       Tpart4=(24*3600-T1)+T2
       endif

       if(T2.ge.TTT1) then
       TSTEP = T2-TTT1
       else
       TSTEP=(24*3600-T1)+T2
       endif

       TCHECK = Tflow1 + Tflow2 + Tflow3
     1   +Tpart1 + Tpart2 + Tpart3 + Tpart4+Tpertvel+Tcollision

       write(*,*)'Tflow1,Tflow2,Tflow3,Tpart1,Tpart2,Tpart3,Tpart4,
     1   Tcollision,Tpertvel,TSTEP,TCHECK'

       write(*,501)Tflow1,Tflow2,Tflow3,Tpart1,Tpart2,
     1  Tpart3,Tpart4,Tcollision,Tpertvel,TSTEP,TCHECK

       write(*,501)Tflow1/tstep,Tflow2/tstep,Tflow3/tstep,
     1  Tpart1/tstep,Tpart2/tstep,Tpart3/tstep,
     1  Tpart4/tstep,Tcollision/tstep,Tpertvel/tstep,TCHECK/tstep
501    format(11f10.6)
c
c
c   Advance time
c   ------------
c
      istep=istep+1
      ttime=ttime+dt
c
      if(mod(istep,nshort).eq.0) then
       write(80,303)ttime,nR11,nR12,nR22,vpmean1x,vpmean1y,
     1            vpmean1z,vpmean2x,vpmean2y,vpmean2z

c       write(8,301)ttime,numb,itocol11,itocol12,itocol22,
c     1    vpmean1x,vpmean1y,vpmean1z,vpmean2x,vpmean2y,vpmean2z,
c     1    icc1,icc2,icc3
c       write(*,301)ttime,numb,itocol11,itocol12,itocol22,
c     1    vpmean1x,vpmean1y,vpmean1z,vpmean2x,vpmean2y,vpmean2z,
c     1    icc1,icc2,icc3
	write(106,302) ttime,vpmean1x,vpmean1y,vpmean1z,vpmean2x,vpmean2y,
     1   vpmean2z,vpvar1x,vpvar1y,vpvar1z,vpvar2x,vpvar2y,vpvar2z
301    format(2x,f9.6,4i10,6f10.4,3i5)
302    format(2x,7f10.5,6f14.6)
303    format(2x,f10.6,3i10,6f10.4)
      icc1=0
      icc2=0
      icc3=0

       write(111,777)ttime,(nir11(ir),ir=0,59)
       write(112,777)ttime,(nir12(ir),ir=0,59)
       write(113,777)ttime,(nir22(ir),ir=0,59)

       write(211,777)ttime,(nir11(ir),ir=60,119)
       write(212,777)ttime,(nir12(ir),ir=60,119)
       write(213,777)ttime,(nir22(ir),ir=60,119)

       write(311,777)ttime,(nir11(ir),ir=120,179)
       write(312,777)ttime,(nir12(ir),ir=120,179)
       write(313,777)ttime,(nir22(ir),ir=120,179)

777    format(2x,f10.6,60i7)
       write(121,778)ttime,(wrbin11(ir),ir=0,59)
       write(122,778)ttime,(wrbin12(ir),ir=0,59)
       write(123,778)ttime,(wrbin22(ir),ir=0,59)
       write(221,778)ttime,(wrbin11(ir),ir=60,119)
       write(222,778)ttime,(wrbin12(ir),ir=60,119)
       write(223,778)ttime,(wrbin22(ir),ir=60,119)
       write(321,778)ttime,(wrbin11(ir),ir=120,179)
       write(322,778)ttime,(wrbin12(ir),ir=120,179)
       write(323,778)ttime,(wrbin22(ir),ir=120,179)
778    format(2x,f10.6,60(1pe12.3))

        do ir=0,179
        nir11(ir)=0
        nir12(ir)=0
        nir22(ir)=0
        wrbin11(ir)=0.0
        wrbin12(ir)=0.0
        wrbin22(ir)=0.0
        enddo

      endif
c
c   Save a time frame with new data if appropriate
c   ----------------------------------------------
c  Need to set this up as needed for job.
c
c      if ( mod ( istep, nlong ).eq. 0 ) then
c
c        open ( unit=51,file='frame.51', form='unformatted',
c    &       status='unknown' )
c
c        write (51) cur,cvr,cwr
c        write (51) a1r,a2r,a3r
c        write (51) yp,vp
c      end if

c
c   End of time-step
c   ----------------
c
       write(*,*)'ending istep=',istep !, fdate()
c      call flush

      if (istep.lt.nhalt) go to 100
c
c   Format statements for unit=16
c   -----------------------------
c
 101  format( 2x, i6, 9(1pe12.4) )
 102  format(2x,i6,7(2x,1pe12.4) )
 103  format( 2(2x,1pe12.4))
 105  format(2x,15(1pe12.4) )
 1001 format(//' Start of new run with particle code ')
 1010 format(/'  n1,n2,n3 =',3(3x,i5))
 1015 format(/'  parameters rnu,dt,var,tf=',/4(3x,1pe12.4))
 1020 format(/'  time at start of run =',3x,e17.10)
 1025 format(/'  logical values of newflo,newpart,forced are',
     & /3(3x,l7))
 1030 format(/'  nhalt,nlong,nshort=',3(3x,i5))
 1035 format(/'  iseedp,iseedf =',2(3x,I20))
c
 1080 format(///' npart,nset, numb =',3(3x,i8))
 1100 format(/'   error in setting npset')
 1110 format(//'  particle parameters used are tau,wp:')
707   format(2x,8f15.6)
901   format(2x,i6,6f8.3)
178   format(2x,'particle leaving box = ',I7,f12.5)
271    format(2x,2f8.4)
c
c
c   End of run; dump data for a restart
c   -----------------------------------
c
      open ( unit=12, file='./HI2/endrun.flo', 
     & form='unformatted',status='unknown' )
      write (12) cur,cvr,cwr
      write (12) sxp,syp,szp
CRANDOM
      write (12)iseedf,iyf,ivf
      write (12) a1r,a2r,a3r
      write (12) b1r,b2r,b3r
      close(12)

        if(les) then
          open(unit=11, file='./HI2/les_sgs_visc_old.dat')
          do i=1,itrunc
          write(11,*) sgs_visck_old(i)
          enddo
          close(11)
         endif

!      do k=1,128
!        do j=1,128
!          do i=1, 130
!             cur128(i,j,k)=cur(i,j,k)
!             cvr128(i,j,k)=cvr(i,j,k)
!             cwr128(i,j,k)=cwr(i,j,k)
!             sxp128(i,j,k)=sxp(i,j,k)
!             syp128(i,j,k)=syp(i,j,k)
!             szp128(i,j,k)=szp(i,j,k)
!             if(i.le.66.and.j.le.64.and.k.le.64) then
!               cur64(i,j,k)=cur(i,j,k)
!               cvr64(i,j,k)=cvr(i,j,k)
!               cwr64(i,j,k)=cwr(i,j,k)
!               sxp64(i,j,k)=sxp(i,j,k)
!               syp64(i,j,k)=syp(i,j,k)
!               szp64(i,j,k)=szp(i,j,k)
!             endif
!          enddo
!        enddo
!      enddo



!      open(unit=12,file='./HI2/endrunles64.flo',
!     & form='unformatted',status='unknown')

!      write (12) cur64,cvr64,cwr64
!      write (12) sxp64,syp64,szp64
!CRANDOM
!      write (12)iseedf,iyf,ivf
!      write (12) a1r,a2r,a3r
!      write (12) b1r,b2r,b3r      
!      close(12)

!      open(unit=12,file='./HI2/endrunles128.flo',
!     & form='unformatted',status='unknown')

!      write (12) cur128,cvr128,cwr128
!      write (12) sxp128,syp128,szp128
!CRANDOM
!      write (12)iseedf,iyf,ivf
!      write (12) a1r,a2r,a3r
!      write (12) b1r,b2r,b3r
!      close(12)


      open ( unit=13,file='./HI2/endrun.part', 
     & form='unformatted', status='unknown' )
      write (13) yp,vp,fvp,fdp
CRANDOM
      write (13)iseedp,iyp,ivp
      close(13)
c
        do i=0,(nnodes-1)
       call destroy_plan1DF(plan1DF(i))
       call destroy_plan1DB(plan1DB(i))
       call destroy_plan2DRC(plan2DRC(i))
       call destroy_plan2DCR(plan2DCR(i))
        enddo
      Tend=SECNDS(0.0)
      Tduration=Tend-Tbegin
      write(16,*) "program begins at:",Tbegin
      write(16,*) "grogram ends at:",Tend
      write(16,*) "Time cost", Tduration

      stop
      end 

