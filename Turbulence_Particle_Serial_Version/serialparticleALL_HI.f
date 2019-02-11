!            Implement fftw 
!
!    Periodic box code to study the geometric collision rates between two groups
!    of particles in an isotropic turbulence or other periodic flow. 
!    The turbulence is (locally) homogeneous and isotropic with forcing at
!    at low wavenumbers to sustain the flow. An AB2/CN2 scheme is used to advance the flow.
!    Calculations are in Fourier space and pseudospectral.
!
!       The particles move in this version as 
!          dV/dt = ( u(Y,t) - V(t) + W )/ tau .
!    The parameters tau and W characterise the motion. Values of
!    u  are found by Lagrange interpolation. Axes are fixed
!    so that  W = ( wp , 0, 0 ) and wp >0 is the still fluid settling
!    velocity.
!
!    Collision detections are used to count individual collision events.
!
!     VERSION March, 2003
!
!	INCLUDE HYDRODYNAMIC INTERACTIONS  (09/16/03)
!
      implicit real(a-h,o-z)
!
      include 'params.inc'
!
      parameter(n1p=n1+1)
      parameter(n1pp=n1+2,n2d=n2*2,n3d=n3*2)
      parameter(nwds=n1pp*n2*n3,npts=n1*n2*n3)
      parameter(n1h=n1/2,n2h=n2/2,n3h=n3/2)
      parameter(n1hp=n1h+1,n2hp=n2h+1,n3hp=n3h+1)
      parameter(inc2x=n1pp,inc3x=n1pp*n2)
      parameter(inc2y=n1h+1,inc3y=(n1h+1)*n2)
      parameter(naux=65000)
      parameter(scaleCR=1.0,scaleRC=1./float(n1*n2*n3))
      parameter(nparthalf=npart/2)
!RANDOM
      parameter (NTAB=32)
!
      real*4, dimension(:,:,:), allocatable :: cur,cvr,cwr
      real*4, dimension(:,:,:), allocatable :: tur,tvr,twr
      real*4, dimension(:,:,:), allocatable :: sxp,syp,szp
      real*4, dimension(:,:,:), allocatable :: sxr,syr,szr
      real*4, dimension(:,:,:), allocatable :: sqx,sqy,sqz
      real*4, dimension(:,:,:), allocatable :: filter,tempa
      real*8, dimension(:), allocatable :: aux
      integer, dimension(:,:,:,:), allocatable :: headd
      integer, dimension(:,:), allocatable :: nir11n,nir12n
      integer, dimension(:,:), allocatable :: nir22n
      real*4, dimension(:,:), allocatable :: wrbin11n,wrbin12n
      real*4, dimension(:,:), allocatable :: wrbin22n

      real*8, dimension(:,:,:), allocatable :: curt
      real*8, dimension(:,:), allocatable :: location
      real*8, dimension(:,:,:), allocatable :: b1rt,b2rt,b3rt
!
      dimension fac1(n1pp), fac2(n2), fac3(n3)
      dimension fsq1(n1pp), fsq2(n2), fsq3(n3)
!
      dimension probe(8),efx(6),phi1(n1pp),phi2(n1pp),	&
        phi3(n1pp),phi4(n1pp),phi5(n1pp),vdist(200)
!
!   Particle arrays.
!   ----------------
      real*4, dimension(:), allocatable ::tau,wp,dmoved
      real*4, dimension(:,:), allocatable :: yp,ypp,vp,up,dmove
      real*4, dimension(:,:), allocatable :: pertvel
      real*4, dimension(:,:), allocatable :: vortp,pssp,yp0
      real*4, dimension(:,:,:), allocatable :: fvp,fdp,bg
      integer, dimension(:,:), allocatable :: lhnode
      integer, dimension(:), allocatable :: removpart
!      integer, dimension(:,:), allocatable :: ipa,ipb
      integer, dimension(:), allocatable :: list
      integer, dimension(:), allocatable :: slab

      dimension dxyz(3)

!
      integer alpha,forcing
      integer head(ncd,ncd,ncd)
      integer nir11(0:179),nir12(0:179),nir22(0:179),pdf(0:360)
      real wrbin11(0:179),wrbin12(0:179),wrbin22(0:179)
      
!
!   Flow forcing arrays.
!   --------------------
      dimension a1r(6,5,5), a2r(6,5,5), a3r(6,5,5),		&
                b1r(6,5,5), b2r(6,5,5), b3r(6,5,5)
!
!
!     double precision rand
!     external function rand,srand
!RANDOM
       integer iseedp,iyp,ivp(NTAB)
       integer iseedf,iyf,ivf(NTAB)

!     FORCING VARIABLES
      real*4          :: force(2),k2,k2_e
      integer, dimension(:), allocatable :: co1, co2
      real*4, dimension(:), allocatable :: tmp1, enum


!       character*28 fconc,fvort
       character*50 fvort
       character*29 direct
       character*34 directvort
       character*25 directinflow
       character*35 directinpart
       character*34 directvel
       character*36 directpart
       character*33 directpart2
 
!
      logical newflo, newpart, forced, pflow, HYDINT, NOTURB,	&	
              OVERL
!
      common /waveno/fac1,fac2,fac3,fsq1,fsq2,fsq3
      common /noise/a1r,a2r,a3r,b1r,b2r,b3r
      common /parms/xkf,var,tf,dt
!     common /interp/lhnode,bg
!     common /part1/yp,vp,up,pssp,vortp
!     common /part2/tau,wp
!     common /rannum/iseedf
      common /geomrad/rad1,rad2,rcoll11,rcoll12,rcoll22,numb,wpvcd
      common /particle/w0diff
      common /boxsize/hx,hy,hz
      common /tiempo/ttime,istep,nhalt
      common /machine/nproc

!============================================================================

! path to output data
      direct='/ptmp/rosa/evolving_OMP/stat/'

! path to vorticity field
      directvort='/ptmp/rosa/evolving_OMP/vort/vort.'

! path to start velocity flow field
!      directinflow='/endrun.flo'

! path to file with initial particles location
       directinpart='/ptmp/rosa/evolving_OMP/endrun.part'

! path to file with output velocity field
      directvel ='/ptmp/rosa/evolving_OMP/endrun.flo'

! path to file with output particle location
      directpart='/ptmp/rosa/evolving_OMP/endrun.part2'

      directpart2='/ptmp/rosa/evolving_OMP/part/part'

!============================================================================

      pi=4.*atan(1.)
      pi2=2.*pi
      my_thread=1
      nnodes=1
!
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
!      allocate ( ipa(npart,3) )
!      allocate ( ipb(npart,3) )
!
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
      allocate ( aux(naux) )
      allocate ( slab(8) )

      allocate ( sqx(n1pp,n2,n3) )
      allocate ( sqy(n1pp,n2,n3) )
      allocate ( sqz(n1pp,n2,n3) )

      allocate ( curt(n1pp,n2,n3) )
      allocate ( location(3,npart) )
      allocate ( b1rt(6,5,5) )
      allocate ( b2rt(6,5,5) )
      allocate ( b3rt(6,5,5) )
!
      allocate ( headd(ncd,ncd,ncd,1) )
      allocate ( nir11n(0:179,1) )
      allocate ( nir12n(0:179,1) )
      allocate ( nir22n(0:179,1) )
      allocate ( wrbin11n(0:179,1) )
      allocate ( wrbin12n(0:179,1) )
      allocate ( wrbin22n(0:179,1) )

      open(12,file='/ptmp/rosa/evolving_OMP/part.dat',&
                 status='unknown',form='unformatted')

      read(12)location
      close(12)

      Tint=0.0
      Tflow=0.0
      Tadv=0.0
      Tpbc=0.0
      Tcd=0.0
      Tod=0.0

!   Set up Fourier grid
!   -------------------
!   hx, hy, hz = box size
!   fx, fy, fz = fundamental wave numbers
!   If n1 > n2 this corresponds to an elongated box in 
!   the vertical.
      hx = pi2
      hy = pi2
      hz = pi2
      DXYZ(1)=HX/real(N1)
      DXYZ(2)=HY/real(N2)
      DXYZ(3)=HZ/real(N3)
      fx = pi2/hx
      fy = pi2/hy
      fz = pi2/hz
!
!   fac1, fac2, fac3 = arrays of wave numbers
!   -----------------------------------------
! -  x : ( 0 --> n1h ) x fx in real/imag pairs
! -  y : ( 0 --> n2h, -(n2h-1) --> -1 ) x fy
! -  z : ( 0 --> n3h, -(n3h-1) --> -1 ) x fz
!
!   fsq1, fsq2, fsq3 = arrays of squares of wave numbers
!
      do  i=1,n1pp,2
        fac1(i)=float((i-1)/2)*fx
        fsq1(i)=fac1(i)**2
        fac1(i+1)=fac1(i)
        fsq1(i+1)=fsq1(i)
      end do

      do  j=1,n2
        jj=j-1
        if ( j .gt. n2hp ) jj=jj-n2
        fac2(j)=float(jj)*fy
        fsq2(j)=fac2(j)**2
      end do
!
      do  k=1,n3
        kk=k-1
        if ( k .gt. n3hp ) kk=kk-n3
        fac3(k)=float(kk)*fz
        fsq3(k)=fac3(k)**2
      end do
!
!
!   Set up spectral filter
!   ----------------------
!   specify filter(i,j,k) here
!
      fkstar=fac1(n1pp)-1.5
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
!   Set up forcing parameters
!   -------------------------
!   This scheme assumes a cubic box!
!
!   xkf= forcing radius
!   var= variance of white noise
!   tf = time constant
!   xnf= no. of modes forced
      xkf=sqrt(8.0)
      var=447.31
      tf=0.038
      xnf=0.0
!
!   Calculate no. of modes forced
!
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

!   Open data files for regular input/output
!   ----------------------------------------
!  input.dat contains basic information to control the job,
!  output.dat is a monitoring log file for the job.
!  probe.dat contains the short set of flow statistics.
      open ( unit=103, file= direct//'index-check.dat')
      open ( unit=106, file= direct//'vpmean-var.dat')
      open ( unit=15, file='partin.dat')
      open ( unit=16, file=direct//'partout.dat')

      open ( unit=8, file=direct//'collision.dat')
      open ( unit=80,file=direct//'detecpairs.dat',status='unknown')
      open ( unit=71,file=direct//'Wr11.dat', status='unknown')
      open ( unit=81,file=direct//'Vp11.dat', status='unknown')
      open ( unit=72,file=direct//'Wr12.dat', status='unknown')
      open ( unit=82,file=direct//'Vp12.dat', status='unknown')
      open ( unit=73,file=direct//'Wr22.dat', status='unknown')
      open ( unit=83,file=direct//'Vp22.dat', status='unknown')

      open ( unit=74,file=direct//'Wr11c.dat', status='unknown')
      open ( unit=84,file=direct//'Vp11c.dat', status='unknown')
      open ( unit=75,file=direct//'Wr12c.dat', status='unknown')
      open ( unit=85,file=direct//'Vp12c.dat', status='unknown')
      open ( unit=76,file=direct//'Wr22c.dat', status='unknown')
      open ( unit=86,file=direct//'Vp22c.dat', status='unknown')

      open ( unit=111,file=direct//'nir11a.dat',status='unknown')
      open ( unit=112,file=direct//'nir12a.dat',status='unknown')
      open ( unit=113,file=direct//'nir22a.dat',status='unknown')

      open ( unit=211,file=direct//'nir11b.dat',status='unknown')
      open ( unit=212,file=direct//'nir12b.dat',status='unknown')
      open ( unit=213,file=direct//'nir22b.dat',status='unknown')

      open ( unit=311,file=direct//'nir11c.dat',status='unknown')
      open ( unit=312,file=direct//'nir12c.dat',status='unknown')
      open ( unit=313,file=direct//'nir22c.dat',status='unknown')

      open ( unit=121,file=direct//'wrbin11a.dat',status='unknown')
      open ( unit=122,file=direct//'wrbin12a.dat',status='unknown')
      open ( unit=123,file=direct//'wrbin22a.dat',status='unknown')

      open ( unit=221,file=direct//'wrbin11b.dat',status='unknown')
      open ( unit=222,file=direct//'wrbin12b.dat',status='unknown')
      open ( unit=223,file=direct//'wrbin22b.dat',status='unknown')

      open ( unit=321,file=direct//'wrbin11c.dat',status='unknown')
      open ( unit=322,file=direct//'wrbin12c.dat',status='unknown')
      open ( unit=323,file=direct//'wrbin22c.dat',status='unknown')
      open ( unit=324,file=direct//'angle.dat',status='unknown')
!
!   Define physical parameters
!   --------------------------
!
!   ndrag=1 for nonlinear drag force, 0 for Stokes drag
!
!   arad1=particle radius of group 1/dx
!   arad2=particle radius of group 2/dx
!
!   rnu = viscosity
!
!   dt = time step size 
!   dtr, dtp = weighted time steps for adams-bashforth
!   rnuf used for viscosity term used in time-stepping
!
!   newflo=.true. for start of new flow, .false. for restart
!   newpart=.true. for start of new particle system
!   forced=.true. if flow is driven by random forcing
!
!   istep = time step counter
!   nhalt = number of time step to be run
!   nlong= number of steps between long data output
!   nshort= number of steps between short data output
!
      dx=pi2/float(n2)
!
      read(15,*)ndrag
      read(15,*)rnu,dt
!
!
      dtr=1.5*dt
      dtp=0.5*dt
      rnuf=0.5*dt*rnu
      ttime=0.0
      smallx=1.0e-18
!
      read (15,*) newflo
      read (15,*) forced
!
      read (15,*) newpart
      read (15,*) pflow
      read (15,*) HYDINT
      read (15,*) NOTURB
      read (15,*) nfactor
      read (15,*) OVERL


	if(NOTURB) then
	  hx=hx/float(nfactor)
	  hy=hy/float(nfactor)
          hz=hz/float(nfactor)
	endif

 	if(HYDINT) OVERL=.false.

      read (15,*) nhalt,nlong,nshort
      read (15,*) itvort
      read (15,*) ttime

      if(pflow) then
      open ( unit=18, file=direct//'probe.dat')
      open ( unit=19, file=direct//'spectrum.dat')
      open ( unit=17, file=direct//'flowstat.dat')
      open ( unit=13, file=direct//'transfer.dat')
      open ( unit=14, file=direct//'vortpdf.dat')
      endif

!     Set iseeds: iseedp for seeding particles in a new particle
!     system and iseedf for the random forcing.

      read (15,*) iseedp,iseedf

!RANDOM
      if(newpart)then
      iyp = 0
      iseedp = - iseedp
      endif

      if(newflo)then
      iyf = 0
      iseedf = - iseedf
      endif

!
!    Start the output file for job
!    -----------------------------
!
      write (16,1001)
      write (16,1010)  n1,n2,n3
      write (16,1015)  rnu,dt,var,tf
      write (16,1020)  ttime
      write (16,1025)  newflo,newpart,forced
      write (16,1030)  nhalt,nlong,nshort
      write (16,1035)  iseedp,iseedf

!
!    Set up particle parameters
!    ------------------------
!   Particles are grouped into npset equal lots of numb
!   particles. The total should be npart. Each set will have
!   a different tau and wp value.
!   Assume for now that npset is .le. 8 .
!

      read (15,*) etk0,tauk0
        vk0=etk0/tauk0
      write (16,*)'etk0,tauk0,vk0=',etk0,tauk0,vk0

! radius normalized by kolmogorov length scale
      read(15,*)arad1,arad2
      write(16,*)'arad1,arad2= ',arad1,arad2
      rad1 = arad1*etk0
      rad2 = arad2*etk0

      read (15,*) npset
      numb=npart/npset
      write (16,1080) npart,npset,numb
      if ( (numb*npset).ne.npart)  then
        write (16,1100)
        stop
      end if
!
!    Open data files for short particle statistics
!
      if (npset.ge.1) open( unit=21,file=direct//'part1a.dat',status='unknown')
      if (npset.ge.1) open( unit=31,file=direct//'part1b.dat',status='unknown')
      if (npset.ge.2) open( unit=22,file=direct//'part2a.dat',status='unknown')
      if (npset.ge.2) open( unit=32,file=direct//'part2b.dat',status='unknown')
!
      write (16,1110)
!
      do  iset=1,npset
        read (15,*) tauset,wpset
        write (16,*) tauset,wpset

        do  jp=1,numb
          ip=jp+(iset-1)*numb
          tau(ip)=tauset*tauk0
          wp(ip)=wpset*vk0
        end do
      end do

        if(wp(1).eq.wp(npart)) then
         w0diff=wp(1)
        else
	 w0diff=abs(wp(1)-wp(npart))
         if( (w0diff.gt.wp(1)).or.(w0diff.gt.wp(npart)) ) then
         w0diff=wp(1)
         if(w0diff.gt.wp(npart))w0diff=wp(npart)
         endif
        endif

       read(15,*)forcing
       read(15,*)integration_scheme


!    If this is a restart, load data cur,..,szp
!    ------------------------------------------

      if ( newflo ) then
        continue
      else
      open(11,file='/ptmp/rosa/evolving_OMP/velo/vx', &
                       status='unknown',form='unformatted')
        read (11) curt
      close(11)
         cur(:,:,:) = curt(:,:,:)

      open(11,file='/ptmp/rosa/evolving_OMP/velo/vy', &
                       status='unknown',form='unformatted')
        read (11) curt
      close(11)
        cvr(:,:,:) = curt(:,:,:)

      open(11,file='/ptmp/rosa/evolving_OMP/velo/vz', &
                       status='unknown',form='unformatted')
        read (11) curt
      close(11)
         cwr(:,:,:) = curt(:,:,:)

      open(11,file='/ptmp/rosa/evolving_OMP/velo/ox', &
                       status='unknown',form='unformatted')
        read (11) curt
      close(11)
        sxp(:,:,:) = curt(:,:,:)

      open(11,file='/ptmp/rosa/evolving_OMP/velo/oy', &
                       status='unknown',form='unformatted')
        read (11) curt
      close(11)
        syp(:,:,:) = curt(:,:,:)

      open(11,file='/ptmp/rosa/evolving_OMP/velo/oz', &
                       status='unknown',form='unformatted')
        read (11) curt
      close(11)
        szp(:,:,:) = curt(:,:,:)

       deallocate(curt)
       
         write(*,*)'Dane wczytane'

      if (forcing.eq.1)then
        open(11,file='/ptmp/rosa/evolving_OMP/velo/force.1', &
                       status='unknown',form='unformatted')

        read(11)iseedf
        read(11)ivf
        read(11)iyf
        read(11)b1rt,b2rt,b3rt

        close(11)

        do  k=1,5
          do j=1,5
            do i=1,6
              b1r( i,j,k) = b1rt(i,j,k)
              b2r( i,j,k) = b2rt(i,j,k)
              b3r( i,j,k) = b3rt(i,j,k)
            end do
          end do
        end do
      endif

      endif

!--------------------------------------------
! initialization of velocity flow field if deterministic forcing has been chosen.

       if (forced.and.(forcing.eq.2).and.newflo)then

        call gaussian(cur,pi)
        call gaussian(cvr,pi)
        call gaussian(cwr,pi)

        call clean (cur)
        call clean (cvr)
        call clean (cwr)

        allocate ( tmp1(int(fkstar)) )
        allocate ( enum(int(fkstar)) )
        allocate ( co1 (int(fkstar)) )
        allocate ( co2 (int(fkstar)) )

        force(1)=0.555440
        force(2)=0.159843

        do  m=1,int(fkstar)
         tmp1(m)=0.
         co1 (m)=0.
         co2 (m)=0.
        enddo

        do  k=1,n3
         do  i=1,n1pp,2
          do  j=1,n2

           k2 = fac1(i)**2 + fac2(j)**2 + fac3(k)**2
           if(k2.eq.0) k2=1e-5
           ik2 = int( sqrt(k2)+0.5 )

           tmpp = ( fac1(i)*cur(i,j,k) + fac2(j)*cvr(i,j,k) + fac3(k)*cwr(i,j,k) ) / k2
           cur(i,j,k) = cur(i,j,k) - fac1(i)*tmpp
           cvr(i,j,k) = cvr(i,j,k) - fac2(j)*tmpp
           cwr(i,j,k) = cwr(i,j,k) - fac3(k)*tmpp

           tmpp = ( fac1(i)*cur(i+1,j,k) + fac2(j)*cvr(i+1,j,k) + fac3(k)*cwr(i+1,j,k) ) / k2
           cur(i+1,j,k) = cur(i+1,j,k) - fac1(i)*tmpp
           cvr(i+1,j,k) = cvr(i+1,j,k) - fac2(j)*tmpp
           cwr(i+1,j,k) = cwr(i+1,j,k) - fac3(k)*tmpp

           tot = cur(i,j,k)**2 + cur(i+1,j,k)**2 +                &
                 cvr(i,j,k)**2 + cvr(i+1,j,k)**2 +                &
                 cwr(i,j,k)**2 + cwr(i+1,j,k)**2
           if (i.eq.1) tot = tot / 2.


           do is=1,int(fkstar)
            if (ik2.eq.is) tmp1(is) = tmp1(is) + tot
            if((ik2.eq.is).and.(fac1(is).gt.0.1)) co1(is)=co1(is)+1
            if((ik2.eq.is).and.(fac1(is).lt.0.1)) co2(is)=co2(is)+1
           enddo

          enddo
         enddo
        enddo

        do is=1,int(fkstar)
         enum(is) = 2.*co1(is) + co2(is)
        enddo

        do  k=1,n3
         do  i=1,n1pp,2
          do  j=1,n2

           k2=fac1(i)**2+fac2(j)**2+fac3(k)**2
           if(k2.eq.0) k2=1e-5
           ik2 = int( sqrt(k2)+0.5 )

           do is=1, int(fkstar)
            e2=is+0.5
            e3=is-0.5

            vol=4./3.*pi*(e2**3-e3**3)
            cccc=force(1)*4./3.*pi*(1.5**3-0.5**3)/enum(1)
            e1=sqrt(cccc/tmp1(is)*enum(is)/vol)/(real(is))**0.833333

            if (ik2.eq.is) then
             cur(i,j,k) = cur(i,j,k) * e1
             cur(i+1,j,k)  = cur(i+1,j,k) * e1
             cvr(i,j,k) = cvr(i,j,k) * e1
             cvr(i+1,j,k)  = cvr(i+1,j,k) * e1
             cwr(i,j,k) = cwr(i,j,k) * e1
             cwr(i+1,j,k)  = cwr(i+1,j,k) * e1
            endif
            if (ik2.gt.fkstar) then
             cur(i,j,k)  = 0.
             cur(i+1,j,k) = 0.
             cvr(i,j,k)  = 0.
             cvr(i+1,j,k) = 0.
             cwr(i,j,k)  = 0.
             cwr(i+1,j,k) = 0.
            endif
           enddo

          enddo
         enddo
        enddo

        sxp=cur
        syp=cvr
        szp=swr

        deallocate (tmp1)
        deallocate (co1)
        deallocate (co2)
        deallocate (enum)
      endif


!
! parameters for collision detection
!
       wcd=hx/float(ncd)
       wpvcd=hx/float(npvcd)
       rcoll11= rad1 + rad1
       rcoll12= rad1 + rad2
       rcoll22= rad2 + rad2
       write(16,*)'wcd, rcoll11, rcoll12, rcoll22 = '		&
      , wcd, rcoll11, rcoll12, rcoll22

        do alpha=0,360
         pdf(alpha)=0
        enddo


        do ir=0,179
        nir11(ir)=0
        nir12(ir)=0
        nir22(ir)=0
        wrbin11(ir)=0.0
        wrbin12(ir)=0.0
        wrbin22(ir)=0.0
        enddo

       rcoll11a = 0.99*rcoll11
       rcoll11b = 1.01*rcoll11
       rcoll12a = 0.99*rcoll12
       rcoll12b = 1.01*rcoll12
       rcoll22a = 0.99*rcoll22
       rcoll22b = 1.01*rcoll22

      istep=0

      if ( newpart ) then

      do  iset=1,npset
        do  jp=1,numb
        ip = jp + (iset-1)*numb
!RANDOM
          ar1=ranpp(iseedp,iyp,ivp)
          ar2=ranpp(iseedp,iyp,ivp)
          ar3=ranpp(iseedp,iyp,ivp)

          yp(ip,1)= ar1*hx
          yp(ip,2)= ar2*hy
          yp(ip,3)= ar3*hz
        end do
        end do


! If you want use location from MPI

        do ip=1,npart

          yp(ip,1)= pi2*location(1,ip)
          yp(ip,2)= pi2*location(2,ip)
          yp(ip,3)= pi2*location(3,ip)

        enddo

 789  format( 3pe12.4)


         write(*,*)'Before ovelapping dettection'

!
!   Overlapping Detection and Correction
!   ***********************************
!
!   Indentifying the particle location relative to
!   the ovarelapping detection grid
!
1981 	continue
!
!  First ensure that particles remain in box, use periodicity
!
      do  ip=1,npart
        if (  yp(ip,1).ge.hx ) yp(ip,1)=yp(ip,1)-hx
        if (  yp(ip,1).lt.0.0) yp(ip,1)=yp(ip,1)+hx
        if (  yp(ip,2).ge.hy ) yp(ip,2)=yp(ip,2)-hy
        if (  yp(ip,2).lt.0.0) yp(ip,2)=yp(ip,2)+hy
        if (  yp(ip,3).ge.hz ) yp(ip,3)=yp(ip,3)-hz
        if (  yp(ip,3).lt.0.0) yp(ip,3)=yp(ip,3)+hz
      end do
!


	if(OVERL)goto 809
!	write(*,*) 'Before Overl Det'

        do i=1,ncd
        do j=1,ncd
        do k=1,ncd
        headd(i,j,k,my_thread)=0
        enddo
        enddo
        enddo

       do ip=1,npart
       ix=1+int(yp(ip,1)/wcd)
       iy=1+int(yp(ip,2)/wcd)
       iz=1+int(yp(ip,3)/wcd)
!
       if(ix.gt.ncd ) then
       ix = ncd
       endif

       if(iy.gt.ncd)then
       iy = ncd
       endif

       if(iz.gt.ncd)then
       iz = ncd
       endif
!
        list(ip)=headd(ix,iy,iz,my_thread)
        headd(ix,iy,iz,my_thread)=ip
       end do

        do l1=1,ncd
        do l2=1,ncd
        do l3=1,ncd
        head(l1,l2,l3)=headd(l1,l2,l3,1)
        enddo
        enddo
        enddo
!
	noverl=0
!
        DO IP1=1, NPART
         ix=1+int(yp(ip1,1)/wcd)
         iy=1+int(yp(ip1,2)/wcd)
         iz=1+int(yp(ip1,3)/wcd)
!
      if(ix.gt.ncd ) then
      write(103,*)ix,iy,iz,ip1
      write(103,*)ttime,yp(ip1,1),yp(ip1,2),	&
        yp(ip1,3),wcd,hx,hy,hz
      ix = ncd
      endif

      if(iy.gt.ncd)then
      write(103,*)ix,iy,iz,ip1
      write(103,*)ttime,yp(ip1,1),yp(ip1,2),	&
        yp(ip1,3),wcd,hx,hy,hz
      iy = ncd
      endif

      if(iz.gt.ncd)then
      write(103,*)ix,iy,iz,ip1
      write(103,*)ttime,yp(ip1,1),yp(ip1,2),	&
        yp(ip1,3),wcd,hx,hy,hz
      iz = ncd
      endif
!
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
!
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
!
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
!
!       ip2 has to be always greater than ip1 to avoid double counting
!
        if(ip2.eq.0)goto 2001
        if(ip2.le.ip1)goto 911
!
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

       dnij=sqrt((yp(ip1,1)-ypb1)**2		&
      +(yp(ip1,2)-ypb2)**2+(yp(ip1,3)-ypb3)**2)

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

!	write(*,*) 'noverl= ',noverl

	if(noverl.ne.0) goto 1981
!
!  END OF DETECTION OF OVERLAPPING
!

	write(*,*) 'After Overl Det'

809	continue
!
      else
        open ( unit=12, file=directinpart,&
        form='unformatted', status='old')
        read (12) yp,vp,fvp,fdp
!RANDOM
        read(12)iseedp,iyp,ivp
        close(12)

      endif

      do  iset=1,npset
        do  jp=1,numb
        ip = jp + (iset-1)*numb
          yp0(ip,1)=yp(ip,1)
          yp0(ip,2)=yp(ip,2)
          yp0(ip,3)=yp(ip,3)
          ypp(ip,1)=yp(ip,1)
          ypp(ip,2)=yp(ip,2)
          ypp(ip,3)=yp(ip,3)
        end do
        end do

!
!   Start of time step loop
!   -----------------------
!
!     Note at this point the flow data at time level n is
!   available and can be used to advance the particles to time
!   level n+1. We need to interpolate the flow data, so
!   some initial setup for this is done and can be applied
!   as needed below. 
!     Once the particles are moved the flow is advanced.
!
!     Note the particles are reintroduced into the box by
!   periodicity if they move out of it. So we can assume that
!   0<= y1 <= hx, etc.
!
!
      itocol11=0
      itocol12=0
      itocol22=0
      icc1=0
      icc2=0
      icc3=0
      nR11=0
      nR12=0
      nR22=0

100  continue
      write(*,*)'begining istep=',istep


	if(NOTURB) goto 800 

       T1 = TIMEF( )/1000.00

!   Set up interpolation factors bg for the particles and
!   locate grid points. lhnode is an index for the nearest
!   grid point coordinate to the left of the particle. The
!   value of lhnode is between 0 and n-1.
!
      do  ix=1,3
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
          bg(ip,3,ix)=( 12.0 - 4.0*z - 15.0*z2 + 5.0*z3 + 3.0*z4 		&
                      - z5 )/12.0
          bg(ip,4,ix)=( 12.0*z + 8.0*z2 - 7.0*z3 - 2.0*z4 + z5 )/12.0
          bg(ip,5,ix)=( - 6.0*z - z2 + 7.0*z3 + z4 - z5 )/24.0
          bg(ip,6,ix)=( 4.0*z - 5.0*z3 + z5 )/120.0
        end do
      end do

800	continue

       T2 = TIMEF( )/1000.00
       Tint = Tint + (T2-T1)

!Timing

        if(NOTURB) goto 801

       T1 = TIMEF( )/1000.00
!
!   Copy fourier coefficients and obtain velocity
!   in physical space. Store in tur, tvr, twr.
!
            tur = cur
            tvr = cvr
            twr = cwr

      CALL SCRFT3(tur,inc2y,inc3y,tur,                          &
          inc2x,inc3x,n1,n2,n3,-1,scaleCR,aux,naux)
      CALL SCRFT3(tvr,inc2y,inc3y,tvr,                          &
          inc2x,inc3x,n1,n2,n3,-1,scaleCR,aux,naux)
      CALL SCRFT3(twr,inc2y,inc3y,twr,                          &
          inc2x,inc3x,n1,n2,n3,-1,scaleCR,aux,naux)

!
!   Calculate vorticity in fourier space and store in sxr,syr,szr
!   -------------------------------------------------------------
!   wz  = dv/dx - du/dy
!  
!   wy  = du/dz - dw/dx
!  
!   wx  = dw/dy - dv/dz
!  
!
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

            szr(i,j,k)  = -fac1(i)*zvi+fac2(j)*zui
            szr(ip,j,k) = +fac1(i)*zvr-fac2(j)*zur

            syr(i,j,k)  = -fac3(k)*zui+fac1(i)*zwi
            syr(ip,j,k) = +fac3(k)*zur-fac1(i)*zwr

            sxr(i,j,k)  = -fac2(j)*zwi+fac3(k)*zvi
            sxr(ip,j,k) = +fac2(j)*zwr-fac3(k)*zvr
          end do
        end do
      end do

!
!
!   Transform vorticity to physical space
!   -------------------------------------
!
      CALL SCRFT3(sxr,inc2y,inc3y,sxr,                          &
          inc2x,inc3x,n1,n2,n3,-1,scaleCR,aux,naux)
      CALL SCRFT3(syr,inc2y,inc3y,syr,                          &
          inc2x,inc3x,n1,n2,n3,-1,scaleCR,aux,naux)
      CALL SCRFT3(szr,inc2y,inc3y,szr,                          &
          inc2x,inc3x,n1,n2,n3,-1,scaleCR,aux,naux)

       T2 = TIMEF( )/1000.00
       Tflow = Tflow + (T2-T1)

!======================== P O S T P R O C E S S I N G ================================

!-- save vorticity field for visualization
      if(mod(istep,itvort).eq.0) then

!      vrms=0.0

        do k=1,n3
           do j=1,n2
              do i=1,n1
      tempa(i,j,k)=sqrt(sxr(i,j,k)**2+syr(i,j,k)**2+szr(i,j,k)**2)
!      vrms=vrms+tempa(i,j,k)**2
               enddo
            enddo
        enddo
!
      num = istep/itvort+1
!
      if (num.le.10) then
      i1d=num
      fvort = directvort//char(i1d+48)
      endif
!
      if ((num.ge.10).and.(num.le.100)) then
      i1d=num/10
      i2d=num-i1d*10
      fvort = directvort//char(i1d+48)//char(i2d+48)
      endif
!
      if ((num.ge.100).and.(num.le.1000)) then
      i1d=num/100
      i2d=num/10-i1d*10
      i3d=num-100*i1d-10*i2d
      fvort = directvort//char(i1d+48)//char(i2d+48)//char(i3d+48)
      endif
!
      if ((num.ge.1000).and.(num.le.10000)) then
      i1d=num/1000
      i2d=num/100-10*i1d
      i3d=num/10-100*i1d-10*i2d
      i4d=num-1000*i1d-100*i2d-10*i3d
      fvort = directvort//char(i1d+48)//char(i2d+48)    &
      //char(i3d+48)//char(i4d+48)
      endif

      open( unit=60,file=fvort,status='unknown')
      write(60,707)(((tempa(i,j,k),i=1,128,2),j=1,128,2),k=1,128,2)
      close(60)
      endif
!
!   Short flow statistics if needed.
!   -------------------------------
!    Do this in physical space for mean square velocity
!    and vorticity. KE dissipation rate probe(7) is found
!    from the vorticity using the periodic bc's/homogeneity.
!    The rate of working by the random forcing is probe(8).
!

      if(pflow)then
      if ( mod(istep,nshort).eq.0 .and. istep.gt.0) then
!     compute energy and dissipation spectrum
        do i=1,n1pp
           phi1(i)=0.0
           phi2(i)=0.0
           phi3(i)=0.0
           phi4(i)=0.0
           phi5(i)=0.0
        enddo
!
!
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
                 usqr=xmult*(cur(i,j,k)**2+cvr(i,j,k)**2		&
                                            +cwr(i,j,k)**2)
                 phi1(itemp)=phi1(itemp)+usqr
           phi4(itemp)=phi4(itemp)+usqr*(temp+fsq1(i))
!
           phi3(itemp)=phi3(itemp)+xmult*(cur(i,j,k)**2*fsq1(i)		&
                   +cvr(i,j,k)**2*fsq2(j)+cwr(i,j,k)**2*fsq3(k))
           phi2(itemp)=phi2(itemp)+2.*rnu*usqr*(temp+fsq1(i))
                 endif
              enddo
           enddo
        enddo
!
!
      itrunc=int(fkstar)

      do i=1,itrunc
      write(19,102) i,phi1(i),phi2(i)
      enddo
!
        do  ic=1,8
          probe(ic)=0.0
        end do

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
!

        if(forcing.eq.1)then 
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
              efx(i)=a1r(i,j,k)*cur(i,jj,kk)+		&
                    a2r(i,j,k)*cvr(i,jj,kk)+		&
                    a3r(i,j,k)*cwr(i,jj,kk)+efx(i)
            end do
          end do
        end do
        efsum= (efx(1)+efx(2))+2.0*(efx(3)+efx(4)+efx(5)+	&
               efx(6))
        probe(8)=efsum
        endif
!
        write(18,101) istep,ttime,probe

        if ( mod(istep,nlong).eq.0) then
        ediss=probe(7)
        eta=(rnu**3/ediss)**(1./(6.-2.0))
        tk = (rnu/ediss)**(1./(3.-1.0))
        vk = eta/tk
!       write(20,*)'ediss=  ',ediss
!       write(20,*)'eta=    ',eta
!       write(20,*)'vk=     ',vk
!       write(20,*)'tk=     ',tk
        uprime=sqrt( (probe(1)+probe(2)+probe(3))/3.0)
        tmse=sqrt(uprime*uprime/tmse)
        ree=uprime*tmse**(2.-1.0)/rnu
        xl=uprime**3/ediss
        et=xl/uprime
        resol=fkstar*eta

        temp=0.0
        do k=1,n3
           do j=1,n2
              do i=1,n1
               outa=abs(tur(i,j,k))+abs(tvr(i,j,k))	&
                   +abs(twr(i,j,k))
               temp = max(temp,outa)
          enddo
         enddo
        enddo

!----  Computing pdf distribution for vorticity magnitude
        vmax=0.0

        do k=1,n3
           do j=1,n2
              do i=1,n1
           sss=sqrt(sxr(i,j,k)**2+syr(i,j,k)**2+szr(i,j,k)**2)
         vmax = max(vmax,sss)
          enddo
         enddo
        enddo

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
!
        CFL=temp*dt/dx
        cfl2=rnu*dt/dx**2
!       write(20,*)'  tmse=  ',tmse
!       write(20,*)'  ree=  ',ree
!       write(20,*)'  L_f=  ',xl
!       write(20,*)'  T_e=  ',et
!       write(20,*)'uprime=  ',uprime
!       write(20,*)'  resol=  ',resol
!       write(20,*)'  CFL=  ',CFL
!       write(20,*)'  CFL2=  ',cfl2

        write(17,105)ttime,uprime,xl,ediss,et,tmse,		&
            eta,vk,tk,ree,resol,CFL,cfl2
!        write(*,*)istep,resol,CFL
        endif
!
      end if
      end if

!Timing

801	continue


!   Interpolate flow velocity ( and vorticity if needed)
!   ---------------------------------------------------

	if (NOTURB) then
        do  ip=1,npart
          up(ip,1)=0.0
          up(ip,2)=0.0
          up(ip,3)=0.0
        end do
	else

      T1 = TIMEF( )/1000.00

       call value(tur,up,1,lhnode,bg)
       call value(tvr,up,2,lhnode,bg)
       call value(twr,up,3,lhnode,bg)

      T2 = TIMEF( )/1000.00
      Tint = Tint + (T2-T1)
      
      endif

!   For a new run set the initial particle velocity

      if ( newpart ) then
        do  ip=1,npart
          vp(ip,1)=wp(ip)+up(ip,1)
          vp(ip,2)=up(ip,2)
          vp(ip,3)=up(ip,3)
        end do
      end if

! INCLUDE HERE THE CHANGE ON vp for ALREADY COLLIDED PARTICLES!!!!

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

! to get the relative velocity including hydrodynamic interactions

       if (HYDINT) then
         if ( newpart ) then
          do  ip=1,npart
          pertvel(ip,1)=0.0
          pertvel(ip,2)=0.0
          pertvel(ip,3)=0.0
          enddo
         endif
         call PERTURBVELOC(vp,yp,up,pertvel,nnodes)
       else
         do  ip=1,npart
          pertvel(ip,1)=0.0
          pertvel(ip,2)=0.0
          pertvel(ip,3)=0.0
         end do
       endif


	if (NOTURB) goto 802

!Timing
       T1 = TIMEF( )/1000.00

!
!   Calculate cross product: u x w in physical space
!   ------------------------------------------------
!   x-cpt  =  tvr*szr - twr*syr 
!   y-cpt  =  twr*sxr - tur*szr 
!   z-cpt  =  tur*syr - tvr*sxr 
!   then overwrite, saving these as sxr, syr, szr
!
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

      CALL SRCFT3(sxr,inc2x,inc3x,sxr,                        &
          inc2y,inc3y,n1,n2,n3,1,scaleRC,aux,naux)
      CALL SRCFT3(syr,inc2x,inc3x,syr,                        &
          inc2y,inc3y,n1,n2,n3,1,scaleRC,aux,naux)
      CALL SRCFT3(szr,inc2x,inc3x,szr,                        &
          inc2y,inc3y,n1,n2,n3,1,scaleRC,aux,naux)
!
!   filter the nonlinear term
!   -----------------------------
!
      do  k=1,n3
        do  j=1,n2
          do  i=1,n1pp
            sxr(i,j,k)=sxr(i,j,k)*filter(i,j,k)
            syr(i,j,k)=syr(i,j,k)*filter(i,j,k)
            szr(i,j,k)=szr(i,j,k)*filter(i,j,k)
          end do
        end do
      end do


       T2 = TIMEF( )/1000.00
       Tflow = Tflow + (T2-T1)


!====== P O S T P R O C E S S I N G =========================
!
!     compute the energy transfer spectrum
!
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
                 usqr=xmult*(cur(i,j,k)*sxr(i,j,k)+cvr(i,j,k)*syr(i,j,k)	&
                           +cwr(i,j,k)*szr(i,j,k) )
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
!
       go to 802
!   More optional data for the pressure.
!   -----------------------------------
!
      if(pflow)then
      if ( mod(istep,nshort).eq.0 ) then
!
!   Compute now q=(u*u)/2 from tur,tvr,twr still saved
!   in physical space, and save.
!
      do  k=1,n3
        do  j=1,n2
          do  i=1,n1
            tur(i,j,k)=0.5*(tur(i,j,k)**2 + tvr(i,j,k)**2 +		&
                           twr(i,j,k)**2)
          end do
        end do
      end do
!
!
!   Compute the modified pressure P from  div ( s ) in spectral form.
!
      do  k=1,n3
        do  j=1,n2
          tempk=fsq3(k)+fsq2(j)
          do  i=1,n1pp
            tvr(i,j,k)=(fac1(i)*sxr(i,j,k)+fac2(j)*syr(i,j,k)+fac3(k)*	&
                 szr(i,j,k))/(tempk+fsq1(i)+smallx)
          end do
        end do
      end do
!
!   Now transform to physical space and combine with q to get the
!   pressure field in physical space.Save P, and p in twr.
!
      CALL SCRFT3(tvr,inc2y,inc3y,tvr,                          &
          inc2x,inc3x,n1,n2,n3,-1,scaleCR,aux,naux)
      do  k=1,n3
        do  j=1,n2
          do  i=1,n1
            twr(i,j,k)=tvr(i,j,k)-tur(i,j,k)
          end do
        end do
      end do
!
      icpt=1
      call value(tur,pssp,icpt, lhnode,bg)
      jcpt=2
      call value(tvr,pssp,jcpt, lhnode,bg)
      kcpt=3
      call value(twr,pssp,kcpt, lhnode,bg)
!
      end if
      end if


802	continue
!
!   Short particle statistics
!
      if ( mod(istep,nshort).eq.0)		&
        call ptdat(numb,npset,ttime,yp,vp,up,pertvel,pssp,vortp,yp0)
!=========================================================================

       T1 = TIMEF( )/1000.00
!
!   Advance the particle velocity
!   ---------------------------
!    Equations of motion set up here are
!
!        dV/dt= (  u(Y,t) )/tau +  ( W - V )/tau
!
!        dY/dt= V(t)
!
!    The velocity is time-stepped first, the first group of terms is
!   advanced with an AB4 scheme while an AM4 scheme is used for the
!   second group. As W is fixed its contribution is simpler.
!    The position is advanced with an AM4 scheme afterwards.
!    The AB4 scheme has the format for  dz/dt=g  of:
!      z(n+1) = z(n) + (dt/24) * ( 55*g(n) - 59*g(n-1) + 37*g(n-2)
!                                  - 9*g(n-3) )
!    The AM4 scheme has the format 
!      z(n+1) = z(n) + (dt/24) * ( 9*g(n+1) + 19*g(n) - 5*g(n-1)
!                                    + g(n-2) )
!
      do  ic=1,3
        do ip=1,npart
          fvp(ip,ic,0)=vp(ip,ic)
          fdp(ip,ic,0)=up(ip,ic)+pertvel(ip,ic)
        end do

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
          do  ip=1,npart
            fvp(ip,ic,1)=fvp(ip,ic,0)
            fvp(ip,ic,2)=fvp(ip,ic,0)
            fdp(ip,ic,1)=fdp(ip,ic,0)
            fdp(ip,ic,2)=fdp(ip,ic,0)
            fdp(ip,ic,3)=fdp(ip,ic,0)
          end do
        end if
      end do
!
!   Set newpart=.false. now as all initialisation is done.
!
      newpart=.false.
!
!
! - write particle motion for later detailed analysis
!
!    Advance the particle velocity
!
        do  ip=1,npart
        Vrel=sqrt((vp(ip,1)-up(ip,1)-pertvel(ip,1))**2+		&
                 (vp(ip,2)-up(ip,2)-pertvel(ip,2))**2+		&
                 (vp(ip,3)-up(ip,3)-pertvel(ip,3))**2)

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
          tempv= 19.0*fvp(ip,ic,0) - 5.0*fvp(ip,ic,1)		&
                + fvp(ip,ic,2)
          tempd= 55.0*fdp(ip,ic,0) - 59.0*fdp(ip,ic,1)		&
                    +37.0*fdp(ip,ic,2) - 9.0*fdp(ip,ic,3)
          vp(ip,ic)=( vp(ip,ic) + dt*unitv*wp(ip)/tau(ip)*fSN	&
        - dtau*tempv + dtau*tempd)/( 1.0 + 9.0*dtau )
        end do
      end do
!
!   Advance the particle position
!   ---------------------------
!
      dtc=dt/24.0
!
      do  ic=1,3
        do  ip=1,npart
        dmove(ip,ic)=dtc*( 9.0*vp(ip,ic) +19.0*fvp(ip,ic,0)	&
                    -5.0*fvp(ip,ic,1) + fvp(ip,ic,2))
          yp(ip,ic)=yp(ip,ic)+dmove(ip,ic)
        end do
      end do

       T2 = TIMEF( )/1000.00
       Tadv = Tadv + (T2-T1)

       T1 = TIMEF( )/1000.00
!
!   Ensure that particles remain in box, use periodicity
!
      do  ip=1,npart
        if (  yp(ip,1).ge.hx ) yp(ip,1)=yp(ip,1)-hx
        if (  yp(ip,1).lt.0.0) yp(ip,1)=yp(ip,1)+hx
        if (  yp(ip,2).ge.hy ) yp(ip,2)=yp(ip,2)-hy
        if (  yp(ip,2).lt.0.0) yp(ip,2)=yp(ip,2)+hy
        if (  yp(ip,3).ge.hz ) yp(ip,3)=yp(ip,3)-hz
        if (  yp(ip,3).lt.0.0) yp(ip,3)=yp(ip,3)+hz
      end do

      T2 = TIMEF( )/1000.00
      Tpbc = Tpbc + (T2-T1)


!====================== P O S T P R O S E S S I N G===============================
! Save particle location
      if(mod(istep,20).eq.0)then

      num = int(istep/20)

      if (num.le.10) then
      i1d=num
      fvort = directpart2//'.'//char(i1d+48)
      endif

      if ((num.ge.10).and.(num.le.100)) then
      i1d=num/10
      i2d=num-i1d*10
      fvort = directpart2//'.'//char(i1d+48)//char(i2d+48)
      endif

      if ((num.ge.100).and.(num.le.1000)) then
      i1d=num/100
      i2d=num/10-i1d*10
      i3d=num-100*i1d-10*i2d
      fvort = directpart2//'.'//char(i1d+48)//char(i2d+48)//char(i3d+48)
      endif

      if ((num.ge.1000).and.(num.le.10000)) then
      i1d=num/1000
      i2d=num/100-10*i1d
      i3d=num/10-100*i1d-10*i2d
      i4d=num-1000*i1d-100*i2d-10*i3d
      fvort = directpart2//'.'//char(i1d+48)//char(i2d+48)    &
      //char(i3d+48)//char(i4d+48)
      endif

!      open( unit=60,file=fvort,status='unknown')
!      do ih=1,5000
!      write(60,729)yp(ih,1),yp(ih,2),yp(ih,3)
!      enddo
!      close(60)

      endif

!
! --    computing particle mean square velocity
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
      do ip=1,numb
      ip2=numb+ip
      vpmean1x=vpmean1x+vp(ip,1)
      vpmean1y=vpmean1y+vp(ip,2)
      vpmean1z=vpmean1z+vp(ip,3)
      vpmean2x=vpmean2x+vp(ip2,1)
      vpmean2y=vpmean2y+vp(ip2,2)
      vpmean2z=vpmean2z+vp(ip2,3)
      vpvar1x=vpvar1x+vp(ip,1)**2
      vpvar1y=vpvar1y+vp(ip,2)**2
      vpvar1z=vpvar1z+vp(ip,3)**2
      vpvar2x=vpvar2x+vp(ip2,1)**2
      vpvar2y=vpvar2y+vp(ip2,2)**2
      vpvar2z=vpvar2z+vp(ip2,3)**2
      enddo
      vpmean1x=vpmean1x/float(numb)
      vpmean1y=vpmean1y/float(numb)
      vpmean1z=vpmean1z/float(numb)
      vpmean2x=vpmean2x/float(numb)
      vpmean2y=vpmean2y/float(numb)
      vpmean2z=vpmean2z/float(numb)

!	write(106,*) 'vpvar1x,vpvar2x= ',vpvar1x,vpvar2x
!	write(106,*) 'vpvar1x/numb,vpvar2x/numb= ',		&
!      vpvar1x/float(numb),vpvar2x/float(numb)
!	write(106,*) 'vpmean1x,vpmean2x= ',vpmean1x,vpmean2x
!        write(106,*) 'vpmean1x**2,vpmean2x**2= ',vpmean1x**2,	&
!       vpmean2x**2

      vpvar1x=vpvar1x/float(numb)
      vpvar1y=vpvar1y/float(numb)
      vpvar1z=vpvar1z/float(numb)
      vpvar2x=vpvar2x/float(numb)
      vpvar2y=vpvar2y/float(numb)
      vpvar2z=vpvar2z/float(numb)

!	write(106,*) 'vpvar1x,vpvar2x= ',vpvar1x,vpvar2x

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


!Timing
       T1 = TIMEF( )/1000.00
!
!   Collision Detection
!   ********************
!
!   Indentifying the particle location relative to
!   the collision detection grid
!
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
        headd(i,j,k,my_thread)=0
        enddo
        enddo
        enddo

      my_thread=1
      do ip=1,npart
      ix=1+int(ypp(ip,1)/wcd)
      iy=1+int(ypp(ip,2)/wcd)
      iz=1+int(ypp(ip,3)/wcd)
!
      if(ix.gt.ncd ) then
      write(103,*)ix,iy,iz,ip
      write(103,*)ttime,ypp(ip,1),ypp(ip,2),	&
        ypp(ip,3),wcd,hx,hy,hz,'COLLDET'
      ix = ncd
      endif

      if(iy.gt.ncd)then
      write(103,*)ix,iy,iz,ip
      write(103,*)ttime,ypp(ip,1),ypp(ip,2),	&
        ypp(ip,3),wcd,hx,hy,hz,'COLLDET'
      iy = ncd
      endif

      if(iz.gt.ncd)then
      write(103,*)ix,iy,iz,ip
      write(103,*)ttime,ypp(ip,1),ypp(ip,2),	&
        ypp(ip,3),wcd,hx,hy,hz,'COLLDET'
      iz = ncd
      endif
!
        list(ip)=headd(ix,iy,iz,my_thread)
        headd(ix,iy,iz,my_thread)=ip
      end do

        do l1=1,ncd
        do l2=1,ncd
        do l3=1,ncd
        head(l1,l2,l3)=headd(l1,l2,l3,1)
	enddo
        enddo
        enddo

      do ip=1,npart
      dmoved(ip)=sqrt(dmove(ip,1)**2+dmove(ip,2)**2	&
          +dmove(ip,3)**2)
      end do

        do i=0,179
        nir11n(i,my_thread)=0
        nir12n(i,my_thread)=0
        nir22n(i,my_thread)=0

        wrbin11n(i,my_thread)=0.0
        wrbin12n(i,my_thread)=0.0
        wrbin22n(i,my_thread)=0.0
        enddo

        DO IP1=1, NPART

	ncollip1=0

         ix=1+int(ypp(ip1,1)/wcd)
         iy=1+int(ypp(ip1,2)/wcd)
         iz=1+int(ypp(ip1,3)/wcd)

      if(ix.gt.ncd ) then
      write(103,*)ix,iy,iz,ip1
      write(103,*)ttime,ypp(ip1,1),ypp(ip1,2),		&
       ypp(ip1,3),wcd,hx,hy,hz
      ix = ncd
      endif

      if(iy.gt.ncd)then
      write(103,*)ix,iy,iz,ip1
      write(103,*)ttime,ypp(ip1,1),ypp(ip1,2),		&
        ypp(ip1,3),wcd,hx,hy,hz
      iy = ncd
      endif

      if(iz.gt.ncd)then
      write(103,*)ix,iy,iz,ip1
      write(103,*)ttime,ypp(ip1,1),ypp(ip1,2),		&
        ypp(ip1,3),wcd,hx,hy,hz
      iz = ncd
      endif
!
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
!
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
!
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
!
!	ip2 has to be always greater than ip1 to avoid double counting
!
        if(ip2.eq.0)goto 200
	if(ip2.le.ip1)goto 91

!       if(ip1.eq.5720) then
!        ijole=ip2
!528     write(525,*) ttime,istep,ijole,ix,iy,iz,ip2,ix2,iy2,iz2
!        ijole=list(ijole)
!        if(ijole.eq.0) goto 199
!        goto 528
!       endif

!
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

!-add bin statistics
!      if(dnij1WRGR.lt.(10.*rcoll))then
       if(dnij1WRGR.ge.rcoll .and. dnij1WRGR.le.(10.*rcoll))then
!       if(dnij1WRGR.ge.(10.*rcoll) .and. dnij1WRGR.le.(19.*rcoll))then
! modified on Oct. 9, 2003
        xr=(xa1-xb1)/dnij1WRGR
        yr=(xa2-xb2)/dnij1WRGR
        zr=(xa3-xb3)/dnij1WRGR
        Vr1=va1*xr+va2*yr+va3*zr
        Vr2=vb1*xr+vb2*yr+vb3*zr
        Wr=Vr1-Vr2

        ir=INT( (dnij1WRGR-rcoll) /(0.05*rcoll) )
!        ir=INT( (dnij1WRGR-10.0*rcoll) /(0.05*rcoll) )

        if(ir.gt.179)ir=179
        if(ir.lt.0)ir=0
!
        if(ip1.le.numb) then
                if(ip2.le.numb) then
        nir11n(ir,my_thread)=nir11n(ir,my_thread)+1
        wrbin11n(ir,my_thread)=wrbin11n(ir,my_thread)+abs(Wr)
                else
!C I made changes here and below !!!
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

!       if(dnij1WRGR.ge.rcoll .and. dnij1WRGR.le.rcollb )then
!       write(110,*) ttime, 'ip1= ',ip1,' ip2= ',ip2 
!       endif

       if(dnij1WRGR.ge.rcolla .and. dnij1WRGR.le.rcollb )then
! modified on Oct. 22, 1997
        xr=(xa1-xb1)/dnij1WRGR
        yr=(xa2-xb2)/dnij1WRGR
        zr=(xa3-xb3)/dnij1WRGR
        Vr1=va1*xr+va2*yr+va3*zr
        Vr2=vb1*xr+vb2*yr+vb3*zr
        Wr=Vr1-Vr2
!
        if(ip1.le.numb) then
                if(ip2.le.numb) then
                  nR11e=nR11e+1
                  write(81,191)ttime,va1,va2,va3,vb1,vb2,vb3
                  write(71,191)ttime,xr,yr,zr,Wr
                else
                  nR12e=nR12e+1
                  write(82,191)ttime,va1,va2,va3,vb1,vb2,vb3
                  write(72,191)ttime,xr,yr,zr,Wr
                endif
        endif
        if(ip1.gt.numb) then
                if(ip2.gt.numb) then
                  nR22e=nR22e+1
                  write(83,191)ttime,va1,va2,va3,vb1,vb2,vb3
                  write(73,191)ttime,xr,yr,zr,Wr
                else
                  nR12e=nR12e+1
                  write(82,191)ttime,va1,va2,va3,vb1,vb2,vb3
                  write(72,191)ttime,xr,yr,zr,Wr
                endif
        endif

! This part was added by Bogdan Rosa Mar. 02, 2007
! It allows calculate RDF as a function of angle
! Important: ip1 must be bigger than ip2!

       if(dnij1WRGR.ge.rcolla .and. dnij1WRGR.le.rcollb )then
         if(ip1.ge.numb .and. ip2.le.numb)then

          alpha=INT(acos(-(xb1-xa1)/dnij1WRGR))
          pdf(alpha)=pdf(alpha)+1

         endif
        endif

!----------------------------------------------------------


191   format(f10.6,6(1x,f11.5))
       endif

 
       yppb1=ypp(ip2,1)+sx*hx
       yppb2=ypp(ip2,2)+sy*hy
       yppb3=ypp(ip2,3)+sz*hz
!
       ypb1=yppb1+dmove(ip2,1)
       ypb2=yppb2+dmove(ip2,2)
       ypb3=yppb3+dmove(ip2,3)
!---   coarse check
       dnij0=sqrt((ypp(ip1,1)-yppb1)**2			&
      +(ypp(ip1,2)-yppb2)**2+(ypp(ip1,3)-yppb3)**2)
       tmove=rcoll+dmoved(ip1)+dmoved(ip2)
       if(dnij0.gt.tmove)goto 91
       dnij1=sqrt((ypa1-ypb1)**2+(ypa2-ypb2)**2+(ypa3-ypb3)**2)

! Searching for type III collisions
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

       xa=x01+v01*tfac1+(3.*x11-(2.*v01+v11)*dt)*tfac2		&
      +((v11+v01)*dt-2.*x11)*tfac3
       ya=x02+v02*tfac1+(3.*x12-(2.*v02+v12)*dt)*tfac2		&
      +((v12+v02)*dt-2.*x12)*tfac3
       za=x03+v03*tfac1+(3.*x13-(2.*v03+v13)*dt)*tfac2		&
      +((v13+v03)*dt-2.*x13)*tfac3

       xb=bx01+bv01*tfac1+(3.*bx11-(2.*bv01+bv11)*dt)*tfac2	&
      +((bv11+bv01)*dt-2.*bx11)*tfac3
       yb=bx02+bv02*tfac1+(3.*bx12-(2.*bv02+bv12)*dt)*tfac2	&
      +((bv12+bv02)*dt-2.*bx12)*tfac3
       zb=bx03+bv03*tfac1+(3.*bx13-(2.*bv03+bv13)*dt)*tfac2	&
      +((bv13+bv03)*dt-2.*bx13)*tfac3

       dnij=sqrt( (xa-xb)**2+(ya-yb)**2+(za-zb)**2 )

       if(dnij.gt.rcoll )then
       icc3e=icc3e+1
       ncollip1=ncollip1+1
       goto 92
       endif

        enddo

       endif

        if(OVERL)then
        continue
        else
        write(*,*) 'WARNING --> TYPE III COLLISION ',ip1,ip2,istep
!      write(325,*) ttime,yp(ip1,1),yp(ip1,2),yp(ip1,3),yp(ip2,1),	&
!                        yp(ip2,2),yp(ip2,3),istep
!      write(325,*) dnij0,dnij1,dnij,rcoll

        goto 92
        endif

       goto 91

       endif
!
! searching for type I collisions
!       if(dnij0.gt.rcoll) then
       if(dnij1.le.rcoll) then
	 icc1e=icc1e+1
	 ncollip1=ncollip1+1
         write(*,*)'type I collisions'
         goto 92
       endif
!       endif

! searching for type  II collisions

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

       xa=x01+v01*tfac1+(3.*x11-(2.*v01+v11)*dt)*tfac2		&
      +((v11+v01)*dt-2.*x11)*tfac3
       ya=x02+v02*tfac1+(3.*x12-(2.*v02+v12)*dt)*tfac2		&
      +((v12+v02)*dt-2.*x12)*tfac3
       za=x03+v03*tfac1+(3.*x13-(2.*v03+v13)*dt)*tfac2		&
      +((v13+v03)*dt-2.*x13)*tfac3

       xb=bx01+bv01*tfac1+(3.*bx11-(2.*bv01+bv11)*dt)*tfac2	&
      +((bv11+bv01)*dt-2.*bx11)*tfac3
       yb=bx02+bv02*tfac1+(3.*bx12-(2.*bv02+bv12)*dt)*tfac2	&
      +((bv12+bv02)*dt-2.*bx12)*tfac3
       zb=bx03+bv03*tfac1+(3.*bx13-(2.*bv03+bv13)*dt)*tfac2	&
      +((bv13+bv03)*dt-2.*bx13)*tfac3

       dnij=sqrt( (xa-xb)**2+(ya-yb)**2+(za-zb)**2 )
       if(dnij.lt.rcoll )then
       icc2e=icc2e+1
       ncollip1=ncollip1+1
       write(*,*)'type II collisions'
!
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

!       write(109,*) ttime, 'ip1= ',ip1,' ip2= ',ip2

! This section is NEW. We want to get the approaching angle when
! particles are colliding

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


        ncollpart=ncollpart+1
        removpart(ncollpart)=ip2

!  We assume that only one collision is possible for a given particle
!     goto 99   WE WANT TO CHECK EVEN IF MORE THAN ONE COLLISION 
!
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
        ncollpart=ncollpart+1
        removpart(ncollpart)=ip1
       endif

      ENDDO
!
      do i=0,179
      nir11(i) = nir11(i) + nir11n(i,1)
      nir12(i) = nir12(i) + nir12n(i,1)
      nir22(i) = nir22(i) + nir22n(i,1)

      wrbin11(i) = wrbin11(i) + wrbin11n(i,1)
      wrbin12(i) = wrbin12(i) + wrbin12n(i,1)
      wrbin22(i) = wrbin22(i) + wrbin22n(i,1)
      enddo
!
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

! Relocate particles that just collided

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


      T2 = TIMEF( )/1000.00
      Tcd = Tcd + (T2-T1)

! Check for Overlapping of the newly relocated particles

       T1 = TIMEF( )/1000.00

	if(OVERL) goto 806
       IF(newcoltot.ne.0) THEN
!
!   Overlapping Detection and Correction
!   ***********************************
!
!   Indentifying the particle location relative to
!   the ovarelapping detection grid
!
1982    continue
!
        do i=1,ncd
        do j=1,ncd
        do k=1,ncd
        headd(i,j,k,my_thread)=0
        enddo
        enddo
        enddo
!
      my_thread=1
      do ip=1,npart
      ix=1+int(yp(ip,1)/wcd)
      iy=1+int(yp(ip,2)/wcd)
      iz=1+int(yp(ip,3)/wcd)
!
      if(ix.gt.ncd ) then
      ix = ncd
      endif

      if(iy.gt.ncd)then
      iy = ncd
      endif

      if(iz.gt.ncd)then
      iz = ncd
      endif
!
        list(ip)=headd(ix,iy,iz,my_thread)
        headd(ix,iy,iz,my_thread)=ip
      end do

        do l1=1,ncd
        do l2=1,ncd
        do l3=1,ncd
        head(l1,l2,l3)=headd(l1,l2,l3,1)
        enddo
        enddo
        enddo
!
        noverl=0
!
        DO IPOS=1, NCOLLPART

          ip1=removpart(ipos)

         ix=1+int(yp(ip1,1)/wcd)
         iy=1+int(yp(ip1,2)/wcd)
         iz=1+int(yp(ip1,3)/wcd)
!
      if(ix.gt.ncd ) then
      write(103,*)ix,iy,iz,ip1
      write(103,*)ttime,yp(ip1,1),yp(ip1,2),	&
        yp(ip1,3),wcd,hx,hy,hz
      ix = ncd
      endif

      if(iy.gt.ncd)then
      write(103,*)ix,iy,iz,ip1
      write(103,*)ttime,yp(ip1,1),yp(ip1,2),	&
        yp(ip1,3),wcd,hx,hy,hz
      iy = ncd
      endif

      if(iz.gt.ncd)then
      write(103,*)ix,iy,iz,ip1
      write(103,*)ttime,yp(ip1,1),yp(ip1,2),	&
        yp(ip1,3),wcd,hx,hy,hz
      iz = ncd
      endif
!
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
!
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
!
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
!
!
        if(ip2.eq.0)goto 2002
        if(ip2.eq.ip1)goto 912
!
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

       dnij=sqrt((yp(ip1,1)-ypb1)**2			&
      +(yp(ip1,2)-ypb2)**2+(yp(ip1,3)-ypb3)**2)

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

!        write(*,*) 'noverlAFTERCOLL= ',noverl

        if(noverl.ne.0) goto 1982
!
!  END OF DETECTION OF OVERLAPPING
!

       ENDIF
806	continue 
!
!Timing
       T2 = TIMEF( )/1000.00
       Tod =Tod+(T2-T1)


!Timing
       T1 = TIMEF( )/1000.00

	if (NOTURB) goto 803
!-Change 6/31/95
!
!   ADVANCE THE FLOW FIELD


!---------------- Exact integration scheme starts ----------------

      if (integration_scheme.eq.2) then

      if(newflo.eqv..false.) goto 3000

      if(istep.eq.0) then

        do  k=1,n3
         do  i=1,n1pp
          do  j=1,n2
             k2 = fac1(i)**2 + fac2(j)**2 + fac3(k)**2
             cur(i,j,k) = cur(i,j,k) * (1.-rnuf*k2) + dtp * sxr(i,j,k)
             cvr(i,j,k) = cvr(i,j,k) * (1.-rnuf*k2) + dtp * syr(i,j,k)
             cwr(i,j,k) = cwr(i,j,k) * (1.-rnuf*k2) + dtp * szr(i,j,k)
          enddo
         enddo
        enddo

        if (forcing.eq.2) call supfordet(cur,cvr,cwr,fac1,fac2,fac3,force)
        call projection(cur,cvr,cwr,fac1,fac2,fac3)

        sqx (:,:,:) = sxr (:,:,:)
        sqy (:,:,:) = syr (:,:,:)
        sqz (:,:,:) = szr (:,:,:)

        ttime = ttime + dtp
        go to 5643

      end if


      if(istep.eq.1) then

        do  k=1,n3
         do  i=1,n1pp
          do  j=1,n2
             k2 = fac1(i)**2 + fac2(j)**2 + fac3(k)**2
             cur(i,j,k) = sxp(i,j,k) + dt * sxr(i,j,k) - rnu * dt * k2 * cur(i,j,k)
             cvr(i,j,k) = syp(i,j,k) + dt * syr(i,j,k) - rnu * dt * k2 * cvr(i,j,k)
             cwr(i,j,k) = szp(i,j,k) + dt * szr(i,j,k) - rnu * dt * k2 * cwr(i,j,k)
          enddo
         enddo
        enddo

       call projection(cur,cvr,cwr,fac1,fac2,fac3)
       if (forcing.eq.1) call supforstoch(cur,cvr,cwr,iseedf,iyf,ivf,dt)
       if (forcing.eq.2) call supfordet(cur,cvr,cwr,fac1,fac2,fac3,force)

       sxp(:,:,:)= sqx (:,:,:)
       syp(:,:,:)= sqy (:,:,:)
       szp(:,:,:)= sqz (:,:,:)

      ttime = ttime + dtp
      go to 5643

      endif

3000  CONTINUE

      do  k=1,n3
        do  i=1,n1pp
          do  j=1,n2

             k2 = fac1(i)**2 + fac2(j)**2 + fac3(k)**2
             k2_e = exp(-k2*dt*rnu)

             if (k2.lt.0.5) then
             k2 = 1.e-5
             k2_e = 1.0
             endif

             cur(i,j,k) = cur(i,j,k) + dtp * ( 3. * sxr(i,j,k) - k2_e * sxp(i,j,k) )
             cvr(i,j,k) = cvr(i,j,k) + dtp * ( 3. * syr(i,j,k) - k2_e * syp(i,j,k) )
             cwr(i,j,k) = cwr(i,j,k) + dtp * ( 3. * szr(i,j,k) - k2_e * szp(i,j,k) )

          enddo
        enddo
      enddo

      if (forcing.eq.1) call supforstoch(cur,cvr,cwr,iseedf,iyf,ivf,dt)
      call projection(cur,cvr,cwr,fac1,fac2,fac3)

      do  k=1,n3
        do  i=1,n1pp
          do  j=1,n2

             k2 = fac1(i)**2 + fac2(j)**2 + fac3(k)**2
             k2_e = exp(-k2*dt*rnu)

             if (k2.lt.0.5) then
              k2 = 1.e-5
              k2_e = 1.0
             endif

             cur(i,j,k) = cur(i,j,k) * k2_e
             cvr(i,j,k) = cvr(i,j,k) * k2_e
             cwr(i,j,k) = cwr(i,j,k) * k2_e

          enddo
        enddo
      enddo

      if (forcing.eq.2) call supfordet(cur,cvr,cwr,fac1,fac2,fac3,force)

       sxp = sxr
       syp = syr
       szp = szr

      endif
      ttime = ttime + dt

5643  continue

      call clean(cur)
      call clean(cvr)
      call clean(cwr)

!   Force a zero mean flow
!   ----------------------
!
      cur(1,1,1)=0.0
      cur(2,1,1)=0.0
      cvr(1,1,1)=0.0
      cvr(2,1,1)=0.0
      cwr(1,1,1)=0.0
      cwr(2,1,1)=0.0

!   Use the filter here if needed
!   -----------------------------
!
      do  k=1,n3
        do  j=1,n2
          do  i=1,n1pp
            cur(i,j,k)=cur(i,j,k)*filter(i,j,k)
            cvr(i,j,k)=cvr(i,j,k)*filter(i,j,k)
            cwr(i,j,k)=cwr(i,j,k)*filter(i,j,k)
          end do
        end do
      end do

!------------ End of Exact integration scheme ------------


!   If this is a new run replace AB2 by an Euler step
!   -------------------------------------------------

      if (integration_scheme.eq.1) then
!
      if ( newflo ) then
!
!$omp parallel do private(j,i)
      do  k=1,n3
        do j=1,n2
          do i=1,n1pp
            sxp(i,j,k)=sxr(i,j,k)
            syp(i,j,k)=syr(i,j,k)
            szp(i,j,k)=szr(i,j,k)
          end do
        end do
      end do
!$omp end parallel do
!
      newflo=.false.
!
      endif
!
!   Advance the flow field - 1st part non-linear terms advanced
!   -----------------------------------------------------------
!                  2
!              nu k  dt  ^      3 dt     ^           dt      ^
!    bu = (1 - --------) u    + ---- [   s   ]   - ---- [    s   ]
!                 2       (n)     2           1(n)    2           1(n-1)
!
!                  2
!              nu k  dt  ^      3 dt     ^           dt      ^
!    bv = (1 - --------) v    + ---- [   s   ]   - ---- [    s   ]
!                 2       (n)     2           2(n)    2           2(n-1)
!
!
!                  2
!              nu k  dt  ^      3 dt     ^           dt      ^
!    bw = (1 - --------) w    + ---- [   s   ]   - ---- [   s   ]
!                 2       (n)     2           3(n)    2           3(n-1)
!
!$omp parallel do private(j,i,facsq,temp)
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
!$omp end parallel do
!
!   Add the forcing
!   ----------------------
!   Note that we could have combined the forcing a1r,a2,a3r with
!   the terms sxr,syr,szr at each step ( and so sxp,etc. ), but
!   as this is random forcing this Euler step is as effective.
!
      if ( forced ) then

      if ( forcing.eq.1 ) call supforstoch (cur,cvr,cwr,iseedf,iyf,ivf,dt)
      if ( forcing.eq.2 ) call supfordet(cur,cvr,cwr,fac1,fac2,fac3,force)

      endif
!
!   Advance the flow field   -  linear terms
!   ----------------------------------------
!   The presure term is eliminated by continuity eqn. Note
!   the pressure here is a pseudo-pressure not the real one.
!   To get the pressure explicitly one must solve a Poisson eq.
!   as done in the earlier data step.
!
!
!   Finally the velocity field in spectral form at (n+1) step is
!                                        2
!       ^              bu - k1 ( k.bu )/k
!       u       = ----------------------
!        (n+1)                  2
!                  (1 + 0.5 nu k  dt)
!
!   etc.
!
!
!$omp parallel do private(j,i,temp,temp1,temp2,pres)
      do  k=1,n3
        do  j=1,n2
          temp=fsq3(k)+fsq2(j)+smallx
          do  i=1,n1pp
            temp1=1.0/(1.0+rnuf*(fsq1(i)+temp))
            temp2=1.0/(fsq1(i)+temp)
            pres=temp2*(fac1(i)*cur(i,j,k)+fac2(j)*cvr(i,j,k)   &
                +fac3(k)*cwr(i,j,k))
            cur(i,j,k)=temp1*(cur(i,j,k)-fac1(i)*pres)
            cvr(i,j,k)=temp1*(cvr(i,j,k)-fac2(j)*pres)
            cwr(i,j,k)=temp1*(cwr(i,j,k)-fac3(k)*pres)
          end do
        end do
      end do
!$omp end parallel do

!
!   Clean and filter data
!   ---------------------
!   Enforce symmetry of fourier coefficients in k_1=0 plane
!   and truncate the n/2 wavenos.
!
      call clean(cur)
      call clean(cvr)
      call clean(cwr)
!   Force a zero mean flow
!   ----------------------
!
      cur(1,1,1)=0.0
      cur(2,1,1)=0.0
      cvr(1,1,1)=0.0
      cvr(2,1,1)=0.0
      cwr(1,1,1)=0.0
      cwr(2,1,1)=0.0
!
!   Use the filter here if needed
!   -----------------------------
!
!$omp parallel do private(j,i)
      do  k=1,n3
        do  j=1,n2
          do  i=1,n1pp
            cur(i,j,k)=cur(i,j,k)*filter(i,j,k)
            cvr(i,j,k)=cvr(i,j,k)*filter(i,j,k)
            cwr(i,j,k)=cwr(i,j,k)*filter(i,j,k)
          end do
        end do
      end do
!$omp end parallel do
!
!   Update (u x w + f )
!   -------------------
!
!$omp parallel do private(j,i)
      do  k=1,n3
        do   j=1,n2
          do  i =1,n1pp
            sxp(i,j,k)=sxr(i,j,k)
            syp(i,j,k)=syr(i,j,k)
            szp(i,j,k)=szr(i,j,k)
          end do
        end do
      end do
!$omp end parallel do
!Timing

     ttime=ttime+dt
     endif

!!!!!!!!!!!! End of CN integration scheme !!!!!!!!!!!!!!!!!!!

803	continue
       T2 = TIMEF( )/1000.00
       Tflow =Tflow+ (T2-T1)

!Timing
       T1 = TIMEF( )/1000.00

!
!   Update particle history
!   ---------------------
!
      do  ic=1,3
        do  ip=1,npart
          fvp(ip,ic,2)=fvp(ip,ic,1)
          fvp(ip,ic,1)=fvp(ip,ic,0)
          fdp(ip,ic,3)=fdp(ip,ic,2)
          fdp(ip,ic,2)=fdp(ip,ic,1)
          fdp(ip,ic,1)=fdp(ip,ic,0)
          ypp(ip,ic)=yp(ip,ic)
        end do
      end do
!
!Timing
       T2 = TIMEF( )/1000.00
       Tadv = Tadv+(T2-T1)


!       write(*,*)'Tflow1,Tflow2,Tflow3,Tpart1,Tpart2,Tpart3,Tpart4,Tcollision,Tpertvel,TSTEP,TCHECK'
!       write(*,501)Tflow1,Tflow2,Tflow3,Tpart1,Tpart2,		&
!       Tpart3,Tpart4,Tcollision,Tpertvel,TSTEP,TCHECK
!       write(*,501)Tflow1/tstep,Tflow2/tstep,Tflow3/tstep,		&
!       Tpart1/tstep,Tpart2/tstep,Tpart3/tstep,				&
!       Tpart4/tstep,Tcollision/tstep,Tpertvel/tstep,TCHECK/tstep
!501    format(11f10.6)
!
!
!   Advance time
!   ------------
!
      istep=istep+1
!
      if(mod(istep,nshort).eq.0) then
       write(80,303)ttime,nR11,nR12,nR22,vpmean1x,vpmean1y,		&
                 vpmean1z,vpmean2x,vpmean2y,vpmean2z

       write(8,301)ttime,numb,itocol11,itocol12,itocol22,		&
         vpmean1x,vpmean1y,vpmean1z,vpmean2x,vpmean2y,vpmean2z,		&
         icc1,icc2,icc3
!       write(*,301)ttime,numb,itocol11,itocol12,itocol22,		&
!         vpmean1x,vpmean1y,vpmean1z,vpmean2x,vpmean2y,vpmean2z,	&
!         icc1,icc2,icc3
	write(106,302) ttime,vpmean1x,vpmean1y,vpmean1z,vpmean2x,vpmean2y,	&
        vpmean2z,vpvar1x,vpvar1y,vpvar1z,vpvar2x,vpvar2y,vpvar2z
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
!
!   End of time-step
!   ----------------
!
       write(*,*)'ending istep=',istep

      if (istep.lt.nhalt) go to 100



      write(88,*)Tint,Tflow,Tadv,Tpbc,Tcd,Tod

!
!   Format statements for unit=16
!   -----------------------------
!
 101  format( 2x, i6, 9(1pe12.4) )
 102  format(2x,i6,2(2x,1pe12.4) )
 103  format( 2(2x,1pe12.4))
 105  format(2x,13(1pe12.4) )
 1001 format(//' Start of new run with particle code ')
 1010 format(/'  n1,n2,n3 =',3(3x,i5))
 1015 format(/'  parameters rnu,dt,var,tf=',/4(3x,1pe12.4))
 1020 format(/'  time at start of run =',3x,e17.10)
 1025 format(/'  logical values of newflo,newpart,forced are',		&
      /3(3x,l7))
 1030 format(/'  nhalt,nlong,nshort=',3(3x,i5))
 1035 format(/'  iseedp,iseedf =',2(3x,I20))
!
 1080 format(///' npart,nset, numb =',3(3x,i8))
 1100 format(/'   error in setting npset')
 1110 format(//'  particle parameters used are tau,wp:')
707   format(2x,8(1pe9.2))
901   format(2x,i6,6f8.3)
178   format(2x,'particle leaving box = ',I7,f12.5)
271    format(2x,2f8.4)
729   format(2x,3(1pe12.2))
!
!
!   End of run; dump data for a restart
!   -----------------------------------
!

      do alpha=0,360
        write(324,*)alpha,pdf(alpha)
      enddo

      open ( unit=12, file=directvel, 	&
      form='unformatted',status='unknown' )
      write (12) cur,cvr,cwr
      write (12) sxp,syp,szp
!RANDOM
      write (12)iseedf,iyf,ivf
      write (12) a1r,a2r,a3r
      write (12) b1r,b2r,b3r
      close(12)
!
      open ( unit=11,file=directpart, 	&
      form='unformatted', status='unknown' )
      write (11) yp,vp,fvp,fdp
!RANDOM
      write (11)iseedp,iyp,ivp
      close(11)

      deallocate(slab)

!
      stop
      end 

