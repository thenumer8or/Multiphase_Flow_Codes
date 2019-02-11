!==================================================================

      SUBROUTINE gaussian(v)
      use var_inc
      use IFPORT
      implicit none

      real, dimension(lx+2,lly,lz)    :: v 
      real t1, t2
      integer i, j, k

      v = 0.0

      do i = 1,lx
        do j = 1,lly,2
          do k = 1,lz
            t1 = rand()
            t2 = rand()
            if(t1 <= 1.e-10) t1 = 1.e-10
            if(t2 <= 1.e-10) t2 = 1.e-10
            t2 = pi2*t2
            v(i,j,k) = dsqrt(-2.d0*log(t1))*cos(t2)
            v(i,j+1,k) = dsqrt(-2.d0*log(t1))*sin(t2)
          end do
        end do
      end do

      END SUBROUTINE gaussian
!==================================================================

      SUBROUTINE wave
      use var_inc
      implicit none

      INTEGER                       :: i,j,k,ig,jg

      do i = 2,(ly+lyext),2
        ig = (ly*indy) + i
        kx(:,i-1,:) = real((ig/2) - 1)
        kx(:,i,:) = real((ig/2) - 1)
      enddo

      do j = 1, lz
       jg = (lz*indz)+j
       if (jg.lt.ny/2+2) then
        ky(:,:,j) = real(jg - 1)
       else
        ky(:,:,j) = real(-(ny+1-jg))
       endif
      enddo

      do k = 1, lx
       if (k.lt.nz/2+2) then
        kz(k,:,:) = real(k - 1)
       else
        kz(k,:,:) = real(-(nz+1-k))
       endif
      enddo

      k2   = kx*kx + ky*ky + kz*kz
      if ( myid.eq.0 ) then
        k2(1,1:2,1)   = 1.e-5
      endif

! Putting zeros (or ones for k2) at the additional padding on the x_matrix 
! direction. This padding is needed for the 1st step FFT calculation but
! not on the spectral space.
      kx(lx+1:lx+2,:,:) = 0.0
      ky(lx+1:lx+2,:,:) = 0.0
      kz(lx+1:lx+2,:,:) = 0.0
      k2(lx+1:lx+2,:,:) = 1.0

      ik2 = int(sqrt(k2) + 0.5)

      RETURN 

      END SUBROUTINE wave
!==================================================================

      SUBROUTINE symmetrize(c)
! For this subroutine, I'm assuming 2^n processors along Z_matrix direction
      use mpi
      use var_inc
      implicit none

      REAL,    DIMENSION (lx+2,lly,lz)   :: c
      REAL,    DIMENSION (lx+2,2,2:lz)   :: ccS, ccR
      REAL,    DIMENSION (lx+2,2)        :: ccvS, ccvR
      INTEGER                           :: npSR, isize,llx
      INTEGER                           :: i, ii, k, kk
      INTEGER                           :: status(MPI_STATUS_SIZE,2),req(2)

      llx = lx+2
! Force a zero mean flow
      if (myid.eq.0) c(1,1:2,1) = 0.0


! Truncate the n/2 wave number as this is noise usually
! For kx = nx/2
      if(indy.eq.(nprocY-1)) c(:,lly-1:lly,:) = 0.0
! For ky = ny/2
      if(indz.eq.(nprocZ/2)) c(:,:,1) = 0.0
! For kz = nz/2
      c(nz/2+1,:,:) = 0.0


! Round and other errors are corrected to ensure the physical data is real
! valued. For that, symmetry must be checked and enforced at the plane kx=0. 
! For the other planes, kx>0, the symmetry is warrantied by the FFT.
      IF(indy.eq.0) THEN
        IF(nproc.ne.1) THEN

          isize = 2*llx*(lz-1)
          ccS = c(:,1:2,2:lz)
          npSR = nprocZ - indz - 1 
          npSR = npSR * nprocY
          call mpi_isend(ccS,isize,MPI_REAL8,npSR,1,MPI_COMM_WORLD, &
                         req(1),ierr)
          call mpi_irecv(ccR,isize,MPI_REAL8,npSR,1,MPI_COMM_WORLD, &
                         req(2),ierr)
          call mpi_waitall(2,req,status,ierr)
          do i=1,lx
            ii=nz+2-i
            if (i.eq.1) ii=1
            do k=2,lz
              kk = lz+2-k
              c(i,1,k) = 0.5*(c(i,1,k) + ccR(ii,1,kk))
              c(i,2,k) = 0.5*(c(i,2,k) - ccR(ii,2,kk))
            enddo
          enddo

          if ((indz.ne.0).and.(indz.ne.nprocZ/2)) then
            isize = llx*2
            ccvS = c(:,1:2,1)
            npSR = nprocZ - indz
            npSR = npSR * nprocY
            call mpi_isend(ccvS,isize,MPI_REAL8,npSR,2,MPI_COMM_WORLD, &
                         req(1),ierr)
            call mpi_irecv(ccvR,isize,MPI_REAL8,npSR,2,MPI_COMM_WORLD, &
                         req(2),ierr)
            call mpi_waitall(2,req,status,ierr)

            do i=1,lx
              ii=nz+2-i
              if (i.eq.1) ii=1
              c(i,1,1) = 0.5*(c(i,1,1) + ccvR(ii,1))
              c(i,2,1) = 0.5*(c(i,2,1) - ccvR(ii,2))
            enddo
          endif

        ELSE

          do i=1,lx
            ii=nz+2-i
            if (i.eq.1) ii=1
            do k=2,lz/2
              kk = ny+2-k
              c(i,1,k) = 0.5*(c(i,1,k) + c(ii,1,kk))
              c(i,2,k) = 0.5*(c(i,2,k) - c(ii,2,kk))
              c(ii,1,kk) =  c(i,1,k)
              c(ii,2,kk) = -c(i,2,k)
            enddo
          enddo

        ENDIF

          if (indz.eq.0) then
            do i=2,lx/2
              ii=nz+2-i
              c(i,1,1) = 0.5*(c(i,1,1) + c(ii,1,1))
              c(i,2,1) = 0.5*(c(i,2,1) - c(ii,2,1))
              c(ii,1,1) =  c(i,1,1)
              c(ii,2,1) = -c(i,2,1)
            enddo
          endif

      ENDIF

      RETURN
      END SUBROUTINE symmetrize
!==================================================================
      SUBROUTINE FORCING
      use mpi
      use var_inc
      implicit none

      REAL      :: fr(6,5,5,3)
      INTEGER   :: mc, k, j, i, ip, jj, kk
      REAL      :: cont1, rmag, rph, divb, radm
      REAL      :: smallx, xkf, var, tf, dt
      REAL      :: t1, t2, t3
      REAL      :: ranff

      smallx = 1.0e-18
      xkf = sqrt(8.0)
      var = 447.31*(vscale*tscale)**2
      tf  = 0.038/tscale
      dt = 1.0

!     Generate white noise
      do mc=1,3
        do  k=1,5
          do  j=1,5
            do  i=1,6,2
              ip=i+1
              cont1 = ranff(iseedf,iyf,ivf)
              if ( cont1.lt.1.e-6 ) cont1 = 1.e-6
              rmag = sqrt ( -4.0 * var * dt * alog(cont1) / tf )
              rph =  ranff (iseedf,iyf,ivf)
              fr(i ,j,k,mc) = rmag * cos(pi2*rph)
              fr(ip,j,k,mc) = rmag * sin(pi2*rph)
            end do
          end do
        end do
      end do

!     Overwrite components for k1=0 to ensure real forcing
!     in physical space with no spatial mean.

      do mc = 1, 3
        do k = 1, 5
          kk = 7 - k
          if ( k.eq.1 ) kk = 1
          do j = 2, 3
            jj = 7 - j
            fr(1,jj,kk,mc) =  fr(1,j,k,mc)
            fr(2,jj,kk,mc) = -fr(2,j,k,mc)
          end do
        end do
        do k = 2, 3
          kk = 7 - k
          fr(1,1,kk,mc) =  fr(1,1,k,mc)
          fr(2,1,kk,mc) = -fr(2,1,k,mc)
        end do
        fr(1:2,1,1,mc)  = 0.0
      end do

!     Euler step for Langevin type process with the...
!     noise forcing. Then obtain incompressible part...
!     and save as (a1r,a2r,a3r). Possible wave nos. are...
!     k1=0,1,2 and k2,k3 = 0,1,2,-2,-1 .

      DO k = 1, 5
        kk = k
        if ( k.gt.3 ) kk = k + nz - 5
        DO j = 1, 5
          jj = j
          if ( j.gt.3 ) jj = j + ny - 5
          DO i = 1, 6
            radm = sqrt ( kxr(i)**2 + kyr(jj)**2 + kzr(kk)**2 )
            IF ( radm.lt.xkf ) THEN
              b1r(i,j,k) = b1r(i,j,k) * (1.0-dt/tf) + fr(i,j,k,1)
              b2r(i,j,k) = b2r(i,j,k) * (1.0-dt/tf) + fr(i,j,k,2)
              b3r(i,j,k) = b3r(i,j,k) * (1.0-dt/tf) + fr(i,j,k,3)

              t1 = kxr(i)  * b1r(i,j,k)
              t2 = kyr(jj) * b2r(i,j,k)
              t3 = kzr(kk) * b3r(i,j,k)

              divb = (t1+t2+t3) / ( (radm+smallx)**2 )

              if ( radm.eq.0. ) divb = 0.0

              a1r(i,j,k) = b1r(i,j,k) - kxr(i)  * divb
              a2r(i,j,k) = b2r(i,j,k) - kyr(jj) * divb
              a3r(i,j,k) = b3r(i,j,k) - kzr(kk) * divb
            ELSE
              a1r(i,j,k) = 0.0
              a2r(i,j,k) = 0.0
              a3r(i,j,k) = 0.0
            ENDIF
          ENDDO
        ENDDO
      ENDDO

!     Eliminate forcing of a mean flow at zero wave no.
      a1r(1:2,1,1) = 0.0
      a2r(1:2,1,1) = 0.0
      a3r(1:2,1,1) = 0.0

!     Ensure that the forcing is real-valued. This step...
!     should only clean up rounding errors if fr is ok.

      DO k = 1, 5
        kk = 7 - k
        if ( k.eq.1 ) kk = 1
        do j = 2, 3
          jj = 7 - j
          a1r(1,jj,kk) =  a1r(1,j,k)
          a1r(2,jj,kk) = -a1r(2,j,k)
          a2r(1,jj,kk) =  a2r(1,j,k)
          a2r(2,jj,kk) = -a2r(2,j,k)
          a3r(1,jj,kk) =  a3r(1,j,k)
          a3r(2,jj,kk) = -a3r(2,j,k)
        enddo
      ENDDO

      do k = 2, 3
        kk = 7 - k
        a1r(1,1,kk) =  a1r(1,1,k)
        a1r(2,1,kk) = -a1r(2,1,k)
        a2r(1,1,kk) =  a2r(1,1,k)
        a2r(2,1,kk) = -a2r(2,1,k)
        a3r(1,1,kk) =  a3r(1,1,k)
        a3r(2,1,kk) = -a3r(2,1,k)
      enddo

! NOW, let's add the forcing term to advance the velocity
! at the substep using Euler scheme which is effective.
! The forcing is added to the spectral domain velocity only
! in a few mesh nodes where the wave number is lower than
! SQRT(8). Those nodes are 3 at the begining and 2 at the end
! of the kz nd ky directions, and 6 at the begining of kx
! direction (considering that the complex pair are saved
! one next to the other.
! THIS SUBROUTINE WILL ASSUME THAT THE DOMAIN DECOMPOSITION
! WAS DONE SUCH AS THE MENTIONED NODES ARE INSIDE ONLY ONE
! PROCESSOR. This means that the number of processors along
! the Y-direction should be: nprocY > nx/6; and  the number
! of processors along Z-direction should be: nprocZ > nx/3.

! NOTE: a1r,a2r & a3r are based on a1r(kx,ky,kz) while velocity
! is based on u(kz,kx,ky)

      vx = 0.0
      vy = 0.0
      vz = 0.0

      IF (indy.eq.0) THEN
        IF (indz.eq.0) THEN
          do k=1,5
            kk = k
            if(k.gt.3) kk = k+nz-5
            do j=1,3
              jj = j
              do i=1,6
                vx(kk,i,jj) = dt*a1r(i,j,k)
                vy(kk,i,jj) = dt*a2r(i,j,k)
                vz(kk,i,jj) = dt*a3r(i,j,k)
              enddo
            enddo
          enddo
        ENDIF

        IF (indz.eq.(nprocZ-1)) THEN
          do k=1,5
            kk = k
            if(k.gt.3) kk = k+nz-5
            do j=4,5
              jj = j+lz-5
              do i=1,6
                vx(kk,i,jj) = dt*a1r(i,j,k)
                vy(kk,i,jj) = dt*a2r(i,j,k)
                vz(kk,i,jj) = dt*a3r(i,j,k)
              enddo
            enddo
          enddo
        ENDIF
      ENDIF

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      call mpifft3DCR(vx)
      call mpifft3DCR(vy)
      call mpifft3DCR(vz)

      force_realx(:,:,:) = vx(1:lx,1:ly,:)
      force_realy(:,:,:) = vy(1:lx,1:ly,:)
      force_realz(:,:,:) = vz(1:lx,1:ly,:)

      RETURN
      END SUBROUTINE FORCING
!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
        Function ranff(idum,iy,iv)
!
!       Minimal random number generator of Park and Miller with
!       Bays-Durham shuffle and added safegards, see Numerical Recipes
!
        Integer idum, IA,IM,IQ,IR,NTAB,NDIV
        Real ranff,AM,EPS,RNMX
        Parameter (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,   &
             NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
        INTEGER j,k,iv(NTAB),iy

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
!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

