!==================================================================

      SUBROUTINE initpop 
      use var_inc
      implicit none

      real, dimension(lx,ly,lz):: usqr, G
      integer ip

      usqr = ux*ux + uy*uy + uz*uz
      usqr = 1.5*usqr   

! initialise density
      rho = 0.0 

      f(0,:,:,:) = ww0*(rho - usqr) 

      do ip = 1,6
        G = (cix(ip)*ux + ciy(ip)*uy + ciz(ip)*uz)
        f(ip,:,:,:) = ww1*(rho + 3.0*G + 4.5*G*G - usqr) 
      end do

      do ip = 7,npop-1
        G = (cix(ip)*ux + ciy(ip)*uy + ciz(ip)*uz)
        f(ip,:,:,:) = ww2*(rho + 3.0*G + 4.5*G*G - usqr) 
      end do

      END SUBROUTINE initpop
!==================================================================

      SUBROUTINE initrand(inseed)
      use mpi
      use var_inc
      implicit none
 
      integer(kind = 4) irand
      integer i, myseed, inseed  
      integer, dimension(nproc):: iseed 
 
      if(myid == 0)then
        call srand(inseed)
        do i = 1,nproc
          iseed(i) = irand()
        end do 
      end if
 
      call MPI_BCAST(iseed,nproc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      
      myseed = iseed(myid + 1)
      call srand(myseed) 
     
      END SUBROUTINE initrand
!==================================================================

      SUBROUTINE initvel    
      use mpi    
      use var_inc
      implicit none
      
      real, dimension(lx+2,lly,lz) :: tmp
      real ek, e_t, k9, Ek9  
      integer ik    

      call gaussian(vx)
      call gaussian(vy)
      call gaussian(vz)

      tmp = ( kx*vx + ky*vy + kz*vz ) / k2
      vx  = vx - kx*tmp
      vy  = vy - ky*tmp
      vz  = vz - kz*tmp

      call symmetrize(vx)
      call symmetrize(vy)
      call symmetrize(vz)

      where(ik2 > nek)
        vx = 0.0 
        vy = 0.0
        vz = 0.0
      end where

      tmp = vx*vx + vy*vy + vz*vz
      if(indy.eq.0) then
         tmp(:,1,:) = 0.5 * tmp(:,1,:)
         tmp(:,2,:) = 0.5 * tmp(:,2,:)
      endif

      do ik = 1,nek
        k9 = real(ik) / real(kpeak)
        Ek9 = 1.5*u0in**2*k9/real(kpeak)*exp(-k9)

        ek = 0.0
        e_t = 0.0
        ek = sum(tmp, mask=(ik2 == ik))

        call MPI_ALLREDUCE(ek,e_t,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

        e_t = sqrt(Ek9/e_t)*vscale

        where(ik2 == ik)
          vx = vx * e_t
          vy = vy * e_t
          vz = vz * e_t
        end where
      end do

      call mpifft3DCR(vx)
      call mpifft3DCR(vy)
      call mpifft3DCR(vz)

      ux = vx(1:lx,1:ly,:)
      uy = vy(1:lx,1:ly,:)
      uz = vz(1:lx,1:ly,:)

      END SUBROUTINE initvel   
!==================================================================


