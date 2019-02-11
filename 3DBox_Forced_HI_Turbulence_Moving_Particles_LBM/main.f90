! -------------------------------------------------------------------------
!       MPI F90 Code for Simulating 3D Decaying Homogeneous Isotropic 
!       Turbulence (DHIT) with finite-size freely-moving particles 
!       embedded in the cube.
!
!       Using LBM D3Q19 MRT Model.
!
!       MPI version applicable on bluefire.ucar.edu.
!
!       The DHIT code was originally written by Dr. Lian-Ping Wang, June 2000.
!
!       The particle part of the code was implemented by Hui Gao, Jan. 2010.
!       The original DHIT code was re-structured as well, using Dr. Yan Peng's
!       code as a reference.
!
!           x   y   z
!       {   0   0   0   }      rest
!
!       {   +   0   0   }      dir  1
!       {   -   0   0   }      dir  2
!       {   0   +   0   }      dir  3
!       {   0   -   0   }      dir  4
!       {   0   0   +   }      dir  5
!       {   0   0   -   }      dir  6
!
!       {   +   +   0   }      dir  7
!       {   -   +   0   }      dir  8
!       {   +   -   0   }      dir  9
!       {   -   -   0   }      dir  10
!       {   +   0   +   }      dir  11
!       {   -   0   +   }      dir  12
!       {   +   0   -   }      dir  13
!       {   -   0   -   }      dir  14
!       {   0   +   +   }      dir  15
!       {   0   -   +   }      dir  16
!       {   0   +   -   }      dir  17
!       {   0   -   -   }      dir  18
! -------------------------------------------------------------------------
! this code is identical with test23 code, except that new subroutines
! for statistics of particle and fluid rms velocities are added 
      PROGRAM main
      use mpi
      use var_inc
      implicit none
!     real,dimension(4999) :: time_stream_max_array
      integer:: i,j,k
      character(len=100)::fnm,fnm1,fnm2
      fnm = '/ptmp/canderse/debug2D/ilinks.630.dat'
      fnm1 = '/ptmp/canderse/debug2D/ilinks.631.dat'

      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)      


! Create FFTW plans, store in plan1DF, plan1DB, plan2DRC, plan2DCR
      call set_plan1D_RC
      call set_plan1D_CC

      call para

      call allocarray

      call wave

      IF(newrun)THEN   

      if(newinitflow)then  

!        call initrand(iflowseed)
!        call initvel

        ux = 0.0
        uy = 0.0
        uz = 0.0

        call initpop

        call statistc 
    
        call rhoupdat

! pre-relaxation of density field after initial forcing
! Note: during this stage, ONLY update density,  keep velocity unchanged
        istep = 0
        do
          call collision_MRT 

          call streaming

          rhop = rho
          call rhoupdat
          rhoerr = maxval(abs(rho - rhop))        

          call MPI_ALLREDUCE(rhoerr,rhoerrmax,1,MPI_REAL8,MPI_MAX,      &
                             MPI_COMM_WORLD,ierr)
          if(myid == 0 .and. mod(istep,1) == 0)                        &
            write(*,*)istep, rhoerrmax 

          if(rhoerrmax <= rhoepsl .or. istep > 15000)then
            if(myid == 0) write(*,*)'final relaxation => ',            &
                                    istep, rhoerrmax
            exit
          end if

          istep = istep + 1 

        end do 


! save population for input of next run 
        call saveinitflow    
!        goto 101
        call macrovar 

      else     
! readin population from previous saving

        call loadinitflow     
!        goto 101
        call macrovar 
      
      end if     

      ELSE
! load data from same number of processes  
      call loadcntdflow      
      call input_outputf(1)

      if(ipart .and. istpload > irelease)then
        call loadcntdpart    
        call beads_links
      end if

      call macrovar

      END IF

! main loop
      istep0 = 0

      do istep = istep0+1,istep0+nsteps 

        if(myid.eq.0) write(*,*) 'istep= ',istep


! Release partilces only after proper skewness (~ -0.5) has been established
! Initialise particle center position and velocity
      if(ipart .and. istep == irelease)then
        istep00 = 0

        call initpartpos
!        call loadinitpartpos

      time_start = MPI_WTIME()

        call initpartvel

        call initpartomg 

        call beads_links

        istep00 = 1
 
      end if

!        if(istep==2) call writeflowfieldstart

!       if(myid == 0 .and. mod(istep-1,1) == 0)                        &
!         write(*,*) 'istep = ', istep 

        call FORCING

        call collision_MRT

        if(ipart .and. istep >= irelease) call beads_collision

        call streaming

        if(ipart .and. istep >= irelease)then
          call beads_lubforce
          call beads_move
          call beads_redistribute
          call beads_links
          call beads_filling
        end if

        call macrovar

         if(mod(istep-1,ndiag) == 0) call diag 
         if(mod(istep-1,nstat) == 0)  call statistc 
!        if(mod(istep-1,nflowout) == 0) call outputflow 

         if(ipart .and. istep >= irelease .and. mod(istep,npartout) == 0)  call outputpart

! output fiels and profiles from the particle surface
         if(ipart .and. istep >= irelease .and. mod(istep,nmovieout) == 0) then
           call moviedata
           call sijstat03
           go to 101
         end if

         if(mod(istep,nstat) == 0) call rmsstat
         if(mod(istep,nsij) == 0) call sijstat03   
!        if(mod(istep,nsij) == 0) call sijstat 

! stop and save in the middle to meet the 6 hour time limit
        if(time_max > time_bond) exit


!      call MPI_ALLREDUCE(time_stream,time_stream_max,1,   &
!                 MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)

!     time_stream_max_array(istep) = time_stream_max

      end do

! main loop ends here

       time_end = MPI_WTIME()
       time_diff = time_end - time_start

       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call MPI_ALLREDUCE(time_diff,time_max,1,MPI_REAL8,  &
                           MPI_MAX,MPI_COMM_WORLD,ierr)

! save data for continue run
      call savecntdflow      
      call input_outputf(2)
      if(ipart .and. istep > irelease) call savecntdpart    

101   continue

!     if(myid.eq.0)then
!     write(*,*)'time_max = ',time_max
!     write(50,*)' ',time_stream_max_array
!     endif

! Destroy plans
      call destroy_plan1D_RC
      call destroy_plan1D_CC


      call MPI_FINALIZE(ierr)

      END PROGRAM main
