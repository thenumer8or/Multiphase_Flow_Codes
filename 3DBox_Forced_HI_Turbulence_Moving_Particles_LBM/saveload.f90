!===========================================================================

      subroutine loadinitpartpos
      use var_inc
      implicit none

      integer iprc1, iprc2, iprc3
      character (len = 150):: fnm

      real, dimension(3,msize) :: yp_db
      real, dimension(3,npart) :: ypglb_db

      iprc1 = myid / 100
      iprc2 = mod(myid,100) / 10
      iprc3 = mod(myid,10)

      fnm = '/glade/scratch/ayala/LBM2Ddd/cntdpart/'           &
            //'initpartpos_409600/initpartpos.'                        &
            //char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)

      open(14, file = trim(fnm), status = 'unknown',                   &
               form = 'unformatted')

      read(14) nps, ipglb
      read(14) yp_db, ypglb_db

      close(14)

      yp = yp_db
      ypglb = ypglb_db

      end subroutine loadinitpartpos
!===========================================================================

      subroutine saveprerelax
      use var_inc
      implicit none

      integer iprc1, iprc2, iprc3
      character (len = 100):: fnm

      iprc1 = myid / 100
      iprc2 = mod(myid,100) / 10
      iprc3 = mod(myid,10)

      fnm = trim(dirinitflow)//'prerelax_01/finit.'                  &
            //char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)

      open(10, file = trim(fnm), status = 'unknown',                   &
               form = 'unformatted')

      write(10) istep
      write(10) f, rho
      write(10) ux, uy, uz

      close(10)

      end subroutine saveprerelax
!===========================================================================

      subroutine loadprerelax
      use var_inc
      implicit none

      integer iprc1, iprc2, iprc3
      character (len = 100):: fnm

      iprc1 = myid / 100
      iprc2 = mod(myid,100) / 10
      iprc3 = mod(myid,10)

      fnm = trim(dirinitflow)//'prerelax_01/finit.'                  &
            //char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)

      open(10, file = trim(fnm), status = 'unknown',                   &
               form = 'unformatted')

      read(10) istep
      read(10) f, rho
      read(10) ux, uy, uz

      close(10)

      end subroutine loadprerelax
!===========================================================================

      subroutine saveinitflow
      use var_inc
      implicit none

      integer iprc1, iprc2, iprc3
      character (len = 100):: fnm

      iprc1 = myid / 100
      iprc2 = mod(myid,100) / 10
      iprc3 = mod(myid,10) 

      fnm = trim(dirinitflow)//'finit.'                                &
            //char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)

      open(10, file = trim(fnm), status = 'unknown',                   &
               form = 'unformatted')

      write(10) istat 
      write(10) f

      close(10)

      end subroutine saveinitflow
!===========================================================================

      subroutine loadinitflow
      use var_inc
      implicit none

      integer iprc1, iprc2, iprc3
      character (len = 100):: fnm

      real, allocatable, dimension(:,:,:,:):: f9
      allocate (f9(0:npop-1,lx,ly,lz))

      iprc1 = myid / 100
      iprc2 = mod(myid,100) / 10
      iprc3 = mod(myid,10)

! read in double precision f9 in binary format
      fnm = trim(dirinitflow)//'finit.'                                &
            //char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)

!      fnm = '/ptmp//finit.'   &
!            //char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)

      open(10, file = trim(fnm), status = 'unknown',                   &
               form = 'unformatted')

      read(10) istat 
      read(10) f9

      close(10)
      
      f = f9   

! output double precision f9 to ascii format
!      fnm = trim(dirinitflow)//'real4_02/finit.'                       &
!          //char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)

!      open(10, file = trim(fnm), status = 'unknown',                   &
!               form = 'formatted')
!      write(10,105) f9 

!      close(10)

      deallocate (f9)

105   format(2x,8(1pe16.6))

!      istat = 1  
! read in single precision f in ascii format
!      fnm = trim(dirinitflow)//'real4_02/finit.'                       &
!          //char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)

!       open(10, file = trim(fnm), status = 'unknown',                   &
!                form = 'formatted')
!       read(10,105) f

!       close(10)

! output single precision f in ascii format for checking
!       fnm = trim(dirinitflow)//'real4_02chk/finit.'                    &
!          //char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)
 
!       open(10, file = trim(fnm), status = 'unknown',                   &
!               form = 'formatted')
!       write(10,105) f

!       close(10) 
     
      end subroutine loadinitflow
!===========================================================================

      subroutine savecntdflow
      use var_inc
      implicit none

      integer iprc1, iprc2, iprc3
      integer istp1, istp2, istp3, istp4, istp5, istp6  
      character (len = 100):: fnm

      iprc1 = myid / 100
      iprc2 = mod(myid,100) / 10
      iprc3 = mod(myid,10) 

      istp1 = istep / 100000
      istp2 = mod(istep,100000) / 10000
      istp3 = mod(istep,10000) / 1000
      istp4 = mod(istep,1000) / 100
      istp5 = mod(istep,100) / 10
      istp6 = mod(istep,10)

      fnm = trim(dircntdflow)//'endrunflow2D8x8S1.'                   &
            //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)     &
            //char(istp4 + 48)//char(istp5 + 48)//char(istp6 + 48)     &
            //'.'                                                      &
            //char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)

      open(12, file = trim(fnm), status = 'unknown',                   &
               form = 'unformatted')

      write(12) nsteps, istat, imovie 
      write(12) f

      close(12)

      end subroutine savecntdflow      
!!===========================================================================

      SUBROUTINE input_outputf(idirec)

      use mpi
      use var_inc
      IMPLICIT NONE

      integer idump,idirec
      integer istp1, istp2, istp3, istp4, istp5, istp6
      character (len = 100):: fnm

      IF(myid.eq.0) THEN
        if(idirec.eq.1) then
          idump = istpload
        elseif(idirec.eq.2) then
          idump = istep
        endif

        istp1 = idump / 100000
        istp2 = mod(idump,100000) / 10000
        istp3 = mod(idump,10000) / 1000
        istp4 = mod(idump,1000) / 100
        istp5 = mod(idump,100) / 10
        istp6 = mod(idump,10)

        fnm = trim(dircntdflow)//'force.'                              &
              //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)   &
              //char(istp4 + 48)//char(istp5 + 48)//char(istp6 + 48)

        if(idirec.eq.1) then
          open(10,file=trim(fnm),status='unknown',form='unformatted')
          read(10)iseedf
          read(10)ivf
          read(10)iyf
          read(10)b1r,b2r,b3r
          close(10)
        elseif(idirec.eq.2) then
          open(10,file=trim(fnm),status='unknown',form='unformatted')
          write(10)iseedf
          write(10)ivf
          write(10)iyf
          write(10)b1r,b2r,b3r
          close(10)
        endif
      ENDIF

      if(idirec.eq.1) then
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST (b1r,150,MPI_REAL8,0,mpi_comm_world,ierr)
        CALL MPI_BCAST (b2r,150,MPI_REAL8,0,mpi_comm_world,ierr)
        CALL MPI_BCAST (b3r,150,MPI_REAL8,0,mpi_comm_world,ierr)
        CALL MPI_BCAST (iseedf,1,MPI_INTEGER,0,mpi_comm_world,ierr)
        CALL MPI_BCAST (ivf,NTAB,MPI_INTEGER,0,mpi_comm_world,ierr)
        CALL MPI_BCAST (iyf,1,MPI_INTEGER,0,mpi_comm_world,ierr)
      endif


      RETURN
      END SUBROUTINE input_outputf

!===========================================================================
      
      subroutine loadcntdflow
      use var_inc
      implicit none

      integer iprc1, iprc2, iprc3
      integer istp1, istp2, istp3, istp4, istp5, istp6
      character (len = 100):: fnm

      iprc1 = myid / 100
      iprc2 = mod(myid,100) / 10
      iprc3 = mod(myid,10) 

      istp1 = istpload / 100000
      istp2 = mod(istpload,100000) / 10000
      istp3 = mod(istpload,10000) / 1000
      istp4 = mod(istpload,1000) / 100
      istp5 = mod(istpload,100) / 10
      istp6 = mod(istpload,10)

      fnm = trim(dircntdflow)//'endrunflow2D8x8S1.'                 &
            //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)     &
            //char(istp4 + 48)//char(istp5 + 48)//char(istp6 + 48)     &
            //'.'                                                      &
            //char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)

      open(12, file = trim(fnm), status = 'unknown',                   &
               form = 'unformatted')

      read(12) istep0, istat, imovie 
      read(12) f

      close(12)

      end subroutine loadcntdflow      
!===========================================================================
! load data from more number of processes
! e.g., (#_proc_current) = (1/2, 1/4, ..., etc.)*(#_proc_previous)
      subroutine loadcntdflow_frmmore
      use mpi
      use var_inc
      implicit none

      integer iprc1, iprc2, iprc3
      integer istp1, istp2, istp3, istp4, istp5, istp6
      character (len = 100):: fnm
      
      integer ii, myiid, lz9   
      real, allocatable, dimension(:,:,:,:):: f9

      lz9 = lz / iprocrate

      allocate (f9(0:npop-1,lx,ly,lz9))

      istp1 = istpload / 100000
      istp2 = mod(istpload,100000) / 10000
      istp3 = mod(istpload,10000) / 1000
      istp4 = mod(istpload,1000) / 100
      istp5 = mod(istpload,100) / 10
      istp6 = mod(istpload,10)

      do ii = 0, iprocrate-1 
        myiid = myid*iprocrate + ii 

        iprc1 = myiid / 100
        iprc2 = mod(myiid,100) / 10
        iprc3 = mod(myiid,10)

        fnm = trim(dircntdflow)//'endrunflow.'                         &
            //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)     &
            //char(istp4 + 48)//char(istp5 + 48)//char(istp6 + 48)     &
            //'.'                                                      &
            //char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)

        open(12, file = trim(fnm), status = 'unknown',                 &
               form = 'unformatted')

        read(12) istep0, istat, imovie
        read(12) f9
        close(12)

        f(:,:,:,1+ii*lz9:(ii+1)*lz9) = f9
      end do

      deallocate (f9)

      end subroutine loadcntdflow_frmmore
!===========================================================================
! load data from less number of processes
! e.g., (#_proc_current) = (2, 4, ..., etc.)*(#_proc_previous)
      subroutine loadcntdflow_frmless
      use mpi
      use var_inc
      implicit none

      integer iprc1, iprc2, iprc3
      integer istp1, istp2, istp3, istp4, istp5, istp6
      character (len = 100):: fnm

      integer ii, myiid, lz9  
      real, allocatable, dimension(:,:,:,:):: f9

      myiid = myid / iprocrate
      ii = mod(myid,iprocrate)   
      lz9 = lz * iprocrate 

      allocate (f9(0:npop-1,lx,ly,lz9))

      istp1 = istpload / 100000
      istp2 = mod(istpload,100000) / 10000
      istp3 = mod(istpload,10000) / 1000
      istp4 = mod(istpload,1000) / 100
      istp5 = mod(istpload,100) / 10
      istp6 = mod(istpload,10)

      iprc1 = myiid / 100
      iprc2 = mod(myiid,100) / 10
      iprc3 = mod(myiid,10)

      fnm = trim(dircntdflow)//'endrunflow.'                           &
            //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)     &
            //char(istp4 + 48)//char(istp5 + 48)//char(istp6 + 48)     &
            //'.'                                                      &
            //char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)

      open(12, file = trim(fnm), status = 'unknown',                   &
               form = 'unformatted')

      read(12) istep0, istat, imovie
      read(12) f9
      close(12)

      f = f9(:,:,:,1+ii*lz:(ii+1)*lz)

      deallocate (f9)

      end subroutine loadcntdflow_frmless   
!===========================================================================

      subroutine savecntdpart00  
      use var_inc
      implicit none

      integer iprc1, iprc2, iprc3
      integer istp1, istp2, istp3, istp4, istp5, istp6   
      character (len = 100):: fnm

      iprc1 = myid / 100
      iprc2 = mod(myid,100) / 10
      iprc3 = mod(myid,10) 

      istp1 = istep / 100000
      istp2 = mod(istep,100000) / 10000
      istp3 = mod(istep,10000) / 1000
      istp4 = mod(istep,1000) / 100
      istp5 = mod(istep,100) / 10
      istp6 = mod(istep,10)

      fnm = trim(dircntdpart)//'endrunpart.'                           &
            //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)     &
            //char(istp4 + 48)//char(istp5 + 48)//char(istp6 + 48)     &
            //'.'                                                      &
            //char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)

      open(14, file = trim(fnm), status = 'unknown',                   &
               form = 'unformatted')

      write(14) nps, ipglb 
      write(14) ibnodes, ibnodes0
      write(14) isnodes, isnodes0
      write(14) ilink, mlink, alink 
      write(14) ypglb, wp, omgp  
      write(14) yp, thetap    
      write(14) forcepp, torqpp 
      write(14) dwdt, domgdt   

      close(14)

      end subroutine savecntdpart00      
!===========================================================================

      subroutine loadcntdpart00  
      use var_inc
      implicit none

      integer iprc1, iprc2, iprc3
      integer istp1, istp2, istp3, istp4, istp5, istp6
      character (len = 100):: fnm

      iprc1 = myid / 100
      iprc2 = mod(myid,100) / 10
      iprc3 = mod(myid,10)

      istp1 = istpload / 100000
      istp2 = mod(istpload,100000) / 10000
      istp3 = mod(istpload,10000) / 1000
      istp4 = mod(istpload,1000) / 100
      istp5 = mod(istpload,100) / 10
      istp6 = mod(istpload,10)

      fnm = trim(dircntdpart)//'endrunpart.'                           &
            //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)     &
            //char(istp4 + 48)//char(istp5 + 48)//char(istp6 + 48)     &
            //'.'                                                      &
            //char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)

      open(14, file = trim(fnm), status = 'unknown',                   &
               form = 'unformatted')

      read(14) nps, ipglb
      read(14) ibnodes, ibnodes0
      read(14) isnodes, isnodes0
      read(14) ilink, mlink, alink
      read(14) ypglb, wp, omgp
      read(14) yp, thetap
      read(14) forcepp, torqpp
      read(14) dwdt, domgdt

      close(14)

      end subroutine loadcntdpart00  
!===========================================================================

      subroutine savecntdpart01
      use var_inc
      implicit none

      integer iprc1, iprc2, iprc3
      integer istp1, istp2, istp3, istp4, istp5, istp6
      character (len = 100):: fnm

      iprc1 = myid / 100
      iprc2 = mod(myid,100) / 10
      iprc3 = mod(myid,10)

      istp1 = istep / 100000
      istp2 = mod(istep,100000) / 10000
      istp3 = mod(istep,10000) / 1000
      istp4 = mod(istep,1000) / 100
      istp5 = mod(istep,100) / 10
      istp6 = mod(istep,10)

      fnm = trim(dircntdpart)//'endrunpart.'                           &
            //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)     &
            //char(istp4 + 48)//char(istp5 + 48)//char(istp6 + 48)     &
            //'.'                                                      &
            //char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)

      open(14, file = trim(fnm), status = 'unknown',                   &
               form = 'unformatted')

      write(14) nps, ipglb
      write(14) ibnodes0, isnodes0
      write(14) ypglb, wp, omgp
      write(14) yp, thetap
      write(14) forcepp, torqpp
      write(14) dwdt, domgdt

      close(14)

      end subroutine savecntdpart01   
!===========================================================================

      subroutine loadcntdpart01  
      use var_inc
      implicit none

      integer iprc1, iprc2, iprc3
      integer istp1, istp2, istp3, istp4, istp5, istp6
      character (len = 100):: fnm

      iprc1 = myid / 100
      iprc2 = mod(myid,100) / 10
      iprc3 = mod(myid,10)

      istp1 = istpload / 100000
      istp2 = mod(istpload,100000) / 10000
      istp3 = mod(istpload,10000) / 1000
      istp4 = mod(istpload,1000) / 100
      istp5 = mod(istpload,100) / 10
      istp6 = mod(istpload,10)

      fnm = trim(dircntdpart)//'endrunpart.'                           &
            //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)     &
            //char(istp4 + 48)//char(istp5 + 48)//char(istp6 + 48)     &
            //'.'                                                      &
            //char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)

      open(14, file = trim(fnm), status = 'unknown',                   &
               form = 'unformatted')

      read(14) nps, ipglb
      read(14) ibnodes0, isnodes0   
      read(14) ypglb, wp, omgp
      read(14) yp, thetap
      read(14) forcepp, torqpp
      read(14) dwdt, domgdt

      close(14)

      end subroutine loadcntdpart01  
!===========================================================================
! this subroutine is the same as "savecntdpart01" above, except that
! the global variables of "ypglb, wp, omgp" are only stored in process 0,
! which will broadcast the global values to the other processes after loading.
      subroutine savecntdpart
      use var_inc
      implicit none

      integer iprc1, iprc2, iprc3
      integer istp1, istp2, istp3, istp4, istp5, istp6
      character (len = 100):: fnm

      iprc1 = myid / 100
      iprc2 = mod(myid,100) / 10
      iprc3 = mod(myid,10)

      istp1 = istep / 100000
      istp2 = mod(istep,100000) / 10000
      istp3 = mod(istep,10000) / 1000
      istp4 = mod(istep,1000) / 100
      istp5 = mod(istep,100) / 10
      istp6 = mod(istep,10)

      fnm = trim(dircntdpart)//'endrunpart2D.'                           &
            //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)     &
            //char(istp4 + 48)//char(istp5 + 48)//char(istp6 + 48)     &
            //'.'                                                      &
            //char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)

      open(14, file = trim(fnm), status = 'unknown',                   &
               form = 'unformatted')

      write(14) nps, ipglb
      write(14) ibnodes0, isnodes0
      write(14) yp, thetap
      write(14) forcepp, torqpp
      write(14) dwdt, domgdt
      if(myid == 0) write(14) ypglb, wp, omgp

      close(14)

      end subroutine savecntdpart
!===========================================================================

      subroutine loadcntdpart  
      use mpi
      use var_inc
      implicit none

      integer iprc1, iprc2, iprc3
      integer istp1, istp2, istp3, istp4, istp5, istp6
      character (len = 100):: fnm
      
      iprc1 = myid / 100
      iprc2 = mod(myid,100) / 10
      iprc3 = mod(myid,10)

      istp1 = istpload / 100000
      istp2 = mod(istpload,100000) / 10000
      istp3 = mod(istpload,10000) / 1000
      istp4 = mod(istpload,1000) / 100
      istp5 = mod(istpload,100) / 10
      istp6 = mod(istpload,10)

      fnm = trim(dircntdpart)//'endrunpart2D.'                           &
            //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)     &
            //char(istp4 + 48)//char(istp5 + 48)//char(istp6 + 48)     &
            //'.'                                                      &
            //char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)

      open(14, file = trim(fnm), status = 'unknown',                   &
               form = 'unformatted')

      read(14) nps, ipglb
      read(14) ibnodes0, isnodes0   
      read(14) yp, thetap
      read(14) forcepp, torqpp
      read(14) dwdt, domgdt
      if(myid == 0) read(14) ypglb, wp, omgp

      close(14)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)     
      call MPI_BCAST(ypglb,3*npart,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(wp,3*npart,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(omgp,3*npart,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      
      end subroutine loadcntdpart      
!===========================================================================

      subroutine outputflow   
      use var_inc
      implicit none

      integer, dimension(lx,ly,nz):: ibnodes9

      call ibnodeassmbl(ibnodes9)

      call outputux(ibnodes9)

      call outputuy(ibnodes9)

      call outputuz(ibnodes9)

      call outputpress(ibnodes9)

      call outputvort(ibnodes9)

      end subroutine outputflow   
!===========================================================================

      subroutine ibnodeassmbl(ibnodes9) 
      use mpi
      use var_inc
      implicit none

      integer, dimension(lx,ly,nz):: ibnodes9
      integer, dimension(lx,ly,lz):: ibnodes1

      integer ip, ilen
      integer status(MPI_STATUS_SIZE)

      if(ipart)then

        ilen = lx*ly*lz

        if(myid == 0)then
          ibnodes9(:,:,1:lz) = ibnodes 

          do ip = 1,nproc-1
            call MPI_RECV(ibnodes1,ilen,MPI_INTEGER,ip,1,              &
                          MPI_COMM_WORLD,status,ierr)

            ibnodes9(:,:,(ip*lz + 1):((ip + 1)*lz)) = ibnodes1 
          end do

        else
          call MPI_SEND(ibnodes,ilen,MPI_INTEGER,0,1,                  &
                        MPI_COMM_WORLD,ierr)
        end if

      else

        ibnodes9 = -1
      
      end if

      end subroutine ibnodeassmbl      
!===========================================================================


      subroutine outputux(ibnodes9)   
      use mpi 
      use var_inc
      implicit none

      integer, dimension(lx,ly,nz):: ibnodes9

      integer ip, ilen
      integer status(MPI_STATUS_SIZE)
      integer istp1, istp2, istp3, istp4, istp5, istp6
      
      real, dimension(lx,ly,nz):: ux9
      real, dimension(lx,ly,lz):: ux0

      character (len = 100):: fnm

      ilen = lx*ly*lz

      if(myid == 0)then

        ux9(:,:,1:lz) = ux

        do ip = 1,nproc-1
          call MPI_RECV(ux0,ilen,MPI_REAL8,ip,1,MPI_COMM_WORLD,status,ierr)

          ux9(:,:,(ip*lz + 1):((ip + 1)*lz)) = ux0          
        end do

      else
        call MPI_SEND(ux,ilen,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
      end if

      if(myid == 0)then

! zero out the velocity inside particle for plot purpose
!        where(ibnodes9 > 0) ux9 = 0.0

        istp1 = istep / 100000
        istp2 = mod(istep,100000) / 10000
        istp3 = mod(istep,10000) / 1000
        istp4 = mod(istep,1000) / 100
        istp5 = mod(istep,100) / 10
        istp6 = mod(istep,10)    

        fnm = trim(dirflowout)//'ux'                                   &
            //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)     &
            //char(istp4 + 48)//char(istp5 + 48)//char(istp6 + 48)     &
            //'.dat' 

        open(16, file = trim(fnm), status = 'unknown',                 &
                 form = 'formatted')

        write(16,160) ux9 

        close(16)

      end if

160   format(2x,8(1pe16.6))

      end subroutine outputux   
!===========================================================================

      subroutine outputuy(ibnodes9) 
      use mpi 
      use var_inc
      implicit none

      integer, dimension(lx,ly,nz):: ibnodes9

      integer ip, ilen
      integer status(MPI_STATUS_SIZE)
      integer istp1, istp2, istp3, istp4, istp5, istp6
      
      real, dimension(lx,ly,nz):: uy9
      real, dimension(lx,ly,lz):: uy0

      character (len = 100):: fnm

      ilen = lx*ly*lz

      if(myid == 0)then

        uy9(:,:,1:lz) = uy   

        do ip = 1,nproc-1
          call MPI_RECV(uy0,ilen,MPI_REAL8,ip,1,MPI_COMM_WORLD,status,ierr)

          uy9(:,:,(ip*lz + 1):((ip + 1)*lz)) = uy0          
        end do

      else
        call MPI_SEND(uy,ilen,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
      end if

      if(myid == 0)then

! zero out the velocity inside particle for plot purpose
!        where(ibnodes9 > 0) uy9 = 0.0

        istp1 = istep / 100000
        istp2 = mod(istep,100000) / 10000
        istp3 = mod(istep,10000) / 1000
        istp4 = mod(istep,1000) / 100
        istp5 = mod(istep,100) / 10
        istp6 = mod(istep,10)    

        fnm = trim(dirflowout)//'uy'                                   &
            //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)     &
            //char(istp4 + 48)//char(istp5 + 48)//char(istp6 + 48)     &
            //'.dat' 

        open(17, file = trim(fnm), status = 'unknown',                 &
                 form = 'formatted')

        write(17,170) uy9 

        close(17)

      end if

170   format(2x,8(1pe16.6))

      end subroutine outputuy      
!===========================================================================

      subroutine outputuz(ibnodes9)   
      use mpi 
      use var_inc
      implicit none

      integer, dimension(lx,ly,nz):: ibnodes9

      integer ip, ilen
      integer status(MPI_STATUS_SIZE)
      integer istp1, istp2, istp3, istp4, istp5, istp6
      
      real, dimension(lx,ly,nz):: uz9    
      real, dimension(lx,ly,lz):: uz0    

      character (len = 100):: fnm

      ilen = lx*ly*lz

      if(myid == 0)then

        uz9(:,:,1:lz) = uz       

        do ip = 1,nproc-1
          call MPI_RECV(uz0,ilen,MPI_REAL8,ip,1,MPI_COMM_WORLD,status,ierr)

          uz9(:,:,(ip*lz + 1):((ip + 1)*lz)) = uz0          
        end do

      else
        call MPI_SEND(uz,ilen,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
      end if

      if(myid == 0)then

! zero out the velocity inside particle for plot purpose
!        where(ibnodes9 > 0) uz9 = 0.0

        istp1 = istep / 100000
        istp2 = mod(istep,100000) / 10000
        istp3 = mod(istep,10000) / 1000
        istp4 = mod(istep,1000) / 100
        istp5 = mod(istep,100) / 10
        istp6 = mod(istep,10)    

        fnm = trim(dirflowout)//'uz'                                   &
            //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)     &
            //char(istp4 + 48)//char(istp5 + 48)//char(istp6 + 48)     &
            //'.dat' 

        open(18, file = trim(fnm), status = 'unknown',                 &
                 form = 'formatted')

        write(18,180) uz9 

        close(18)

      end if

180   format(2x,8(1pe16.6))

      end subroutine outputuz      
!===========================================================================

      subroutine outputpress(ibnodes9)      
      use mpi 
      use var_inc
      implicit none

      integer, dimension(lx,ly,nz):: ibnodes9

      integer ip, ilen
      integer status(MPI_STATUS_SIZE)
      integer istp1, istp2, istp3, istp4, istp5, istp6
      
      real, dimension(lx,ly,nz):: rho9    
      real, dimension(lx,ly,lz):: rho1    

      character (len = 100):: fnm

      ilen = lx*ly*lz

      if(myid == 0)then

        rho9(:,:,1:lz) = rho       

        do ip = 1,nproc-1
          call MPI_RECV(rho1,ilen,MPI_REAL8,ip,1,MPI_COMM_WORLD,status,ierr)

          rho9(:,:,(ip*lz + 1):((ip + 1)*lz)) = rho1            
        end do

      else
        call MPI_SEND(rho,ilen,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
      end if

      if(myid == 0)then

! zero out the density inside particle for plot purpose
!        where(ibnodes9 > 0) rho9 = 0.0

        istp1 = istep / 100000
        istp2 = mod(istep,100000) / 10000
        istp3 = mod(istep,10000) / 1000
        istp4 = mod(istep,1000) / 100
        istp5 = mod(istep,100) / 10
        istp6 = mod(istep,10)    

        fnm = trim(dirflowout)//'press'                                &
            //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)     &
            //char(istp4 + 48)//char(istp5 + 48)//char(istp6 + 48)     &
            //'.dat' 

        open(19, file = trim(fnm), status = 'unknown',                 &
                 form = 'formatted')

        write(19,190) rho9/3.0 

        close(19)

      end if

190   format(2x,8(1pe16.6))

      end subroutine outputpress      
!===========================================================================

      subroutine outputvort(ibnodes9)      
      use mpi 
      use var_inc
      implicit none

      integer, dimension(lx,ly,nz):: ibnodes9

      integer ip, ilen
      integer status(MPI_STATUS_SIZE)
      integer istp1, istp2, istp3, istp4, istp5, istp6
      
      real, dimension(lx,ly,nz):: vort9    
      real, dimension(lx,ly,lz):: vort, vort0    

      character (len = 100):: fnm

! prepare vorticity field ox, oy, oz, and magnitude field vort
      call vortcalc

      vort = sqrt(ox*ox + oy*oy + oz*oz)
   
      ilen = lx*ly*lz

      if(myid == 0)then

        vort9(:,:,1:lz) = vort         

        do ip = 1,nproc-1
          call MPI_RECV(vort0,ilen,MPI_REAL8,ip,1,MPI_COMM_WORLD,status,ierr)

          vort9(:,:,(ip*lz + 1):((ip + 1)*lz)) = vort0            
        end do

      else
        call MPI_SEND(vort,ilen,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
      end if

      if(myid == 0)then

! zero out the vorticity inside particle for plot purpose
!        where(ibnodes9 > 0) vort9 = 0.0

        istp1 = istep / 100000
        istp2 = mod(istep,100000) / 10000
        istp3 = mod(istep,10000) / 1000
        istp4 = mod(istep,1000) / 100
        istp5 = mod(istep,100) / 10
        istp6 = mod(istep,10)    

        fnm = trim(dirflowout)//'vort'                                 &
            //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)     &
            //char(istp4 + 48)//char(istp5 + 48)//char(istp6 + 48)     &
            //'.dat' 

!       open(20, file = trim(fnm), status = 'unknown',                 &
!                form = 'formatted')

!       write(20,200) vort9      

!       close(20)

      end if

200   format(2x,8(1pe16.6))

      end subroutine outputvort      
!===========================================================================

      subroutine outputpart       
      use mpi
      use var_inc
      implicit none

      integer ip, id, ilen 
      real, dimension(3,npart):: thglb0, thglb

      integer istp1, istp2, istp3, istp4, istp5, istp6   
      character (len = 100):: fnm

! prepare for thglb, the global thetap
      thglb0 = 0.0

      do ip = 1,nps
        id = ipglb(ip)

        thglb0(:,id) = thetap(:,ip)
      end do 

      ilen = 3*npart
      call MPI_ALLREDUCE(thglb0,thglb,ilen,MPI_REAL8,MPI_SUM,           &
                         MPI_COMM_WORLD,ierr)

      if(myid == 0)then

        istp1 = istep / 100000
        istp2 = mod(istep,100000) / 10000
        istp3 = mod(istep,10000) / 1000
        istp4 = mod(istep,1000) / 100
        istp5 = mod(istep,100) / 10
        istp6 = mod(istep,10)    

        fnm = trim(dirpartout)//'partS1'                                 &
            //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)     &
            //char(istp4 + 48)//char(istp5 + 48)//char(istp6 + 48)     &
            //'.dat' 

        open(22, file = trim(fnm), status = 'unknown',                 &
                 form = 'formatted')

        do id = 1,npart
          write(22,220) id, ypglb(1,id), ypglb(2,id), ypglb(3,id),     &
                            thglb(1,id), thglb(2,id), thglb(3,id),     &
                            wp(1,id), wp(2,id), wp(3,id),              &
                            omgp(1,id), omgp(2,id), omgp(3,id),        &
                            fHIp(1,id), fHIp(2,id), fHIp(3,id),        &
                            flubp(1,id), flubp(2,id), flubp(3,id),     &
                            torqp(1,id), torqp(2,id), torqp(3,id)
        end do

        close(22)
      
      end if

220   format(2x,i5,21(1pe16.6))

      end subroutine outputpart    
!===========================================================================
! to compute kinetic energy spectrum, dissipation rate spectrum
! also the skewness and flatness

      subroutine statistc
      use mpi
      use var_inc
      implicit none

      integer ik, nxyz,k,kglb,idZ
      character (len = 100):: fnm1, fnm2, fnm3, fnm4
      real, dimension(lx+2,lly,lz) :: tmph,tmp1h
      real ek, e_t, eks, e_ts, dissp, qq
      real ttt, qqq, dissppp, eta, vk, tk, uprm
      real tmse, Re, xl, et, kmxeta, xintls, vskew, vflat
      real cc2, cc2t, cc3, cc3t, cc4, cc4t

      REAL,DIMENSION     (lx,ly)  ::  tmp2D
      REAL,DIMENSION  (nz) :: vxave,vyave,vzave,vxsq,vysq,vzsq,stress_xz,stress_xy,stress_yz
      REAL,DIMENSION  (nz) :: vxavet,vyavet,vzavet,vxsqt,vysqt,vzsqt,stress_xzt,stress_xyt,stress_yzt

      vx = 0.0
      vy = 0.0
      vz = 0.0

      vx(1:lx,1:ly,:) = ux
      vy(1:lx,1:ly,:) = uy
      vz(1:lx,1:ly,:) = uz

      call mpifft3DRC(vx)
      call mpifft3DRC(vy)
      call mpifft3DRC(vz)


      if(myid == 0)then
        qq = 0.0
        dissp = 0.0
        xintls = 0.0

        fnm1 = trim(dirstat)//'spectrum.dat'
        fnm2 = trim(dirstat)//'monitora.dat'
        fnm3 = trim(dirstat)//'monitorb.dat'
        fnm4 = trim(dirstat)//'profiles.dat'

        open(24, file = trim(fnm1), status = 'unknown',                &
                 form = 'formatted', position = 'append')
        open(25, file = trim(fnm2), status = 'unknown',                &
                 form = 'formatted', position = 'append')
        open(26, file = trim(fnm3), status = 'unknown',                &
                 form = 'formatted', position = 'append')
        open(27, file = trim(fnm4), status = 'unknown',                &
                 form = 'formatted', position = 'append')

        if(mod(istep-1,nspec) == 0) write(24,*) istat
      end if

! skewness and flatness
      wx = kx * vx
      wy = ky * vy
      wz = kz * vz
      tmph = wx
      wx(1:lx , 1:lly-1:2 , :) = -tmph(1:lx , 2:lly:2   , :)
      wx(1:lx , 2:lly:2   , :) =  tmph(1:lx , 1:lly-1:2 , :)
      tmph = wy
      wy(1:lx , 1:lly-1:2 , :) = -tmph(1:lx , 2:lly:2   , :)
      wy(1:lx , 2:lly:2   , :) =  tmph(1:lx , 1:lly-1:2 , :)
      tmph = wz
      wz(1:lx , 1:lly-1:2 , :) = -tmph(1:lx , 2:lly:2   , :)
      wz(1:lx , 2:lly:2   , :) =  tmph(1:lx , 1:lly-1:2 , :)

      call mpifft3DCR(wx)
      call mpifft3DCR(wy)
      call mpifft3DCR(wz)

      tmph = wx**2 + wy**2 + wx**2
      tmph = tmph / 3.0

      cc2 = sum(tmph(1:lx,1:ly,:))
      call MPI_REDUCE(cc2,cc2t,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      tmph = wx**3 + wy**3 + wx**3
      tmph = tmph / 3.0

      cc3 = sum(tmph(1:lx,1:ly,:))
      call MPI_REDUCE(cc3,cc3t,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      tmph = wx**4 + wy**4 + wx**4
      tmph = tmph / 3.0

      cc4 = sum(tmph(1:lx,1:ly,:))
      call MPI_REDUCE(cc4,cc4t,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      if(myid == 0)then
        nxyz = nx*ny*nz
        vskew = cc3t / real(nxyz) / (cc2t/real(nxyz))**1.5
        vflat = cc4t / real(nxyz) / (cc2t/real(nxyz))**2
      end if

! Profiles
       IF ( mod(istep,nspec).eq.0 ) THEN
        vxave = 0.0
        vyave = 0.0
        vzave = 0.0
        vxsq = 0.0
        vysq = 0.0
        vzsq = 0.0
        stress_xz = 0.0
        stress_xy = 0.0
        stress_yz = 0.0

       idZ = int(myid/nprocY)

       do k=1,lz
       kglb = lz*idZ+k
       tmp2D = ux(:,:,k)
       vxave(kglb) = sum (tmp2D(:,:) )
       tmp2D = uy(:,:,k)
       vyave(kglb) = sum (tmp2D(:,:) )
       tmp2D = uz(:,:,k)
       vzave(kglb) = sum (tmp2D(:,:) )

       tmp2D = ux(:,:,k)*uz(:,:,k)
       stress_xz(kglb) = sum ( tmp2D(:,:) )
       tmp2D = ux(:,:,k)*uy(:,:,k)
       stress_xy(kglb) = sum ( tmp2D(:,:) )
       tmp2D = uy(:,:,k)*uz(:,:,k)
       stress_yz(kglb) = sum ( tmp2D(:,:) )

       tmp2D = (ux(:,:,k))**2
       vxsq(kglb) = sum ( tmp2D(:,:) )
       tmp2D = (uy(:,:,k))**2
       vysq(kglb) = sum ( tmp2D(:,:) )
       tmp2D = (uz(:,:,k))**2
       vzsq(kglb) = sum ( tmp2D(:,:) )

       end do

       call MPI_BARRIER(MPI_COMM_WORLD,ierr)

       CALL MPI_ALLREDUCE (vxave,vxavet,nz,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_ALLREDUCE (vyave,vyavet,nz,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_ALLREDUCE (vzave,vzavet,nz,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_ALLREDUCE (vxsq,vxsqt,nz,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_ALLREDUCE (vysq,vysqt,nz,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_ALLREDUCE (vzsq,vzsqt,nz,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_ALLREDUCE (stress_xz,stress_xzt,nz,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_ALLREDUCE (stress_xy,stress_xyt,nz,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_ALLREDUCE (stress_yz,stress_yzt,nz,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)

       if (myid.eq.0) then
       vxavet = vxavet / (nx*ny) / vscale
       vyavet = vyavet / (nx*ny) /vscale
       vzavet = vzavet / (nx*ny) /vscale
       vxsqt = vxsqt / (nx*ny) /vscale**2
       vysqt = vysqt / (nx*ny) /vscale**2
       vzsqt = vzsqt / (nx*ny) /vscale**2
       stress_xzt = stress_xzt / (nx*ny) /vscale**2
       stress_xyt = stress_xyt / (nx*ny) /vscale**2
       stress_yzt = stress_yzt / (nx*ny) /vscale**2

       do k=1,nz
       write(27,460)k-1,vxavet(k),vyavet(k),vzavet(k), &
             vxsqt(k),vysqt(k),vzsqt(k), stress_xzt(k),stress_xyt(k),stress_yzt(k)
       end do
460    format(2x,i5,9(1pe15.6))
        close(27)
       end if

       ENDIF

! kinetic energy and dissipation rate
     tmph = vx*vx + vy*vy + vz*vz
      if(indy.eq.0) then
        tmph(:,1,:) = 0.5 * tmph(:,1,:) ! note: this is already1/2*(u')^2 = tke
        tmph(:,2,:) = 0.5 * tmph(:,2,:)
      endif
      tmp1h = 2. * visc * tmph * k2

      do ik = 1,nek
        ek = 0.0
        e_t = 0.0
        ek = sum(tmph(1:lx,:,:), mask=(ik2(1:lx,:,:) == ik))

        call MPI_REDUCE(ek,e_t,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

        eks = 0.0
        e_ts = 0.0
        eks = sum(tmp1h(1:lx,:,:), mask=(ik2(1:lx,:,:) == ik))

        call MPI_REDUCE(eks,e_ts,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

        if(myid == 0 )then
          qq = qq + e_t
          dissp = dissp + e_ts
          xintls = xintls + e_t / real(ik)
          if( mod(istep-1,nspec)==0) write(24,240) ik, e_t*escale, e_ts*dscale
        end if
      end do

      if(myid == 0)then
!       ttt = real(istat*nstat)*tscale
        ttt = real(istep)*tscale
        qqq = qq*escale
        xintls = xintls*escale
        dissppp = dissp*dscale

        eta    = (anu**3 / dissppp)**0.25 ! Kolmogorov lengh scale
        vk     = (anu * dissppp)**0.25 ! Kolmogorov velocity scale
        tk     = (anu / dissppp)**0.5 ! Kolmogorov time scale
        uprm   = sqrt(2.0 * qqq / 3.0) ! u_prime = u_rms
        tmse   = sqrt(15.0 * anu * uprm**2 / dissppp) ! Taylor microscale
        Re     = uprm * tmse / anu ! Taylor microscale Reynolds number
        xl     = uprm**3 / dissppp ! large eddy lengthscale
        et     = xl / uprm ! large eddy turn-over time
        kmxeta = (real(lxh) - 1.5) * eta ! kmax * Kolmogorov lengthscale
        xintls = xintls * pi / 2.0 / uprm**2 !longitudinal integral lengthscale

        write(25,250) ttt, eta, tmse, xintls, vk, tk, et, xl
        write(26,250) ttt, uprm, qqq, dissppp,      &
                      Re, kmxeta, vskew, vflat
        write(*,*)'istep,uprm,dissppp,kmxeta,vskew, vflat=',istep,uprm,dissppp,kmxeta,vskew, vflat
        close(24)
        close(25)
        close(26)
      end if

      istat = istat + 1

240   format(2x,i5,2(1pe16.6))
250   format(2x,15(1pe16.6))

      end subroutine statistc 
!===========================================================================
! this is to monitor the maximum flow velocity and particle velocity

      subroutine diag 
      use mpi
      use var_inc
      implicit none

!********THIS IS CHANGED*******************************
!      integer, dimension(2):: idwp, idomgp  

      integer, dimension(1):: idwp, idomgp

      real ttt, vmax0, vmax, wpmax, omgpmax
      real, dimension(lx,ly,lz):: vel 
      real, dimension(npart):: wpmag, omgpmag   
     
      character (len = 100):: fnm   

      vel = sqrt(ux*ux + uy*uy + uz*uz)
      where(ibnodes > 0) vel = 0.0

      vmax0 = maxval(vel) 
      call MPI_REDUCE(vmax0,vmax,1,MPI_REAL8,MPI_MAX,0,MPI_COMM_WORLD,ierr)

      if(myid == 0)then
        wpmax = 0.0
        omgpmax = 0.0
        idwp = 0
        idomgp = 0 

        if(ipart)then
          wpmag = sqrt(wp(1,:)**2 + wp(2,:)**2 + wp(3,:)**2)
          wpmax = maxval(wpmag)
          idwp = maxloc(wpmag)
 
          omgpmag = sqrt(omgp(1,:)**2 + omgp(2,:)**2 + omgp(3,:)**2)
          omgpmax = maxval(omgpmag)
          idomgp = maxloc(omgpmag)
        end if

      end if
 
      if(myid == 0)then

        fnm = trim(dirdiag)//'diag.dat'

        open(26, file = trim(fnm), status = 'unknown',                 &
                 form = 'formatted', position = 'append')

        ttt = real(istep)*tscale
        write(26,260) ttt, vmax,                                       &
                      maxval(idwp), wpmax, maxval(idomgp), omgpmax    
!       write(*,260) ttt, vmax,                                       &
!                     maxval(idwp), wpmax, maxval(idomgp), omgpmax
        close(26)

      end if

260   format(2x,2(1pe16.6), i8, 1pe16.6, i8, 1pe16.6) 
 
      end subroutine diag    
!===========================================================================
!===========================================================================
! this is to output vorticity magnitude and center location of particles,
! and radius projected on the xy-plane of z=256.5 (512^3 domain) in the
! Tecplot format. The data will be used for movie production

      subroutine moviedata
      use mpi
      use var_inc
      implicit none

      integer,parameter:: npart9 = npart/10
! eddy turnover time at t = 0
!      real,parameter:: et0 = 2.356480E+01

      integer id, ii, ilen, mpart, i, j, iroot, sly
      integer,dimension(npart9):: idpart

      real xx, yy, zz9, zdist, ttt
      real, dimension(lx,ly,lz):: vort
      real, dimension(nx,ny):: vort9,vort9c
      real, dimension(npart9):: rad9

      character (len = 100):: fnm1, fnm2

      zz9 = 256.5

      iroot = nprocZ/2
      ilen = lx*ly

! prepare vorticity field ox, oy, oz, and vort the magnitude
      call vortcalc

      vort = sqrt(ox*ox + oy*oy + oz*oz)/tscale*pi2/real(ny)
!Set vortivity to zero inside particles
      where(ibnodes > 0) vort = 0.0

      vort9 = 0.0d0

      if(indz == iroot) then
      sly = indy*ly
      vort9(1:lx,sly+1:sly+ly)=vort(1:lx,1:ly,1)
      end if

! Merge the field
! update ypglb
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(vort9,vort9c,nx*ny,MPI_REAL8,MPI_SUM,           &
                         MPI_COMM_WORLD,ierr)

      if(myid == 0)then

        ttt = real(istep)*tscale
! normalized by eddy turnover time at t = 0
!       ttt = ttt / et0

        mpart = 0
        do id = 1,npart
          zdist = abs(ypglb(3,id) - zz9)
          if(zdist < rad)then
            mpart = mpart + 1

            if(mpart > npart9)then
              write(*,*) 'Too many particles on z = 256.5 plane.'
              stop
            end if

            idpart(mpart) = id
            rad9(mpart) = sqrt(rad*rad - zdist*zdist)
          end if
        end do

      end if

      if(myid == 0)then

        fnm1 = trim(dirmoviedata)//'vort_z2565.dat'
        fnm2 = trim(dirmoviedata)//'part_z2565.dat'

        open(28, file = trim(fnm1), status = 'unknown',                &
                 form = 'formatted', position = 'append')
        open(29, file = trim(fnm2), status = 'unknown',                &
                 form = 'formatted', position = 'append')

!       if(imovie == 0)then
!         write(28,280) ' '
!         write(29,290) ' '
!       end if

!       write(28,281) istep, lx/2, ly/2
!       do j = 1,ny
!       do i = 1,nx
!         xx = real(i) - 0.5
!         yy = real(j) - 0.5
          write(28,282) ((vort9c(i,j),i=1,nx),j=1,ny)
!       end do
!       end do

        write(29,291) istep, mpart
        do ii = 1,mpart
          id = idpart(ii)
          write(29,292) id, ypglb(1,id), ypglb(2,id), rad9(ii)
        end do

        close(28)
        close(29)

      end if

      imovie = imovie + 1

280   format(1x,'title = "',a2,'"',/,1x,                               &
             'variables = "x", "y", "vorticity"')
281   format(/,1x,'zone t = "',i5,'", i = ',i5, ', j = ',i5,           &
             ', f = point')
282   format(2x,8(1pe12.4))

290   format(1x,'title = "',a2,'"',/,1x,                               &
             'variables = "num", "x", "y", "diam", "time"')
291   format(/,1x,'istep = "',i8,'", mpart = ',i5,', f = point')
292   format(2x,i5,3(1pe16.5))

      end subroutine moviedata

!===========================================================================
! this is to claculate the temporal evolution of rms particle velocity 
! with rms fluid velocity, and rms particle angular velocity with rms fluid
! vorticity

      subroutine rmsstat 
      use mpi
      use var_inc
      implicit none
! eddy turnover time at t = 0
!      real,parameter:: et0 = 2.356480E+01

      integer nf0, nf 
      character (len = 100):: fnm
      real ttt, portionf     

      real uxt0, uyt0, uzt0, uxt, uyt, uzt      
      real uxmn, uymn, uzmn     
      real uxrms0, uyrms0, uzrms0, uxrms, uyrms, uzrms      
      real velfrms, velfmn, velprms, velpmn   
      real wpxmn, wpymn, wpzmn, wpxrms, wpyrms, wpzrms    

      real oxt0, oyt0, ozt0, oxt, oyt, ozt  
      real oxmn, oymn, ozmn 
      real oxrms0, oyrms0, ozrms0, oxrms, oyrms, ozrms    
      real omgfrms, omgfmn, omgprms, omgpmn  
      real omgpxmn, omgpymn, omgpzmn 
      real omgpxrms, omgpyrms, omgpzrms  

! first calculate fluid rms velocity
      nf0 = count(ibnodes < 0) 
      call MPI_REDUCE(nf0,nf,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      uxt0 = sum(ux, MASK = (ibnodes < 0))
      call MPI_REDUCE(uxt0,uxt,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      uyt0 = sum(uy, MASK = (ibnodes < 0))
      call MPI_REDUCE(uyt0,uyt,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      uzt0 = sum(uz, MASK = (ibnodes < 0))
      call MPI_REDUCE(uzt0,uzt,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      if(myid == 0)then
        uxmn = uxt / real(nf)        
        uymn = uyt / real(nf)        
        uzmn = uzt / real(nf)        
      end if

      call MPI_BCAST(uxmn,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(uymn,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(uzmn,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

      uxrms0 = sum((ux - uxmn)**2, MASK = (ibnodes < 0))
      call MPI_REDUCE(uxrms0,uxrms,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      uyrms0 = sum((uy - uymn)**2, MASK = (ibnodes < 0))
      call MPI_REDUCE(uyrms0,uyrms,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      uzrms0 = sum((uz - uzmn)**2, MASK = (ibnodes < 0))
      call MPI_REDUCE(uzrms0,uzrms,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      if(myid == 0)then
        uxrms = sqrt(uxrms / real(nf))
        uyrms = sqrt(uyrms / real(nf))
        uzrms = sqrt(uzrms / real(nf))

        velfrms = sqrt(uxrms**2 + uyrms**2 + uzrms**2)
        velfmn = sqrt(uxmn**2 + uymn**2 + uzmn**2) 
      end if


! then calculate particle rms velocity
      if(ipart .and. istep >= irelease)then

      if(myid == 0)then
        wpxmn = sum(wp(1,:)) / real(npart)
        wpymn = sum(wp(2,:)) / real(npart)
        wpzmn = sum(wp(3,:)) / real(npart)

        wpxrms = sqrt(sum((wp(1,:) - wpxmn)**2) / real(npart))
        wpyrms = sqrt(sum((wp(2,:) - wpymn)**2) / real(npart))
        wpzrms = sqrt(sum((wp(3,:) - wpzmn)**2) / real(npart))

        velprms = sqrt(wpxrms**2 + wpyrms**2 + wpzrms**2)
        velpmn = sqrt(wpxmn**2 + wpymn**2 + wpzmn**2)
      end if 

      else

      if(myid == 0)then
        wpxmn = 0.0
        wpymn = 0.0
        wpzmn = 0.0
        wpxrms = 0.0
        wpyrms = 0.0
        wpzrms = 0.0
        velprms = 0.0
        velpmn = 0.0
      end if

      end if 
      

! compute fluid rms vorticity
! prepare vorticity field ox, oy, and oz
      call vortcalc

      oxt0 = sum(ox, MASK = (ibnodes < 0))
      call MPI_REDUCE(oxt0,oxt,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      oyt0 = sum(oy, MASK = (ibnodes < 0))
      call MPI_REDUCE(oyt0,oyt,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      ozt0 = sum(oz, MASK = (ibnodes < 0))
      call MPI_REDUCE(ozt0,ozt,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      if(myid == 0)then
        oxmn = oxt / real(nf)
        oymn = oyt / real(nf)
        ozmn = ozt / real(nf)
      end if

      call MPI_BCAST(oxmn,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(oymn,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ozmn,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

      oxrms0 = sum((ox - oxmn)**2, MASK = (ibnodes < 0))
      call MPI_REDUCE(oxrms0,oxrms,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      oyrms0 = sum((oy - oymn)**2, MASK = (ibnodes < 0))
      call MPI_REDUCE(oyrms0,oyrms,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      ozrms0 = sum((oz - ozmn)**2, MASK = (ibnodes < 0))
      call MPI_REDUCE(ozrms0,ozrms,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      if(myid == 0)then
        oxrms = sqrt(oxrms / real(nf))
        oyrms = sqrt(oyrms / real(nf))
        ozrms = sqrt(ozrms / real(nf))

        omgfrms = sqrt(oxrms**2 + oyrms**2 + ozrms**2)
        omgfmn = sqrt(oxmn**2 + oymn**2 + ozmn**2)
      end if

! compute particle rms angular velocity. Note that vorticity is twice 
! the rotation rate of a small fluid line segment oriented along a 
! principal direction of rate-of-strain tensor

      if(ipart .and. istep >= irelease)then

      if(myid == 0)then
        omgpxmn = sum(omgp(1,:)) / real(npart)
        omgpymn = sum(omgp(2,:)) / real(npart)
        omgpzmn = sum(omgp(3,:)) / real(npart)

        omgpxrms = sqrt(sum((omgp(1,:) - omgpxmn)**2) / real(npart))
        omgpyrms = sqrt(sum((omgp(2,:) - omgpymn)**2) / real(npart))
        omgpzrms = sqrt(sum((omgp(3,:) - omgpzmn)**2) / real(npart))

        omgprms = sqrt(omgpxrms**2 + omgpyrms**2 + omgpzrms**2)
        omgpmn = sqrt(omgpxmn**2 + omgpymn**2 + omgpxmn**2)
      end if

      else

      if(myid == 0)then
        omgpxmn = 0.0
        omgpymn = 0.0
        omgpzmn = 0.0
        omgpxrms = 0.0
        omgpyrms = 0.0
        omgpzrms = 0.0
        omgprms = 0.0
        omgpmn = 0.0
      end if

      end if

      if(myid == 0)then
        ttt = real(istep)*tscale
! normalized by eddy turnover time at t = 0
        ttt = ttt / et0

        portionf = real(nf) / real(nx**3)

        fnm = trim(dirstat)//'rms.dat'
        open(30, file = trim(fnm), status = 'unknown',                 &
                 form = 'formatted', position = 'append')

        write(30,300) ttt, nf, portionf,                               &
                      uxmn, uymn, uzmn, uxrms, uyrms, uzrms,           &
                      velfrms, velfmn,                                 &
                      wpxmn, wpymn, wpzmn, wpxrms,wpyrms, wpzrms,      &
                      velprms, velpmn,                                 &
                      oxmn, oymn, ozmn, oxrms, oyrms, ozrms,           &
                      omgfrms/2.0, omgfmn,                             &
             omgpxmn, omgpymn, omgpzmn, omgpxrms, omgpyrms, omgpzrms,  &
                      omgprms, omgpmn,                                 &
                      omgfrms**2, velfrms**2          
        close(30)

      end if

300   format(2x,1(1pe16.6),i12,35(1pe16.6))
       
      end subroutine rmsstat 
!===========================================================================
! This is to compute local strain rate, Sij, and then to get local dissipation 
! rate epsilon = 2*visc*Sij*Sij, and the local kinetic energy 
! tke = 3/2*urms*urms .

      subroutine sijstat00 
      use mpi
      use var_inc
      implicit none
! eddy turnover time at t = 0
!      real,parameter:: et0 = 2.356480E+01

      integer ix, iy, iz
      real, dimension(0:npop-1) :: f9
      real, dimension(lx,ly,lz) :: Sij2
      real rho9, ux9, uy9, uz9, ux9s, uy9s, uz9s  
      real eqm1, eqm6, eqm8, eqm10, eqm11, eqm12   
      real sum1, sum2, sum6, sum7, sum8, sum9, sum10, sum11
      real evlm1, evlm6, evlm8, evlm10, evlm11, evlm12  
      real neqm1, neqm9, neqm11, neqm13, neqm14, neqm15 
      real Sxx, Syy, Szz, Sxy, Syz, Szx 

      integer, allocatable, dimension(:,:):: icnt0, icnt
      real, allocatable, dimension(:,:):: Sij2p0, Sij2p  
      real, allocatable, dimension(:,:):: uprm2p0, uprm2p  

      integer ip, ilen, i, j, k, is  
      integer irrbnd, irr0, icnttt         
      integer istp1, istp2, istp3, istp4, istp5, istp6  
      real xc, yc, zc, xpnt, ypnt, zpnt, rrbnd
      real xx0, yy0, zz0, rr0, Sij2tt, Sij2mn, rrmn, ttt    

      real, dimension(lx,ly,lz) :: uxprm2, uyprm2, uzprm2 
      real, dimension(lx,ly,lz) :: uprm2

      integer nf0, nf
      real uxt0, uyt0, uzt0, uxt, uyt, uzt
      real uxmn, uymn, uzmn
      real uprm2tt, uprm2mn 
      
      character (len = 100):: fnm    

! to calculate Sij*Sij as a local array of (lx,ly,lz), 
! see Yu H. et al. Computers & Fluids 35, pp. 957-965, 2006, Appendix.
      do iz = 1,lz
      do iy = 1,ly
      do ix = 1,lx
      if(ibnodes(ix,iy,iz) < 0)then

        f9 = f(:,ix,iy,iz)

        rho9 = rho(ix,iy,iz)
        ux9 = ux(ix,iy,iz)
        uy9 = uy(ix,iy,iz)
        uz9 = uz(ix,iy,iz)
        ux9s = ux9*ux9
        uy9s = uy9*uy9
        uz9s = uz9*uz9

        eqm1 = -11.0*rho9 + 19.0*(ux9s + uy9s + uz9s)   
        eqm6 = 2.0*ux9s - uy9s - uz9s
        eqm8 = uy9s - uz9s
        eqm10 = ux9*uy9
        eqm11 = uy9*uz9
        eqm12 = ux9*uz9

        sum1 = f9(1) + f9(2) + f9(3) + f9(4) + f9(5) + f9(6)
        sum2 = f9(7) + f9(8) + f9(9) + f9(10) + f9(11) + f9(12)        &
             + f9(13) + f9(14) + f9(15) + f9(16) + f9(17) + f9(18)
        sum6 = f9(1) + f9(2)
        sum7 = f9(3) + f9(4) + f9(5) + f9(6)
        sum8 = f9(7) + f9(8) + f9(9) + f9(10) + f9(11) + f9(12)        &
             + f9(13) + f9(14)
        sum9 = f9(15) + f9(16) + f9(17) + f9(18)
        sum10 = f9(3) + f9(4) - f9(5) - f9(6)
        sum11 = f9(7) + f9(8) + f9(9) + f9(10) - f9(11) - f9(12)       &
              - f9(13) - f9(14)

        evlm1 = -30.0*f9(0) + coef2*sum1 + coef3*sum2
        evlm6 = coef5*sum6 - sum7 + sum8 - coef5*sum9
        evlm8 = sum10 + sum11
        evlm10 = f9(7) - f9(8) - f9(9) + f9(10)
        evlm11 = f9(15) - f9(16) - f9(17) + f9(18)
        evlm12 = f9(11) - f9(12) - f9(13) + f9(14)

        neqm1 = s1*(evlm1 - eqm1)  
        neqm9 = s9*(evlm6 - eqm6)
        neqm11 = s9*(evlm8 - eqm8)
        neqm13 = s9*(evlm10 - eqm10)
        neqm14 = s9*(evlm11 - eqm11)
        neqm15 = s9*(evlm12 - eqm12)   

        Sxx = -(neqm1 + 19.0*neqm9) / 38.0
        Syy = -(2.0*neqm1 - 19.0*(neqm9 - 3.0*neqm11)) / 76.0
        Szz = -(2.0*neqm1 - 19.0*(neqm9 + 3.0*neqm11)) / 76.0
        Sxy = -1.5*neqm13    
        Syz = -1.5*neqm14    
        Szx = -1.5*neqm15    

        Sij2(ix,iy,iz) = Sxx*Sxx + Syy*Syy + Szz*Szz                   &
                 + 2.0*(Sxy*Sxy + Syz*Syz + Szx*Szx)
      end if
      end do
      end do
      end do

! to calculate urms*urms as a local array of (lx,ly,lz)
! first calculate fluid rms velocity
      nf0 = count(ibnodes < 0)
      call MPI_REDUCE(nf0,nf,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      uxt0 = sum(ux, MASK = (ibnodes < 0))
      call MPI_REDUCE(uxt0,uxt,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      uyt0 = sum(uy, MASK = (ibnodes < 0))
      call MPI_REDUCE(uyt0,uyt,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      uzt0 = sum(uz, MASK = (ibnodes < 0))
      call MPI_REDUCE(uzt0,uzt,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      if(myid == 0)then
        uxmn = uxt / real(nf)
        uymn = uyt / real(nf)
        uzmn = uzt / real(nf)
      end if

      call MPI_BCAST(uxmn,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(uymn,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(uzmn,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

      uxprm2 = (ux - uxmn)**2
      uyprm2 = (uy - uymn)**2
      uzprm2 = (uz - uzmn)**2
      
      uprm2 = (uxprm2 + uyprm2 + uzprm2) / 3.0


! to calculate SijSij and urms*urms as a function of (r/a), need to 
! perform volume-averaging in sperical shells around a particle  
      if(ipart .and. istep >= irelease)then

      rrbnd = real(lx) / (real(npart)**(1.0/3.0)) / 2.0 
      irrbnd = int(rrbnd - rad) + 1 

      allocate (Sij2p0(irrbnd,npart))
      allocate (uprm2p0(irrbnd,npart))
      allocate (icnt0(irrbnd,npart))
      allocate (Sij2p(irrbnd,npart))
      allocate (uprm2p(irrbnd,npart))
      allocate (icnt(irrbnd,npart))

      Sij2p0 = 0.0
      uprm2p0 = 0.0
      icnt0 = 0

      do ip = 1,npart
        xc = ypglb(1,ip)
        yc = ypglb(2,ip)
        zc = ypglb(3,ip)
        
        do k = 1,lz
          zpnt = real(k) - 0.5 + real(myid*lz)

! use the nearest particle center instead of the real center
          if((zc - zpnt) > real(nzh)) zc = zc - real(nz)
          if((zc - zpnt) < -real(nzh)) zc = zc + real(nz)

          zz0 = zpnt - zc
        IF(abs(zz0) <= rrbnd)THEN

        do j = 1,ly
          ypnt = real(j) - 0.5

! use the nearest particle center instead of the real center
          if((yc - ypnt) > real(nyh)) yc = yc - real(ny)
          if((yc - ypnt) < -real(nyh)) yc = yc + real(ny)

          yy0 = ypnt - yc
        IF(abs(yy0) <= rrbnd)THEN

        do i = 1,lx
        IF(ibnodes(i,j,k) < 0)THEN
          xpnt = real(i) - 0.5

! use the nearest particle center instead of the real center
          if((xc - xpnt) > real(lxh)) xc = xc - real(lx)
          if((xc - xpnt) < -real(lxh)) xc = xc + real(lx)
!          if(abs(xc -xpnt) <= real(lxh)) xc = xc

          xx0 = xpnt - xc
        IF(abs(xx0) <= rrbnd)THEN

          rr0 = sqrt(xx0*xx0 + yy0*yy0 + zz0*zz0)  
          if(rr0 <= rrbnd)then

            irr0 = int(rr0 - rad) + 1  
            if(rr0 < rad)then
              write(*,*) 'warning: counting node inside particles'
              stop
            end if

            Sij2p0(irr0,ip) = Sij2p0(irr0,ip) + Sij2(i,j,k)   
            uprm2p0(irr0,ip) = uprm2p0(irr0,ip) + uprm2(i,j,k)   
            icnt0(irr0,ip) = icnt0(irr0,ip) + 1
          end if  
        END IF
        END IF
        end do
        END IF
        end do
        END IF
        end do
      end do

! now collect info. for Sij2p and icnt
      ilen = irrbnd*npart

      call MPI_ALLREDUCE(Sij2p0,Sij2p,ilen,MPI_REAL8,MPI_SUM,           &
                         MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(uprm2p0,uprm2p,ilen,MPI_REAL8,MPI_SUM,         &
                         MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(icnt0,icnt,ilen,MPI_INTEGER,MPI_SUM,          &
                         MPI_COMM_WORLD,ierr)

      deallocate (Sij2p0)
      deallocate (uprm2p0)
      deallocate (icnt0)

! then average the Sij2 for each particle in each spherical shell layer
      if(myid == 0)then
        ttt = real(istep)*tscale
! normalized by eddy turnover time at t = 0
        ttt = ttt / et0

        istp1 = istep / 100000
        istp2 = mod(istep,100000) / 10000
        istp3 = mod(istep,10000) / 1000
        istp4 = mod(istep,1000) / 100
        istp5 = mod(istep,100) / 10
        istp6 = mod(istep,10)

        fnm = trim(dirstat)//'sij'                                     &
            //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)     &
            //char(istp4 + 48)//char(istp5 + 48)//char(istp6 + 48)     &
            //'.dat'

        open(31, file = trim(fnm), status = 'unknown',                 &
                 form = 'formatted', position = 'append')

        do is = 1,irrbnd
          Sij2tt = sum(Sij2p(is,:))
          uprm2tt = sum(uprm2p(is,:))
          icnttt = sum(icnt(is,:))  
          if(icnttt == 0)then
            Sij2mn = 0.0
            uprm2mn = 0.0
          else
            Sij2mn = Sij2tt / real(icnttt)
            uprm2mn = uprm2tt / real(icnttt)
          end if
          rrmn = (real(is) - 0.5) / rad + 1.0  
          write(31,310) ttt, rrmn, Sij2mn, uprm2mn    
        end do  

      end if
310   format(2x,4(1pe16.6))  

      deallocate (Sij2p)
      deallocate (uprm2p)
      deallocate (icnt)

      end if 

      end subroutine sijstat00 
!===========================================================================
! This is to compute local strain rate, Sij, and then to get local dissipation 
! rate epsilon = 2*visc*Sij*Sij, and the local kinetic energy 
! tke = 3/2*urms*urms .
! 2nd edition modification: 
! (1) calculate the Sij using both LBGK and MRT
! (2) spherical shell width size variation: from 1*dx to different drr's setup,
!     note it's better to set multiple different drr's simultaneously
! (3) change the range of profile: from 3*dx to r/a_p = 10 
! (4) Note: each fluid node to be counted only once for the compuation of 
!     SijSij and Uprime*Uprime, i.e., only consider the nearest particle
 
      subroutine sijstat01 
      use mpi
      use var_inc
      implicit none
! eddy turnover time at t = 0
!      real,parameter:: et0 = 2.356480E+01

      integer ix, iy, iz
      real, dimension(0:npop-1) :: f9
      real, dimension(lx,ly,lz) :: Sij2_m1, Sij2_m2
      real rho9, ux9, uy9, uz9, ux9s, uy9s, uz9s  
      real eqm1, eqm6, eqm8, eqm10, eqm11, eqm12   
      real sum1, sum2, sum6, sum7, sum8, sum9, sum10, sum11
      real evlm1, evlm6, evlm8, evlm10, evlm11, evlm12  
      real neqm1, neqm9, neqm11, neqm13, neqm14, neqm15 
      real Sxx, Syy, Szz, Sxy, Syz, Szx 

      real, dimension(lx,ly,lz):: usqr, edtu 
      real, dimension(0:npop-1,lx,ly,lz):: feq, fneq 
      real fneq9 

      integer, allocatable, dimension(:,:):: icnt0, icnt
      real, allocatable, dimension(:,:)::Sij2p0_m1, Sij2p_m1
      real, allocatable, dimension(:,:)::Sij2p0_m2, Sij2p_m2
      real, allocatable, dimension(:,:):: uprm2p0, uprm2p  

      integer ip, ilen, i, j, k, is  
      integer irrbnd, irr0, icnttt, iprrmin         
      integer istp1, istp2, istp3, istp4, istp5, istp6  
      real xc, yc, zc, xpnt, ypnt, zpnt, rrbnd, drr, rrmin  
      real xx0, yy0, zz0, rr0, ttt    
      real Sij2tt_m1, Sij2tt_m2, Sij2mn_m1, Sij2mn_m2, rrmn
      real Sij2_m1t0, Sij2_m2t0, Sij2_m1t, Sij2_m2t  
      real Sij2_m1_avg, Sij2_m2_avg, Sij2m2_tmp  
      real Sij2_m1_norm, Sij2_m2_norm  

      real, dimension(lx,ly,lz) :: uxprm2, uyprm2, uzprm2 
      real, dimension(lx,ly,lz) :: uprm2

      integer nf0, nf
      real uxt0, uyt0, uzt0, uxt, uyt, uzt
      real uxmn, uymn, uzmn
      real uprm2tt, uprm2mn 
      real uprm2t0, uprm2t, uprm2_avg, uprm2_norm  
      
      character (len = 100):: fnm    

! to calculate Sij*Sij as a local array of (lx,ly,lz), 
! Method 1: calculate in the moment space. 
! see Yu H. et al. Computers & Fluids 35, pp. 957-965, 2006, Appendix.
      do iz = 1,lz
      do iy = 1,ly
      do ix = 1,lx
      if(ibnodes(ix,iy,iz) < 0)then

        f9 = f(:,ix,iy,iz)

        rho9 = rho(ix,iy,iz)
        ux9 = ux(ix,iy,iz)
        uy9 = uy(ix,iy,iz)
        uz9 = uz(ix,iy,iz)
        ux9s = ux9*ux9
        uy9s = uy9*uy9
        uz9s = uz9*uz9

        eqm1 = -11.0*rho9 + 19.0*(ux9s + uy9s + uz9s)   
        eqm6 = 2.0*ux9s - uy9s - uz9s
        eqm8 = uy9s - uz9s
        eqm10 = ux9*uy9
        eqm11 = uy9*uz9
        eqm12 = ux9*uz9

        sum1 = f9(1) + f9(2) + f9(3) + f9(4) + f9(5) + f9(6)
        sum2 = f9(7) + f9(8) + f9(9) + f9(10) + f9(11) + f9(12)        &
             + f9(13) + f9(14) + f9(15) + f9(16) + f9(17) + f9(18)
        sum6 = f9(1) + f9(2)
        sum7 = f9(3) + f9(4) + f9(5) + f9(6)
        sum8 = f9(7) + f9(8) + f9(9) + f9(10) + f9(11) + f9(12)        &
             + f9(13) + f9(14)
        sum9 = f9(15) + f9(16) + f9(17) + f9(18)
        sum10 = f9(3) + f9(4) - f9(5) - f9(6)
        sum11 = f9(7) + f9(8) + f9(9) + f9(10) - f9(11) - f9(12)       &
              - f9(13) - f9(14)

        evlm1 = -30.0*f9(0) + coef2*sum1 + coef3*sum2
        evlm6 = coef5*sum6 - sum7 + sum8 - coef5*sum9
        evlm8 = sum10 + sum11
        evlm10 = f9(7) - f9(8) - f9(9) + f9(10)
        evlm11 = f9(15) - f9(16) - f9(17) + f9(18)
        evlm12 = f9(11) - f9(12) - f9(13) + f9(14)

        neqm1 = s1*(evlm1 - eqm1)  
        neqm9 = s9*(evlm6 - eqm6)
        neqm11 = s9*(evlm8 - eqm8)
        neqm13 = s9*(evlm10 - eqm10)
        neqm14 = s9*(evlm11 - eqm11)
        neqm15 = s9*(evlm12 - eqm12)   

        Sxx = -(neqm1 + 19.0*neqm9) / 38.0
        Syy = -(2.0*neqm1 - 19.0*(neqm9 - 3.0*neqm11)) / 76.0
        Szz = -(2.0*neqm1 - 19.0*(neqm9 + 3.0*neqm11)) / 76.0
        Sxy = -1.5*neqm13    
        Syz = -1.5*neqm14    
        Szx = -1.5*neqm15    

        Sij2_m1(ix,iy,iz) = Sxx*Sxx + Syy*Syy + Szz*Szz                &
                    + 2.0*(Sxy*Sxy + Syz*Syz + Szx*Szx)
      end if
      end do
      end do
      end do

! to calculate Sij*Sij as a local array of (lx,ly,lz), 
! Method 2: calculate in the discrete velocity space, f-space. 
! Sij = -3/(2*rho0*tau)*Sigma_alpha{f_alpha(x,t)-f^(0)_alpha(x,t)}*e_alpha_i*
! e_alpha_j, where alpha = 0, 1, 2, ..., 18, f^(0) = f^(eq), the equilibrium 
! distribution function, i,j = x, y, z 

! first, calculate the equilibrium distribution function 
      usqr = ux*ux + uy*uy + uz*uz
      usqr = 1.5*usqr
    
      feq(0,:,:,:) = ww0*(rho - usqr) 

      do ip = 1,6
        edtu = (cix(ip)*ux + ciy(ip)*uy + ciz(ip)*uz)
        feq(ip,:,:,:) = ww1*(rho + 3.0*edtu + 4.5*edtu**2 - usqr)
      end do

      do ip = 7,npop-1
        edtu = (cix(ip)*ux + ciy(ip)*uy + ciz(ip)*uz)
        feq(ip,:,:,:) = ww2*(rho + 3.0*edtu + 4.5*edtu**2 - usqr)
      end do

      fneq = f - feq

      do iz = 1,lz
      do iy = 1,ly
      do ix = 1,lx
      if(ibnodes(ix,iy,iz) < 0)then
        Sxx = 0.0
        Syy = 0.0
        Szz = 0.0
        Sxy = 0.0
        Syz = 0.0
        Szx = 0.0
        do ip = 1,npop-1
          fneq9 = fneq(ip,ix,iy,iz)
          Sxx = Sxx + fneq9*real(cix(ip)*cix(ip))
          Syy = Syy + fneq9*real(ciy(ip)*ciy(ip))
          Szz = Szz + fneq9*real(ciz(ip)*ciz(ip))
          Sxy = Sxy + fneq9*real(cix(ip)*ciy(ip))
          Syz = Syz + fneq9*real(ciy(ip)*ciz(ip))
          Szx = Szx + fneq9*real(ciz(ip)*cix(ip))
        end do

        Sij2m2_tmp = Sxx*Sxx + Syy*Syy + Szz*Szz                &
                    + 2.0*(Sxy*Sxy + Syz*Syz + Szx*Szx)

        Sij2_m2(ix,iy,iz) = Sij2m2_tmp*(1.5 / tau)**2
      end if
      end do
      end do
      end do

! to calculate urms*urms as a local array of (lx,ly,lz)
! first calculate fluid rms velocity
      nf0 = count(ibnodes < 0)
      call MPI_REDUCE(nf0,nf,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      uxt0 = sum(ux, MASK = (ibnodes < 0))
      call MPI_REDUCE(uxt0,uxt,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      uyt0 = sum(uy, MASK = (ibnodes < 0))
      call MPI_REDUCE(uyt0,uyt,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      uzt0 = sum(uz, MASK = (ibnodes < 0))
      call MPI_REDUCE(uzt0,uzt,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      if(myid == 0)then
        uxmn = uxt / real(nf)
        uymn = uyt / real(nf)
        uzmn = uzt / real(nf)
      end if

      call MPI_BCAST(uxmn,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(uymn,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(uzmn,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

      uxprm2 = (ux - uxmn)**2
      uyprm2 = (uy - uymn)**2
      uzprm2 = (uz - uzmn)**2
      
      uprm2 = (uxprm2 + uyprm2 + uzprm2) / 3.0

! to calculate the averaged value of Sij2_m1, Sij2_m2, and uprm2 
! over all the fluid nodes in the domain
      Sij2_m1t0 = sum(Sij2_m1, MASK = (ibnodes < 0))
      call MPI_REDUCE(Sij2_m1t0,Sij2_m1t,1,MPI_REAL8,MPI_SUM,0,         &
                      MPI_COMM_WORLD,ierr)

      Sij2_m2t0 = sum(Sij2_m2, MASK = (ibnodes < 0))
      call MPI_REDUCE(Sij2_m2t0,Sij2_m2t,1,MPI_REAL8,MPI_SUM,0,         &
                      MPI_COMM_WORLD,ierr)

      uprm2t0 = sum(uprm2, MASK = (ibnodes < 0))
      call MPI_REDUCE(uprm2t0,uprm2t,1,MPI_REAL8,MPI_SUM,0,             &
                      MPI_COMM_WORLD,ierr)

      if(myid == 0)then
        Sij2_m1_avg = Sij2_m1t / real(nf)
        Sij2_m2_avg = Sij2_m2t / real(nf)
        uprm2_avg = uprm2t / real(nf)
      end if

! to calculate SijSij and urms*urms as a function of (r/a), need to 
! perform volume-averaging in sperical shells around a particle  
      if(ipart .and. istep >= irelease)then

      rrbnd = 10.0*rad
      drr =0.05*rad
      irrbnd = int((rrbnd - rad) / drr) + 1 

      allocate (Sij2p0_m1(irrbnd,npart))
      allocate (Sij2p0_m2(irrbnd,npart))
      allocate (uprm2p0(irrbnd,npart))
      allocate (icnt0(irrbnd,npart))
      allocate (Sij2p_m1(irrbnd,npart))
      allocate (Sij2p_m2(irrbnd,npart))
      allocate (uprm2p(irrbnd,npart))
      allocate (icnt(irrbnd,npart))

      Sij2p0_m1 = 0.0
      Sij2p0_m2 = 0.0
      uprm2p0 = 0.0
      icnt0 = 0

      do k = 1,lz
      do j = 1,ly
      do i = 1,lx
      if(ibnodes(i,j,k) < 0)then

        xpnt = real(i) - 0.5
        ypnt = real(j) - 0.5
        zpnt = real(k) - 0.5 + real(myid*lz)

        rrmin = real(lx) 
        iprrmin = -1
        do ip = 1,npart
          xc = ypglb(1,ip)
          yc = ypglb(2,ip)
          zc = ypglb(3,ip)     

! use the nearest particle center instead of the real center
          if((xc - xpnt) > real(lxh)) xc = xc - real(lx)
          if((xc - xpnt) < -real(lxh)) xc = xc + real(lx)

          if((yc - ypnt) > real(nyh)) yc = yc - real(ny)
          if((yc - ypnt) < -real(nyh)) yc = yc + real(ny)

          if((zc - zpnt) > real(nzh)) zc = zc - real(nz)
          if((zc - zpnt) < -real(nzh)) zc = zc + real(nz)

          xx0 = xpnt - xc
          yy0 = ypnt - yc
          zz0 = zpnt - zc

          rr0 = sqrt(xx0*xx0 + yy0*yy0 + zz0*zz0)
          if(rr0 <= rrbnd .and. rr0 < rrmin)then
            rrmin = rr0
            iprrmin = ip
          end if
        end do

        if(iprrmin > 0)then
          irr0 = int((rrmin - rad) / drr) + 1          
          Sij2p0_m1(irr0,iprrmin) = Sij2p0_m1(irr0,iprrmin) + Sij2_m1(i,j,k)
          Sij2p0_m2(irr0,iprrmin) = Sij2p0_m2(irr0,iprrmin) + Sij2_m2(i,j,k)
          uprm2p0(irr0,iprrmin) = uprm2p0(irr0,iprrmin) + uprm2(i,j,k)
          icnt0(irr0,iprrmin) = icnt0(irr0,iprrmin) + 1
        end if

      end if
      end do
      end do
      end do

! now collect info. for Sij2p and icnt
      ilen = irrbnd*npart

      call MPI_ALLREDUCE(Sij2p0_m1,Sij2p_m1,ilen,MPI_REAL8,MPI_SUM,     &
                         MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(Sij2p0_m2,Sij2p_m2,ilen,MPI_REAL8,MPI_SUM,     &
                         MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(uprm2p0,uprm2p,ilen,MPI_REAL8,MPI_SUM,         &
                         MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(icnt0,icnt,ilen,MPI_INTEGER,MPI_SUM,          &
                         MPI_COMM_WORLD,ierr)

      deallocate (Sij2p0_m1)
      deallocate (Sij2p0_m2)
      deallocate (uprm2p0)
      deallocate (icnt0)

! then average the Sij2 for each particle in each spherical shell layer
      if(myid == 0)then
        ttt = real(istep)*tscale
! normalized by eddy turnover time at t = 0
        ttt = ttt / et0

        istp1 = istep / 100000
        istp2 = mod(istep,100000) / 10000
        istp3 = mod(istep,10000) / 1000
        istp4 = mod(istep,1000) / 100
        istp5 = mod(istep,100) / 10
        istp6 = mod(istep,10)

        fnm = trim(dirstat)//'sij'                                     &
            //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)     &
            //char(istp4 + 48)//char(istp5 + 48)//char(istp6 + 48)     &
            //'.dat'

        open(31, file = trim(fnm), status = 'unknown',                 &
                 form = 'formatted', position = 'append')

        do is = 1,irrbnd
          Sij2tt_m1 = sum(Sij2p_m1(is,:))
          Sij2tt_m2 = sum(Sij2p_m2(is,:))
          uprm2tt = sum(uprm2p(is,:))
          icnttt = sum(icnt(is,:))  
          if(icnttt == 0)then
            Sij2mn_m1 = 0.0
            Sij2mn_m2 = 0.0
            uprm2mn = 0.0
          else
            Sij2mn_m1 = Sij2tt_m1 / real(icnttt)
            Sij2mn_m2 = Sij2tt_m2 / real(icnttt)
            uprm2mn = uprm2tt / real(icnttt)
          end if
          rrmn = real(is)*drr / rad + 1.0 
          Sij2_m1_norm = Sij2mn_m1 / Sij2_m1_avg 
          Sij2_m2_norm = Sij2mn_m2 / Sij2_m2_avg 
          uprm2_norm = uprm2mn / uprm2_avg 
          write(31,310) ttt, rrmn, Sij2mn_m1, Sij2mn_m2, uprm2mn,      &
                        Sij2_m1_avg, Sij2_m2_avg, uprm2_avg,           &
                        Sij2_m1_norm, Sij2_m2_norm, uprm2_norm   
        end do  

      end if
310   format(2x,11(1pe16.6))  

      deallocate (Sij2p_m1)
      deallocate (Sij2p_m2)
      deallocate (uprm2p)
      deallocate (icnt)

      end if 

      end subroutine sijstat01 
!===========================================================================
! this subroutine "sijstat" is identical with "sijstat01" above, the difference
! is two arrays of dimension(0:npop-1,lx,ly,lz):: feq, fneq are replaced. 
! the goal is to optimize memory usage.

      subroutine sijstat02
      use mpi
      use var_inc
      implicit none

      integer ix, iy, iz
      real, dimension(0:npop-1) :: f9
      real, dimension(lx,ly,lz) :: Sij2_m1, Sij2_m2
      real rho9, ux9, uy9, uz9, ux9s, uy9s, uz9s
      real eqm1, eqm6, eqm8, eqm10, eqm11, eqm12
      real sum1, sum2, sum6, sum7, sum8, sum9, sum10, sum11
      real evlm1, evlm6, evlm8, evlm10, evlm11, evlm12
      real neqm1, neqm9, neqm11, neqm13, neqm14, neqm15
      real Sxx, Syy, Szz, Sxy, Syz, Szx
      real usqr, edtu, wght, feq9, fneq9

      integer, allocatable, dimension(:,:):: icnt0, icnt
      real, allocatable, dimension(:,:)::Sij2p0_m1, Sij2p_m1
      real, allocatable, dimension(:,:)::Sij2p0_m2, Sij2p_m2
      real, allocatable, dimension(:,:):: uprm2p0, uprm2p
      real, allocatable, dimension(:,:,:) :: uprm2

      integer ip, ilen, i, j, k, is
      integer irrbnd, irr0, icnttt, iprrmin
      integer istp1, istp2, istp3, istp4, istp5, istp6
      real xc, yc, zc, xpnt, ypnt, zpnt, rrbnd, drr, rrmin
      real xx0, yy0, zz0, rr0, ttt
      real Sij2tt_m1, Sij2tt_m2, Sij2mn_m1, Sij2mn_m2, rrmn
      real Sij2_m1t0, Sij2_m2t0, Sij2_m1t, Sij2_m2t
      real Sij2_m1_avg, Sij2_m2_avg, Sij2m2_tmp
      real Sij2_m1_norm, Sij2_m2_norm


      integer nf0, nf
      real uxt0, uyt0, uzt0, uxt, uyt, uzt
      real uxmn, uymn, uzmn
      real uprm2tt, uprm2mn
      real uprm2t0, uprm2t, uprm2_avg, uprm2_norm

      character (len = 100):: fnm

! to calculate Sij*Sij as a local array of (lx,ly,lz),
! Method 1: calculate in the moment space.
! see Yu H. et al. Computers & Fluids 35, pp. 957-965, 2006, Appendix.
      do iz = 1,lz
      do iy = 1,ly
      do ix = 1,lx
      if(ibnodes(ix,iy,iz) < 0)then

        f9 = f(:,ix,iy,iz)

        rho9 = rho(ix,iy,iz)
        ux9 = ux(ix,iy,iz)
        uy9 = uy(ix,iy,iz)
        uz9 = uz(ix,iy,iz)
        ux9s = ux9*ux9
        uy9s = uy9*uy9
        uz9s = uz9*uz9

        eqm1 = -11.0*rho9 + 19.0*(ux9s + uy9s + uz9s)
        eqm6 = 2.0*ux9s - uy9s - uz9s
        eqm8 = uy9s - uz9s
        eqm10 = ux9*uy9
        eqm11 = uy9*uz9
        eqm12 = ux9*uz9

        sum1 = f9(1) + f9(2) + f9(3) + f9(4) + f9(5) + f9(6)
        sum2 = f9(7) + f9(8) + f9(9) + f9(10) + f9(11) + f9(12)        &
             + f9(13) + f9(14) + f9(15) + f9(16) + f9(17) + f9(18)
        sum6 = f9(1) + f9(2)
        sum7 = f9(3) + f9(4) + f9(5) + f9(6)
        sum8 = f9(7) + f9(8) + f9(9) + f9(10) + f9(11) + f9(12)        &
             + f9(13) + f9(14)
        sum9 = f9(15) + f9(16) + f9(17) + f9(18)
        sum10 = f9(3) + f9(4) - f9(5) - f9(6)
        sum11 = f9(7) + f9(8) + f9(9) + f9(10) - f9(11) - f9(12)       &
              - f9(13) - f9(14)

        evlm1 = -30.0*f9(0) + coef2*sum1 + coef3*sum2
        evlm6 = coef5*sum6 - sum7 + sum8 - coef5*sum9
        evlm8 = sum10 + sum11
        evlm10 = f9(7) - f9(8) - f9(9) + f9(10)
        evlm11 = f9(15) - f9(16) - f9(17) + f9(18)
        evlm12 = f9(11) - f9(12) - f9(13) + f9(14)

        neqm1 = s1*(evlm1 - eqm1)
        neqm9 = s9*(evlm6 - eqm6)
        neqm11 = s9*(evlm8 - eqm8)
        neqm13 = s9*(evlm10 - eqm10)
        neqm14 = s9*(evlm11 - eqm11)
        neqm15 = s9*(evlm12 - eqm12)

        Sxx = -(neqm1 + 19.0*neqm9) / 38.0
        Syy = -(2.0*neqm1 - 19.0*(neqm9 - 3.0*neqm11)) / 76.0
        Szz = -(2.0*neqm1 - 19.0*(neqm9 + 3.0*neqm11)) / 76.0
        Sxy = -1.5*neqm13
        Syz = -1.5*neqm14
        Szx = -1.5*neqm15

        Sij2_m1(ix,iy,iz) = Sxx*Sxx + Syy*Syy + Szz*Szz                &
                    + 2.0*(Sxy*Sxy + Syz*Syz + Szx*Szx)
      end if
      end do
      end do
      end do

! to calculate Sij*Sij as a local array of (lx,ly,lz),
! Method 2: calculate in the discrete velocity space, f-space.
! Sij = -3/(2*rho0*tau)*Sigma_alpha{f_alpha(x,t)-f^(0)_alpha(x,t)}*e_alpha_i*
! e_alpha_j, where alpha = 0, 1, 2, ..., 18, f^(0) = f^(eq), the equilibrium
! distribution function, i,j = x, y, z

      do iz = 1,lz
      do iy = 1,ly
      do ix = 1,lx
      if(ibnodes(ix,iy,iz) < 0)then

        Sxx = 0.0
        Syy = 0.0
        Szz = 0.0
        Sxy = 0.0
        Syz = 0.0
        Szx = 0.0

        f9 = f(:,ix,iy,iz)
        rho9 = rho(ix,iy,iz)
        ux9 = ux(ix,iy,iz)
        uy9 = uy(ix,iy,iz)
        uz9 = uz(ix,iy,iz)
       
        usqr = 1.5*(ux9*ux9 + uy9*uy9 + uz9*uz9)   

        do ip = 1,npop-1
! first, calculate the equilibrium distribution function
          if(ip <= 6) wght = ww1  
          if(ip > 6) wght = ww2

          edtu = cix(ip)*ux9 + ciy(ip)*uy9 + ciz(ip)*uz9   
          feq9 = wght*(rho9 + 3.0*edtu + 4.5*edtu**2 - usqr)

          fneq9 = f9(ip) - feq9

          Sxx = Sxx + fneq9*real(cix(ip)*cix(ip))
          Syy = Syy + fneq9*real(ciy(ip)*ciy(ip))
          Szz = Szz + fneq9*real(ciz(ip)*ciz(ip))
          Sxy = Sxy + fneq9*real(cix(ip)*ciy(ip))
          Syz = Syz + fneq9*real(ciy(ip)*ciz(ip))
          Szx = Szx + fneq9*real(ciz(ip)*cix(ip))
        end do

        Sij2m2_tmp = Sxx*Sxx + Syy*Syy + Szz*Szz                &
                    + 2.0*(Sxy*Sxy + Syz*Syz + Szx*Szx)

        Sij2_m2(ix,iy,iz) = Sij2m2_tmp*(1.5 / tau)**2
      end if
      end do
      end do
      end do

! to calculate urms*urms as a local array of (lx,ly,lz)
! first calculate fluid rms velocity
      nf0 = count(ibnodes < 0)
      call MPI_REDUCE(nf0,nf,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      uxt0 = sum(ux, MASK = (ibnodes < 0))
      call MPI_REDUCE(uxt0,uxt,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      uyt0 = sum(uy, MASK = (ibnodes < 0))
      call MPI_REDUCE(uyt0,uyt,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      uzt0 = sum(uz, MASK = (ibnodes < 0))
      call MPI_REDUCE(uzt0,uzt,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      if(myid == 0)then
        uxmn = uxt / real(nf)
        uymn = uyt / real(nf)
        uzmn = uzt / real(nf)
      end if

      call MPI_BCAST(uxmn,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(uymn,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(uzmn,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)


      allocate (uprm2(lx,ly,lz))      
      uprm2 = ((ux - uxmn)**2 + (uy - uymn)**2 + (uz - uzmn)**2) / 3.0

! to calculate the averaged value of Sij2_m1, Sij2_m2, and uprm2
! over all the fluid nodes in the domain
      Sij2_m1t0 = sum(Sij2_m1, MASK = (ibnodes < 0))
      call MPI_REDUCE(Sij2_m1t0,Sij2_m1t,1,MPI_REAL8,MPI_SUM,0,         &
                      MPI_COMM_WORLD,ierr)

      Sij2_m2t0 = sum(Sij2_m2, MASK = (ibnodes < 0))
      call MPI_REDUCE(Sij2_m2t0,Sij2_m2t,1,MPI_REAL8,MPI_SUM,0,         &
                      MPI_COMM_WORLD,ierr)

      uprm2t0 = sum(uprm2, MASK = (ibnodes < 0))
      call MPI_REDUCE(uprm2t0,uprm2t,1,MPI_REAL8,MPI_SUM,0,             &
                      MPI_COMM_WORLD,ierr)

      if(myid == 0)then
        Sij2_m1_avg = Sij2_m1t / real(nf)
        Sij2_m2_avg = Sij2_m2t / real(nf)
        uprm2_avg = uprm2t / real(nf)
      end if

! to calculate SijSij and urms*urms as a function of (r/a), need to
! perform volume-averaging in sperical shells around a particle
      if(ipart .and. istep >= irelease)then

      rrbnd = 10.0*rad
      drr =0.05*rad
      irrbnd = int((rrbnd - rad) / drr) + 1

      allocate (Sij2p0_m1(irrbnd,npart))
      allocate (Sij2p0_m2(irrbnd,npart))
      allocate (uprm2p0(irrbnd,npart))
      allocate (icnt0(irrbnd,npart))
      allocate (Sij2p_m1(irrbnd,npart))
      allocate (Sij2p_m2(irrbnd,npart))
      allocate (uprm2p(irrbnd,npart))
      allocate (icnt(irrbnd,npart))

      Sij2p0_m1 = 0.0
      Sij2p0_m2 = 0.0
      uprm2p0 = 0.0
      icnt0 = 0

      do k = 1,lz
      do j = 1,ly
      do i = 1,lx
      if(ibnodes(i,j,k) < 0)then

        xpnt = real(i) - 0.5
        ypnt = real(j) - 0.5
        zpnt = real(k) - 0.5 + real(myid*lz)

        rrmin = real(lx)
        iprrmin = -1
        do ip = 1,npart
          xc = ypglb(1,ip)
          yc = ypglb(2,ip)
          zc = ypglb(3,ip)

! use the nearest particle center instead of the real center
          if((xc - xpnt) > real(lxh)) xc = xc - real(lx)
          if((xc - xpnt) < -real(lxh)) xc = xc + real(lx)

          if((yc - ypnt) > real(nyh)) yc = yc - real(ny)
          if((yc - ypnt) < -real(nyh)) yc = yc + real(ny)

          if((zc - zpnt) > real(nzh)) zc = zc - real(nz)
          if((zc - zpnt) < -real(nzh)) zc = zc + real(nz)

          xx0 = xpnt - xc
          yy0 = ypnt - yc
          zz0 = zpnt - zc

          rr0 = sqrt(xx0*xx0 + yy0*yy0 + zz0*zz0)
          if(rr0 <= rrbnd .and. rr0 < rrmin)then
            rrmin = rr0
            iprrmin = ip
          end if
        end do

        if(iprrmin > 0)then
          irr0 = int((rrmin - rad) / drr) + 1
          Sij2p0_m1(irr0,iprrmin) = Sij2p0_m1(irr0,iprrmin) + Sij2_m1(i,j,k)
          Sij2p0_m2(irr0,iprrmin) = Sij2p0_m2(irr0,iprrmin) + Sij2_m2(i,j,k)
          uprm2p0(irr0,iprrmin) = uprm2p0(irr0,iprrmin) + uprm2(i,j,k)
          icnt0(irr0,iprrmin) = icnt0(irr0,iprrmin) + 1
        end if

      end if
      end do
      end do
      end do

      deallocate (uprm2)         

! now collect info. for Sij2p and icnt
      ilen = irrbnd*npart

      call MPI_ALLREDUCE(Sij2p0_m1,Sij2p_m1,ilen,MPI_REAL8,MPI_SUM,     &
                         MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(Sij2p0_m2,Sij2p_m2,ilen,MPI_REAL8,MPI_SUM,     &
                         MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(uprm2p0,uprm2p,ilen,MPI_REAL8,MPI_SUM,         &
                         MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(icnt0,icnt,ilen,MPI_INTEGER,MPI_SUM,          &
                         MPI_COMM_WORLD,ierr)

      deallocate (Sij2p0_m1)
      deallocate (Sij2p0_m2)
      deallocate (uprm2p0)
      deallocate (icnt0)

! then average the Sij2 for each particle in each spherical shell layer
      if(myid == 0)then
        ttt = real(istep)*tscale
! normalized by eddy turnover time at t = 0
        ttt = ttt / et0

        istp1 = istep / 100000
        istp2 = mod(istep,100000) / 10000
        istp3 = mod(istep,10000) / 1000
        istp4 = mod(istep,1000) / 100
        istp5 = mod(istep,100) / 10
        istp6 = mod(istep,10)

        fnm = trim(dirstat)//'sij'                                     &
            //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)     &
            //char(istp4 + 48)//char(istp5 + 48)//char(istp6 + 48)     &
            //'.dat'

        open(31, file = trim(fnm), status = 'unknown',                 &
                 form = 'formatted', position = 'append')

        do is = 1,irrbnd
          Sij2tt_m1 = sum(Sij2p_m1(is,:))
          Sij2tt_m2 = sum(Sij2p_m2(is,:))
          uprm2tt = sum(uprm2p(is,:))
          icnttt = sum(icnt(is,:))
          if(icnttt == 0)then
            Sij2mn_m1 = 0.0
            Sij2mn_m2 = 0.0
            uprm2mn = 0.0
          else
            Sij2mn_m1 = Sij2tt_m1 / real(icnttt)
            Sij2mn_m2 = Sij2tt_m2 / real(icnttt)
            uprm2mn = uprm2tt / real(icnttt)
          end if
          rrmn = real(is)*drr / rad + 1.0
          Sij2_m1_norm = Sij2mn_m1 / Sij2_m1_avg
          Sij2_m2_norm = Sij2mn_m2 / Sij2_m2_avg
          uprm2_norm = uprm2mn / uprm2_avg
          write(31,310) ttt, rrmn, Sij2mn_m1, Sij2mn_m2, uprm2mn,      &
                        Sij2_m1_avg, Sij2_m2_avg, uprm2_avg,           &
                        Sij2_m1_norm, Sij2_m2_norm, uprm2_norm
        end do

      end if
310   format(2x,11(1pe16.6))

      deallocate (Sij2p_m1)
      deallocate (Sij2p_m2)
      deallocate (uprm2p)
      deallocate (icnt)

      end if

      end subroutine sijstat02
!===========================================================================
! This subroutine "sijstat" is essentially the same as "sijstat02" above.
! The goal is to optimize memory usage. In the case of 409600 particles, 
! the particle realted arrays of dimension(180,npart) as follows
! can be of huge size with employment of 512 processors.  
! They require the memory access of more than 1.1 TB! 
!      integer, allocatable, dimension(:,:):: icnt0, icnt
!      real, allocatable, dimension(:,:)::Sij2p0_m1, Sij2p_m1
!      real, allocatable, dimension(:,:)::Sij2p0_m2, Sij2p_m2
!      real, allocatable, dimension(:,:):: uprm2p0, uprm2p

! In this subroutine, these 2D arrays are replaced by their 1D conuterparts
! of dimension (180), as the Sij profile of each particle is not necessary.
! Also, the calculation of Sij2_m1(lx,ly,lz) and Sij2_m2(lx,ly,lz) are 
! merged in a same loop instead of two loops, hence saving some CPU time.
! Moreover, the searching of the nearest particle center to a specific 
! fluid node was improved, in the way that this fluid node should be located
! within 10*rad_particle distance to the particle center being considered.

      subroutine sijstat03
      use mpi
      use var_inc
      implicit none

      integer ix, iy, iz
      real, dimension(0:npop-1) :: f9
      real, dimension(lx,ly,lz) :: Sij2_m1, Sij2_m2, uprm2, vort2
      real rho9, ux9, uy9, uz9, ux9s, uy9s, uz9s
      real eqm1, eqm6, eqm8, eqm10, eqm11, eqm12
      real sum1, sum2, sum6, sum7, sum8, sum9, sum10, sum11
      real evlm1, evlm6, evlm8, evlm10, evlm11, evlm12
      real neqm1, neqm9, neqm11, neqm13, neqm14, neqm15
      real Sxx, Syy, Szz, Sxy, Syz, Szx
      real usqr, edtu, wght, feq9, fneq9

      real, dimension(nx,ny):: tke9,diss_m1_9,diss_m2_9
      real, dimension(nx,ny):: tke9c,diss_m1_9c,diss_m2_9c
      character (len = 100):: fnm1, fnm2, fnm3

      integer iroot,sly,i,j
      integer, allocatable, dimension(:):: icnt0, icnt
      real, allocatable, dimension(:):: Sij2p0_m1, Sij2p_m1
      real, allocatable, dimension(:):: Sij2p0_m2, Sij2p_m2
      real, allocatable, dimension(:):: uprm2p0, uprm2p
      real, allocatable, dimension(:):: vort2p0, vort2p

      integer ip, ilen, is
      integer irrbnd, irr0, icnttt, iprrmin
      integer istp1, istp2, istp3, istp4, istp5, istp6
      real xc, yc, zc, xpnt, ypnt, zpnt, rrbnd, drr, rrmin, rrbnds
      real xx0, yy0, zz0, rr0, ttt
      real Sij2tt_m1, Sij2tt_m2, Sij2mn_m1, Sij2mn_m2, rrmn
      real Sij2_m1t0, Sij2_m2t0, Sij2_m1t, Sij2_m2t
      real Sij2_m1_avg, Sij2_m2_avg, Sij2m2_tmp
      real Sij2_m1_norm, Sij2_m2_norm


      integer nf0, nf
      real uxt0, uyt0, uzt0, uxt, uyt, uzt
      real uxmn, uymn, uzmn
      real uprm2tt, uprm2mn
      real uprm2t0, uprm2t, uprm2_avg, uprm2_norm
      real vort2tt, vort2mn
      real vort2t0, vort2t, vort2_avg, vort2_norm

      character (len = 100):: fnm

      do iz = 1,lz
      do iy = 1,ly
      do ix = 1,lx
      if(ibnodes(ix,iy,iz) < 0)then
! to calculate Sij*Sij as a local array of (lx,ly,lz),

! Method 1: calculate in the moment space.
! see Yu H. et al. Computers & Fluids 35, pp. 957-965, 2006, Appendix.

        f9 = f(:,ix,iy,iz)

        rho9 = rho(ix,iy,iz)
        ux9 = ux(ix,iy,iz)
        uy9 = uy(ix,iy,iz)
        uz9 = uz(ix,iy,iz)
        ux9s = ux9*ux9
        uy9s = uy9*uy9
        uz9s = uz9*uz9

        eqm1 = -11.0*rho9 + 19.0*(ux9s + uy9s + uz9s)
        eqm6 = 2.0*ux9s - uy9s - uz9s
        eqm8 = uy9s - uz9s
        eqm10 = ux9*uy9
        eqm11 = uy9*uz9
        eqm12 = ux9*uz9

        sum1 = f9(1) + f9(2) + f9(3) + f9(4) + f9(5) + f9(6)
        sum2 = f9(7) + f9(8) + f9(9) + f9(10) + f9(11) + f9(12)        &
             + f9(13) + f9(14) + f9(15) + f9(16) + f9(17) + f9(18)
        sum6 = f9(1) + f9(2)
        sum7 = f9(3) + f9(4) + f9(5) + f9(6)
        sum8 = f9(7) + f9(8) + f9(9) + f9(10) + f9(11) + f9(12)        &
             + f9(13) + f9(14)
        sum9 = f9(15) + f9(16) + f9(17) + f9(18)
        sum10 = f9(3) + f9(4) - f9(5) - f9(6)
        sum11 = f9(7) + f9(8) + f9(9) + f9(10) - f9(11) - f9(12)       &
              - f9(13) - f9(14)

        evlm1 = -30.0*f9(0) + coef2*sum1 + coef3*sum2
        evlm6 = coef5*sum6 - sum7 + sum8 - coef5*sum9
        evlm8 = sum10 + sum11
        evlm10 = f9(7) - f9(8) - f9(9) + f9(10)
        evlm11 = f9(15) - f9(16) - f9(17) + f9(18)
        evlm12 = f9(11) - f9(12) - f9(13) + f9(14)

        neqm1 = s1*(evlm1 - eqm1)
        neqm9 = s9*(evlm6 - eqm6)
        neqm11 = s9*(evlm8 - eqm8)
        neqm13 = s9*(evlm10 - eqm10)
        neqm14 = s9*(evlm11 - eqm11)
        neqm15 = s9*(evlm12 - eqm12)

        Sxx = -(neqm1 + 19.0*neqm9) / 38.0
        Syy = -(2.0*neqm1 - 19.0*(neqm9 - 3.0*neqm11)) / 76.0
        Szz = -(2.0*neqm1 - 19.0*(neqm9 + 3.0*neqm11)) / 76.0
        Sxy = -1.5*neqm13
        Syz = -1.5*neqm14
        Szx = -1.5*neqm15

        Sij2_m1(ix,iy,iz) = Sxx*Sxx + Syy*Syy + Szz*Szz                &
                    + 2.0*(Sxy*Sxy + Syz*Syz + Szx*Szx)

! Method 2: calculate in the discrete velocity space, f-space.
! Sij = -3/(2*rho0*tau)*Sigma_alpha{f_alpha(x,t)-f^(0)_alpha(x,t)}*e_alpha_i*
! e_alpha_j, where alpha = 0, 1, 2, ..., 18, f^(0) = f^(eq), the equilibrium
! distribution function, i,j = x, y, z

        Sxx = 0.0
        Syy = 0.0
        Szz = 0.0
        Sxy = 0.0
        Syz = 0.0
        Szx = 0.0

        usqr = 1.5*(ux9s + uy9s + uz9s)

        do ip = 1,npop-1
! first, calculate the equilibrium distribution function
          if(ip <= 6) wght = ww1
          if(ip > 6) wght = ww2

          edtu = cix(ip)*ux9 + ciy(ip)*uy9 + ciz(ip)*uz9
          feq9 = wght*(rho9 + 3.0*edtu + 4.5*edtu**2 - usqr)

          fneq9 = f9(ip) - feq9

          Sxx = Sxx + fneq9*real(cix(ip)*cix(ip))
          Syy = Syy + fneq9*real(ciy(ip)*ciy(ip))
          Szz = Szz + fneq9*real(ciz(ip)*ciz(ip))
          Sxy = Sxy + fneq9*real(cix(ip)*ciy(ip))
          Syz = Syz + fneq9*real(ciy(ip)*ciz(ip))
          Szx = Szx + fneq9*real(ciz(ip)*cix(ip))
        end do

        Sij2m2_tmp = Sxx*Sxx + Syy*Syy + Szz*Szz                &
                    + 2.0*(Sxy*Sxy + Syz*Syz + Szx*Szx)
        Sij2_m2(ix,iy,iz) = Sij2m2_tmp*(1.5 / tau)**2
      end if
      end do
      end do
      end do
      if(myid == 0)write(*,*)'pass sijstat03 step 1'

! calculate vorticity squared as local array
      call vortcalc
      vort2 = ox*ox + oy*oy + oz*oz

! to calculate urms*urms as a local array of (lx,ly,lz)
! first calculate fluid rms velocity
      nf0 = count(ibnodes < 0)
      call MPI_REDUCE(nf0,nf,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      uxt0 = sum(ux, MASK = (ibnodes < 0))
      call MPI_REDUCE(uxt0,uxt,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      uyt0 = sum(uy, MASK = (ibnodes < 0))
      call MPI_REDUCE(uyt0,uyt,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      uzt0 = sum(uz, MASK = (ibnodes < 0))
      call MPI_REDUCE(uzt0,uzt,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      if(myid == 0)then
        uxmn = uxt / real(nf)
        uymn = uyt / real(nf)
        uzmn = uzt / real(nf)
      end if

      call MPI_BCAST(uxmn,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(uymn,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(uzmn,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

      uprm2 = ((ux - uxmn)**2 + (uy - uymn)**2 + (uz - uzmn)**2) / 3.0

! to calculate the averaged value of Sij2_m1, Sij2_m2, and uprm2
! over all the fluid nodes in the domain
      Sij2_m1t0 = sum(Sij2_m1, MASK = (ibnodes < 0))
      call MPI_REDUCE(Sij2_m1t0,Sij2_m1t,1,MPI_REAL8,MPI_SUM,0,         &
                      MPI_COMM_WORLD,ierr)

      Sij2_m2t0 = sum(Sij2_m2, MASK = (ibnodes < 0))
      call MPI_REDUCE(Sij2_m2t0,Sij2_m2t,1,MPI_REAL8,MPI_SUM,0,         &
                      MPI_COMM_WORLD,ierr)

      uprm2t0 = sum(uprm2, MASK = (ibnodes < 0))
      call MPI_REDUCE(uprm2t0,uprm2t,1,MPI_REAL8,MPI_SUM,0,             &
                      MPI_COMM_WORLD,ierr)

      vort2t0 = sum(vort2, MASK = (ibnodes < 0))
      call MPI_REDUCE(vort2t0,vort2t,1,MPI_REAL8,MPI_SUM,0,             &
                      MPI_COMM_WORLD,ierr)

      if(myid == 0)then
        Sij2_m1_avg = Sij2_m1t / real(nf)
        Sij2_m2_avg = Sij2_m2t / real(nf)
        uprm2_avg = uprm2t / real(nf)
        vort2_avg = vort2t / real(nf)
      end if

      if(myid == 0)write(*,*)'pass sijstat03 step 2'

! to calculate SijSij and urms*urms as a function of (r/a), need to
! perform volume-averaging in sperical shells around a particle
      if(ipart .and. istep >= irelease)then

      rrbnd = 10.0*rad
!      rrbnd = 6.0*rad
      rrbnds = 3.0*rad
      drr =0.05*rad
      irrbnd = int((rrbnd - rad) / drr) + 1

      allocate (Sij2p0_m1(irrbnd))
      allocate (Sij2p0_m2(irrbnd))
      allocate (uprm2p0(irrbnd))
      allocate (vort2p0(irrbnd))
      allocate (icnt0(irrbnd))
      allocate (Sij2p_m1(irrbnd))
      allocate (Sij2p_m2(irrbnd))
      allocate (uprm2p(irrbnd))
      allocate (vort2p(irrbnd))
      allocate (icnt(irrbnd))

      Sij2p0_m1 = 0.0
      Sij2p0_m2 = 0.0
      uprm2p0 = 0.0
      vort2p0 = 0.0
      icnt0 = 0

      do iz = 1,lz
      do iy = 1,ly
      do ix = 1,lx
      if(ibnodes(ix,iy,iz) < 0)then

        xpnt = real(ix) - 0.5
        ypnt = real(iy) - 0.5 + real(indy*ly)
        zpnt = real(iz) - 0.5 + real(indz*lz)

        rrmin = real(lx)
        iprrmin = -1

        DO
        do ip = 1,npart
          xc = ypglb(1,ip)

! use the nearest particle center instead of the real center
          if((xc - xpnt) > real(nxh)) xc = xc - real(lx)
          if((xc - xpnt) < -real(nxh)) xc = xc + real(lx)

          xx0 = xpnt - xc
          if(abs(xx0) <= rrbnds)then

          yc = ypglb(2,ip)
          if((yc - ypnt) > real(nyh)) yc = yc - real(ny)
          if((yc - ypnt) < -real(nyh)) yc = yc + real(ny)

          yy0 = ypnt - yc
          if(abs(yy0) <= rrbnds)then

          zc = ypglb(3,ip)
          if((zc - zpnt) > real(nzh)) zc = zc - real(nz)
          if((zc - zpnt) < -real(nzh)) zc = zc + real(nz)

          zz0 = zpnt - zc
          if(abs(zz0) <= rrbnds)then

          rr0 = sqrt(xx0*xx0 + yy0*yy0 + zz0*zz0)
          if(rr0 <= rrbnds .and. rr0 < rrmin)then
            rrmin = rr0
            iprrmin = ip
          end if

          end if
          end if
          end if
        end do

        if(iprrmin > 0 .or. (rrbnds - rrbnd) >= 0.0) then
             exit
        else
          rrbnds = rrbnds + rad
        end if

        END DO

        if(iprrmin > 0)then
          irr0 = int((rrmin - rad) / drr) + 1
          Sij2p0_m1(irr0) = Sij2p0_m1(irr0) + Sij2_m1(ix,iy,iz)
          Sij2p0_m2(irr0) = Sij2p0_m2(irr0) + Sij2_m2(ix,iy,iz)
          uprm2p0(irr0) = uprm2p0(irr0) + uprm2(ix,iy,iz)
          vort2p0(irr0) = vort2p0(irr0) + vort2(ix,iy,iz)
          icnt0(irr0) = icnt0(irr0) + 1
        end if

      end if
      end do
      end do
      end do

! now collect info. for Sij2p and icnt
      ilen = irrbnd

      call MPI_ALLREDUCE(Sij2p0_m1,Sij2p_m1,ilen,MPI_REAL8,MPI_SUM,     &
                         MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(Sij2p0_m2,Sij2p_m2,ilen,MPI_REAL8,MPI_SUM,     &
                         MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(uprm2p0,uprm2p,ilen,MPI_REAL8,MPI_SUM,         &
                         MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(vort2p0,vort2p,ilen,MPI_REAL8,MPI_SUM,         &
                         MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(icnt0,icnt,ilen,MPI_INTEGER,MPI_SUM,          &
                         MPI_COMM_WORLD,ierr)

      deallocate (Sij2p0_m1)
      deallocate (Sij2p0_m2)
      deallocate (uprm2p0)
      deallocate (vort2p0)
      deallocate (icnt0)

! then average the Sij2 for each particle in each spherical shell layer
      if(myid == 0)then
        ttt = real(istep)*tscale
! normalized by eddy turnover time at t = 0
        ttt = ttt / et0

        istp1 = istep / 100000
        istp2 = mod(istep,100000) / 10000
        istp3 = mod(istep,10000) / 1000
        istp4 = mod(istep,1000) / 100
        istp5 = mod(istep,100) / 10
        istp6 = mod(istep,10)

        fnm = trim(dirstat)//'sij'                                     &
            //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)     &
            //char(istp4 + 48)//char(istp5 + 48)//char(istp6 + 48)     &
            //'.dat'

        open(31, file = trim(fnm), status = 'unknown',                 &
                 form = 'formatted', position = 'append')

        do is = 1,irrbnd
          Sij2tt_m1 = Sij2p_m1(is)
          Sij2tt_m2 = Sij2p_m2(is)
          uprm2tt = uprm2p(is)
          vort2tt = vort2p(is)
          icnttt = icnt(is)
          if(icnttt == 0)then
            Sij2mn_m1 = 0.0
            Sij2mn_m2 = 0.0
            uprm2mn = 0.0
            vort2mn = 0.0
          else
            Sij2mn_m1 = Sij2tt_m1 / real(icnttt)
            Sij2mn_m2 = Sij2tt_m2 / real(icnttt)
            uprm2mn = uprm2tt / real(icnttt)
            vort2mn = vort2tt / real(icnttt)
          end if
          rrmn = (real(is)-0.5)*drr / rad + 1.0
          Sij2_m1_norm = Sij2mn_m1 / Sij2_m1_avg
          Sij2_m2_norm = Sij2mn_m2 / Sij2_m2_avg
          uprm2_norm = uprm2mn / uprm2_avg
          vort2_norm = vort2mn / vort2_avg
          write(31,310) ttt, rrmn, real(icnttt), Sij2mn_m1, Sij2mn_m2, uprm2mn,vort2mn, &
                        Sij2_m1_avg, Sij2_m2_avg, uprm2_avg, vort2_avg,   &
                   Sij2_m1_norm, Sij2_m2_norm, uprm2_norm, vort2_norm
        end do

      close(31)

      end if
310   format(2x,15(1pe16.6))

      deallocate (Sij2p_m1)
      deallocate (Sij2p_m2)
      deallocate (uprm2p)
      deallocate (vort2p)
      deallocate (icnt)

      end if

      if(myid == 0)write(*,*)'pass sijstat03 step 3'

! Print out for visualization
      iroot = nprocZ/2

!Set vortivity to zero inside particles
      where(ibnodes > 0) uprm2 = 0.0
      where(ibnodes > 0) sij2_m1 = 0.0
      where(ibnodes > 0) sij2_m2 = 0.0

      tke9 = 0.0d0
      diss_m1_9 = 0.0d0
      diss_m2_9 = 0.0d0

      if(indz == iroot) then
      sly = indy*ly
      tke9(1:lx,sly+1:sly+ly)=sqrt(uprm2(1:lx,1:ly,1))/vscale
      diss_m1_9(1:lx,sly+1:sly+ly)=sqrt(sij2_m1(1:lx,1:ly,1))/tscale
      diss_m2_9(1:lx,sly+1:sly+ly)=sqrt(sij2_m2(1:lx,1:ly,1))/tscale
      end if

! Merge the field
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(tke9,tke9c,nx*ny,MPI_REAL8,MPI_SUM,           &
                         MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(diss_m1_9,diss_m1_9c,nx*ny,MPI_REAL8,MPI_SUM, &
                         MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(diss_m2_9,diss_m2_9c,nx*ny,MPI_REAL8,MPI_SUM,     &
                         MPI_COMM_WORLD,ierr)

      if(myid == 0)then

        fnm1 = trim(dirmoviedata)//'Sij_m1_z2565.dat'
        fnm2 = trim(dirmoviedata)//'Sij_m2_z2565.dat'
        fnm3 = trim(dirmoviedata)//'urms_z2565.dat'

        open(28, file = trim(fnm1), status = 'unknown',                &
                 form = 'formatted', position = 'append')
        open(29, file = trim(fnm2), status = 'unknown',                &
                 form = 'formatted', position = 'append')
        open(30, file = trim(fnm3), status = 'unknown',                &
                 form = 'formatted', position = 'append')

          write(28,282) ((diss_m1_9c(i,j),i=1,nx),j=1,ny)
          write(29,282) ((diss_m2_9c(i,j),i=1,nx),j=1,ny)
          write(30,282) ((tke9c(i,j),i=1,nx),j=1,ny)
        close(28)
        close(29)
        close(30)

      end if
      if(myid == 0)write(*,*)'pass sijstat03 step 4'

282   format(2x,8(1pe12.4))

      end subroutine sijstat03
!===========================================================================
! This subroutine is the same as "sijstat03" above, except that it only
! calculate till the Sij2_m1, Sij2_m2, uprm2, and then save these arrays
! for further processing.
      subroutine sijstat
      use mpi
      use var_inc
      implicit none

      integer ix, iy, iz
      real, dimension(0:npop-1) :: f9
      real, dimension(lx,ly,lz) :: Sij2_m1, Sij2_m2, uprm2
      real rho9, ux9, uy9, uz9, ux9s, uy9s, uz9s
      real eqm1, eqm6, eqm8, eqm10, eqm11, eqm12
      real sum1, sum2, sum6, sum7, sum8, sum9, sum10, sum11
      real evlm1, evlm6, evlm8, evlm10, evlm11, evlm12
      real neqm1, neqm9, neqm11, neqm13, neqm14, neqm15
      real Sxx, Syy, Szz, Sxy, Syz, Szx
      real usqr, edtu, wght, feq9, fneq9

      integer iprc1, iprc2, iprc3
      integer istp1, istp2, istp3, istp4, istp5, istp6

      real Sij2_m1t0, Sij2_m2t0, Sij2_m1t, Sij2_m2t
      real Sij2_m1_avg, Sij2_m2_avg, Sij2m2_tmp

      integer nf0, nf, ip
      real uxt0, uyt0, uzt0, uxt, uyt, uzt
      real uxmn, uymn, uzmn
      real uprm2t0, uprm2t, uprm2_avg

      character (len = 100):: fnm

      Sij2_m1 = 0.0  
      Sij2_m2 = 0.0  

      do iz = 1,lz
      do iy = 1,ly
      do ix = 1,lx
      if(ibnodes(ix,iy,iz) < 0)then
! to calculate Sij*Sij as a local array of (lx,ly,lz),

! Method 1: calculate in the moment space.
! see Yu H. et al. Computers & Fluids 35, pp. 957-965, 2006, Appendix.

        f9 = f(:,ix,iy,iz)

        rho9 = rho(ix,iy,iz)
        ux9 = ux(ix,iy,iz)
        uy9 = uy(ix,iy,iz)
        uz9 = uz(ix,iy,iz)
        ux9s = ux9*ux9
        uy9s = uy9*uy9
        uz9s = uz9*uz9

        eqm1 = -11.0*rho9 + 19.0*(ux9s + uy9s + uz9s)
        eqm6 = 2.0*ux9s - uy9s - uz9s
        eqm8 = uy9s - uz9s
        eqm10 = ux9*uy9
        eqm11 = uy9*uz9
        eqm12 = ux9*uz9

        sum1 = f9(1) + f9(2) + f9(3) + f9(4) + f9(5) + f9(6)
        sum2 = f9(7) + f9(8) + f9(9) + f9(10) + f9(11) + f9(12)        &
             + f9(13) + f9(14) + f9(15) + f9(16) + f9(17) + f9(18)
        sum6 = f9(1) + f9(2)
        sum7 = f9(3) + f9(4) + f9(5) + f9(6)
        sum8 = f9(7) + f9(8) + f9(9) + f9(10) + f9(11) + f9(12)        &
             + f9(13) + f9(14)
        sum9 = f9(15) + f9(16) + f9(17) + f9(18)
        sum10 = f9(3) + f9(4) - f9(5) - f9(6)
        sum11 = f9(7) + f9(8) + f9(9) + f9(10) - f9(11) - f9(12)       &
              - f9(13) - f9(14)

        evlm1 = -30.0*f9(0) + coef2*sum1 + coef3*sum2
        evlm6 = coef5*sum6 - sum7 + sum8 - coef5*sum9
        evlm8 = sum10 + sum11
        evlm10 = f9(7) - f9(8) - f9(9) + f9(10)
        evlm11 = f9(15) - f9(16) - f9(17) + f9(18)
        evlm12 = f9(11) - f9(12) - f9(13) + f9(14)

        neqm1 = s1*(evlm1 - eqm1)
        neqm9 = s9*(evlm6 - eqm6)
        neqm11 = s9*(evlm8 - eqm8)
        neqm13 = s9*(evlm10 - eqm10)
        neqm14 = s9*(evlm11 - eqm11)
        neqm15 = s9*(evlm12 - eqm12)

        Sxx = -(neqm1 + 19.0*neqm9) / 38.0
        Syy = -(2.0*neqm1 - 19.0*(neqm9 - 3.0*neqm11)) / 76.0
        Szz = -(2.0*neqm1 - 19.0*(neqm9 + 3.0*neqm11)) / 76.0
        Sxy = -1.5*neqm13
        Syz = -1.5*neqm14
        Szx = -1.5*neqm15

        Sij2_m1(ix,iy,iz) = Sxx*Sxx + Syy*Syy + Szz*Szz                &
                    + 2.0*(Sxy*Sxy + Syz*Syz + Szx*Szx)

! Method 2: calculate in the discrete velocity space, f-space.
! Sij = -3/(2*rho0*tau)*Sigma_alpha{f_alpha(x,t)-f^(0)_alpha(x,t)}*e_alpha_i*
! e_alpha_j, where alpha = 0, 1, 2, ..., 18, f^(0) = f^(eq), the equilibrium
! distribution function, i,j = x, y, z

        Sxx = 0.0
        Syy = 0.0
        Szz = 0.0
        Sxy = 0.0
        Syz = 0.0
        Szx = 0.0

        usqr = 1.5*(ux9s + uy9s + uz9s)

        do ip = 1,npop-1
! first, calculate the equilibrium distribution function
          if(ip <= 6) wght = ww1
          if(ip > 6) wght = ww2

          edtu = cix(ip)*ux9 + ciy(ip)*uy9 + ciz(ip)*uz9
          feq9 = wght*(rho9 + 3.0*edtu + 4.5*edtu**2 - usqr)

          fneq9 = f9(ip) - feq9

          Sxx = Sxx + fneq9*real(cix(ip)*cix(ip))
          Syy = Syy + fneq9*real(ciy(ip)*ciy(ip))
          Szz = Szz + fneq9*real(ciz(ip)*ciz(ip))
          Sxy = Sxy + fneq9*real(cix(ip)*ciy(ip))
          Syz = Syz + fneq9*real(ciy(ip)*ciz(ip))
          Szx = Szx + fneq9*real(ciz(ip)*cix(ip))
        end do

        Sij2m2_tmp = Sxx*Sxx + Syy*Syy + Szz*Szz                       &
                    + 2.0*(Sxy*Sxy + Syz*Syz + Szx*Szx)

        Sij2_m2(ix,iy,iz) = Sij2m2_tmp*(1.5 / tau)**2
      end if
      end do
      end do
      end do

! to calculate urms*urms as a local array of (lx,ly,lz)
! first calculate fluid rms velocity
      nf0 = count(ibnodes < 0)
      call MPI_REDUCE(nf0,nf,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      uxt0 = sum(ux, MASK = (ibnodes < 0))
      call MPI_REDUCE(uxt0,uxt,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      uyt0 = sum(uy, MASK = (ibnodes < 0))
      call MPI_REDUCE(uyt0,uyt,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      uzt0 = sum(uz, MASK = (ibnodes < 0))
      call MPI_REDUCE(uzt0,uzt,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      if(myid == 0)then
        uxmn = uxt / real(nf)
        uymn = uyt / real(nf)
        uzmn = uzt / real(nf)
      end if

      call MPI_BCAST(uxmn,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(uymn,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(uzmn,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

      uprm2 = ((ux - uxmn)**2 + (uy - uymn)**2 + (uz - uzmn)**2) / 3.0

! to calculate the averaged value of Sij2_m1, Sij2_m2, and uprm2
! over all the fluid nodes in the domain
      Sij2_m1t0 = sum(Sij2_m1, MASK = (ibnodes < 0))
      call MPI_REDUCE(Sij2_m1t0,Sij2_m1t,1,MPI_REAL8,MPI_SUM,0,         &
                      MPI_COMM_WORLD,ierr)

      Sij2_m2t0 = sum(Sij2_m2, MASK = (ibnodes < 0))
      call MPI_REDUCE(Sij2_m2t0,Sij2_m2t,1,MPI_REAL8,MPI_SUM,0,         &
                      MPI_COMM_WORLD,ierr)

      uprm2t0 = sum(uprm2, MASK = (ibnodes < 0))
      call MPI_REDUCE(uprm2t0,uprm2t,1,MPI_REAL8,MPI_SUM,0,             &
                      MPI_COMM_WORLD,ierr)

      if(myid == 0)then
        Sij2_m1_avg = Sij2_m1t / real(nf)
        Sij2_m2_avg = Sij2_m2t / real(nf)
        uprm2_avg = uprm2t / real(nf)
      end if


      IF(ipart .and. istep >= irelease) THEN

      iprc1 = myid / 100
      iprc2 = mod(myid,100) / 10
      iprc3 = mod(myid,10)

      istp1 = istep / 100000
      istp2 = mod(istep,100000) / 10000
      istp3 = mod(istep,10000) / 1000
      istp4 = mod(istep,1000) / 100
      istp5 = mod(istep,100) / 10
      istp6 = mod(istep,10)

      fnm = trim(dirstat)//'sijdata/sij.'                            &
          //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)     &
          //char(istp4 + 48)//char(istp5 + 48)//char(istp6 + 48)     &
          //'.'                                                      &
          //char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)

      open(31, file = trim(fnm), status = 'unknown',                 &
               form = 'unformatted')

      write(31) Sij2_m1, Sij2_m2, uprm2, ibnodes  
      if(myid == 0) write(31) ypglb, Sij2_m1_avg, Sij2_m2_avg, uprm2_avg


      close(31)

      END IF

      end subroutine sijstat
!==========================================================================
      subroutine vortcalc
      use mpi
      use var_inc
      implicit none

      real, dimension(lx+2,lly,lz) :: tmp

      vx = 0.0
      vy = 0.0
      vz = 0.0

      vx(1:lx,1:ly,:) = ux
      vy(1:lx,1:ly,:) = uy
      vz(1:lx,1:ly,:) = uz

      call mpifft3DRC(vx)
      call mpifft3DRC(vy)
      call mpifft3DRC(vz)

      wx = ky*vz - kz*vy
      wy = kz*vx - kx*vz
      wz = kx*vy - ky*vx
      tmp = wx
      wx(1:lx , 1:lly-1:2 , :) = -tmp(1:lx , 2:lly:2   , :)
      wx(1:lx , 2:lly:2   , :) =  tmp(1:lx , 1:lly-1:2 , :)
      tmp = wy
      wy(1:lx , 1:lly-1:2 , :) = -tmp(1:lx , 2:lly:2   , :)
      wy(1:lx , 2:lly:2   , :) =  tmp(1:lx , 1:lly-1:2 , :)
      tmp = wz
      wz(1:lx , 1:lly-1:2 , :) = -tmp(1:lx , 2:lly:2   , :)
      wz(1:lx , 2:lly:2   , :) =  tmp(1:lx , 1:lly-1:2 , :)

      call mpifft3DCR(wx)
      call mpifft3DCR(wy)
      call mpifft3DCR(wz)

      ox = wx(1:lx,1:ly,:)
      oy = wy(1:lx,1:ly,:)
      oz = wz(1:lx,1:ly,:)

      end subroutine vortcalc
!===========================================================================

!===========================================================================
! This subroutine is for debugging overlap problem and can be deleted after
      subroutine writepartpair 
      use mpi
      use var_inc
      implicit none


      end subroutine writepartpair
!===========================================================================


!==========================================================================
! This subroutine is for writing out a slice of the instantaneous flow field
      subroutine writeflowfieldstart
      use mpi
      use var_inc
      implicit none

      real,dimension(ly,lz)::uxplane
      integer i,k,j
      integer iprc1, iprc2, iprc3
      character (len = 100):: fnm1

        iprc1 = myid / 100
        iprc2 = mod(myid,100) / 10
        iprc3 = mod(myid,10)

        fnm1 = trim(dirstat)//'flowfield2Dstart.dat.'         &
              //char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)

        open(44, file = trim(fnm1), status = 'unknown',                 &
                 form = 'unformatted')

          uxplane = ux(lx/2,:,:)

          write(44) uxplane

       close(44)

      end subroutine writeflowfieldstart
!===========================================================================

!==========================================================================
! This subroutine is for writing out a slice of the instantaneous flow field
      subroutine writeflowfield
      use mpi
      use var_inc
      implicit none

      real,dimension(ly,lz)::uxplane
      integer i,k,j
      integer iprc1, iprc2, iprc3
      character (len = 100):: fnm1

        iprc1 = myid / 100
        iprc2 = mod(myid,100) / 10
        iprc3 = mod(myid,10)

        fnm1 = trim(dirstat)//'flowfield2D.dat.'         &
              //char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)

        open(66, file = trim(fnm1), status = 'unknown',                 &
                 form = 'unformatted')

          uxplane = ux(lx/2,:,:)

          write(66) uxplane

          close(66)

      end subroutine writeflowfield
!===========================================================================

!==========================================================================
! This subroutine is for writing out a slice of the instantaneous flow field
      subroutine writeflowfieldend
      use mpi
      use var_inc
      implicit none

      real,dimension(ly,lz)::uxplane
      integer i,k,j
      integer iprc1, iprc2, iprc3
      character (len = 100):: fnm1

        iprc1 = myid / 100
        iprc2 = mod(myid,100) / 10
        iprc3 = mod(myid,10)

        fnm1 = trim(dirstat)//'flowfield2Dend.dat.'         &
              //char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)

        open(55, file = trim(fnm1), status = 'unknown',                 &
                 form = 'unformatted')

          uxplane = ux(lx/2,:,:)

          write(55) uxplane

          close(55)

      end subroutine writeflowfieldend
!===========================================================================

!===========================================================================
      subroutine colldata(fnm2) 
      use mpi
      use var_inc
      implicit none

      integer ip, id, ilen
      integer istp1, istp2, istp3, istp4, istp5, istp6
      character (len = 100):: fnm,fnm2

!     if(myid == 0)then

        istp1 = istep / 100000
        istp2 = mod(istep,100000) / 10000
        istp3 = mod(istep,10000) / 1000
        istp4 = mod(istep,1000) / 100
        istp5 = mod(istep,100) / 10
        istp6 = mod(istep,10)

        fnm = trim(dirpartout)//trim(fnm2)         &
            //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)     &
            //char(istp4 + 48)//char(istp5 + 48)//char(istp6 + 48)     &
            //'.dat'

        open(60, file = trim(fnm), status = 'unknown',                 &
                 form = 'formatted', position = 'append')

        do id = 1,npart
          write(60,600) id, ypglb(1,id), ypglb(2,id), ypglb(3,id),     &
                            wp(1,id), wp(2,id), wp(3,id),              &
                            fHIp(1,id), fHIp(2,id), fHIp(3,id)
        end do

        close(60)

!     end if

600   format(2x,i5,9(1pe16.6))
      end subroutine colldata 
!============================================================================

!===========================================================================
      subroutine checkredisbefore 
      use mpi
      use var_inc
      implicit none

      integer ip, id, ilen
      integer istp1, istp2, istp3, istp4, istp5, istp6
      character (len = 100):: fnm

      if(myid == 0)then

        istp1 = istep / 100000
        istp2 = mod(istep,100000) / 10000
        istp3 = mod(istep,10000) / 1000
        istp4 = mod(istep,1000) / 100
        istp5 = mod(istep,100) / 10
        istp6 = mod(istep,10)

        fnm = trim(dirpartout)//'checkRedisBefore2D'                                 &
            //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)     &
            //char(istp4 + 48)//char(istp5 + 48)//char(istp6 + 48)     &
            //'.dat'

        open(70, file = trim(fnm), status = 'unknown',                 &
                 form = 'formatted', position = 'append')

        do id = 1,npart
          write(70,700) id, ypglb(1,id), ypglb(2,id), ypglb(3,id),     &
                            wp(1,id), wp(2,id), wp(3,id),              &
                            fHIp(1,id), fHIp(2,id), fHIp(3,id)
        end do

        close(70)

      end if

700   format(2x,i5,9(1pe16.6))
      end subroutine checkredisbefore
!============================================================================

!===========================================================================
      subroutine checkredisafter
      use mpi
      use var_inc
      implicit none

      integer ip, id, ilen
      integer istp1, istp2, istp3, istp4, istp5, istp6
      character (len = 100):: fnm

      if(myid == 0)then

        istp1 = istep / 100000
        istp2 = mod(istep,100000) / 10000
        istp3 = mod(istep,10000) / 1000
        istp4 = mod(istep,1000) / 100
        istp5 = mod(istep,100) / 10
        istp6 = mod(istep,10)

        fnm = trim(dirpartout)//'checkRedisAfter2D'                                 &
            //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)     &
            //char(istp4 + 48)//char(istp5 + 48)//char(istp6 + 48)     &
            //'.dat'

        open(80, file = trim(fnm), status = 'unknown',                 &
                 form = 'formatted', position = 'append')

        do id = 1,npart
          write(80,800) id, ypglb(1,id), ypglb(2,id), ypglb(3,id),     &
!                           forcep(1,id), forcep(2,id), forcep(3,id),              &
                            wp(1,id), wp(2,id), wp(3,id),              &
                            fHIp(1,id), fHIp(2,id), fHIp(3,id)
        end do

        close(80)

      end if

800   format(2x,i5,9(1pe16.6))
      end subroutine checkredisafter
!============================================================================

      subroutine dataf(fnm2)
      use mpi
      use var_inc
      implicit none

      integer ip,i9,j9,k9,id, ilen
      integer istp1, istp2, istp3, istp4, istp5, istp6
      character (len = 100):: fnm,fnm2

      if(myid == 0)then

        istp1 = istep / 100000
        istp2 = mod(istep,100000) / 10000
        istp3 = mod(istep,10000) / 1000
        istp4 = mod(istep,1000) / 100
        istp5 = mod(istep,100) / 10
        istp6 = mod(istep,10)

        fnm = trim(dirpartout)//trim(fnm2)                           &
            //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)     &
            //char(istp4 + 48)//char(istp5 + 48)//char(istp6 + 48)     &
            //'.dat'

        open(60, file = trim(fnm), status = 'unknown',                 &
                 form = 'formatted', position = 'append')

       do k9 = 2,31
       do j9 = 2,31
       do i9 = 1,64
       write(60,600)i9,j9,k9,ibnodes(i9,j9,k9),f(1,i9,j9,k9),f(2,i9,j9,k9) &
                   ,f(3,i9,j9,k9),f(4,i9,j9,k9)
       end do
       end do
       end do

        close(60)

      end if

600   format(2x,4i5,4(1pe16.6))
      end subroutine dataf

!===========================================================================

