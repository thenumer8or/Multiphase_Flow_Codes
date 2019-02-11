      program divide_density
      implicit none
      character (len = 100):: fnm
      character(len=80):: dirgenr, dircntdflow
      integer iprc1, iprc2, iprc3
      integer istp1, istp2, istp3, istp4, istp5, istp6
      integer myid, lx, ly, lz, llz, lly, nx, ny, nz, i, j
      integer nsteps, istat, imovie, istpload,nsplit,msplit
      double precision,allocatable,dimension(:,:,:,:):: f

      dirgenr = '/ptmp/canderse/64_LBM_part_1DFFT/'
      dircntdflow = trim(dirgenr)//'cntdflow/'

      nsplit = 8     !The number of procs in the z
      msplit = 4     !The number of procs in the y

      nx = 256  
      lx = nx

      ny = 256 
      nz = 256 

      ly = ny       
      lz = nz/nsplit   

      lly = ly/msplit

      istpload = 6

      istp1 = istpload / 100000
      istp2 = mod(istpload,100000) / 10000
      istp3 = mod(istpload,10000) / 1000
      istp4 = mod(istpload,1000) / 100
      istp5 = mod(istpload,100) / 10
      istp6 = mod(istpload,10)


      allocate (f(0:18,lx,ly,lz))

!*******Read in this***********************  

      DO i = 1,nsplit

      myid = (i-1)

      iprc1 = myid / 100
      iprc2 = mod(myid,100) / 10
      iprc3 = mod(myid,10)

      fnm = trim(dircntdflow)//'endrunflow1D8.'                     &
            //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)  &
            //char(istp4 + 48)//char(istp5 + 48)//char(istp6 + 48)  &
            //'.'                                                   &
            //char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)


       open(12, file = trim(fnm), status = 'unknown',               &
               form = 'unformatted')
       read(12) nsteps, istat, imovie
       read(12) f

       close(12)


       do  j = 1,msplit

        myid = ((i-1)*msplit) + (j-1)

        iprc1 = myid / 100
        iprc2 = mod(myid,100) / 10
        iprc3 = mod(myid,10)

!jxk means j procs in the y k procs in the z 
        fnm = trim(dircntdflow)//'endrunflow2D4x8.'                  &
              //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)   &
              //char(istp4 + 48)//char(istp5 + 48)//char(istp6 + 48)   &
              //'.'                                                    &
              //char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)

         open(12, file = trim(fnm), status = 'unknown',                &
                 form = 'unformatted')
         write(12) nsteps, istat, imovie
         write(12) f(:,:,(j-1)*lly+1:(j-1)*lly+lly,:)

         close(12)

        write(*,*) i,j,myid
 
        enddo


      ENDDO


      deallocate(f)
   
      end program divide_density
