C
c       Moving particles in a 2D Couette flow
C       Multiple-relaxarion times
c       x - flow direction, y - cross channel
c       Domain x=0 to Lx, periodic BC in outer boundary
c             y=0 to Ly, no-slip in outer boundary
c       Details follow P. Lallemand and L.S. Luo, 2003,
C                      J. Comp. Phys. 184: 406 - 421.
c
c ================================
	program lbe2D
c ================================
	include'lbe.par'
        logical newflo
c
        data icx/0,1,0,-1,0,1,-1,-1,1/
        data icy/0,0,1,0,-1,1,1,-1,-1/
C Note in Fortran, the order is transm(0,0),transm(1,0),transm(2,0),...
        data transm/1.d0,-4.d0,4.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,
     .    1.d0,-1.d0,-2.d0, 1.d0,-2.d0, 0.d0, 0.d0, 1.d0, 0.d0,
     .    1.d0,-1.d0,-2.d0, 0.d0, 0.d0, 1.d0,-2.d0,-1.d0, 0.d0,
     .    1.d0,-1.d0,-2.d0,-1.d0, 2.d0, 0.d0, 0.d0, 1.d0, 0.d0,
     .    1.d0,-1.d0,-2.d0, 0.d0, 0.d0,-1.d0, 2.d0,-1.d0, 0.d0,
     .    1.d0, 2.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 0.d0, 1.d0,
     .    1.d0, 2.d0, 1.d0,-1.d0,-1.d0, 1.d0, 1.d0, 0.d0,-1.d0,
     .    1.d0, 2.d0, 1.d0,-1.d0,-1.d0,-1.d0,-1.d0, 0.d0, 1.d0,
     .    1.d0, 2.d0, 1.d0, 1.d0, 1.d0,-1.d0,-1.d0, 0.d0,-1.d0/
c The insverse is the transpose plus normalization below
        data transmi/1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,
     .   -4.d0,-1.d0,-1.d0,-1.d0,-1.d0, 2.d0, 2.d0, 2.d0, 2.d0,
     .    4.d0,-2.d0,-2.d0,-2.d0,-2.d0, 1.d0, 1.d0, 1.d0, 1.d0,
     .    0.d0, 1.d0, 0.d0,-1.d0, 0.d0, 1.d0,-1.d0,-1.d0, 1.d0,
     .    0.d0,-2.d0, 0.d0, 2.d0, 0.d0, 1.d0,-1.d0,-1.d0, 1.d0,
     .    0.d0, 0.d0, 1.d0, 0.d0,-1.d0, 1.d0, 1.d0,-1.d0,-1.d0,
     .    0.d0, 0.d0,-2.d0, 0.d0, 2.d0, 1.d0, 1.d0,-1.d0,-1.d0,
     .    0.d0, 1.d0,-1.d0, 1.d0,-1.d0, 0.d0, 0.d0, 0.d0, 0.d0,
     .    0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0,-1.d0, 1.d0,-1.d0/
C note the last two sss will be reset by omega or viscosity
c       data sss/0.d0,-1.64d0,-1.54d0,0.d0,-1.9d0,0.d0,-1.9d0,
c    .         0.d0,0.d0/
         data sss/0.d0,-1.5d0,-1.4d0,0.d0,-1.5d0,0.d0,-1.5d0,
     .         0.d0,0.d0/
C
          xcenter(1) = 100.d0
          ycenter(1) = 20.d0
           write(*,*)'xcenter,ycenter=',xcenter(1),ycenter(1)
c  The pairing link
        data ipopp/0,3,4,1,2,7,8,5,6/
C
c       write(*,*)'newflo =? .true. or .flase.'
       newflo = .true.
c        read(*,*)newflo
        pi = datan(1.0d0)*4.0d0
C
        if(newflo)then
c --- input parameters
c
C parameter matching Lallemand and Luo, JCP, 2003, 184: 406-421
       omega = 1.2d0
       visc = (1.d0/omega - 0.5d0)/3.0d0
c      uf = 1.5     
       ub1 = 0.d0
       ub2 = 0.d0
       u0 = 0.d0
       v0 = 0.d0
       umax = 0.05d0

       sss(7) = - omega
       sss(8) = - omega

c normalization
        do i=0,npop-1
          transmi(i,0)=transmi(i,0)/9.d0*sss(0)
          transmi(i,1)=transmi(i,1)/36.d0*sss(1)
          transmi(i,2)=transmi(i,2)/36.d0*sss(2)
          transmi(i,3)=transmi(i,3)/6.d0*sss(3)
          transmi(i,4)=transmi(i,4)/12.d0*sss(4)
          transmi(i,5)=transmi(i,5)/6.d0*sss(5)
          transmi(i,6)=transmi(i,6)/12.d0*sss(6)
          transmi(i,7)=transmi(i,7)/4.d0*sss(7)
          transmi(i,8)=transmi(i,8)/4.d0*sss(8)
        enddo
C Parameter for MRT
          ccc1=-1.d0
          alpha2=-2.d0
          alpha3=1.d0
          gamma1=1.0
          gamma2=3.d0
          gamma3=gamma1
          gamma4=-gamma2

	call input
        call beads_links
	
c --- initialisation

	call inithydro
	call equili
        call initpop

        else
        open ( unit=21, file='restart.flo.16000',
     .    form='unformatted', status='old' )
        read (21) f
        read (21) cs2,cs22,cssq,omega,den,visc,
     .    rho0,rt0,rt1,rt2,ub1,ub2,u0
        read(21)istep0
        read(21)fpois
        read(21)sss,transmi
        read(21)ccc1,alpha2,alpha3
        read(21)gamma1,gamma2,gamma3,gamma4
        read(21)ilink,alink
        read(21)nobst0,nobst1,ibnodes
        read(21)feqzero
        close(21)
C
        call hydrovar
C
        print*,' Number of time steps'
        read(5,*)nsteps

        print*,' Number of steps between printing profile'
        read(5,*)nout

        print*,' Number of steps between performing diagnostics'
        read(5,*)ndiag

        print*,' File for output: 5 chars'
        read(5,'(A)')fileout

        open(10,file=fileout//'.uy')
        open(11,file=fileout//'.vy')
        open(14,file=fileout//'.uvx')
        open(17,file=fileout//'.flux')
        open(18,file=fileout//'.fill')

        endif

c ------- MAIN LOOP
        kout = nout

        iconf=0
	do 10 istep = istep0,istep0 + nsteps

           time = istep

           if (mod(istep,ndiag).eq.0) then
              call diag0D
C maximum fluid velocity
          vmax = 0.d0
          do j=1,ny
           do i=1,nx
           v9 = dsqrt(u(i,j)**2 + v(i,j)**2)
           if(v9.gt.vmax) vmax = v9
           enddo
          enddo
C  total volume flux
            flux = 0.0d0
            do j=1,ny
C i=1 is at 0.5dx, not at x = 0
            flux = flux + (u(1,j)+u(nx,j))/2.0d0
            end do
           write(17,171)time,flux,vmax,xmom,ymom,torq,
     .        nobst0,nobst1,(nobst0+nobst1),ntlink
           write(18,181)time,xmom,ymom,torq,xmom_gain,ymom_gain,
     .         xmom_loss,ymom_loss,xmom-xmom_gain+xmom_loss,
     .        ymom-ymom_gain+ymom_loss,ifill,icover

171        format(2x,6(1pe16.5),4I10)
181        format(2x,10(1pe16.5),2I10)
           endif

           if (mod(istep,nout).eq.0) then
              call profil(istep,frce)
           endif

	   call collis_MRT
           call mixbc
c Then on circular beads
           call beads(istep)
           call move
	   call hydrovar
c
           if (mod(istep,kout).eq.0) then
              call config(istep,iconf)
              iconf=iconf+1
           endif
C
        do j = 1,ny
        do i = 1,nx
         ibnodes0(i,j)=ibnodes(i,j)
          do ip=1,npop-1
          ilink0(ip,i,j) = ilink(ip,i,j)
          end do
         end do
        enddo

        grav = 8*visc*umax/(dfloat(ny)*dfloat(ny)) 
        pmass = 1./6. * pi * (2*rad)**3 *2*rho0  
        du0 = 0.5d0 * (xmom+xmomp)/pmass + grav 
        xcenter(1) = xcenter(1) + 0.5d0 * (2.0*u0 + du0)
        u0 = u0 + du0

        if ( (ycenter(1).ge.(ny-rad)) .or. (ycenter(1).lt.(rad)) ) then
           v0 = -v0
        endif

        dv0 = 0.5d0 * (ymom+ymomp)/pmass 
        ycenter(1) = ycenter(1) + 0.5d0 * (2.0*v0 + dv0) 
        v0 = v0 + dv0

        xmomp = xmom
        ymomp = ymom

        angvel = torq / (2.0/5.0 * pmass * rad**2)

        write (23,*) istep, xcenter(1), ycenter(1), u0, v0, angvel

        call beads_links
        call filling

c-------- end of main loop

10	continue

C Save information for continued run
      open ( unit=22, file='endrun.flo', form='unformatted',
     &       status='unknown' )
        write (22) f
        write (22) cs2,cs22,cssq,omega,den,visc,
     .    rho0,rt0,rt1,rt2,ub1,ub2,u0
        write(22)(istep0+nsteps)
        write(22)fpois
        write(22)sss,transmi
        write(22)ccc1,alpha2,alpha3
        write(22)gamma1,gamma2,gamma3,gamma4
        write(22)ilink,alink
        write(22)nobst0,nobst1,ibnodes
        write(22)feqzero
        close(22)

        Re_ave = flux/visc
        write(*,*)'omega,visc,uf, flux, Re_ave, vmax = ',
     >     omega,visc,uf,flux,Re_ave,vmax

        stop
	end
C -------------------------------------------------
        subroutine filling
        include'lbe.par'
      real cc(8),cc2
      parameter(cc2=1.41421356)
      data cc /cc2,cc2,cc2,cc2, 1., 1., 1., 1./
      xmom_gain = 0.0
      ymom_gain = 0.0
      ifill = 0
      icover = 0

      xmom_loss = 0.0
      ymom_loss = 0.0

        do j = 1, ny
          do i = 1, nx
          IF(ibnodes0(i,j).gt.0 .and. ibnodes(i,j).lt.0)then
          ifill = ifill + 1
C use extrapolation to fill new fluid nodes
C
c first identify the direction most close to the local normal
        prod = -100.
        do ip=1,npop-1
        alpha = alink(ip,i,j)
        ix = icx(ip)
        iy = icy(ip)
        xx0=dfloat(i)+ix*alpha-0.5d0-xcenter(1)
        yy0=dfloat(j)+iy*alpha-0.5d0-ycenter(1)
        prod1 = cc(ip)*(xx0*ix+yy0*iy)
        if(prod1.gt.prod)then
        ip9 = ip
        prod = prod1
        endif
        enddo
C filling
c       if(ip9.le.0) then
c       write(*,*)'ip9=',ip9
c       stop
c       endif 
        ix = icx(ip9)
        iy = icy(ip9)
C ip = 0  previously I did not handle this part correctly ....
        ip = 0
        f(ip,i,j) = 3.d0*f(ip,i+ix,j+iy)
     .     -3.d0*f(ip,i+2*ix,j+2*iy) + f(ip,i+3*ix,j+3*iy)

        do ip=1,npop-1

        ix1 = icx(ip)
        iy1 = icy(ip)
C  Only apply extropolation for unknown f
        if(ibnodes0(i-ix1,j-iy1).gt.0)then
        f(ip,i,j) = 3.d0*f(ip,i+ix,j+iy)
     .     -3.d0*f(ip,i+2*ix,j+2*iy) + f(ip,i+3*ix,j+3*iy)
        endif

C compute added momentum to the fluid domain due to filling
C All links are included here
        xmom_gain = xmom_gain + f(ip,i,j)*float(ix1)
        ymom_gain = ymom_gain + f(ip,i,j)*float(iy1)
        enddo

          ENDIF
          end do
        end do

C Covered fluid node: compute the loss of momentum
        do j = 1, ny
          do i = 1, nx
          IF(ibnodes0(i,j).lt.0 .and. ibnodes(i,j).gt.0)then
          icover = icover + 1
        do ip=1,npop-1
        ix1 = icx(ip)
        iy1 = icy(ip)
C compute lost momentum to the fluid domain due to covering by solid particle
C All links are included here
        xmom_loss = xmom_loss + f(ip,i,j)*float(ix1)
        ymom_loss = ymom_loss + f(ip,i,j)*float(iy1)
        enddo

          ENDIF
          end do
        end do

        return
        end
c--------------------------------------------------
        subroutine beads_links
        include'lbe.par'
c 
        do j = 1, ny
          do i = 1, nx
            do ip=1,npop-1
              ilink(ip,i,j) = -1          
            end do
          end do
        end do
c
C  bead interior nodes not on wall
         nobst0 = 0
C  bead nodes at interface
         nobst1 = 0

        do j = 1,ny
        do i = 1,nx
         ibnodes(i,j)=-1
         end do
        enddo
C
        DO IPART = 1,MPART
        do j = 1,ny
        yy0 = dfloat(j) - 0.5d0 - ycenter(ipart)
        do i = 1,nx
        xx0 = dfloat(i) - 0.5d0 - xcenter(ipart)

        rr0 = dsqrt(xx0*xx0+yy0*yy0)

        IF(rr0.lt.rad)nobst0 = nobst0 + 1 
        IF(rr0.eq.rad)nobst1 = nobst1 + 1 
         if(rr0.le.rad)then
         ibnodes(i,j)=1
         endif
         end do
        enddo

       END DO
C-----------------
c
        ntlink = 0

        DO IPART = 1,MPART
        do j = 1,ny
        yy0 = dfloat(j) - 0.5d0 - ycenter(ipart)
        do i = 1,nx
        xx0 = dfloat(i) - 0.5d0 - xcenter(ipart)

        rr0 = dsqrt(xx0*xx0+yy0*yy0)

        IF(rr0.gt.rad)then

         do ip=1,npop-1
         imove = i + icx(ip)
         jmove = j + icy(ip)

         xx1 = dfloat(imove) - 0.5d0 - xcenter(ipart)
         yy1 = dfloat(jmove) - 0.5d0 - ycenter(ipart)

         rr1 = dsqrt(yy1**2+xx1**2)

         if(rr1.le.rad)then

C note (i,j) is the flow-domain point, and (imove,jmove) the wall node
C alpha is the percentage of link from wall measured from (i,j)
C  alpha = percentage of distance from the flow-domain node to actual wall 
C
C   [x1-alpha(x0-x1)]^2 + [y1-alpha(y0-y1)]^2 = rad^2
C    then alpha = 1 - alpha
C
         ilink(ip,i,j) = 1
         ntlink = ntlink + 1

         aa = rr1*rr1 + rr0*rr0 - 2.0d0*(yy1*yy0+xx1*xx0)
         bb = 2.d0*(yy1*(yy0-yy1)+xx1*(xx0-xx1))
         cc = rr1*rr1 - rad*rad
         alpha = dsqrt(bb*bb-4.d0*aa*cc)
         alpha = (-bb + alpha)/(2.d0*aa)
C
           alpha = 1.d0 - alpha

           alink(ip,i,j) = alpha
            endif
           end do
        ENDIF
         end do
        enddo

       END DO

        return
        end

c---------------------------------------------------
	subroutine input

	include'lbe.par'
c---------------------------------------------------
        istep0 = 0

	print*,' Number of steps'
        read(5,*)nsteps

	print*,' Number of steps between printing profile'
        read(5,*)nout

	print*,' Number of steps between performing diagnostics'
	read(5,*)ndiag

c print*,' Initial density for the Poiseuille force'
c read(5,*)rho0
        rho0=1.0d0

	cs2  = 1.0d0 / 3.0d0
	cs22 = 2.0d0 * cs2
	cssq = 2.0d0 / 9.0d0

        write(*,*)'omega,visc,uf = ',omega,visc,uf

        print*,' File for output: 5 chars'
        read(5,'(A)')fileout

        open(10,file=fileout//'.uy')
        open(11,file=fileout//'.vy')
        open(14,file=fileout//'.uvx')
        open(17,file=fileout//'.flux')
        open(18,file=fileout//'.fill')

	print*,'*****************************************'
	print*,' Lattice BGK model, 2D with 9 velocities'
	print*,'*****************************************'
	print*,'Number of cells :',nx,'*',ny
	print*,'Nsteps :',nsteps
	print*,'Relaxation frequency :',omega
c       print*,'Initial velocity for this Poiseuille force :',u0
	print*,'Initial density :',rho0
	write(6,'(A)')'Output file :',fileout
c       pause

c reduced density

	den = rho0/float(npop) 

	print*,' Viscosity :',visc

c calculation of the constant applied force

	fpois = 8.0d0 * visc * uf / dfloat(ny) / dfloat(ny)
        fpois = rho0*fpois/6.  ! # of biased populations
	print*,' Intensity of the applied force ',fpois
	
	return
	end
c--------------------------------------------------
	subroutine inithydro
	
	include'lbe.par'
c---------------------------------------------------
	do j = 0, ny+1
	  do i = 0, nx+1
             rho(i,j)  = 0.0d0
             u(i,j) = u0
             v(i,j) = 0.0d0
	  enddo
        enddo

        do j = 1, ny
           u(nx,j)=umax*(1 - (j-0.5d0-ny/2.0)**2/(ny/2.0)**2) 
           rho(nx,j) = 2 * ( f(8,nx,j) + f(1,nx,j) + f(5,nx,j) )
     .        + f(0,nx,j) + f(2,nx,j) + f(4,nx,j) - rho0*u(nx,j)
        enddo

	rt0 = rho0* 4.0d0 / 9.0d0
	rt1 = rho0/ 9.0d0
	rt2 = rho0/ 36.0d0

        wwp(1) = rho0*2.d0/9.d0
        wwp(2) = rho0*2.d0/9.d0
        wwp(3) = rho0*2.d0/9.d0
        wwp(4) = rho0*2.d0/9.d0
        wwp(5) = rho0*2.d0/36.d0
        wwp(6) = rho0*2.d0/36.d0
        wwp(7) = rho0*2.d0/36.d0
        wwp(8) = rho0*2.d0/36.d0

	return	
	end
c --------------------------------------------------
	subroutine initpop
	
	include'lbe.par'
c---------------------------------------------------
	do j = 0, ny+1
	  do i = 0, nx+1
            do ip=0,npop-1
        if(ibnodes(i,j).lt.0)then
            f(ip,i,j)=feq(ip,i,j)
        endif
            end do
          end do
        end do

        return
        end
c----------------------------------------------
	subroutine move
c----------------------------------------------
	include'lbe.par'
c---------------------------------------------

        do j = ny,1,-1
           do i = 1,nx
              f(2,i,j) = f(2,i,j-1)
              f(6,i,j) = f(6,i+1,j-1)
           enddo
        enddo

        do j = ny,1,-1
           do i = nx,1,-1
              f(1,i,j) = f(1,i-1,j)
              f(5,i,j) = f(5,i-1,j-1)
           enddo
        enddo

        do j = 1,ny
           do i = nx,1,-1
           f(4,i,j) = f(4,i,j+1)
           f(8,i,j) = f(8,i-1,j+1)
           enddo
        enddo

        do j = 1,ny
           do i = 1,nx
           f(3,i,j) = f(3,i+1,j)
           f(7,i,j) = f(7,i+1,j+1)
           enddo
        enddo

	return
	end
c---------------------------------------------
	subroutine hydrovar

	include'lbe.par'
c--------------------------------------------
c Calculation of density / pressure and velocities

	do j = 1,ny
	 do i = 1,nx
        if(ibnodes(i,j).lt.0)then
          rho(i,j) = f(0,i,j)+f(1,i,j)+f(2,i,j)
     .     + f(3,i,j)+ f(4,i,j)+f(5,i,j)
     .     + f(6,i,j)+f(7,i,j)+f(8,i,j)

          rhoi=1.d0/rho0

	  u(i,j)=(f(1,i,j)-f(3,i,j)+f(5,i,j)- 
     .            f(6,i,j)-f(7,i,j)+f(8,i,j))*rhoi 

	  v(i,j)=(f(5,i,j)+f(2,i,j)+f(6,i,j)
     .          - f(7,i,j)-f(4,i,j)-f(8,i,j))*rhoi

        endif
	  enddo
	enddo

	return
	end
c-------------------------------------------------
	subroutine equili

	include'lbe.par'
c-------------------------------------------------
	do j = 1,ny
	  do i = 1,nx
c
        if(ibnodes(i,j).lt.0)then
            rl=rho(i,j)/rho0
c
	    usq = u(i,j) * u(i,j) 
	    vsq = v(i,j) * v(i,j)
	    sumsq = (usq + vsq) / cs22
	    sumsq2 = sumsq * (1.0d0 - cs2) / cs2
	    u2 = usq / cssq 
            v2 = vsq / cssq
	    ui = u(i,j) / cs2
	    vi = v(i,j) / cs2
	    uv = ui * vi
c   Here I took the suggestion of He and Luo, 1997, lattice Boltzmann
C   model for incompressible N-S eqn. J. Statistical Phys.
C   88: 927-944.
	    feq(0,i,j) = rt0*(rl - sumsq)
	    feq(1,i,j) = rt1*(rl - sumsq + u2 + ui)
	    feq(2,i,j) = rt1*(rl - sumsq + v2 + vi)
	    feq(3,i,j) = rt1*(rl - sumsq + u2 - ui)
	    feq(4,i,j) = rt1*(rl - sumsq + v2 - vi)
	    feq(5,i,j) = rt2*(rl + sumsq2 + ui + vi + uv)
	    feq(6,i,j) = rt2*(rl + sumsq2 - ui + vi - uv)
	    feq(7,i,j) = rt2*(rl + sumsq2 - ui - vi + uv)
	    feq(8,i,j) = rt2*(rl + sumsq2 + ui - vi - uv)
        endif
	
	   enddo
	enddo

c       write(*,*)(ip,feq(ip,10,10),ip=0,npop-1)
		
	return
	end
c----------------------------------------------------------
	subroutine collis_MRT

	include'lbe.par'
        dimension f9(0:npop-1),fchange(0:npop-1)
c----------------------------------------------------------
	 do j = 1,ny
	  do i = 1,nx
        if(ibnodes(i,j).lt.0)then

C Hui corrected the problem for me
        aden = 0.d0
        do ip=0,npop-1
        f9(ip)=f(ip,i,j)
        aden = aden + f(ip,i,j)
        enddo
C Part 1: Compute moments
C Improved computational efficiency, although the code
C looks longer
CC       aden = rho(i,j)
        ajx = rho0*u(i,j)
        ajy = rho0*v(i,j)

        p13 = f9(1)+f9(3)
        q13 = f9(1)-f9(3)
        p24 = f9(2)+f9(4)
        q24 = f9(2)-f9(4)
        p56 = f9(5)+f9(6)
        q56 = f9(5)-f9(6)
        p78 = f9(7)+f9(8)
        q78 = f9(7)-f9(8)

         p1to4 = p13+p24
         p5to8 = p56+p78
         rmom(1) = -4.d0*f9(0)-p1to4+2.d0*p5to8
         rmom(2) = 4.d0*f9(0)-2.d0*p1to4+p5to8
         rmom(4) = -2.d0*q13+q56-q78
         comb1 = p56-p78
         rmom(6) = -2.d0*q24+comb1
         rmom(7) = p13-p24
         rmom(8) = q56+q78

C Part 2: Compute (rmom - rmome)
        ajx2 = ajx*ajx
        ajy2 = ajy*ajy
        ajxy = ajx*ajy
        ajx2y2 = ajx2 + ajy2
        
        drmom1 = rmom(1) - alpha2*aden - gamma2*ajx2y2
        drmom2 = rmom(2) - alpha3*aden - gamma4*ajx2y2
        drmom4 = rmom(4) - ccc1*ajx
        drmom6 = rmom(6) - ccc1*ajy
        drmom7 = rmom(7) - gamma1*(ajx2-ajy2)
        drmom8 = rmom(8) - gamma3*ajxy
C Part 3: Inverse
        fchange(0) = transmi(0,1)*drmom1+transmi(0,2)*drmom2

        fchange(1) = transmi(1,1)*drmom1+transmi(1,2)*drmom2
     .    + transmi(1,4)*drmom4+transmi(1,7)*drmom7
        fchange(3) = transmi(3,1)*drmom1+transmi(3,2)*drmom2
     .    + transmi(3,4)*drmom4+transmi(3,7)*drmom7

        fchange(2) = transmi(2,1)*drmom1+transmi(2,2)*drmom2
     .    + transmi(2,6)*drmom6+transmi(2,7)*drmom7
        fchange(4) = transmi(4,1)*drmom1+transmi(4,2)*drmom2
     .    + transmi(4,6)*drmom6+transmi(4,7)*drmom7

        fchange(5) = transmi(5,1)*drmom1+transmi(5,2)*drmom2
     .    + transmi(5,4)*drmom4+transmi(5,6)*drmom6
     .    + transmi(5,8)*drmom8
        fchange(6) = transmi(6,1)*drmom1+transmi(6,2)*drmom2
     .    + transmi(6,4)*drmom4+transmi(6,6)*drmom6
     .    + transmi(6,8)*drmom8
        fchange(7) = transmi(7,1)*drmom1+transmi(7,2)*drmom2
     .    + transmi(7,4)*drmom4+transmi(7,6)*drmom6
     .    + transmi(7,8)*drmom8
        fchange(8) = transmi(8,1)*drmom1+transmi(8,2)*drmom2
     .    + transmi(8,4)*drmom4+transmi(8,6)*drmom6
     .    + transmi(8,8)*drmom8

	do k = 0,npop-1
	   f(k,i,j)=f9(k) + fchange(k)
	  enddo

        endif
	 enddo
	enddo
	    
	return 
	end  
c-------------------------------------------------------------
	subroutine mixbc
	
	include'lbe.par'
c-------------------------------------------------------------

!        call equili

!    Inlet

	do j = 1,ny

           u(1,j) = 2./3. * umax
! (f(1,1,j)-f(3,1,j)+f(5,1,j)-
!     .               f(6,1,j)-f(7,1,j)+f(8,1,j))/rho0
           rho(0,j) = rho0*u(1,j) + 2*(f(3,0,j)+f(6,0,j)+f(7,0,j))
     .     + f(0,0,j) + f(2,0,j) + f(4,0,j) 
           f(1,0,j) = feq(1,0,j) + f(3,0,j) - feq(3,0,j)
           f(5,0,j) = 0.5d0 * ( rho0*u(1,j) + feq(3,0,j) - feq(1,0,j)
     .                 - f(2,0,j) + f(4,0,j) ) + f(7,0,j)
           f(8,0,j) = 0.5d0 * ( rho0*u(1,j) + feq(3,0,j) - feq(1,0,j)
     .     + f(2,0,j) - f(4,0,j) ) + f(6,0,j)

        enddo

!     Outlet 

	do j = 1,ny

           u(nx,j) = umax*(1-(j-0.5d0-ny/2.0)**2/(ny/2.0)**2) 
           rho(nx+1,j) = 2*(f(8,nx+1,j) + f(1,nx+1,j) + f(5,nx+1,j))
     .       + f(0,nx+1,j) + f(2,nx+1,j) + f(4,nx+1,j) - rho0*u(nx,j) 

           f(3,nx+1,j)= f(1,nx+1,j)-feq(1,nx+1,j)+feq(3,nx+1,j)

           f(6,nx+1,j)=0.5d0*(f(4,nx+1,j)-f(2,nx+1,j)+feq(1,nx+1,j)
     .          - feq(3,nx+1,j) - rho0*u(nx,j)) + f(8,nx+1,j)

           f(7,nx+1,j)=0.5d0*(f(2,nx+1,j)-f(4,nx+1,j)+feq(1,nx+1,j)
     .              - feq(3,nx+1,j) - rho0*u(nx,j)) + f(5,nx+1,j)

        enddo

c NORTH case
clpw we use the buffer layer as temporary storage
c
	do i = 1,nx
        f(4,i,ny+1) = f(2,i,ny) 
        f(8,i-1,ny+1) = f(6,i,ny) + rho0/6.d0*ub1
        f(7,i+1,ny+1) = f(5,i,ny) - rho0/6.d0*ub1
	enddo

c SOUTH case

	do i = 1,nx
           f(2,i,0) = f(4,i,1)
           f(5,i-1,0) = f(7,i,1) + rho0/6.d0*ub2
           f(6,i+1,0) = f(8,i,1) - rho0/6.d0*ub2
        enddo

	return
	end
c ==========================
	subroutine beads(istep)
c ==========================
	include'lbe.par'
c--------------------------------------------------------
c
C  Also compute the force and torque acting on the particle
c
        xmom = 0.0d0
        ymom = 0.0d0
        torq = 0.0d0
        do j=1,ny
        do i=1,nx
        do ip=1,npop-1

        IF(ilink(ip,i,j).gt.0)THEN

        alpha = alink(ip,i,j)

        ipp = ipopp(ip)

        ix = icx(ip)
        iy = icy(ip)

        xx0=dfloat(i)+ix*alpha-0.5d0-xcenter(1)
        yy0=dfloat(j)+iy*alpha-0.5d0-ycenter(1)
        
        ff1 = f(ip,i,j)


        xdist= dfloat(i)-0.5d0-xcenter(1)
        ydist= dfloat(j)-0.5d0-ycenter(1)
        dist = dsqrt (xdist**2 + ydist**2) 
        xvel = ydist * angvel * rad/ dist
        yvel = - xdist * angvel * rad / dist

        if(alpha.le.0.5d0)then
C use the 3-point interpolation scheme of Lallemand and Luo (2003) JCP
        ff2 = f(ip,i-ix,j-iy)
        ff3 = f(ip,i-2*ix,j-2*iy)
        c1 = alpha*(1.d0+2.d0*alpha)
        c2 = 1.d0-4.d0*alpha*alpha
        c3 = - alpha*(1.d0-2.d0*alpha)
        f(ipp,i+ix,j+iy) = c1*ff1+c2*ff2+c3*ff3 - 3.d0*wwp(ip)
     .            *((u0+xvel)*dfloat(ix)+(v0+yvel)*dfloat(iy))
        else
        ff2 = f(ipp,i,j)
        ff3 = f(ipp,i-ix,j-iy)
        c1 = 1.d0/alpha/(2.d0*alpha+1.d0)
        c2 = (2.d0*alpha-1.d0)/alpha
        c3 = -(2.d0*alpha-1.d0)/(2.d0*alpha+1.d0)
        f(ipp,i+ix,j+iy) = c1*ff1+c2*ff2+c3*ff3 - 3.d0*wwp(ip)
     .            *((u0+xvel)*dfloat(ix)+(v0+yvel)*dfloat(iy))
     .            /alpha/(2.d0*alpha+1.d0)
        endif

        dxmom = (ff1 + f(ipp,i+ix,j+iy))*dfloat(ix)
        dymom = (ff1 + f(ipp,i+ix,j+iy))*dfloat(iy)

        xmom = xmom + dxmom
        ymom = ymom + dymom
        torq = torq + dymom*xx0 - dxmom*yy0

        ENDIF
        end do
        end do
        end do

c       write(90,*)'istep,xmom,ymom,torq=',istep,xmom,ymom,torq
c900     format(2x,4i4,10(1pe13.3))

	return
	end
c-----------------------------------------------------------
	subroutine profil(it,frce)
c-----------------------------------------------------------
	include'lbe.par'
c----------------------------------------------------------

        do j = 1,ny
          y9 = dfloat(j) - 0.5d0
          uu1 = u(nx-10,j)
          uu25 = u(10,j)
          vv1 = v(nx-10,j)
          vv25 = v(10,j)
          write(10,*) y9,uu1,uu25
          write(11,*) y9,vv1,vv25
        enddo

        write(10,'(bn)')
        write(11,'(bn)')
        write(14,'(bn)')

        if (it .gt. 1.0) then
        write(28,197) ((u(i,j),i=1,nx),j=1,ny)
        write(29,197) ((v(i,j),i=1,nx),j=1,ny)
        endif

197     format(2x,1000(1pe16.5))

	return
	end
c---------------------------------------------------------
	subroutine diag0D
c---------------------------------------------------------
        include'lbe.par'
c----------------------------------------------------------
	densit = 0.0d0

	   do j= 1,ny
	      do i = 1,nx
        if(ibnodes(i,j).lt.0)then
	do k = 0,npop-1
	         densit = densit + f(k,i,j)
	      enddo
        endif
	   enddo
	enddo

        numb = nx*ny - nobst0 - nobst1
	densit = densit/float(numb)

	umoy = 0.0d0
	vmoy = 0.0d0

	do j = 1,ny
	   do i = 1,nx
        if(ibnodes(i,j).lt.0)then
	      umoy = umoy + dabs(u(i,j))
	      vmoy = vmoy + dabs(v(i,j))
        endif
	   enddo
	enddo
	
	umoy = umoy / dfloat(numb)
	vmoy = vmoy / dfloat(numb)

        stime = time*visc/(float(ny)*float(ny))
        press = rho(2,1)/3.0

	print*,'diagnostic 0D : istep density umoy and vmoy',
     .          time,densit,umoy,vmoy

	return
	end
c =========================================================
	subroutine config(istep,iconf) 
c =========================================================
        include'lbe.par'
        dimension vort(nx,ny)
c -------------------------------------------
c Calculation of density / pressure and velocities
c
        iout=60+iconf
c
          do j=1,ny
          do i=1,nx
C Vorticity field ...........
          IF(ibnodes(i,j).lt.0)THEN
          axp = 0.5d0
          axm = 0.5d0
          byp = 0.5d0
          bym = 0.5d0

          vxp = 0.5d0*(v(i,j) + v(i+1,j))
          if(i.eq.nx)vxp = 0.5d0*(v(i,j) + v(1,j))
          vxm = 0.5d0*(v(i,j) + v(i-1,j))
          if(i.eq.1)vxm = 0.5d0*(v(i,j) + v(nx,j))
          uyp = 0.5d0*(u(i,j) + u(i,j+1))
          if(j.eq.ny) uyp = 0.0d0
          uym = 0.5d0*(u(i,j) + u(i,j-1))
          if(j.eq.1) uym = 0.0d0

          if(ilink(1,i,j).gt.0)then
            axp = alink(1,i,j)
            vxp = 0.d0
          endif
          if(ilink(3,i,j).gt.0)then
            axm = alink(3,i,j)
            vxm = 0.d0
          endif
          if(ilink(2,i,j).gt.0)then
            byp = alink(2,i,j)
            uyp = 0.d0
          endif
          if(ilink(4,i,j).gt.0)then
            bym = alink(4,i,j)
            uym = 0.d0
          endif

          vort(i,j)=(vxp-vxm)/(axp+axm)-(uyp-uym)/(byp+bym)

          ELSE
           vort(i,j) = 0.0d0
          ENDIF
          enddo
          enddo

        do j = 1,ny
	  do i = 1,nx
            write(iout,*)u(i,j)
            write(iout+10,*)v(i,j)
            write(iout+20,*)vort(i,j)
	   enddo
	enddo

         close(iout)
         close(iout+10)
         close(iout+20)

	return
	end


