C
C 2D LBE code for isothermal flow in a square cavity
C using structured nonuniform mesh and the least square procedure of 
C Shu C, Niu XD, Chew YT, 2002, Taylor-series expansion and 
C  least-squares-based lattice Boltzmann method: Two-dimensional 
C  formulation and its applications. Phys. Rev. E, 65, 036708.
C
C The lid is at y=L. The domain is 0<x<L and 0<y<L
C The lid velocity is ulip
C -------------------------------------------------
c
c ================================
	program lbe2D
c ================================
	implicit double precision(a-h,o-z)
	include'lbe.par'

        icx(0)=0
        icx(1)=1
        icx(2)=0
        icx(3)=-1
        icx(4)=0
        icx(5)=1
        icx(6)=-1
        icx(7)=-1
        icx(8)=1

        icy(0)=0
        icy(1)=0
        icy(2)=1
        icy(3)=0
        icy(4)=-1
        icy(5)=1
        icy(6)=1
        icy(7)=-1
        icy(8)=-1

c --- Matrix preparation
        call prepmatrix
c       stop

c --- input parameters
	call input

c --- initialisation

	call inithydro
	call equili
        call initpop

        Umax0 = ulip

c ------- MAIN LOOP
        iconf=0
	do 10 istep = 0,nsteps

           time = istep
           if (mod(istep,ndiag).eq.0) then
              call diag0D(istep)
           endif

           if (mod(istep,nout).eq.0) then
              call profil(istep)
           endif

           if (mod(istep,nout).eq.0) then
              call config(istep,iconf)
              iconf=iconf+1
           endif

C We do collision first to compute all g.
C
	   call collis
C
c now we do the matrix multiplications
c
           call transform
C
c  The following takes care of BCs.
c
           call mbc

	   call hydrovar
	   call equili

c   check convergence ...
       Umax1 = 0.d0
        do j = 1,ny
           do i = 1,nx
c          sss = dabs(u(i,j)**2 + v(i,j)**2)
           sss = - u(i,j)
           if(sss.gt.Umax1)Umax1=sss
           enddo
        enddo
        relerr =dabs(Umax1-Umax0)/Umax0
        write(19,*)time,relerr,Umax0,Umax1
        Umax0 = Umax1

c-------- end of main loop

10	continue
	end
C-----------------------------------------------
        subroutine prepmatrix

        implicit double precision(a-h,o-z)
        include'lbe.par'

       dimension ssij(9,6),ssij0(9,6),
     .    sinv(7,7),ssc(7,7),ssa(9)
       integer indx(7)
        n = 6
        np = n+1
        n2 = 9

c define the nonuniform mesh
c in lattice units
       xx(1)=0.5d0
       xx(2)=1.5d0
c
C Can do non-uniform mesh if needed
c      iss = (nx-1)/2-1
c      beta = 3.35d0**(1.d0/dfloat(iss) )

c      ddx = beta
c      do i=3,(nx+1)/2
c      xx(i) = xx(i-1) + ddx
c      ddx = ddx*beta
c      enddo

c      WL = 2.d0*xx((nx+1)/2)

c      do i=1,(nx-1)/2
c      xx(nx+1-i) = WL - xx(i)
c      enddo

C For comparison with the finite volume code, use uniform mesh
       do i=3,nx
       xx(i) = xx(i-1) + 1.0
       enddo
       WL = xx(nx) + 0.5d0

       if(ny.ne.nx)stop
       do j=1,ny
       yy(j)=xx(j)
       enddo

       do i=2,nx
       dx(i)=xx(i)-xx(i-1)
       dy(i)=dx(i)
       enddo

       dx(1)=1.0d0
       dy(1)=1.0d0
       dx(nx+1)=1.0d0
       dy(ny+1)=1.0d0

       do i=1,nx
       write(21,*)i,xx(i),yy(i)
       enddo
       do i=1,nx+1
       write(22,*)i,dx(i),dy(i)
       enddo

C  maximum ratio
       rmax = dx((nx+1)/2)/dx(2)
       write(*,*)'rmax=',rmax

C Prepare the matrix for the least-square procedure
C Part 1: On near-boundary lines, using 3x3 Lagrangian
C interpolation
       jy = 1
       ddy1 = dy(jy)
       ddy2 = dy(jy+1)
       do 42 ix=1,nx
       ddx1 = dx(ix)
       ddx2 = dx(ix+1)
       do 41 id=1,8
       y0 = -icy(id)
       ym = -icy(id)+ddy1
       yp = -icy(id)-ddy2
       x0 = -icx(id)
       xm = -icx(id)+ddx1
       xp = -icx(id)-ddx2
       ssa(1) = xm*xp*ym*yp/(ddx1*ddx2*ddy1*ddy2)
       ssa(2) = x0*xp*y0*yp/((ddx1+ddx2)*ddx1*(ddy1+ddy2)*ddy1)
       ssa(3) = - x0*xp*ym*yp/((ddx1+ddx2)*ddx1*ddy1*ddy2)
       ssa(4) = x0*xp*y0*ym/((ddx1+ddx2)*ddx1*(ddy1+ddy2)*ddy2)
       ssa(5) = - xm*xp*y0*yp/(ddx2*ddx1*(ddy1+ddy2)*ddy1)
       ssa(6) = - xm*xp*y0*ym/(ddx2*ddx1*(ddy1+ddy2)*ddy2)
       ssa(7) = x0*xm*y0*yp/((ddx1+ddx2)*ddx2*(ddy1+ddy2)*ddy1)
       ssa(8) = - x0*xm*ym*yp/((ddx1+ddx2)*ddx2*ddy1*ddy2)
       ssa(9) = x0*xm*y0*ym/((ddx1+ddx2)*ddx2*(ddy1+ddy2)*ddy2)
      do j=1,n2
       fam(id,ix,jy,j)=ssa(j)
       enddo
41     enddo
42      ENDDO

       jy = ny
       ddy1 = dy(jy)
       ddy2 = dy(jy+1)
       do 52 ix=1,nx
       ddx1 = dx(ix)
       ddx2 = dx(ix+1)
       do 51 id=1,8
       y0 = -icy(id)
       ym = -icy(id)+ddy1
       yp = -icy(id)-ddy2
       x0 = -icx(id)
       xm = -icx(id)+ddx1
       xp = -icx(id)-ddx2
       ssa(1) = xm*xp*ym*yp/(ddx1*ddx2*ddy1*ddy2)
       ssa(2) = x0*xp*y0*yp/((ddx1+ddx2)*ddx1*(ddy1+ddy2)*ddy1)
       ssa(3) = - x0*xp*ym*yp/((ddx1+ddx2)*ddx1*ddy1*ddy2)
       ssa(4) = x0*xp*y0*ym/((ddx1+ddx2)*ddx1*(ddy1+ddy2)*ddy2)
       ssa(5) = - xm*xp*y0*yp/(ddx2*ddx1*(ddy1+ddy2)*ddy1)
       ssa(6) = - xm*xp*y0*ym/(ddx2*ddx1*(ddy1+ddy2)*ddy2)
       ssa(7) = x0*xm*y0*yp/((ddx1+ddx2)*ddx2*(ddy1+ddy2)*ddy1)
       ssa(8) = - x0*xm*ym*yp/((ddx1+ddx2)*ddx2*ddy1*ddy2)
       ssa(9) = x0*xm*y0*ym/((ddx1+ddx2)*ddx2*(ddy1+ddy2)*ddy2)
      do j=1,n2
       fam(id,ix,jy,j)=ssa(j)
       enddo
51     enddo
52      ENDDO

       do 63 jy=2,ny-1
       ddy1 = dy(jy)
       ddy2 = dy(jy+1)
c      do 62 ix=1,nx
       ix = 1
       ddx1 = dx(ix)
       ddx2 = dx(ix+1)
       do 61 id=1,8
       y0 = -icy(id)
       ym = -icy(id)+ddy1
       yp = -icy(id)-ddy2
       x0 = -icx(id)
       xm = -icx(id)+ddx1
       xp = -icx(id)-ddx2
       ssa(1) = xm*xp*ym*yp/(ddx1*ddx2*ddy1*ddy2)
       ssa(2) = x0*xp*y0*yp/((ddx1+ddx2)*ddx1*(ddy1+ddy2)*ddy1)
       ssa(3) = - x0*xp*ym*yp/((ddx1+ddx2)*ddx1*ddy1*ddy2)
       ssa(4) = x0*xp*y0*ym/((ddx1+ddx2)*ddx1*(ddy1+ddy2)*ddy2)
       ssa(5) = - xm*xp*y0*yp/(ddx2*ddx1*(ddy1+ddy2)*ddy1)
       ssa(6) = - xm*xp*y0*ym/(ddx2*ddx1*(ddy1+ddy2)*ddy2)
       ssa(7) = x0*xm*y0*yp/((ddx1+ddx2)*ddx2*(ddy1+ddy2)*ddy1)
       ssa(8) = - x0*xm*ym*yp/((ddx1+ddx2)*ddx2*ddy1*ddy2)
       ssa(9) = x0*xm*y0*ym/((ddx1+ddx2)*ddx2*(ddy1+ddy2)*ddy2)
      do j=1,n2
       fam(id,ix,jy,j)=ssa(j)
       enddo
61     enddo
c62      ENDDO
63      ENDDO

       do 73 jy=2,ny-1
       ddy1 = dy(jy)
       ddy2 = dy(jy+1)
c      do 72 ix=1,nx
       ix = nx
       ddx1 = dx(ix)
       ddx2 = dx(ix+1)
       do 71 id=1,8
       y0 = -icy(id)
       ym = -icy(id)+ddy1
       yp = -icy(id)-ddy2
       x0 = -icx(id)
       xm = -icx(id)+ddx1
       xp = -icx(id)-ddx2
       ssa(1) = xm*xp*ym*yp/(ddx1*ddx2*ddy1*ddy2)
       ssa(2) = x0*xp*y0*yp/((ddx1+ddx2)*ddx1*(ddy1+ddy2)*ddy1)
       ssa(3) = - x0*xp*ym*yp/((ddx1+ddx2)*ddx1*ddy1*ddy2)
       ssa(4) = x0*xp*y0*ym/((ddx1+ddx2)*ddx1*(ddy1+ddy2)*ddy2)
       ssa(5) = - xm*xp*y0*yp/(ddx2*ddx1*(ddy1+ddy2)*ddy1)
       ssa(6) = - xm*xp*y0*ym/(ddx2*ddx1*(ddy1+ddy2)*ddy2)
       ssa(7) = x0*xm*y0*yp/((ddx1+ddx2)*ddx2*(ddy1+ddy2)*ddy1)
       ssa(8) = - x0*xm*ym*yp/((ddx1+ddx2)*ddx2*ddy1*ddy2)
       ssa(9) = x0*xm*y0*ym/((ddx1+ddx2)*ddx2*(ddy1+ddy2)*ddy2)
      do j=1,n2
       fam(id,ix,jy,j)=ssa(j)
       enddo
71     enddo
c72      ENDDO
73      ENDDO

C Prepare the matrix for the least-square procedure
C Part 2: using the least-square procedure
C for interier nodes
       do 33 jy=2,ny-1
c      jy = 45
       ddy1 = dy(jy)
       ddy2 = dy(jy+1)

       do 32 ix=2,nx-1
c      ix = 1
       ddx1 = dx(ix)
       ddx2 = dx(ix+1)

C i: P,A,B,C,D,E,F,G,H
C       C4      E6     H9
C       
C       B3      P1     G8
C
C       A2      D5     F7
C
C id is the velocity direction
       do 31 id=1,8
       do i=1,9
       ssij(i,1) = 1.0d0
       enddo
       ssij(1,2) = icx(id)
       ssij(2,2) = -ddx1 + icx(id)
       ssij(3,2) = -ddx1 + icx(id)
       ssij(4,2) = -ddx1 + icx(id)
       ssij(5,2) = icx(id)
       ssij(6,2) = icx(id)
       ssij(7,2) = ddx2 + icx(id)
       ssij(8,2) = ddx2 + icx(id)
       ssij(9,2) = ddx2 + icx(id)

       ssij(1,3) = icy(id)
       ssij(2,3) = -ddy1 + icy(id)
       ssij(3,3) = icy(id)
       ssij(4,3) = ddy2 + icy(id)
       ssij(5,3) = -ddy1 + icy(id)
       ssij(6,3) = ddy2 + icy(id)
       ssij(7,3) = -ddy1 + icy(id)
       ssij(8,3) = icy(id)
       ssij(9,3) = ddy2 + icy(id)

       do i=1,n2
       ssij(i,4) = ssij(i,2)**2/2.d0
       ssij(i,5) = ssij(i,3)**2/2.d0
       ssij(i,6) = ssij(i,2)*ssij(i,3)
       enddo

       do j=1,n
       do i=1,n2
       ssij0(i,j)=ssij(i,j)
       enddo
       enddo 

C Compute [S^T S]

       do j=1,n
       do i=1,n
       ssc(i,j)=0.d0
       do k=1,n2
       ssc(i,j)=ssc(i,j) + ssij(k,i)*ssij(k,j)
       enddo

       enddo
       enddo
C Compute the inverse
       do 12 i=1,n
       do 11 j=1,n
       sinv(i,j)=0.d0
11     enddo
       sinv(i,i)=1.d0
12     enddo
       call ludcmp(ssc,n,np,indx,d)
       do 13 j=1,n
       call lubksb(ssc,n,np,indx,sinv(1,j))
13     enddo

C compute  [S^T S]^(-1) X S^T, but only the first row
      do j=1,n2
       ssa(j)=0.d0
       do k=1,n
       ssa(j)=ssa(j) + sinv(1,k)*ssij0(j,k)
       enddo
       fam(id,ix,jy,j)=ssa(j)
       enddo

31     enddo

32     ENDDO
33     ENDDO

       return
       end
c--------------------------------------------------
	subroutine input

	implicit double precision(a-h,o-z)
	include'lbe.par'
c---------------------------------------------------
	print*,' Number of steps'
        read(5,*)nsteps

	print*,' Number of steps between printing profile'
        read(5,*)nout

	print*,' Number of steps between performing diagnostics'
	read(5,*)ndiag

	print*,' Initial density'
	read(5,*)rho0

	print*,'Lip moving speed'
	read(5,*) ulip

        print*,' Reynolds number'
        read(5,*)Reynolds

        visc = ulip*WL/Reynolds

c constants
        cs2  = 1.0d0 / 3.0d0
        cs22 = 2.0d0 * cs2
        cssq = 2.0d0 / 9.0d0

        omega = 1.d0/(visc/cs2 + 0.5d0)

c       print*,' Relaxation frequency omega'
c       read(5,*)omega

        print*,' File for output: 5 chars'
        read(5,'(A)')fileout

        open(10,file='velprof1.dat')
        open(11,file='velprof2.dat')
        open(12,file=fileout//'.uc')
        open(16,file=fileout//'.pop')
        open(17,file=fileout//'.probe')

	print*,'*****************************************'
	print*,' Lattice BGK model, 2D with 9 velocities'
	print*,'*****************************************'
	print*,'Number of cells :',nx,'*',ny
	print*,'Nsteps :',nsteps
	print*,'Initial density :',rho0
	if (iobst) then
	    print*,' Linear Obstacle with length :',nobst
	endif
	write(6,'(A)')'Output file :',fileout

c reduced density
	den = rho0/float(npop) 

	print*,'Relaxation frequency :',omega
	print*,' Viscosity :',visc
	print*,' driving speed ulip :', ulip
	print*,' Cavity Dimension:', WL
        print*,'Reynolds number:',Reynolds

	return
	end
c--------------------------------------------------
	subroutine inithydro
	
	implicit double precision(a-h,o-z)
	include'lbe.par'
c---------------------------------------------------
c       write(6,*) 'u0',u0
	do j = 0, ny+1
	  do i = 0, nx+1
           rho(i,j)  = rho0
c      u(i,j) = u0
             u(i,j) = 0.0
	      v(i,j) = 0.0d0
c         write(6,*) i,j,rho(i,j)
	  enddo
        enddo

	rt0 = rho0* 4.0d0 / 9.0d0
	rt1 = rho0/ 9.0d0
	rt2 = rho0/ 36.0d0

	return	
	end
c --------------------------------------------------
	subroutine initpop
	
	implicit double precision(a-h,o-z)
	include'lbe.par'
c---------------------------------------------------
	do j = 0, ny+1
	  do i = 0, nx+1
            do ip=0,npop-1
              f(ip,i,j)=feq(ip,i,j)
            end do
          end do
        end do

        return
        end
c----------------------------------------------
	subroutine transform
c----------------------------------------------
	implicit double precision(a-h,o-z)
	include'lbe.par'
c---------------------------------------------
	do j = 1,ny
	   do i = 1,nx

        do id=1,8
           f(id,i,j) = gg(id,i,j)*fam(id,i,j,1) 
     .    +  gg(id,i-1,j-1)*fam(id,i,j,2)
     .  +  gg(id,i-1,j)*fam(id,i,j,3) +  gg(id,i-1,j+1)*fam(id,i,j,4)
     .   +  gg(id,i,j-1)*fam(id,i,j,5) +  gg(id,i,j+1)*fam(id,i,j,6)
     .   +  gg(id,i+1,j-1)*fam(id,i,j,7) +  gg(id,i+1,j)*fam(id,i,j,8)
     .   +  gg(id,i+1,j+1)*fam(id,i,j,9)
        enddo

          f(0,i,j) = gg(0,i,j)

           enddo
	enddo

	return
	end
c---------------------------------------------
	subroutine hydrovar

	implicit double precision(a-h,o-z)
	include'lbe.par'
c--------------------------------------------
c Calculation of density / pressure and velocities

	do j = 1,ny
	 do i = 1,nx
          rho(i,j) = f(0,i,j)+f(1,i,j)+f(2,i,j)
     .     + f(3,i,j)+ f(4,i,j)+f(5,i,j)
     .     +f(6,i,j)+f(7,i,j)+f(8,i,j)

          rhoi=1./rho(i,j)

	  u(i,j)=(f(1,i,j)-f(3,i,j)+f(5,i,j)- 
     .            f(6,i,j)-f(7,i,j)+f(8,i,j))*rhoi 

	  v(i,j)=(f(5,i,j)+f(2,i,j)+f(6,i,j)
     .          - f(7,i,j)-f(4,i,j)-f(8,i,j))*rhoi

	  enddo
	enddo

	return
	end
c-------------------------------------------------
	subroutine equili

	implicit double precision(a-h,o-z)
	include'lbe.par'
c-------------------------------------------------
	do j = 0,ny+1
	  do i = 0,nx+1
c
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

	    feq(0,i,j) = rl*rt0*(1.0d0 - sumsq)
	    feq(1,i,j) = rl*rt1*(1.0d0 - sumsq + u2 + ui)
	    feq(2,i,j) = rl*rt1*(1.0d0 - sumsq + v2 + vi)
	    feq(3,i,j) = rl*rt1*(1.0d0 - sumsq + u2 - ui)
	    feq(4,i,j) = rl*rt1*(1.0d0 - sumsq + v2 - vi)
	    feq(5,i,j) = rl*rt2*(1.0d0 + sumsq2 + ui + vi + uv)
	    feq(6,i,j) = rl*rt2*(1.0d0 + sumsq2 - ui + vi - uv)
	    feq(7,i,j) = rl*rt2*(1.0d0 + sumsq2 - ui - vi + uv)
	    feq(8,i,j) = rl*rt2*(1.0d0 + sumsq2 + ui - vi - uv)
	
	   enddo
	enddo
	return
	end
c----------------------------------------------------------
	subroutine collis

	implicit double precision(a-h,o-z)
	include'lbe.par'
c----------------------------------------------------------
	do k = 0,npop-1
	 do j = 1,ny
	  do i = 1,nx
	   gg(k,i,j)=f(k,i,j)*(1.0d0-omega)+omega*feq(k,i,j)
	  enddo
	 enddo
	enddo
	    
	return 
	end  
c =========================
	subroutine pbc
c =========================
	implicit double precision(a-h,o-z)
	include'lbe.par'
c-----------------------------------------------------------
c EAST case
	
	do j = 1,ny
	   f(1,0,j) = f(1,nx,j)
	   f(5,0,j) = f(5,nx,j)
	   f(8,0,j) = f(8,nx,j)
	enddo

c WEST case
	do j = 1,ny
	   f(3,nx+1,j) = f(3,1,j)
	   f(6,nx+1,j) = f(6,1,j)
	   f(7,nx+1,j) = f(7,1,j)
	enddo

c NORTH case
	do i = 1,nx
	   f(2,i,0) = f(2,i,ny)
	   f(5,i,0) = f(5,i,ny)
           f(6,i,0) = f(6,i,ny)
	enddo

c SOUTH case
	do i = 1,nx
	   f(4,i,ny+1) = f(4,i,1)
           f(7,i,ny+1) = f(7,i,1)
           f(8,i,ny+1) = f(8,i,1)
	enddo

	return
	end
c-------------------------------------------------------------
	subroutine mbc
	
	implicit double precision(a-h,o-z)
	include'lbe.par'
c-------------------------------------------------------------
c Top 

	do i = 1,nx
        f(4,i,ny) = gg(2,i,ny)
        f(8,i,ny) = gg(6,i,ny) + rho(i,ny)/6.*ulip
        f(7,i,ny) = gg(5,i,ny) - rho(i,ny)/6.*ulip
	enddo

C Note here that the coner points on top uses no-slip
c LEFT 

	 do j = 1,ny
           f(1,1,j) = gg(3,1,j)
           f(5,1,j) = gg(7,1,j)
           f(8,1,j) = gg(6,1,j)
        enddo

c Right    

	do j = 1,ny
           f(3,nx,j) = gg(1,nx,j)
           f(7,nx,j) = gg(5,nx,j)
           f(6,nx,j) = gg(8,nx,j)
        enddo

c Bottom     
	do i = 1,nx
           f(2,i,1) = gg(4,i,1)
           f(6,i,1) = gg(8,i,1)
           f(5,i,1) = gg(7,i,1)
        enddo

	return
	end
c ==========================
	subroutine obst1
c ==========================
	implicit double precision(a-h,o-z)
	include'lbe.par'
c--------------------------------------------------------
	k = nx / 4	
	

        nt = ny/2+nobst/2
        nb = ny/2-nobst/2+1

        do j = nb,nt
           f(1,k,j) = f(3,k+1,j)
           f(5,k,j) = f(7,k+1,j+1)
           f(8,k,j) = f(6,k+1,j-1)

           f(3,k,j) = f(1,k-1,j)
           f(7,k,j) = f(5,k-1,j-1)
           f(6,k,j) = f(8,k-1,j+1)
        enddo

           temp = f(8,k,nt+1)
           f(8,k,nt+1) = f(6,k+1,nt)
           f(6,k+1,nt) = temp

           temp = f(7,k,nt+1)
           f(7,k,nt+1) = f(5,k-1,nt)
           f(5,k-1,nt) = temp

           temp = f(6,k,nb-1)
           f(6,k,nb-1) = f(8,k-1,nb)
           f(8,k-1,nb) = temp

           temp = f(5,k,nb-1)
           f(5,k,nb-1) = f(7,k+1,nb)
           f(7,k+1,nb) = temp

           f(2,k,nt) = f(4,k,nt+1)
           f(4,k,nb) = f(2,k,nb-1)

	return
	end
c ==========================
	subroutine obst2
c ==========================
	implicit double precision(a-h,o-z)
	include'lbe.par'
c--------------------------------------------------------
	k = nx / 4	
	
        do j = ny/2-nobst/2+1,ny/2+nobst/2
           f(1,k+1,j) = f(3,k+1,j)
	   f(3,k  ,j) = f(1,k,j)
        enddo

        do j = ny/2-nobst/2,ny/2+nobst/2+1
	   f(5,k+1,j) = f(7,k+1,j)
	   f(8,k+1,j) = f(6,k+1,j)
	   f(7,k,  j) = f(5,k,  j)
	   f(6,k,  j) = f(8,k,  j)
        enddo

	return
	end
c ==========================
	subroutine obstequil
c ==========================
	implicit double precision(a-h,o-z)
	include'lbe.par'
        k=nx/4
c zero-speed equilibria
        do j = ny/2-nobst/2+1,ny/2+nobst/2
          f(0,k,j)     = rt0
          f(0,k+1,j)   = rt0
          do ip=1,4
           f(ip,k,j)   = rt1
           f(ip,k+1,j) = rt1
          enddo
          do ip=5,8
           f(ip,k,j)   = rt2
           f(ip,k+1,j) = rt2
          enddo
        enddo

	return
	end
c-----------------------------------------------------------
	subroutine profil(it)
c-----------------------------------------------------------
	implicit double precision(a-h,o-z)
	include'lbe.par'
c----------------------------------------------------------
        write(6,*) 'ucenter,force',u(nx/2,ny/2),fpois
	do j = 1,ny
	  write(10,301) yy(j)/WL,0.5*(u(nx/2,j)+u(nx/2+1,j))/ulip,
     >        0.5*(v(nx/2,j)+v(nx/2+1,j))/ulip
	enddo
	write(10,'(bn)')

	do i = 1,nx
	  write(11,301) xx(i)/WL,0.5*(u(i,ny/2)+u(i,ny/2+1))/ulip,
     >        0.5*(v(i,ny/2)+v(i,ny/2+1))/ulip
	enddo
	write(11,'(bn)')

        write(17,*) it,u(nx/2,ny/2)

301     format(2x,3(1pe15.4))
   
	return
	end
c---------------------------------------------------------
	subroutine diag0D(istep)
c---------------------------------------------------------
	implicit double precision(a-h,o-z)
        include'lbe.par'
c----------------------------------------------------------
	densit = 0.0d0

	do k = 0,npop-1
	   do j= 1,ny
	      do i = 1,nx
	         densit = densit + f(k,i,j)
	      enddo
	   enddo
	enddo

        numb = dfloat(nx*ny)

        if (iobst) then
        numb = numb -nobst
        endif

	densit = densit/float(numb)

	umoy = 0.0d0
	vmoy = 0.0d0

	do j = 1,ny
	   do i = 1,nx
	      umoy = umoy + u(i,j)
	      vmoy = vmoy + v(i,j)
	   enddo
	enddo
	
	umoy = umoy / dfloat(nx*ny)
	vmoy = vmoy / dfloat(nx*ny)

        ijc = nx/2
        write(12,*)time*ulip/float(nx),
     >0.25*(u(ijc,ijc)+u(ijc+1,ijc)+u(ijc,ijc+1)+u(ijc+1,ijc+1))/ulip,
     >0.25*(v(ijc,ijc)+v(ijc+1,ijc)+v(ijc,ijc+1)+v(ijc+1,ijc+1))/ulip

	print*,'diagnostic 0D : istep density umoy and vmoy',
     .          istep,densit,umoy,vmoy


	return
	end
c =========================================================
	subroutine config(istep,iconf) 
c =========================================================
	implicit double precision(a-h,o-z)
        include'lbe.par'
c -------------------------------------------
C streamfunction psi
        dimension psi(nx,ny),psi2(nx,ny)
c Calculation of density / pressure and velocities

        do j = 1,ny
         do i = 1,nx
          rho(i,j) = f(0,i,j)+f(1,i,j)+f(2,i,j)
     .     + f(3,i,j)+ f(4,i,j)+f(5,i,j)
     .     +f(6,i,j)+f(7,i,j)+f(8,i,j)

          rhoi=1./rho(i,j)

          u(i,j)=(f(1,i,j)-f(3,i,j)+f(5,i,j)-
     .            f(6,i,j)-f(7,i,j)+f(8,i,j))*rhoi

          v(i,j)=(f(5,i,j)+f(2,i,j)+f(6,i,j)
     .          - f(7,i,j)-f(4,i,j)-f(8,i,j))*rhoi

          enddo
        enddo
c
        do i=1,nx
        psi(i,1)=0.5*u(i,1)*yy(1)
        do j = 2,ny
        psi(i,j)=psi(i,j-1)+0.5*(u(i,j)+u(i,j-1))*(yy(j)-yy(j-1))
        enddo
        enddo
c
        do j=1,ny
 
        psi2(1,j)= - 0.5*v(1,j)*xx(1)

        do i = 2,nx
        psi2(i,j)=psi2(i-1,j)-0.5*(v(i,j)+v(i-1,j))*(xx(i)-xx(i-1))
        enddo
        enddo
c
        iout=60+iconf
c
	do j = 1,ny
           
           ycc = (j-ny/2) - 0.5d0
           ycc = ycc/(ny/2)
           theory = uf*(1.d0- ycc**2)
           error = u(nx/2,j) - theory

           write(iout-10,*) j,ycc,u(nx/2,j)/uf,v(nx/2,j),
     1      theory/uf,error

	  do i = 1,nx
c           write(iout,*) i,j,u(i,j),v(i,j)
c           write(iout,*) i,j,u(i,j),v(i,j)
            write(iout,*)psi(i,j)
            write(iout+10,*)psi2(i,j)
            write(iout+20,*)rho(i,j)
	   enddo
	enddo

         close(iout)
         close(iout+10)
         close(iout+20)

        do i=1,nx
        write(iout-20,*)i,rho(i,1),rho(i,ny/4),rho(i,ny/2),rho(i,ny)
        enddo
        write(6,*) 'configuration at time and file>>>>>> ',istep,iout
	
	return
	end






