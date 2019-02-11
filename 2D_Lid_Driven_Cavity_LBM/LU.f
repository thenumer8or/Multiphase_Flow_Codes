	SUBROUTINE ludcmp(a,n,np,indx,d)
        implicit double precision(a-h,o-z)

        dimension indx(n),a(np,np)
c	Integer n,np,indx(n),NMAX
c	Real d,a(np,np),TINY
	parameter (NMAX = 5000,TINY = 1.0d-20)
        dimension vv(NMAX)
c	Ineteger i,imax,j,k
c	Real aamax,dum,sum,vv(NMAX)

	d = 1.d0

	do 12  i = 1, n
        aamax = 0.d0
        do 11 j = 1, n

        if (dabs(a(i,j)).gt.aamax) aamax = dabs(a(i,j))
11      enddo 
c       if (aamax.eq.0) pause 'singular matrix in ludcmp'
  
        vv(i) = 1.d0 / aamax
12      enddo 

        do 19 j = 1, n

        do 14 i = 1, j - 1 

        sum = a(i,j)

        do 13 k = 1, i - 1

        sum = sum - a(i,k) * a(k,j)
13      enddo 
        a(i,j) = sum
14      enddo 
        aamax = 0.d0

        do 16 i = j, n

        sum = a(i,j) 
        do 15 k = 1, j-1
        sum = sum - a(i,k) * a(k,j) 
15      enddo 

        a(i,j) = sum 
        dum = vv(i) * dabs(sum)

        if (dum.ge.aamax) then

        imax = i
        aamax = dum
        endif
16      enddo

        if (j.ne.imax) then
        do 17 k=1, n
        dum = a(imax, k)
        a(imax,k) = a(j,k)
        a(j,k) = dum
17      enddo 

       d = -d
       vv(imax) = vv(j)
       endif

       indx(j) = imax

       if (a(j,j).eq.0.d0)a(j,j) = TINY 
       if (j.ne.n) then 
       dum = 1.d0 / a(j,j) 

       do 18 i = j + 1, n

       a(i,j) = a(i,j) * dum
18     enddo 

       endif
19     enddo 
       return
       end


      SUBROUTINE lubksb(a,n,np,indx,b)
      implicit double precision(a-h,o-z)
        dimension indx(n),a(np,np),b(n)
c     Integer n,np,indx(n)
c     Real a(np,np),b(n)
c     Integer i,ii,j,ll
c     Real sum
C *************************

      ii= 0
      do 12 i = 1, n 
      ll = indx(i)
      sum = b(ll)
       b(ll) = b(i)
      if   (ii.ne.0) then
      do 11 j = ii, i - 1 
      sum = sum - a(i,j) * b(j)
 11   enddo 

      else if (sum.ne.0.d0) then

       ii = i
       endif
       b(i) = sum
12    enddo 

       do 14 i = n, 1, -1
       sum= b(i)
        do 13 j = i + 1, n
        sum = sum - a(i,j) * b (j)
13      enddo 
        b(i) = sum / a(i,i)
 
14      enddo  
        return
        end

