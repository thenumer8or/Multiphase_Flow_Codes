C	-------------------------
C
C	Fourier to real 
C
C	-------------------------
C
	SUBROUTINE F2RC1(ERSATZ)
C
C
	include 'param.inc'
        PARAMETER ( n1h=n1/2 , n2h=n2/2 )
        PARAMETER ( INCX=1,JUMPX=n1p2 )
        INTEGER IFAXX(13)
        REAL TRIGSX(3*n1h+1)
C
        REAL ERSATZ(n1,n2),TEMPW(n1p2,n2)
        REAL WORKX(n2*n1p1)
        double precision theta,wi,wi1,wpi,wpr,wr,wr1,wtemp,pi
C
        COMMON/TRIG/ TRIGSX
        COMMON/FAX/ IFAXX
C
        pi=4.0d0*datan(1.0d0)
        theta=0.5d0*pi/n1
        wr=1.0d0
        wi=0.0d0
        wr1=cos(theta)
        wi1=sin(theta)
        wpr=-2.0d0*wi1**2
        wpi=sin(2.d0*theta)

        do j=1,n2
        tempw(1,j)=ERSATZ(1,j)
        tempw(2,j)=0.0
        do i=n1,4,-2
        tempw(i,j) = ERSATZ(i-2,J) - ERSATZ(i,J)
        enddo
        tempw(n1p1,j)=ERSATZ(n1,j)
        tempw(n1p2,j)=0.0
        enddo

        DO I=3,n1,2
        wtemp=wr
        wr=wr*wpr-wi*wpi+wr
        wi=wi*wpr+wtemp*wpi+wi

        DO J=1,N2
        y1=ERSATZ(i,j)*wr+TEMPW(i+1,j)*wi
        y2=TEMPW(i+1,j)*wr-ERSATZ(i,j)*wi
        tempw(I,J)=y1/2.0
        tempw(I+1,J)=-y2/2.0
        ENDDO

        ENDDO

        CALL FFT991(TEMPW(1,1),WORKX,TRIGSX,IFAXX,
     1                   INCX,JUMPX,n1,n2,+1)

 
         DO I=1,n1h
         IA=n1p1-I

         do j=1,n2
         y1=tempw(I,J)+tempw(IA,J)
         y2=0.5/wi1*(tempw(I,J)-tempw(IA,J))
         ERSATZ(i,j)=0.5*(y1+y2)
         ERSATZ(ia,j)=0.5*(y1-y2)
         enddo

         wtemp=wr1
         wr1=wr1*wpr-wi1*wpi+wr1
         wi1=wi1*wpr+wtemp*wpi+wi1
         ENDDO

	RETURN
	END
C
C
C       -------------------------
C
C       real to Fourier
C
C       -------------------------
C
        SUBROUTINE R2FC1(ERSATZ)
C
C
        include 'param.inc'
        PARAMETER ( n1h=n1/2 , n2h=n2/2 )
        PARAMETER ( INCX=1,JUMPX=n1p2 )
        INTEGER IFAXX(13)
        REAL TRIGSX(3*n1h+1)
C
        REAL ERSATZ(n1,n2),TEMPW(n1p2,n2)
        REAL WORKX(n2*n1p1)
        double precision theta,wi,wi1,wpi,wpr,wr,wr1,wtemp,pi
C
        COMMON/TRIG/ TRIGSX
        COMMON/FAX/ IFAXX
C
        pi=4.0d0*datan(1.0d0)
        theta=0.5d0*pi/n1
        wr=1.0d0
        wi=0.0d0
        wr1=cos(theta)
        wi1=sin(theta)
        wpr=-2.0d0*wi1**2
        wpi=sin(2.d0*theta)
 
         DO I=1,n1h
         IA=n1p1-I

         do j=1,n2
         y1=0.5*(ERSATZ(I,J)+ERSATZ(IA,J))
         y2=wi1*(ERSATZ(I,J)-ERSATZ(IA,J))
         TEMPW(i,j)=y1+y2
         TEMPW(ia,j)=y1-y2
         enddo

         wtemp=wr1
         wr1=wr1*wpr-wi1*wpi+wr1
         wi1=wi1*wpr+wtemp*wpi+wi1
         ENDDO

         do j=1,n2
         TEMPW(n1p1,J)=0.0
         TEMPW(n1p2,J)=0.0
         ENDDO

        CALL FFT991(TEMPW(1,1),WORKX,TRIGSX,IFAXX,
     1                   INCX,JUMPX,n1,N2,-1)

c - even modes -
        do j=1,n2
        ERSATZ(1,J)=tempw(1,j)
        enddo
        DO I=3,n1,2
        wtemp=wr
        wr=wr*wpr-wi*wpi+wr
        wi=wi*wpr+wtemp*wpi+wi

        DO J=1,N2
        y1=TEMPW(i,j)*wr+TEMPW(i+1,j)*wi
        y2=-TEMPW(i+1,j)*wr+TEMPW(i,j)*wi
        ERSATZ(I,J)=2.0*y1
        ERSATZ(I+1,J)=y2
        ENDDO

        ENDDO

c  - odd modes -
        DO J=1,N2
        sum=0.5*TEMPW(n1p1,J)
        DO I=n1,2,-2
        sum1=sum
        sum=sum+ERSATZ(i,j)
        ERSATZ(i,j)=2.0*sum1
        ENDDO
        ENDDO
C
        return
        end
c
