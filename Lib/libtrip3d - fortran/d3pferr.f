

c Contents of this file:
c
c	subroutine D3PFERR
c	  entry BD3PFERR
c	subroutine D3PFGEOM
c	subroutine D3PFERRB

c=======================================================================



      subroutine D3PFERR

c-----------------------------------------------------------------------
c                                            By: M. Schaffer 2006 jun 10
c                                 Last Modified: M. Schaffer 2007 oct 25
c-----------------------------------------------------------------------
c D3PFERR defines the geometry of and calculates the magnetic field 
c from a model of the shifted and tilted DIII-D PF coils (F-coils).
c D3PFERR receives user variables defining PF coil tilts and shifts 
c through common /d3pfer/, which is defined in D3PFGEOM.
c D3PFERR handles all the tasks needed to calculate this field.
c
c This subroutine has several main parts:
c  1) Set flags to signal the components to use.
c  2) Call routines to set up coil models to be used.
c  3) Write parameters to appropriate output files.
c  4) Upon separate ENTRY at BD3PFERR, call subroutines to calculate
c     the F-coil error magnetic field at a point (x,y,z).
c 
c Units: SI (mks)
c Fully double precision.
c-----------------------------------------------------------------------
c Parameters are passed by include 'D3PFERR.i':
c	nfc    = number of D3D F-coils = 18
c	nflps  = 2*nfc = total number of current loops in model
c-----------------------------------------------------------------------

      IMPLICIT NONE
      INCLUDE 'd3pferr.i'	! included parameter & variable types
      				! are declared in the include files.

c-----------------------------------------------------------------------

      common /consts/ pi, twopi, cir, rtd, dtr
      common /flags2/ iPROBE, iSURFMN, iTRIP3D,
     &                iD3Dcoils, iNSTXcoils, iITERcoils, iFDFcoils
      common /d3pfer/ fcshft, fashft, fctilt, fatilt, fcur
      common /eqpol/  dsbp,alfsbp,dthbp,alfthbp,iplasbp,ipdir,lbpol


!  Common variables:
      real*8  pi, twopi, cir, rtd, dtr
      integer iPROBE, iSURFMN, iTRIP3D
      integer iD3Dcoils, iNSTXcoils, iITERcoils, iFDFcoils
      real*8  fcshft(nfc), fashft(nfc), fcur(nfc)
      real*8  fctilt(nfc), fatilt(nfc)
      real*8  dsbtpl, alfsbtpl, dthbtpl, alfthbtpl
      real*8  dsbp, alfsbp, dthbp, alfthbp
      integer iplasbp, ipdir, lbpol

!  Arguments at ENTRY BD3PFERR:
      real*8  x, y, z, bx, by, bz

!  Local and Namelist Variables:
      logical ipfs
      integer kusepfs(nflps)
      integer i, k

      SAVE

c-----------------------------------------------------------------------
c Set some default values
c-----------------------------------------------------------------------
c Variables initialized by data statement acquire the SAVE attribute.

      data kusepfs /nflps*0/


c-----------------------------------------------------------------------
c Set flags
c-----------------------------------------------------------------------

      ipfs = .false. 
      DO k=1, nfc
       if(fcur(k) .ne. 0) then
        ipfs = .true.		! at least 1 PF loop has current
        iD3Dcoils = 1		! at least 1 D3D coil model is used
       end if
      END DO
      
c-----------------------------------------------------------------------
c Set up geometries of PF errors that will be used
c-----------------------------------------------------------------------

      if (ipfs)  CALL D3PFGEOM (kusepfs)

c-----------------------------------------------------------------------
c Write some information about this model to file CASE
c-----------------------------------------------------------------------

      if(ipfs) then			! Write F-coil currents
       write(80,*) ''
       write(80,*) 'DIII-D F-coil currents (A):'
       write(80,1010) (fcur(i),i=01,05)
       write(80,1010) (fcur(i),i=06,09)
       write(80,1010) (fcur(i),i=10,14)
       write(80,1010) (fcur(i),i=15,18)
      end if


      RETURN		! End of D3D PF Error Model setup



c=======================================================================
c Calculate magnetic field for all requested PF error.
c
c ENTRY BD3PFERR receives the field point (x,y,z) at which the field is
c  to be calculated. It returns the Cartesian magentic field components
c  (bx,by,bz).
c-----------------------------------------------------------------------

      ENTRY BD3PFERR (x, y, z, bx, by, bz)

      bx = 0.0d0		! initialize magnetic field components
      by = 0.0d0		! that will be calculated and returned
      bz = 0.0d0

       if (ipfs)  CALL D3PFERRB 
     &            (kusepfs, x, y, z, bx, by, bz)

      RETURN


c-----------------------------------------------------------------------
c Format statements
c-----------------------------------------------------------------------
 1003 format(1x,30i3)
 1010 format(1x,8f9.0)
 1014 format(1x,8f9.4)


      END  	SUBROUTINE	D3PFERR 
c=======================================================================







      subroutine D3PFGEOM (kuse)

c-----------------------------------------------------------------------
c                                            By: M. Schaffer 2006 jun 10
c                                 Last Modified: M. Schaffer 2007 oct 25
c-----------------------------------------------------------------------
c This code is a DIII-D adaptation of LPS.f for shifted and/or tilted 
c PF-coils written by Todd Evans. 
c
c Added tilt capability to PF coils (MS 2002jul05)
c-----------------------------------------------------------------------
c This subroutine is called by D3PFERR once to build two single-filament
c circular loops for each of the PF-coils.
c 
c Each DIII-D PF-coil consists of two parts, a 'physical' filamentary loop
c and a reversed one. The physical loop carries the actual PF-coil current
c and is user-shifted with respect to the unperturbed location by 
c magnitude fcshft(nfc=nflps/2) in the toroidal angle direction 
c fashft(nflps). The reversed current loop is at the location of the
c corresponding PF-coil in the equilibrium (from BEQPOL), which may also
c be shifted and tilted. The pair gives the difference between the actual
c and ideally placed loops. The net effect is to simulate the PF-coil 
c shift error B.
c
c Positive fcur(nfc) current means in the direction of code's positive
c toroidal angle (i.e., the same direction as positive Ip in DIII-D). 
c The current direction is reversed if variable ipdir .lt. 0.
c
c D3PFGEOM multiplies fcur(i) by turns(i), so that the user does not 
c have the chore of multiplying ptdata F-coil currents by turn numbers.
c
c Note: See subroutines BDLOOP.f and TRLOOP.f for the orientation specs.
c-----------------------------------------------------------------------
c Set up common blocks, arrays and F-coil data (this data is consistent
c with efit data located in /link/efit/dprobe.dat on 4/22/02).
c F-coil order is:    F1A,F2A,F3A,F4A,F5A,F6A,F7A,F8A,F9A,
c                     F1B,F2B,F3B,F4B,F5B,F6B,F7B,F8B,F9B
c
c Parameters are passed by include 'D3PFERR.i':
c	nfc    = number of D3D F-coils = 18
c	nflps  = 2*nfc = total number of current loops in model
c
c-----------------------------------------------------------------------

      IMPLICIT NONE
      INCLUDE 'd3pferr.i'	! included parameter & variable types
				! are declared in the include file.

c-----------------------------------------------------------------------

      common /consts/  pi, twopi, cir, rtd, dtr
      common /d3pfer/  fcshft, fashft, fctilt, fatilt, fcur
      common /d3pflps/ amat, origin, rcoil, curlps
      common /eqpol/   dsbp,alfsbp,dthbp,alfthbp,iplasbp,ipdir,lbpol


!  Common variables:
      real*8  pi, twopi, cir, rtd, dtr
      real*8  fcshft(nfc), fashft(nfc), fcur(nfc)
      real*8  fctilt(nfc), fatilt(nfc)
      real*8  amat(3,3,nflps), origin(3,nflps)
      real*8  rcoil(nflps), curlps(nflps)
      real*8  dsbp, alfsbp, dthbp, alfthbp
      integer iplasbp, ipdir, lbpol

!  Arguments
      integer kuse(nflps)

!  Local Variables
      real*8  dp0, dp1
      real*8  zcoil(nflps), turns(nflps)
      real*8  costh(nflps), sinth(nflps), cosphi(nflps), sinphi(nflps)
      integer i, j, k, indx

      SAVE

c-----------------------------------------------------------------------

! save shift & tilt transform arrays that are only calculated once
CC      SAVE amat, origin	! THIS IS ILLEGAL


      data dp0, dp1 / 0.0d0, 1.0d0 /
      

c F-coil order is:    F1A,F2A,F3A,F4A,F5A,F6A,F7A,F8A,F9A,
c                     F1B,F2B,F3B,F4B,F5B,F6B,F7B,F8B,F9B

c Data for the 18 DIII-D F-coils and their reversed-current twins

      data curlps /nflps*0.0/	! default to zero current

      data rcoil /2*0.8608, 2*0.8614, 2*0.8628, 2*0.8611, 2*1.0041,
     &               2*2.6124, 2*2.3733, 2*1.2518, 2*1.6890, 
     &            2*0.8608, 2*0.8607, 2*0.8611, 2*0.8630, 2*1.0025,
     &               2*2.6124, 2*2.3834, 2*1.2524, 2*1.6889/

      data zcoil/ 2*0.1683,  2*0.5081,  2*0.8491,  2*1.1899,  2*1.5169,
     &              2*0.4376,  2*1.1171,  2*1.6019,  2*1.5874,
     &           2*-0.1737, 2*-0.5135, 2*-0.8543, 2*-1.1957, 2*-1.5169,
     &              2*-0.4376, 2*-1.1171, 2*-1.6027, 2*-1.5780/

      data turns /2*58.0, 2*58.0, 2*58.0, 2*58.0, 2*58.0,
     &               2*55.0, 2*55.0, 2*58.0, 2*55.0, 
     &            2*58.0, 2*58.0, 2*58.0, 2*58.0, 2*58.0,
     &               2*55.0, 2*55.0, 2*58.0, 2*55.0/

! Perturbed PF coil positions specified with respect to code axis.
! Actual D3D coils wrt ERROR MEASUREMENT AXIS:

      data fcshft /.0049,.0046,.0060,.0064,.0020, !shift (m) wrt code axis
     & 	              .0052,.0082,.0008,.0012,
     &             .0027,.0033,.0029,.0032,.0019,
     & 	              .0031,.0049,.0063,.0018/

      data fashft /071, 063, 060, 066, 063, 	  !shifted toward (deg)
     & 	              166, 247, 072, 120,
     &             077, 051, 064, 072, 072,
     & 	              224, 180, 152, 174/

      data fctilt /.010, .062, .033, .132, .137,  !tilt(deg) wrt code axis
     & 	              .029,.044,.185,.068,
     &             .099, .178, .106, .110, .046,
     & 	              .016, .061, .145, .098/

      data fatilt /029, 293, 053, 031, 277,       !axis tilts toward (deg)
     & 	              308, 168, 274, 316,
     &             343, 021, 098, 017, 198,
     & 	              314, 104, 175, 178/
c------------------------------------------------------------------------


c Initialize trigonometric terms
      DO  i=1,nflps
       costh(i)=  dp1
       sinth(i)=  dp0
       cosphi(i)= dp1
       sinphi(i)= dp0
       kuse = 0
      END DO

c Precalculate terms needed to to transform the field point to loops'
c coordinate frame and Magnetic field from loop frame back to code frame.
c This code is from T Evans' subroutine D3LPS, 2002.

c Initialize indices
      indx = 0		! indx = 0 signals a physical PF coil
      			! indx = 1 signals a reversed current loop
			!           to subtract from the equilibrium
      j = 1		! here i is loop number, j is pair number

c Calculate shifted origins (Cartesian) of all PF-coils.
c Calculate trigonometric functions of their tilts.
c Assign the proper total current curlps(i) to each loop.


CC      WRITE(0,*)'D3PFGEOM  ORIGIN, CURLPS / costh,sinth / cosphi,sinphi'

      DO  i=1,nflps
      
      IF (indx.eq.0) THEN 		! indx=0 means a physical coil
       costh(i)=cos(dtr*fctilt(j))	! fctilt= coil tilt (deg)
       sinth(i)=sin(dtr*fctilt(j))
       cosphi(i)=cos(dtr*fatilt(j))	! fatilt= toroidal angle toward
       sinphi(i)=sin(dtr*fatilt(j))	!  which coil axis tilts (deg)
       
       origin(1,i)=fcshft(j)*cos(dtr*fashft(j))	 ! fashft are in degrees
       origin(2,i)=fcshft(j)*sin(dtr*fashft(j))
       origin(3,i)=zcoil(i)
       curlps(i) = fcur(j)*turns(i)	! for shifted PF coil loop
       if (ipdir .lt. 0) curlps(i) = -curlps(i)	! reversed plasma current
       if (curlps(i) .ne. dp0)  kuse(i) = 1	! flag to use loop i

C       IF (kuse(i) .ne. 0) THEN		! ***** DEBUG
C       WRITE(0,929) i,indx,turns(i),rcoil(i),zcoil(i),fcur(j)
CC       WRITE(0,929) i,j,(origin(k,i),k=1,3),curlps(i)
CC       WRITE(0,929) i,j,costh(i),sinth(i)
CC       WRITE(0,929) i,j,cosphi(i),sinphi(i)
C       END IF

      indx = 1
       
      ELSE 	! reversed, shift & tilt with the plasma equilibrium
      
c Shift and tilt of equilibrium poloidal fields are specified by 
c  dsbp, alfsbp, dthbp, alfthbp, respectively. Put the second coil of 
c  each pair in its perturbed equilibrium position, to cancel its
c  field contribution from the equilibrium field.

       costh(i)=cos(dthbp)	! dthbp  = equilibrium tilt (rad)
       sinth(i)=sin(dthbp)
       cosphi(i)=cos(alfthbp)	! alfthbp= toroidal angle toward which
       sinphi(i)=sin(alfthbp)	! magnetic axis tilts (rad)

       origin(1,i)=dsbp*cos(alfsbp)	! alfsbp= tor angle toward which
       origin(2,i)=dsbp*sin(alfsbp) 	! equilibrium shifts (rad)
       origin(3,i)=zcoil(i)
       curlps(i) = -fcur(j)*turns(i) 	! opposite to coil currnt
       if (ipdir .lt. 0) curlps(i) = -curlps(i)	! reversed plasma current
       if (curlps(i) .ne. dp0)  kuse(i) = 1	! flag to use loop i

C       IF (kuse(i) .ne. 0) THEN
C       WRITE(0,929) i,indx,turns(i),rcoil(i),zcoil(i),fcur(j)
CC       WRITE(0,929) i,j,(origin(k,i),k=1,3),curlps(i)
CC       WRITE(0,929) i,j,costh(i),sinth(i)
CC       WRITE(0,929) i,j,cosphi(i),sinphi(i)
C       END IF

       j = j + 1
       indx = 0

      END IF

      END DO


c Precalculate matrix elements used to transform between the code's
c coordinate system and the coordinate system tied to a rotated loop.
c This code is from T Evans' subroutine D3LPS.

      DO 200 I=1,NFLPS
       AMAT(1,1,I)=COSTH(I)*COSPHI(I)
       AMAT(1,2,I)=COSTH(I)*SINPHI(I)
       AMAT(1,3,I)=-SINTH(I)
       AMAT(2,1,I)=-SINPHI(I)
       AMAT(2,2,I)=COSPHI(I)
       AMAT(2,3,I)=0.
       AMAT(3,1,I)=SINTH(I)*COSPHI(I)
       AMAT(3,2,I)=SINTH(I)*SINPHI(I)
       AMAT(3,3,I)=COSTH(I)
200   CONTINUE

C      WRITE (15,*) 'd3pfgeom:'
C      WRITE (15,*) nflps,nxtra,fcur
C      WRITE (15,*) nflps,nxtra,turns
C      WRITE (15,*) ' '
      
 929  format(1x,2i3,8f22.10)	
 1010 format(1x,8f9.0)
 1014 format(1x,8f9.4)
      
      RETURN

      END 	SUBROUTINE	D3PFGEOM
c=======================================================================







      subroutine D3PFERRB (kuse, x, y, z, bxf, byf, bzf)

c-----------------------------------------------------------------------
c                                            By: M. Schaffer 2006 jun 10
c                                 Last Modified: M. Schaffer 2007 oct 25
c-----------------------------------------------------------------------
c D3PFERRB calculates magnetic field from the DIII-D PF-coils that may
c be shifted and/or tilted with respect to the code coordinate system.
c
c The magnetic field due to these loops is obtained using an elliptic
c integral expression. The shifts and tilts use code written by 
c Todd Evans.
c
c Input variable kuse signals loops with no current (by kuse = 0).
c Input variables x,y,z = field point in Cartesian coordinates.
c Output bxf,byf,bzf = magnetic field in Cartesian coordinates, summed
c over all nfc coils and the corresponding nfc equal and opposite coils
c shifted with the equilibrium.
c
c Parameters are passed by include 'D3PFERR.i':
c	nfc    = number of D3D F-coils = 18
c	nflps  = 2*nfc = total number of current loops in model
c-----------------------------------------------------------------------

      IMPLICIT NONE
      INCLUDE 'd3pferr.i'	! included parameter & variable types
				! are declared in the include file.

c-----------------------------------------------------------------------

      common /d3pflps/ amat, origin, rcoil, curlps


!  Common variables:
      real*8  amat(3,3,nflps), origin(3,nflps)
      real*8  rcoil(nflps), curlps(nflps)

!  Arguments
      integer kuse(nflps)
      real*8  x, y, z, bxf, byf, bzf

!  Local Variables
      real*8  dp0, dp1
      real*8  C(3), CP(3,nflps) 
      real*8  cpk(3), Bpk(3), BP(3,nflps)
      real*8  BX(3,nflps), BXSUM(3)
      integer i, j, k, m, NLOOPS

      SAVE

c-----------------------------------------------------------------------

      data dp0, dp1 /0.0d0, 1.0d0/

      NLOOPS = nflps

c-----------------------------------------------------------------------
c Initialize field point and magnetic field component values
c-----------------------------------------------------------------------

      C(1) = x		! x of field point
      C(2) = y		! y
      C(3) = z		! z

      do  i=1,3
       BXSUM(i) = dp0
       do j=1,nflps
        BP(i,j) = dp0
       end do
      end do

c-----------------------------------------------------------------------
c Transform field point (C) to loop coordinates (CP) for all loops.
c This code is from the old TRLOOP routine written by Todd Evans (2002).
C
C   LET X,Y,Z BE SOME CARTESIAN COORDINATE SYSTEM.                        
C   LET XP,YP,ZP BE A CARTESIAN COORDINATE SYSTEM ATTACHED TO A LOOP.    
C   GIVEN A FIELD POINT C=(X,Y,Z), TRANSFORM IT TO THE COORDINATES     
C   CP=(XP,YP,ZP) OF A LOOP WHOSE PERPENDICULAR     
C   (ZP-AXIS) IS SPECIFIED BY THE SPHERICAL COORDINATES (THETA,PHI)
C   (I.E. TWO ROTATIONS) AND WHOSE ORIGIN IS X0=ORIGIN(1,K),       
C   Y0=ORIGIN(2,K),Z0=ORIGIN(3,K),K=1,NLOOPS.                      
C   THIS IS DONE FOR NLOOPS.                                       
c-----------------------------------------------------------------------

CC      WRITE(0,*)'D3PFERRB  CP C / Cp' 	! **********  DEBUG

        DO 10 I=1,3
         DO 20 K=1,NLOOPS     		! reset CP(I,K)
          CP(I,K)=0.0d0
20       CONTINUE

         DO 30 J=1,3
          DO 40 K=1,NLOOPS
          if (kuse(K).ne.0) then	! do only loops with current

CC           WRITE(0,930) i,j,k,CP(I,K),C(i)

            CP(I,K) = CP(I,K) + AMAT(I,J,K)*(C(J) - ORIGIN(J,K))

CC           WRITE(0,930) i,j,k,CP(I,K)

	  end if
40        CONTINUE
30       CONTINUE
10      CONTINUE


c-----------------------------------------------------------------------
c Calculate magnetic field components of BP at field point CP,
c both in loop coordinates. The vector B is calculated for each loop by 
c subroutine CIRCLEB, which calls ELLIPINTS for the elliptic functions.
c  rcoil(k) = radius of loop k
c  currl(k) = current of loop k
c  cpk      = field point    CP=(XP,YP,ZP) of loop k in k's frame
c  Bpk      = magnetic field BP=(BX,BY,BZ) of loop k in k's frame
c-----------------------------------------------------------------------


CC      WRITE(0,*)'D3PFERRB  cpk / Bpk' 	! **********  DEBUG


      DO k=1, nloops
       if (kuse(k).ne.0) then	! do only loops with current
        DO i=1,3
         cpk(i) = CP(i,k)	! field coordinates in loop k's frame
        END DO

CC      WRITE(0,929) k,k,(cpk(m),m=1,3)


         CALL CIRCLEB (rcoil(k), curlps(k), cpk, Bpk)

CC      WRITE(0,929) k,k,(Bpk(m),m=1,3)

        DO i=1,3
         BP(i,k) = Bpk(i) 	! add loop k's B to total field BP 
        END DO
       end if
      END DO

c-----------------------------------------------------------------------
c Transform BP back to the main code's coordinate system.
c This code is from the old LOPLAB routine written by Todd Evans (2002).
C
C  GIVEN THE VECTOR BP AND ITS DERIVATIVES DBP IN THE (XP,YP,ZP) 
C  COORDINATE SYSTEM OF THE LOOP (OR ANY PLANE), TRANSFORM       
C  BP TO BX IN THE MAIN CODE'S (X,Y,Z) COORDINATE SYSTEM         
C  THE ORIGIN OF THE LOOP IS AT (X0,X0,Z0) AND                   
C  ITS ORIENTATION IS SPECIFIED BY (THETA,PHI) I.E. THE          
C  PERPENDICULAR TO THE PLANE OF THE LOOP (ZP-AXIS).             
C  THIS IS DONE FOR NLOOPS.                                      
c-----------------------------------------------------------------------

      DO  I=1,3
       DO  J=1,NLOOPS
        if (kuse(J).ne.0) then 	! do only loops with current

         BX(I,J) = AMAT(1,I,J)*BP(1,J) + AMAT(2,I,J)*BP(2,J)
     &           + AMAT(3,I,J)*BP(3,J)

        end if
       END DO
      END DO


CC      WRITE(0,*)'D3PFERRB  BX' 	! **********  DEBUG
CC      do j=1,nloops	! ******** DEBUG **********
CC       if(kuse(j).ne.0) then
CC        WRITE(0,929) j,j,(BX(m,j),m=1,3)
CC       end if
CC      end do

c-----------------------------------------------------------------------
c Sum magnetic fields of all current loops, transform to the main code's
c cylindrical coordinate system.
c-----------------------------------------------------------------------

      DO  J=1,NLOOPS
       if(kuse(J).ne.0) then 	! do only loops with current
        DO  I=1,3
         BXSUM(I) = BXSUM(I) + BX(I,J)
        END DO
       end if
      END DO

      bxf = BXSUM(1)
      byf = BXSUM(2)
      bzf = BXSUM(3)



C      WRITE (15,*) 'D3PFERRB:'
C      WRITE (15,*) 'nfc,nxtra,nlcur,nld:', nfc,nxtra,nlcur,nld
C      WRITE (15,*) 'curlsp:', curlps
C      WRITE (15,*) 'brsum:', brsum(1),brsum(2),brsum(3)
      
 929  format(1x,2i3,8f22.10)	
 930  format(1x,3i3,8f22.10)	

      RETURN

      END 	SUBROUTINE	D3PFERRB
c=======================================================================
