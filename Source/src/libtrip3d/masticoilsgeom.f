c Contents of file masticoilsgeom.f:
c
c	subroutine MASTIGEOM
c
c
c=======================================================================
 
      subroutine MASTIGEOM (kuse,nbands,nloops,nsegs,xs,dvs,curntw)
 
c=======================================================================
c                                      Started by: DM Orlov, 2009 jul 21
c                                      Last Modified:   DMO, 2009 jul 21
c                                      Last Modified:    MJS 2010 jul 24
c-----------------------------------------------------------------------
c MASTIGEOM reads geometry-defining information, then writes arrays of
c MAST I-coil geometry parameters in the form required by MJS's
c routines MASTICOILS, POLYGONB, BIOTLOOP.
c
c MAST I-coils can be defined with HELICAL SPIRAL PATHS on conical
c surfaces.
c
c The INTERNAL coils (I-coils) are represented as sets of polygon
c current loops, each comprised of straight-line segments. To offer the
c user considerable flexibility in defining arrays of coils, the total
c array is divided into toroidal "bands" or "rows".
c  A band is a set of IDENTICAL coils spaced circumferentially around
c   the torus on an imaginary CONICAL surface. Conical includes cylin-
c   drical and planar surfaces as limiting cases.
c  All coils in a band have same defining values of radii R, elevations
c   (heights) Z, coil TOROIDAL ANGLE SPAN phispan, coil-to-coil toroidal
c   angle spacing dphi, and HELICAL ARC TOROIDAL SPAN phiHel (deg). 
c   The 1st coil in a band is CENTERED at toroidal angle phicentr1.
c  An individual coil is a single polygon loop, constructed of
c   straight-line segments. Each coil has its own specified current.
c  All coils in a band array can be rotated toroidally by an optional 
c   'addangle' (deg), and its current distribution can be multiplied by
c   a constant current 'scale' factor (dimensionless).
c
c Even though the MAST internal coil geometries are fixed,
c  here we use the flexible-input P-coil geometry subroutine as  
c  a template. This makes the program simpler and maybe also easier to
c  maintain.
c
c  Maximum array dimensions set by INCLUDE MASTICOILS.i:
c   mxIbands, maximum # of "bands"
c   mxIloops, maximum # of "loops"
c   mxIsegs,  maximum # of "segments"
c   ncturns, number of electrical turns in an I-coil
c
c INPUTS:
c  MASTIGEOM gets parameters from the NAMELIST input file mast.in.
c  Icur(1:N,L) = currents in Iloops(L) coils (single turn, Amp)
c  Iadj(1:4,L) = addangle, arcmax, arcHmax = max. helical arc, (deg) 
c                current scale factor (dimless)   for row L
c
c OUTPUTS:
c  kuse(k,L), k=1...nloops(L) is integer flag to tell which loops are
c             active (1) or not (0) in band L.
c  nbands,    number of bands actually used in the set of coils.
c  nloops(L), number of identical loops actually used in band(L).
c  nsegs(k,l),  total number of segments used in loop(k) of band(L).
c     Though identical, need explicit value for each loop for POLYGONB.
c  xs(i,j,k,L), i=1,2,3 j=1...nsegs are (x,y,z) Cartesian position 
c     vectors of nsegs coil-defining points of loop(k) in band(L).
c  dvs(i,j,k,L), i=1,2,3,4 j=1...nsegs are Cartesian direction vectors  
c     (3 elements) and lengths (4th element) of nsegs straignt segments.
c     The start point of segment j is row j of xs, its end point is 
c     row j+1 of xs, except end point of segment nisegs is row 1 of xs.
c  curntw(l,k), total (Amp*turns) current of kth coil of lth band.
c     Number of turns in coil, ncturns, is set by INCLUDE "*.i' file
c
c Coordinate system: Geometry information is all written in right-
c   handed Cartesian coordinates, as used by POLYGONB, BIOTLOOP.
c
c This version orders array indices as (i,j,k,L). Previously (L,k,j,i).
c
c Units: SI (mks)
c Fully double precision.
c-----------------------------------------------------------------------
 
      IMPLICIT NONE
 
      INCLUDE 'masticoils.i'	! gives mxIbands, mxIloops, mxIsegs
				! niturns (integer)
 
      INTEGER    mxbands, mxloops, mxsegs	!local parameters
      PARAMETER (mxbands = mxIbands)		! to dimension arrays
      PARAMETER (mxloops = mxIloops)		! to dimension arrays
      PARAMETER (mxsegs  = mxIsegs) 		! to dimension arrays
 
      INTEGER npolsegs
      PARAMETER (npolsegs = 1)	!*** ONLY 1 STRAIGHT POLOIDAL SEGMENT
				!*** ALLOWED PER I-COIL END AT PRESENT
	! npolsegs does not apply to helical paths.
c-----------------------------------------------------------------------
 
 
c      INCLUDE 'constsi.i'
c	   Change by Wingen: replace INCLUDE by common block
      common /consts/  pi, twopi, cir, rtd, dtr
      common /Icurrents/ Iadj, Icur 
 
 
!  Common variables:
      real*8  pi, twopi, cir, rtd, dtr
 
!  Arguments:
      integer kuse(mxloops,mxbands)
      integer nbands, nloops(mxbands), nsegs(mxloops,mxbands)
      real*8  xs(3,mxsegs,mxloops,mxbands),dvs(4,mxsegs,mxloops,mxbands)
      real*8  curntw(mxloops,mxbands)
      intent(out):: kuse, nbands, nloops, nsegs, xs, dvs, curntw
 
!  Local and Namelist Variables:
      real*8  dp0, dp1
      logical useMASTIcoil, livecoil
      integer Iloops(mxbands)
      real*8  Igeom(8,mxbands), Iadj(4,mxbands)
      real*8  Icur(mxloops,mxbands)
      real*8  curnt(mxloops,mxbands)
      integer nbandsd, nloopsd(mxbands) ! duplicates for internal use
      real*8  curntd(mxloops,mxbands) 	! duplicates
      integer ntorsegs, nhelsegs
      real*8  phianchor, phispan, dphi, phiHel
      real*8  arcmax, segarc, arcHmax, segarcH
      real*8  phicent, phistart, phiend, phipoint
      real*8  phipoint1, phipoint2, phipoint3, phipoint4
      real*8  addangl, scalec
      real*8  R1, Z1, R2, Z2, R, Z, dRdphi, dZdphi
      integer i, j, k, L, m, n, j1, j2, j3, j4, jend

      REAL*8  CHECKR	! *******************

      SAVE

c-----------------------------------------------------------------------
c  Initialize
c-----------------------------------------------------------------------

      data arcmax /05.0d0/		! default I coil chord = 5 (deg)

      dp0 = 0.0d0			! double precision 0
      dp1 = 1.0d0

 
  ! Initialize namelist variables
      useMASTIcoil = .FALSE.
 
      Iloops = 0		! using fortran90 array assignment
      Igeom  = dp0
      Icur   = dp0
      Iadj   = dp0
 
  ! Initialize output arguments
      kuse   = 0
      nbands = 0
      nloops = 0
      nsegs  = 0
      xs     = dp0
      dvs    = dp0
      curnt  = dp0
 
! Set the FIXED geometry parameters of the 2 bands of MAST internal
! nonaxisymmetric coils:
 ! Using a fortran90 array initialization method
 ! This is a simple 4-sided model to roughly approximate MAST I coils

 ! Iloops(L) = N = number of coils (loops) in band (row) L
      Iloops(1) = 6 			! 1 upper band of 6 coils
      Iloops(2) = 12  			! 1 lower band of 6 coils

 ! Igeom(1:8,L)= R1,Z1,R2,Z2 (m), phicentr1,phispan,dphi, phiHel (deg) 
 !               for row L
      Igeom(1:8,1)= 
     &   (/ 1.452, +0.591, 1.339, +0.791,  15.0, 22.85, 60.0,  0.0 /)
      Igeom(1:8,2)= 
     &   (/ 1.339, -0.791, 1.452, -0.591,  15.0, 22.85, 30.0,  0.0 /)

 ! Iadj(1:4,L) = addangle, arcmax, arcHmax = max. helical arc, (deg) 
 !               current scale factor (dimless)   for row L
      Iadj(2,1) = arcmax
      Iadj(2,2) = arcmax

c-----------------------------------------------------------------------
c Open namelist input file 'mast.in' and read input data
c-----------------------------------------------------------------------
 
c changed by Wingen: Iadj and Icur are read elsewhere and included here 
c					 through common block
c      CALL MASTIREAD (useMASTIcoil, Iadj, Icur)
 
 
c      if(.not. useMASTIcoil)   RETURN	! don't need MAST I-coils
 
c-----------------------------------------------------------------------
c More initial stuff
c-----------------------------------------------------------------------
 
      nbands = 0
      DO L=1, mxbands			! find number of coil bands
       if(Iloops(L) .ne. 0) then
        nbands = nbands + 1
       end if
      END DO
 
c      write(81,*) ''		! write more details to file 81 = case2
c      write(81,*) nbands,' MAST I-coil bands'
 
      if(nbands .eq. 0)    RETURN   	! no I-coils are defined


      nloops = Iloops			! using f90 array assign
      curnt  = Icur

c adjust currents by scale factors from Iadj(4)
      DO  L=1,nbands  	! do all bands of I-coils
       DO  k=1,nloops(L)
        scalec = Iadj(4,L)
        curnt(k,L) = scalec*curnt(k,L)	! current for loop(k) of band(L)
       END DO
      END DO	! end DO  l=1,nbands
 
c set array of flags, kuse, to identify loops with non-zero current
      livecoil = .false.
      DO L=1,nbands
       DO k=1,nloops(L)
        if(curnt(k,L) .ne. dp0) then
	 kuse(k,L) = 1	! use loop(k) of band(L)
	 livecoil = .true.
	end if
       END DO
      END DO
 
  ! write input information to file 81 = case2
c      write(81,'(a,16i4)') ' nloops: ',(nloops(L),L=1,nbands)
c      write(81,'(a,16i4)') ' npoints:',(npoints(L),L=1,nbands)
C      WRITE(81,*) ''
c      write(81,*) 'Igeom:'
c      do L=1,nbands
c       write(81,'(8f10.3 )') (Igeom(k,L),k=1,8)
c      end do
c      write(81,*) 'Iadj:'
c      do L=1,nbands
c       write(81,'(8f10.3 )') (Iadj(k,L),k=1,4)
c      end do
c      write(81,*) 'I Currents (Amp*turns):'
c      do L=1,nbands
c       write(81,'(10f8.0)') (curnt(k,L)*niturns,k=1,nloops(L))
c      end do

C      WRITE(*,*) 'kuse:'
C      do L=1,nbands
C       WRITE(*,'(24i3)') (kuse(k,L),k=1,nloops(L))
C      end do
 
 
c=======================================================================
c            Fill Position Vector Array	
c=======================================================================
C To write points to fort.15, usually named OUT.DEBUG:
C        WRITE(15,*) ' '
C        WRITE(15,*) '***************  MAST I-COIL  *******************'
C        WRITE(15,*) 'nbands:', nbands
C        WRITE(15,*) 'nloops:', nloops
C        WRITE(15,*) ' '
 
c-----------------------------------------------------------------------
c Loop through all requested toroidal bands of I-coils
c-----------------------------------------------------------------------
c Toroidal angle arithmetic is done in degrees; 
c  convert to radians to evaluate trigonometric functions
c Here we consisently use indices i,j,k,L to identify (respectively)
c  Cartesian component, segment, loop or coil, band or row.  
 
      DO  L=1, nbands  		! do all bands
 
       R1        = Igeom(1,L)                
       Z1        = Igeom(2,L) 
       R2        = Igeom(3,L) 
       Z2        = Igeom(4,L) 
       phianchor = Igeom(5,L)
       phispan   = abs(Igeom(6,L))	! must be .ge. 0
       dphi      = Igeom(7,L)
       phiHel    = Igeom(8,L)
       addangl   = Iadj(1,L)
       arcmax    = abs(Iadj(2,L))	! must be .gt. 0
       arcHmax   = abs(Iadj(3,L))	! must be .ge. 0

       phianchor = phianchor + addangl

       IF(arcmax .eq. 0. .and. Iadj(4,L) .ne. 0.) THEN
        write(*,*)'MASTECGEOM: Iadj(2,*) = arcmax = 0 not allowed.'
	write(*,*)'MASTECGEOM: Stopping.'
	STOP
       END IF


c       write(81,*)''		! write geometry parameters to case2 file
c       write(81,*)'band',L
c       write(81,*) '       ',
c     &  'R1        Z1        R2        Z2 phianchor   phispan      dphi'
c     & ,'    phiHel'
c       write(81,'(8f10.3)') 
c     &         R1, Z1, R2, Z2, phianchor, phispan, dphi, phiHel


c     set sufficient number of toroidal segments for coils in band(L) 
       segarc   = phispan
       ntorsegs = 0
       n = 0
       DO
	 ntorsegs = ntorsegs + 1
	 segarc = abs(phispan/float(ntorsegs))
	 IF (segarc .le. arcmax)  EXIT		! satisfied arcmax

	 IF(2*ntorsegs + 2 .gt. mxsegs) THEN 	! too many segments
	  write(*,*)'MAST I-coil toroidal segments are too small'
c	  write(*,'( '' in band (row)'', I3)') L
c	  write(*,'(I4,'' total segments allowed per loop.'')') mxsegs
c          write(81,'(''segarc,  ntorsegs for this band:'', f10.3,i5)')
c     &               segarc, ntorsegs
	  write(*,*)'STOPPING'
	  STOP
	 END IF
       END DO				! end of 'WHILE'

C*       write(81,*) 'segarc, ntorsegs for this band'
C*       write(81,'(f10.3,i5)') segarc, ntorsegs
c       write(81,'(''segarc,  ntorsegs for this band:'', f10.3,i5)')
c     &            segarc, ntorsegs


c    set a sufficient number of HELICAL SEGMENTS for coils in band(L)
       segarcH  = phiHel	! phiHel and segarch can be + or -
       nhelsegs = 0
       n = 0
       DO				! this is a 'DO WHILE' loop
	 nhelsegs = nhelsegs + 1
	 segarcH = phiHel/float(nhelsegs)	!phiHel is + or -
	 IF (abs(segarcH) .le. arcHmax)  EXIT	! satisfied arcHmax

	 IF(2*nhelsegs + 2 .gt. mxsegs) THEN	! too many segments
	  write(*,*)'MAST I-coil helical segments are too small'
c	  write(*,'( '' in band (row)'', I3)') L
c	  write(*,'(I4,'' total segments allowed per loop.'')') mxsegs
c          write(81,'(''segarcH, nhelsegs for this band:'', f10.3,i5)')
c     &               segarcH, nhelsegs
	  write(*,*)'STOPPING'
	  STOP
	 END IF
       END DO				! end of 'WHILE'

c       write(81,'(''segarcH, nhelsegs for this band:'', f10.3,i5)')
c     &            segarcH, nhelsegs


c     check total number of loop segments vs. array dimension
       IF(2*(ntorsegs + nhelsegs) .gt. mxsegs) THEN
       	write(*,*)'There are too many total segments in MASTicoil loops'
c	write(*,'( '' in band (row)'', I3)') L
c	write(*,'(I4,'' total segments allowed per loop.'')') mxsegs
	write(*,*)'STOPPING'	  
	STOP
       END IF


    !-------------------------------------------------------------------
    ! Loop through all coils (loops) in band L of I-coils
    !-------------------------------------------------------------------
 
       DO 100 k=1,nloops(L) 	! do all loops in I-coil band L
 
        phicent = phianchor + (k - 1)*dphi	! center of loop k	
C*        phistart = phicent - phispan/2.0d0	! coil start angle
C*        phiend   = phicent + phispan/2.0d0	! coil end angle

c	! Define "Trapezoidal" coil corner angles (toroidal degrees)
c       ! phiHel = toroidal "advance" of slanted or helical side
        phipoint1 = phicent + (-phispan - phiHel)/2.0d0	
	phipoint2 = phicent + (+phispan - phiHel)/2.0d0 
        phipoint3 = phicent + (+phispan + phiHel)/2.0d0 
	phipoint4 = phicent + (-phispan + phiHel)/2.0d0 

c        write(81,*) ' L, k =', L, k
C*        write(81,*) '   phicent,  phistart,    phiend:'
C*        write(81,'(7f10.3)') phicent, phistart, phiend
c        write(81,*) 
c     &   '   phicent  phipoint1  phipoint2  phipoint3  phipoint4'
c        write(81,'(7f11.3)') 
c     &      phicent, phipoint1, phipoint2, phipoint3, phipoint4
 
C To write points to fort.15, usually named OUT.DEBUG:
C        WRITE(15,*) 'band =', L,'loop =', k
C        WRITE(15,*) 
C     &   ' phipoint(deg)           Rpoint(m)                 point #'
C        WRITE(15,*) 
C     &   ' Xpoint(m)               Ypoint(m)               Zpoint(m)'
 
 
c Go around loop k
 
        j = 0		! j is segment index counter for this loop
       
        ! 1st TOROIDAL conductor, from point1 to point2
	! It is always defined in direction of increasing phi
        ! USER SHOULD CHOOSE 1ST CONDUCTOR FOR DESIRED B VS CURRENT SIGN
        ! e.g. bottom conductor for outward B with positive coil current

        DO j1=1, ntorsegs+1 	! ntorsegs arcs, ntorsegs+1 points
         j = j+1 		! j is cumulative point counter
         phipoint = phipoint1 + (j1-1)*segarc	! increasing phi
         R = R1
	 Z = Z1
         xs(1,j,k,L) = R*cos(phipoint*dtr)	! Cartesian X
         xs(2,j,k,L) = R*sin(phipoint*dtr)	! Cartesian Y
         xs(3,j,k,L) = Z
	 CHECKR = SQRT(XS(1,j,k,L)**2 + XS(2,j,k,L)**2)
 
C         WRITE(*,'(3f10.3,I10)') phipoint, R, Z, j
C         WRITE(*,'(7f10.3)') xs(1,j,k,L), xs(2,j,k,L), xs(3,j,k,L),
C     &                      CHECKR
C         WRITE(*,*) ' '
        END DO			! end of DO j1=1, ...
 
C         WRITE(*,*) '----------------------------------------'

 
        ! 1st HELICAL conductor, from point2 toward point3
        ! If only 1 helical segment is needed, then jump directly
	!  to point 3; otherwise calculate points with nhelsegs to 
	!  approximate a curved surface

        IF (abs(nhelsegs) .ge. 1) THEN
         dRdphi = (R2 - R1)/(phipoint3 - phipoint2)
 	 dZdphi = (Z2 - Z1)/(phipoint3 - phipoint2)
         DO j2=1, nhelsegs-1	! nhelsegs sections, nhelsegs-1 pts
          j = j+1		! j is cumulative point counter
          phipoint = phipoint2 + j2*segarcH	! step toward point3
					! segarcH has sign of phiHel
          R = R1 + dRdphi*j2*segarch
	  Z = Z1 + dZdphi*j2*segarch
          xs(1,j,k,L) = R*cos(phipoint*dtr)
          xs(2,j,k,L) = R*sin(phipoint*dtr)
          xs(3,j,k,L) = Z
	  CHECKR = SQRT(XS(1,j,k,L)**2 + XS(2,j,k,L)**2)

C         WRITE(*,'(3f10.3,I10)') phipoint, R, Z, j
C         WRITE(*,'(7f10.3)') xs(1,j,k,L), xs(2,j,k,L), xs(3,j,k,L),
C     &                      CHECKR
C         WRITE(*,*) ' '
         END DO			!end of DO j2=1, ...
        END IF

C         WRITE(*,*) '----------------------------------------'

 
        ! 2nd TOROIDAL conductor, from point3 to point4
	! It is always defined in direction of decreasing phi
         
         DO j3=1, ntorsegs+1 	! ntorsegs arcs, ntorsegs+1 points
         j = j+1 		! j is cumulative point counter
         phipoint = phipoint3 - (j3-1)*segarc	! decreasing phi
         R = R2
	 Z = Z2
         xs(1,j,k,L) = R*cos(phipoint*dtr)	! Cartesian X
         xs(2,j,k,L) = R*sin(phipoint*dtr)	! Cartesian Y
         xs(3,j,k,L) = Z
	 CHECKR = SQRT(XS(1,j,k,L)**2 + XS(2,j,k,L)**2)

C         WRITE(*,'(3f10.3,I10)') phipoint, R, Z, j
C         WRITE(*,'(7f10.3)') xs(1,j,k,L), xs(2,j,k,L), xs(3,j,k,L),
C     &                      CHECKR
C         WRITE(*,*) ' '
        END DO			! end of DO j3=1, ...

C         WRITE(*,*) '----------------------------------------'
	 
 
        ! 2nd HELICAL conductor, from point4 toward point1
        ! If only 1 helical segment is needed, then jump directly
	!  to point 3; otherwise calculate points with nhelsegs to 
	!  approximate a curved surface

        IF (abs(nhelsegs) .ge. 1) THEN
         dRdphi = (R1 - R2)/(phipoint1 - phipoint4)
 	 dZdphi = (Z1 - Z2)/(phipoint1 - phipoint4)
         DO j4=1, nhelsegs-1	! nhelsegs sections, nhelsegs-1 pts
          j = j+1		! j is cumulative point counter
          phipoint = phipoint4 - j4*segarch	! step toward point3
          R = R2 - dRdphi*j4*segarch
	  Z = Z2 - dZdphi*j4*segarch
          xs(1,j,k,L) = R*cos(phipoint*dtr)
          xs(2,j,k,L) = R*sin(phipoint*dtr)
          xs(3,j,k,L) = Z
	 CHECKR = SQRT(XS(1,j,k,L)**2 + XS(2,j,k,L)**2)

C         WRITE(*,'(3f10.3,I10)') phipoint, R, Z, j
C         WRITE(*,'(7f10.3)') xs(1,j,k,L), xs(2,j,k,L), xs(3,j,k,L),
C     &                      CHECKR
C         WRITE(*,*) ' '
         END DO			!end of DO j4=1, ...
        END IF

C         WRITE(*,*) '****************************************'


      nsegs(k,L) = j	! number of segments in loop k of band L
 
C         WRITE(15,*)' End of loop',k,' nsegs =', nsegs(k,L)
 
  100 CONTINUE		! end DO k=1,nloops(L) ! do all loops in band L
 
      END DO	! end DO  L=1,nbands  	! do all bands of I-coils
 
 
c=======================================================================
c            Fill Direction Vector Array	
c=======================================================================
 
      DO  l=1,nbands  	! do all bands of I-coils
 
      DO 200 k = 1, nloops(L)
       DO 202 j = 1, nsegs(k,L)
          
        IF (j.ne.nsegs(k,L)) THEN  
         jend = j + 1		! end point for all segments but last
        ELSE
         jend = 1 		! end point for last segment
        ENDIF
 
        ! segment length vector:
        dvs(1,j,k,L) = xs(1,jend,k,L) - xs(1,j,k,L)
        dvs(2,j,k,L) = xs(2,jend,k,L) - xs(2,j,k,L)
        dvs(3,j,k,L) = xs(3,jend,k,L) - xs(3,j,k,L)
 
        ! segment length:
        dvs(4,j,k,L) = 
     &       sqrt(dvs(1,j,k,L)**2 + dvs(2,j,k,L)**2 + dvs(3,j,k,L)**2)
 
        ! segment direction unit vector:
        DO 203 i=1,3
         dvs(i,j,k,L) = dvs(i,j,k,L)/dvs(4,j,k,L)
  203   CONTINUE
 
  202  CONTINUE
  200 CONTINUE
 
      END DO	! end DO  L=1,nbands  	! do all bands of I-coils
 
 

!-----------------------------------------------------------------------
!       Write data file for plotgeom
!-----------------------------------------------------------------------

c      write(93,'(I3,A)'),nbands,'   bands. I coils'
c      DO  L=1,nbands
c	write(93,'(I3,A)'),nloops(L),'   loops'
c	DO  k = 1, nloops(L)
c	  write(93,'(I3,A)'),nsegs(k,L),'   segs'
c	  DO  j = 1, nsegs(k,L)
c	    write(93,'(3F10.4)'),xs(1,j,k,L),xs(2,j,k,L),xs(3,j,k,L)
c	  enddo
c	enddo
c      enddo 

c-----------------------------------------------------------------------
c Dummy argument variables returned to the program calling subroutine 
c MASTIGEOM become unavailable for execution by following procedures 
c (like ENTRY ...). This restriction does not apply to all other 
c variables, which are still available after the return.
c Therefore, before returning, here we define duplicates of those dummy
c variables that will be needed in the rest of this subroutine.
c
c We also multiply here the input current array by number of turns.
 
      nbandsd = nbands
      nloopsd = nloops
      curntd  = curnt

      curntw  = curnt*niturns		! multiply by number of turns
 
      
C      WRITE(*,*) 'LEAVING MASTIGEOM'	! for DEBUG
 
      RETURN
      
      END 	SUBROUTINE	MASTIGEOM
c=======================================================================
