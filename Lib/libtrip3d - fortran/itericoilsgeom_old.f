

c Contents of file itereccoils.f:
c
c	subroutine ITERIGEOM
c
c
c=======================================================================

      subroutine ITERIGEOM (kuse,nbands,nloops,nsegs,xs,dvs,curntw)

c=======================================================================
c                                   Started by: MJ Schaffer, 2006 dec 15
c                                   Last Modified: D. Orlov, 2010 jul 28
c-----------------------------------------------------------------------
c ITERIGEOM reads geometry-defining information, then writes arrays of
c ITER I-coil geometry parameters in the form required by MJS's 
c routines ITERICOILS, POLYGONB, BIOTLOOP.
c
c The Internal RMP (I-coils) are represented as sets of polygon
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
c Even though the ITER I coil geometries are fixed (per 2010 ITER 
c   Design), here we use the flexible-input P-coil geometry
c  subroutine as a template. This makes the program simpler and maybe
c  also easier to maintain.
c
c  Maximum array dimensions set by INCLUDE ITERICOILS.i:
c   mxIbands, maximum # of "bands"
c   mxIloops, maximum # of "loops"
c   mxIegs,  maximum # of "segments"
c   niturns, number of electrical turns in an I-coil
c
c INPUTS:
c  ITERIGEOM gets parameters from the NAMELIST input file iter.in.
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
c     Number of turns in coil, niturns, is set by INCLUDE "*.i' file
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
 
      INCLUDE 'itericoils.i'	! gives mxIbands, mxIloops, mxIsegs
				! niturns (integer)
	  
      INTEGER    mxbands, mxloops, mxsegs	!local parameters
      PARAMETER (mxbands = mxIbands)		! to dimension arrays
      PARAMETER (mxloops = mxIloops)		! to dimension arrays
      PARAMETER (mxsegs  = mxIsegs) 		! to dimension arrays

 
      INTEGER npolsegs
      PARAMETER (npolsegs = 1)	!*** ONLY 1 STRAIGHT POLOIDAL SEGMENT
				!*** ALLOWED PER I-COIL SIDE AT PRESENT
c-----------------------------------------------------------------------


c      INCLUDE 'constsi.i'
c	   Change by Wingen: replace INCLUDE by common block
      common /consts/  pi, twopi, cir, rtd, dtr
      common /currents/ Iadj, Icur 
 
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
      logical useITERIcoil, livecoil
      integer Iloops(mxbands)
      real*8  Igeom(7,mxbands), Iadj(3,mxbands)
      real*8  Icur(mxloops,mxbands)
      real*8  curnt(mxloops,mxbands)
      integer nbandsd, nloopsd(mxbands) ! duplicates for internal use
      real*8  curntd(mxloops,mxbands) 	! duplicates
      integer ntorsegs
      real*8  phianchor, phispan, dphi
      real*8  arcmax, segarc
      real*8  phicent, phistart, phiend, phipoint
      real*8  addangl, scalec
      real*8  R1, Z1, R2, Z2, R, Z
      integer i, j, k, l, m, n, j1, j2, j3, j4, jend
 
      SAVE
 
c-----------------------------------------------------------------------
c  Initialize
c-----------------------------------------------------------------------
 
      data arcmax /15.0d0/		! default 15 (deg)
 
      dp0 = 0.0d0			! double precision 0
      dp1 = 1.0d0			!    "       "     1

 ! Initialize geometry and namelist variables
      useITERIcoil = .FALSE.
 
      Iloops = 0		! using fortran90 array assignment
      Igeom  = dp0
c      Icur   = dp0
c      Iadj   = dp0
 
 ! Initialize output arguments
      kuse   = 0
      nbands = 0
      nloops = 0
      nsegs  = 0
      xs     = dp0
      dvs    = dp0
      curnt  = dp0
 
! Set the FIXED geometry parameters of the 3 bands of ITER
! 2010 Design Internal RMP Coils:
 ! Using a fortran90 array initialization method
 
 ! Iloops(1:L) = number of coils (loops) in band (row) L  
      Iloops(1:3) = (/9, 9, 9/) 	! 3 rows, Upper, Midplane, Lower 
 
 ! Igeom(1:7,L)= R1,Z1,R2,Z2, phicentr1, phispan, dphi of band L (m,deg)

      Igeom(1:7,1)=(/7.735, 3.380, 8.262, 2.626, 0., 28.5, 40./)
C     9 top loops

      Igeom(1:7,2)=(/8.618, 1.790, 8.661, -0.550, -3.3, 20.2, 40./)
C     9 middle loops

      Igeom(1:7,3)=(/8.230, -1.546, 7.771, -2.381, 0., 30.0, 40./)
C     9 bottom loops      
 
 ! Iadj(1:3:L) = addangl, arcmax (deg), scalec (dimless)
c      Iadj(2,1) = arcmax
 
c-----------------------------------------------------------------------
c Open namelist input file 'iter.in' and read input data relevant to
c ITER Error Correction coils
c-----------------------------------------------------------------------
c Namelist file iter.in was read by subroutine ITERCOILSREAD.
c We get the needed EC error information from entry point ITERECREAD.

c changed by Wingen: Iadj and Icur are read elsewhere and included here 
c					 through common block
 
c      CALL ITERIREAD (useITERIcoil, Iadj, Icur)
 
c      if(.not. useITERIcoil)   RETURN	! don't need ITER I-coils
 
c-----------------------------------------------------------------------
c More initial stuff
c-----------------------------------------------------------------------
 
      nbands = 0
      DO L=1, mxbands			! find number of coil bands
       if(Iloops(L) .ne. 0) then
        nbands = nbands + 1
       end if
      END DO
 
      if(nbands .eq. 0)    RETURN   	! no I-coils are defined
 
 
      nloops = Iloops			! f90 array assign
      curnt  = Icur
  
c adjust currents by scale factors from Iajd(3)
      DO  L=1,nbands  	! do all bands of I-coils
       DO  k=1,nloops(L)
        scalec = Iadj(3,L)
        curnt(k,L) = scalec*curnt(k,L)	! current of loop(k) of band(L)
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

c=======================================================================
c            Fill Position Vector Array	
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
       phispan   = abs(Igeom(6,L))	! must be .gt. 0
       dphi      = Igeom(7,L)
       addangl   = Iadj(1,L)
       arcmax    = abs(Iadj(2,L))	! must be .gt. 0

       phianchor = phianchor + addangl
       
       IF(arcmax .eq. 0. .and. Iadj(3,L) .ne. 0.) THEN
        write(*,*)'ITERIGEOM: Iadj(2,*) = arcmax = 0 not allowed.'
	write(*,*)'ITERIGEOM: Stopping.'
	STOP
       END IF
 
c     set sufficient number of toroidal segments for coils in band(L) 
       segarc   = phispan
       ntorsegs = 0
       n = 0
       DO				! this is a 'DO WHILE' loop
	 ntorsegs = ntorsegs + 1
	 segarc = phispan/ntorsegs
	 IF (segarc .le. arcmax)  EXIT		! satisfied segarc
 
	 IF(ntorsegs .gt. 100) THEN		! stop runaway loop
	  write(*,*)'ITERIGEOM: Coil segment arc size setup ran away.'
	  write(*,*)'ITERIGEOM: Stopping.'
	  STOP
	 END IF
       END DO				! end of 'WHILE'
 
    !-------------------------------------------------------------------
    ! Loop through all coils (loops) in band L of I-coils
    !-------------------------------------------------------------------
 
       DO 100 k=1,nloops(L) 	! do all loops in I-coil band L
 
        phicent = phianchor + (k - 1)*dphi	! center of loop k	
        phistart = phicent - phispan/2.0d0	! coil start angle
        phiend   = phicent + phispan/2.0d0	! coil end angle
 
c Go around loop k
 
        j = 0		! j is segment index counter for this loop
       
        ! 1st toroidal conductor, increasing coordinate phi 
        ! USER SHOULD CHOOSE 1ST CONDUCTOR FOR DESIRED B VS CURRENT SIGN
        ! e.g. bottom conductor for outward B with positive coil current
 
        DO 20 j1=1,ntorsegs+1   ! ntorsegs arcs, ntorsegs+1 points
         j = j+1 		! j is cumulative segment index counter
         R = R1
	 Z = Z1
         phipoint = phistart + (j1-1)*segarc
         xs(1,j,k,L) = R*cos(phipoint*dtr)
         xs(2,j,k,L) = R*sin(phipoint*dtr)
         xs(3,j,k,L) = Z
 
   20   CONTINUE
 
        ! 1st poloidal conductor, only 1 segment, no new point
C ! NOT IMPLEMENTED, SKIP FOR NOW
C        DO 30 j2=1,npolsegs-1	! npolsegs sections, npolsegs-1 pts
C         j = j+1 		! j is cumulative segment index counter
C         phipoint = phiend
C         xs(1,j,k,L) = Ru*cos(phipoint*dtr)
C         xs(2,j,k,L) = Ru*sin(phipoint*dtr)
C         xs(3,j,k,L) = Zu - j2*segz
C         
   30   CONTINUE
 
         ! 2nd toroidal conductor, decreasing coordinate phi
         
         DO 40 j3=1,ntorsegs+1  ! ntorsegs arcs, ntorsegs+1 points
         j = j+1 		! j is cumulative segment index counter
         R = R2
	 Z = Z2
         phipoint = phiend - (j3-1)*segarc
         xs(1,j,k,L) = R*cos(phipoint*dtr)
         xs(2,j,k,L) = R*sin(phipoint*dtr)
         xs(3,j,k,L) = Z
 
   40   CONTINUE
 
        ! 2nd poloidal conductor, only 1 segment, no new point
C ! NOT IMPLEMENTED, SKIP FOR NOW
C        DO 50 j4=1,npolsegs-1 	! npolsegs sections, npolsegs-1 pts
C         j = j+1 		! j is cumulative segment index counter
C         phipoint = phistart
C         xs(1,j,k,L) = Ru*cos(phipoint*dtr)
C         xs(2,j,k,L) = Ru*sin(phipoint*dtr)
C         xs(3,j,k,L) = Zu - j4*segz
C         
   50   CONTINUE
 
      nsegs(k,L) = j	! number of segments in loop k of band L
 
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
 
c-----------------------------------------------------------------------
c Dummy argument variables returned to the program calling subroutine 
c ITERIGEOM become unavailable for execution by following procedures 
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

      RETURN
 
      END 	SUBROUTINE	ITERIGEOM
c=======================================================================
