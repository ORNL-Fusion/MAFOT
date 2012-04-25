

c Contents of file itereccoils.f:
c
c	subroutine ITERICOILS
c	  entry BITERICOILS
c	subroutine ITERIGEOM
c         entry ITERIGEOMCASE
c
c
c=======================================================================
c
c      subroutine ITERICOILS
c
c-----------------------------------------------------------------------
c                                            By: D. Orlov    2010 jul 28
c                                 Last Modified: D. Orlov    2010 jul 28
c-----------------------------------------------------------------------
c ITERICOILS defines the geometry of and calculates the magnetic field
c from one or more specified models of the ITER nonaxisymmetric or
c Internal RMP coils (I-coils).
c This subroutine handles all the tasks needed to calculate this field.
c
c This subroutine has several main purposes:
c  1) Call routines to set up  models to be used.
c  2) Set flags for components to use.
c  3) Upon separate entry at BITERICOILS, call subroutines to calculate
c     the magnetic field at a point (x,y,z) made by the coils, using the
c     generic polynomial loop B-field evaluator POLYGONB.
c
c  Maximum array dimensions set by INCLUDE ITERECCOILS.i:
c   mxIbands, maximum # of "bands"
c   mxIloops, maximum # of "loops"
c   mxIsegs,  maximum # of "segments"
c
c This version orders array indices as (i,j,k,L). Previously (L,k,j,i).
c   i for x,y or z; j for segment; k for loop; L for band or row
c 
c Units: SI (mks)
c Fully double precision.
c-----------------------------------------------------------------------
c 
c     IMPLICIT NONE
c 
c     INCLUDE 'itericoils.i'	! gives mxIbands, mxIloops, mxIsegs
c				! niturns (integer)
c	  
c      INTEGER    mxbands, mxloops, mxsegs	!local parameters
c      PARAMETER (mxbands = mxIbands)		! to dimension arrays
c      PARAMETER (mxloops = mxIloops)		! to dimension arrays
c      PARAMETER (mxsegs  = mxIsegs) 		! to dimension arrays
c-----------------------------------------------------------------------
c 
c      INCLUDE 'constsi.i'
c      common /consts/  pi, twopi, cir, rtd, dtr
c
c      INCLUDE 'flags2i.i'
c      common /flags2/ iPROBE, iSURFMN, iTRIP3D,
c     &                iD3Dcoils, iNSTXcoils, iJETcoils, iITERcoils, 
c     &                iFDFcoils, iASDEXcoils, iMASTcoils
c  
c 
!  Common variables:
c      real*8  pi, twopi, cir, rtd, dtr
c      integer iPROBE, iSURFMN, iTRIP3D
c      integer iD3Dcoils, iNSTXcoils, iJETcoils, iITERcoils, iFDFcoils
c      integer iASDEXcoils, iMASTcoils
c 
!  Arguments at ENTRY BITERECCOILS:
c      real*8  x, y, z, bx, by, bz
c      intent(in)  ::  x,  y,  z
c      intent(out) :: bx, by, bz
c
!  Local and Namelist Variables:
c      logical idone, iIcoil
c     integer kuse(mxloops,mxbands)
c      integer nbands, nloops(mxbands), nsegs(mxloops,mxbands) 
c      real*8  xs(3,mxsegs,mxloops,mxbands),dvs(4,mxsegs,mxloops,mxbands)
c      real*8  curnt(mxloops,mxbands)
c      real*8  bbx, bby, bbz
c      integer i,j,k,l
c
c      SAVE
c
c-----------------------------------------------------------------------
c Set some default values
c-----------------------------------------------------------------------
c Variables initialized by data statement acquire the SAVE attribute.
c
c      idone = .false.
c
c-----------------------------------------------------------------------
c Set up geometries of EC-coil filaments that will be used; do only once
c-----------------------------------------------------------------------
c 
c     if(.not. idone) then
c 
c       CALL  ITERIGEOM (kuse,nbands,nloops,nsegs,xs,dvs,curnt)
c 
c       idone = .true.		! to avoid repeating geometry setup
c 
c      end if
c 
c-----------------------------------------------------------------------
c Set flags
c-----------------------------------------------------------------------
c global flag iITERcoils is initialized at 0 in ITER_MODEL
c
c      iIcoil = .false.
c      DO L=1, nbands
c       DO k=1, nloops(L)
c        if(kuse(k,L) .ne. 0) then
c         iIcoil = .true.	! at least 1 I-coil loop has current
c        iITERcoils = 1		! at least 1 ITER coil model is used
c         GO TO 5		! no need to test further
c        end if
c 
c       END DO
c      END DO
c    5 CONTINUE
c
c      RETURN		! End of ITER I Coils Model setup
c
c
c
c
c
c=======================================================================
c Calculate magnetic field for all requested ITER I-coils.
c
c ENTRY BITERICOILS
c  Receives the point (x,y,z) at which the field is to be calculated;
c  Calls subroutine POLYGONB to evaluate the field at (x,y,z);
c  Returns the Cartesian magnetic field components (bx,by,bz).
c
c POLYGONB calculates the magnetic field from a generic group of polygon
c  loops using subroutine BIOTLOOP.
c-----------------------------------------------------------------------
c
c
c
c      ENTRY BITERICOILS (x, y, z, bx, by, bz)
c
C      WRITE(*,*) 'ENTERED BITERICOILS'
c
c      bx = 0.0d0		! initialize magnetic field components
c      by = 0.0d0		! that will be calculated and returned
c      bz = 0.0d0
c
c      IF (.not. iIcoil) THEN	! return 0 if no I-coil is active
c
c       RETURN
c
c      ELSE			! calculate B-field from I-coil model
c				! Need to loop through nbands bands
c				! because POLYGONB only accepts a
c				! single group of loops.
c       DO L=1,nbands
c      using 'NEW' METHOD with order of indices reversed to i,j,k,L
c
c        CALL POLYGONB(mxloops, mxsegs, nloops(L), nsegs(:,L), kuse(:,L),
c     &                xs(:,:,:,L), dvs(:,:,:,L), curnt(:,L), 
c     &                x, y, z, bbx, bby, bbz)
c
c        bx = bx + bbx		! sum each band's contribution
c        by = by + bby
c        bz = bz + bbz
c       END DO
c
c      END IF			! end of IF (.not. iIcoil)
c
c      RETURN		! END OF ENTRY BITERICOILS
c
c
c-----------------------------------------------------------------------
c Format statements
c-----------------------------------------------------------------------
c 1003 format(1x,30i3)
c 1010 format(1x,8f10.0)
c 1014 format(1x,8f9.4)
c 
c 
c      END 	SUBROUTINE 		ITERICOILS
c=======================================================================
 
 
 
 
 
 
 
 
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
 
c      write(81,*) ''		! write more details to file 81 = case2
c      write(81,*) nbands,' ITER I-coil bands'
 
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

  ! write input information to file 81 = case2
c       write(81,'(a,16i4)') ' nloops: ',(nloops(L),L=1,nbands)
c      write(81,'(a,16i4)') ' npoints:',(npoints(L),L=1,nbands)
C      WRITE(81,*) ''
c      write(81,*) 'Igeom:'
c      do L=1,nbands
c       write(81,'(8f10.3 )') (Igeom(k,L),k=1,7)
c      end do
c      write(81,*) 'Iadj:'
c      do L=1,nbands
c       write(81,'(8f10.3 )') (Iadj(k,L),k=1,3)
c      end do
c      WRITE(81,*) 'I Currents (Amp):'
c      do L=1,nbands
c       WRITE(81,'(10f8.0)') (curnt(k,L),k=1,nloops(L))
c      end do

C      WRITE(*,*) 'kuse:'
C      do L=1,nbands
C       WRITE(*,'(24i3)') (kuse(k,L),k=1,nloops(L))
C      end do
 
 
c=======================================================================
c            Fill Position Vector Array	
c-----------------------------------------------------------------------
C To write points to fort.15, usually named OUT.DEBUG:
C        WRITE(15,*) ' '
C        WRITE(15,*) '**************   ITER I-COIL   *****************'
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
 
 
c       write(81,*)''
c       write(81,*)'band',L
c       write(81,*) '       ',
c     & 'R1,       Z1,       R2,       Z2, phianchor, phispan,    dphi:'
c       write(81,'(7f10.3)') 
c     &         R1, Z1, R2, Z2, phianchor, phispan, dphi
 
c     set sufficient number of toroidal segments for coils in band(L) 
       segarc   = phispan
       ntorsegs = 0
       n = 0
       DO				! this is a 'DO WHILE' loop
	 ntorsegs = ntorsegs + 1
	 segarc = phispan/ntorsegs
	 IF (segarc .le. arcmax)  EXIT		! satisfied segarc
 
	 IF(ntorsegs .gt. 100) THEN		! stop runaway loop
      write(*,'(f10.3,i5,f10.3)') segarc, ntorsegs, arcmax
	  write(*,*)'ITERIGEOM: Coil segment arc size setup ran away.'
	  write(*,*)'ITERIGEOM: Stopping.'
	  STOP
	 END IF
       END DO				! end of 'WHILE'
 
c       write(81,*) 'segarc, ntorsegs for this band'
c       write(81,'(f10.3,i5)') segarc, ntorsegs
 
 
    !-------------------------------------------------------------------
    ! Loop through all coils (loops) in band L of I-coils
    !-------------------------------------------------------------------
 
       DO 100 k=1,nloops(L) 	! do all loops in I-coil band L
 
        phicent = phianchor + (k - 1)*dphi	! center of loop k	
        phistart = phicent - phispan/2.0d0	! coil start angle
        phiend   = phicent + phispan/2.0d0	! coil end angle
 
c        WRITE(81,*) ' L, k =', L, k
c        WRITE(81,*) '   phicent,  phistart,    phiend:'
c        WRITE(81,'(7f10.3)') phicent, phistart, phiend
 
C To write points to fort.15, usually named OUT.DEBUG:
C        WRITE(15,*) 'band =', L,'loop =', k
C        WRITE(15,*) 
C     &   ' phipoint(deg)           Rpoint(m)                 point #'
C        WRITE(15,*) 
C     &   ' Xpoint(m)               Ypoint(m)               Zpoint(m)'
 
 
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
 
C         WRITE(15,*) phipoint, R, j
C         WRITE(15,*) xs(1,j,k,L),xs(2,j,k,L),xs(3,j,k,L)
C         WRITE(15,*) ' '
 
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
 
C         WRITE(15,*) phipoint, R, j
C         WRITE(15,*) xs(1,j,k,L),xs(2,j,k,L),xs(3,j,k,L)
C         WRITE(15,*) ' '
 
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
c     enddo

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

      
C      WRITE(*,*) 'LEAVING ITERIGEOM'		! for DEBUG
 
      RETURN
      
 
 
 
 
c=======================================================================
c Write some information about this model to file 80 = CASE 
c-----------------------------------------------------------------------

c      ENTRY ITERIGEOMCASE
 
c      if(.not. livecoil) RETURN		! no current-carrying I-coil
 
C      WRITE(*,*) 'ENTRY ITERIGEOMCASE'
C      WRITE(*,*) 'nbands'
C      WRITE(*,*) 'nloops', nloops
C      WRITE(*,*) 'nbands', nbands
 
c Write I-coil information to 'CASE' file:
c      write(80,*) ''
c      write(80,'('' I-coil information for'',i3,'' bands:'')') nbandsd
c      write(80,'(a7,16i4)') ' nloops:',(nloopsd(l),l=1,nbandsd)
c      write(80,*)' I-coil parameters:',
c     &           ' R1, Z1, R2, Z2, phicent1, phispan, dphi, addangl'
c      do L=1, nbandsd
c       write(80,1006) (Igeom(k,L),k=1,7), ECadj(1,L)
c      end do
 
c      write(80,*)' I-coil parameters',
c     &           ' R1, Z1, R2, Z2, coil1phi, torwdth, cntr-cntr:'
c      do L=1, nbandsd
c       write(80,1006) (Igeom(k,L),k=1,4), Igeom(5,L)+Iadj(1,L),
c     &                (Igeom(k,L),k=6,7)
c      end do
c 
c      write(80,*) ' I-coil Currents (A):'
c      do L=1, nbandsd
c       write(80,1005) (curntd(k,L),k=1,nloopsd(L))
c      end do
c 
c      RETURN
 
c 1003 format(1x,30i3)
c 1005 format(9f9.0)
c 1006 format(9f9.3)
 
 
      END 	SUBROUTINE	ITERIGEOM
c=======================================================================
