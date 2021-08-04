c Contents of this file:
c
c	subroutine TCABRCPGEOM
c
c
c=======================================================================
 
      subroutine TCABRCPGEOM (kuse,nbands,nloops,nsegs,xs,dvs,curntw)
 
c=======================================================================
c                                   Started by: MJ Schaffer, 2006 dec 15
c                                      Last Modified:   MJS, 2010 jul 23
c-----------------------------------------------------------------------
c TCABRCPGEOM reads geometry-defining information, then writes arrays of
c TCABR CP-coil geometry parameters in the form required by MJS's 
c routines TCABRCPCOILS, POLYGONB, BIOTLOOP.
c
c The Center-post (CP-coils) are represented as sets of polygon
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
c Even though the TCABR Center-post coil geometries are 
c  fixed, here we use the flexible-input P-coil geometry subroutine as  
c  a template. This makes the program simpler and maybe also easier to
c  maintain.
c
c  Maximum array dimensions set by INCLUDE TCABRCPCOILS.i:
c   mxECbands, maximum # of "bands"
c   mxECloops, maximum # of "loops"
c   mxECsegs,  maximum # of "segments"
c   ncturns, number of electrical turns in an CP-coil
c
c INPUTS:
c  TCABRCPGEOM gets parameters from the NAMELIST input file tcabr.in.
c  ECcur(1:N,L) = currents in ECloops(L) coils (single turn, Amp)
c  ECadj(1:4,L) = addangle, arcmax, arcHmax = max. helical arc, (deg) 
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
 
      INCLUDE 'tcabrcpcoils.i'	! gives mxECbands, mxECloops, mxECsegs,
				! ncturns (integer)
 
      INTEGER    mxbands, mxloops, mxsegs	!local parameters
      PARAMETER (mxbands = mxECbands)		! to dimension arrays
      PARAMETER (mxloops = mxECloops)		! to dimension arrays
      PARAMETER (mxsegs  = mxECsegs) 		! to dimension arrays

 
      INTEGER npolsegs
      PARAMETER (npolsegs = 1)	!*** ONLY 1 STRAIGHT POLOIDAL SEGMENT
				!*** ALLOWED PER CP-COIL END AT PRESENT
c-----------------------------------------------------------------------
 
 
c	   Change by Wingen: replace INCLUDE by common block 
c      INCLUDE 'constsi.i'
      common /consts/  pi, twopi, cir, rtd, dtr
      common /ECcurrents/ ECadj, ECcur 
 
 
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
      logical useTCABRCPcoil, livecoil
      integer ECloops(mxbands)
      real*8  ECgeom(7,mxbands), ECadj(3,mxbands)
      real*8  ECcur(mxloops,mxbands)
      real*8  curnt(mxloops,mxbands)
      integer nbandsd, nloopsd(mxbands) ! duplicates for internal use
      real*8  curntd(mxloops,mxbands) 	! duplicates
      integer ntorsegs
      real*8  phianchor, phispan, dphi
      real*8  arcmax, segarc
      real*8  phicent, phistart, phiend, phipoint
      real*8  addangl, scalec
      real*8  R1, Z1, R2, Z2, R, Z
      integer i, j, k, L, m, n, j1, j2, j3, j4, jend
 
      SAVE
 
c-----------------------------------------------------------------------
c  Initialize
c-----------------------------------------------------------------------
 
      data arcmax /15.0d0/		! default = 15 (deg)
 
      dp0 = 0.0d0			! double precision 0
      dp1 = 1.0d0


 ! Initialize geometry and namelist variables
      useTCABRCPcoil = .FALSE.
 
      ECloops = 0		! using fortran90 array assignment
      ECgeom  = dp0
c      ECcur   = dp0
c      ECadj   = dp0
 
 ! Initialize output arguments
      kuse   = 0
      nbands = 0
      nloops = 0
      nsegs  = 0
      xs     = dp0
      dvs    = dp0
      curnt  = dp0
 
! Set the FIXED geometry parameters of the 1 band of TCABR external
! Error Correction Coils:
 ! Using a fortran90 array initialization method
 ! This is a simple 4-sided model to roughly approximate TCABR CP coils

 ! ECloops(1:L) = number of coils (loops) in band (row) L  
      ECloops(1:3) = (/18,18,18/) 	! 3 bands of 18 coils

! Lower CP Coil
      ECgeom(1:7,1) = (/ 0.417, -0.027, 0.417, -0.063, 20., 16., 20. /) 
! Mid(Central) CP coil
      ECgeom(1:7,2) = (/ 0.417, 0.018, 0.417, -0.018, 20., 16., 20. /)
! Upper CP Coil
      ECgeom(1:7,3) = (/ 0.417, 0.063, 0.417, 0.027, 20., 16., 20. /)

! ECadj(1:3,L) = addangle, arcmax (deg), scalec (dimless)
c      ECadj(2,1) = arcmax

c-----------------------------------------------------------------------
c Open namelist input file 'tcabr.in' and read input data
c-----------------------------------------------------------------------
c changed by Wingen: Iadj and Icur are read elsewhere and included here 
c					 through common block
 
c      CALL TCABRCPREAD (useTCABRCPcoil, ECadj, ECcur)
 
 
c      if(.not. useTCABRCPcoil)   RETURN	! don't need TCABR CP-coils
 
c-----------------------------------------------------------------------
c More initial stuff
c-----------------------------------------------------------------------
 
      nbands = 0
      DO L=1, mxbands			! find number of coil bands
       if(ECloops(L) .ne. 0) then
        nbands = nbands + 1
       end if
      END DO
 
c      write(81,*) ''		! write more details to file 81 = case2
c      write(81,*) nbands,' TCABR CP-coil bands'
 
      if(nbands .eq. 0)    RETURN   	! no CP-coils are defined

      DO L=1, nbands
       ECadj(2,L)  = abs(ECadj(2,L))	! make arcmax     "
      END DO 
 
      nloops = ECloops			! f90 array assign
      curnt  = ECcur
 
  ! write input information to file 81 = case2
c      write(81,'(a,16i4)') ' nloops: ',(nloops(L),L=1,nbands)
c      write(81,'(a,16i4)') ' npoints:',(npoints(L),L=1,nbands)
C      WRITE(81,*) ''
c      write(81,*) 'ECgeom:'
c      do L=1,nbands
c       write(81,'(8f10.3 )') (ECgeom(k,L),k=1,7)
c      end do
C      WRITE(81,*) ''
c      WRITE(0,*) 'ECcur:'
c      do L=1,nbands
c       WRITE(0,'(8f10.0 )') (ECcur(k,L),k=1,nloops(L))
c      end do
c      WRITE(0,*) ''
c      write(0,*) 'ECadj:'
c      do L=1,nbands
c       write(0,'(8f10.3 )') (ECadj(k,L),k=1,3)
c      end do
 
c adjust currents by scale factors from ECajd(3)
      DO  L=1,nbands  	! do all bands of CP-coils
       DO  k=1,nloops(L)
        scalec = ECadj(3,L)
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
 
c      WRITE(81,*) 'CURNT:'
c      do L=1,nbands
c       WRITE(81,'(10f8.0)') (curnt(k,L),k=1,nloops(L))
c      end do
C      write(0,*) 'kuse:'
C      do L=1,nbands
C       write(0,'(24i3)') (kuse(k,L),k=1,nloops(L))
C      end do
 
 
c=======================================================================
c            Fill Position Vector Array	
c=======================================================================
C To write points to fort.15, usually named OUT.DEBUG:
C        WRITE(15,*) ' '
C        WRITE(15,*) '**************   TCABR CP-COIL   *****************'
C        WRITE(15,*) 'nbands:', nbands
C        WRITE(15,*) 'nloops:', nloops
C        WRITE(15,*) ' '

c-----------------------------------------------------------------------
c Loop through all requested toroidal bands of CP-coils
c-----------------------------------------------------------------------
c Toroidal angle arithmetic is done in degrees; 
c  convert to radians to evaluate trigonometric functions
c Here we consisently use indices i,j,k,L to identify (respectively)
c  Cartesian component, segment, loop or coil, band or row.  
 
      DO  L=1, nbands  		! do all bands
 
       R1        = ECgeom(1,L)                
       Z1        = ECgeom(2,L) 
       R2        = ECgeom(3,L) 
       Z2        = ECgeom(4,L) 
       phianchor = ECgeom(5,L)
       phispan   = abs(ECgeom(6,L))	! must be .ge. 0
       dphi      = ECgeom(7,L)
       addangl   = ECadj(1,L)
       arcmax    = abs(ECadj(2,L))	! must be .ge. 0

       phianchor = phianchor + addangl
       
       IF(ECadj(2,L) .eq. 0. .and. ECadj(3,L) .ne. 0.) THEN
        write(0,*)'TCABRCPGEOM: ECadj(2,*) = arcmax = 0 not allowed.'
	write(0,*)'TCABRCPGEOM: Stopping.'
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
       DO
	 ntorsegs = ntorsegs + 1
	 segarc = phispan/ntorsegs
	 IF (segarc .le. arcmax)  EXIT		! satisfied segarc
 
	 IF(ntorsegs .gt. 100) THEN		! stop runaway loop
	  write(0,*)'TCABRCPGEOM: Coil segment arc size setup ran away.'
	  write(0,*)'TCABRCPGEOM: Stopping.'
	  STOP
	 END IF
       END DO				! end of 'WHILE'
 
c       write(81,*) 'segarc, ntorsegs for this band'
c       write(81,'(f10.3,i5)') segarc, ntorsegs
 
 
    !-------------------------------------------------------------------
    ! Loop through all coils (loops) in band L of CP-coils
    !-------------------------------------------------------------------
 
       DO 100 k=1,nloops(L) 	! do all loops in CP-coil band L
 
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
 
      END DO	! end DO  L=1,nbands  	! do all bands of CP-coils
 
 
c=======================================================================
c            Fill Direction Vector Array	
c=======================================================================
 
      DO  l=1,nbands  	! do all bands of CP-coils
 
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
 
      END DO	! end DO  L=1,nbands  	! do all bands of CP-coils
 

!       Write data for plotgeom
!-----------------------------------------------------------------------

c      write(93,'(I3,A)'),nbands,'   bands. CP coils'
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
c TCABRCPGEOM become unavailable for execution by following procedures 
c (like ENTRY ...). This restriction does not apply to all other 
c variables, which are still available after the return.
c Therefore, before returning, here we define duplicates of those dummy
c variables that will be needed in the rest of this subroutine.
c
c We also multiply here the input current array by number of turns.
 
      nbandsd = nbands
      nloopsd = nloops
      curntd  = curnt
      
      curntw  = curnt*ncturns		! multiply by number of turns
 
      
C      write(0,*) 'LEAVING TCABRCPGEOM'		! for DEBUG
 
      RETURN 
 
      END 	SUBROUTINE	TCABRCPGEOM
c=======================================================================
