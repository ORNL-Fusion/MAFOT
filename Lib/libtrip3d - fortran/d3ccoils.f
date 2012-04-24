

c Contents of this file:
c
c	subroutine D3CCOILS
c	  entry BD3CCOILS
c	subroutine D3CGEOM
c
c
c=======================================================================

      subroutine D3CCOILS

c=======================================================================
c                                            By: M. Schaffer 2006 jul 25
c                                 Last Modified: M. Schaffer 2007 oct 25
c-----------------------------------------------------------------------
c D3DCCOILS defines the geometry of and calculates the magnetic field
c from one or more specified models of the DIII-D correction coils 
c (C-coils).
c This subroutine handles all the tasks needed to calculate this field.
c
c This subroutine has several main parts:
c  1) Set flags for components to use.
c  2) Call routines to set up models to be used.
c  3) Write important parameters to appropriate output files.
c  4) Upon separate entry at BD3CCOILS, call subroutines to calculate
c     the C-coil magnetic field at a point (x,y,z).
c
c  Array dimensions set by INCLUDE CCOILSD3.i:
c   nCloops, # of "loops" in C-coil array
c   nCsegs,  # of "segments" in each loop (same for all)
c
c This version orders array indices as (i,j,k). Previously (k,j,i).
c	i for x,y,z; j for segment; k for loop
c
c Units: SI (mks)
c Fully double precision.
c-----------------------------------------------------------------------

      IMPLICIT NONE
      INCLUDE 'ccoilsd3.i'	! included parameter & variable types
				! are declared in the include file.

c-----------------------------------------------------------------------

      common /consts/  pi, twopi, cir, rtd, dtr
      common /flags2/ iPROBE, iSURFMN, iTRIP3D,
     &                iD3Dcoils, iNSTXcoils, iITERcoils, iFDFcoils
      common /d3ccoil/ curntC, curntw, addanglC, scaleC
      common /d3cloops/ xcs, dcvs


!  Common variables:
      real*8  pi, twopi, cir, rtd, dtr
      integer iPROBE, iSURFMN, iTRIP3D
      integer iD3Dcoils, iNSTXcoils, iITERcoils, iFDFcoils
      real*8  curntC(nCloops), curntw(nCloops), addanglC, scaleC
      real*8  xcs(3,nCsegs,nCloops), dcvs(4,nCsegs,nCloops)

!  Arguments at ENTRY BD3BUS:
      real*8  x, y, z, bx, by, bz

!  Local and Namelist Variables:
      integer nccsegs(nCloops)		! permits use of POLYGONB
      logical iCcoil
      integer kuseC(ncloops)
      integer i,j,k

      SAVE

c-----------------------------------------------------------------------
c Set some default values
c-----------------------------------------------------------------------
c Variables initialized by data statement acquire the SAVE attribute.

      data kusec /nCloops*0/

      DO k=1, nCloops
       nccsegs(k) = nCsegs	! nsegs from ccoilsd3.i, same for all
      END DO

c-----------------------------------------------------------------------
c Set flags
c-----------------------------------------------------------------------

      iCcoil = .false.
      DO k=1, nCloops
       if(curntC(k) .ne. 0) then
        iCcoil = .true.		! at least 1 C-coil loop has current
        iD3Dcoils = 1		! at least 1 D3D coil model is used
	kuseC(k) = 1		! flag for this active C-coil loop
       end if
      END DO

c-----------------------------------------------------------------------
c Set up C-coil geometry
c-----------------------------------------------------------------------


      if (iCcoil)  CALL D3CGEOM (kusec)


c-----------------------------------------------------------------------
c Write some information about this model to file CASE
c-----------------------------------------------------------------------

      if(iCcoil) then			! Write C-coil information
       write(80,*) ''
       write(80,*) 'C-coil, addanglC (deg) and currents (A):'
       write(80,1010) addanglC
       write(80,1010) (curntC(i),i=1,6)
      end if


      RETURN		! End of D3D C-coil Model setup





c=======================================================================
c Calculate magnetic field for all requested C-coils.
c
c ENTRY BD3CCOILS
c  Receives the point (x,y,z) at which the field is to be calculated;
c  Calls POLYGONB (or D3CB) to evaluate the field at (x,y,z);
c  Returns the Cartesian magnetic field components (bx,by,bz).
c-----------------------------------------------------------------------



      ENTRY BD3CCOILS (x, y, z, bx, by, bz)

      bx = 0.0d0		! initialize magnetic field components
      by = 0.0d0		! that will be calculated and returned
      bz = 0.0d0


      IF (.not. iCcoil) THEN	! return 0 if no C-coil is active

       RETURN

      ELSE

       CALL POLYGONB(nCloops, nCsegs, nCloops, nccsegs, kuseC,
     &                xcs, dcvs, curntw, x, y, z, bx, by, bz)

      END IF


      RETURN

c-----------------------------------------------------------------------
c Format statements
c-----------------------------------------------------------------------
 1003 format(1x,30i3)
 1010 format(1x,8f9.0)
 1014 format(1x,8f9.4)


      END 	SUBROUTINE 		D3CCOILS
c=======================================================================








c=======================================================================

      subroutine D3CGEOM (kuse)

c=======================================================================
c			       Based on CGEOM: MJ Schaffer, 2004 aug 09
c                                  Started by: MJ Schaffer, 2006 jul 23
c                                     Last Modified:   MJS, 2007 oct 25
c----------------------------------------------------------------------
c Based on previous subroutine CGEOM.
c D3CGEOM writes arrays of DIII-D C-coil geometry parameters 
c in the form required by MJS's routines D3CB,, BIOTLOOP.
c
c Input variable kuse(nCloops) is a flag signaling active loops.
c
c The C-coil is represented as a set of 6 closed polygon current loops
c each comprised of straight-line segments.
c
c  ncloops= number of loops (ncloops = 6 for DIII-D C-coil)
c  ncsegs= number of segments in a loop
c  xcs(k,j,i), i=1,2,3 j=1...ncsegs are (x,y,z) Cartesian position 
c    vectors of ncsegs coil-defining points of loop(k), k=1,...,niloops.
c  dcvs(j,i), i=1,2,3,4 j=1...ncsegs are Cartesian direction vectors  
c     (3 elements) and lengths (4th element) of ncsegs straignt segments.
c     The start point of segment j is row j of xcs, its end point is 
c     row j+1 of xcs, except end point of segment ncsegs is row 1 of xcs.
c  curntC(k), current (A/turn) in kth C-coil, k=1,...,ncloops.
c  ncturns= number of turns in each C-coil = 4
c
c The C-coil array can be rotated toroidally by addanglc (deg) and its 
c current distribution can be multiplied by constant factor scaleC 
c (dimensionless). Addanglec and scaleC are provided for convenience 
c when doing scans of the C-coil B-field.
c
c Coordinate system: Coil geometry information is all written in 
c    right-handed Cartesian coordinates, as used by D3CB, BIOTLOOP.
c    
c Note: For convenience in relating to DIII-D machine "geographical" 
c   C-coil nomenclature and prior TRIP3D practice, the 6 coils are 
c   indexed in the order 1,...6  <---> C79, C139, C199, C259, C319, C19.
c   However, the coil toroidal angles here are all written in
c   right-handed cylindrical coordinates, so coordinate toroidal angle
c   negative of machine geographical toroidal angle.
c
c This version orders array indices as (i,j,k). Previously (k,j,i).
c    
c Units: SI (mks)
c Fully double precision.
c-----------------------------------------------------------------------

      IMPLICIT NONE
      INCLUDE 'ccoilsd3.i'	! included parameter & variable types
				! are declared in the include file.

c ======================================================================
c SIMPLE "PICTURE FRAME" C-COIL SET
c
c Cannot assign values to subscripted variables in parameter statements.
c Use mcpolsegs=1; extra pts on vertical runs are unnecessary, wasteful.
c-----------------------------------------------------------------------
c
c      integer ncloops, nctorsegs, mcpolsegs, ncsegs
c      integer ncturns            !number of electrical turns in a C-coil
c
c      parameter (ncloops= 6, nctorsegs= 4, mcpolsegs= 1)
c      parameter (ncturns= 4)
c      parameter (ncsegs= 2*nctorsegs + 2*mcpolsegs)
c
c      real*8 torarcc 			! toroidal span of one loop, deg
c      data   torarcc /58.0/
c
c      real*8 phictrc(ncloops)   !tor angles of ncloops loop centers, deg
c      data   phictrc  / -079., -139., -199., -259., -319., -019. /
c
c      real*8 Rc, Zcu, Zcl 		!radius, heights of loop arcs, m
c      data   Rc, Zcu, Zcl /3.23, 0.8, -0.8/
c
c ======================================================================


      common /consts/   pi, twopi, cir, rtd, dtr
      common /d3ccoil/  curntC, curntw, addanglC, scaleC
      common /d3cloops/ xcs, dcvs


!  Common variables:
      real*8  pi, twopi, cir, rtd, dtr
      real*8  curntC(nCloops), curntw(nCloops), addanglC, scaleC
      real*8  xcs(3,nCsegs,nCloops), dcvs(4,nCsegs,nCloops)

!  Arguments:
      integer kuse(nCloops)
      intent(in) :: kuse

!  Local Variables:
      real*8  segarc, segz
      real*8  phicent, phistart, phiend, phipoint
      integer i, j, k, j1, j2, j3, j4, jend

      SAVE

c-----------------------------------------------------------------------
c  Initialize
c-----------------------------------------------------------------------
      segarc= torarcc/float(nctorsegs)    !arc of one toroidal segment
      segarc= abs(segarc)
      segz  = (Zcu - Zcl)/float(mcpolsegs)  !length of one vert. segment
      segz   = abs(segz)

      DO 10 k=1,ncloops  	! initialize all C-coil loops
       DO 11 j=1,ncsegs 	! and all straight line segments
        xcs(1,j,k) = 0.  	! (x,y,z) point arrays
        xcs(2,j,k) = 0.
        xcs(3,j,k) = 0.

        dcvs(1,j,k) = 0. 	! initialize all C-coil 
        dcvs(2,j,k) = 0. 	! segment length vectors
        dcvs(3,j,k) = 0.
        dcvs(4,j,k) = 0.
   11  CONTINUE
   10 CONTINUE

c=======================================================================
c            Fill Position Vector Array	
c=======================================================================
C To write points to fort.15, usually named OUT.DEBUG:
C        WRITE(15,*) ' '
C        WRITE(15,*) '******************   C-COIL   ******************'
C        WRITE(15,*) 'ncloops, ncsegs:', ncloops, ncsegs
C        WRITE(15,*) 'segarc(deg), segz(m):', segarc, segz
C        WRITE(15,*) ' '

c Toroidal angle arithmetic is done in degrees; 
c  convert to radians to evaluate trigonometric functions
c Here we consisently use indices i,j,k to identify (respectively)
c  Cartesian component, segment, loop or coil.


      DO 100 k=1,ncloops			! do all C-coil loops
        
        phicent  = phictrc(k) + addanglc	! center of loop k
        phistart = phicent - torarcc/2.
        phiend   = phicent + torarcc/2.
        
CC To write points to fort.11:
C        WRITE(11,*) ' '
C        WRITE(11,*) 
C     &   '      phicent                  start                    end'
C        WRITE(11,*) phicent, phistart, phiend

c Go around loop k

        j = 0                 ! j is segment index counter for this loop
       
        ! upper conductor, increasing coordinate phi 
        ! is DIII-D direction of positive C-coil current
        
        DO 20 j1=1,nctorsegs+1     ! nctorsegs arcs, nctorsegs+1 points
         j = j+1                   ! j is segment index counter for loop
         phipoint = phistart + (j1-1)*segarc
         xcs(1,j,k) = Rc*cos(phipoint*dtr)
         xcs(2,j,k) = Rc*sin(phipoint*dtr)
         xcs(3,j,k) = Zcu
         
CC To write points to fort.15:
C         WRITE(15,*) 
C     &   ' phipoint(deg)           Rpoint(m)'
C         WRITE(15,*) 
C     &   ' Xpoint(m)               Ypoint(m)               Zpoint(m)'
C         WRITE(15,*) phipoint, Rc 
C         WRITE(15,*) xcs(1,j,k),xcs(2,j,k),xcs(3,j,k)

   20   CONTINUE

        ! first vertical conductor, going downward (decreasing Z)
        ! is DIII-D direction of positive C-coil current
        
        DO 30 j2=1,mcpolsegs-1  ! mcpols straight sections, mcpols-1 pts
         j = j+1                   ! j is segment index counter for loop
         phipoint = phiend
         xcs(1,j,k) = Rc*cos(phipoint*dtr)
         xcs(2,j,k) = Rc*sin(phipoint*dtr)
         xcs(3,j,k) = Zcu - j2*segz

C         write(11,*) 'going down'
   30   CONTINUE

         ! do lower conductor, decreasing coordinate phi
         ! is DIII-D direction of positive C-coil current
         
         DO 40 j3=1,nctorsegs+1     ! nctorsegs arcs, nctorsegs+1 points
         j = j+1                   ! j is segment index counter for loop
         phipoint = phiend - (j3-1)*segarc
         xcs(1,j,k) = Rc*cos(phipoint*dtr)
         xcs(2,j,k) = Rc*sin(phipoint*dtr)
         xcs(3,j,k) = Zcl
         
CC To write points to fort.15:
C         WRITE(15,*) 
C     &   ' phipoint(deg)           Rpoint(m)'
C         WRITE(15,*) 
C     &   ' Xpoint(m)               Ypoint(m)               Zpoint(m)'
C         WRITE(15,*) phipoint , Rc
C         WRITE(15,*) xcs(1,j,k),xcs(2,j,k),xcs(3,j,k)
C         WRITE(15,*) ' '

   40   continue

        ! do last vertical conductor, going upward (increasing Z)
        ! is DIII-D direction of positive C-coil current
        DO 50 j4=1,mcpolsegs-1  ! mcpols straight sections, mcpols-1 pts
         j = j+1                   ! j is segment index counter for loop
         phipoint = phistart
         xcs(1,j,k) = Rc*cos(phipoint*dtr)
         xcs(2,j,k) = Rc*sin(phipoint*dtr)
         xcs(3,j,k) = Zcl + j4*segz
         
C         write(*,*) 'going up'
   50   CONTINUE

  100 CONTINUE


c-----------------------------------------------------------------------
c            Fill Direction Vector Array	
c-----------------------------------------------------------------------

      DO 200 k=1,ncloops
       DO 202 j=1,ncsegs
          
        if (j.ne.ncsegs) then  
         jend = j+1              ! end point for all segments but last
        else
         jend = 1                ! end point for last segment
        endif

        ! segment length vector:
        dcvs(1,j,k) = xcs(1,jend,k) - xcs(1,j,k)
        dcvs(2,j,k) = xcs(2,jend,k) - xcs(2,j,k)
        dcvs(3,j,k) = xcs(3,jend,k) - xcs(3,j,k)

        ! segment length:
        dcvs(4,j,k) = 
     &       sqrt(dcvs(1,j,k)**2 + dcvs(2,j,k)**2 + dcvs(3,j,k)**2)

        ! segment direction vector:
        DO 203 i=1,3
         dcvs(i,j,k) = dcvs(i,j,k)/dcvs(4,j,k)
  203   CONTINUE

  202  CONTINUE
  200 CONTINUE


c-----------------------------------------------------------------------
c       Adjust Wire Current to Scale Factor and Number of Turns
c-----------------------------------------------------------------------

      DO 300 k=1,ncloops
       curntc(k) = scaleC*curntc(k)
       curntw(k) = ncturns*curntc(k)
  300 CONTINUE

CCCCC
C      WRITE(0,*) 'ncloops',ncloops
C      WRITE(0,*) 'ncsegs', ncsegs
C      do j=1,ncsegs
C      WRITE(0,*) 'xcs',j, (xcs(1,j,i), i=1,3)
C      end do
C      do j=1,ncsegs
C      WRITE(0,*) 'dcvs',j,(dcvs(1,j,i), i=1,4)
C      end do
C      WRITE(0,*) 'curntw', curntw(1)
C
C      STOP
CCCCC


      RETURN
      

      END 	SUBROUTINE	D3CGEOM
c=======================================================================
