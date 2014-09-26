

c Contents of this file:
c
c	subroutine D3ICOILS
c	  entry BD3ICOILS
c	subroutine D3IGEOM
c
c
c=======================================================================

      subroutine D3ICOILS

c=======================================================================
c                                            By: M. Schaffer 2006 jul 31
c                                 Last Modified: M. Schaffer 2007 oct 25
c-----------------------------------------------------------------------
c D3DICOILS defines the geometry of and calculates the magnetic field
c from one or more specified models of the DIII-D internal coils 
c (I-coils).
c This subroutine handles all the tasks needed to calculate this field.
c
c This subroutine has several main parts:
c  1) Set flags for components to use.
c  2) Call routines to set up models to be used.
c  3) Write important parameters to appropriate output files.
c  4) Upon separate entry at BD3ICOILS, call subroutines to calculate
c     the I-coil magnetic field at a point (x,y,z).
c
c  Array dimensions set by INCLUDE ICOILSD3.i:
c   nIloops, # of "loops" in I-coil array
c   nIsegs,  # of "segments" in each loop (same for all)
c
c This version orders array indices as (i,j,k). Previously (k,j,i).
c	i for x,y,z; j for segment; k for loop
c    
c Units: SI (mks)
c Fully double precision.
c-----------------------------------------------------------------------

      IMPLICIT NONE
      INCLUDE 'icoilsd3.i'	! included parameter & variable types
				! are declared in the include file.

c-----------------------------------------------------------------------

      common /consts/  pi, twopi, cir, rtd, dtr
      common /flags2/ iPROBE, iSURFMN, iTRIP3D,
     &                iD3Dcoils, iNSTXcoils, iITERcoils, iFDFcoils
      common /d3icoil/ curntIc, addanglIU, addanglIL, scaleIU, scaleIL
      common /d3iloops/ xis, divs


!  Common variables:
      real*8  pi, twopi, cir, rtd, dtr
      integer iPROBE, iSURFMN, iTRIP3D
      integer iD3Dcoils, iNSTXcoils, iITERcoils, iFDFcoils
      real*8  curntIc(nIloops), addanglIU, addanglIL, scaleIU, scaleIL
      real*8  xis(3,nIsegs,nIloops), divs(4,nIsegs,nIloops)

!  Arguments at ENTRY BD3BUS:
      real*8  x, y, z, bx, by, bz

!  Local and Namelist Variables:
      integer nicsegs(nIloops)		! permits use of POLYGONB
      logical iIcoil
      integer kuseI(nIloops)
      real*8  curnt(nIloops)
      integer i,j,k	! i for x,y,z; j for segment; k for loop

      SAVE

c-----------------------------------------------------------------------
c Set some default values
c-----------------------------------------------------------------------
c Variables initialized by data statement acquire the SAVE attribute.

      data kuseI /nIloops*0/

      DO k=1, nIloops
       nicsegs(k) = nIsegs	! nsegs from icoilsd3.i, same for all
      END DO

c-----------------------------------------------------------------------
c Set flags
c-----------------------------------------------------------------------

      iIcoil = .false.
      DO k=1, nIloops
       if(curntIc(k) .ne. 0) then
        iIcoil = .true.		! at least 1 i-coil loop has current
        iD3Dcoils = 1		! at least 1 D3D coil model is used
	kuseI(k) = 1		! flag for this active I-coil loop
       end if
      END DO

c-----------------------------------------------------------------------
c Set up I-coil geometry
c-----------------------------------------------------------------------


      if (iIcoil)  CALL D3IGEOM (kuseI)
      

c-----------------------------------------------------------------------
c Write some information about this model to file CASE
c-----------------------------------------------------------------------

      if(iIcoil) then			! Write I-coil information
       write(80,*) ''
       write(80,*) 'I-coil, addanglIU, addanglIL (deg) currents (A):'
       write(80,1010) addanglIU,addanglIL
       write(80,1010) (curntIc(i),i=1,6)
       write(80,1010) (curntIc(i),i=7,12)
      end if


      RETURN		! End of D3D I-coil Model setup




c=======================================================================
c Calculate magnetic field for all requested I-coils.
c
c ENTRY BD3ICOILS 
c  Receives the point (x,y,z) at which the field is to be calculated;
c  Calls POLYGONB (or D3IB) to evaluate the field at (x,y,z);
c  Returns the Cartesian magnetic field components (bx,by,bz).
c-----------------------------------------------------------------------



      ENTRY BD3ICOILS (x, y, z, bx, by, bz)

      bx = 0.0d0		! initialize magnetic field components
      by = 0.0d0		! that will be calculated and returned
      bz = 0.0d0


      IF (.not. iIcoil) THEN	! return 0 if no I-coil is active

       RETURN

      ELSE

       CALL POLYGONB(nIloops, nIsegs, nIloops, nicsegs, kuseI,
     &               xis, divs, curntIc, x, y, z, bx, by, bz)

      END IF


      RETURN

c-----------------------------------------------------------------------
c Format statements
c-----------------------------------------------------------------------
 1003 format(1x,30i3)
 1010 format(1x,8f9.0)
 1014 format(1x,8f9.4)


      END 	SUBROUTINE 		D3ICOILS
c=======================================================================








c=======================================================================

      subroutine D3IGEOM (kuse)

c=======================================================================
c			        Based on IGEOM: MJ Schaffer, 2005 mar 05
c                                   Started by: MJ Schaffer, 2006 jul 23
c                                      Last Modified:   MJS, 2007 oct 25
c-----------------------------------------------------------------------
c Based on previous subroutine CGEOM.
c D3IGEOM writes arrays of DIII-D C-coil geometry parameters 
c in the form required by MJS's routines D3IB, BIOTLOOP.
c
c Input variable kuse(nIloops) is a flag signaling active loops.
c
c The I-coil is represented as a set of 12 closed polygon current loops
c each comprised of straight-line segments.
c
c  niloops= number of loops (niloops = 12 for DIII-D I-coil)
c  nisegs= number of segments in an I-coil loop
c  xis(i,j,k), i=1,2,3 j=1...nisegs are (x,y,z) Cartesian position 
c    vectors of nisegs coil-defining points of loop(k), k=1,...,niloops.
c  divs(i,j,k), i=1,2,3,4 j=1...nisegs are Cartesian direction vectors  
c    (3 elements) and lengths (4th element) of nisegs straignt segments.
c    The start point of segment j is row j of xis, its end point is 
c    row j+1 of xis, except end point of segment nisegs is row 1 of xis.
c  curntIc(k), current in kth I-coil, k=1,...,niloops.
c
c The I-coil arrays can be rotated toroidally by addanglIU and addanglIL
c (deg), and their current distributions can be multiplied by constant
c factors scaleIU and scaleIL (dimensionless) respectively. Addangl.. 
c and scale.. are provided for convenience when doing scans of the 
c I-coil B-field.
c
c Coordinate system: All geometry information is written in 
c    right-handed Cartesian coordinates, as also used by BIOTLOOP.
c
c Note: For convenience in relating to DIII-D machine "geographical" 
c   I-coil nomenclature and prior TRIP3D practice, the 12 coils are 
c   indexed in the order 1,...12 
c          <--->  IU30, IU90, IU150, IU210, IU270, IU330,
c                 IL30, IL90, IL150, IL210, IL270, IL330.
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
      INCLUDE 'icoilsd3.i'	! included parameter & variable types
				! are declared in the include file.

c ======================================================================
c SIMPLE "PICTURE FRAME" I-COIL SET
c
c Cannot assign values to subscripted variables in parameter statements.
c Use mipolsegs=1; extra pts on vertical runs are unnecessary, wasteful.
c NOTE:  D3IGEOM can presently do ONLY mipolsegs=1  !!!
c-----------------------------------------------------------------------
c
c      integer niloops, nuploops, nloloops
c      integer nitorsegs, mipolsegs
c      integer niturns           !number of electrical turns in an I-coil
c      integer nisegs
c
c      parameter (niloops= 12, nuploops= 6, nloloops= 6)
c      parameter (nitorsegs= 6, mipolsegs= 1)
c      parameter (niturns= 1)
c      parameter (nisegs= 2*nitorsegs + 2*mipolsegs)
c
c
c SPECIFIC PARAMETERS FOR DIII-D: * * * * * * * * * * * * * * * * * *  *
c
c      real*8 torarci  			! toroidal span of one loop, deg
c      data   torarci /51.72/
c
c      real*8 phictri(niloops)   !tor angles of niloops loop centers, deg
c      data   phictri /-032.7,-087.3,-152.7,-207.3,-272.7,-327.3,
c     &                -032.7,-087.3,-152.7,-207.3,-272.7,-327.3 /
c
c Original 2003 Mar geometric parameters:
c      real*8 Riu1, Ril1, Ziu1, Zil1   !radii, heights of upper loop arcs
c      data   Riu1, Ril1, Ziu1, Zil1  /2.184, 2.394, +1.012, +0.504/
c      real*8 Riu2, Ril2, Ziu2, Zil2   !radii, heights of lower loop arcs
c      data   Riu2, Ril2, Ziu2, Zil2  /2.394, 2.184, -0.504, -1.012/
c
c Revised 2006 Dec parameters from 4:1 "large" drawing:
C      real*8 Riu1, Ril1, Ziu1, Zil1   !radii, heights of upper loop arcs
C      data   Riu1, Ril1, Ziu1, Zil1  /2.164, 2.373, +1.016, +0.504/
C      real*8 Riu2, Ril2, Ziu2, Zil2   !radii, heights of lower loop arcs
C      data   Riu2, Ril2, Ziu2, Zil2  /2.373, 2.164, -0.504, -1.016/
c
C Specific parameters for DIII-D:
C      parameter (torarci = 51.72)   !toroidal extent of a loop, deg
C      parameter (Riu1= 2.184, Ril1= 2.394, Riu2= 2.394, Ril2= 2.184) !m
C      parameter (Ziu1= 1.012, Zil1= 0.504, Ziu2=-0.504, Zil2=-1.012) !m
c ======================================================================


      common /consts/   pi, twopi, cir, rtd, dtr
      common /d3icoil/  curntIc, addanglIU, addanglIL, scaleIU, scaleIL
      common /d3iloops/ xis, divs


!  Common variables:
      real*8  pi, twopi, cir, rtd, dtr
      real*8  curntIc(nIloops), addanglIU, addanglIL, scaleIU, scaleIL
      real*8  xis(3,nIsegs,nIloops), divs(4,nIsegs,nIloops)

!  Arguments:
      integer kuse(nIloops)
      intent(in) :: kuse

!  Local Variables:
      real*8  segarc, segz
      real*8  phicent, phistart, phiend, phipoint
      real*8  Ru, Zu, Rl, Zl
      integer i, j, k, j1, j2, j3, j4, jend

      SAVE

c-----------------------------------------------------------------------
c  Initialize
c-----------------------------------------------------------------------
      segarc = torarci/float(nitorsegs)     !arc of one toroidal segment
      segarc = abs(segarc)
      segz   = (Ziu1 - Zil1)/float(mipolsegs) !length of a vert. segment
      segz   = abs(segz)
      
      DO 10 k=1,niloops        ! initialize all I-coil loops
       DO 11 j=1,nisegs        ! and all straight line segments
        xis(1,j,k) = 0.        ! (x,y,z) point arrays
        xis(2,j,k) = 0.
        xis(3,j,k) = 0.

        divs(1,j,k) = 0.       ! segment length vectors:
        divs(2,j,k) = 0.
        divs(3,j,k) = 0.
        divs(4,j,k) = 0.
   11  CONTINUE
   10 CONTINUE


c=======================================================================
c            Fill Position Vector Array	
c=======================================================================
c Subroutine BIOT assumes that segment currents flow in the same direc-
c  tion as increasing index 'j' of its defining points in xics(i,j,k). 
c  Therefore, here we define the I-coil geometry in the order of a 
c  positive loop current according to DIII-D convention and PTDATA. 
c  A positive I-coil current makes a B-field directed outward from the
c  plasma surface.
c-----------------------------------------------------------------------
C To write points to fort.15, usually named OUT.DEBUG:
C        WRITE(15,*) ' '
C        WRITE(15,*) '*******************   I-COIL   *******************'
C        WRITE(15,*) 'niloops, nisegs:', niloops, nisegs
C        WRITE(15,*) 'segarc(deg), segz(m):', segarc, segz
C        WRITE(15,*) ' '

c Toroidal angle arithmetic is done in degrees; 
c  convert to radians to evaluate trigonometric functions
c Here we consisently use indices i,j,k to identify (respectively)
c  Cartesian component, segment, loop (coil).


      DO 100 k=1,niloops        ! do all I-coil loops
        
        if (k .le. nuploops) then		! upper I-coil set
	 phicent  = phictri(k) + addangliu	! center of loop k
         Ru = Riu1                  
         Zu = Ziu1
         RL = RiL1
         ZL = ZiL1
        else					! lower I-coil set
	 phicent  = phictri(k) + addanglil	! center of loop k
         Ru = Riu2                  
         Zu = Ziu2
         RL = RiL2 
         ZL = ZiL2
        end if
        phistart = phicent - torarci/2.
        phiend   = phicent + torarci/2.

C To write points to fort.15:
C        WRITE(15,*) ' '
C        WRITE(15,*) 
C     &   '      phicent                  start                    end'
C        WRITE(15,*) phicent, phistart, phiend
C        WRITE(15,*) 
C     &   '      Ru                       RL                       ZL'
C        WRITE(15,*) Ru,RL,ZL
C        WRITE(15,*) ' '

c Go around loop k

        j = 0                 ! j is segment index counter for this loop
       
        ! lower conductor, increasing coordinate phi 
        
        DO 20 j1=1,nitorsegs+1   ! nitorsegs arcs, nitorsegs+1 points
         j = j+1                 ! j is cumulative segment index counter
         phipoint = phistart + (j1-1)*segarc
         xis(1,j,k) = RL*cos(phipoint*dtr)
         xis(2,j,k) = RL*sin(phipoint*dtr)
         xis(3,j,k) = ZL

C To write points to fort.15:
C         WRITE(15,*) 
C     &   ' phipoint(deg)           Rpoint(m)      k        j'
C         WRITE(15,*) 
C     &   ' Xpoint(m)               Ypoint(m)               Zpoint(m)'
C         WRITE(15,*) phipoint, RL, k, j
C         WRITE(15,*) xis(1,j,k),xis(2,j,k),xis(3,j,k)
C         WRITE(15,*) ' '

   20   CONTINUE

        ! first vertical conductor, going upward (increasing Z)
        
        DO 30 j2=1,mipolsegs-1   ! mipolsegs sections, mipolsegs-1 pts
         j = j+1                 ! j is cumulative segment index counter
         phipoint = phiend
         xis(1,j,k) = Ru*cos(phipoint*dtr)
         xis(2,j,k) = Ru*sin(phipoint*dtr)
         xis(3,j,k) = ZL + j2*segz

C To write points to fort.11:
C         WRITE(*,*) 'going down'

   30   CONTINUE

         ! upper conductor, decreasing coordinate phi
         
         DO 40 j3=1,nitorsegs+1  ! nitorsegs arcs, nitorsegs+1 points
         j = j+1                 ! j is cumulative segment index counter
         phipoint = phiend - (j3-1)*segarc
         xis(1,j,k) = Ru*cos(phipoint*dtr)
         xis(2,j,k) = Ru*sin(phipoint*dtr)
         xis(3,j,k) = Zu

C To write points to fort.15:
C         WRITE(15,*) 
C     &   ' phipoint(deg)           Rpoint(m)'
C         WRITE(15,*) 
C     &   ' Xpoint(m)               Ypoint(m)               Zpoint(m)'
C         WRITE(15,*) phipoint, Ru, k, j
C         WRITE(15,*) xis(1,j,k),xis(2,j,k),xis(3,j,k)
C         WRITE(15,*) ' '

   40   CONTINUE

        ! last vertical conductor, going downward (decreasing Z)
        
        DO 50 j4=1,mipolsegs-1   ! mipolsegs sections, mipolsegs-1 pts
         j = j+1                 ! j is cumulative segment index counter
         phipoint = phistart
         xis(1,j,k) = RL*cos(phipoint*dtr)
         xis(2,j,k) = RL*sin(phipoint*dtr)
         xis(3,j,k) = Zu - j4*segz

C         write(*,*) 'going up'
   50   CONTINUE

  100 CONTINUE


c-----------------------------------------------------------------------
c            Fill Direction Vector Array	
c-----------------------------------------------------------------------

      DO 200 k=1,niloops
       DO 202 j=1,nisegs
          
        if (j.ne.nisegs) then  
         jend = j+1              ! end point for all segments but last
        else
         jend = 1                ! end point for last segment
        endif

        ! segment length vector:
        divs(1,j,k) = xis(1,jend,k) - xis(1,j,k)
        divs(2,j,k) = xis(2,jend,k) - xis(2,j,k)
        divs(3,j,k) = xis(3,jend,k) - xis(3,j,k)

        ! segment length:
        divs(4,j,k) = 
     &       sqrt(divs(1,j,k)**2 + divs(2,j,k)**2 + divs(3,j,k)**2)

        ! segment direction unit vector:
        DO 203 i=1,3
         divs(i,j,k) = divs(i,j,k)/divs(4,j,k)
  203   CONTINUE

  202  CONTINUE
  200 CONTINUE


c-----------------------------------------------------------------------
c       Adjust Wire Current by Scale Factors
c-----------------------------------------------------------------------

      DO 300 k=1,nuploops
       curntIc(k) = scaleIU*curntIc(k)
  300 CONTINUE
  
      DO 301 k=nuploops + 1, nuploops + nloloops
       curntIc(k) = scaleIL*curntIc(k)
  301 CONTINUE


      RETURN
      

      END 	SUBROUTINE	D3IGEOM
c=======================================================================
