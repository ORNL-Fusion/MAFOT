

c=======================================================================

      subroutine POLYGONB (loopsdim, segsdim, nloops, nsegs, kuse,
     &                      xs, dvs, curnt, x, y, z, bx, by, bz)

c=======================================================================
c                                            By: M. Schaffer 2006 aug 15
c                                 Last Modified: M. Schaffer 2007 oct 25
c-----------------------------------------------------------------------
c POLYGONB calculates the magnetic field from a generic group of polygon
c loops (usually constituting a magnetic field source model) using the 
c BIOTLOOP subroutine. The model geometry is defined elsewhere, and is 
c embodied for POLYGONB in the input variables xs and dvs.
c
c Input arrays kuse,...,curnt, have adjustable dimensions.
c
c In the following: 
c    i=1,2, or 3 selects the x-, y-, or z-component of the point;
c    j is the index identifying a segment; 
c    k is the index identifying a loop.
c
c INPUTS:
c loopsdim = full array dimension for loops, maximum allowable k.
c segsdim  = full array dimension for segments, max allowable j.
c nloops   = integer, number of loops actually used in this group.
c nsegs(k) = integer array specifying number of polygon segments in
c            loop k of the group of loops.
c kuse(k)  = integer array used as a flag to use (.ne.0) or skip (.eq.0)
c            loop k of the group of loops.
c xs(i,j,k)= Cartesian position-defining  point i=1,2,3 = x,y,z of 
c            segment j of loop k.
c   The start point of segment j is the point (i=1,2,3) j of xs; its end
c   point is the point j+1 of xs, except end point of segment nsegs is
c   j=1 of xs.
c dvs(i,j,k), i=1,2,3,4, j=1...nsegs(k), are the precomputed Cartesian
c   direction vectors (3 elements) and lengths (4th element) of nsegs(k)
c   straignt segments.
c curnt(k) is the array of currents (amp) in polygon loop(k).
c x,y,z    = Cartesian coordinates of the field evaluation point
c
c OUTPUTS:
c bx,by,bz = Cartesian components of magnetic field at the field point
c
c Coordinate system: Geometry information is all written in right-
c   handed Cartesian coordinates, as used by BIOTLOOP.
c
c Units: SI (mks)
c Fully double precision.
c=======================================================================


      IMPLICIT NONE

!  Common variables:

!  Arguments:  Adjustable dimensions used here
      integer loopsdim, segsdim
      integer nloops, nsegs(loopsdim), kuse(loopsdim)
      real*8  xs(3,segsdim,loopsdim), dvs(4,segsdim,loopsdim)
      real*8  curnt(loopsdim)
      real*8  x, y, z, bx, by, bz 
      intent(in) :: loopsdim, segsdim, nloops, nsegs, kuse
      intent(in) :: xs, dvs, curnt, x, y, z
      intent(out):: bx, by, bz

!  Local Variables:
      real*8  xbk(3,segsdim), dbvk(4,segsdim)
      real*8  xc(3), bc(3)
      integer i, j, k

      SAVE

c-----------------------------------------------------------------------

C      WRITE(0,*) 'POLYGONB: dims1=', loopsdim, segsdim
C      WRITE(0,*) 'POLYGONB: nloops1=', nloops
C      WRITE(0,*) 'POLYGONB: nsegs1 =', nsegs

c Initialize magnetic field components to zero

      bx = 0.d0
      by = 0.d0
      bz = 0.d0
      
      xc(1) = x   	! BIOTLOOP input uses a 3-component Cartesian
      xc(2) = y   	! position vector, here xc, specifying the point 
      xc(3) = z  	! at which the magnetic field is calculated.

c BIOTLOOP calculates the magnetic field of a single current loop. It
c needs 2D arrays, here xbk and dbvk, specifying that polygon loop's 
c coordinate points and direction vectors. The external code calling 
c POLYGONB passes 3D arrays, here xs and dvs, in which the 3rd dimension
c is loop number. Here, we construct the necessary 2D arrays xbk and dbvk 
c and then call BIOTLOOP.


      DO k=1,nloops 				! do all the loops
       IF(kuse(k) .ne. 0) THEN			! flag to do this loop


			! NEW METHOOD using f90
       CALL BIOTLOOP (segsdim, nsegs(k), 
     &                xs(:,:,k), dvs(:,:,k), curnt(k), xc, bc)  


       bx = bx + bc(1)  	! BIOTLOOP returns the magnetic field
       by = by + bc(2)  	! as a 3-component Cartesian vector,
       bz = bz + bc(3)  	! here bc().
	
C       WRITE(15,*)'POLYGONB: bc(1, 2, 3): ', bc(1),bc(2),bc(3)

       END IF		! end of IF(kuse(k) .ne. 0)
      END DO


      RETURN
      
      END 	SUBROUTINE 	POLYGONB
c=======================================================================

