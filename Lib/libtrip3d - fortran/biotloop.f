

c=======================================================================

      subroutine BIOTLOOP (narray, ns, xc, dv, cc, x, b)

c=======================================================================
c                                Started by:     MJSchaffer 2002 feb 21
c                                Last Modified:  MJSchaffer 2007 oct 25
c----------------------------------------------------------------------
c BIOTLOOP computes the magnetic field 3-vector b(3) of an arbitrary
c closed polygon current loop comprised of straight-line segments. 
c The numerical method of Hanson & Hirshman, Phys. Plasmas vol.9 (2002)
c p.4410 is used for efficiency and because it has no problem when the 
c field point lies on the straight line extension of a segment.
c BIOTLOOP is written to be very general, making no assumptions about 
c the nature of the loop.
c BIOTLOOP fails if the field point lies so close to a segment as to 
c cause numerical overflow.
c
c Input arrays xc,dv,cc, have adjustable dimensions.
c
c Coordinate system: Cartesian
c
c Units: SI (mksA)
c
c Inputs: 
c   narray= Size of the 2nd or 'segments' dimension of arrays xc(i,j) 
c      and dv(i,j). This may be larger than ns (see below) when the 
c      actual number of loop segments is less than the array size. 
c      This information lets BIOTLOOP handle loop specification arrays
c      that have empty elements. Only the first ns rows of these arrays
c      are read and must contain the loop data. Thus, ns < = narray.
c   ns= Number of actual segments in the given loop.
c
c   xc(i,j), i=1,2,3 j=1...narray are (x,y,z) Cartesian position 
c      vectors. There must be j=1,ns path-defining loop-point 
c      position vectors.
c   dv(i,j), i=1,2,3,4 j=1...narray are Cartesian direction vectors  
c      (3 elements) and lengths (4th element). There must be j=1,ns
c      path-defining loop-point position vectors.
c      contain ns segment direction vectors and lengths. 
c      The start point of segment j is row j of xc, its end point is 
c      row j+1 of xc, except the end point of segment ns is row 1 of xc.
c   cc= Electric current (amp) in the closed polygon loop.
c   x(i), i=1,2,3 is (x,y,z) Cartesian position vector of field point.
c
c Outputs:
c   b(i), i=1,2,3 is (Bx,By,Bz) Cartesian magnetic induction vector 
c      (tesla) of the closed current polygon at x.
c
c Note:  Positive electric current is assigned to the direction of 
c   increasing segment and path-defining point identification indices, 
c   j, from j=1 to j=ns.
c
c Note: Although segments of zero length, i.e., dv(j,4)=0, are skipped
c   if detected, it is best practice to delete unused elements in the 
c   range 1 <= j <= ns from the loop specification files xc and dv.
c
c The coding nomenclature follows Hanson & Hirshman, except the use of 
c dv for their L, the length of a segment.
c The programming style and structure are copied from Wayne D. Bard's 
c subroutine BSEGS, developed by him at least as early as 1979.
c
c Fully double precision.
c-----------------------------------------------------------------------

      IMPLICIT NONE

c-----------------------------------------------------------------------

!  Arguments
      integer   narray, ns
      real*8    xc(3,narray), dv(4,narray), cc, x(3), b(3)
      intent(in)  :: narray, ns
      intent(in)  :: xc, dv, cc, x
      intent(out) :: b

!  Local Variables

      real*8    dp0, Ri(3), Rf(3), cross(3)
      real*8	Risq, Rfsq, Rimag, Rfmag, Rsum
      real*8    SIunit, dbf, dbx, dby, dbz, fac
      integer   i, j, jend

      SAVE

!-----------------------------------------------------------------------

      data SIunit/ 2.0d-7 /       ! mu_0/(2 pi)

c-----------------------------------------------------------------------

c Initialize
      dp0  = 0.0d0

      b(1) = dp0
      b(2) = dp0
      b(3) = dp0
      
      if (cc .eq. dp0) RETURN		! no current in this loop


c Sum contributions of all segments
      DO j=1,ns
       IF(dv(4,j) .ne. dp0) THEN	! non trivial segment, use it

        if (j .ne. ns) then  
         jend = j+1	! segment end point for all segments but last
        else
         jend = 1	! segment end point for last segment
        endif

        Ri(1) = x(1) - xc(1,j)      	! start point of segment j
        Ri(2) = x(2) - xc(2,j)
        Ri(3) = x(3) - xc(3,j)

        Rf(1) = x(1) - xc(1,jend)	! end point of segment j
        Rf(2) = x(2) - xc(2,jend)
        Rf(3) = x(3) - xc(3,jend)

        Risq  = Ri(1)**2 + Ri(2)**2 + Ri(3)**2
        Rfsq  = Rf(1)**2 + Rf(2)**2 + Rf(3)**2

        Rimag = sqrt(Risq)          ! magnitudes of Ri & Rf
        Rfmag = sqrt(Rfsq)

        Rsum  = Rimag + Rfmag

c   calculate cross product between segment direction vector and Ri
        cross(1) = dv(2,j)*Ri(3) - dv(3,j)*Ri(2)
        cross(2) = dv(3,j)*Ri(1) - dv(1,j)*Ri(3)
        cross(3) = dv(1,j)*Ri(2) - dv(2,j)*Ri(1)

c   calculate Hansen & Hirshman's factor
        dbf = (dv(4,j)*Rsum/Rimag/Rfmag)/(Rsum**2 - dv(4,j)**2)

c   the "geometric B" (the vector components of B without the common
c   scale factors) of segment j is:
        dbx = cross(1)*dbf
        dby = cross(2)*dbf
        dbz = cross(3)*dbf

c   sum to get total "geometric B"
        b(1) = b(1) + dbx
        b(2) = b(2) + dby
        b(3) = b(3) + dbz

       END IF		! end of IF(dv(4,j) .ne. dp0)
      END DO

CC      WRITE(0,*) 'BIOT: b(1)=',b(1)		! for debug
CC      WRITE(0,*) 'BIOT: b(2)=',b(2)
CC      WRITE(0,*) 'BIOT: b(3)=',b(3)

c Convert B to physical unit
      fac = SIunit*cc        ! converts "geometric B" to tesla
      b(1) = b(1)*fac
      b(2) = b(2)*fac
      b(3) = b(3)*fac


      RETURN

      END 	SUBROUTINE 	BIOTLOOP
c=======================================================================
