

c include file CCOILSD3.i

c                                            By: M. Schaffer 2003 feb 25
c                                 Last Modified: M. Schaffer 2004 apr 02
c    MS: varied Rc data for a C-coil geometry scan, for Todd 2007 may 26
c-----------------------------------------------------------------------
c Parameters for the D3D C-coils
c-----------------------------------------------------------------------
c The C-coil is represented as a set of 6 closed polygon current loops
c each comprised of straight-line segments.
c
c   ncloops= number of loops (nloops = 6 for DIII-D C-coil)
c   nctorsegs= number of segments in one toroidal side
c   mcpolsegs= number of segments in one poloidal side
c   ncsegs= number of segments in a single C-coil loop
c   xcs(k,j,i), i=1,2,3 j=1...ncsegs are (x,y,z) Cartesian position 
c      vectors of ncsegs coil-defining loop points, for k loops.
c
c Coordinate system: Coil geometry information is all written in 
c    right-handed Cartesian coordinates, as used by BD3CCOIL, BIOTLOOP.
c    
c Units: SI (mksA)
c
c Note: For convenience in relating to DIII-D machine "geographical" 
c   C-coil nomenclature and prior TRIP3D practice, the 6 coils are 
c   indexed in the order 1,...6  
c          <--->  C79, C139, C199, C259, C319, C19.
c   However, the coil toroidal angles here are all written in
c   right-handed cylindrical coordinates, so coordinate toroidal angle
c   negative of machine geographical toroidal angle.
c
c-----------------------------------------------------------------------
c
c Include file can have type, data, parameter declarations.
c Cannot assign values to subscripted variables in PARAMETER statements.
c
c Do not declare 'IMPLICIT NONE' in this include file; it will conflict
c  with 'IMPLICIT NONE' declarations in routines that use the variables
c  that are defined here.
c Explicitly declare here the types of all variables used here. They do
c  not have to be declared again in routines that INCLUDE 'CCOILSD3.I'
c
c-----------------------------------------------------------------------


c ======================================================================
c SIMPLE "PICTURE FRAME" C-COIL SET
c
c Cannot assign values to subscripted variables in parameter statements.
c Use mcpolsegs=1; extra pts on vertical runs are unnecessary, wasteful.
c-----------------------------------------------------------------------

      integer ncloops, nctorsegs, mcpolsegs, ncsegs
      integer ncturns  		!number of electrical turns in a C-coil

      parameter (ncloops= 6, nctorsegs= 4, mcpolsegs= 1)
      parameter (ncturns= 4)   
      parameter (ncsegs= 2*nctorsegs + 2*mcpolsegs)


      real*8 torarcc 			! toroidal span of one loop, deg
      data   torarcc /58.0/

      real*8 phictrc(ncloops)   !tor angles of ncloops loop centers, deg
      data   phictrc  / -079., -139., -199., -259., -319., -019. /

      real*8 Rc, Zcu, Zcl 		!radius, heights of loop arcs, m
      data   Rc, Zcu, Zcl /3.23, 0.8, -0.8/
C      data   Rc, Zcu, Zcl /2.40, 0.25, -0.25/	! for Todd's scan

c ======================================================================
