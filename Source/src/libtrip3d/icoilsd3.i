

c include file ICOILSD3.i

c                                            By: M. Schaffer 2003 mar 03
c                                 Last Modified: M. Schaffer 2007 jan 23
c-----------------------------------------------------------------------
c Parameters for scaled-down ITER I-coils in an otherwise DIII-D model.
c This include file uses the same variables as the DIII-D I-coil model.
c  In order to calculate with ITER I-coils, rename this file ICOILSD3.i
c  and recompile the SURFMN, TRIP3D, and/or PROBE codes as needed. Then
c  run codes with this new I-coil.
c
c ======================================================================
c SIMPLE "PICTURE FRAME" I-COIL SET
c
c Cannot assign values to subscripted variables in parameter statements.
c Use mipolsegs=1; extra pts on vertical runs are unnecessary, wasteful.
c NOTE:  IGEOM can presently do ONLY mipolsegs=1  !!!
c-----------------------------------------------------------------------

      integer niloops, nuploops, nloloops
      integer nitorsegs, mipolsegs
      integer niturns           !number of electrical turns in an I-coil
      integer nisegs

      parameter (niloops= 12, nuploops= 6, nloloops= 6)
      parameter (nitorsegs= 6, mipolsegs= 1)
      parameter (niturns= 1)
      parameter (nisegs= 2*nitorsegs + 2*mipolsegs)


c SPECIFIC PARAMETERS FOR DIII-D: * * * * * * * * * * * * * * * * * *  *

      real*8 torarci  			! toroidal span of one loop, deg
      data   torarci /51.72/

      real*8 phictri(niloops)   !tor angles of niloops loop centers, deg
      data   phictri /-032.7,-087.3,-152.7,-207.3,-272.7,-327.3,
     &                -032.7,-087.3,-152.7,-207.3,-272.7,-327.3 /

c Original 2003 Mar geometric parameters:
      real*8 Riu1, Ril1, Ziu1, Zil1   !radii, heights of upper loop arcs
      data   Riu1, Ril1, Ziu1, Zil1  /2.184, 2.394, +1.012, +0.504/
      real*8 Riu2, Ril2, Ziu2, Zil2   !radii, heights of lower loop arcs
      data   Riu2, Ril2, Ziu2, Zil2  /2.394, 2.184, -0.504, -1.012/

c Revised 2006 Dec parameters from 4:1 "large" drawing:
C      real*8 Riu1, Ril1, Ziu1, Zil1   !radii, heights of upper loop arcs
C      data   Riu1, Ril1, Ziu1, Zil1  /2.164, 2.373, +1.016, +0.504/
C      real*8 Riu2, Ril2, Ziu2, Zil2   !radii, heights of lower loop arcs
C      data   Riu2, Ril2, Ziu2, Zil2  /2.373, 2.164, -0.504, -1.016/

C Specific parameters for DIII-D:
C      parameter (torarci = 51.72)   !toroidal extent of a loop, deg
C      parameter (Riu1= 2.184, Ril1= 2.394, Riu2= 2.394, Ril2= 2.184) !m
C      parameter (Ziu1= 1.012, Zil1= 0.504, Ziu2=-0.504, Zil2=-1.012) !m
c ======================================================================
