 
 
c=======================================================================
 
      subroutine D3BUSNEW2GEOM (kbus, nloops, nsegs, xs, dvs, curnt)
 
c=======================================================================
c                                            By: M. Schaffer 2005 aug 15
c                                 Last Modified: M. Schaffer 2011 jun 05
c-----------------------------------------------------------------------
c D3BUSNEW2GEOM writes geometry parameter and current arrays for the new
c DIII-D B-coil buswork at 210 deg and updated 30 deg bus information.
c The 30 deg geometry is based on Big DEE drawings and on measurements
c from 2006 July, 2006 Nov, 2007 Oct.  A presumed ferromagnetic contri-
c bution from neutral beam support structure, fit originally with a 
c separate point dipole, is now added as an extra loop in D3BUSNEW2GEOM.
c The 210 deg geometry comes from the 2005 AUGUST Final Design Review
c (FDR) design, except:
c  1) The bends in the return belt bus are now each moved 1.5 inch
c     toward the 210 deg gap.
c  2) The quadrupole bus was not part of this FDR. The geometry 
c     presently in this subroutine is extrapolated from FDR figures.
c
c The buswork is represented as a set of nloops closed polygon current
c loops, each comprised of nsegs straight-line segments.
c
c This subroutine replaces the older D3BUS...GEOM routine, to conform
c to the new (2007 Jan) code architecture, wherein all DIII-D bus
c routines are grouped together and only interact locally, and xs and
c dvs array incides appear in order (i,j,k).
c
c INPUTS:
c  kbus(k), k=1...nloops is integer flag to tell which loops to activate.
c
c OUTPUTS:
c  nloops = number of polygon current loops used in this model
c  nsegs(k), k=1...nloops are numbers of segments in each buswork
c   loop. k is used consistently here as the index identifying a current
c   loop. Buswork model loops are not all the same size.
c  xs(i,j,k), i=1,2,3, j=1...nsegs, are the (x,y,z) Cartesian position 
c   vectors of nsegs(k) coil-defining loop points of one loop. j is
c   used consistently here as the index identifying a segment. i=1,2,
c   or 3 select the x-, y-, or z-component of the point.
c  dvs(i,j,k), i=1,2,3,4 j=1...nsegs are Cartesian direction vectors  
c   (3 elements) and lengths (4th element) of nsegs(k) straight segments.
c   The start point of segment j is row j of xs, its end point is 
c   row j+1 of xs, except end point of segment nsegs is row 1 
c   of xs.
c  curdist(k), k=1...nloops, specifies relative bus current in
c   loop k. In general, different loops comprising the buswork model
c   may have different currents, but each loop has only one current.
c   Example of curdist for nloops=4 is:
c         (1., 0.5, 1., 0.5) for (full,half,full,half), respectively. 
c  curnt(k) is the array of actual bus currents (amp) in loop k of the
c   various loops comprising the buswork model. Different loops in the
c   model may have different currents, per curdist(k).  Curnt is used
c   in the Biot-Savart calculation of B by BIOTLOOP. 
c
c Coordinate system: Bus geometry information is all written in 
c   right-handed Cartesian coordinates, as used by POLYGON, BIOTLOOP.
c
c This version orders array indices as (i,j,k). Previously (k,j,i).
c
c Units: SI (mks)
c Fully double precision.
c-----------------------------------------------------------------------
c  Parameter variable types are declared in PARAMETER stmt before using,
c  if also using IMPLICIT NONE:
 
      IMPLICIT NONE
 
      INCLUDE 'd3bus.i'		! included parameter & variable types
				! are declared in the include file.
      INTEGER    mxpgloops, mxpgsegs		! local parameters
      PARAMETER (mxpgloops = maxpolygloops)	! to dimension arrays
      PARAMETER (mxpgsegs  = maxpolygsegs) 	! to dimension arrays
c-----------------------------------------------------------------------
 
 
      common /consts/ pi, twopi, cir, rtd, dtr
      common /tfcoil/ tflimits,rma,bma,dsbtc,alfsbtc,dthbtc,alfthbtc,
     &                ntfturns,ibcoil,lbtc,iripple,bcoil,usebcoilmds
      
!  Common variables:
      real*8  pi, twopi, cir, rtd, dtr
      real*8  tflimits(4),rma, bma, dsbtc, alfsbtc, dthbtc, alfthbtc
      integer ntfturns, ibcoil, lbtc, iripple
      real*8  bcoil
      integer usebcoilmds
 
!  Arguments:
      integer kbus(mxpgloops), nloops, nsegs(mxpgloops)
      real*8  xs(3,mxpgsegs,mxpgloops), dvs(4,mxpgsegs,mxpgloops)
      real*8  curnt(mxpgloops)
      intent(in) :: kbus
      intent(out):: nloops, nsegs, xs, dvs, curnt
 
!  Local Variables:
      real*8  curdist(mxpgloops)
      real*8  phi, z, dumbr, bphi, dumbz, xpsi, curturn
      real*8  dpolA, dpolD, dpolX, dpolY, dpolZ	!mjs for ferromag model
      integer i, j, k, jend
 
      SAVE
 
c-----------------------------------------------------------------------
c            Get Toroidal Field Coil conductor current
c-----------------------------------------------------------------------
c TF-coil buswork current is proportional to the electric current in one
c  TF-coil turn. Therefore, we call BTFCOIL, which returns the coil's 
c  vacuum B at the reference coil radius, here called bphi and rma,
c  respectively; then use Ampere's law to get the total (AMP-turns)
c  current.  The single-turn current is the 'unit current' in the 
c  TF buswork, here called curturn. The physical current curnt(k) in 
c  buswork loop k is curturn times that loop's current distribution
c  factor, curdist(k).
c
cms 2008 Aug 05:
c Added a feature, to calculate the Bus error field ALONE that would be
c  generated by the DIII-D Bcoil, even when the user does not want a 
c  B-coil toroidal field to be present. This user intention is signaled
c  if the flag IBCOIL has a negative value.
C      write(0,*)'D3BUSNEW2GEOM:'
C      write(0,*)' ibcoil_1:', ibcoil

c      phi = 0.0d0	! a dummy toroidal angle for BTFCOIL routine
c      z   = 0.0d0	! a dummy elevation for BTFCOIL routine

c      if(ibcoil .ge. 0) then

c       CALL BTFCOIL(rma,phi,z,dumbr,bphi,dumbz,xpsi) !get vacuum Bt0r

cms 2008 Aug 06 New code to calculate TF bus error alone:
c     To accomodate calculation of TF-bus error, even when user has
c      specified vacuum Btor = 0:  This feature is flagged by negative
c      values of ibcoil. However, subroutine BTFCOIL needs ibcoil > 0.
c      Therefore, if ibcoil is negative, we temporarily change its sign
c      to positive, then call BTFCOIL to get a returned Bphi that we use
c      to calculate the B-coil current = the bus current. Then we return
c      ibcoil to its original sign.

c      else		! ibcoil .lt. 0; change sign in to call BTFCOIL
c       ibcoil = -ibcoil
C       write(0,*)' ibcoil_2:', ibcoil

c        CALL BTFCOIL(rma,phi,z,dumbr,bphi,dumbz,xpsi) !get vacuum Btor

C       write(0,*)' rma, bphi:', rma, bphi
c       ibcoil = -ibcoil	! return ibcoil to negative, for no vacuum Btor

C       write(0,*)' ibcoil_3:', ibcoil
C       write(0,*)' kbus:', kbus
c      endif
cms End of 2008 Aug 06 changes.


      curturn = bcoil       	! external input via common block
c      write(0,*)' curturn, bcoil:', curturn, bcoil
c      if (usebcoilmds.eq.0) then
c        curturn = 0.5d7 *rma*Bphi/ntfturns      ! Ampere's law applied
c        write(80,*) 'Bcoil current is given as 0.5d7*rma*Bphi/ntfturns='
c        write(80,*) curturn
c      else
c        curturn = bcoil					! bcoil from MDSplus
c        write(80,*) 'B-coil current is from MDSplus signal bcoil = '
c        write(80,*) curturn
c      end if
 
C      write(0,*)'d3busnew2geom: curturn, nloops=', curturn,nloops
 
 
c-----------------------------------------------------------------------
c            Fill Cartesian Point Position Vector Array
c-----------------------------------------------------------------------
c Tables of relative loop currents and loop-defining points xs(i,j,k)
cms 2008 Sept 05:  MODELS UPDATED WITH LATEST INFORMATION (from 2006 Oct)
c Fortran 90 implied-do constructor is used to assign values to this 
c  array, which appear in common.
 
 
      nloops = 11	!mjs 2011 Jun 06 increased to 11 polygonal loops
			! previous D3BUSNEW2GEOM bus model had 10 loops
 
 
c Fill current array with unit TF current for nloops loops
 
      DO 300 k=1,nloops
       if(kbus(k) .ne. 0) then
        curnt(k) = curturn
       else
        curnt(k) = 0.0d0
       endif
C       write(0,*) 'd3busnew2geom: k, curnt(k)=', k,curnt(k)
  300 CONTINUE
 
c 		BEGIN 210 deg BUS MODEL: 
c k=1  MISSING IDEAL LOOP at 210 machine deg
c 4 points, 4 segments, carries full current
c xs(i,j,k)  example statement:  xs(1:3,j,k) = (/+3.05,-1.63,-0.45/)
      k = 1
      nsegs(k)   = 4
      curdist(k) = 1.0d0 			!full current
      xs(1:3,1,k) = (/-2.800, +1.460, -0.454/)  !-208 deg, large R
      xs(1:3,2,k) = (/-2.660, +1.700, -0.454/)  !-212 deg, large R
      xs(1:3,3,k) = (/-2.390, +1.485, -0.454/)  !-212 deg, small R
      xs(1:3,4,k) = (/-2.480, +1.320, -0.454/)  !-208 deg, small R
c     there is an implicit radial segment outward to starting point
      curnt(k) = curnt(k)*curdist(k)
 
c k=2  B14 (-202.5 deg) RETURN BARS WITH VERTICAL JOG
c 6 points, 6 segments, carries full current
c xs(k,j,i)  example statement:  xs(2,1,1:3) = (/+2.37,-1.54,-0.37/)
      k = 2
      nsegs(k)   = 6
      curdist(k) = 1.0d0 			!full current
      xs(1:3,1,k) = (/-2.540, +1.345, -0.454/)	!above B14 tab
      xs(1:3,2,k) = (/-2.680, +1.405, -0.454/)	!went radially out
 
c FDR design, bend at middle of bar
c      xs(1:3,3,k) = (/-2.780, +1.150, -0.454/)	!went across B14 face
c      xs(1:3,4,k) = (/-2.780, +1.150, -0.483/)	!jogged down
c Final, bend 1.5 inch closer to 210 deg end than FDR figure, 2005aug21
      xs(1:3,3,k) = (/-2.770, +1.190, -0.454/)	!went across B14 face
      xs(1:3,4,k) = (/-2.770, +1.190, -0.483/)	!jogged down
c Modified FDR, bend 10 cm closer to 210 deg end
c      xs(1:3,3,k) = (/-2.750, +1.240, -0.454/)	!went across B14 face
c      xs(1:3,4,k) = (/-2.750, +1.240, -0.483/)	!down to tab level
 
      xs(1:3,5,k) = (/-2.680, +1.405, -0.483/)	!back across B14
      xs(1:3,6,k) = (/-2.540, +1.345, -0.483/)	!B14 tab point
c     there is an implicit vertical segment up to starting point
      curnt(k) = curnt(k)*curdist(k)
 
c k=3  B15 (-217.5 deg) RETURN BARS WITH VERTICAL JOG
c 6 points, 6 segments, carries full current
c xs(k,j,i)  example statement:  xs(3,1,1:3) = (/-3.65, +1.30, -0.78/)
      k = 3
      nsegs(k)   = 6
      curdist(k) = 1.0d0 			!full current
      xs(1:3,1,k) = (/-2.435, +1.525, -0.454/)	!above B15 tab
      xs(1:3,2,k) = (/-2.435, +1.525, -0.489/)	!B15 tab point
      xs(1:3,3,k) = (/-2.550, +1.615, -0.489/)	!went radially out
 
c FDR design, bend at middle of bar
c      xs(1:3,4,k) = (/-2.390, +1.830, -0.489/)	!went across B15 face
c      xs(1:3,5,k) = (/-2.390, +1.830, -0.454/)	!jogged up to ideal level
c Final, bend 1.5 inch closer to 210 deg end than FDR figure, 2005aug21
      xs(1:3,4,k) = (/-2.410, +1.800, -0.489/)	!went across B15 face
      xs(1:3,5,k) = (/-2.410, +1.800, -0.454/)	!jogged up to ideal level
c Modified FDR, bend 10 cm closer to 210 deg end
c      xs(1:3,4,k) = (/-2.450, +1.750, -0.489/)	!went across B15 face
c      xs(1:3,5,k) = (/-2.450, +1.750, -0.454/)	!up to ideal loop level
 
      xs(1:3,6,k) = (/-2.550, +1.615, -0.454/)	!back across B15
c     there is an implicit radial segment inward to starting point
      curnt(k) = curnt(k)*curdist(k)
 
c k=4  Loop from 210 deg FEED CONDUCTORS + INNER INTERCONNECT
c 10 points, 10 segments, carries full current
c xs(k,j,i)  example statement:  xs(4,1,1:3) = (/-2.86,+1.45,-0.45/)
      k = 4
      nsegs(k)   = 10
      curdist(k) = 1.0d0 			!full current
      xs(1:3,1,k) = (/-2.92, +1.600, -0.598/)	!feed & connect to quad bus)
      xs(1:3,2,k) = (/-2.72, +1.510, -0.598/)	!feed went radially in
      xs(1:3,3,k) = (/-2.72, +1.510, -0.483/)	!jogged up to tab level
      xs(1:3,4,k) = (/-2.54, +1.425, -0.483/)	!B14 tab, feed connect pt.
      xs(1:3,5,k) = (/-2.54, +1.345, -0.483/)	!B14 tab point
      xs(1:3,6,k) = (/-2.48, +1.320, -0.483/)	!innermost pt. near B14
      xs(1:3,7,k) = (/-2.39, +1.485, -0.489/)	!innermost pt. near B15
      xs(1:3,8,k) = (/-2.66, +1.700, -0.489/)	!return went radialy out
      xs(1:3,9,k) = (/-2.66, +1.700, -0.598/)	!return, jogged down 
      xs(1:3,10,k)= (/-2.80, +1.810, -0.598/)	!return out at bus radius
c     there is an implicit segment toroidally back to starting point
      curnt(k) = curnt(k)*curdist(k)
 
c k=5  210 deg QUADRUPOLE BUS from deep in pit to point (4,1), 2005sep13
c 12 points, 12 segments, carries HALF current
c xs(k,j,i)  example statement:  xs(4,1,1:3) = (/-2.86,+1.45,-0.45/)
      k = 5
      nsegs(k)   = 12
      curdist(k) = 0.5d0 			!HALF current
      xs(1:3,1,k) = (/-2.920, +1.600, -0.598/)	!end quad center bar
      xs(1:3,2,k) = (/-2.920, +1.600, -0.448/)	!jumped to end quad top bar
      xs(1:3,3,k) = (/-3.280, +0.980, -0.448/)	!top bar corner
      xs(1:3,4,k) = (/-3.280, +0.980, -3.400/)	!top bar end in pit
      xs(1:3,5,k) = (/-3.250, +1.040, -3.400/)	!jumped to center bar in pit
      xs(1:3,6,k) = (/-3.250, +1.040, -0.598/)	!center bar corner
      xs(1:3,7,k) = (/-2.920, +1.600, -0.598/)	!end quad center bar
      xs(1:3,8,k) = (/-2.920, +1.600, -0.748/)	!jumped to end quad bot bar
      xs(1:3,9,k) = (/-3.220, +1.100, -0.748/)	!bottom bar corner
      xs(1:3,10,k)= (/-3.220, +1.100, -3.400/)	!bottom bar end in pit
      xs(1:3,11,k)= (/-3.250, +1.040, -3.400/)	!jumped to center bar in pit
      xs(1:3,12,k)= (/-3.250, +1.040, -0.598/)	!center bar corner
c     there is an implicit segment back to starting point
      curnt(k) = curnt(k)*curdist(k)

 
c 		BEGIN 30 deg BUS MODEL: 
c k=6  30 deg(machine) RETURN BELT CROSSOVERS, 2006 July 21 model
c 4 points, 4 segments, carries HALF current
c xs(k,j,i)  model statement:  xs(6,1,1:3) = (/+2.370,-1.540,-0.37/)
      k = 6
      nsegs(k)   = 4
      curdist(k) = 0.5d0 			!HALF current
      xs(1:3,1,k) = (/+2.571, -1.360, -0.359/)   !j=1
      xs(1:3,2,k) = (/+2.464, -1.548, -0.547/)   !j=2
      xs(1:3,3,k) = (/+2.390, -1.487, -0.547/)   !j=2
      xs(1:3,4,k) = (/+2.485, -1.322, -0.359/)   !j=4
c     there is an implicit segment back to starting point
      curnt(k) = curnt(k)*curdist(k)


c k=7  30 deg(machine) VERTICAL BUS PAIR, 2006 July 21 model
c 4 points, 4 segments, carries full current
c xs(k,j,i) example statement:  xs(5,j,1:3) = (/+3.052,-1.631,-0.45/)
      k = 7
      nsegs(k)   = 4
      curdist(k) = 1.0d0			!full current
      xs(1:3,1,k) = (/+2.79 , -1.47 , -0.453/)	!missing outer belt piece
      xs(1:3,2,k) = (/+2.673, -1.672, -0.453/)	!other end of missing belt
      xs(1:3,3,k) = (/+2.673, -1.673, -2.61/)	!went down 1st vertical bus
      xs(1:3,4,k) = (/+2.79 , -1.47 , -2.61/)	!at bottom, jumper to 2nd bus
c     there is an implicit segment up along 2nd vert. bus to starting point
      curnt(k) = curnt(k)*curdist(k)


c k=8  30 deg(machine) HORIZONTAL LOWER LOOP
c 8 points, 8 segments, carries full current
c xs(k,j,i) example statement:  xs(5,j,1:3) = (/+3.052,-1.631,-0.45/)
      k = 8
      nsegs(k)   = 8
      curdist(k) = 1.0d0			!full current
      xs(1:3,1,k) = (/+2.673, -1.673, -2.61/)	!bottom of 1st vertical bus
      xs(1:3,2,k) = (/+3.415, -2.10 , -2.61/)	!went radially outward
      xs(1:3,3,k) = (/+3.75 , -2.00 , -2.61/)	!went out toward north wall
      xs(1:3,4,k) = (/+3.75 , -1.80 , -2.61/)	!went parallel to north wall
      xs(1:3,5,k) = (/+3.53 , -1.80 , -2.61/)	!jumper, outer-to-inner feed
      xs(1:3,6,k) = (/+3.40 , -1.61 , -2.61/)	!away from inner feed
      xs(1:3,7,k) = (/+3.028, -1.61 , -2.61/)	!went radially inward
      xs(1:3,8,k)= (/+2.79 , -1.47 , -2.61/)	!went radially in to 2nd bus
c     there is an implicit jumper back to starting point, bottom of 1st bus
      curnt(k) = curnt(k)*curdist(k)


c k=9  30 deg(machine) HORIZONTAL LOOP REPRESENTING FERROMAGNETIC SOURCES
c 4 points, 4 segments, carries full current
c xs(k,j,i) example statement:  xs(5,j,1:3) = (/+3.052,-1.631,-0.45/)
      k = 9
      nsegs(k)   = 4
      curdist(k) = -1.0d0		!OPPOSITE full current
c Define a square loop of area dpolA, such that full bus current
c in it equals the dipole fitted (2007 Oct 26) to B_ext measured at
c code phi = -8 and -15 deg, z = -.26, -.04, +.21 meter.
      dpolA = 23700.0/58871.0	! area of ferromagnetic effect dipole
      dpolD = 0.5*sqrt(dpolA)	! half-edge of square dipole loop
      dpolX = +3.595		! x,y,z coords of dipole loop Center
      dpolY = -1.753
      dpolZ = -2.600
c
      xs(1:3,1,k) = (/dpolX+dpolD, dpolY+dpolD, dpolZ/)
      xs(1:3,2,k) = (/dpolX-dpolD, dpolY+dpolD, dpolZ/)
      xs(1:3,3,k) = (/dpolX-dpolD, dpolY-dpolD, dpolZ/)
      xs(1:3,4,k) = (/dpolX+dpolD, dpolY-dpolD, dpolZ/)
c     there is an implicit jumper back to starting point, bottom of 1st bus
      curnt(k) = curnt(k)*curdist(k)


c k=10  30 deg(machine) LOWEST LOOPS, 2006 Nov + 2007 Oct information.
c TWO physical loops connected together at points j=1 & j=6.
c 10 points, 10 segments, carries full current; SIGN DOES NOT CHANGE!!!
      k = 10
      nsegs(k)   = 10
      curdist(k) = 1.0d0			!full current
c  wpper loop, in east-west plane at x = +3.64 m.
      xs(1:3,1,k) = (/+3.64, -1.80, -3.09/)	!+ feed pt. to E-W loop
      xs(1:3,2,k) = (/+3.64, -2.03, -3.09/)	!loop lower E corner
      xs(1:3,3,k) = (/+3.64, -2.03, -2.88/) 	!uppr E corner; switch bar
      xs(1:3,4,k) = (/+3.64, -1.54, -2.88/)	!uppr W corner; switch bar
      xs(1:3,5,k) = (/+3.64, -1.54, -3.28/)	!loop lowr W corner
      xs(1:3,6,k) = (/+3.64, -1.80, -3.28/)	!- feed pt. to E-W loop
c  lower loop, in north-south plane at y = -1.80 m.
      xs(1:3,7,k) = (/+3.84, -1.80, -3.28/)	!- feed corner
      xs(1:3,8,k) = (/+3.84, -1.80, -3.61/)	!- feed lowest corner
      xs(1:3,9,k) = (/+4.07, -1.80, -3.61/)	!connection to quadrupole
      xs(1:3,10,k)= (/+4.07, -1.80, -3.09/)	!+ feed upper N corner
c     there is an implicit segment back to starting point
      curnt(k) = abs(curnt(k)*curdist(k)) 	!CURRENT SIGN FIXED HERE
 
 
c k=11  165 deg(machine) MODIFIED BELT-BUS JUMPER (to clear off-axis NB)
c       For all shots taken in 2011 and later.
c 4 points, 4 segments, carries full current
c xs(i,j,k)  example statement:  xs(1:3,j,k) = (/+3.05,-1.63,-0.45/)
      k = 11
      nsegs(k)   = 4
      curdist(k) = 1.0d0 			!full current
      xs(1:3,1,k) = (/-3.022, -0.924, -0.454/)  !-163 deg, large R
      xs(1:3,2,k) = (/-3.079, -0.711, -0.454/)  !-167 deg, large R
      xs(1:3,3,k) = (/-2.982, -0.688, -0.454/)  !-167 deg, small R
      xs(1:3,4,k) = (/-2.926, -0.895, -0.454/)  !-163 deg, small R
c     there is an implicit radial segment outward to starting point
      curnt(k) = curnt(k)*curdist(k)
 
 

c k=NN  30 deg(machine) LOWEST LOOP, 2006 Nov information; superceded.
c 4 points, 4 segments, carries full current
c      k = 8
c      nsegs(k)   = 4
c      curdist(k) = 1.0d0			!full current
c      xs(1:3,1,k) = (/+3.55 , -1.45 , -2.88/) 	!top west
c      xs(1:3,2,k) = (/+3.55 , -1.95 , -2.88/) 	!top east
c      xs(1:3,3,k) = (/+3.55 , -1.95 , -3.10/) 	!bottom east
c      xs(1:3,4,k) = (/+3.55 , -1.45 , -3.28/) 	!bottom west
c      curnt(k) = -abs(curnt(k)*curdist(k)) 	!current ALWAYS CCW here
 
 
 
C To write points to fort.11:
C      write(11,*) ' '
C      write(11,*) '***************   D3BUSNEW2GEOM   ***************'
C      write(11,*) 'nloops:', nloops
C      write(11,*) 'nsegs: ',
C     &             nsegs(1),nsegs(2),nsegs(3),nsegs(4)
C	
C      DO 900 k=1,nloops			! do all model loops
C       write(11,*) ' '
C       write(11,*) 'loop #', k
C       write(11,*) 
C     &   ' Xpoint(m)               Ypoint(m)               Zpoint(m)'
C       DO 902 j=1,nsegs(k)
C        write(11,*) xs(k,j,1),xs(k,j,2),xs(k,j,3)
C  902  CONTINUE
C  900 CONTINUE
 
 
c-----------------------------------------------------------------------
c            Fill Direction Vector Array
c-----------------------------------------------------------------------
 
      DO 200 k=1,nloops
c      write(0,*) 'd3busnew2geom: k, curnt(k)=', k,curnt(k)
       DO 202 j=1,nsegs(k)
       
        if (j.ne.nsegs(k)) then  
         jend = j+1              ! end point for all segments but last
        else
         jend = 1                ! end point for last segment
        endif
 
        ! segment length vector:
        dvs(1,j,k) = xs(1,jend,k) - xs(1,j,k)
        dvs(2,j,k) = xs(2,jend,k) - xs(2,j,k)
        dvs(3,j,k) = xs(3,jend,k) - xs(3,j,k)
        
        ! segment length:
        dvs(4,j,k) = sqrt(dvs(1,j,k)**2 + dvs(2,j,k)**2 + dvs(3,j,k)**2)
 
        ! segment direction vector:
         DO 203 i=1,3
          dvs(i,j,k) = dvs(i,j,k)/dvs(4,j,k)
  203    CONTINUE
 
  202  CONTINUE
 
  200 CONTINUE
 
 
 
C      WRITE(15,*) 'd3busnew2geom: Bphi, rma, curturn', bphi,rma,curturn
C      WRITE(15,*) 'd3busnew2geom: curnt(1)', curnt(1)
 
      RETURN
 
      END	SUBROUTINE 	D3BUSNEW2GEOM
c=======================================================================
 
 
