

c include file TCABRCPCOILS.i

c-----------------------------------------------------------------------
c                                            By: D.M.Orlov 2015 jun 22
c                                 Last Modified: D.M.Orlov 2015 jun 22
c-----------------------------------------------------------------------
c Parameters for TCABR Centerpost-coil (CP-coil) subroutines
c-----------------------------------------------------------------------
c Include file can have type, data, parameter declarations.
c
c Cannot assign values to subscripted variables in PARAMETER statements.
c
c Do not declare 'IMPLICIT NONE' in this include file; it will conflict
c with 'IMPLICIT NONE' declarations in routines that use the variables
c that are defined here.
c
c Explicitly declare here the types of all variables used in this file. 
c-----------------------------------------------------------------------

      INTEGER   mxECbands, mxECloops, mxECsegs	! set array dimensions
      PARAMETER (mxECbands =  3)	! max # toroidal bands of loops
      PARAMETER (mxECloops =  24)	! max loops in a toroidal band
      PARAMETER	(mxECsegs  = 52)	! max polygon segments in a loop

      INTEGER ncturns  	!number of electrical turns in an EC-coil
      PARAMETER (ncturns = 10)   

c-----------------------------------------------------------------------
