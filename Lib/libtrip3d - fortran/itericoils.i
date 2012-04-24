

c include file ITERECCOILS.i

c-----------------------------------------------------------------------
c                                            By: M. Schaffer 2006 dec 22
c                                 Last Modified: D. Orlov    2010 jul 28
c-----------------------------------------------------------------------
c Parameters for ITER Internal RMP coils (I-coil) subroutines
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

      INTEGER   mxIbands, mxIloops, mxIsegs	! set array dimensions
      PARAMETER (mxIbands =  3)	! max # toroidal bands of loops
      PARAMETER (mxIloops =  9)	! max loops in a toroidal band
      PARAMETER	(mxIsegs  = 52)	! max polygon segments in a loop

      INTEGER niturns  	!number of electrical turns in an EC-coil
      PARAMETER (niturns = 1) 		!placeholder; ITER not yet built

c-----------------------------------------------------------------------
