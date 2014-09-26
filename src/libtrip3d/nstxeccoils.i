

c include file NSTXECCOILS.i

c-----------------------------------------------------------------------
c                                            By: D. Orlov 2009 jun 12
c
c-----------------------------------------------------------------------
c Parameters for NSTX Error Correction-coil (EC-coil) subroutines
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
      PARAMETER (mxECbands =  1)	! max # toroidal bands of loops
      PARAMETER (mxECloops =  6)	! max loops in a toroidal band
      PARAMETER	(mxECsegs  = 26)	! max polygon segments in a loop

      INTEGER ncturns  	!number of electrical turns in an EC-coil
      PARAMETER (ncturns = 2)   

c-----------------------------------------------------------------------
