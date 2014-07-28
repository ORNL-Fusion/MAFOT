

c include file MASTICOILS.i

c-----------------------------------------------------------------------
c                                            By: D. Orlov 2009 jul 21
c                                 Last Modified: D. Orlov 2009 oct 13
c-----------------------------------------------------------------------
c Parameters for MAST Internal-coil (I-coil) subroutines
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
      PARAMETER (mxIbands =  20) 	! max # toroidal bands of loops
      PARAMETER (mxIloops =  60) 	! max loops in a toroidal band
      PARAMETER	(mxIsegs  = 100) 	! max polygon segments in a loop

      INTEGER niturns  	!number of electrical turns in an I-coil
      PARAMETER (niturns = 4)   

c-----------------------------------------------------------------------
