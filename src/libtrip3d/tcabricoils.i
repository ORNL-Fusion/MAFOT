

c include file TCABRICOILS.i

c-----------------------------------------------------------------------
c                                            By: D.M.Orlov    2015 jun 22
c                                 Last Modified: D.M.Orlov    2015 jun 22
c-----------------------------------------------------------------------
c Parameters for TCABR Internal ELM coils subroutines
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
      PARAMETER (mxIloops =  24)! max loops in a toroidal band
      PARAMETER	(mxIsegs  = 52)	! max polygon segments in a loop

      INTEGER niturns  	!number of electrical turns in an ELM-coil
      PARAMETER (niturns = 10) 		

c-----------------------------------------------------------------------
