

c include file D3BUS.i

c-----------------------------------------------------------------------
c                                            By: M. Schaffer 2006 aug 15
c                                 Last Modified: M. Schaffer 2006 dec 18
c-----------------------------------------------------------------------
c General parameters to set array sizes in d3bus subroutines that work
c  with polygon loop coils. 
c This include file is like the original BBUSD3.i file
c-----------------------------------------------------------------------
c Include file can have type, data, parameter declarations.
c
c Cannot assign values to subscripted variables in PARAMETER statements.
c
c Do not declare 'IMPLICIT NONE' in this include file; it will conflict
c with 'IMPLICIT NONE' declarations in routines that use the variables
c that are defined here.
c Explicitly declare here the types of all variables used here.
c-----------------------------------------------------------------------

      INTEGER maxpolygloops, maxpolygsegs

      PARAMETER (maxpolygloops =  30)	! array dimension and maximum 
					! number of polygon loops 
					! allowed per model.
      PARAMETER (maxpolygsegs  =  30)	! array dimension and max number
					! of polygon segments allowed 
					! allowed per loop.

c ======================================================================
