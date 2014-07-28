

c include file D3PFERR.i

c                                            By: M. Schaffer 2006 jun 10
c                                 Last Modified: M. Schaffer 2006 jun 12
c-----------------------------------------------------------------------
c Parameters for D3D poloidal field coils (F-coils)
c
c Parameters nfc and nflps are used to dimension arrays, etc.
c  in the subroutines used to calculate magnetic fields
c  from shifted and tilted poloidal field coils.
c
c D3D has 18 poloidal field shaping coils (F-coils).
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
c  not have to be declared again in routines that INCLUDE 'D3PFERR.i'
c
c-----------------------------------------------------------------------

      integer nfc, nflps

      parameter (nfc    = 18) 		!number of F-coils
      parameter (nflps  = 2*nfc)	!number of loops in model

