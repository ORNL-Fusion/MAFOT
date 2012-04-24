

c=======================================================================

      subroutine ELLIPINTS(ksq, Kint, Eint, S, Dk4, err)

c-----------------------------------------------------------------------
c                                    Started By: M. Schaffer 2006 feb 27
c                                 Last Modified: M. Schaffer 2007 oct 25
c-----------------------------------------------------------------------
c ELLIPINTS calculates the complete elliptic integrals of the first and
c second kinds as functions of argument ksq. The integrals are defined as
c
c              pi/2
c   K(k^2) = Integral [1/sqrt(1 - (k^2 * sin^2 phi))] dphi  (1st kind)
c               0
c
c              pi/2
c   E(k^2) = Integral [sqrt(1 - (k^2 * sin^2 phi))] dphi    (2nd kind)
c               0
c
c where k^2 = k*k is also known as (sin(alpha))^2, p and m in the 
c literature (refs 1 & 2).
c
c ELLIPINTS also returns difference functions, here called S and Dk4,
c
c   S   = [2(K - E) - ksq*K]/K = Sum of [ (2^n)*(c_n)^2 ],  (n=1,2,...)
c          where c_n are terms calculated during the iterative 
c          calculation of K and E
c   Dk4 = [2(K - E) - ksq*K]/k^4 = K*S/k^4
c
c which are very useful for numerical calculation of the magnetic field
c of a current loop. These extensions are documented in notes by
c M. Schaffer dated 2006 Apr 4, 6, 7. NOTE that S is the sum denoted by
c S in most of Schaffer's notes, but also by Dint. Unfortunately, both 
c 'sum' and 'dint' are fortran intrinsic functions, so they shouldn't be
c used as variables.
c
c INPUT --
c   ksq  = argument k^2 in the integrand; 0 <= ksq <= 1
c	   other notations:  k^2 = m = p;  k = sin(alpha)
c
c OUTPUT --
c   Kint  = value of K(k^2)
c   Eint  = value of E(k^2)
c   S     = value of S (k^2)
c   Dk4   = value of Dk4(k^2)
c   err   = 0 if no error
c         = 1 if ksq is out of range
c         = 2 if there was runaway in iteration loop
c
c K(k^2) is calculated by the convergence of the 'arithmetic-geometric'
c or 'common' mean upon their repeated evaluation.
c E(k^2) is calculated by summing a series of terms constructed at
c each step of the common-mean iteration. See (refs 1,2).
c
c The algorithms of the two references are identical for K but differ
c for E. Both algorithms for K converge to double precision numerical 
c accuracy limit (about 15 decimal digits) in .le. 10 iterations when
c the argument is at numerical limit, ksq = 0.99999 99999 99999 d00. 
c The computed K(ksq) stops converging to its ksq --> 1 limit,
c K(ksq) --> (1/2)*ln(16/(1-ksq)) at about ksq = 15 nines. The ref 1
c algorithm is slightly more accurate and is better behaved in this 
c limit, so it is the one used in this code. 
c
c This code iterates the common mean to convergence using Abramowitz &
c Stegun's method.
c  It is a modification of code by J.A. Leuer. The main change is
c  in the way the terms in the series for E(ksq) evaluation are
c  calculated, to avoid several unnecessary multiplications.
c  Another change here that the primary iteration control is arrival
c  at the numerical resolution limit. My tests of the iteration show
c  that variable 'term', which is the square of a numerical difference,
c  adds numerical noise to the iteration after reaching the
c  resolution limit, which makes subsequent computed 'term' values
c  increase. Therefore, maximum accuracy is attained from the last 
c  iteration for which the value of 'term' decreased. Note that 'term'
c  often reaches exact numerical zero, from which it never increases,
c  so the convergence test must consider both possibilities.
c  An absolute convergence number, 'tolitr' may be implemented, 
c  particularly if one wants to shorten the iteration as ksq --> 1
c  and can accept less accuracy. However, since convergence is fast,
c  there is probably no need for this option.
c
c
c ref 1. M. Abramowitz & I.A. Stegun, Handbook of Mathemetical Functions,
c        (Dover, 1965) pp 598-9.
c ref 2. J. Spanier & K.B. Oldham, Atlas of Functions, (Springer, 1987)
c        pp 613-4.
c
c Fully double precision.
c-----------------------------------------------------------------------

      IMPLICIT NONE

!  Arguments:
      integer err
      real*8  ksq, Kint, Eint, S, Dk4

!  Local Variables:

      real*8, PARAMETER:: pi = 3.14159265358979324d0, zero= 0.0d0
      real*8, PARAMETER:: one = 1.0d0, two = 2.0d0, thre= 3.0d0
      real*8, PARAMETER:: four= 4.0d0, fiv = 5.0d0, six = 6.0d0
      real*8, PARAMETER:: svn = 7.0d0, eght= 8.0d0, nin = 9.0d0
      real*8, PARAMETER:: half = one/two,  third= one/thre 
      real*8, PARAMETER:: forth= one/four, fifth= one/fiv
      real*8, PARAMETER:: hpi= pi/two, pi16th= pi/(two*eght)
      real*8, PARAMETER:: threfor = thre/four, threforsq= (thre/four)**2
      real*8, PARAMETER:: fivsixsq = (fiv/six)**2
      real*8, PARAMETER:: svneghtsq= (svn/eght)**2
      real*8, PARAMETER:: df2= 25.d0/32.d0, df3= 49.d0/60.d0
      real*8, PARAMETER:: df4= 27.d0/32.d0

      real*8, PARAMETER:: tolitr= 1.0d-33

      integer i, n
      real*8  ao, an, bo, bn, tn, summ, oldterm, term
      real*8  gsq, g, a, t, e, ssq
      real*8  Kser, Eser, Dser

      SAVE

c------------------------------------------------------------------------
c tolitr is convergence tolerance, 1.d-33 to 1.d-32 for double precision
c or it can be made larger for shorter iteration if less accuracy is
c acceptable.
c------------------------------------------------------------------------


c check that argument is in allowed range
      err = 0
      if((ksq.lt.zero) .or. (ksq.gt.one)) then
       Kint = 0.0d0
       Eint = 0.0d0
       err = 1
       RETURN
      end if

c values at limits

      if(ksq.eq.zero) then		! case of ksq = 0
       Kint = hpi
       Eint = hpi
       S    = zero
       Dk4  = pi16th			! pi/16
C      WRITE(0,*) 'ELLIPINTS K, E:', Kint, Eint		! DEBUG
       RETURN

      else if(ksq.eq.one) then 		! case of ksq = 1
       Kint = 1.0d99			! large, to approximate infinity
       Eint = one
       S    = one
       Dk4  = Kint
C      WRITE(0,*) 'ELLIPINTS K, E:', Kint, Eint		! DEBUG
       RETURN
      end if

c Main calculation, using Abramowitz & Stegun's iterative method
c initialize variables for iteration
      ao  = one			! for arithmetic mean
      bo  = sqrt(one - ksq)	! for geometric mean
				! Difference (one-ksq) is a source   
				! of error as ksq --> 1.
      tn   = one  		! initialize tn
      summ = zero		! initialize sum
      		! Don't use 'sum'; it's a fortran intrinsic function
      oldterm = one  		! initialize; used in convergence test

c		Iteration loop
      DO n=1,21			! using Fortran 90/95 WHILE syntax
       an  = half*(ao + bo) 		! new arithmetic mean
       bn  = sqrt(ao * bo)		! new geometric mean
       term= tn*(ao - bo)**2		! term = 2*[2^n (c_n)^2]

       if(term.ge.oldterm) EXIT  	! converged to best numerical
       					! resolution, exit DO loop
       oldterm = term

       summ = summ + term			! update series for E(ksq)
       tn  = two*tn			! update tn for next iteration
       ao  = an				! store new arith. mean as old
       bo  = bn				! store new geom. mean as old
C       WRITE(0,*) 'term=', term   	! *** DEBUG
       if(term.eq.zero) EXIT		! converged, can exit
C       if(term.lt.tolitr) EXIT		! alternate convergence criterion
      END DO

      if(n .gt. 20) then		! did not converge
        err = 2  		
        Kint = hpi			! return K(0) = pi/2 
        Eint = hpi			! return E(0) = pi/2
	RETURN
      end if

      Kint = hpi/ao			! K(ksq)
      S    = half*summ			! sum above has extra factor 2
      Eint = Kint*(one - half*(ksq + S))	! E(ksq)

C      WRITE(0,*) 'iterations=', n			! DEBUG
C      WRITE(0,903) 'ELLIPINTS:', Kint, Eint, S 	! DEBUG


c S varies as k^4 as k^2 --> 0, and S must be divided by k^4 to
c calculate the radial component of the magnetic field of a circular
c current loop. Here we calculate Dk4 = K*S/k^4 as a power series
c for ksq << 1 to avoid this numerical problem and make it possible to
c calculate B accurately when ksq is small.

      Dk4  = pi16th*(one + threfor*ksq*(one + df2*ksq*
     &              (one +    df3*ksq*(one + df4*ksq))) )


  901 format(1x, 5f15.9)	! some formats useful for
  902 format(1x, a, 3f20.15)	! test and debug
  903 format(1x, a, 1p3e22.15)

      RETURN
      END


c------------------    END OF SUBROUTINE ELLIPINTS    ------------------

C The following power series were useful to check the calculation of 
C K, E, and DK4 for ksq << 1. These series are written as I derived
C them. They are intended just to test the series and are not optimized
C for numerical efficiency.
C
C      Kser = hpi*(one +half**2*ksq*(one +(thre/four)**2*ksq*
C     &           (one +(fiv/six)**2*ksq*(one +(sevn/eght)**2*ksq))) )
C      Eser = hpi*(one -half**2*ksq*(one +(thre/four)**2*ksq*
C     &  (third +(fiv/six)**2*ksq*(fifth +(sevn/eght)**2*ksq/sevn))) )
C      Dk4  = (hpi/eght)*(one +(thre/four)*ksq*(one +(25.d0/32.d0)*ksq*
C     &          (one +(49.d0/60.d0)*ksq*(one +(27.d0/32.d0)*ksq))) )
C      Dser = Dk4*ksq*ksq/Kser


C These are more efficient versions of the above:
C
C      Kser = hpi*(one + forth*ksq*(one + threforsq*ksq*
C     &           (one + fivsixsq*ksq*(one + svneghtsq**2*ksq))) )
C      Eser = hpi*(one - forth*ksq*(one + threforsq*ksq*
C     &      (third + fivsixsq*ksq*(fifth + svneghtsq*ksq/sevn))) )
C      Dk4  = pi16th*(one + threfor*ksq*(one + df2*ksq*
C     &              (one +     df3*ksq*(one + df4*ksq))) )
C      Dser = Dk4*ksq*ksq/Kser


CC       WRITE(0,903) 'ELSERIES: ', Kser, Eser, Dser	! DEBUG
CC       WRITE(0,903) 'ELSERIES: ', S - Dser, Dk4	! DEBUG

c------------------    END OF SUBROUTINE ELLIPINTS    ------------------
