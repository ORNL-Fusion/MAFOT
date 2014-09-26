

c=======================================================================

      subroutine CIRCLEB(a, cur, p, B)

c-----------------------------------------------------------------------
c                                    Started By: M. Schaffer 2006 feb 23
c                                 Last Modified: M. Schaffer 2007 oct 25
c-----------------------------------------------------------------------
c Given a point p = (x,y,z) in a right-handed Cartesean coordinate
c system, CIRCLEB calculates the magnetic field vector B of a single 
c circular current filament of radius = a, centered on the origin 
c and lying in the xy plane, carryng a current = cur amperes  
c counter-clockwise around the z-axis. 
c
c INPUT --
c  a   = radius (m) of circular current filament, real
c  cur = electric current (A) in loop, real
c  p   = CARTESIAN coordinates (m) of the point at which B-field is 
c        to be calculated, as a 3-element real array, p = (x,y,z)
c
c OUTPUT --
c  B   = magnetic field (T) at p, as a 3-element real array, (Bx,By,Bz)
c
c
c In this subroutine, R = sqrt(x^2 + y^2) is the cylindrical radius
c coordinate of field point P from the z-axis.
c
c Kint(ksq) and Eint(ksq) are complete elliptic integrals of 1st and
c 2nd kind, respectively, calculated here by subroutine ELLIPINTS,
c and ksq = k^2 = 4aR/[(a+R)^2 + z^2] for calculation of loop's B.
c The range of ksq is [0,1).
c
c The expression used to calculate B numerically is rearranged from
c the common one (e.g. J.A. Stratton, Electromagnetic Theory, 1941)
c that is used in the original TRIP3D code, in order to avoid subtrac-
c tion of large numbers in the limit of large distance from the loop, 
c which greatly compromises numerical accuracy. The expressions, here, 
c developed by M. Schaffer, notes dated 2006 Apr 4, 6, 7, give much 
c better accuracy and are used heres. They are:
c
c Bx= (2mu/pi)(cur*a^2)(x*z)(Kint/(f3*f1^(3/2)))*(1 -2(f2/f1)(S/k^4))
c By= replace x by y in the Bx expression
c Bz= (mu/2pi)cur(Kint/(sqrtf1*f3))*(2*(a^2 - aR) - f4*(k^2 + S)/2)
c
c S is a function of ksq and (Kint - Eint) and is used to eliminate
c the explicit appearance of Eint in expressions for B. Subroutine
c ELLIPINTS, that calculates Kint(ksq), Eint(ksq) and S(ksq), was
c specially written for this purpose. To further improve numerical
c accuracy of B,  When ksq is smaller than about 0.001 (this value
c for double precision with a 15 or 16 decimal digit mantissa), the
c direct calculation of S is replaced by a power series in k^2.
c ELLIPINTS returns the power series result as a variable called Dk4,
c and S = ((k^4)/Kint)Dk4 is calculated in CIRCLEB. In some of Schaffer's
c notes, S is called Dint. However, dint is a fortran intrinsic fortran
c function and should not be used as a variable name here.
c
c These techniques return B accurate to 15 decimal digits almost
c everywhere in (R/a, z/a) space, degrading possibly to 14 decimal
c digits near ksq = 0.001.   - M.J. Schaffer, 2006 Apr 10.
c
c Units: SI (mksA)
c Fully double precision.
c-----------------------------------------------------------------------

      IMPLICIT NONE
 
!  Common variables:

!  Arguments:
      real*8  a, cur, p(3), B(3)

!  Local Variables:

      real*8, PARAMETER:: pi    = 3.14159265358979324d0
      real*8, PARAMETER:: zero= 0.0d0, one= 1.0d0, two= 2.0d0
      real*8, PARAMETER:: four= 4.0d0, half = one/two
      real*8, PARAMETER:: mu    = pi*4.0d-7		! mu = 4e-7 pi
      real*8, PARAMETER:: muo2  = mu/two		! mu/2
      real*8, PARAMETER:: muo4  = mu/four		! mu/4
      real*8, PARAMETER:: muopi = 4.0d-7		! mu/pi
      real*8, PARAMETER:: muo2pi= 2.0d-7		! mu/(2*pi)
      real*8, PARAMETER:: hilow = 1.0d-3

      real*8  asq, Rsq, R, zsq, rhosq, aR
      real*8  f1, sqrtf1, f2, f3, f4, ksq
      integer err
      real*8  Kint, Eint, S, Dk4
      real*8  Bfac, Bror
      real*8  Roasq, fd, sqrtfd, Bfak, Br, Brn

      SAVE

c-----------------------------------------------------------------------
c Use hilow = 1.d-03 (approximately) with 64 bit arithmetic 
c as the divide between the two methods (higher ksq and lower ksq) 
c of calculating B.
c-----------------------------------------------------------------------


      if (cur.eq.0) then		! no need to compute
       B(1) = zero
       B(2) = zero
       B(3) = zero
       RETURN
      end if

c Calculate argument for elliptic integrals

      Rsq    = p(1)**2 + p(2)**2	! R^2 = x^2 + y^2
      zsq    = p(3)**2			! z^2
      R      = sqrt(Rsq)		! R
      aR     = a*R			! aR
      f1     = (a + R)**2 + zsq
      ksq    = four*aR/f1		! k^2
C      WRITE(0,*) ''
C      WRITE(0,*) 'ksq =', ksq		! for DEBUG
      if (ksq .ge. one) then		! error condition
       write(0,*) ''
       write(0,*) 'Attempt to calculate a magnetic field at the'
       write(0,*) 'exact location of a circular current element'
       write(0,*) 'in subroutine CIRCLEB.'
      end if

c Calculate elliptic integrals of 1st and 2nd kind

      CALL ELLIPINTS(ksq, Kint, Eint, S, Dk4, err)

      if (err.ne.0) then
       WRITE(0,*) 'Subroutine ELLIPINTS returned err', err, 'to CIRCLEB'
      end if

c Calculate magnetic field

      asq    = a**2
      rhosq  = Rsq + zsq		! R^2 + z^2
      sqrtf1 = sqrt(f1)
      f2     = asq + rhosq		! a^2 + R^2 + z^2
      f3     = (a - R)**2 + zsq
      f4     = asq - rhosq		! a^2 - R^2 - z^2

      Bfac = cur*Kint 
      IF (ksq .gt. hilow) THEN		! calculate B using S

        if (R.eq.zero .or. zsq.eq.zero) then	! point on an axis
         B(1) = zero				! Bx
         B(2) = zero				! By

        else					! point is general
         Bror = two*muopi*Bfac*asq*p(3)/(f1*sqrtf1*f3)
     &                   *(one - 2*(f2/f1)*(S/(ksq*ksq)))
         B(1) = Bror*p(1)			! Bx
         B(2) = Bror*p(2)			! By
        end if  	! end if(R.eq.zero .or. z.eq.zero) else
        B(3)  = muo2pi*Bfac/(sqrtf1*f3)
     &                   *(two*(asq - aR) - half*f4*(ksq + S))
C       WRITE(0,*) 'Numerical with S:'
C       WRITE(0,903) B(1), B(2), B(3)

      ELSE			! ksq is small, calculate using Dk4
        if (R.eq.zero .or. zsq.eq.zero) then	! point on an axis
         B(1) = zero				! Bx
         B(2) = zero				! By

        else					! general point
         Bror = two*muopi*Bfac*asq*p(3)/(f1*sqrtf1*f3)
     &                   *(one - two*(f2/f1)*Dk4/Kint)
         B(1) = Bror*p(1)			! Bx
         B(2) = Bror*p(2)			! By
        end if  	! end if(R.eq.zero .or. z.eq.zero) else
        B(3)  = muo2pi*bfac/(sqrtf1*f3)
     &            *(two*(asq - ar) - half*f4*(ksq + ksq*ksq*Dk4/Kint))
C       WRITE(0,*) 'Numerical with Dk4:'
C       WRITE(0,903) B(1), B(2), B(3)

      END IF		! end if (ksq .gt. hilow) else


  901 format(1x, 5f15.9)	! some formats useful for 
  902 format(1x, 3f20.15)	! test and debug
  903 format(1x, 1p3e22.14)

      RETURN
      END


c------------------    END OF SUBROUTINE CIRCLEB    --------------------

CC The following are useful to check the calculations of B in the 
CC active part of this subroutine in the ksq ranges where they overlap
CC with high accuracy.
C
C
CC This is the "Standard Numerical" calculation of B that was in the 
CC original TRIP3D. It is double precision accurate everywhere for all  
CC three Bcomponents only when 0.3 .lt. ksq .lt. 1 (approximately).
C
C      Br    = muo2pi*cur*p(3)/(R*sqrtf1)*(f2/f3*Eint - Kint)
C      B(1)  = Br*p(1)/R			! Bx
C      B(2)  = Br*p(2)/R			! By
C      B(3)   = muo2pi*cur/sqrtf1*(f4/f3*Eint + Kint)
CC       WRITE(0,*) 'Standard numerical:'
C      WRITE(0,903) B(1), B(2), B(3)	! for DEBUG
C
C
C
CC This is the calculation of B from a series valid for "Small R/a"
CC developed by M.J. Schaffer specifically to test this code in 
CC that limit. Div B and curl B are zero to the order of the expansion.
CC (M. Schaffer's notes dated 2006 Mar 4)
C
C      Roasq  = Rsq/asq			! (a/R)^2
C      fd     = one + zsq/asq
C      sqrtfd = sqrt(fd)
C      Bfak = muo2*cur/(a*fd*sqrtfd)
C      Brn  = 1.5d0*(one - 2.5d0*(Roasq/fd)*(one - 1.75d0/fd) )
C      Brn  = Bfak*Brn/fd
C      B(1) = Brn*p(1)*p(3)/asq 		! Bx
C      B(2) = Brn*p(2)*p(3)/asq 		! By
C      B(3) = Bfak*(one - 3.0d0*(Roasq/fd)*(one - 1.25d0/fd))
C       WRITE(0,*) 'Small R/a series:'
C      WRITE(0,903) B(1), B(2), B(3)	! for DEBUG
C
C
C      RETURN
C      END

c------------------    END OF SUBROUTINE CIRCLEB    --------------------
