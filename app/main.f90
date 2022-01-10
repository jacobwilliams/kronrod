    program main

    use mpmodule
    use iso_fortran_env, only: wp => real64

    implicit none

    integer, parameter :: ndp1 = 500
    integer, parameter :: ndp2 = 1000
    integer, parameter :: neps1 = -ndp1
    integer, parameter :: neps2 = -ndp2
    integer, parameter :: nq1 = 11
    integer, parameter :: nq2 = 12 * 2 ** nq1
    integer, parameter :: nwds1 = int (ndp1 / mpdpw + 2)
    integer, parameter :: nwds2 = int (ndp2 / mpdpw + 2)

    type (mp_real) :: zero,one_half,one,two,three,four,pi,half_pi,eight
    type (mp_real) :: epsilon    ! epsilon ( x )

    integer :: i !! counter
    integer :: icase
    integer :: n

    ! kronrod variables:
    !integer,dimension(*),parameter :: cases = [7, 10, 87]
    integer,dimension(*),parameter :: cases = [(i, i = 7, 87)]
    type (mp_real),dimension(:),allocatable :: w1
    type (mp_real),dimension(:),allocatable :: w2
    type (mp_real),dimension(:),allocatable :: x
    type (mp_real) :: eps

    if (ndp2 > mpipl) then
        write (6, '("Increase default precision in module MPFUNF.")')
        stop
      endif

    call init()

    do icase = 1, size(cases)

        write(*,*) ''
        write(*,*) '----------------------'

        n = cases(icase)

        if (allocated(w1)) deallocate(w1)
        if (allocated(w2)) deallocate(w2)
        if (allocated(x)) deallocate(x)

        allocate(w1(n+1))
        allocate(w2(n+1))
        allocate(x(n+1))

        call kronrod ( n, eps, x, w1, w2 )

        write(*,*) ''
        write(*,*) 'n = ', n
        write(*,*) ''
        write(*,*) 'abscissas:'
        do i = 1, n+1
            call mpwrite (6, 60, 40, x(i))
        end do
        write(*,*) ''
        write(*,*) 'Gauss-Kronrod weights:'
        do i = 1, n+1
            call mpwrite (6, 60, 40, W1(i))
        end do
        write(*,*) ''
        write(*,*) 'Gauss weights:'
        do i = 1, n+1
            call mpwrite (6, 60, 40, W2(i))
        end do

    end do

    contains
!*********************************************************************************

    subroutine init()

        !! initialize constants

        implicit none

        zero     = mpreal(0.0_wp, nwds2)
        one_half = mpreal(0.5_wp, nwds2)
        one      = mpreal(1.0_wp, nwds2)
        two      = mpreal(2.0_wp, nwds2)
        three    = mpreal(3.0_wp, nwds2)
        four     = mpreal(4.0_wp, nwds2)
        eight    = mpreal(8.0_wp, nwds2)

        pi = mppi (nwds2)
        half_pi = pi / two

        epsilon  = mpreal (10.0_wp, nwds2) ** (-100)
        eps      = mpreal (10.0_wp, nwds2) ** (-100)

    end subroutine init

    subroutine kronrod ( n, eps, x, w1, w2 )

        !*****************************************************************************80
        !
        !! KRONROD adds N+1 points to an N-point Gaussian rule.
        !
        !  Discussion:
        !
        !    This subroutine calculates the abscissas and weights of the 2N+1
        !    point Gauss Kronrod quadrature formula which is obtained from the
        !    N point Gauss quadrature formula by the optimal addition of N+1 points.
        !
        !    The optimally added points are called Kronrod abscissas.  The
        !    abscissas and weights for both the Gauss and Gauss Kronrod rules
        !    are calculated for integration over the interval [-1,+1].
        !
        !    Since the quadrature formula is symmetric with respect to the origin,
        !    only the nonnegative abscissas are calculated.
        !
        !    Note that the code published in Mathematics of Computation
        !    omitted the definition of the variable which is here called COEF2.
        !
        !  Storage:
        !
        !    Given N, let M = ( N + 1 ) / 2.
        !
        !    The Gauss-Kronrod rule will include 2*N+1 points.  However, by symmetry,
        !    only N + 1 of them need to be listed.
        !
        !    The arrays X, W1 and W2 contain the nonnegative abscissas in decreasing
        !    order, and the weights of each abscissa in the Gauss-Kronrod and
        !    Gauss rules respectively.  This means that about half the entries
        !    in W2 are zero.
        !
        !    For instance, if N = 3, the output is:
        !
        !    I      X               W1              W2
        !
        !    1    0.960491        0.104656         0.000000
        !    2    0.774597        0.268488         0.555556
        !    3    0.434244        0.401397         0.000000
        !    4    0.000000        0.450917         0.888889
        !
        !    and if N = 4, (notice that 0 is now a Kronrod abscissa)
        !    the output is
        !
        !    I      X               W1              W2
        !
        !    1    0.976560        0.062977        0.000000
        !    2    0.861136        0.170054        0.347855
        !    3    0.640286        0.266798        0.000000
        !    4    0.339981        0.326949        0.652145
        !    5    0.000000        0.346443        0.000000
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license.
        !
        !  Modified:
        !
        !    03 August 2010
        !
        !  Author:
        !
        !    Original FORTRAN77 version by Robert Piessens, Maria Branders.
        !    FORTRAN90 version by John Burkardt.
        !
        !  Reference:
        !
        !    Robert Piessens, Maria Branders,
        !    A Note on the Optimal Addition of Abscissas to Quadrature Formulas
        !    of Gauss and Lobatto,
        !    Mathematics of Computation,
        !    Volume 28, Number 125, January 1974, pages 135-139.
        !
        !  Parameters:
        !
        !    Input, integer N, the order of the Gauss rule.
        !
        !    Input, type (mp_real) EPS, the requested absolute accuracy of the
        !    abscissas.
        !
        !    Output, type (mp_real) X(N+1), the abscissas.
        !
        !    Output, type (mp_real) W1(N+1), the weights for
        !    the Gauss-Kronrod rule.
        !
        !    Output, type (mp_real) W2(N+1), the weights for
        !    the Gauss rule.
        !
          implicit none

          integer n
          type (mp_real) eps
          type (mp_real) w1(n+1)
          type (mp_real) w2(n+1)
          type (mp_real) x(n+1)

          type (mp_real) ak
          type (mp_real) an
          type (mp_real) b(((n+1)/2)+1)
          type (mp_real) bb
          type (mp_real) c
          type (mp_real) coef
          type (mp_real) coef2
          type (mp_real) d
          logical even
          integer i
          integer k
          integer l
          integer ll
          integer m
          type (mp_real) s
          type (mp_real) tau(( n + 1 ) / 2 )
          type (mp_real) x1
          type (mp_real) xx
          type (mp_real) y

          m = ( n + 1 ) / 2
          even = ( 2 * m == n )

          d = two
          an = zero
          do k = 1, n
            an = an + one
            d = d * an / ( an + one_half )
          end do
        !
        !  Calculation of the Chebyshev coefficients of the orthogonal polynomial.
        !
          tau(1) = ( an + two ) / ( an + an + three )
          b(m) = tau(1) - one
          ak = an

          do l = 1, m - 1

            ak = ak + two
            tau(l+1) = ( ( ak - one ) * ak &
              - an * ( an + one ) ) * ( ak + two ) * tau(l) &
              / ( ak * ( ( ak + three ) * ( ak + two ) &
              - an * ( an + one ) ) )
            b(m-l) = tau(l+1)

            do ll = 1, l
              b(m-l) = b(m-l) + tau(ll) * b(m-l+ll)
            end do

          end do

          b(m+1) = one
        !
        !  Calculation of approximate values for the abscissas.
        !
          bb = sin ( half_pi / ( an + an + one ) )
          x1 = sqrt ( one - bb * bb )
          s = two * bb * x1
          c = sqrt ( one - s * s )
          coef = one - ( one - one / an ) / ( eight * an * an )
          xx = coef * x1
        !
        !  Coefficient needed for weights.
        !
        !  COEF2 = 2^(2*n+1) * n! * n! / (2n+1)!
        !

          coef2 = two / real ( 2 * n + 1, kind = wp )
          do i = 1, n
            coef2 = coef2 * four * real ( i, kind = wp ) / real ( n + i, kind = wp )
          end do
        !
        !  Calculation of the K-th abscissa (a Kronrod abscissa) and the
        !  corresponding weight.
        !
          do k = 1, n, 2

            call abwe1 ( n, m, eps, coef2, even, b, xx, w1(k) )
            w2(k) = zero

            x(k) = xx
            y = x1
            x1 = y * c - bb * s
            bb = y * s + bb * c

            if ( k == n ) then
              xx = zero
            else
              xx = coef * x1
            end if
        !
        !  Calculation of the K+1 abscissa (a Gaussian abscissa) and the
        !  corresponding weights.
        !
            call abwe2 ( n, m, eps, coef2, even, b, xx, w1(k+1), w2(k+1) )

            x(k+1) = xx
            y = x1
            x1 = y * c - bb * s
            bb = y * s + bb * c
            xx = coef * x1

          end do
        !
        !  If N is even, we have one more Kronrod abscissa to compute,
        !  namely the origin.
        !
          if ( even ) then
            xx = zero
            call abwe1 ( n, m, eps, coef2, even, b, xx, w1(n+1) )
            w2(n+1) = zero
            x(n+1) = xx
          end if

        end

        subroutine abwe1 ( n, m, eps, coef2, even, b, x, w )

        !*****************************************************************************80
        !
        !! ABWE1 calculates a Kronrod abscissa and weight.
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license.
        !
        !  Modified:
        !
        !    03 August 2010
        !
        !  Author:
        !
        !    Original FORTRAN77 version by Robert Piessens, Maria Branders.
        !    FORTRAN90 version by John Burkardt.
        !
        !  Reference:
        !
        !    Robert Piessens, Maria Branders,
        !    A Note on the Optimal Addition of Abscissas to Quadrature Formulas
        !    of Gauss and Lobatto,
        !    Mathematics of Computation,
        !    Volume 28, Number 125, January 1974, pages 135-139.
        !
        !  Parameters:
        !
        !    Input, integer N, the order of the Gauss rule.
        !
        !    Input, integer M, the value of ( N + 1 ) / 2.
        !
        !    Input, type (mp_real) EPS, the requested absolute accuracy of the
        !    abscissas.
        !
        !    Input, type (mp_real) COEF2, a value needed to compute weights.
        !
        !    Input, logical EVEN, is TRUE if N is even.
        !
        !    Input, type (mp_real) B(M+1), the Chebyshev coefficients.
        !
        !    Input/output, type (mp_real) X; on input, an estimate for
        !    the abscissa, and on output, the computed abscissa.
        !
        !    Output, type (mp_real) W, the weight.
        !
          implicit none

          integer m

          type (mp_real) ai
          type (mp_real) b(m+1)
          type (mp_real) b0
          type (mp_real) b1
          type (mp_real) b2
          type (mp_real) coef2
          type (mp_real) d0
          type (mp_real) d1
          type (mp_real) d2
          type (mp_real) delta
          type (mp_real) dif
          type (mp_real) eps
          logical even
          type (mp_real) f
          type (mp_real) fd
          integer i
          integer iter
          integer k
          integer ka
          integer n
          type (mp_real) w
          type (mp_real) x
          type (mp_real) yy

          if ( x == zero ) then
            ka = 1
          else
            ka = 0
          end if
        !
        !  Iterative process for the computation of a Kronrod abscissa.
        !
          do iter = 1, 50

            b1 = zero
            b2 = b(m+1)
            yy = four * x * x - two
            d1 = zero

            if ( even ) then
              ai = m + m + 1
              d2 = ai * b(m+1)
              dif = two
            else
              ai = m + 1
              d2 = zero
              dif = one
            end if

            do k = 1, m
              ai = ai - dif
              i = m - k + 1
              b0 = b1
              b1 = b2
              d0 = d1
              d1 = d2
              b2 = yy * b1 - b0 + b(i)
              if ( .not. even ) then
                i = i + 1
              end if
              d2 = yy * d1 - d0 + ai * b(i)
            end do

            if ( even ) then
              f = x * ( b2 - b1 )
              fd = d2 + d1
            else
              f = one_half * ( b2 - b0 )
              fd = four * x * d2
            end if
        !
        !  Newton correction.
        !
            delta = f / fd
            x = x - delta

            if ( ka == 1 ) then
              exit
            end if

            if ( abs ( delta ) <= eps ) then
              ka = 1
            end if

          end do
        !
        !  Catch non-convergence.
        !
          if ( ka /= 1 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'ABWE1 - Fatal error!'
            write ( *, '(a)' ) '  Iteration limit reached.'
            write ( *, '(a,g14.6)' ) '  EPS is ', eps
            write ( *, '(a,g14.6)' ) '  Last DELTA was ', delta
            stop 1
          end if
        !
        !  Computation of the weight.
        !
          d0 = one
          d1 = x
          ai = zero
          do k = 2, n
            ai = ai + one
            d2 = ( ( ai + ai + one ) * x * d1 - ai * d0 ) / ( ai + one )
            d0 = d1
            d1 = d2
          end do

          w = coef2 / ( fd * d2 )

        end

        subroutine abwe2 ( n, m, eps, coef2, even, b, x, w1, w2 )

        !*****************************************************************************80
        !
        !! ABWE2 calculates a Gaussian abscissa and 2 weights.
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license.
        !
        !  Modified:
        !
        !    30 April 2013
        !
        !  Author:
        !
        !    Original FORTRAN77 version by Robert Piessens, Maria Branders.
        !    FORTRAN90 version by John Burkardt.
        !
        !  Reference:
        !
        !    Robert Piessens, Maria Branders,
        !    A Note on the Optimal Addition of Abscissas to Quadrature Formulas
        !    of Gauss and Lobatto,
        !    Mathematics of Computation,
        !    Volume 28, Number 125, January 1974, pages 135-139.
        !
        !  Parameters:
        !
        !    Input, integer N, the order of the Gauss rule.
        !
        !    Input, integer M, the value of ( N + 1 ) / 2.
        !
        !    Input, type (mp_real) EPS, the requested absolute accuracy of the
        !    abscissas.
        !
        !    Input, type (mp_real) COEF2, a value needed to compute weights.
        !
        !    Input, logical EVEN, is TRUE if N is even.
        !
        !    Input, type (mp_real) B(M+1), the Chebyshev coefficients.
        !
        !    Input/output, type (mp_real) X; on input, an estimate for
        !    the abscissa, and on output, the computed abscissa.
        !
        !    Output, type (mp_real) W1, the Gauss-Kronrod weight.
        !
        !    Output, type (mp_real) W2, the Gauss weight.
        !
          implicit none

          integer m

          type (mp_real) ai
          type (mp_real) an
          type (mp_real) b(m+1)
          type (mp_real) coef2
          type (mp_real) delta
          type (mp_real) eps
          logical even
          integer i
          integer iter
          integer k
          integer ka
          integer n
          type (mp_real) p0
          type (mp_real) p1
          type (mp_real) p2
          type (mp_real) pd0
          type (mp_real) pd1
          type (mp_real) pd2
          type (mp_real) w1
          type (mp_real) w2
          type (mp_real) x
          type (mp_real) yy

          if ( x == zero ) then
            ka = 1
          else
            ka = 0
          end if
        !
        !  Iterative process for the computation of a Gaussian abscissa.
        !
          do iter = 1, 50

            p0 = one
            p1 = x
            pd0 = zero
            pd1 = one
        !
        !  When N is 1, we need to initialize P2 and PD2 to avoid problems with DELTA.
        !
            if ( n <= 1 ) then
              if ( epsilon < abs ( x ) ) then
                p2 = ( three * x * x - one ) / two
                pd2 = three * x
              else
                p2 = three * x
                pd2 = three
              end if
            end if

            ai = zero
            do k = 2, n
              ai = ai + one
              p2 = ( ( ai + ai + one ) * x * p1 - ai * p0 ) / ( ai + one )
              pd2 = ( ( ai + ai + one ) * ( p1 + x * pd1 ) - ai * pd0 ) &
                / ( ai + one )
              p0 = p1
              p1 = p2
              pd0 = pd1
              pd1 = pd2
            end do
        !
        !  Newton correction.
        !
            delta = p2 / pd2
            x = x - delta

            if ( ka == 1 ) then
              exit
            end if

            if ( abs ( delta ) <= eps ) then
              ka = 1
            end if

          end do
        !
        !  Catch non-convergence.
        !
          if ( ka /= 1 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'ABWE2 - Fatal error!'
            write ( *, '(a)' ) '  Iteration limit reached.'
            write ( *, '(a,g14.6)' ) '  EPS is ', eps
            write ( *, '(a,g14.6)' ) '  Last DELTA was ', delta
            stop 1
          end if
        !
        !  Computation of the weight.
        !
          an = n

          w2 = two / ( an * pd2 * p0 )

          p1 = zero
          p2 = b(m+1)
          yy = four * x * x - two
          do k = 1, m
            i = m - k + 1
            p0 = p1
            p1 = p2
            p2 = yy * p1 - p0 + b(i)
          end do

          if ( even ) then
            w1 = w2 + coef2 / ( pd2 * x * ( p2 - p1 ) )
          else
            w1 = w2 + two * coef2 / ( pd2 * ( p2 - p0 ) )
          end if

          return
        end

    end program main