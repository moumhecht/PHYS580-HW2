!-----------------------------------------------------------------------
!Module: quadrature
!-----------------------------------------------------------------------
!! By: Malida Hecht
!!
!! This module contains subroutines and functions that computes Boole's rule
!! and Monte Carlo integration.
!!
!! Boole's quadrature is a numerical solution for integration over a 5 point
!! interval that approximates the function to a 4th degree polynomial.
!!
!! Monte Carlo Quadrature is a numerical solution for integration that
!! calls on random points to evaluate the integrand.
!!
!! booles_quadrature integrates the function using booles_rule and only passes
!!          5 points in the array at a time
!!
!! booles_rule calculates Boole's five point rule
!!
!!
!! monte_carlo_quad calculates the integral on the interval [a,b] using 
!!          monte carlo quadrature by computing the volume and multiplying
!!          by the probability that the points are successfully in the 
!!          area under the curve and computes the error from monte carlo
!!          integration
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! 
!! monte_carlo_quad
!!----------------------------------------------------------------------
!! Included functions:
!!
!! booles_quadrature
!! booles_rule
!-----------------------------------------------------------------------
module quadrature
use types

implicit none

private
public :: booles_quadrature, monte_carlo_quad

!-----------------------------------------------------------------------
!Interface: func
!-----------------------------------------------------------------------
!! This defines a new type of procedure in order to allow callbacks
!! in the Monte Carlo quadrature subroutine of an arbitrary function that is given
!! as input and declared as a procedure
!!
!! The arbitrary function receives two rank 1 arrays of arbitrary size.
!! The first array contains an n-dimensional vector representing the
!! point sampled by the Monte Carlo method. The second is a "work array"
!! that contains parameters  necessary to calculate the function to be
!! integrated.
!!----------------------------------------------------------------------
interface
    real(dp) function func(x, data)
        use types, only : dp
        implicit none
        real(dp), intent(in) :: x(:), data(:)
        ! This is the interface we saw in class that allows callbacks
    end function func
end interface

contains

!-----------------------------------------------------------------------
!! Function: booles_quadrature
!-----------------------------------------------------------------------
!! By: Malida Hecht
!!
!! This function evaluates the integral of the function fx with spacing
!! delta_x using booles_rule by passing slices of array fx.
!! ----------------------------------------------------------------------
!! Input:
!!
!! fx           real        Array containing the evaluated function
!! delta_x      real        Distance between the evaluation points
!!----------------------------------------------------------------------
!! Output:
!!
!! s            real        Result of the Boole's quadrature
!-----------------------------------------------------------------------
real(dp) function booles_quadrature(fx, delta_x) result(s)
    implicit none
    real(dp), intent(in) :: fx(1:), delta_x

    integer :: fx_size, i

    fx_size = size(fx)

    ! As the diagram below shows, only certain number of grid points
    ! fit the scheme of Boole's quadrature. Implement a test 
    ! to make sure that the number of evaluated points in the fx array
    ! is the correct one

    ! |--interval 1---|--interval 2---|--interval 3---|
    ! 1   2   3   4   5   6   7   8   9   10  11  12  13
    ! |---|---|---|---|---|---|---|---|---|---|---|---|
    ! x0  x1  x2  x3  x4
    !                 x0  x1  x2  x3  x4
    !                                 x0  x1  x2  x3  x4

    !  5 mod 4 = 1
    !  9 mod 4 = 1
    ! 13 mod 4 = 1


     if  (MOD(fx_size-1,4) .NE. 0) then
         print *, 'fx array size in booles_quadrature has to be multiple of 4'
         stop
     endif

    ! We could implement the full integration here, however to make a cleaner,
    ! easy to read (and debug or maintain) code we will define a smaller
    ! function that returns Boole's five point rule and pass slices (1:5), (5:9),
    ! (9:13), ... of fx to such function to then add all the results. 

    s = 0._dp

     do i = 1, (fx_size-1), 4
         s = s + booles_rule(fx(i:i+4),delta_x)
         
     enddo
end function booles_quadrature

!-----------------------------------------------------------------------
!! Function: booles_rule
!-----------------------------------------------------------------------
!! By: Malida Hecht
!!
!! This function computes Boole's five point rule for only 5 points,
!! Boole's five point rule is given by:
!! dx * (2/45) * (7*f(x4)+32*f(x3)+12f(x2)+32f(x1)+7f(x0))
!! ----------------------------------------------------------------------
!! Input:
!!
!! fx           real        Array containing the evaluated function
!! delta_x      real        Distance between the evaluation points
!!----------------------------------------------------------------------
!! Output:
!!
!! s            real        Result of the Boole's quadrature
!-----------------------------------------------------------------------
real(dp) function booles_rule(fx, delta_x) result(s)
    implicit none
    real(dp), intent(in) :: fx(1:), delta_x

    integer :: fx_size
    real(dp) :: fx0, fx1, fx2, fx3, fx4

    fx_size = size(fx)

    ! Let's make an additional test to make sure that the array
    ! received has 5 and only 5 points 

     if (fx_size .NE. 5) then
        print * , 'The array does not have exaclty 5 points'
        stop
     endif
    
    ! define f_0, f_1, f_2, .. f_4 from array fx
     fx0 = fx(1)
     fx1 = fx(2)
    ! print*, fx1
     fx2 = fx(3)
     fx3 = fx(4)
     fx4 = fx(5)
    
    ! compute sum from Boole's quadrature
     s = delta_x * (2._dp / 45._dp) * (7._dp * fx4 + 32._dp * fx3 + 12._dp * fx2 + 32._dp * fx1 + 7._dp*fx0 )
end function booles_rule

!-----------------------------------------------------------------------
!! Subroutine: monte_carlo_quad
!-----------------------------------------------------------------------
!! By: Malida Hecht
!!
!! This subroutine computes the integral of a function on the interval [a,b]
!! for a generalized vector szie.
!! The volume of the area of integration is computed using a do loop
!! 
!! Random numbers on the interval [0,1] are generated in an array x_vector
!! and is rescaled to the integration volume
!! The function, and its square are evaluated at that point
!!
!! Monte carlo integration is then computed by V/M * SUM(fx)
!!
!! The monte carlo error is given by  V/SQRT(M) * sigma_f
!! Where sigma_f^^2 is given by  SUM(fx^2) / M - (SUM(fx)/M)^2
!! ----------------------------------------------------------------------
!! Input:
!!
!! f            procedure   function to be integrated
!! a            real        array containing the lower limits of the integral
!! b            real        array containing the upper limits of the integral
!! data         real        array containing parameters necessary to calculate the function f
!! n_samples    integer     number of sample points in the Monte Carlo integration
!!----------------------------------------------------------------------
!! Output:
!!
!! s            real        Result of the Monte Carlo integral
!! sigma_s      real        Estimate of the uncertainty in the Monte Carlo integral
!-----------------------------------------------------------------------
subroutine monte_carlo_quad(f, a, b, data, n_samples, s, sigma_s)
    implicit none
    procedure(func) :: f
    real(dp), intent(in) :: a(:), b(:), data(:)
    integer, intent(in) :: n_samples
    real(dp), intent(out) :: s, sigma_s
    real(dp) :: volume, f_x, f2_x, s_squared

    integer :: i, vector_size
    real(dp), allocatable :: x_vector(:), fx(:), fx_square(:)! ...you might need to declare other arrays here


    vector_size = size(a)

    ! We're defining a Monte Carlo routine that works for an arbitrary number of 
    ! dimensions in the integral (Remember, that's the advantage of Monte Carlo integration,
    ! it's very efficient for high dimensional integrals)

    ! Since a and b give the lower and upper limits they need to have the same size.
    ! Make a check to see if they do have the same size

     if (size(a) .NE. size(b) )then
          print *, 'a and b arrays in monte_carlo_quad have to be the same size'
         stop       
     endif
     !print*, a, b

    ! Here we allocate memory for the vector containing the sample points and 
    ! for a vector that contains the evaluated function
    allocate(x_vector(1:vector_size))
    allocate(fx(1:n_samples))
    allocate(fx_square(1:n_samples))

    ! Initializing Volume

    volume = 1._dp

    ! Computing volume
     ! V = (b1 -a1) * (b2-a2) *(b3-a3) = x* y * z 
    do i = 1, vector_size
    volume =  (b(i) - a(i)) * volume  ! V = (b1 -a1) * (b2-a2)
    enddo 

    do i= 1, n_samples
        ! NEVER USE THE INTRINSIC FUNCTION RAND() to generate 'random' numbers
        ! It's NOT reliable when a large number of random numbers are necessary.
        ! The only reason it is still in some compilers is for backwards compatibility with old code. 

        ! Instead use the intrinsic function random_number() it's a more modern version in the
        ! current standard with a period of 2^{1024} - 1 (that's a HUGE number)

        ! However, caution should be used in parallel applications. Other pseudo random number
        ! generators are more appropriate for such cases. But that's beyond the scope of this class

        call random_number(x_vector) !generates an array with random numbers in the [0,1) interval
        x_vector = a + x_vector*(b-a) !rescaling to the integration volume [a,b)
        fx(i) = f(x_vector,data)
        fx_square(i) = f(x_vector,data) * f(x_vector,data)

    enddo

    f_x  = sum(fx)
    f2_x = sum(fx_square)
    ! I'll leave the calculation of the integral and it's uncertainty to you
    s = (volume / n_samples) * f_x ! integral
    s_squared = (f2_x / n_samples ) - (f_x / n_samples)**2._dp !sigma f squared
    sigma_s = volume * sqrt(s_squared/ n_samples)
end subroutine monte_carlo_quad

end module quadrature