!-----------------------------------------------------------------------
!Module: read_write
!-----------------------------------------------------------------------
!! By: Malida Hecht
!!
!! The 'read' subroutines in this module takes user input and stores the 
!! values for the dimensions of the box and sphere for the nuclear reactor,
!! the coordinate positions for the detector location, step size for integration,
!! sample points for integration, and minimum and maximum radius sizes for a nuclear
!! reactor with a box that is spherically hollow.
!!
!! The 'write' subroutines write the results of Boole's integration and 
!! Monte Carlo integration to data files.
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! read_input
!! read_advanced_input
!! write_neutron_flux
!! write_advanced_flux
!!----------------------------------------------------------------------
!! Included functions:
!!
!! read_real
!! read_integer
!-----------------------------------------------------------------------
module read_write
use types
use neutron_flux, only : box_flux_booles, large_x0_flux, box_flux_monte_carlo, total_flux_booles, hollow_box_flux_mc

implicit none

private
public :: read_input, write_neutron_flux, read_advanced_input, write_advanced_flux

contains

!-----------------------------------------------------------------------
!! Subroutine: read_input
!-----------------------------------------------------------------------
!! By: Malida Hecht
!!
!! This subroutine calls the functions read_real and read_integer that takes
!! and stores user input  to variables.
!!----------------------------------------------------------------------
!! Output:
!!
!! depth        real        Depth of the rectangular nuclear reactor
!! width        real        Width of the rectangular nuclear reactor
!! height       real        Height of the rectangular nuclear reactor
!! y_zero       real        y coordinate of the detector's position
!! x_min        real        minimum x coordinate of the detector's position
!! x_max        real        maximum x coordinate of the detector's position
!! x_step       real        increment size for the x coordinate of the detector's position
!! n_grid       integer     number of grid points in each dimension of Boole's integration
!! m_samples    integer     number of sample points in the Monte Carlo integration
!-----------------------------------------------------------------------
subroutine read_input(depth, width, height, y_zero, x_min, x_max, x_step, n_grid, m_samples)
    implicit none
    real(dp), intent(out) :: depth, width, height, y_zero, x_min, x_max, x_step
    integer, intent(out) ::  n_grid, m_samples

    ! use the print statement to give a message to the user describing
    ! what the program does and what input should the user give to the
    ! program  
    print *, 'This program calculates neutron flux of a reactor using'
    print *, "Boole's Quadrature and Monte Carlo Quadrature"
    print *, "For varying distances x_0 of a detector location away from the reactor."
    print *, "Please provide the following inputs:"
    print *, ''

    ! In assignment 01 there was a single input and we used a do loop 
    ! to ask for such input. The loop was only exited once the input given
    ! was the correct one (the user had to give an actual number)
    ! and we put that loop directly in this function.

    ! However, here we need several inputs from the user, it wouldn't 
    ! make sense to write a do loop for each input needed. Instead 
    ! it's more efficient to put such loop inside a function and have
    ! such function return the input given by the user. 

    ! In order to ask for a different input in the screen, the new function
    ! we'll define below will take a string as input with the name of the
    ! variable we want the user to give us.
    depth = read_real('depth D')
    width = read_real('width W')
    height = read_real('height H')
    y_zero = read_real("y coordinate of the detector's position y_0")
    x_min = read_real("minimum x coordinate of detector's position x_min")
    x_max = read_real("maximum x coordinate of detector's position x_max")
    x_step = read_real("increment size for x coordinate of detector's position x_step")

    ! ...
    ! make sure to get a value for every real number the user needs to give
    ! The function read_real used above returns a double precision real,
    ! however n_grid and m_samples are integers, that means that we need
    ! another function to get integers. For that we'll define read_integer
    ! as well
    n_grid = read_integer('number of lattice points N')
    m_samples = read_integer('number Monte Carlo samples M')

    print *, ''
end subroutine read_input

!-----------------------------------------------------------------------
!! Function: read_real
!-----------------------------------------------------------------------
!! By: Malida Hecht
!!
!! This function takes a string, and prompts the user to input a value for 
!! each variable stored.
!! The input is then passed to check if it is a real number, and then is checked
!! to make sure that it is nonzero and positive. 
!! If the input is not a positive number, the user will be prompted to provide another
!! input.
!!----------------------------------------------------------------------
!! Input:
!!
!! name     character   A string with a brief description of the value being asked for
!!----------------------------------------------------------------------
!! Output:
!!
!! x        real        A positive non negative number given by the user
!-----------------------------------------------------------------------
real(dp) function read_real(name) result(x)
    implicit none
    character(len=*), intent(in) :: name
    character(len=120) :: string
    integer :: ierror

    print *,  ''
    print *, 'Provide a nonzero positive value for the '//trim(name)//':'

    ! Use the do loop from assignment 01 to get an input from the user.
    ! Also, modify accordingly to make sure that the input is a non zero 
    ! positive number (there's no box with height H = -4.5 or width W = 0)

    ! Remember to declare above the type for `string` and `ierror`



     do
        read(*,'(a)', iostat=ierror) string
        if(string /= '') then
            read(string, *, iostat=ierror) x
            if (ierror == 0) then
                if (x .GT. 0) exit !checks for positive number, will send prompt to provide a new input
                print *,  "'"//trim(string)//"'"//' is not a positive number, please provide a positive number'
            else 
                print *, "'"//trim(string)//"'"//' is not a number, please provide a number'
            endif
        else
            print *, 'that was an empty input, please provide a number'
        endif
    enddo

    

end function read_real

!-----------------------------------------------------------------------
!! Function: read_integer
!-----------------------------------------------------------------------
!! By: Malida Hecht
!!
!! This function prompts the user to provide a positive integer. The function 
!! takes the input and verifies that it is nonzero and positive and returns the
!! value to a stored input.
!!----------------------------------------------------------------------
!! Input:
!!
!! name     character   A string with a brief description of the value being asked for
!!----------------------------------------------------------------------
!! Output:
!!
!! x        integer     A positive non negative number given by the user
!-----------------------------------------------------------------------
integer function read_integer(name) result(x)
    implicit none
    character(len=*), intent(in) :: name
    character(len=120) :: string
    integer :: ierror

    print *, ''
    print *, 'Provide a nonzero positive value for the '//trim(name)//':'

    ! Here you can use the same structure that you used in read_real.
    ! The fact that we declared x as an integer will take care of 
    ! checking that the number is an actual integer

    ! Also, put checks to make sure that the user gives positive non zero integers
    ! You can't integrate with -17 grid points! 

    ! Remember to declare above what type of variables string and ierror are
    
     do
        read(*,'(a)', iostat=ierror) string
        if(string /= '') then
            read(string, *, iostat=ierror) x
            if  (ierror == 0 ) then
                if (x .GT. 0 ) exit !checks for positive number, if it is not positive user will get prompt
                print *,  "'"//trim(string)//"'"//' is not a nonzero positive integer, please provide a nonzero positive integer'
            else 
                print *, "'"//trim(string)//"'"//' is not a number, please provide a number'
            endif
        else
            print *, 'that was an empty input, please provide a number'
        endif
    enddo
end function read_integer

!-----------------------------------------------------------------------
!Subroutine: write_neutron_flux
!-----------------------------------------------------------------------
!! By: Malida Hecht
!!
!! This subroutine takes inputs provided by the user and computes the neutron 
!! flux using boole's quadrature and monte carlo quadrature for increasing
!! values of x_0 (distance between detector and reactor) by calling on 
!! box_flux_booles and box_flux_monte_carlo. The subroutine also computes
!! the large x_0 approximation for flux, when the detector is further away 
!! from the reactor and treated as a point source. The values for each x_0
!! are then written onto the file 'results_basic.dat'
!!----------------------------------------------------------------------
!! Input:
!!
!! depth        real        Depth of the rectangular nuclear reactor
!! width        real        Width of the rectangular nuclear reactor
!! height       real        Height of the rectangular nuclear reactor
!! y_zero       real        y coordinate of the detector's position
!! x_min        real        minimum x coordinate of the detector's position
!! x_max        real        maximum x coordinate of the detector's position
!! x_step       real        increment size for the x coordinate of the detector's position
!! n_grid       integer     number of grid points in each dimension of Boole's integration
!! m_samples    integer     number of sample points in the Monte Carlo integration
!-----------------------------------------------------------------------
subroutine write_neutron_flux(depth, width, height, y_zero, x_min, x_max, x_step, n_grid, m_samples)
    implicit none
    real(dp), intent(in) :: depth, width, height, y_zero, x_min, x_max, x_step
    integer, intent(in) :: n_grid, m_samples

    real(dp) :: x_zero, box_booles, box_mc, box_large_x0, sigma_box
    character(len=*), parameter :: file_name = 'results_basic.dat'
    integer :: unit

    open(newunit=unit,file=file_name)
    write(unit,'(5a28)') 'x_0', 'booles', 'large x_0', 'monte carlo', 'MC uncertainty'
    x_zero = x_min
    do 
        if(x_zero > x_max) exit
        box_booles = box_flux_booles(depth, width, height, x_zero, y_zero, n_grid)
        box_large_x0 = large_x0_flux(depth, width, height, x_zero, y_zero)
        call box_flux_monte_carlo(depth, width, height, x_zero, y_zero, m_samples, box_mc, sigma_box)
        write(unit,'(5e28.16)') x_zero, box_booles, box_large_x0, box_mc, sigma_box
        x_zero = x_zero + x_step
    enddo
    close(unit)
    print *, 'The fluxes were written in the '//file_name//' file'
    !print *,
end subroutine write_neutron_flux

! Below are the read and write subroutines for the advanced part the project.
! Remember to make them public at the top of the module. And to also `use` them
! in the main program

!-----------------------------------------------------------------------
!! Subroutine: read_advanced_input
!-----------------------------------------------------------------------
!! By: Malida Hecht
!!
!! This subroutine takes the additional values needed by the program to 
!! compute the neutron flux of a box with a hollow sphere. 
!!----------------------------------------------------------------------
!! Output:
!!
!! x_zero       real        x coordinate of the detector's position
!! r_min        real        minimum radius of the hollow sphere
!! r_max        real        maximum radius of the hollow sphere
!! r_step       real        increment size for the radius of the hollow sphere
!-----------------------------------------------------------------------
 subroutine read_advanced_input(x_zero, r_min, r_max, r_step)
     implicit none
     real(dp), intent(out) :: x_zero, r_min, r_max, r_step

!     ! Base this new subroutine on the read_input one
    print *, "This part of the project uses Boole's and Monte Carlo Quadratures to compute"
    print *, "The neutron flux of a rectangular box with a spherical hollow inside"
    print *, "Please provide the following inputs:"

    x_zero = read_real("x coordinate of the detector's position")
    r_min = read_real("Minimum radius of the hollow sphere")
    r_max = read_real("Maximum radius of the hollow sphere")
    r_step = read_real("Increment size for the radius of the hollow sphere")
 end subroutine read_advanced_input

!-----------------------------------------------------------------------
!Subroutine: write_advanced_flux
!-----------------------------------------------------------------------
!! By: Malida Hecht
!!
!! This subroutine takes inputs provided by the user and computes the neutron 
!! flux  of a box with a hollow sphere using boole's quadrature and monte carlo 
!! quadrature for various values of radius by calling on 
!! total_flux_booles and hollow_box_flux_mc. The subroutine also computes
!! the flux for the box when it is full. The values for each radius
!! are then written onto the file 'results_advanced.dat'
!!----------------------------------------------------------------------
!! Input:
!! 
!! depth        real        Depth of the rectangular nuclear reactor
!! width        real        Width of the rectangular nuclear reactor
!! height       real        Height of the rectangular nuclear reactor
!! x_zero       real        x coordinate of the detector's position
!! y_zero       real        y coordinate of the detector's position
!! r_min        real        minimum radius of the hollow sphere
!! r_max        real        maximum radius of the hollow sphere
!! r_step       real        increment size for the radius of the hollow sphere
!! n_grid       integer     number of grid points in each dimension of Boole's integration
!! m_samples    integer     number of sample points in the Monte Carlo integration
!-----------------------------------------------------------------------
subroutine write_advanced_flux(depth, width, height, x_zero, y_zero, r_min, r_max, r_step, n_grid, m_samples)
     implicit none
     real(dp), intent(in) :: depth, width, height, x_zero, y_zero, r_min, r_max, r_step
     integer, intent(in) :: n_grid, m_samples

     real(dp) :: radius, box_booles, hollow_booles, hollow_mc, sigma_hollow
     character(len=*), parameter :: file_name = 'results_advanced.dat'
     integer :: unit

     open(newunit=unit, file=file_name)
     write(unit,'(5a28)') 'radius', 'box booles', 'hollow booles', 'hollow monte carlo', 'MC uncertainty'
     radius = r_min


    do
        if(radius .GT. r_max) exit ! Stops do loop when radius is greater than r_max

        box_booles = box_flux_booles(depth, width, height, x_zero, y_zero, n_grid)!box when full
        hollow_booles = total_flux_booles(depth, width, height, radius, x_zero, y_zero, n_grid)
        call hollow_box_flux_mc(depth, width, height, radius, x_zero, y_zero, m_samples, hollow_mc, sigma_hollow) !compute montecarlo
        write(unit,'(5e28.16)') radius, box_booles, hollow_booles, hollow_mc, sigma_hollow !write to file

        !update radius
        radius = radius + r_step
    enddo
    close(unit)

    print *, 'The fluxes were written in the '//file_name//' file'


!     ! The goal is to compare Boole's and Monte Carlo integration when there's a hollow 
!     ! sphere inside the reactor to the calculation of the solid reactor using Boole's method
!     ! Base the rest of the subroutine on the one from the basic part of the project.
 end subroutine write_advanced_flux

end module read_write
