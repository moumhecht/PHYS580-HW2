! Program: nuclear_reactor
! By: Malida Hecht
!-----------------------------------------------------------------------------
! This program computes the flux of a nuclear reactor.
!
! The basic part of the program computes the flux of the nuclear reactor where
! the position of the dector is varied over distances.
! 
! The program takes user inputs for depth, width, height, y0, xmin, xmax, ngrid,
! and samples.
! 
! These values are then used to compute the flux using Boole's quadrature and
! Monte Carlo quadrature techniques. The fluxes are then written to 'results.dat'
! The basic part also uses Boole's quadrature to compute the large x0 approximation
! For flux where the reactor is treated as a point source, at varying distances for
! the detector position.
!
! The advanced program computes the flux of the nuclear reactor computes the flux
! of a nuclear reactor with a hollow sphere in its center using Booles quadrature and 
! Monte Carlo integration. The user is prompted to provide input for x0, rmin, rmax,
! and rstep. The program also computes the flux for a solid box. The fluxes are recorded
! in the file 'results_advanced.dat'
!
! Basic part of the project
! depth = 40.0_dp
! width = 100.0_dp
! height = 60.0_dp
! y_zero = 30._dp
! x_min = 5_dp
! x_max = 200._dp
! x_step = 5._dp
! n_grid = 25
! m_samples = 10000

! Advanced part of the project
! x_zero = 80._dp
! r_min = 1.0_dp
! r_max = 19._dp
! r_step = 1.0_dp
!-----------------------------------------------------------------------------
program nuclear_reactor
use types 

use read_write, only : read_input, write_neutron_flux, read_advanced_input, write_advanced_flux

implicit none

real(dp) :: depth, width, height, y_zero, x_min, x_max, x_step
integer :: n_grid, m_samples
real(dp) :: x_zero, r_min, r_max, r_step



! Basic part of the project
call read_input(depth, width, height, y_zero, x_min, x_max, x_step, n_grid, m_samples)
call write_neutron_flux(depth, width, height, y_zero, x_min, x_max, x_step, n_grid, m_samples)

! Advanced part of the project
call read_advanced_input(x_zero, r_min, r_max, r_step)
call write_advanced_flux(depth, width, height, x_zero, y_zero, r_min, r_max, r_step, n_grid, m_samples)

end program nuclear_reactor