# Nuclear Reactor Program
This program computes the flux of a nuclear reactor using Boole's
quadrature and Monte Carlo quadratures.

The basic part of the program determines the flux of the nuclear 
reactor for varying values of detector position in the x direction.
The large x0 approximation flux is also computed for the same x0 values.
The fluxes are written in the 'results_basic.dat' file

The advanced part of the program computes the flux of a hollow
spherical nuclear reactor for a detector at a fixed position.
In this section, the radius of the sphere is varied. The fluxes
are written in the 'results_advanced.dat' file.

## File Descriptions
main.f90 calls read_input and write_output subroutines from
    read_write.f90 for both basic and advanced parts.

read_write.f90 prompts the user to provide inputs for depth, width
    height, xmin, xmax, y0, rmin, rmax, number of steps and number 
    of sample points for monte carlo integration. The code 
    determines if the values are positive and nonzero, before
    saving to the variables. the code then computes the fluxes 
    using booles and monte carlo quadratures and writes them to
    files 'results_basic.dat' and 'results_advanced.dat'

neutron_flux.f90 computes the flux of the nuclear reactor using 
    boole's and monte carlo quadratures. The program calls 
    functions and subroutines from quadrature.f90 to compute flux
    kernels for the integration are included in this file

quadrature.f90 contains the functions and subroutines for boole's
    and monte carlo integration 

makefile compiles main.f90, read_write.f90, neutron_flux.f90, and 
    quadrature.f90 and writes program 'nuclear-reactor'

plots.ipynb calls values from files 'results_advanced.dat' and 
    'results_basic.dat' and plots the data and error

## Running the Code

To compile, type in terminal:

```$ make ```

The code will then compile and create the program 'nuclear-reactor'.

To run the program, type in terminal:
```$ ./nuclear-reactor```

The program will then write the fluxes to 'results_basic.dat' and 
'results_advanced.dat'

To view the plots, open 'plots.ipynb' and select 'Run All'.

To clear output files, type in terminal:

```$make clean```
