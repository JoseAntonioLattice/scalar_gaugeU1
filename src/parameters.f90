module parameters

  use precision
  implicit none

  integer(i4) :: Lx, Ly

  integer(i4) :: Nthermalization, Nmeasurements, Nskip
  
  real(dp) :: lambda
  real(dp) :: epsilon
  real(dp) :: beta

  character(100) :: algorithm, start

  namelist /lattice/ Lx, Ly
  namelist /simulation_params/ start,Nthermalization,Nmeasurements, Nskip, algorithm
  namelist /pars/ lambda, epsilon,beta
  
  
contains

  subroutine read_input()

    character(100) :: parameters_file
    integer(i4) :: inunit
    write(*,"(a)") "Please enter the parameters file:"
    read(*,"(a)") parameters_file
    write(*,"(2a)") "User typed: ", trim(parameters_file)
    open(newunit = inunit, file = trim(parameters_file))

    read(inunit, nml = lattice)
    read(inunit, nml = simulation_params)
    read(inunit, nml = pars)

    close(inunit)

    write(*, nml = lattice)
    write(*, nml = simulation_params)
    write(*, nml = pars)
    
  end subroutine read_input

    
end module parameters
