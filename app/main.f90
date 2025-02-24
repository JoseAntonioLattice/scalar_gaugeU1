program lambda_phi4_2d

  use parameters
  use field
  use pbc
  use dynamics
  use observables
  implicit none

  integer :: i, save_unit, term_unit, obs_unit, conf_unit, avr_unit
  real(dp) :: msq(1), msqi = -10.0_dp, msqf = 10.0_dp

  print"(a)", "-----------------------------------"
  print"(a)", "|2 dimensional lambda phi^4 theory|"
  print"(a)", "-----------------------------------"
  msq = [1.0_dp]![(msqi + ((msqf - msqi)/(size(msq) - 1))*i, i = 0, size(msq) - 1 )]
  
  call read_input()
  call allocate_pbc(Lx)
  allocate(magnetization(Nmeasurements),mag1(Nmeasurements))
  allocate(phi(Lx,Lx))
  allocate(correlation(Nmeasurements,Lx))
  
  call hot_start(phi)
  !call cold_start(phi)
  
  open(newunit = save_unit, file = "data/unthermalized_configurations.dat")
  open(newunit = conf_unit, file = "data/thermalized_configurations.dat")
  open(newunit = term_unit,  file = "data/thermalization.dat")
  open(newunit = obs_unit,  file = "data/data.dat")
  open(newunit = avr_unit, file = "data/avr.dat")
  
  do i = 1, size(msq)
     !call save_configuration(phi,save_unit)
     call thermalization(phi,Nthermalization,algorithm,[msq(i),lambda],save_unit,term_unit)
     call take_measurements(phi, Nmeasurements, Nskip, algorithm, [m7878sq(i),lambda], obs_unit, conf_unit)
     print*, msq(i), sum(abs(magnetization))/size(magnetization)
  end do

  !do i = 1, Lx
  !   write(avr_unit,*) avr(correlation(:,i)), std_error(correlation(:,i)) !- (sum(mag1)/Nmeasurements)**2
  !end do

contains

  function avr(x)
    real(dp), intent(in) :: x(:)
    real(dp) :: avr

    avr = sum(x)/size(x)
    
  end function avr

  function std_error(x)
    real(dp), intent(in) :: x(:)
    real(dp) :: std_error, s, mean

    mean = avr(x)
    s = sqrt(sum((x - mean)**2)/(size(x) - 1))

    std_error = s/sqrt(real(size(x),dp))
    
  end function std_error
  
end program lambda_phi4_2d
