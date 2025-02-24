module functions

  use precision
  use pbc, only: ip, im
  implicit none

  real(dp), parameter :: pi = acos(-1.0_dp)

contains

  function lagrangian(phi,x,msq,lambda)
    complex(dp), intent(in) :: phi(:,:)
    integer(i4), intent(in) :: x(2)
    real(dp), intent(in) :: msq, lambda

    real(dp) :: lagrangian

    lagrangian = 0.5 * ( (abs(phi(ip(x(1)),x(2)) - phi(x(1),x(2))) )**2  &
                       + (abs(phi(x(1),ip(x(2))) - phi(x(1),x(2))) )**2 )&
               + 0.5 * msq * (abs(phi(x(1),x(2))))**2 &
               + 0.25 * lambda * (abs(phi(x(1),x(2))))**4 
    
  end function lagrangian
  
  function action(phi,msq,lambda)
    complex(dp), intent(in) :: phi(:,:)
    real(dp), intent(in) :: msq, lambda
    integer(i4) :: x1, x2
    real(dp) :: action

    action  =  0.0_dp

    do x1 = 1, size(phi(:,1))
       do x2 = 1, size(phi(1,:))
          action = action + lagrangian(phi,[x1,x2],msq,lambda)
       end do
    end do

  end function action

  function DS(phi,phi_p,x,msq,lambda)
    use parameters, only : epsilon
    complex(dp), intent(in) :: phi(:,:), phi_p
    integer(i4), intent(in) :: x(2)
    real(dp), intent(in) ::  msq, lambda
    real(dp) :: r
    real(dp) :: DS

    DS = 2 * (abs(phi_p)**2 - abs(phi(x(1),x(2)))**2) &
    - real((phi_p - phi(x(1),x(2)))*conjg(neighbors(phi,x)))&
    + 0.5 * msq * ( abs(phi_p)**2 - abs(phi(x(1),x(2)))**2 ) &
    + 0.25 * lambda * ( abs(phi_p)**4 - abs(phi(x(1),x(2)))**4 )
    
  end function DS

  function neighbors(phi,x)
    complex(dp) :: neighbors
    complex(dp), dimension(:,:), intent(in) :: phi
    integer, dimension(2), intent(in) :: x
    integer :: mu, xp(2)

    neighbors = phi(ip(x(1)),x(2)) + phi(x(1),ip(x(2))) &
              + phi(im(x(1)),x(2)) + phi(x(1),im(x(2))) 
    
  end function neighbors
  
  function DS2(phi,phi_p,x,msq,lambda)
    use parameters, only : epsilon
    complex(dp), intent(in) :: phi(:,:), phi_p
    integer(i4), intent(in) :: x(2)
    real(dp), intent(in) :: msq, lambda
    real(dp) :: r, DS2
    complex(dp) ::  phip(size(phi(:,1)),size(phi(1,:)))

    phip = phi
    phip(x(1),x(2)) = phi_p
    DS2 = action(phip,msq,lambda) - action(phi,msq,lambda)
  end function DS2


end module functions
