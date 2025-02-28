module functions

  use precision
  use pbc, only: ip, im, ipf, imf
  implicit none

  real(dp), parameter :: pi = acos(-1.0_dp)

contains

  function lagrangian(phi,u,x,msq,lambda,beta)
    complex(dp), intent(in) :: phi(:,:), u(:,:,:)
    integer(i4), intent(in) :: x(2)
    real(dp), intent(in) :: msq, lambda, beta

    real(dp) :: lagrangian

    lagrangian = 2*abs(phi(x(1),x(2)))**2 - &
           real( &
                conjg(phi(x(1),x(2)))*  &
                ( &
                conjg(U(1,im(x(1)),x(2))) * phi(im(x(1)),x(2)) + U(1,x(1),x(2)) * phi(ip(x(1)),x(2)) +  &
                conjg(U(2,x(1),im(x(2)))) * phi(x(1),im(x(2))) + U(2,x(1),x(2)) * phi(x(1),ip(x(2)))  &
                ) &
               ) + &
               -beta*(1.0_dp - real(plaquette(u,x))) & 
               + 0.5 * msq * (abs(phi(x(1),x(2))))**2 &
               + 0.25 * lambda * (abs(phi(x(1),x(2))))**4 
    
  end function lagrangian

  function plaquette(u,x)

    complex(dp) :: plaquette
    
    complex(dp), dimension(:,:,:), intent(in) :: u
    integer(i4), intent(in) :: x(2)
    integer(i4), dimension(2) :: x2, x3


    x2 = ipf(x,1)
    x3 = ipf(x,2)
    
    plaquette = U(1,x(1),x(2)) * U(2,x2(1),x2(2)) * &
          conjg(U(1,x3(1),x3(2))) * conjg(U(2,x(1),x(2)))
    
  end function plaquette

  function action(phi,u,msq,lambda,beta)
    complex(dp), intent(in) :: phi(:,:), u(:,:,:)
    real(dp), intent(in) :: msq, lambda, beta
    integer(i4) :: x1, x2
    real(dp) :: action

    action  =  0.0_dp

    do x1 = 1, size(phi(:,1))
       do x2 = 1, size(phi(1,:))
          action = action + lagrangian(phi,u,[x1,x2],msq,lambda,beta)
       end do
    end do

  end function action

  function DS_scalar(phi,phi_p,u,x,msq,lambda)
    use parameters, only : epsilon
    complex(dp), intent(in) :: phi(:,:), phi_p, u(:,:,:)
    integer(i4), intent(in) :: x(2)
    real(dp), intent(in) ::  msq, lambda
    real(dp) :: r
    real(dp) :: DS_scalar

    DS_scalar = (2.0_dp + msq/2)*(abs(phi_p)**2 - abs(phi(x(1),x(2)))**2) &
         + 0.25*lambda*(abs(phi_p)**4 - abs(phi(x(1),x(2)))**4) &
         -2*real( &
                conjg(phi_p - phi(x(1),x(2)))*  &
                ( &
                conjg(U(1,im(x(1)),x(2))) * phi(im(x(1)),x(2)) + U(1,x(1),x(2)) * phi(ip(x(1)),x(2)) +  &
                conjg(U(2,x(1),im(x(2)))) * phi(x(1),im(x(2))) + U(2,x(1),x(2)) * phi(x(1),ip(x(2)))  &
                ) &
               )
  end function DS_SCALAR

  
  function DS_GAUGE(phi,u,unew,x,mu,beta)
    use parameters, only : epsilon
    complex(dp), intent(in) :: phi(:,:), u(:,:,:), unew
    integer(i4), intent(in) :: x(2), mu
    real(dp), intent(in) :: beta
    real(dp) :: r
    real(dp) :: DS_gauge
    integer(i4) :: x2(2)
    
    x2 = ipf(x,mu)
    DS_gauge = -beta * real( (unew - u(mu,x(1),x(2))) * conjg(staples(u,x,mu)) ) &
               -2*real((unew - u(mu,x(1),x(2)))*(conjg(phi(x(1),x(2)))*phi(x2(1),x2(2))))
    
  end function DS_GAUGE

  function staples(u,x,mu)

    complex(dp) :: staples
    complex(dp), dimension(:,:,:), intent(in) :: u
    integer(i4), intent(in) :: x(2), mu
    integer(i4), dimension(2) :: x2, x3, x4, x5, x6

    integer(i4) :: nu
    
    if ( mu == 1 ) then
       nu = 2
    elseif( mu == 2)then
       nu = 1
    end if

    x2 = ipf(x,nu)
    x3 = ipf(x,mu)
    x4 = imf(x,nu)
    x5 = x4
    x6 = imf(x3,nu)
    
    staples = u(nu,x(1),x(2)) * u(mu,x2(1), x2(2)) * conjg( u(nu,x3(1), x3(2)) ) + &
         conjg( u(nu,x4(1),x4(2)) ) * u(mu,x5(1), x5(2)) * u(nu,x6(1), x6(2))
    
  end function staples
  
  function neighbors(phi,x)
    complex(dp) :: neighbors
    complex(dp), dimension(:,:), intent(in) :: phi
    integer, dimension(2), intent(in) :: x
    integer :: mu, xp(2)

    neighbors = phi(ip(x(1)),x(2)) + phi(x(1),ip(x(2))) &
              + phi(im(x(1)),x(2)) + phi(x(1),im(x(2))) 
    
  end function neighbors
  
  !function DS2(phi,phi_p,x,msq,lambda)
  !  use parameters, only : epsilon
  !  complex(dp), intent(in) :: phi(:,:), phi_p
  !  integer(i4), intent(in) :: x(2)
  !  real(dp), intent(in) :: msq, lambda
  !  real(dp) :: r, DS2
  !  complex(dp) ::  phip(size(phi(:,1)),size(phi(1,:)))

   ! phip = phi
   ! phip(x(1),x(2)) = phi_p
   ! DS2 = action(phip,msq,lambda,beta) - action(phi,msq,lambda)
  !end function DS2


end module functions
