module pbc

  use precision
  implicit none

  integer(i4), allocatable, dimension(:) :: ip, im

contains

  subroutine allocate_pbc(L)
    integer(i4), intent(in) :: L

    allocate(ip(L),im(L))
    call set_pbc(L)
  end subroutine allocate_pbc

  subroutine set_pbc(L)
    integer(i4), intent(in) :: L
    integer(i4) :: i

    do i = 1, L
       ip(i) = i + 1
       im(i) = i - 1
    end do
    ip(L) = 1
    im(1) = L
  end subroutine set_pbc

   function ipf(vector, mu)

    integer(i4), dimension(:), intent(in) :: vector
    integer(i4) :: mu

    integer(i4), dimension(size(vector)) :: ipf
    
    ipf = vector

    ipf(mu) = ip(vector(mu))
    
  end function ipf

  function imf(vector, mu)

    integer(i4), dimension(:), intent(in) :: vector
    integer(i4) :: mu

    integer(i4), dimension(size(vector)) :: imf
    
    imf = vector

    imf(mu) = im(vector(mu))
    
  end function imf
 
  
end module pbc
