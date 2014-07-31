module m_phys_domain
  
  implicit none
  public
  
  integer, parameter, private :: dp = kind(0.0d0)

  real(dp) :: PD_r_max(3)
  real(dp) :: PD_dr(3)
  integer  :: PD_size(3)
  logical :: PD_use_elec

contains
  
  subroutine PD_set(rmax, npoints, use_elec)
    real(dp), intent(in) :: rmax(3)
    integer              :: npoints(3)
    logical :: use_elec
    
    PD_r_max = rmax
    PD_size  = npoints
    PD_dr    = rmax / (npoints-1)

    PD_use_elec = use_elec
  end subroutine PD_set

  ! Checks whether pos is outside the computational domain or inside the
  ! electrode, and if so returns .TRUE., otherwise returns .FALSE.
  logical function PD_outside_domain(pos)
    use m_electrode
    real(dp), intent(in) :: pos(3)

    ! Check whether any component of the position is outside of the domain
    PD_outside_domain = any(pos <= 0 .or. pos >= PD_r_max)

    ! Now if any component is outside the domain or inside the electrode,
    ! PD_outside_domain = .TRUE.
    if (PD_use_elec) then
       if (EL_inside_elec(pos)) then
          PD_outside_domain = .true.
       end if
    end if
  end function PD_outside_domain

end module m_phys_domain