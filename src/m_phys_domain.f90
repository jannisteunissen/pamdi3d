module m_phys_domain
  
  implicit none
  public
  
  integer, parameter, private :: dp = kind(0.0d0)

  real(dp) :: PD_r_max(3)
  real(dp) :: PD_dr(3)
  integer  :: PD_size(3)
  logical :: PD_use_elec

contains
  
  subroutine PD_set(cfg)
    use m_config
    type(CFG_t), intent(in) :: cfg
    call CFG_get(cfg, "grid_size", PD_size)
    call CFG_get(cfg, "grid_delta", PD_dr)
    call CFG_get(cfg, "sim_use_electrode", PD_use_elec)
    PD_r_max = PD_size * PD_dr
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