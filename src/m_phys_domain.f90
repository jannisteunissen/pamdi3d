module m_phys_domain

  implicit none
  public

  integer, parameter, private :: dp = kind(0.0d0)

  real(dp) :: PD_r_max(3)
  real(dp) :: PD_plasma_r_min(3)
  real(dp) :: PD_plasma_r_max(3)
  real(dp) :: PD_dr(3)
  integer  :: PD_size(3)
  logical :: PD_use_elec

contains

  subroutine PD_set(cfg)
    use m_config
    type(CFG_t), intent(inout) :: cfg
    call CFG_get(cfg, "grid_size", PD_size)
    call CFG_get(cfg, "grid_delta", PD_dr)
    call CFG_get(cfg, "elec_enabled", PD_use_elec)
    call CFG_get(cfg, "grid_plasma_min_rel_pos", PD_plasma_r_min)
    call CFG_get(cfg, "grid_plasma_max_rel_pos", PD_plasma_r_max)
    
    PD_r_max        = (PD_size-1) * PD_dr
    PD_plasma_r_min = PD_r_max * PD_plasma_r_min
    PD_plasma_r_max = PD_r_max * PD_plasma_r_max

  end subroutine PD_set

  ! Checks whether pos is outside the computational domain or inside the
  ! electrode, and if so returns .TRUE., otherwise returns .FALSE.
  logical function PD_outside_domain(pos)
    use m_electrode
    real(dp), intent(in) :: pos(3)

    ! Check whether any component of the position is outside of the domain
    PD_outside_domain = any(pos <= PD_plasma_r_min &
         .or. pos >= PD_plasma_r_max)

    ! Now if any component is outside the domain or inside the electrode,
    ! PD_outside_domain = .TRUE.
    if (PD_use_elec) then
       if (EL_inside_elec(pos)) then
          PD_outside_domain = .true.
       end if
    end if
  end function PD_outside_domain

end module m_phys_domain
