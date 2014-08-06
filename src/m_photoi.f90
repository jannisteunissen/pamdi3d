module m_photoi
  use m_random

  implicit none
  private

  integer, parameter                     :: dp = kind(0.0d0)

  real(dp)                               :: pi_quench_fac
  real(dp)                               :: pi_min_inv_abs_len, pi_max_inv_abs_len
  real(dp), dimension(:, :), allocatable :: pi_photo_eff_table
  type(RNG_t)                            :: pi_rng

  ! Public methods
  public :: pi_initialize
  public :: pi_from_ionization

contains

  subroutine pi_initialize(cfg)
    use m_gas
    use m_units_constants
    use m_config
    type(CFG_t), intent(in) :: cfg
    integer                 :: t_size, t_size_2
    real(dp)                :: frac_O2, temp_vec(2)

    frac_O2 = GAS_get_fraction("O2")
    if (frac_O2 <= epsilon(1.0_dp)) then
       print *, "There is no oxygen, you should disable photoionzation"
       stop
    end if

    call CFG_get(cfg, "photoi_absorp_inv_lengths", temp_vec)
    pi_min_inv_abs_len = temp_vec(1) * frac_O2 * GAS_get_pressure()
    pi_max_inv_abs_len = temp_vec(2) * frac_O2 * GAS_get_pressure()

    ! print *, "Max abs. length photoi.", 1.0d3 / pi_min_inv_abs_len, "mm"
    ! print *, "Min abs. length photoi.", 1.0d3 / pi_max_inv_abs_len, "mm"

    pi_quench_fac = (30.0D0 * UC_torr_to_bar) / &
         (GAS_get_pressure() + (30.0D0 * UC_torr_to_bar))

    call CFG_get_size(cfg, "photoi_efficiency_table", t_size)
    call CFG_get_size(cfg, "photoi_efield_table", t_size_2)
    if (t_size_2 /= t_size) then
       print *, "size(photoi_efield_table) /= size(photoi_efficiency_table)"
       stop
    end if

    allocate(pi_photo_eff_table(2, t_size))
    call CFG_get(cfg, "photoi_efield_table", pi_photo_eff_table(1,:))
    call CFG_get(cfg, "photoi_efficiency_table", pi_photo_eff_table(2,:))
  end subroutine pi_initialize

  subroutine pi_from_ionization(my_part, photons)
    use m_particle_core
    use m_units_constants
    type(PC_part_t), intent(in)          :: my_part
    real(dp), allocatable, intent(inout) :: photons(:,:)
    real(dp)                             :: mean_gammas, en_frac, fly_len
    real(dp)                             :: fld, psi, chi, x_end(3)
    integer                              :: n, n_photons

    fld         = norm2(my_part%a / UC_elec_q_over_m)
    mean_gammas = get_photoi_eff(fld) * my_part%w * pi_quench_fac
    n_photons   = pi_rng%poisson(mean_gammas)

    if (allocated(photons)) deallocate(photons)
    allocate(photons(3, n_photons))

    do n = 1, n_photons
       ! Select random direction and absorption length
       en_frac  = pi_rng%uni_01()
       fly_len  = -log(1.0_dp - pi_rng%uni_01()) / get_photoi_lambda(en_frac)
       psi      = 2 * UC_pi * pi_rng%uni_01()
       chi      = acos(1.0_dp - 2 * pi_rng%uni_01())

       x_end(1) = my_part%x(1) + fly_len * sin(chi) * cos(psi)
       x_end(2) = my_part%x(2) + fly_len * sin(chi) * sin(psi)
       x_end(3) = my_part%x(3) + fly_len * cos(chi)
       photons(:, n) = x_end
    end do
  end subroutine pi_from_ionization

  ! Returns the photo-efficiency coefficient corresponding to an electric
  ! field of strength fld
  real(dp) function get_photoi_eff(fld)
    use m_lookup_table
    real(dp), intent(in) :: fld
    call LT_lin_interp_list(pi_photo_eff_table(1,:), &
         pi_photo_eff_table(2,:), fld, get_photoi_eff)
  end function get_photoi_eff

   ! Returns the inverse mean free path for a photon.
   real(dp) function get_photoi_lambda(en_frac)
      real(dp), intent(in) :: en_frac
      get_photoi_lambda = pi_min_inv_abs_len * &
           (pi_max_inv_abs_len/pi_min_inv_abs_len)**en_frac
   end function
end module m_photoi