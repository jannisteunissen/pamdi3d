! Jannis Teunissen, Casper Rutjes
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

!> Module that adds photo ionization functionallity to particle core
module m_photoi
  use m_random
  use m_particle_core
  use m_units_constants
  use m_gas

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  ! photo ionization settings
  real(dp), protected              :: pi_quench_fac
  real(dp), protected              :: pi_min_inv_abs_len, pi_max_inv_abs_len
  real(dp), allocatable, protected :: pi_photo_eff_table(:, :)

  public :: PI_initialize
  public :: PI_do_photoi

contains

  !> After the first argument, if argument is present then overwrite default
  subroutine PI_initialize(quench_fac_shift,     &
                           min_inv_abs_len_resc, &
                           max_inv_abs_len_resc, &
                           size_photo_eff_table, &
                           photo_eff_table)
   real(dp), intent(in), optional              :: quench_fac_shift
   real(dp), intent(in), optional              :: min_inv_abs_len_resc
   real(dp), intent(in), optional              :: max_inv_abs_len_resc
   integer, intent(in), optional               :: size_photo_eff_table
   real(dp), intent(in), allocatable,optional  :: photo_eff_table(:,:)

   !> Check if GAS is initialized
   if (GAS_initialized .eqv. .false.) then
     print *, "Initialize GAS before PI"
     stop
   end if

    !> Check oxygen concentration
    if (GAS_get_fraction("O2") <= epsilon(1.0_dp)) then
       print *, "There is no oxygen, you should disable photoionization"
       stop
    end if

    pi_quench_fac = (30.0D0 * UC_torr_to_bar) / &
         (GAS_pressure + (30.0D0 * UC_torr_to_bar))
    if (present(quench_fac_shift)) then
      pi_quench_fac = (quench_fac_shift * UC_torr_to_bar) / &
           (GAS_pressure + (quench_fac_shift * UC_torr_to_bar))
    end if

    if (present(min_inv_abs_len_resc)) then
      pi_min_inv_abs_len = min_inv_abs_len_resc * GAS_get_fraction("O2") * GAS_pressure
    else
      pi_min_inv_abs_len = 3.5D0 * GAS_get_fraction("O2") * GAS_pressure
    end if

    if (present(max_inv_abs_len_resc)) then
      pi_max_inv_abs_len = max_inv_abs_len_resc * GAS_get_fraction("O2") * GAS_pressure
    else
      pi_max_inv_abs_len = 200D0 * GAS_get_fraction("O2") * GAS_pressure
    end if

    if (present(size_photo_eff_table) .and. present(photo_eff_table)) then
      allocate(pi_photo_eff_table(2,size_photo_eff_table))
      pi_photo_eff_table = photo_eff_table
    else
        allocate(pi_photo_eff_table(2,6))
        pi_photo_eff_table(1,:) = [0.0D0, 0.25D7,  0.4D7, 0.75D7,  1.5D7, 3.75D7]
        pi_photo_eff_table(2,:) = [0.0D0, 0.05D0, 0.12D0, 0.08D0, 0.06D0, 0.04D0]
    end if

    write(*,'(A,F8.1,A,F8.1,A)') "** Photo ionization is initialized **"
    write(*,'(A,F8.1,A,F8.1,A)') "The photon mean free path (non-uniformly distributed) ranges between:", &
    1.0/get_photoi_lambda(1.0d0)*1.0d3," mm",&
    1.0/get_photoi_lambda(0.0d0)*1.0d3," mm"

    write(*,'(A,F9.3,A,F9.3,A)') "The mean photon ionization conversion rate (efficiency times quench_fac) ranges between:", &
    minval(pi_photo_eff_table(2,:)) * pi_quench_fac * 100,"% and",&
    maxval(pi_photo_eff_table(2,:)) * pi_quench_fac * 100,"%"
  end subroutine PI_initialize

  !> Ionization callback procedure
  subroutine PI_do_photoi(pc, my_part, c_ix, c_type)
    class(PC_t), intent(inout)  :: pc
    type(PC_part_t), intent(in) :: my_part
    integer, intent(in)         :: c_ix, c_type

    ! dummy
    type(PC_part_t)                      :: my_new_part
    real(dp)                             :: mean_gammas, en_frac, fly_len
    real(dp)                             :: fld, psi, chi
    integer                              :: n, n_photons

    fld         = norm2(my_part%a / UC_elec_q_over_m)
    mean_gammas = get_photoi_eff(fld) * my_part%w * pi_quench_fac
    n_photons   = pc%rng%poisson(mean_gammas)

    do n = 1, n_photons
       ! Select random direction and absorption length
       en_frac  = pc%rng%unif_01()
       fly_len  = -log(1.0_dp - pc%rng%unif_01()) / get_photoi_lambda(en_frac)
       psi      = 2 * UC_pi * pc%rng%unif_01()
       chi      = acos(1.0_dp - 2 * pc%rng%unif_01())

       my_new_part%x(1)   = my_part%x(1) + fly_len * sin(chi) * cos(psi)
       my_new_part%x(2)   = my_part%x(2) + fly_len * sin(chi) * sin(psi)
       my_new_part%x(3)   = my_part%x(3) + fly_len * cos(chi)
       my_new_part%v      = [0.0D0, 0.0D0, 0.0D0]
       my_new_part%a      = my_part%a
       my_new_part%w      = 1.0D0
       my_new_part%t_left = 0.0D0

       if (pc%outside_check(my_new_part) .eqv. .false.) then
         call pc%add_part(my_new_part)
       end if
    end do
  end subroutine PI_do_photoi

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
   end function get_photoi_lambda

end module m_photoi
