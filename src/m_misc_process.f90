! Copyright 2005-2012, Chao Li, Margreet Nool, Anbang Sun, Jannis Teunissen
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

module m_misc_process

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)
  integer  :: MISC_sum_detach

  real(dp) :: MISC_N2_bgdens
  real(dp) :: MISC_O2_bgdens
  real(dp) :: MISC_dt
  real(dp) :: MISC_gamma_decay_time

  public :: MISC_initialize
  public :: MISC_detachment
  public :: MISC_gamma

contains

  subroutine MISC_initialize(cfg)
    use m_config
    use m_gas
    type(CFG_t), intent(inout) :: cfg
    logical :: use_detach

    call CFG_get(cfg, "sim_use_o2m_detach", use_detach)
    call CFG_get(cfg, "gamma_decay_time", MISC_gamma_decay_time)

    if (use_detach) then
       MISC_O2_bgdens = GAS_number_dens * GAS_get_fraction("O2")
       MISC_N2_bgdens = GAS_number_dens * GAS_get_fraction("N2")
    end if
  end subroutine MISC_initialize

  subroutine MISC_detachment(pc, rng, dt, myrank, root)
    use m_efield_amr
    use m_particle_core
    use m_random
    real(dp), intent(in) :: dt
    integer, intent(in)  :: myrank, root
    class(PC_t), intent(inout) :: pc
    type(RNG_t), intent(inout) :: rng

    ! Set global vars
    MISC_dt         = dt
    MISC_sum_detach = 0

    call E_collect_mpi((/E_i_O2m/), myrank, root)

    if (myrank == root) then
       call E_loop_over_grids(pc, rng, update_detachment)
       print *, "Number of detachments", MISC_sum_detach
    end if
  end subroutine MISC_detachment

  subroutine update_detachment(pc, rng, amr_grid)
    use m_efield_amr
    use m_particle_core
    use m_particle
    use m_random
    use m_phys_domain
    class(PC_t), intent(inout)       :: pc
    type(RNG_t), intent(inout)      :: rng
    type(amr_grid_t), intent(inout) :: amr_grid
    real(dp), allocatable           :: loss(:,:,:)
    logical, allocatable            :: child_region(:, :, :)
    integer                         :: i, j, k, Nx, Ny, Nz
    integer                         :: n, nc, i_min(3), i_max(3), n_detach
    real(dp)                        :: xyz(3), e_str

    Nx = amr_grid%Nr(1)
    Ny = amr_grid%Nr(2)
    Nz = amr_grid%Nr(3)

    ! Create mask where children are
    allocate(child_region(Nx, Ny, Nz))
    allocate(loss(Nx, Ny, Nz))
    child_region = .false.

    do nc = 1, amr_grid%n_child
       i_min = E_xyz_to_ix(amr_grid, amr_grid%children(nc)%r_min)
       i_max = E_xyz_to_ix(amr_grid, amr_grid%children(nc)%r_max)
       child_region(i_min(1):i_max(1), i_min(2):i_max(2), i_min(3):i_max(3)) = .true.
    end do

    do k = 1, Nz
       do j = 1, Ny
          do i = 1, Nx
             xyz           = E_ix_to_xyz(amr_grid, (/i, j, k/))
             e_str         = norm2(E_get_field(xyz))
             loss(i, j, k) = MISC_get_O2m_loss(amr_grid%vars(i,j,k, E_i_O2m), e_str)
          end do
       end do
    end do

    amr_grid%vars(:,:,:, E_i_O2m) = amr_grid%vars(:,:,:, E_i_O2m) - loss

    where (child_region)
       loss = 0.0_dp
    elsewhere
       loss = loss * product(amr_grid%dr)
    end where

    do k = 1, Nz
       do j = 1, Ny
          do i = 1, Nx
             if (loss(i,j,k) <= epsilon(1.0_dp)) cycle

             xyz             = E_ix_to_xyz(amr_grid, (/i, j, k/))
             n_detach        = rng%poisson(loss(i,j,k))
             MISC_sum_detach = MISC_sum_detach + n_detach ! Global variable :(

             do n = 1, n_detach
                if (.not. PD_outside_domain(xyz)) &
                     call PM_create_ei_pair(pc, xyz)
             end do
          end do
       end do
    end do

  end subroutine update_detachment

  real(dp) function MISC_get_O2m_loss(O2min_dens, e_str)
    use m_gas
    use m_units_constants
    real(dp), intent(in) :: O2min_dens, e_str
    real(dp) :: T_eff

    T_eff = GAS_temperature + UC_N2_mass / (3.0d0 * UC_boltzmann_const) * (e_str * 2.0d-4)**2
    MISC_get_O2m_loss = 1.9d-12 * sqrt(T_eff/3.0d2) * exp(-4.990e3/T_eff) * MISC_N2_bgdens * 1.0d-6 + &
         2.7d-10 * sqrt(T_eff/3.0d2) * exp(-5.590e3/T_eff) * MISC_O2_bgdens * 1.0d-6
    MISC_get_O2m_loss = MISC_get_O2m_loss * O2min_dens * MISC_dt
  end function MISC_get_O2m_loss

  subroutine MISC_gamma(dt)
    use m_efield_amr
    real(dp), intent(in) :: dt
    real(dp)             :: factor

    factor = exp(-dt / MISC_gamma_decay_time)

    call E_multiply_grids(E_i_gamma, factor)

  end subroutine MISC_gamma

end module m_misc_process
