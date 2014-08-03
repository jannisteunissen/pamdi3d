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

   real(dp) :: MISC_tau_excited, MISC_quench_fac
   real(dp) :: MISC_min_inv_abs_len, MISC_max_inv_abs_len, MISC_frac_O2

   integer :: MISC_table_size
   integer :: MISC_sum_detach
   integer :: MISC_sum_photons

   ! Global variables
   real(dp) :: MISC_loss_frac
   real(dp) :: MISC_N2_bgdens
   real(dp) :: MISC_O2_bgdens
   real(dp) :: MISC_dt

   real(dp), dimension(:, :), allocatable :: MISC_photo_eff_table

   public :: MISC_photoionization
   public :: MISC_initialize
   public :: MISC_detachment

contains

   subroutine MISC_initialize(cfg)
      use m_gas
      use m_units_constants
      use m_config
      type(CFG_t), intent(in) :: cfg
      real(dp) :: temp_vec(2)
      MISC_tau_excited = CFG_get_real(cfg, "PI_meanLifeTimeExcited")

      MISC_frac_O2 = GAS_get_fraction("O2")

      if (MISC_frac_O2 <= epsilon(1.0_dp)) then
         print *, "There is no oxygen, you should disable photoionzation"
         stop
      end if

      MISC_O2_bgdens = GAS_get_number_dens() * GAS_get_fraction("O2")
      MISC_N2_bgdens = GAS_get_number_dens() * GAS_get_fraction("N2")

      call CFG_get_array(cfg, "PI_absorpInvLengths", temp_vec)
      MISC_min_inv_abs_len = temp_vec(1) * MISC_frac_O2 * GAS_get_pressure()
      MISC_max_inv_abs_len = temp_vec(2) * MISC_frac_O2 * GAS_get_pressure()

      print *, "Max absorbp. length for photoionization ", 1.0d3 / MISC_min_inv_abs_len, "mm"
      print *, "Min absorbp. length for photoionization ", 1.0d3 / MISC_max_inv_abs_len, "mm"

      MISC_quench_fac = (30.0D0 * UC_torr_to_bar) / (GAS_get_pressure() + (30.0D0 * UC_torr_to_bar))
      MISC_table_size = CFG_get_size(cfg, "PI_EfieldTable")

      if (MISC_table_size /= CFG_get_size(cfg, "PI_efficiencyTable")) then
         print *, "Make sure MISC_efficiencyTable and MISC_EfieldTable have the same size"
         stop
      end if

      allocate( MISC_photo_eff_table(2, MISC_table_size) )

      call CFG_get_array(cfg, "PI_EfieldTable", MISC_photo_eff_table(1,:) )
      call CFG_get_array(cfg, "PI_efficiencyTable", MISC_photo_eff_table(2,:) )

   end subroutine MISC_initialize

   real(dp) function findPhotoEff(E_f)
     use m_lookup_table
      ! Returns the photo-efficiency coefficient corresponding to an electric
      ! field of strength E_f
      real(dp), intent(IN) :: E_f
      call LT_lin_interp_list(MISC_photo_eff_table(1,:), &
           MISC_photo_eff_table(2,:), E_f, findPhotoEff)
   end function findPhotoEff

   real(dp) function Pho_kf(energy_frac)
      ! Returns the inverse mean free path for a photon.
      real(dp), intent(IN) :: energy_frac
      Pho_kf = MISC_min_inv_abs_len * (MISC_max_inv_abs_len/MISC_min_inv_abs_len)**energy_frac
   end function

   subroutine MISC_photoionization(pc, rng, dt, myrank, root)
      use m_efield_amr
      use m_particle_core
      use m_random
      type(PC_t), intent(inout) :: pc
      type(RNG_t), intent(inout) :: rng
      real(dp), intent(in) :: dt
      integer, intent(in)  :: myrank, root

      ! Set global vars
      MISC_loss_frac   = 1.0_dp - exp(-dt / MISC_tau_excited)
      MISC_sum_photons = 0

      call E_collect_mpi((/E_i_exc/), myrank, root)

      if (myrank == root) then
         call E_loop_over_grids(pc, rng, update_excited)
         print *, "Number of ionizing photons:", MISC_sum_photons
      end if
   end subroutine MISC_photoionization

   subroutine update_excited(pc, rng, amr_grid)
      use m_efield_amr
      use m_phys_domain
      use m_random
      use m_particle_core
      type(PC_t), intent(inout) :: pc
      type(RNG_t), intent(inout) :: rng
      type(amr_grid_t), intent(inout) :: amr_grid
      real(dp), allocatable           :: loss(:,:,:)
      logical, allocatable            :: child_region(:, :, :)
      integer                         :: i, j, k, Nx, Ny, Nz
      integer                         :: n, nc, i_min(3), i_max(3), n_photons
      real(dp)                        :: xyz(3), x_end(3), flylen, chi, psi, e_str, energy_frac
      real(dp) :: v(3), a(3), w, t_left

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

      loss = amr_grid%vars(:,:,:, E_i_exc) * MISC_loss_frac ! Global variable :(
      amr_grid%vars(:,:,:, E_i_exc) = amr_grid%vars(:,:,:, E_i_exc) - loss

      where (child_region)
         loss = 0.0_dp
      elsewhere
         loss = loss * product(amr_grid%dr)
      end where

      do k = 1, Nz
         do j = 1, Ny
            do i = 1, Nx
               if (loss(i,j,k) <= epsilon(1.0_dp)) cycle

               xyz            = E_ix_to_xyz(amr_grid, (/i, j, k/))
               e_str          = norm2(E_get_field(xyz))
               n_photons      = rng%poisson(findPhotoEff(e_str) * loss(i,j,k) * MISC_quench_fac)
               MISC_sum_photons = MISC_sum_photons + n_photons ! Global variable :(

               do n = 1, n_photons
                  energy_frac = rng%uni_01()
                  flylen      = -log(1.0_dp - rng%uni_01()) / Pho_kf(energy_frac)
                  psi         = 2 * acos(-1.0_dp) * rng%uni_01()
                  chi         = acos(1.0_dp - 2 * rng%uni_01())

                  x_end(1)    = xyz(1) + flylen * sin(chi) * cos(psi)
                  x_end(2)    = xyz(2) + flylen * sin(chi) * sin(psi)
                  x_end(3)    = xyz(3) + flylen * cos(chi)

                  v      = 0
                  a      = 0
                  w      = 1
                  t_left = 0

                   ! Create electron-ion pair with excess enery of photon
                  if (.not. PD_outside_domain(x_end)) &
                       call pc%create_part(x_end, v, a, w, t_left)
               end do
            end do
         end do
      end do
   end subroutine update_excited

   subroutine MISC_detachment(pc, rng, dt, myrank, root)
     use m_efield_amr
     use m_particle_core
     use m_random
     real(dp), intent(in) :: dt
     integer, intent(in)  :: myrank, root
     type(PC_t), intent(inout) :: pc
     type(RNG_t), intent(inout) :: rng
     ! Set global vars
     MISC_dt = dt
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
      use m_random
      use m_phys_domain
      type(PC_t), intent(inout) :: pc
      type(RNG_t), intent(inout) :: rng
      type(amr_grid_t), intent(inout) :: amr_grid
      real(dp), allocatable           :: loss(:,:,:)
      logical, allocatable            :: child_region(:, :, :)
      integer                         :: i, j, k, Nx, Ny, Nz
      integer                         :: n, nc, i_min(3), i_max(3), n_detach
      real(dp)                        :: xyz(3), e_str
      real(dp) :: v(3), a(3), w, t_left

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

      amr_grid%vars(:,:,:, E_i_exc) = amr_grid%vars(:,:,:, E_i_exc) - loss

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
                  ! Give a small energy up to 1 eV
                  v = 0
                  a = 0
                  w = 1
                  t_left = 0
                  if (.not. PD_outside_domain(xyz)) &
                       call pc%create_part(xyz, v, a, w, t_left)
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

      T_eff = GAS_get_temperature() + UC_N2_mass / (3.0d0 * UC_boltzmann_const) * (e_str * 2.0d-4)**2
      MISC_get_O2m_loss = 1.9d-12 * sqrt(T_eff/3.0d2) * exp(-4.990e3/T_eff) * MISC_N2_bgdens * 1.0d-6 + &
           2.7d-10 * sqrt(T_eff/3.0d2) * exp(-5.590e3/T_eff) * MISC_O2_bgdens * 1.0d-6
      MISC_get_O2m_loss = MISC_get_O2m_loss * O2min_dens * MISC_dt
   end function MISC_get_O2m_loss

 end module m_misc_process
