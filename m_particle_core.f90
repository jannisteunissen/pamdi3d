! Authors: Jannis Teunissen, concepts based on work of Chao Li, Margreet Nool
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

!> Core particle module. The user has to supply his own routines for customization.
module m_particle_core
   use m_lookup_table
   use m_linked_list
   use m_random

   implicit none
   private

   integer, parameter  :: dp               = kind(0.0d0)
   integer, parameter  :: PC_max_num_colls = 100
   real(dp), parameter :: PC_dead_weight   = -HUGE(1.0_dp)

   !> The particle type
   type PC_part_t
      real(dp) :: t_left   ! The time until the next timestep
      real(dp) :: x(3)
      real(dp) :: v(3)
      real(dp) :: a(3)
      real(dp) :: weight
   end type PC_part_t

   ! Collision type. TODO: rethink this
   type PC_coll_t
      integer                      :: num                    ! The number of different collisions
      real(dp)                     :: max_rate, inv_max_rate ! Maximum collision rate and inverse
      integer, allocatable         :: types(:)               ! The types of the collisions that can occur
      real(dp), allocatable        :: special_val(:)         ! Special value used by this collision
      type(LT_table_t)             :: rate_lt                ! Lookup table with collision rates
   end type PC_coll_t

   ! We then also need multiple RNG states, that should be separated by
   ! at least a few cache lines in memory
   type PC_rng_t
      type(RNG_state_t) :: state
      integer           :: separator(32) ! Large enough to be sure
   end type PC_rng_t

   type PC_list_t
      type(PC_part_t), allocatable :: parts(:)
      type(LL_int_head_t)          :: clean_list
      type(PC_coll_t)              :: colls
      real(dp)                     :: mass
      type(RNG_state_t)            :: rng
      integer                      :: n_part
   end type PC_list_t

   ! We can have multiple lists of particles, for OpenMP
   type(PC_list_t), allocatable :: PC_PL(:)

   interface
      subroutine if_ipart(my_part)
         import
         type(PC_part_t), intent(in) :: my_part
      end subroutine if_ipart

      subroutine if_ipart_oreals(my_part, my_reals)
         import
         type(PC_part_t), intent(in) :: my_part
         real(dp), intent(out)       :: my_reals(:)
      end subroutine if_ipart_oreals

      subroutine if_ipart_oreal(my_part, my_real)
         import
         type(PC_part_t), intent(in) :: my_part
         real(dp), intent(out)       :: my_real
      end subroutine if_ipart_oreal

      subroutine if_iopart(my_part)
         import
         type(PC_part_t), intent(inout) :: my_part
      end subroutine if_iopart

      subroutine if_iint2_ipart(my_part, my_int_1, my_int_2)
         import
         type(PC_part_t), intent(in) :: my_part
         integer, intent(in)         :: my_int_1, my_int_2
      end subroutine if_iint2_ipart

      subroutine if_ovec3(my_vec)
         import
         real(dp), intent(out) :: my_vec(3)
      end subroutine if_ovec3

      subroutine if_ipart_oreal3(my_part, my_vec)
         import
         type(PC_part_t), intent(in) :: my_part
         real(dp), intent(out)       :: my_vec(3)
      end subroutine if_ipart_oreal3

      real(dp) function if_freal_ipart(my_part)
         import
         type(PC_part_t), intent(in) :: my_part
      end function if_freal_ipart

      logical function if_filter_func(my_part, real_args)
         import
         type(PC_part_t), intent(in) :: my_part
         real(dp), intent(in) :: real_args(:)
      end function if_filter_func

      subroutine io_2part_rng(part_a, part_b, rng_state)
         import
         type(PC_part_t), intent(inout) :: part_a, part_b
         type(RNG_state_t), intent(inout) :: rng_state
      end subroutine io_2part_rng
   end interface

   procedure(if_iint2_ipart), pointer :: PC_pptr_coll_callback => null()
   procedure(if_ovec3), pointer :: PC_pptr_bg_vel_sampler => null()

   ! Types
   public :: PC_part_t
   public :: PC_coll_t

   ! Procedures
   public :: PC_initialize
   public :: PC_set_pptrs
   public :: PC_advance
   public :: PC_create_part
   public :: PC_add_part
   public :: PC_periodify
   public :: PC_translate
   public :: PC_get_part_mass
   public :: PC_get_num_sim_part
   public :: PC_reset
   public :: PC_get_num_real_part
   public :: PC_set_accel
   public :: PC_correct_new_accel
   public :: PC_get_num_colls
   public :: PC_get_colls
   public :: PC_get_max_coll_rate
   public :: PC_loop_iopart
   public :: PC_compute_scalar_sum
   public :: PC_compute_vector_sum
   public :: PC_resize_part_list
   public :: PC_share_particles
   public :: sort_pl
   public :: PC_merge_and_split
   public :: PC_clean_up
   public :: histogram_pl
   public :: PC_v_to_en
   public :: PC_speed_to_en
   public :: PC_en_to_vel
   public :: PC_merge_part_rxv
   public :: PC_split_part
   public :: PC_get_coeffs

contains

   !> Initialization routine for the particle module
   subroutine PC_initialize(mass, cross_secs, lookup_table_size, max_en_eV, &
        n_part_max, n_part_lists)
      use m_cross_sec
      use m_units_constants
      type(CS_type), intent(in) :: cross_secs(:)
      integer, intent(in)       :: lookup_table_size
      real(dp), intent(in)      :: mass, max_en_eV
      integer, intent(in)       :: n_part_max, n_part_lists
      integer                   :: ix, n_part_per_list

      if (size(cross_secs) < 1) then
         print *, "No cross sections given, will abort"
         stop
      else if (size(cross_secs) > PC_max_num_colls) then
         print *, "Increase PC_max_num_colls, current buffer sizes"
         print *, "are too small to hold all the collisions"
         stop
      end if

      ! Allocate particle lists
      n_part_per_list = n_part_max/n_part_lists
      allocate(PC_PL(n_part_lists))

      ! For now, all particles get the same cross sections and mass,
      ! but this could be changed easily in the future
      do ix = 1, n_part_lists
         PC_PL(ix)%mass   = mass
         PC_PL(ix)%n_part = 0

         allocate(PC_PL(ix)%parts(n_part_per_list))
         call RNG_set_state(PC_PL(ix)%rng, (/ix, 1337, 31337, 89234/))
         call create_coll_rate_table(PC_PL(ix)%colls, cross_secs, &
              mass, 0.0_dp, max_en_eV, lookup_table_size)
      end do

   end subroutine PC_initialize

   subroutine PC_reset()
      integer :: ix
      do ix = 1, size(PC_PL)
         PC_PL(ix)%n_part = 0
         call LL_clear(PC_PL(ix)%clean_list)
      end do
   end subroutine PC_reset

   subroutine PC_resize_part_list(ix, new_size)
      integer, intent(in)          :: ix, new_size
      type(PC_part_t), allocatable :: parts_copy(:)

      allocate(parts_copy(PC_PL(ix)%n_part))
      deallocate(PC_PL(ix)%parts)
      allocate(PC_PL(ix)%parts(new_size))
      PC_PL(ix)%parts(1:PC_PL(ix)%n_part) = parts_copy
   end subroutine PC_resize_part_list

   ! Share particles between particle lists
   subroutine PC_share_particles(arg_ixs)
      integer, intent(in), optional :: arg_ixs(:)
      integer, allocatable :: ixs(:)
      integer :: i, i_temp(1)
      integer :: n_avg, i_min, i_max, n_min, n_max, n_send

      if (present(arg_ixs)) then
         allocate(ixs(size(arg_ixs)))
         ixs = arg_ixs
      else
         allocate(ixs(size(PC_PL)))
         ixs = (/ (i, i=1,size(PC_PL)) /)
      end if

      n_avg = ceiling(sum(PC_PL(ixs)%n_part) / (1.0_dp * size(ixs)))

      do
         i_temp = maxloc(PC_PL(ixs)%n_part)
         i_max  = i_temp(1)
         n_max  = PC_PL(i_max)%n_part
         i_temp = minloc(PC_PL(ixs)%n_part)
         i_min  = i_temp(1)
         n_min  = PC_PL(i_min)%n_part

         ! Difference it at most size(ixs) - 1, if all lists get one more particle
         ! than the last list
         if (n_max - n_min < size(ixs)) exit

         ! Send particles from i_max to i_min
         n_send = min(n_max - n_avg, n_avg - n_min)
         ! print *, n_avg, n_send, i_min, n_min, i_max, n_max
         PC_PL(i_min)%parts(n_min+1:n_min+n_send) = &
              PC_PL(i_max)%parts(n_max-n_send+1:n_max)

         ! Always at the end of a list, so do not need to clean up later
         PC_PL(i_min)%n_part = PC_PL(i_min)%n_part + n_send
         PC_PL(i_max)%n_part = PC_PL(i_max)%n_part - n_send
      end do

   end subroutine PC_share_particles

   subroutine PC_set_pptrs(pptr_bg_vel_sampler, pptr_coll_callback)
      procedure(if_ovec3), optional      :: pptr_bg_vel_sampler
      procedure(if_iint2_ipart), optional :: pptr_coll_callback
      if (present(pptr_bg_vel_sampler)) PC_pptr_bg_vel_sampler => pptr_bg_vel_sampler
      if (present(pptr_coll_callback)) PC_pptr_coll_callback => pptr_coll_callback
   end subroutine PC_set_pptrs

   !> Loop over all the particles, and for each set the time left 'dt' for this step.
   !! Then they are fed to the move_and_collide() routine, which advances them in time
   subroutine PC_advance(dt)
      use omp_lib
      real(dp), intent(in)      :: dt
      integer                   :: ix

      !$omp parallel do
      do ix = 1, size(PC_PL)
         ! print *, ix, PC_PL(ix)%n_part, omp_get_thread_num()
         call advance_pl(PC_PL(ix), dt)
      end do
      !$omp end parallel do
   end subroutine PC_advance

   subroutine advance_pl(pl, dt)
      type(PC_list_t), intent(inout) :: pl
      real(dp), intent(in)           :: dt
      integer                        :: ll

      pl%parts(1:pl%n_part)%t_left = dt
      ll = 1

      do while (ll <= pl%n_part)
         call move_and_collide(pl, ll)
         ll = ll + 1
      end do

      call PC_clean_up(pl)
   end subroutine advance_pl

   !> Perform a collision for an electron, either elastic, excitation, ionizationCollision,
   !! attachment or null.
   subroutine move_and_collide(pl, ll)
      use m_cross_sec

      type(PC_list_t), intent(inout) :: pl
      integer, intent(in)            :: ll
      integer                        :: cIx, cType
      real(dp)                       :: coll_time, new_vel

      do
         coll_time = sample_coll_time(pl%colls, pl%rng) ! Set the next collision time
         if (coll_time > pl%parts(ll)%t_left) exit

         call advance_particle(pl%parts(ll), coll_time)   ! Set x,v at the collision time
         new_vel = norm2(pl%parts(ll)%v)
         cIx     = get_coll_index(pl%colls, new_vel, pl%rng)

         if (cIx > 0) then
            cType = pl%colls%types(cIx)

            if (associated(PC_pptr_coll_callback)) then
               ! TODO: should be thread safe
               ! TODO: add check here whether the particle is still
               ! in a valid region of the domain
               call PC_pptr_coll_callback(pl%parts(ll), cType, cIx)
            end if

            ! Perform the corresponding collision
            select case (cType)
            case (CS_attach_t)
               call attach_collision(pl, ll)
               go to 100 ! Particle is removed, so exit
            case (CS_elastic_t)
               call elastic_collision(pl, ll, cIx)
            case (CS_excite_t)
               call excite_collision(pl, ll, cIx, new_vel)
            case (CS_ionize_t)
               call ionization_collision(pl, ll, cIx, new_vel)
            end select
         end if
      end do

      ! Update the particle position and velocity to the next timestep
      call advance_particle(pl%parts(ll), pl%parts(ll)%t_left)
100   continue
   end subroutine move_and_collide

   !> Returns a sample from the exponential distribution of the collision times
   ! RNG_uniform() is uniform on [0,1), but log(0) = nan, so we take 1 - RNG_uniform()
   real(dp) function sample_coll_time(colls, rng_state)
      type(PC_coll_t), intent(in) :: colls
      type(RNG_state_t), intent(inout) :: rng_state
      sample_coll_time = -log(1.0D0 - RNG_U01(rng_state)) * colls%inv_max_rate
   end function sample_coll_time

   !> From the list crosssec(:) select the index of the process that will occur,
   !! or set colIndex = 0 if there is a null collision
   ! TODO: improve as in paper
   integer function get_coll_index(colls, velocity, rng_state)
      type(PC_coll_t), intent(in) :: colls
      real(dp), intent(IN) :: velocity
      type(RNG_state_t), intent(inout) :: rng_state
      integer              :: j
      real(dp)             :: rate, rand_rate

      ! Making this an automatic array leads to a big slowdown with OpenMP
      real(dp) :: buffer(PC_max_num_colls)

      get_coll_index = 0
      rand_rate = RNG_U01(rng_state) * colls%max_rate ! Random collision frequency

      ! Fill an array with interpolated rates
      buffer = LT_get_mcol(colls%rate_lt, velocity)

      ! Determine the type of collision
      ! This uses that process j has a probability coll_rate(j) / PC_coll%max_rate
      ! TODO: binary search
      rate = 0.0D0
      do j = 1, colls%num
         rate = rate + buffer(j)
         if (rand_rate <= rate) then
            get_coll_index = j
            exit
         end if
      end do

   end function get_coll_index

   real(dp) function PC_get_max_coll_rate(ix)
      integer :: ix
      PC_get_max_coll_rate = PC_PL(ix)%colls%max_rate
   end function PC_get_max_coll_rate

   !> Perform an elastic collision for particle 'll'
   subroutine elastic_collision(pl, ll, coll_ix)
      type(PC_list_t), intent(inout) :: pl
      integer, intent(IN)  :: ll, coll_ix
      real(dp)             :: bg_vel(3), com_vel(3)

      if (associated(PC_pptr_bg_vel_sampler)) then
         call PC_pptr_bg_vel_sampler(bg_vel)
      else
         bg_vel = 0.0_dp
      end if

      ! Compute center of mass velocity
      com_vel = (pl%colls%special_val(coll_ix) * pl%parts(ll)%v + bg_vel) / &
           (1 + pl%colls%special_val(coll_ix))

      ! Scatter in center of mass coordinates
      pl%parts(ll)%v = pl%parts(ll)%v - com_vel
      call scatter_isotropic(pl%parts(ll), norm2(pl%parts(ll)%v), pl%rng)
      pl%parts(ll)%v = pl%parts(ll)%v + com_vel
   end subroutine elastic_collision

   !> Perform an excitation-collision for particle 'll'
   subroutine excite_collision(pl, ll, coll_ix, velocity)
      use m_units_constants
      type(PC_list_t), intent(inout) :: pl
      integer, intent(IN)  :: coll_ix, ll
      real(dp), intent(IN) :: velocity
      real(dp)             :: new_vel, energy, old_en

      old_en  = PC_speed_to_en(velocity, pl%mass)
      energy  = max(0.0_dp, old_en - pl%colls%special_val(coll_ix))
      new_vel = PC_en_to_vel(energy, pl%mass)

      call scatter_isotropic(pl%parts(ll), new_vel, pl%rng)
   end subroutine excite_collision

   !> Perform an ionizing collision for particle 'll'
   subroutine ionization_collision(pl, ll, coll_ix, velocity)
      use m_units_constants
      type(PC_list_t), intent(inout) :: pl
      integer, intent(in)            :: coll_ix, ll
      real(dp), intent(IN)           :: velocity
      type(PC_part_t)                :: new_part
      real(dp)                       :: energy, old_en, first_vel, second_vel
      real(dp)                       :: en_s1, en_s2

      old_en     = PC_speed_to_en(velocity, pl%mass)
      energy     = max(0.0_dp, old_en - pl%colls%special_val(coll_ix))
      en_s1      = 0.5D0 * energy
      en_s2      = 0.5D0 * energy
      first_vel  = PC_en_to_vel(en_s1, pl%mass)
      second_vel = PC_en_to_vel(en_s2, pl%mass)

      new_part = pl%parts(ll)
      call scatter_isotropic(pl%parts(ll), first_vel, pl%rng)
      call scatter_isotropic(new_part, second_vel, pl%rng)
      call PC_add_part(pl, new_part)
   end subroutine ionization_collision

   !> Perform attachment of electron 'll'
   subroutine attach_collision(pl, ll)
      type(PC_list_t), intent(inout) :: pl
      integer, intent(in)            :: ll
      call PC_remove_part(pl, ll)
   end subroutine attach_collision

   subroutine scatter_isotropic(part, vel_norm, rng_state)
      use m_units_constants
      type(PC_part_t), intent(inout) :: part
      real(dp), intent(in)           :: vel_norm
      type(RNG_state_t), intent(inout) :: rng_state
      real(dp)                       :: sum_sq, tmp_sqrt, rands(2)

      ! Marsaglia method for uniform sampling on sphere
      do
         rands(1) = RNG_U01(rng_state)
         rands(2) = RNG_U01(rng_state)
         rands = rands * 2 - 1
         sum_sq = rands(1)**2 + rands(2)**2
         if (sum_sq <= 1) exit
      end do

      tmp_sqrt = sqrt(1 - sum_sq)
      part%v(1) = 2 * rands(1) * tmp_sqrt
      part%v(2) = 2 * rands(2) * tmp_sqrt
      part%v(3) = 1 - 2 * sum_sq
      part%v = part%v * vel_norm ! Normalization

   end subroutine scatter_isotropic

   !> Advance the particle position and velocity over time tt
   subroutine advance_particle(part, tt)
      type(PC_part_t), intent(inout) :: part
      real(dp), intent(IN) :: tt

      part%x      = part%x + part%v * tt + &
           0.5_dp * part%a * tt**2
      part%v      = part%v + part%a * tt
      part%t_left = part%t_left - tt
   end subroutine advance_particle

   subroutine PC_set_accel(accel_func)
      procedure(if_ipart_oreal3) :: accel_func
      integer                   :: ix, ll
      real(dp)                  :: new_accel(3)

      do ix = 1, size(PC_PL)
         do ll = 1, PC_PL(ix)%n_part
            call accel_func(PC_PL(ix)%parts(ll), new_accel)
            PC_PL(ix)%parts(ll)%a = new_accel
         end do
      end do
   end subroutine PC_set_accel

   !> Correct particle velocities for the previous timestep of 'dt'
   !!
   !! During the timestep x,v have been advanced to:
   !! x(t+1) = x(t) + v(t)*dt + 0.5*a(t)*dt^2,
   !! v(t+1) = v(t) + a(t)*dt
   !! But the velocity at t+1 should be v(t+1) = v(t) + 0.5*(a(t) + a(t+1))*dt,
   !! to have a second order leapfrog scheme, so here we set it to that value.
   subroutine PC_correct_new_accel(dt, accel_func)
      use m_units_constants
      real(dp), intent(IN)      :: dt
      procedure(if_ipart_oreal3) :: accel_func
      integer                   :: ix, ll
      real(dp)                  :: new_accel(3)

      do ix = 1, size(PC_PL)
         do ll = 1, PC_PL(ix)%n_part
            call accel_func(PC_PL(ix)%parts(ll), new_accel)
            PC_PL(ix)%parts(ll)%v = PC_PL(ix)%parts(ll)%v + &
                 0.5_dp * (new_accel - PC_PL(ix)%parts(ll)%a) * dt
            PC_PL(ix)%parts(ll)%a = new_accel
         end do
      end do
   end subroutine PC_correct_new_accel

   subroutine PC_clean_up(pl)
      type(PC_list_t), intent(inout) :: pl
      integer :: ix_end, ix_clean
      logical :: success

      do
         ! Get an index that has to be cleaned from the list
         call LL_pop(pl%clean_list, ix_clean, success)
         if (.not. success) exit

         ! Find the last "alive" particle in the list
         do ix_end = pl%n_part, 1, -1
            if (pl%parts(ix_end)%weight /= PC_dead_weight) then
               exit
            else
               pl%n_part = pl%n_part - 1
            end if
         end do

         ! Fill in empty spot ix_clean, if it lies before n_part
         if (ix_clean < pl%n_part) then
            pl%parts(ix_clean) = pl%parts(pl%n_part)
            pl%n_part = pl%n_part - 1
         end if
      end do
   end subroutine PC_clean_up

   subroutine PC_add_part(pl, part)
      type(PC_list_t), intent(inout) :: pl
      type(PC_part_t), intent(in)    :: part
      integer                        :: ix
      ix                   = get_ix_new_particle(pl)
      pl%parts(ix) = part
   end subroutine PC_add_part

   function get_ix_new_particle(pl) result(ix)
      type(PC_list_t), intent(inout) :: pl
      integer                        :: ix
      pl%n_part = pl%n_part + 1
      ix            = pl%n_part
   end function get_ix_new_particle

   subroutine PC_create_part(pos, vel, accel, weight, t_left, list_ix)
      real(dp), intent(IN) :: pos(3), vel(3), accel(3), weight, t_left
      integer, optional    :: list_ix
      integer              :: ix, n_part

      if (present(list_ix)) then
         ix = list_ix
      else
         ix = 1
      end if

      n_part                             = PC_PL(ix)%n_part
      PC_PL(ix)%parts(n_part + 1)%weight = weight
      PC_PL(ix)%parts(n_part + 1)%t_left = t_left
      PC_PL(ix)%parts(n_part + 1)%x      = pos
      PC_PL(ix)%parts(n_part + 1)%v      = vel
      PC_PL(ix)%parts(n_part + 1)%a      = accel
      PC_PL(ix)%n_part                   = n_part + 1
   end subroutine PC_create_part

   !> Mark particle for removal
   subroutine PC_remove_part(pl, ix_to_remove)
      type(PC_list_t), intent(inout) :: pl
      integer, intent(in)            :: ix_to_remove

      call LL_add(pl%clean_list, ix_to_remove)
      pl%parts(ix_to_remove)%weight = PC_dead_weight
   end subroutine PC_remove_part

   subroutine PC_periodify(is_periodic, lengths)
      logical, intent(in) :: is_periodic(3)
      real(dp), intent(in) :: lengths(3)
      integer :: ix, i_dim, n_part

      do ix = 1, size(PC_PL)
         n_part = PC_PL(ix)%n_part
         do i_dim = 1, 3
            if (is_periodic(i_dim)) then
               PC_PL(ix)%parts(1:n_part)%x(i_dim) = &
                    modulo(PC_PL(ix)%parts(1:n_part)%x(i_dim), lengths(i_dim))
            end if
         end do
      end do
   end subroutine PC_periodify

   subroutine PC_translate(delta_x)
      real(dp), intent(in) :: delta_x(3)
      integer              :: ix, ll

      do ix = 1, size(PC_PL)
         do ll = 1, PC_PL(ix)%n_part
            PC_PL(ix)%parts(ll)%x = PC_PL(ix)%parts(ll)%x + delta_x
         end do
      end do
   end subroutine PC_translate

   real(dp) function PC_get_part_mass(ix)
      integer, intent(in) :: ix
      PC_get_part_mass = PC_PL(ix)%mass
   end function PC_get_part_mass

   !> Return the number of real particles
   real(dp) function PC_get_num_real_part()
      integer :: ix
      PC_get_num_real_part = 0.0_dp
      do ix = 1, size(PC_PL)
         PC_get_num_real_part = PC_get_num_real_part + &
              sum(PC_PL(ix)%parts(1:PC_PL(ix)%n_part)%weight)
      end do
   end function PC_get_num_real_part

   !> Return the number of simulation particles
   integer function PC_get_num_sim_part()
      PC_get_num_sim_part = sum(PC_PL(:)%n_part)
   end function PC_get_num_sim_part

   !> Loop over all the particles and call pptr for each of them
   subroutine PC_loop_iopart(pptr)
      procedure(if_iopart) :: pptr
      integer              :: ix, ll
      do ix = 1, size(PC_PL)
         do ll = 1, PC_PL(ix)%n_part
            call pptr(PC_PL(ix)%parts(ll))
         end do
      end do
   end subroutine PC_loop_iopart

   subroutine PC_compute_vector_sum(pptr, my_sum)
      procedure(if_ipart_oreals) :: pptr
      real(dp), intent(out)      :: my_sum(:)

      integer                    :: ix, ll
      real(dp)                   :: temp(size(my_sum))

      my_sum = 0.0_dp
      do ix = 1, size(PC_PL)
         do ll = 1, PC_PL(ix)%n_part
            call pptr(PC_PL(ix)%parts(ll), temp)
            my_sum = my_sum + temp
         end do
      end do
   end subroutine PC_compute_vector_sum

   subroutine PC_compute_scalar_sum(pptr, my_sum)
      procedure(if_ipart_oreal) :: pptr
      real(dp), intent(out)     :: my_sum

      integer                   :: ix, ll
      real(dp)                  :: temp

      my_sum = 0.0_dp
      do ix = 1, size(PC_PL)
         do ll = 1, PC_PL(ix)%n_part
            call pptr(PC_PL(ix)%parts(ll), temp)
            my_sum = my_sum + temp
         end do
      end do
   end subroutine PC_compute_scalar_sum

   real(dp) function PC_v_to_en(v, mass)
      real(dp), intent(in) :: v(3), mass
      PC_v_to_en = 0.5_dp * mass * sum(v**2)
   end function PC_v_to_en

   real(dp) elemental function PC_speed_to_en(vel, mass)
      real(dp), intent(in) :: vel, mass
      PC_speed_to_en = 0.5_dp * mass * vel**2
   end function PC_speed_to_en

   real(dp) elemental function PC_en_to_vel(en, mass)
      real(dp), intent(in) :: en, mass
      PC_en_to_vel = sqrt(2 * en / mass)
   end function PC_en_to_vel

   !> Create a lookup table with cross sections for a number of energies
   subroutine create_coll_rate_table(colls, cross_secs, mass, &
        min_eV, max_eV, table_size)
      use m_units_constants
      use m_cross_sec
      use m_lookup_table
      type(PC_coll_t), intent(out) :: colls
      type(CS_type), intent(in) :: cross_secs(:)
      integer, intent(IN)       :: table_size
      real(dp), intent(in)      :: mass, min_eV, max_eV

      real(dp)                  :: vel_list(table_size), rate_list(table_size)
      real(dp)                  :: sum_rate_list(table_size)
      integer                   :: ix, i_coll, i_row
      integer, allocatable      :: collIxs(:)
      real(dp)                  :: min_vel, max_vel, en_eV, temp

      min_vel = PC_en_to_vel(min_eV * UC_elec_volt, mass)
      max_vel = PC_en_to_vel(max_eV * UC_elec_volt, mass)

      colls%num = size(cross_secs)
      allocate( collIxs(colls%num) )
      allocate( colls%types(colls%num) )
      allocate( colls%special_val(colls%num) )

      ! Set a range of velocities
      do ix = 1, table_size
         temp = real(ix-1, dp) / (table_size-1)
         vel_list(ix) = min_vel + temp * (max_vel - min_vel)
      end do

      sum_rate_list = 0.0_dp

      ! Create collision rate table
      colls%rate_lt = LT_create(min_vel, max_vel, table_size, colls%num)

      do i_coll = 1, colls%num
         colls%types(i_coll) = cross_secs(i_coll)%col_type

         select case (colls%types(i_coll))
         case (CS_ionize_t, CS_excite_t)
            ! Ionizations and excitations use a threshold energy which is lost
            colls%special_val(i_coll) = cross_secs(i_coll)%spec_value * UC_elec_volt
         case (CS_elastic_t)
            ! Mass ration electron : neutral
            colls%special_val(i_coll) = cross_secs(i_coll)%spec_value
         case DEFAULT
            colls%special_val(i_coll) = 0.0_dp
         end select

         ! Linear interpolate cross sections by energy
         do i_row = 1, table_size
            en_eV = PC_speed_to_en(vel_list(i_row), mass) / UC_elec_volt
            call LT_lin_interp_list(cross_secs(i_coll)%en_cs(1, :), &
                 cross_secs(i_coll)%en_cs(2, :), en_eV, rate_list(i_row))
            rate_list(i_row) = rate_list(i_row) * vel_list(i_row)
         end do

         call LT_set_col(colls%rate_lt, i_coll, vel_list, rate_list)
         sum_rate_list = sum_rate_list + rate_list
      end do

      colls%max_rate = maxval(sum_rate_list)
      colls%inv_max_rate = 1 / colls%max_rate

   end subroutine create_coll_rate_table

   subroutine sort_pl(pl, sort_func)
      use m_mrgrnk
      type(PC_list_t), intent(inout) :: pl
      procedure(if_freal_ipart)      :: sort_func

      integer                        :: ix, n_part
      integer, allocatable           :: part_ixs(:)
      real(dp), allocatable          :: part_values(:)
      type(PC_part_t), allocatable   :: part_copies(:)

      n_part = pl%n_part
      allocate(part_ixs(n_part))
      allocate(part_values(n_part))
      allocate(part_copies(n_part))

      do ix = 1, n_part
         part_copies(ix) = pl%parts(ix)
         part_values(ix) = sort_func(part_copies(ix))
      end do

      call mrgrnk(part_values, part_ixs)

      do ix = 1, n_part
         pl%parts(ix) = part_copies(part_ixs(ix))
      end do
   end subroutine sort_pl

   subroutine PC_get_histogram(hist_func, filter_func, &
        filter_args, x_values, y_values)
      procedure(if_freal_ipart) :: hist_func
      real(dp), intent(in)      :: x_values(:)
      real(dp), intent(out)     :: y_values(:)
      procedure(if_filter_func) :: filter_func
      real(dp), intent(in)      :: filter_args(:)

      integer                   :: ix
      real(dp)                  :: temp(size(y_values))

      y_values = 0.0_dp
      do ix = 1, size(PC_PL)
         call histogram_pl(PC_PL(ix), hist_func, filter_func, &
              filter_args, x_values, temp)
         y_values = y_values + temp
      end do
   end subroutine PC_get_histogram

   subroutine histogram_pl(pl, hist_func, filter_func, &
        filter_args, x_values, y_values)
      use m_mrgrnk
      type(PC_list_t), intent(in) :: pl
      procedure(if_freal_ipart)   :: hist_func
      real(dp), intent(in)        :: x_values(:)
      real(dp), intent(out)       :: y_values(:)
      procedure(if_filter_func)   :: filter_func
      real(dp), intent(in)        :: filter_args(:)

      integer                     :: ix, p_ix, o_ix, n_used, num_bins, n_part
      integer, allocatable        :: part_ixs(:)
      real(dp)                    :: boundary_value
      real(dp), allocatable       :: part_values(:)
      logical, allocatable        :: part_mask(:)

      n_part = pl%n_part
      allocate(part_mask(n_part))
      do ix = 1, n_part
         part_mask(ix) = filter_func(pl%parts(ix), filter_args)
      end do

      n_used = count(part_mask)

      allocate(part_values(n_used))
      allocate(part_ixs(n_used))

      p_ix = 0
      do ix = 1, pl%n_part
         if (part_mask(ix)) then
            p_ix = p_ix + 1
            part_values(p_ix) = hist_func(pl%parts(p_ix))
         end if
      end do

      call mrgrnk(part_values, part_ixs)

      num_bins = size(x_values)
      p_ix     = 1
      y_values = 0

      outer: do ix = 1, num_bins-1
         boundary_value = 0.5_dp * (x_values(ix) + x_values(ix+1))
         do
            if (p_ix == n_used + 1) exit outer
            o_ix = part_ixs(p_ix) ! Index in the 'old' list
            if (part_values(o_ix) > boundary_value) exit

            y_values(ix) = y_values(ix) + pl%parts(o_ix)%weight
            p_ix         = p_ix + 1
         end do
      end do outer

      ! Fill last bin
      y_values(num_bins) = sum(pl%parts(part_ixs(p_ix:n_used))%weight)
   end subroutine histogram_pl

   subroutine PC_merge_and_split(coord_weights, max_distance, weight_func, &
        pptr_merge, pptr_split)
      real(dp), intent(in)           :: coord_weights(6), max_distance
      procedure(if_freal_ipart)      :: weight_func
      procedure(io_2part_rng)            :: pptr_merge, pptr_split

      integer :: ix
      do ix = 1, size(PC_PL)
         call merge_and_split_pl(PC_PL(ix), coord_weights, max_distance, &
              weight_func, pptr_merge, pptr_split)
      end do
   end subroutine PC_merge_and_split

   ! Routine to merge and split particles. Input arguments are the coordinate weights, used
   ! to determine the 'distance' between particles. The first three elements of the array are
   ! the weights of the xyz position coordinates, the next three the weights of the xyz
   ! velocity coordinates. Max_distance is the maxium Euclidean distance between particles
   ! to be merged. The weight_func returns the desired weight for a particle, whereas the
   ! pptr_merge and pptr_split procedures merge and split particles.
   subroutine merge_and_split_pl(pl, coord_weights, max_distance, weight_func, &
        pptr_merge, pptr_split)
      use m_mrgrnk
      use kdtree2_module
      type(PC_list_t), intent(inout) :: pl
      real(dp), intent(in)           :: coord_weights(6), max_distance
      procedure(if_freal_ipart)      :: weight_func
      procedure(io_2part_rng)            :: pptr_merge, pptr_split

      integer, parameter           :: num_neighbors = 1
      real(dp), parameter          :: large_ratio   = 1.5_dp, small_ratio = 1 / large_ratio
      real(dp)                     :: distance
      type(kdtree2), pointer       :: kd_tree
      type(kdtree2_result)         :: kd_results(num_neighbors)

      integer                      :: num_part, num_merge, num_split
      integer                      :: p_min, p_max, n_too_far
      integer                      :: o_ix, o_nn_ix, new_ix
      integer                      :: ix, neighbor_ix, cntr, num_coords
      logical, allocatable         :: already_merged(:)
      integer, allocatable         :: part_ixs(:)
      real(dp), allocatable        :: coord_data(:, :), weight_ratios(:)
      type(PC_part_t), allocatable :: part_copy(:)

      p_min    = 1
      p_max    = pl%n_part
      num_part = p_max - p_min + 1

      allocate(weight_ratios(num_part))
      allocate(part_ixs(num_part))
      allocate(part_copy(num_part))

      do ix = 1, num_part
         weight_ratios(ix) = pl%parts(p_min+ix-1)%weight / &
              weight_func(pl%parts(p_min+ix-1))
      end do

      num_merge      = count(weight_ratios <= small_ratio)
      num_coords     = count(coord_weights /= 0.0_dp)
      allocate(coord_data(num_coords, num_merge))
      allocate(already_merged(num_merge))
      already_merged = .false.
      n_too_far      = 0

      ! Sort particles by their relative weight and store them in part_copy
      ! so that particles to be merged are at the beginning of the list
      call mrgrnk(weight_ratios, part_ixs)
      part_copy = pl%parts(p_min + part_ixs - 1)

      ! Only create a k-d tree if there are enough particles to be merged
      if (num_merge > num_coords) then
         ! Store the coordinates of the particles to be merged in coord_data
         cntr = 0
         do ix = 1, 6
            if (coord_weights(ix) /= 0.0_dp) then
               cntr = cntr + 1
               if (ix <= 3) then ! Spatial coordinates
                  coord_data(cntr, :) = part_copy(1:num_merge)%x(ix) * coord_weights(ix)
               else              ! Velocity coordinates
                  coord_data(cntr, :) = part_copy(1:num_merge)%v(ix-3) * coord_weights(ix)
               end if
            end if
         end do

         ! Create k-d tree
         kd_tree => kdtree2_create(coord_data)

         ! Merge particles
         do ix = 1, num_merge
            if (already_merged(ix)) cycle

            call kdtree2_n_nearest_around_point(kd_tree, idxin=ix, &
                 nn=num_neighbors, correltime=1, results=kd_results)
            neighbor_ix = kd_results(1)%idx

            if (already_merged(neighbor_ix)) cycle

            distance = norm2(coord_data(:, ix) - coord_data(:, neighbor_ix))
            if (distance > max_distance) cycle

            ! Get indices in the original particle list
            o_ix = part_ixs(ix)
            o_nn_ix = part_ixs(neighbor_ix)

            ! Merge, then remove neighbor
            call pptr_merge(pl%parts(o_ix), pl%parts(o_nn_ix), pl%rng)
            call PC_remove_part(pl, o_nn_ix)
            already_merged((/ix, neighbor_ix/)) = .true.
         end do
      end if

      ! Split particles. These are at the end of part_copy
      num_split = count(weight_ratios >= large_ratio)

      do ix = num_part - num_split + 1, num_part
         ! Change part_copy(ix), then add an extra particle at then end of pl
         o_ix = part_ixs(ix)
         new_ix = get_ix_new_particle(pl)
         call pptr_split(pl%parts(o_ix), pl%parts(new_ix), pl%rng)
      end do

      call PC_clean_up(pl)
   end subroutine merge_and_split_pl

   ! Merge two particles into part_a, should remove part_b afterwards
   subroutine PC_merge_part_rxv(part_a, part_b, rng_state)
      type(PC_part_t), intent(inout) :: part_a, part_b
      type(RNG_state_t), intent(inout) :: rng_state
      real(dp)                       :: rr

      rr = RNG_U01(rng_state)

      if (rr > part_a%weight / (part_a%weight + part_b%weight)) then
         ! Keep particle b's position and velocity
         part_a%v      = part_b%v
         part_a%x      = part_b%x
      end if
      part_a%weight = part_a%weight + part_b%weight
   end subroutine PC_merge_part_rxv

   subroutine PC_split_part(part_a, part_b, rng_state)
      type(PC_part_t), intent(inout) :: part_a, part_b
      type(RNG_state_t), intent(inout) :: rng_state
      part_a%weight = part_a%weight * 0.5_dp
      part_b        = part_a
   end subroutine PC_split_part

   integer function PC_get_num_colls(ix)
      integer, intent(in) :: ix
      PC_get_num_colls = PC_PL(ix)%colls%num
   end function PC_get_num_colls

   subroutine PC_get_colls(ix, out_colls)
      integer, intent(in) :: ix
      type(PC_coll_t), intent(out) :: out_colls
      out_colls = PC_PL(ix)%colls
   end subroutine PC_get_colls

   subroutine PC_get_coeffs(coeff_data, coeff_names, n_coeffs)
      use m_cross_sec
      use m_lookup_table
      real(dp), intent(out), allocatable          :: coeff_data(:,:)
      character(len=*), intent(out), allocatable :: coeff_names(:)
      integer, intent(out)                        :: n_coeffs
      type(PC_coll_t)                             :: coll_data
      integer                                     :: nn, n_rows

      call PC_get_colls(1, coll_data)
      n_coeffs = coll_data%num + 2
      n_rows = LT_get_num_rows(coll_data%rate_lt)
      allocate(coeff_data(n_coeffs, n_rows))
      allocate(coeff_names(n_coeffs))

      call LT_get_data(coll_data%rate_lt, coeff_data(1, :), coeff_data(3:,:))
      coeff_names(1) = "velocity (m/s)"
      coeff_names(2) = "sum coll_rate (1/s)"
      do nn = 1, coll_data%num
         select case (coll_data%types(nn))
         case (CS_ionize_t)
            coeff_names(2+nn) = "ionization"
         case (CS_attach_t)
            coeff_names(2+nn) = "attachment"
         case (CS_elastic_t)
            coeff_names(2+nn) = "elastic"
         case (CS_excite_t)
            coeff_names(2+nn) = "excitation"
         case default
            coeff_names(2+nn) = "unknown"
         end select
      end do
   end subroutine PC_get_coeffs

end module m_particle_core
