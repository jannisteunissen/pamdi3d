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

   implicit none
   private

   integer, parameter :: dp = kind(0.0d0)

   !> The number of active particles in the simulation
   integer     :: PC_num_part

   !> The electron type
   type PC_part_t
      real(dp) :: T_left   ! The time until the next timestep
      real(dp) :: x(3)     ! The position
      real(dp) :: v(3)     ! The velocity
      real(dp) :: a(3)     ! The electric field's acceleration
      real(dp) :: weight   ! The weight factor of the (super)particle
      logical  :: live     ! If .false. the particle will not be used (anymore)
   end type PC_part_t

   !> The mass of the particles in the simulation
   real(dp) :: PC_particle_mass

   type PC_coll_t
      integer               :: num                    ! The number of different collisions
      real(dp)              :: max_rate, inv_max_rate ! Maximum collision rate and inverse
      integer, allocatable  :: types(:)               ! The types of the collisions that can occur
      real(dp), allocatable :: special_val(:)         ! Special value used by this collision
      real(dp), allocatable :: counts(:)              ! To collect statistics during runtime
      type(LT_col_t)        :: sum_rate_lt
      type(LT_mcol_t)       :: rate_lt                ! Lookup table with collision rates
   end type PC_coll_t

   type(PC_coll_t) :: PC_coll

   !> The list that contains all the particles in the simulation
   type(PC_part_t), allocatable :: PC_particles(:)

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
   public :: PC_set_part
   public :: PC_get_part
   public :: PC_get_part_mass
   public :: PC_get_num_sim_part
   public :: PC_reset
   public :: PC_get_num_real_part
   public :: PC_set_accel
   public :: PC_correct_new_accel
   public :: PC_get_num_colls
   public :: PC_get_colls
   public :: PC_get_coll_count
   public :: PC_get_coll_rates
   public :: PC_get_max_coll_rate
   public :: PC_reset_coll_count
   public :: PC_loop_ipart
   public :: PC_loop_iopart
   public :: PC_compute_sum
   public :: PC_compute_vector_sum
   public :: PC_sort_particles
   public :: PC_merge_and_split
   public :: PC_clean_up_dead
   public :: PC_get_histogram
   public :: PC_vel_to_en
   public :: PC_en_to_vel
   public :: PC_merge_part_rxv
   public :: PC_split_part
   public :: PC_get_coeffs

contains

   !> Initialization routine for the particle module
   subroutine PC_initialize(mass, cross_secs, lookup_table_size, max_en_eV)
      use m_cross_sec
      use m_units_constants
      type(CS_type), intent(in) :: cross_secs(:)
      integer, intent(in)       :: lookup_table_size
      real(dp), intent(in)      :: mass, max_en_eV

      PC_particle_mass = mass
      allocate(PC_particles(1)) ! Storage for first particle
      PC_num_part      = 0

      if (size(cross_secs) > 0) then
         ! Create the lookup table for the collision rates, minimum velocity is 0.0
         call create_coll_rate_table(cross_secs, 0.0_dp, max_en_eV, lookup_table_size)
         allocate(PC_coll%counts(PC_coll%num))
         PC_coll%counts(:) = 0
      else
         PC_coll%num = 0
      end if

   end subroutine PC_initialize

   subroutine PC_reset()
      PC_num_part = 0
      PC_coll%counts(:) = 0
   end subroutine PC_reset

   subroutine PC_set_pptrs(pptr_bg_vel_sampler, pptr_coll_callback)
      procedure(if_ovec3), optional      :: pptr_bg_vel_sampler
      procedure(if_iint2_ipart), optional :: pptr_coll_callback
      if (present(pptr_bg_vel_sampler)) PC_pptr_bg_vel_sampler => pptr_bg_vel_sampler
      if (present(pptr_coll_callback)) PC_pptr_coll_callback => pptr_coll_callback
   end subroutine PC_set_pptrs

   !> Loop over all the particles, and for each set the time left 'dt' for this step.
   !! Then they are fed to the move_and_collide() routine, which advances them in time
   subroutine PC_advance(dt)
      real(dp), intent(IN)      :: dt
      integer                   :: ll

      call ensure_storage_particles(int(PC_num_part * 2))

      if (PC_coll%num > 0) then ! There are collisions
         PC_particles(1:PC_num_part)%t_left = dt
         ll = 1
         do while (ll <= PC_num_part)
            call move_and_collide(PC_particles(ll))
            ll = ll + 1
         end do
         call PC_clean_up_dead()
      else                      ! There are no collisions
         do ll = 1, PC_num_part
            call advance_particle(PC_particles(ll), dt)
         end do
      end if

   end subroutine PC_advance

   !> Perform a collision for an electron, either elastic, excitation, ionizationCollision,
   !! attachment or null.
   subroutine move_and_collide(part)
      use m_cross_sec

      type(PC_part_t), intent(inout) :: part
      integer                        :: cIx, cType
      real(dp)                       :: coll_time, new_vel

      do
         coll_time = sample_coll_time() ! Set the next collision time
         if (coll_time > part%t_left .or. .not. part%live) exit

         call advance_particle(part, coll_time)   ! Set x,v at the collision time
         new_vel = norm2(part%v)
         cIx     = get_coll_index(new_vel)

         if (cIx > 0) then
            cType               = PC_coll%types(cIx)
            PC_coll%counts(cIx) = PC_coll%counts(cIx) + part%weight

            if (associated(PC_pptr_coll_callback)) then
               call PC_pptr_coll_callback(part, cType, cIx)
            end if

            ! Perform the corresponding collision
            select case (cType)
            case (CS_attach_t)
               call attach_collision(part)
            case (CS_elastic_t)
               call elastic_collision(part, cIx)
            case (CS_excite_t)
               call excite_collision(part, cIx, new_vel)
            case (CS_ionize_t)
               call ionization_collision(part, cIx, new_vel)
            end select
         end if
      end do

      ! Update the particle position and velocity to the next timestep
      call advance_particle(part, part%t_left)

   end subroutine move_and_collide

   !> Returns a sample from the exponential distribution of the collision times
   ! RNG_uniform() is uniform on [0,1), but log(0) = nan, so we take 1 - RNG_uniform()
   real(dp) function sample_coll_time()
      use m_random
      real(dp) :: rand_val
      call random_number(rand_val)
      sample_coll_time = -log(1.0D0 - rand_val) * PC_coll%inv_max_rate
   end function sample_coll_time

   !> From the list crosssec(:) select the index of the process that will occur,
   !! or set colIndex = 0 if there is a null collision
   ! TODO: improve as in paper
   integer function get_coll_index(velocity)
      use m_random
      real(dp), intent(IN) :: velocity
      integer              :: j
      real(dp)             :: rate, rand_rate, sum_rate
      real(dp)             :: rate_buffer(PC_coll%num)

      get_coll_index = 0
      sum_rate       = LT_get_col(PC_coll%sum_rate_lt, velocity)
      call random_number(rand_rate)
      rand_rate      = rand_rate * PC_coll%max_rate ! Random collision frequency

      if (rand_rate > sum_rate) return ! If we have a null collision we are done

      ! Otherwise, fill an array with interpolated rates
      call PC_get_coll_rates(velocity, rate_buffer)

      ! Determine the type of collision
      ! This uses that process j has a probability coll_rate(j) / PC_coll%max_rate
      rate = 0.0D0
      do j = 1, PC_coll%num
         rate = rate + rate_buffer(j)
         if (rand_rate <= rate) then
            get_coll_index = j
            exit
         end if
      end do

   end function get_coll_index

   subroutine PC_get_coll_rates(velocity, coll_rates)
      real(dp), intent(in) :: velocity
      real(dp), intent(out) :: coll_rates(:)
      coll_rates = LT_get_mcol(PC_coll%rate_lt, velocity)
   end subroutine PC_get_coll_rates

   real(dp) function PC_get_max_coll_rate()
      PC_get_max_coll_rate = PC_coll%max_rate
   end function PC_get_max_coll_rate

   !> Perform an elastic collision for particle 'll'
   subroutine elastic_collision(part, coll_ix)
      type(PC_part_t), intent(inout) :: part
      integer, intent(IN)  :: coll_ix
      real(dp)             :: bg_vel(3), com_vel(3)

      if (associated(PC_pptr_bg_vel_sampler)) then
         call PC_pptr_bg_vel_sampler(bg_vel)
      else
         bg_vel = 0.0_dp
      end if

      ! Compute center of mass velocity
      com_vel = (PC_coll%special_val(coll_ix) * part%v + bg_vel) / &
           (1 + PC_coll%special_val(coll_ix))

      ! Scatter in center of mass coordinates
      ! part%v = part%v - com_vel
      call scatter_isotropic(part, norm2(part%v))
      ! part%v = part%v + com_vel
   end subroutine elastic_collision

   !> Perform an excitation-collision for particle 'll'
   subroutine excite_collision(part, coll_ix, velocity)
      use m_units_constants
      type(PC_part_t), intent(inout) :: part
      integer, intent(IN)  :: coll_ix
      real(dp), intent(IN) :: velocity
      real(dp)             :: new_vel, energy, old_en

      old_en  = PC_vel_to_en(velocity)
      energy  = max(0.0_dp, old_en - PC_coll%special_val(coll_ix))
      new_vel = PC_en_to_vel(energy)

      call scatter_isotropic(part, new_vel)
   end subroutine excite_collision

   !> Perform an ionizing collision for particle 'll'
   subroutine ionization_collision(part, coll_ix, velocity)
      use m_units_constants
      type(PC_part_t), intent(inout) :: part
      integer, intent(in)            :: coll_ix
      real(dp), intent(IN)           :: velocity
      type(PC_part_t)                :: new_part
      real(dp)                       :: energy, old_en, first_vel, second_vel
      real(dp)                       :: en_s1, en_s2

      old_en     = PC_vel_to_en(velocity)
      energy     = max(0.0_dp, old_en - PC_coll%special_val(coll_ix))
      en_s1      = 0.5D0 * energy
      en_s2      = 0.5D0 * energy
      first_vel  = PC_en_to_vel(en_s1)
      second_vel = PC_en_to_vel(en_s2)

      new_part = part
      call scatter_isotropic(part, first_vel)
      call scatter_isotropic(new_part, second_vel)
      call PC_add_part(new_part)
   end subroutine ionization_collision

   !> Perform attachment of electron 'll'
   subroutine attach_collision(part)
      type(PC_part_t), intent(inout) :: part
      call PC_remove_part(part)
   end subroutine attach_collision

   subroutine scatter_isotropic(part, vel_norm)
      use m_random
      use m_units_constants
      type(PC_part_t), intent(inout) :: part
      real(dp), intent(in)           :: vel_norm
      real(dp)                       :: sum_sq, tmp_sqrt, rands(2)

      ! Marsaglia method for uniform sampling on sphere
      do
         call random_number(rands)
         rands = rands * 2 - 1
         sum_sq = rands(1)**2 + rands(2)**2
         if (sum_sq > 1) cycle

         tmp_sqrt = sqrt(1 - sum_sq)
         part%v(1) = 2 * rands(1) * tmp_sqrt
         part%v(2) = 2 * rands(2) * tmp_sqrt
         part%v(3) = 1 - 2 * sum_sq
         part%v = part%v * vel_norm ! Normalization
         exit
      end do

   end subroutine scatter_isotropic

   !> Rotate vel by chi and psi (spherical angles)
   subroutine scatter_velocity(vel, chi, psi)
      use m_units_constants
      real(dp), intent(IN)    :: chi, psi
      real(dp), intent(inout) :: vel(3)

      real(dp) :: vel_norm, theta, phi
      real(dp) :: costheta, cosphi, sintheta, sinphi, coschi
      real(dp) :: cospsi, sinchi, sinpsi

      ! Determine incoming angles
      call UC_xyz_to_spherical(vel, vel_norm, theta, phi)

      costheta = cos(theta)
      sintheta = sin(theta)
      cosphi   = cos(phi)
      sinphi   = sin(phi)
      coschi   = cos(chi)
      sinchi   = sin(chi)
      cospsi   = cos(psi)
      sinpsi   = sin(psi)

      ! The matrix on the right hand side is obtained by doing the following elemental rotations:
      ! Rotate around z-axis over -phi, rotate around y-axis over -theta,
      ! now vel lies on the z-axis, so only take z-component into next step
      ! rotate around y-axis over xhi, rotate around z-axis over psi
      ! Do the inverse of the first two rotations (y-axis over theta, z-axis over phi).

      vel(1) = ( sintheta*cosphi*coschi + sinchi*(costheta*cosphi*cospsi-sinphi*sinpsi) )
      vel(2) = ( sintheta*sinphi*coschi + sinchi*(costheta*sinphi*cospsi+cosphi*sinpsi) )
      vel(3) = ( costheta*coschi - sintheta*sinchi*cospsi )
      vel    = vel * vel_norm
   end subroutine scatter_velocity

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
      integer                   :: ll
      real(dp)                  :: new_accel(3)

      do ll = 1, PC_num_part
         call accel_func(PC_particles(ll), new_accel)
         PC_particles(ll)%a = new_accel
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
      integer                   :: ll
      real(dp)                  :: new_accel(3)

      do ll = 1, PC_num_part
         call accel_func(PC_particles(ll), new_accel)
         PC_particles(ll)%v = PC_particles(ll)%v + 0.5_dp * (new_accel - PC_particles(ll)%a) * dt
         PC_particles(ll)%a = new_accel
      end do
   end subroutine PC_correct_new_accel

   !> Remove all particles for which live == .FALSE. from the list, which can happen
   !! due moving outside the computational domain or attachment.
   subroutine PC_clean_up_dead()
      integer :: ll, ix_end
      ix_end = PC_num_part

      do ll = 1, PC_num_part
         if (.not. PC_particles(ll)%live) then
            ! Find the first alive particle from the end of the list, and place it at index ll
            do while (.not. PC_particles(ix_end)%live .and. ix_end > ll)
               ix_end = ix_end - 1
            end do

            ! If there is no alive particle available in the range [ll+1 : PC_num_part] we have
            ! ll >= ix_end and we exit the routine. Else we put the alive particle at index ll
            if (ll >= ix_end) then
               PC_num_part = ll - 1
               exit
            else
               PC_particles(ll)          = PC_particles(ix_end)
               PC_particles(ix_end)%live = .false.
               ix_end                    = ix_end - 1
            end if
         end if
      end do
   end subroutine PC_clean_up_dead

   subroutine PC_add_part(part)
      type(PC_part_t), intent(in) :: part

      PC_num_part = PC_num_part + 1
      call ensure_storage_particles(PC_num_part)
      PC_particles(PC_num_part) = part
   end subroutine PC_add_part

   subroutine PC_create_part(pos, vel, accel, weight, t_left)
      real(dp), intent(IN) :: pos(3), vel(3), accel(3), weight, t_left

      PC_num_part = PC_num_part + 1
      call ensure_storage_particles(PC_num_part)

      PC_particles(PC_num_part)%weight = weight
      PC_particles(PC_num_part)%live   = .true.
      PC_particles(PC_num_part)%t_left = t_left
      PC_particles(PC_num_part)%x      = pos
      PC_particles(PC_num_part)%v      = vel
      PC_particles(PC_num_part)%a      = accel
   end subroutine PC_create_part

   !> Mark particle 'll' as inactive, so that it will be removed from the particle list
   subroutine PC_remove_part(part)
      type(PC_part_t), intent(inout) :: part
      part%live   = .false.
      part%weight = 0
   end subroutine PC_remove_part

   subroutine PC_periodify(is_periodic, lengths)
      logical, intent(in) :: is_periodic(3)
      real(dp), intent(in) :: lengths(3)
      integer :: ix

      do ix = 1, 3
         if (is_periodic(ix)) then
            PC_particles(1:PC_num_part)%x(ix) = &
                 modulo(PC_particles(1:PC_num_part)%x(ix), lengths(ix))
         end if
      end do
   end subroutine PC_periodify

   subroutine PC_translate(delta_x)
      real(dp), intent(in) :: delta_x(3)
      integer :: ix
      do ix = 1, PC_num_part
         PC_particles(ix)%x = PC_particles(ix)%x + delta_x
      end do
   end subroutine PC_translate

   subroutine PC_get_part(ll, my_part)
      integer, intent(in)          :: ll
      type(PC_part_t), intent(out) :: my_part
      my_part = PC_particles(ll)
   end subroutine PC_get_part

   subroutine PC_set_part(ll, my_part)
      type(PC_part_t), intent(in) :: my_part
      integer, intent(in)         :: ll
      PC_particles(ll) = my_part
   end subroutine PC_set_part

   real(dp) function PC_get_part_mass()
      PC_get_part_mass = PC_particle_mass
   end function PC_get_part_mass

   !> Return the number of real particles
   real(dp) function PC_get_num_real_part()
      PC_get_num_real_part = sum(PC_particles(1:PC_num_part)%weight)
   end function PC_get_num_real_part

   !> Return the number of simulation particles
   integer function PC_get_num_sim_part()
      PC_get_num_sim_part = PC_num_part
   end function PC_get_num_sim_part

   !> Loop over all the particles and call pptr_loop for each of them
   subroutine PC_loop_ipart(pptr)
      procedure(if_ipart) :: pptr
      integer             :: n
      do n = 1, PC_num_part
         call pptr(PC_particles(n))
      end do
   end subroutine PC_loop_ipart

   !> Loop over all the particles and call pptr for each of them
   subroutine PC_loop_iopart(pptr)
      procedure(if_iopart) :: pptr
      integer              :: n
      do n = 1, PC_num_part
         call pptr(PC_particles(n))
      end do
   end subroutine PC_loop_iopart

   subroutine PC_compute_vector_sum(pptr, my_sum)
      procedure(if_ipart_oreals) :: pptr
      real(dp), intent(out)      :: my_sum(:)

      integer                    :: n
      real(dp)                   :: temp(size(my_sum))

      my_sum = 0.0_dp
      do n = 1, PC_num_part
         call pptr(PC_particles(n), temp)
         my_sum = my_sum + temp
      end do
   end subroutine PC_compute_vector_sum

   subroutine PC_compute_sum(pptr, my_sum)
      procedure(if_ipart_oreal) :: pptr
      real(dp), intent(out)     :: my_sum

      integer                   :: n
      real(dp)                  :: temp

      my_sum = 0.0_dp
      do n = 1, PC_num_part
         call pptr(PC_particles(n), temp)
         my_sum = my_sum + temp
      end do
   end subroutine PC_compute_sum

   subroutine ensure_storage_particles(req_size)
      integer, intent(in)         :: req_size
      integer                     :: new_size, old_size
      type(PC_part_t), allocatable :: copy_particles(:)

      ! These parameters control how the array grows
      real(dp), parameter :: increase_factor = 1.2_dp
      integer, parameter :: min_increase = 100*1000

      old_size = size(PC_particles)

      if (req_size > old_size) then
         new_size = max(int(increase_factor * req_size), req_size + min_increase)

         allocate(copy_particles(old_size))
         copy_particles = PC_particles
         deallocate(PC_particles)
         allocate(PC_particles(new_size))
         PC_particles(1:old_size) = copy_particles
         print *, "new size", size(PC_particles), req_size, old_size
      end if
   end subroutine ensure_storage_particles

   real(dp) elemental function PC_vel_to_en(vel)
      real(dp), intent(in) :: vel
      PC_vel_to_en = 0.5_dp * PC_particle_mass * vel**2
   end function PC_vel_to_en

   real(dp) elemental function PC_en_to_vel(en)
      real(dp), intent(in) :: en
      PC_en_to_vel = sqrt(2 * en / PC_particle_mass)
   end function PC_en_to_vel

   !> Create a lookup table with cross sections for a number of energies
   subroutine create_coll_rate_table(cross_secs, min_eV, max_eV, table_size)
      use m_units_constants
      use m_cross_sec
      use m_lookup_table
      type(CS_type), intent(in) :: cross_secs(:)
      integer, intent(IN)       :: table_size
      real(dp), intent(in)      :: min_eV, max_eV

      real(dp)                  :: vel_list(table_size), rate_list(table_size), sum_rate_list(table_size)
      integer                   :: ix, i_coll, i_row
      integer, allocatable      :: collIxs(:)
      real(dp)                  :: min_vel, max_vel, en_eV, temp

      min_vel = PC_en_to_vel(min_eV * UC_elec_volt)
      max_vel = PC_en_to_vel(max_eV * UC_elec_volt)

      PC_coll%num = size(cross_secs)
      allocate( collIxs(PC_coll%num) )
      allocate( PC_coll%types(PC_coll%num) )
      allocate( PC_coll%special_val(PC_coll%num) )

      ! Set a range of velocities
      do ix = 1, table_size
         temp = real(ix-1, dp) / (table_size-1)
         vel_list(ix) = min_vel + temp * (max_vel - min_vel)
      end do

      sum_rate_list = 0.0_dp

      ! Create collision rate table
      PC_coll%rate_lt = LT_create_mcol(min_vel, max_vel, table_size, PC_coll%num)

      do i_coll = 1, PC_coll%num
         PC_coll%types(i_coll) = cross_secs(i_coll)%col_type

         select case (PC_coll%types(i_coll))
         case (CS_ionize_t, CS_excite_t)
            ! Ionizations and excitations use a threshold energy which is lost
            PC_coll%special_val(i_coll) = cross_secs(i_coll)%spec_value * UC_elec_volt
         case (CS_elastic_t)
            ! Mass ration electron : neutral
            PC_coll%special_val(i_coll) = cross_secs(i_coll)%spec_value
         case DEFAULT
            PC_coll%special_val(i_coll) = 0.0_dp
         end select

         ! Linear interpolate cross sections by energy
         do i_row = 1, table_size
            en_eV = PC_vel_to_en(vel_list(i_row)) / UC_elec_volt
            call LT_lin_interp_list(cross_secs(i_coll)%en_cs(1, :), &
                 cross_secs(i_coll)%en_cs(2, :), en_eV, rate_list(i_row))
            rate_list(i_row) = rate_list(i_row) * vel_list(i_row)
         end do

         call LT_set_mcol(PC_coll%rate_lt, i_coll, vel_list, rate_list)
         sum_rate_list = sum_rate_list + rate_list
      end do

      ! Create a new table with the sums of the rates. This can speed up the null collision method.
      PC_coll%sum_rate_lt = LT_create_col(minval(vel_list), maxval(vel_list), table_size)
      call LT_set_col(PC_coll%sum_rate_lt, vel_list, sum_rate_list)

      PC_coll%max_rate = maxval(sum_rate_list)
      PC_coll%inv_max_rate = 1 / PC_coll%max_rate

   end subroutine create_coll_rate_table

   subroutine PC_sort_particles(sort_func)
      use m_mrgrnk
      procedure(if_freal_ipart) :: sort_func

      integer                      :: ix
      integer, allocatable         :: part_ixs(:)
      real(dp), allocatable        :: part_values(:)
      type(PC_part_t), allocatable :: part_copies(:)

      allocate(part_ixs(PC_num_part))
      allocate(part_values(PC_num_part))
      allocate(part_copies(PC_num_part))

      do ix = 1, PC_num_part
         part_values(ix) = sort_func(PC_particles(ix))
         part_copies(ix) = PC_particles(ix)
      end do

      call mrgrnk(part_values, part_ixs)

      do ix = 1, PC_num_part
         PC_particles(ix) = part_copies(part_ixs(ix))
      end do
   end subroutine PC_sort_particles

   subroutine PC_get_histogram(hist_func, filter_func, filter_args, x_values, y_values)
      use m_mrgrnk
      procedure(if_freal_ipart) :: hist_func
      real(dp), intent(in)      :: x_values(:)
      real(dp), intent(out)     :: y_values(:)
      procedure(if_filter_func) :: filter_func
      real(dp), intent(in)      :: filter_args(:)

      integer                   :: ix, p_ix, num_part, num_bins
      integer, allocatable      :: part_ixs(:)
      real(dp)                  :: boundary_value
      real(dp), allocatable     :: part_values(:)
      logical, allocatable      :: part_filter(:)

      allocate(part_filter(PC_num_part))
      do ix = 1, PC_num_part
         part_filter(ix) = filter_func(PC_particles(ix), filter_args)
      end do

      num_part = count(part_filter)

      allocate(part_values(num_part))
      allocate(part_ixs(num_part))

      p_ix = 0
      do ix = 1, PC_num_part
         if (part_filter(ix)) then
            p_ix = p_ix + 1
            part_values(p_ix) = hist_func(PC_particles(p_ix))
         end if
      end do

      call mrgrnk(part_values, part_ixs)

      num_bins = size(x_values)
      p_ix     = 1
      y_values = 0

      outer: do ix = 1, num_bins-1
         boundary_value = 0.5_dp * (x_values(ix) + x_values(ix+1))
         do
            if (p_ix == num_part + 1) exit outer
            if (part_values(part_ixs(p_ix)) > boundary_value) exit
            y_values(ix) = y_values(ix) + PC_particles(part_ixs(p_ix))%weight
            p_ix         = p_ix + 1
         end do
      end do outer

      ! Fill last bin
      y_values(num_bins) = sum(PC_particles(part_ixs(p_ix:num_part))%weight)
   end subroutine PC_get_histogram

   ! Routine to merge and split particles. Input arguments are the coordinate weights, used
   ! to determine the 'distance' between particles. The first three elements of the array are
   ! the weights of the xyz position coordinates, the next three the weights of the xyz
   ! velocity coordinates. Max_distance is the maxium Euclidean distance between particles
   ! to be merged. The weight_func returns the desired weight for a particle, whereas the
   ! pptr_merge and pptr_split procedures merge and split particles.
   subroutine PC_merge_and_split(coord_weights, max_distance, weight_func, pptr_merge, &
        pptr_split, ix_range, periodic_lengths, remove_dead)
      use m_mrgrnk
      use kdtree2_module
      real(dp), intent(in)           :: coord_weights(6), max_distance
      integer, intent(in), optional  :: ix_range(2)
      real(dp), intent(in), optional :: periodic_lengths(3)
      logical, intent(in), optional  :: remove_dead

      interface
         real(dp) function weight_func(my_part)
            import
            type(PC_part_t), intent(in) :: my_part
         end function weight_func

         subroutine pptr_merge(part_a, part_b)
            ! Takes in two particles, and marks one of them as inactive,
            ! while increasing the weight of the other.
            import
            type(PC_part_t), intent(inout) :: part_a, part_b
         end subroutine pptr_merge

         subroutine pptr_split(part_a, part_b)
            ! Split part_a into two, with the other 'half' put into part_b
            import
            type(PC_part_t), intent(inout) :: part_a, part_b
         end subroutine pptr_split
      end interface

      integer, parameter           :: num_neighbors = 1
      real(dp), parameter          :: large_ratio   = 1.5_dp, small_ratio = 1 / large_ratio
      real(dp)                     :: distance, per_fac(3)
      logical                      :: is_periodic(3), clean_up
      type(kdtree2), pointer       :: kd_tree
      type(kdtree2_result)         :: kd_results(num_neighbors)

      integer                      :: num_part, num_merge, num_split
      integer                      :: p_min, p_max, n_too_far
      integer                      :: ix, neighbor_ix, cntr, num_coords
      logical, allocatable         :: already_merged(:)
      integer, allocatable         :: part_ixs(:)
      real(dp), allocatable        :: coord_data(:, :), weight_ratios(:)
      type(PC_part_t), allocatable :: part_copy(:)

      if (.not. all(PC_particles(1:PC_num_part)%live)) then
         print *, "ERROR: There are dead particles when merging!"
         stop
      end if

      if (present(ix_range)) then
         p_min = ix_range(1)
         p_max = ix_range(2)
      else
         p_min = 1
         p_max = PC_num_part
      end if

      if (present(periodic_lengths)) then
         do ix = 1, 3
            is_periodic(ix) = (periodic_lengths(ix) /= 0.0_dp .and. coord_weights(ix) /= 0.0_dp)
            per_fac(ix) = periodic_lengths(ix) / (2 * acos(-1.0_dp))
         end do
      else
         is_periodic = .false.
      end if

      if (present(remove_dead)) then
         clean_up = remove_dead
      else
         clean_up = .true.
      end if

      num_part = p_max - p_min + 1
      allocate(weight_ratios(num_part))
      allocate(part_ixs(num_part))
      allocate(part_copy(num_part))

      do ix = 1, num_part
         weight_ratios(ix) = PC_particles(p_min+ix-1)%weight / weight_func(PC_particles(p_min+ix-1))
      end do

      num_merge      = count(weight_ratios <= small_ratio)
      num_coords     = count(coord_weights /= 0.0_dp) + count(is_periodic)
      allocate(coord_data(num_coords, num_merge))
      allocate(already_merged(num_merge))
      already_merged = .false.
      n_too_far      = 0

      ! Sort particles by their relative weight and store them in part_copy
      ! so that particles to be merged are at the beginning of the list
      call mrgrnk(weight_ratios, part_ixs)
      part_copy = PC_particles(p_min + part_ixs - 1)

      ! Only create a k-d tree if there are enough particles to be merged
      if (num_merge > num_coords) then
         ! Store the coordinates of the particles to be merged in coord_data
         cntr = 0
         do ix = 1, 6
            if (coord_weights(ix) /= 0.0_dp) then
               cntr = cntr + 1

               if (ix <= 3) then ! Spatial coordinates
                  if (is_periodic(ix)) then
                     ! Periodic coordinates are doubled
                     coord_data(cntr, :) = per_fac(ix) * &
                          cos(part_copy(1:num_merge)%x(ix) / per_fac(ix)) * coord_weights(ix)
                     cntr = cntr + 1
                     coord_data(cntr, :) = per_fac(ix) * &
                          sin(part_copy(1:num_merge)%x(ix) / per_fac(ix)) * coord_weights(ix)
                  else
                     coord_data(cntr, :) = part_copy(1:num_merge)%x(ix) * coord_weights(ix)
                  end if
               else             ! Velocity coordinates
                  coord_data(cntr, :) = part_copy(1:num_merge)%v(ix-3) * coord_weights(ix)
               end if
            end if
         end do

         ! Create k-d tree
         kd_tree => kdtree2_create(coord_data)

         ! Merge particles
         do ix = 1, num_merge
            if (already_merged(ix)) cycle

            call kdtree2_n_nearest_around_point(kd_tree, idxin=ix, nn=num_neighbors, correltime=1, results=kd_results)
            neighbor_ix = kd_results(1)%idx

            if (already_merged(neighbor_ix)) cycle
            distance = sqrt(sum((coord_data(:, ix) - coord_data(:, neighbor_ix))**2))
            if (distance > max_distance) then
               n_too_far = n_too_far + 1
               cycle
            end if
            call pptr_merge(part_copy(ix), part_copy(neighbor_ix))
            already_merged((/ix, neighbor_ix/)) = .true.
         end do
      end if

      ! Split particles. These are at the end of part_copy
      num_split = count(weight_ratios >= large_ratio)

      do ix = num_part - num_split + 1, num_part
         PC_num_part = PC_num_part + 1 ! An extra particle is created at the end of the particle list
         call ensure_storage_particles(PC_num_part)
         call pptr_split(part_copy(ix), PC_particles(PC_num_part))
      end do

      PC_particles(p_min:p_max) = part_copy ! Copy back the merged and split particles
      if (clean_up) call PC_clean_up_dead()
   end subroutine PC_merge_and_split

   subroutine PC_merge_part_rxv(part_a, part_b)
      use m_random
      type(PC_part_t), intent(inout) :: part_a, part_b

      real(dp) :: w1, w2, rr

      w1 = part_a%weight
      w2 = part_b%weight

      call random_number(rr)
      if (rr > w1 / (w1+w2)) then
         ! Keep particle b's position and velocity
         part_a%v      = part_b%v
         part_a%x      = part_b%x
      end if

      part_a%weight = w1 + w2
      part_b%live   = .false.
   end subroutine PC_merge_part_rxv

   subroutine PC_split_part(part_a, part_b)
      type(PC_part_t), intent(inout) :: part_a, part_b
      part_a%weight = part_a%weight * 0.5_dp
      part_b        = part_a
   end subroutine PC_split_part

   integer function PC_get_num_colls()
      PC_get_num_colls = PC_coll%num
   end function PC_get_num_colls

   subroutine PC_get_colls(out_colls)
      type(PC_coll_t), intent(out) :: out_colls
      out_colls = PC_coll
   end subroutine PC_get_colls

   !> Returns the number of collisions as real, to avoid overflow with standard integers
   subroutine PC_get_coll_count(collCount)
      real(dp), intent(out) :: collCount(:)
      collCount = PC_coll%counts
   end subroutine PC_get_coll_count

   !> Reset the number of collisions
   subroutine PC_reset_coll_count()
      PC_coll%counts(:) = 0
   end subroutine PC_reset_coll_count

   subroutine PC_get_coeffs(coeff_data, coeff_names, n_coeffs)
      use m_cross_sec
      use m_lookup_table
      real(dp), intent(out), allocatable          :: coeff_data(:,:)
      character(len=*), intent(out), allocatable :: coeff_names(:)
      integer, intent(out)                        :: n_coeffs
      type(PC_coll_t)                             :: coll_data
      integer                                     :: nn, n_rows

      call PC_get_colls(coll_data)
      n_coeffs = coll_data%num + 2
      n_rows = LT_get_num_rows_mcol(coll_data%rate_lt)
      allocate(coeff_data(n_coeffs, n_rows))
      allocate(coeff_names(n_coeffs))

      call LT_get_data_col(coll_data%sum_rate_lt, coeff_data(1, :), coeff_data(2,:))
      call LT_get_data_mcol(coll_data%rate_lt, coeff_data(1, :), coeff_data(3:,:))
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
