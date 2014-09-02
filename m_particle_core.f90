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
  use m_cross_sec

  implicit none
  private

  integer, parameter  :: dp               = kind(0.0d0)
  real(dp), parameter :: PC_dead_weight   = -HUGE(1.0_dp)
  ! This has to do with openmp. It's quite interesting. Ask me about it.
  integer, parameter, public :: PC_max_num_coll = 100

  !> The particle type
  type, public :: PC_part_t
     real(dp) :: x(3)   = 0
     real(dp) :: v(3)   = 0
     real(dp) :: a(3)   = 0
     real(dp) :: w      = 0
     real(dp) :: t_left = 0
  end type PC_part_t

  type, public :: PC_t
     type(PC_part_t), allocatable :: particles(:)
     integer                      :: n_part
     type(CS_coll_t), allocatable :: colls(:)
     integer                      :: n_colls
     integer, allocatable         :: ionization_colls(:)
     integer, allocatable         :: attachment_colls(:)
     type(LT_table_t)             :: rate_lt
     real(dp)                     :: max_rate, inv_max_rate
     type(LL_int_head_t)          :: clean_list
     real(dp)                     :: mass
     type(RNG_t)                  :: rng
     integer                      :: separator(100)  ! Separate rng data

     procedure(p_to_logic_f), pointer, nopass :: outside_check => null()
     procedure(coll_callback_p), pointer      :: coll_callback => null()

   contains
     procedure, non_overridable :: initialize
     procedure, non_overridable :: set_coll_callback
     procedure, non_overridable :: set_outside_check
     procedure, non_overridable :: resize_part_list
     procedure, non_overridable :: remove_particles
     procedure, non_overridable :: advance
     procedure, non_overridable :: clean_up
     procedure, non_overridable :: create_part
     procedure, non_overridable :: add_part
     procedure, non_overridable :: remove_part
     procedure, non_overridable :: periodify
     procedure, non_overridable :: translate
     procedure, non_overridable :: get_mass
     procedure, non_overridable :: get_part
     procedure, non_overridable :: get_num_sim_part
     procedure, non_overridable :: get_num_real_part
     procedure, non_overridable :: set_accel
     procedure, non_overridable :: correct_new_accel
     procedure, non_overridable :: get_max_coll_rate
     procedure, non_overridable :: loop_iopart
     procedure, non_overridable :: compute_scalar_sum
     procedure, non_overridable :: compute_vector_sum
     procedure, non_overridable :: move_and_collide
     procedure, non_overridable :: set_coll_rates
     procedure, non_overridable :: get_mean_energy
     procedure, non_overridable :: get_coll_rates
     procedure, non_overridable :: check_space

     procedure, non_overridable :: merge_and_split
     procedure, non_overridable :: histogram

     procedure, non_overridable :: get_num_colls
     procedure, non_overridable :: get_colls
     procedure, non_overridable :: get_coeffs
  end type PC_t

  type, abstract, public :: PC_bin_t
     integer :: n_bins
   contains
     procedure(bin_f), deferred :: bin_func
  end type PC_bin_t

  abstract interface
     subroutine coll_callback_p(self, my_part, c_ix, c_type)
       import
       class(PC_t), intent(inout)  :: self
       type(PC_part_t), intent(in) :: my_part
       integer, intent(in)         :: c_ix, c_type
     end subroutine coll_callback_p

     subroutine p_to_r3_p(my_part, my_vec)
       import
       type(PC_part_t), intent(in) :: my_part
       real(dp), intent(out)       :: my_vec(3)
     end subroutine p_to_r3_p

     function p_to_r3_f(my_part) result(my_vec)
       import
       type(PC_part_t), intent(in) :: my_part
       real(dp)                    :: my_vec(3)
     end function p_to_r3_f

     real(dp) function p_to_r_f(my_part)
       import
       type(PC_part_t), intent(in) :: my_part
     end function p_to_r_f

     logical function p_to_logic_f(my_part)
       import
       type(PC_part_t), intent(in) :: my_part
     end function p_to_logic_f

     integer function bin_f(binner, my_part)
       import
       class(PC_bin_t), intent(in) :: binner
       type(PC_part_t), intent(in) :: my_part
     end function bin_f
  end interface

  ! Public procedures
  public :: PC_merge_part_rxv
  public :: PC_split_part
  public :: PC_v_to_en
  public :: PC_share
  public :: PC_reorder_by_bins

contains

  !> Initialization routine for the particle module
  subroutine initialize(self, mass, cross_secs, lookup_table_size, &
       max_en_eV, n_part_max)
    use m_cross_sec
    use m_units_constants
    class(PC_t), intent(inout) :: self
    type(CS_t), intent(in)     :: cross_secs(:)
    integer, intent(in)        :: lookup_table_size
    real(dp), intent(in)       :: mass, max_en_eV
    integer, intent(in)        :: n_part_max

    if (size(cross_secs) < 1) then
       print *, "No cross sections given, will abort"
       stop
    end if

    allocate(self%particles(n_part_max))
    self%mass   = mass
    self%n_part = 0

    call self%rng%set_seed((/1, 3, 3, 1337/))
    call self%set_coll_rates(cross_secs, mass, max_en_eV, lookup_table_size)
    print *, "Max coll rate", self%max_rate, self%inv_max_rate

    call get_colls_of_type(self, CS_ionize_t, self%ionization_colls)
    call get_colls_of_type(self, CS_attach_t, self%attachment_colls)
  end subroutine initialize

  subroutine set_coll_callback(self, pptr)
    class(PC_t), intent(inout) :: self
    procedure(coll_callback_p) :: pptr
    self%coll_callback => pptr
  end subroutine set_coll_callback

  subroutine set_outside_check(self, pptr)
    class(PC_t), intent(inout) :: self
    procedure(p_to_logic_f) :: pptr
    self%outside_check => pptr
  end subroutine set_outside_check

  subroutine get_colls_of_type(pc, ctype, ixs)
    class(PC_t), intent(in) :: pc
    integer, intent(in) :: ctype
    integer, intent(inout), allocatable :: ixs(:)
    integer :: nn, i

    nn = count(pc%colls(:)%type == ctype)
    allocate(ixs(nn))

    nn = 0
    do i = 1, pc%n_colls
       if (pc%colls(i)%type == ctype) then
          nn = nn + 1
          ixs(nn) = i
       end if
    end do
  end subroutine get_colls_of_type

  function get_mass(self) result(mass)
    class(PC_t), intent(in) :: self
    real(dp) :: mass
    mass = self%mass
  end function get_mass

  subroutine remove_particles(self)
    class(PC_t), intent(inout) :: self
    self%n_part = 0
    call LL_clear(self%clean_list)
  end subroutine remove_particles

  subroutine resize_part_list(self, new_size)
    class(PC_t), intent(inout) :: self
    integer, intent(in)          :: new_size
    type(PC_part_t), allocatable :: parts_copy(:)

    allocate(parts_copy(self%n_part))
    deallocate(self%particles)
    allocate(self%particles(new_size))
    self%particles(1:self%n_part) = parts_copy
  end subroutine resize_part_list

  subroutine advance(self, dt)
    class(PC_t), intent(inout) :: self
    real(dp), intent(in)           :: dt
    integer                        :: ll

    self%particles(1:self%n_part)%t_left = dt
    ll = 1

    do while (ll <= self%n_part)
       call self%move_and_collide(ll)
       ll = ll + 1
    end do

    call self%clean_up()
  end subroutine advance

  subroutine move_and_collide(self, ll)
    use m_cross_sec
    class(PC_t), intent(inout) :: self
    integer, intent(in)        :: ll
    integer                    :: cIx, cType, n_part_out
    real(dp)                   :: coll_time, new_vel
    integer, parameter         :: max_num_part_out = 2
    type(PC_part_t)            :: part_out(max_num_part_out)

    do
       ! Get the next collision time
       coll_time = sample_coll_time(self%rng%uni_01(), self%inv_max_rate)
       if (coll_time > self%particles(ll)%t_left) exit

       ! Set x,v at the collision time
       call advance_particle(self%particles(ll), coll_time)

       if (associated(self%outside_check)) then
          if (self%outside_check(self%particles(ll))) then
             call self%remove_part(ll)
             go to 100
          end if
       end if

       new_vel = norm2(self%particles(ll)%v)
       cIx     = get_coll_index(self%rate_lt, self%n_colls, self%max_rate, &
            new_vel, self%rng%uni_01())

       if (cIx > 0) then
          ! Perform the corresponding collision
          cType   = self%colls(cIx)%type

          if (associated(self%coll_callback)) &
               call self%coll_callback(self%particles(ll), cIx, cType)

          select case (cType)
          case (CS_attach_t)
             call attach_collision(self%particles(ll), part_out, &
                  n_part_out, self%colls(cIx), self%rng)
             go to 100 ! Particle is removed, so exit
          case (CS_elastic_t)
             call elastic_collision(self%particles(ll), part_out, &
                  n_part_out, self%colls(cIx), self%rng)
          case (CS_excite_t)
             call excite_collision(self%particles(ll), part_out, &
                  n_part_out, self%colls(cIx), self%rng)
          case (CS_ionize_t)
             call ionization_collision(self%particles(ll), part_out, &
                  n_part_out, self%colls(cIx), self%rng)
          case default
             stop "Wrong collision type"
          end select

          ! Store the particles returned. The first one goes at location ll, the
          ! rest at the end of the list.
          if (n_part_out == 0) then
             call self%remove_part(ll)
          else if (n_part_out == 1) then
             self%particles(ll) = part_out(1)
          else
             self%particles(ll) = part_out(1)
             self%particles(self%n_part+1:self%n_part+n_part_out-1) = &
                  part_out(2:n_part_out)
             self%n_part = self%n_part + n_part_out - 1
          end if
       end if
    end do

    ! Update the particle position and velocity to the next timestep
    call advance_particle(self%particles(ll), self%particles(ll)%t_left)
100 continue
  end subroutine move_and_collide

  !> Returns a sample from the exponential distribution of the collision times
  ! RNG_uniform() is uniform on [0,1), but log(0) = nan, so we take 1 - RNG_uniform()
  real(dp) function sample_coll_time(uni_01, inv_max_rate)
    real(dp), intent(in) :: uni_01, inv_max_rate
    sample_coll_time = -log(1 - uni_01) * inv_max_rate
  end function sample_coll_time

  integer function get_coll_index(rate_lt, n_colls, max_rate, &
       velocity, rand_unif)
    use m_find_index
    type(LT_table_t), intent(in) :: rate_lt
    real(dp), intent(IN)         :: velocity, rand_unif, max_rate
    integer, intent(in)          :: n_colls
    real(dp)                     :: rand_rate
    real(dp)                     :: buffer(PC_max_num_coll)

    ! Fill an array with interpolated rates
    buffer(1:n_colls) = LT_get_mcol(rate_lt, velocity)
    rand_rate         = rand_unif * max_rate
    get_coll_index    = FI_adaptive_r(buffer(1:n_colls), rand_rate)

    ! If there was no collision, the index exceeds the list and is set to 0
    if (get_coll_index == n_colls+1) get_coll_index = 0
  end function get_coll_index

  real(dp) function get_max_coll_rate(self)
    class(PC_t), intent(in) :: self
    get_max_coll_rate = self%max_rate
  end function get_max_coll_rate

  !> Perform an elastic collision for particle 'll'
  subroutine elastic_collision(part_in, part_out, n_part_out, coll, rng)
    type(PC_part_t), intent(in)        :: part_in
    type(PC_part_t), intent(inout)     :: part_out(:)
    integer, intent(out)               :: n_part_out
    type(CS_coll_t), intent(in)        :: coll
    type(RNG_t), intent(inout) :: rng

    real(dp)                           :: bg_vel(3), com_vel(3)

    ! TODO: implement random bg velocity
    bg_vel      = 0.0_dp
    n_part_out  = 1
    part_out(1) = part_in

    ! Compute center of mass velocity
    com_vel = (coll%rel_mass * part_out(1)%v + bg_vel) / (1 + coll%rel_mass)

    ! Scatter in center of mass coordinates
    part_out(1)%v = part_out(1)%v - com_vel
    call scatter_isotropic(part_out(1), norm2(part_out(1)%v), rng)
    part_out(1)%v = part_out(1)%v + com_vel
  end subroutine elastic_collision

  !> Perform an excitation-collision for particle 'll'
  subroutine excite_collision(part_in, part_out, n_part_out, coll, rng)
    use m_units_constants
    type(PC_part_t), intent(in)        :: part_in
    type(PC_part_t), intent(inout)     :: part_out(:)
    integer, intent(out)               :: n_part_out
    type(CS_coll_t), intent(in)        :: coll
    type(RNG_t), intent(inout) :: rng

    real(dp)             :: energy, old_en, new_vel

    old_en  = PC_v_to_en(part_in%v, coll%part_mass)
    energy  = max(0.0_dp, old_en - coll%en_loss)
    new_vel = PC_en_to_vel(energy, coll%part_mass)

    n_part_out  = 1
    part_out(1) = part_in
    call scatter_isotropic(part_out(1), new_vel, rng)
  end subroutine excite_collision

  !> Perform an ionizing collision for particle 'll'
  subroutine ionization_collision(part_in, part_out, n_part_out, coll, rng)
    use m_units_constants
    type(PC_part_t), intent(in)        :: part_in
    type(PC_part_t), intent(inout)     :: part_out(:)
    integer, intent(out)               :: n_part_out
    type(CS_coll_t), intent(in)        :: coll
    type(RNG_t), intent(inout) :: rng

    real(dp)                       :: energy, old_en, velocity

    old_en      = PC_v_to_en(part_in%v, coll%part_mass)
    energy      = max(0.0_dp, old_en - coll%en_loss)
    velocity    = PC_en_to_vel(0.5_dp * energy, coll%part_mass)

    n_part_out  = 2
    part_out(1) = part_in
    part_out(2) = part_in
    call scatter_isotropic(part_out(1), velocity, rng)
    call scatter_isotropic(part_out(2), velocity, rng)
  end subroutine ionization_collision

  !> Perform attachment of electron 'll'
  subroutine attach_collision(part_in, part_out, n_part_out, coll, rng)
    use m_units_constants
    type(PC_part_t), intent(in)        :: part_in
    type(PC_part_t), intent(inout)     :: part_out(:)
    integer, intent(out)               :: n_part_out
    type(CS_coll_t), intent(in)        :: coll
    type(RNG_t), intent(inout) :: rng
    n_part_out = 0
  end subroutine attach_collision

  subroutine scatter_isotropic(part, vel_norm, rng)
    type(PC_part_t), intent(inout) :: part
    real(dp), intent(in)           :: vel_norm
    type(RNG_t), intent(inout)     :: rng
    real(dp)                       :: sum_sq, tmp_sqrt, rands(2)

    ! Marsaglia method for uniform sampling on sphere
    do
       rands(1) = rng%uni_ab(-1.0_dp, 1.0_dp)
       rands(2) = rng%uni_ab(-1.0_dp, 1.0_dp)
       sum_sq   = rands(1)**2 + rands(2)**2
       if (sum_sq <= 1) exit
    end do

    tmp_sqrt    = sqrt(1 - sum_sq)
    part%v(1:2) = 2 * rands(1:2) * tmp_sqrt
    part%v(3)   = 1 - 2 * sum_sq
    part%v      = part%v * vel_norm ! Normalization
  end subroutine scatter_isotropic

  !> Advance the particle position and velocity over time dt
  subroutine advance_particle(part, dt)
    type(PC_part_t), intent(inout) :: part
    real(dp), intent(in) :: dt

    part%x      = part%x + part%v * dt + &
         0.5_dp * part%a * dt**2
    part%v      = part%v + part%a * dt
    part%t_left = part%t_left - dt
  end subroutine advance_particle

  subroutine set_accel(self, accel_func)
    class(PC_t), intent(inout) :: self
    procedure(p_to_r3_f)    :: accel_func
    integer                    :: ll

    do ll = 1, self%n_part
       self%particles(ll)%a = accel_func(self%particles(ll))
    end do
  end subroutine set_accel

  !> Correct particle velocities for the previous timestep of 'dt'
  !!
  !! During the timestep x,v have been advanced to:
  !! x(t+1) = x(t) + v(t)*dt + 0.5*a(t)*dt^2,
  !! v(t+1) = v(t) + a(t)*dt
  !! But the velocity at t+1 should be v(t+1) = v(t) + 0.5*(a(t) + a(t+1))*dt,
  !! to have a second order leapfrog scheme, so here we set it to that value.
  subroutine correct_new_accel(self, dt, accel_func)
    use m_units_constants
    class(PC_t), intent(inout) :: self
    real(dp), intent(IN)      :: dt
    procedure(p_to_r3_p) :: accel_func
    integer                   :: ll
    real(dp)                  :: new_accel(3)

    do ll = 1, self%n_part
       call accel_func(self%particles(ll), new_accel)
       self%particles(ll)%v = self%particles(ll)%v + &
            0.5_dp * (new_accel - self%particles(ll)%a) * dt
       self%particles(ll)%a = new_accel
    end do
  end subroutine correct_new_accel

  subroutine clean_up(self)
    class(PC_t), intent(inout) :: self
    integer :: ix_end, ix_clean, n_part
    logical :: success

    do
       ! Get an index that has to be cleaned from the list
       call LL_pop(self%clean_list, ix_clean, success)
       if (.not. success) exit

       n_part      = self%n_part
       ! This is overridden if a replacement is found
       self%n_part = min(self%n_part, ix_clean-1)

       ! Find the last "alive" particle in the list
       do ix_end = n_part, ix_clean+1, -1
          if (self%particles(ix_end)%w /= PC_dead_weight) then
             self%particles(ix_clean) = self%particles(ix_end)
             self%n_part              = ix_end-1
             exit
          end if
       end do
    end do
  end subroutine clean_up

  subroutine add_part(self, part)
    class(PC_t), intent(inout)  :: self
    type(PC_part_t), intent(in) :: part
    self%n_part                 = self%n_part + 1
    self%particles(self%n_part) = part
  end subroutine add_part

  subroutine create_part(self, x, v, a, weight, t_left)
    class(PC_t), intent(inout) :: self
    real(dp), intent(IN)       :: x(3), v(3), a(3), weight, t_left
    type(PC_part_t)            :: my_part
    my_part%x      = x
    my_part%v      = v
    my_part%a      = a
    my_part%w      = weight
    my_part%t_left = t_left
    call self%add_part(my_part)
  end subroutine create_part

  !> Mark particle for removal
  subroutine remove_part(self, ix_to_remove)
    class(PC_t), intent(inout) :: self
    integer, intent(in) :: ix_to_remove

    call LL_add(self%clean_list, ix_to_remove)
    self%particles(ix_to_remove)%w = PC_dead_weight
  end subroutine remove_part

  subroutine periodify(self, is_periodic, lengths)
    class(PC_t), intent(inout) :: self
    logical, intent(in) :: is_periodic(3)
    real(dp), intent(in) :: lengths(3)
    integer :: i_dim, n_part

    n_part = self%n_part
    do i_dim = 1, 3
       if (is_periodic(i_dim)) then
          self%particles(1:n_part)%x(i_dim) = &
               modulo(self%particles(1:n_part)%x(i_dim), lengths(i_dim))
       end if
    end do
  end subroutine periodify

  subroutine translate(self, delta_x)
    class(PC_t), intent(inout) :: self
    real(dp), intent(in) :: delta_x(3)
    integer              :: ll
    do ll = 1, self%n_part
       self%particles(ll)%x = self%particles(ll)%x + delta_x
    end do
  end subroutine translate

  real(dp) function get_part_mass(self)
    class(PC_t), intent(inout) :: self
    get_part_mass = self%mass
  end function get_part_mass

  function get_part(self, ll) result(part)
    class(PC_t), intent(inout) :: self
    integer, intent(in) :: ll
    type(PC_part_t) :: part
    part = self%particles(ll)
  end function get_part

  !> Return the number of real particles
  real(dp) function get_num_real_part(self)
    class(PC_t), intent(in) :: self
    get_num_real_part = sum(self%particles(1:self%n_part)%w)
  end function get_num_real_part

  !> Return the number of simulation particles
  integer function get_num_sim_part(self)
    class(PC_t), intent(in) :: self
    get_num_sim_part = self%n_part
  end function get_num_sim_part

  !> Loop over all the particles and call pptr for each of them
  subroutine loop_iopart(self, pptr)
    class(PC_t), intent(inout) :: self
    interface
       subroutine pptr(my_part)
         import
         type(PC_part_t), intent(inout) :: my_part
       end subroutine pptr
    end interface

    integer :: ll

    do ll = 1, self%n_part
       call pptr(self%particles(ll))
    end do
  end subroutine loop_iopart

  subroutine compute_vector_sum(self, pptr, my_sum)
    class(PC_t), intent(inout) :: self
    interface
       subroutine pptr(my_part, my_reals)
         import
         type(PC_part_t), intent(in) :: my_part
         real(dp), intent(out)       :: my_reals(:)
       end subroutine pptr
    end interface
    real(dp), intent(out)      :: my_sum(:)

    integer                    :: ll
    real(dp)                   :: temp(size(my_sum))

    my_sum = 0.0_dp
    do ll = 1, self%n_part
       call pptr(self%particles(ll), temp)
       my_sum = my_sum + temp
    end do
  end subroutine compute_vector_sum

  subroutine compute_scalar_sum(self, pptr, my_sum)
    class(PC_t), intent(in) :: self
    interface
       subroutine pptr(my_part, my_real)
         import
         type(PC_part_t), intent(in) :: my_part
         real(dp), intent(out)       :: my_real
       end subroutine pptr
    end interface
    real(dp), intent(out)     :: my_sum

    integer                   :: ll
    real(dp)                  :: temp

    my_sum = 0.0_dp
    do ll = 1, self%n_part
       call pptr(self%particles(ll), temp)
       my_sum = my_sum + temp
    end do
  end subroutine compute_scalar_sum

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

  subroutine check_space(self, n_req)
    class(PC_t), intent(in) :: self
    integer, intent(in) :: n_req
    if (n_req > size(self%particles)) then
       print *, "Particle list too small!"
       stop
    end if
  end subroutine check_space

  !> Create a lookup table with cross sections for a number of energies
  subroutine set_coll_rates(self, cross_secs, mass, max_ev, table_size)
    use m_units_constants
    use m_cross_sec
    use m_lookup_table
    class(PC_t), intent(inout) :: self
    type(CS_t), intent(in)     :: cross_secs(:)
    integer, intent(in)        :: table_size
    real(dp), intent(in)       :: mass, max_ev

    real(dp)                  :: vel_list(table_size), rate_list(table_size)
    real(dp)                  :: sum_rate_list(table_size)
    integer                   :: ix, i_c, i_row, n_colls
    real(dp)                  :: max_vel, en_eV

    max_vel      = PC_en_to_vel(max_ev * UC_elec_volt, mass)
    n_colls      = size(cross_secs)
    self%n_colls = n_colls
    allocate(self%colls(n_colls))

    ! Set a range of velocities
    do ix = 1, table_size
       vel_list(ix) = (ix-1) * max_vel / (table_size-1)
       sum_rate_list(ix) = 0
    end do

    ! Create collision rate table
    self%rate_lt = LT_create(0.0_dp, max_vel, table_size, n_colls)

    do i_c = 1, n_colls
       self%colls(i_c) = cross_secs(i_c)%coll

       ! Linear interpolate cross sections by energy
       do i_row = 1, table_size
          en_eV = PC_speed_to_en(vel_list(i_row), mass) / UC_elec_volt
          call LT_lin_interp_list(cross_secs(i_c)%en_cs(1, :), &
               cross_secs(i_c)%en_cs(2, :), en_eV, rate_list(i_row))
          rate_list(i_row) = rate_list(i_row) * vel_list(i_row)
       end do

       sum_rate_list = sum_rate_list + rate_list
       call LT_set_col(self%rate_lt, i_c, vel_list, sum_rate_list)
    end do

    self%max_rate = maxval(sum_rate_list)
    self%inv_max_rate = 1 / self%max_rate
  end subroutine set_coll_rates

  subroutine sort(self, sort_func)
    use m_mrgrnk
    class(PC_t), intent(inout)   :: self
    procedure(p_to_r_f)          :: sort_func

    integer                      :: ix, n_part
    integer, allocatable         :: sorted_ixs(:)
    real(dp), allocatable        :: part_values(:)
    type(PC_part_t), allocatable :: part_copies(:)

    n_part = self%n_part
    allocate(sorted_ixs(n_part))
    allocate(part_values(n_part))
    allocate(part_copies(n_part))

    do ix = 1, n_part
       part_copies(ix) = self%particles(ix)
       part_values(ix) = sort_func(part_copies(ix))
    end do

    call mrgrnk(part_values, sorted_ixs)

    do ix = 1, n_part
       self%particles(ix) = part_copies(sorted_ixs(ix))
    end do
  end subroutine sort

  subroutine histogram(self, hist_func, filter_f, &
       filter_args, x_values, y_values)
    use m_mrgrnk
    class(PC_t), intent(in) :: self
    procedure(p_to_r_f)   :: hist_func
    real(dp), intent(in)        :: x_values(:)
    real(dp), intent(out)       :: y_values(:)
    interface
       logical function filter_f(my_part, real_args)
         import
         type(PC_part_t), intent(in) :: my_part
         real(dp), intent(in) :: real_args(:)
       end function filter_f
    end interface
    real(dp), intent(in)        :: filter_args(:)

    integer                     :: ix, p_ix, o_ix, n_used, num_bins, n_part
    integer, allocatable        :: sorted_ixs(:)
    real(dp)                    :: boundary_value
    real(dp), allocatable       :: part_values(:)
    logical, allocatable        :: part_mask(:)

    n_part = self%n_part
    allocate(part_mask(n_part))
    do ix = 1, n_part
       part_mask(ix) = filter_f(self%particles(ix), filter_args)
    end do

    n_used = count(part_mask)

    allocate(part_values(n_used))
    allocate(sorted_ixs(n_used))

    p_ix = 0
    do ix = 1, self%n_part
       if (part_mask(ix)) then
          p_ix = p_ix + 1
          part_values(p_ix) = hist_func(self%particles(p_ix))
       end if
    end do

    call mrgrnk(part_values, sorted_ixs)

    num_bins = size(x_values)
    p_ix     = 1
    y_values = 0

    outer: do ix = 1, num_bins-1
       boundary_value = 0.5_dp * (x_values(ix) + x_values(ix+1))
       do
          if (p_ix == n_used + 1) exit outer
          o_ix = sorted_ixs(p_ix) ! Index in the 'old' list
          if (part_values(o_ix) > boundary_value) exit

          y_values(ix) = y_values(ix) + self%particles(o_ix)%w
          p_ix         = p_ix + 1
       end do
    end do outer

    ! Fill last bin
    y_values(num_bins) = sum(self%particles(sorted_ixs(p_ix:n_used))%w)
  end subroutine histogram

  ! Routine to merge and split particles. Input arguments are the coordinate
  ! weights, used to determine the 'distance' between particles. The first three
  ! elements of the array are the weights of the xyz position coordinates, the
  ! next three the weights of the xyz velocity coordinates. Max_distance is the
  ! maxium Euclidean distance between particles to be merged. The weight_func
  ! returns the desired weight for a particle, whereas the pptr_merge and
  ! pptr_split procedures merge and split particles.
  subroutine merge_and_split(self, x_mask, v_fac, use_v_norm, weight_func, &
       pptr_merge, pptr_split)
    use m_mrgrnk
    use kdtree2_module
    class(PC_t), intent(inout) :: self
    real(dp), intent(in)       :: v_fac
    logical, intent(in)        :: x_mask(3), use_v_norm

    interface
       subroutine pptr_merge(part_a, part_b, part_out, rng)
         import
         type(PC_part_t), intent(in)  :: part_a, part_b
         type(PC_part_t), intent(out) :: part_out
         type(RNG_t), intent(inout)   :: rng
       end subroutine pptr_merge

       subroutine pptr_split(part_a, w_ratio, part_out, n_part_out, rng)
         import
         type(PC_part_t), intent(in)    :: part_a
         real(dp), intent(in)           :: w_ratio
         type(PC_part_t), intent(inout) :: part_out(:)
         integer, intent(inout)         :: n_part_out
         type(RNG_t), intent(inout)     :: rng
       end subroutine pptr_split
    end interface

    procedure(p_to_r_f)    :: weight_func

    integer, parameter     :: num_neighbors  = 1
    integer, parameter     :: n_part_out_max = 16
    real(dp), parameter    :: large_ratio    = 1.5_dp
    real(dp), parameter    :: small_ratio    = 1 / large_ratio
    type(kdtree2), pointer :: kd_tree
    type(kdtree2_result)   :: kd_results(num_neighbors)

    integer               :: n_x_coord, n_coords
    integer               :: num_part, num_merge, num_split
    integer               :: p_min, p_max, n_too_far
    integer               :: o_ix, o_nn_ix
    integer               :: i, ix, neighbor_ix
    integer               :: n_part_out
    logical, allocatable  :: already_merged(:)
    integer, allocatable  :: sorted_ixs(:), coord_ixs(:)
    real(dp), allocatable :: coord_data(:, :), weight_ratios(:)
    type(PC_part_t)       :: part_out(n_part_out_max)

    p_min                                    = 1
    p_max                                    = self%n_part
    num_part                                 = p_max - p_min + 1

    allocate(weight_ratios(num_part))
    allocate(sorted_ixs(num_part))

    do ix = 1, num_part
       weight_ratios(ix) = self%particles(p_min+ix-1)%w / &
            weight_func(self%particles(p_min+ix-1))
    end do

    num_merge = count(weight_ratios <= small_ratio)
    n_x_coord = count(x_mask)
    n_coords  = n_x_coord + 3
    if (use_v_norm) n_coords = n_coords - 2

    allocate(coord_data(n_coords, num_merge))
    allocate(coord_ixs(n_coords))
    allocate(already_merged(num_merge))
    already_merged = .false.
    n_too_far      = 0

    ! print *, "num_merge", num_merge
    ! Sort particles by their relative weight
    call mrgrnk(weight_ratios, sorted_ixs)

    ! Only create a k-d tree if there are enough particles to be merged
    if (num_merge > n_coords) then
       do ix = 1, num_merge
          o_ix = sorted_ixs(ix)
          coord_data(1:n_x_coord, ix) = pack(self%particles(o_ix)%x, x_mask)
          if (use_v_norm) then
             coord_data(n_x_coord+1, ix) = v_fac * norm2(self%particles(o_ix)%v)
          else
             coord_data(n_x_coord+1:, ix) = v_fac * self%particles(o_ix)%v
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

          ! Get indices in the original particle list
          o_ix = sorted_ixs(ix)
          o_nn_ix = sorted_ixs(neighbor_ix)

          ! Merge, then remove neighbor
          call pptr_merge(self%particles(o_ix), self%particles(o_nn_ix), &
               part_out(1), self%rng)
          self%particles(o_ix) = part_out(1)
          call self%remove_part(o_nn_ix)
          already_merged((/ix, neighbor_ix/)) = .true.
       end do

       call kdtree2_destroy(kd_tree)
    end if

    ! Split particles.
    num_split = count(weight_ratios >= large_ratio)
    ! print *, "num_split:", num_split

    do ix = num_part - num_split + 1, num_part
       o_ix = sorted_ixs(ix)
       call pptr_split(self%particles(o_ix), weight_ratios(o_ix), part_out, &
            n_part_out, self%rng)
       self%particles(o_ix) = part_out(1)
       do i = 2, n_part_out
          call self%add_part(part_out(i))
       end do
    end do

    ! print *, "clean up"
    call self%clean_up()
  end subroutine merge_and_split

  ! Merge two particles into part_a, should remove part_b afterwards
  subroutine PC_merge_part_rxv(part_a, part_b, part_out, rng)
    type(PC_part_t), intent(in)    :: part_a, part_b
    type(PC_part_t), intent(out) :: part_out
    type(RNG_t), intent(inout) :: rng

    if (rng%uni_01() > part_a%w / (part_a%w + part_b%w)) then
       part_out%x      = part_b%x
       part_out%v      = part_b%v
       part_out%a      = part_b%a
       part_out%t_left = part_b%t_left
    else
       part_out%x      = part_a%x
       part_out%v      = part_a%v
       part_out%a      = part_a%a
       part_out%t_left = part_a%t_left
    end if
    part_out%w = part_a%w + part_b%w
  end subroutine PC_merge_part_rxv

  subroutine PC_split_part(part_a, w_ratio, part_out, n_part_out, rng)
    type(PC_part_t), intent(in)    :: part_a
    real(dp), intent(in)           :: w_ratio
    type(PC_part_t), intent(inout) :: part_out(:)
    integer, intent(inout)         :: n_part_out
    type(RNG_t), intent(inout)     :: rng

    n_part_out    = 2
    part_out(1)   = part_a
    part_out(1)%w = 0.5_dp * part_out(1)%w
    part_out(2)   = part_out(1)
  end subroutine PC_split_part

  integer function get_num_colls(self)
    class(PC_t), intent(in) :: self
    get_num_colls = self%n_colls
  end function get_num_colls

  function get_mean_energy(self) result(mean_en)
    class(PC_t), intent(in) :: self
    real(dp)                :: mean_en, weight
    integer                 :: ll
    mean_en = 0
    weight  = 0
    do ll = 1, self%n_part
       weight  = weight + self%particles(ll)%w
       mean_en = mean_en + self%particles(ll)%w * &
            PC_v_to_en(self%particles(ll)%v, self%mass)
    end do
    mean_en = mean_en / weight
  end function get_mean_energy

  subroutine get_coll_rates(self, velocity, coll_rates)
    class(PC_t), intent(in) :: self
    real(dp), intent(in) :: velocity
    real(dp), intent(inout) :: coll_rates(:)
    integer :: i
    coll_rates = LT_get_mcol(self%rate_lt, velocity)

    do i = size(coll_rates), 2, -1
       coll_rates(i) = coll_rates(i) - coll_rates(i-1)
    end do
  end subroutine get_coll_rates

  subroutine get_colls(self, out_colls)
    class(PC_t), intent(in) :: self
    type(CS_coll_t), intent(inout), allocatable :: out_colls(:)
    allocate(out_colls(self%n_colls))
    out_colls = self%colls
  end subroutine get_colls

  subroutine get_coeffs(self, coeff_data, coeff_names, n_coeffs)
    use m_cross_sec
    use m_lookup_table
    class(PC_t), intent(in)                    :: self
    real(dp), intent(out), allocatable         :: coeff_data(:,:)
    character(len=*), intent(out), allocatable :: coeff_names(:)
    integer, intent(out)                       :: n_coeffs
    type(CS_coll_t), allocatable               :: coll_data(:)
    integer                                    :: nn, n_rows

    call self%get_colls(coll_data)
    n_coeffs = self%n_colls + 1
    n_rows = LT_get_num_rows(self%rate_lt)
    allocate(coeff_data(n_coeffs, n_rows))
    allocate(coeff_names(n_coeffs))

    call LT_get_data(self%rate_lt, coeff_data(1, :), coeff_data(2:,:))
    coeff_names(1) = "velocity (m/s)"
    do nn = 1, self%n_colls
       select case (self%colls(nn)%type)
       case (CS_ionize_t)
          coeff_names(1+nn) = "ionization"
       case (CS_attach_t)
          coeff_names(1+nn) = "attachment"
       case (CS_elastic_t)
          coeff_names(1+nn) = "elastic"
       case (CS_excite_t)
          coeff_names(1+nn) = "excitation"
       case default
          coeff_names(1+nn) = "unknown"
       end select
    end do
  end subroutine get_coeffs

  ! Share particles between PC_t objects
  subroutine PC_share(pcs)
    type(PC_t), intent(inout) :: pcs(:)

    integer :: i_temp(1), n_pc
    integer :: n_avg, i_min, i_max, n_min, n_max, n_send

    n_pc = size(pcs)
    n_avg = ceiling(sum(pcs(:)%n_part) / real(n_pc, dp))

    do
       i_temp = maxloc(pcs(:)%n_part)
       i_max  = i_temp(1)
       n_max  = pcs(i_max)%n_part
       i_temp = minloc(pcs(:)%n_part)
       i_min  = i_temp(1)
       n_min  = pcs(i_min)%n_part

       ! Difference it at most n_pc - 1, if all lists get one more particle
       ! than the last list
       if (n_max - n_min < n_pc) exit

       ! Send particles from i_max to i_min
       n_send = min(n_max - n_avg, n_avg - n_min)
       pcs(i_min)%particles(n_min+1:n_min+n_send) = &
            pcs(i_max)%particles(n_max-n_send+1:n_max)

       ! Always at the end of a list, so do not need to clean up later
       pcs(i_min)%n_part = pcs(i_min)%n_part + n_send
       pcs(i_max)%n_part = pcs(i_max)%n_part - n_send
    end do
  end subroutine PC_share

  subroutine PC_reorder_by_bins(pcs, binner)
    use m_mrgrnk
    type(PC_t), intent(inout) :: pcs(:)
    class(PC_bin_t), intent(in) :: binner
    integer, allocatable :: bin_counts(:,:)
    integer, allocatable :: bin_counts_sum(:)
    integer, allocatable :: bin_owner(:)
    integer :: prev_num_part(size(pcs))
    integer :: n_pc, n_avg, tmp_sum
    integer :: ll, io, ib, ip, new_loc

    type tmp_t
       integer, allocatable :: ib(:)
    end type tmp_t
    type(tmp_t), allocatable :: pcs_bins(:)

    n_pc = size(pcs)
    n_avg = ceiling(sum(pcs(:)%n_part) / real(n_pc, dp))

    allocate(bin_counts(binner%n_bins, n_pc))
    allocate(bin_counts_sum(binner%n_bins))
    allocate(bin_owner(binner%n_bins))
    allocate(pcs_bins(n_pc))
    bin_counts(:, :) = 0

    ! Get the counts in the bins for each pcs(ip)
    do ip = 1, n_pc
       allocate(pcs_bins(ip)%ib(pcs(ip)%n_part))
       do ll = 1, pcs(ip)%n_part
          ib                  = binner%bin_func(pcs(ip)%particles(ll))
          pcs_bins(ip)%ib(ll) = ib
          bin_counts(ib, ip)  = bin_counts(ib, ip) + 1
       end do
    end do

    bin_counts_sum = sum(bin_counts, dim=2)
    tmp_sum        = 0
    ip             = 1

    ! Set the owners of the bins
    do ib = 1, binner%n_bins
       tmp_sum = tmp_sum + bin_counts_sum(ib)
       bin_owner(ib) = ip
       if (tmp_sum >= n_avg) then
          ip = ip + 1
          tmp_sum = 0
       end if
    end do

    prev_num_part(:) = pcs(:)%n_part

    do ip = 1, n_pc
       do ll = 1, prev_num_part(ip)
          ib = pcs_bins(ip)%ib(ll)
          io = bin_owner(ib)
          if (io /= ip) then
             ! Insert at owner
             new_loc = pcs(io)%n_part + 1
             pcs(io)%particles(new_loc) = pcs(ip)%particles(ll)
             pcs(io)%n_part = new_loc
             call remove_part(pcs(ip), ll)
          end if
       end do
    end do

    do ip = 1, n_pc
       call pcs(ip)%clean_up()
    end do
  end subroutine PC_reorder_by_bins

end module m_particle_core
