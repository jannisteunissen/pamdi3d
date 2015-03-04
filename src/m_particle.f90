module m_particle
  use m_particle_core

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)
  real(dp)           :: PM_max_weight
  real(dp)           :: PM_part_per_cell
  real(dp)           :: PM_v_rel_weight
  logical            :: PM_use_photoi

  abstract interface
     integer function p_to_int(my_part)
       import
       type(PC_part_t), intent(in) :: my_part
     end function p_to_int
  end interface

  type, extends(PC_bin_t) :: bin_t
    real(dp) :: x_min
    real(dp) :: inv_dx
  contains
    procedure :: bin_func => PM_bin_func
  end type bin_t

  type(bin_t), public :: PM_binner

  ! Public methods
  public :: PM_initialize
  public :: PM_create_ei_pair
  public :: PM_particles_to_density
  public :: PM_fld_error
  public :: PM_get_max_dt
  public :: PM_bin_func
  public :: PM_adjust_weights
  public :: PM_limit_dens

contains

  subroutine PM_initialize(pc, cross_secs, cfg, myrank, ntasks)
    use m_units_constants
    use m_cross_sec
    use m_config
    use m_phys_domain
    class(PC_t), intent(inout) :: pc
    type(CS_t), intent(in)    :: cross_secs(:)
    type(CFG_t), intent(in)   :: cfg
    integer, intent(in)       :: myrank, ntasks

    integer                   :: n_part_max, part_per_task, tbl_size
    real(dp)                  :: max_ev
    type(PC_part_t)           :: temp_part

    call CFG_get(cfg, "part_max_weight", PM_max_weight)
    call CFG_get(cfg, "part_per_cell", PM_part_per_cell)
    call CFG_get(cfg, "part_v_rel_weight", PM_v_rel_weight)
    call CFG_get(cfg, "part_max_num", n_part_max)
    call CFG_get(cfg, "part_lkp_tbl_size", tbl_size)
    call CFG_get(cfg, "part_max_ev", max_ev)
    call CFG_get(cfg, "photoi_enabled", PM_use_photoi)

    part_per_task = n_part_max / ntasks
    call pc%initialize(UC_elec_mass, cross_secs, tbl_size, &
         max_ev, part_per_task)
    call pc%set_coll_callback(coll_callback)
    call pc%set_outside_check(outside_check)

    PM_binner%n_bins = 5000
    PM_binner%inv_dx = (PM_binner%n_bins-1) / PD_r_max(1)
    if (myrank == 0) print *, "Particle model uses", &
         (storage_size(temp_part) / 2.0_dp**33) * n_part_max, "GB"
  end subroutine PM_initialize

  subroutine coll_callback(pc, my_part, c_ix, c_type)
    use m_cross_sec
    use m_efield_amr
    use m_photoi
    use m_phys_domain
    class(PC_t), intent(inout)  :: pc
    type(PC_part_t), intent(in) :: my_part
    integer, intent(in)         :: c_ix, c_type
    integer :: n, n_photons
    real(dp), allocatable :: photons(:,:)

    select case (c_type)
    case (CS_ionize_t)
       call E_add_to_var(E_i_pion, my_part%x, my_part%w)

       if (PM_use_photoi) then
          call pi_from_ionization(my_part, photons)

          n_photons = size(photons, 2)
          do n = 1, n_photons
             if (.not. PD_outside_domain(photons(:, n))) &
                  call PM_create_ei_pair(pc, photons(:, n))
          end do
       end if

    case (CS_attach_t)
       call E_add_to_var(E_i_O2m, my_part%x, my_part%w)
    end select
  end subroutine coll_callback

  logical function outside_check(my_part)
    use m_phys_domain
    type(PC_part_t), intent(in) :: my_part
    outside_check = PD_outside_domain(my_part%x)
  end function outside_check

  subroutine PM_create_ei_pair(pc, lx, v, a, w, t_left)
    use m_efield_amr
    class(PC_t), intent(inout)      :: pc
    real(dp), intent(in)           :: lx(3)
    real(dp), intent(in), optional :: v(3), a(3), w, t_left
    real(dp)                       :: lv(3), la(3), lw, lt_left

    lv      = 0; if (present(v)) lv = v
    la      = 0; if (present(a)) la = a
    lw      = 1; if (present(w)) lw = w
    lt_left = 0; if (present(t_left)) lt_left = t_left

    call pc%create_part(lx, lv, la, lw, lt_left)
    call E_add_to_var(E_i_pion, lx, lw)
  end subroutine PM_create_ei_pair

  subroutine PM_adjust_weights(pc)
    class(PC_t), intent(inout) :: pc
    call pc%merge_and_split((/.true., .true., .true./), PM_v_rel_weight, &
         .true., weight_func, PC_merge_part_rxv, split_part)
  end subroutine PM_adjust_weights

  function weight_func(my_part) result(weight)
    use m_efield_amr
    type(PC_part_t), intent(in) :: my_part
    real(dp)                    :: weight, n_elec
    n_elec = E_get_var_x_dr(E_i_elec, my_part%x)
    weight = n_elec / PM_part_per_cell
    weight = max(1.0_dp, min(PM_max_weight, weight))
  end function weight_func

  ! This routine can be used to estimate the error in the electric field. It
  ! should first be called with store_samples = .true., then a later call with
  ! store_samples = .false. sets the error in fld_err.
  subroutine PM_fld_error(pc, rng, n_samples, fld_err, store_samples)
    use m_efield_amr
    use m_random
    use mpi
    class(PC_t), intent(in) :: pc
    type(RNG_t), intent(inout) :: rng
    integer, intent(in) :: n_samples
    real(dp), intent(out) :: fld_err
    logical, intent(in) :: store_samples

    integer :: i, n, ierr
    real(dp) :: fld(3), this_err, avg_norm
    real(dp), allocatable, save :: pos_samples(:, :)
    real(dp), allocatable, save :: fld_samples(:, :)
    logical, save :: have_samples = .false.

    if (store_samples) then
       ! Store samples of the electric field

       if (pc%n_part < 1) then
          ! Not enough particles to store samples
          have_samples = .false.
       else
          have_samples = .true.
          if (allocated(fld_samples)) deallocate(fld_samples)
          if (allocated(pos_samples)) deallocate(pos_samples)
          allocate(fld_samples(3, n_samples))
          allocate(pos_samples(3, n_samples))

          do i = 1, n_samples
             ! Randomly select a particle
             n = rng%int_ab(1, pc%n_part)
             pos_samples(:,i) = pc%particles(n)%x
             fld_samples(:,i) = E_get_field(pos_samples(:,i))
          end do
       end if

    else if (.not. store_samples) then
       ! Compute the error
       fld_err = 0

       if (have_samples) then
          avg_norm = norm2(fld_samples) / sqrt(1.0_dp * n_samples)

          do i = 1, n_samples
             fld      = E_get_field(pos_samples(:,i))
             this_err = norm2(fld - fld_samples(:,i)) / &
                  max(norm2(fld), norm2(fld_samples(:,i)), avg_norm)
             fld_err  = max(fld_err, this_err)
          end do
       end if
    end if

    call MPI_ALLREDUCE(MPI_IN_PLACE, fld_err, 1, MPI_DOUBLE_PRECISION, &
         MPI_MAX, MPI_COMM_WORLD, ierr)
  end subroutine PM_fld_error

  subroutine PM_particles_to_density(pc)
    use m_efield_amr
    class(PC_t), intent(in) :: pc
    integer                :: ll

    call E_set_vars((/E_i_elec/), (/0.0_dp/))

    do ll = 1, pc%n_part
       call E_add_to_var(E_i_elec, pc%particles(ll)%x, pc%particles(ll)%w)
    end do
  end subroutine PM_particles_to_density

  function PM_get_max_dt(pc, rng, n_samples, cfl_num) result(dt_max)
    use m_random
    use m_efield_amr
    use mpi
    use m_mrgrnk
    class(PC_t), intent(in) :: pc
    type(RNG_t), intent(inout) :: rng
    integer, intent(in) :: n_samples
    real(dp), intent(in) :: cfl_num
    real(dp) :: dt_max

    real(dp) :: vel_est, min_dr
    real(dp), allocatable :: velocities(:)
    integer, allocatable :: ix_list(:)
    integer :: n, ix, ierr

    allocate(velocities(n_samples))
    allocate(ix_list(n_samples))

    if (pc%n_part == 0) then
       print *, "No particles, dt_max = 1.0e-12"
       dt_max = 1.0e-12_dp
       return
    end if

    ! Estimate maximum velocity of particles
    do n = 1, n_samples
       ix = rng%int_ab(1, pc%n_part)
       velocities(n) = norm2(pc%particles(ix)%v)
    end do

    call mrgrnk(velocities, ix_list)
    velocities = velocities(ix_list)

    vel_est = velocities(nint(n_samples * 0.9_dp))

    ! Get smallest grid delta
    min_dr = minval(E_get_smallest_dr())

    dt_max = cfl_num * min_dr / vel_est
    call MPI_ALLREDUCE(MPI_IN_PLACE, dt_max, 1, MPI_DOUBLE_PRECISION, &
         MPI_MIN, MPI_COMM_WORLD, ierr)
  end function PM_get_max_dt

  function PM_bin_func(binner, my_part) result(i_bin)
    class(bin_t), intent(in)    :: binner
    type(PC_part_t), intent(in) :: my_part
    integer                     :: i_bin

    i_bin  = 1 + int((my_part%x(1) - binner%x_min) * binner%inv_dx)
    i_bin  = max(1, i_bin)
    i_bin  = min(binner%n_bins, i_bin)
  end function PM_bin_func

  ! Perform attachment of electrons to stabilize the simulation in regions with
  ! a very high positive ion density.
  subroutine PM_limit_dens(pc, rng, max_dens)
    use m_units_constants
    use m_random
    use m_efield_amr
    class(PC_t), intent(inout)  :: pc
    type(RNG_t), intent(inout) :: rng
    real(dp), intent(in)       :: max_dens
    integer                    :: ll
    real(dp)                   :: local_dens, convert_prob, cur_max_dens(1)
    real(dp), parameter        :: convert_frac = 0.1_dp ! 10% converted per call

    call E_get_max_of_vars((/E_i_elec/), cur_max_dens)
    if (cur_max_dens(1) < max_dens) return

    print *, "Artificially attaching electrons for stability"
    print *, "Max elec dens", cur_max_dens, max_dens
    do ll = 1, pc%n_part
       local_dens = E_get_var(E_i_elec, pc%particles(ll)%x)

       if (local_dens > max_dens) then
          convert_prob = 1 - max_dens / local_dens
          convert_prob = convert_prob * convert_frac
          if (rng%uni_01() < convert_prob) then
             call E_add_to_var(E_i_nion, pc%particles(ll)%x, &
                  dble(pc%particles(ll)%w))
             call pc%remove_part(ll)
          end if
       end if
    end do

    call pc%clean_up()
  end subroutine PM_limit_dens

  subroutine split_part(part_a, w_ratio, part_out, n_part_out, rng)
    use m_random
    use m_efield_amr
    type(PC_part_t), intent(in)    :: part_a
    real(dp), intent(in)           :: w_ratio ! Current w / desired w
    type(PC_part_t), intent(inout) :: part_out(:)
    integer, intent(inout)         :: n_part_out
    type(RNG_t), intent(inout)     :: rng

    integer                        :: ix
    real(dp)                       :: dr(3), diff_dr

    ! The idea here is that each particle occupies a volume product(dr) /
    ! part_per_cell. The split-up particles are distributed over this volume.
    dr         = E_get_dr(part_a%x)
    diff_dr    = dr(1) / PM_part_per_cell**(1.0_dp/3.0_dp)
    n_part_out = min(nint(w_ratio), size(part_out))

    do ix = 1, n_part_out
       part_out(ix)   = part_a
       part_out(ix)%w = part_a%w / n_part_out
       part_out(ix)%x = part_out(ix)%x + &
            (/rng%uni_ab(-0.5_dp, 0.5_dp), rng%uni_ab(-0.5_dp, 0.5_dp), &
            rng%uni_ab(-0.5_dp, 0.5_dp)/) * diff_dr
    end do
  end subroutine split_part

end module m_particle
