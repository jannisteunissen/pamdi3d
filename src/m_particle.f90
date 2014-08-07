module m_particle
  use m_particle_core

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)
  real(dp) :: PM_max_weight
  real(dp) :: PM_part_per_cell
  real(dp) :: PM_coord_weights(6)
  real(dp) :: PM_max_distance
  logical :: PM_use_photoi

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

contains

  subroutine PM_initialize(pc, cross_secs, cfg, myrank, ntasks)
    use m_units_constants
    use m_cross_sec
    use m_config
    use m_phys_domain
    type(PC_t), intent(inout) :: pc
    type(CS_t), intent(in)    :: cross_secs(:)
    type(CFG_t), intent(in)   :: cfg
    integer, intent(in)       :: myrank, ntasks

    integer                   :: n_part_max, part_per_task, tbl_size
    real(dp)                  :: max_ev
    type(PC_part_t)           :: temp_part

    call CFG_get(cfg, "part_max_weight", PM_max_weight)
    call CFG_get(cfg, "part_per_cell", PM_part_per_cell)
    call CFG_get(cfg, "part_merge_coord_weights", PM_coord_weights)
    call CFG_get(cfg, "part_merge_max_distance", PM_max_distance)
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
    type(PC_t), intent(inout)      :: pc
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
    type(PC_t), intent(inout) :: pc
    call pc%merge_and_split(PM_coord_weights, PM_max_distance, &
         weight_func, PC_merge_part_rxv, PC_split_part)
  end subroutine PM_adjust_weights

  function weight_func(my_part) result(weight)
    use m_efield_amr
    type(PC_part_t), intent(in) :: my_part
    real(dp) :: weight, n_elec
    n_elec = E_get_var_x_dr(E_i_elec, my_part%x)
    weight = max(1.0_dp, min(PM_max_weight, n_elec / PM_part_per_cell))
  end function weight_func

  subroutine PM_fld_error(pc, rng, n_samples, fld_err, only_store)
    use m_efield_amr
    use m_random
    type(PC_t), intent(in) :: pc
    type(RNG_t), intent(inout) :: rng
    integer, intent(in) :: n_samples
    real(dp), intent(out) :: fld_err
    logical, intent(in) :: only_store

    integer :: i, n
    real(dp) :: fld(3), fld_sum
    real(dp), allocatable, save :: pos_samples(:, :)
    real(dp), allocatable, save :: fld_samples(:, :)

    if (pc%n_part < 1) then
       fld_err = 0
       return
    end if

    if (only_store) then
       fld_err = 0
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
    else
       fld_sum = 0
       do i = 1, n_samples
          fld = E_get_field(pos_samples(:,i))
          fld_err = max(fld_err, norm2(fld - fld_samples(:,i)))
          fld_sum = fld_sum + norm2(fld)
       end do
       fld_err = fld_err * n_samples / fld_sum
    end if
  end subroutine PM_fld_error

  subroutine PM_particles_to_density(pc)
    use m_efield_amr
    type(PC_t), intent(in) :: pc
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
    type(PC_t), intent(in) :: pc
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
         MPI_MAX, MPI_COMM_WORLD, ierr)
  end function PM_get_max_dt

  function PM_bin_func(binner, my_part) result(i_bin)
    class(bin_t), intent(in)    :: binner
    type(PC_part_t), intent(in) :: my_part
    integer                     :: i_bin

    i_bin  = 1 + int((my_part%x(1) - binner%x_min) * binner%inv_dx)
    i_bin  = max(1, i_bin)
    i_bin  = min(binner%n_bins, i_bin)
  end function PM_bin_func
end module m_particle
