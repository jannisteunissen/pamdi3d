module m_particle
  use m_particle_core

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)
  real(dp) :: PM_max_weight
  real(dp) :: PM_part_per_cell
  real(dp) :: PM_coord_weights(6)
  real(dp) :: PM_max_distance

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

  subroutine PM_initialize(pc, cross_secs, cfg)
    use m_units_constants
    use m_cross_sec
    use m_config
    use m_phys_domain
    type(PC_t), intent(inout) :: pc
    type(CS_t), intent(in)    :: cross_secs(:)
    type(CFG_t), intent(in)   :: cfg

    integer                   :: n_part_max, tbl_size
    real(dp) :: max_ev

    call CFG_get(cfg, "part_max_weight", PM_max_weight)
    call CFG_get(cfg, "part_per_cell", PM_part_per_cell)
    call CFG_get(cfg, "part_merge_coord_weights", PM_coord_weights)
    call CFG_get(cfg, "part_merge_max_distance", PM_max_distance)
    call CFG_get(cfg, "part_max_num", n_part_max)
    call CFG_get(cfg, "part_lkp_tbl_size", tbl_size)
    call CFG_get(cfg, "part_max_ev", max_ev)

    call pc%initialize(UC_elec_mass, cross_secs, tbl_size, max_ev, n_part_max)
    call pc%set_coll_callback(coll_callback)
    call pc%set_outside_check(outside_check)

    PM_binner%n_bins = PD_size(1)
    PM_binner%inv_dx = 1/PD_dr(1)
  end subroutine PM_initialize

  subroutine coll_callback(my_part, c_ix, c_type)
    use m_cross_sec
    use m_efield_amr
    type(PC_part_t), intent(in) :: my_part
    integer, intent(in) :: c_ix, c_type
    select case (c_type)
    case (CS_ionize_t)
       call E_add_to_var(E_i_pion, my_part%x, my_part%w)
       ! This is for photoionization
       call E_add_to_var(E_i_exc, my_part%x, my_part%w)
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
    real(dp) :: fld(3)
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
          n = int(rng%uni_01() * pc%n_part) + 1
          pos_samples(:,i) = pc%particles(n)%x
          fld_samples(:,i) = E_get_field(pos_samples(:,i))
       end do
    else
       do i = 1, size(pos_samples, 2)
          fld = E_get_field(pos_samples(:,i))
          fld_err  = max(fld_err, norm2(fld - fld_samples(:,i)))
       end do
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

  function PM_get_max_dt(pc, myrank, root) result(dt_max)
    type(PC_t), intent(in) :: pc
    integer, intent(in) :: myrank, root
    real(dp) :: dt_max
    ! Do something with cfl..
    print *, "TODO CFL"
    dt_max = 1.0e-12_dp
  end function PM_get_max_dt

  function PM_bin_func(binner, my_part) result(i_bin)
    class(bin_t), intent(in)     :: binner
    type(PC_part_t), intent(in) :: my_part

    integer                     :: i_bin
    real(dp)                    :: x, x_min, inv_dx
    x      = my_part%x(1)
    i_bin  = 1 + int((x - binner%x_min) * binner%inv_dx)
    i_bin  = max(1, i_bin)
    i_bin  = min(binner%n_bins, i_bin)
  end function PM_bin_func
end module m_particle
