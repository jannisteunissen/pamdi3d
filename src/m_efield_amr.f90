module m_efield_amr

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  type amr_grid_t
     integer                   :: Nr(3)
     integer                   :: lvl
     real(dp)                  :: dr(3), inv_dr(3)
     real(dp)                  :: r_min(3), r_max(3)
     real(dp), allocatable     :: vars(:,:,:,:)
     integer, allocatable      :: child_ix(:,:,:)
     type(amr_grid_t), pointer :: parent
     type(amr_grid_t), pointer :: children(:)
     integer                   :: n_child
  end type amr_grid_t

  type amr_grid_p
     type(amr_grid_t), pointer :: ptr
  end type amr_grid_p

  type box_t
     integer :: i_min(3)
     integer :: i_max(3)
  end type box_t

  integer, parameter :: E_i_pot  = 1, E_i_pion = 2, E_i_elec = 3, &
       E_i_nion = 4, E_i_O2m  = 5, E_i_src  = 6, &
       E_i_Ex   = 7, E_i_Ey   = 8, E_i_Ez   = 9, E_i_tmp = 10, E_n_vars = 10
  character(len=*), parameter :: E_var_names(*) = (/"pot ", "pion", "elec", &
       "nion", "O2m ", "src ", "Ex  ", "Ey  ", "Ez  ", "E   "/)
  character(len=*), parameter :: E_var_units(*) = (/"V   ", "1_m3", "1_m3", &
       "1_m3", "1_m3", "1_m3", "V_m ", "V_m ", "V_m ", "V_m "/)

  integer  :: E_min_grid_size
  integer  :: E_bc_type
  integer  :: E_min_grid_separation
  integer  :: E_grid_buffer_width

  integer  :: E_ref_max_lvl
  integer  :: E_ref_lvl_electrode
  real(dp) :: E_ref_min_xyz(3)
  real(dp) :: E_ref_max_xyz(3)
  real(dp) :: E_ref_min_elec_dens
  real(dp) :: E_prev_max_diff

  !> Electric field values
  real(dp), allocatable :: E_efield_values(:)

  !> Electric field times (at which the values given above are used)
  real(dp), allocatable :: E_efield_times(:)

  !> List of the max. dr vs. electric field
  real(dp), allocatable :: E_dr_vs_efield(:,:)

  !> The root grid
  type(amr_grid_t), pointer :: root_grid

  public :: amr_grid_t
  public :: E_i_elec, E_i_pion, E_i_nion, E_i_O2m

  public :: E_initialize
  public :: E_benchmark
  public :: E_compute_field
  public :: E_write_grids
  public :: E_update_grids
  public :: E_get_vars
  public :: E_get_var
  public :: E_get_var_x_dr
  public :: E_get_field
  public :: E_get_accel_part
  public :: E_get_dr
  public :: E_get_smallest_dr
  public :: E_get_max_of_vars
  public :: E_collect_mpi
  public :: E_loop_over_grids
  public :: E_add_to_var
  public :: E_xyz_to_ix
  public :: E_ix_to_xyz
  public :: E_share_vars
  public :: E_set_vars
  public :: E_adjust_electrode
  public :: E_readjust_elec_new_grid

contains

  subroutine E_initialize(cfg)
    use m_config
    use m_phys_domain
    type(CFG_t), intent(in) :: cfg
    integer                 :: dyn_size, ref_size

    call CFG_get(cfg, "ref_min_grid_size", E_min_grid_size)
    call CFG_get(cfg, "ref_min_grid_separation", E_min_grid_separation)
    call CFG_get(cfg, "ref_min_lvl_electrode", E_ref_lvl_electrode)

    call CFG_get_size(cfg, "sim_efield_times", dyn_size)
    allocate(E_efield_times(dyn_size))
    allocate(E_efield_values(dyn_size))
    call CFG_get(cfg, "sim_efield_values", E_efield_values)
    call CFG_get(cfg, "sim_efield_times", E_efield_times)

    call CFG_get(cfg, "elec_bc_type", E_bc_type)
    call CFG_get_size(cfg, "ref_delta_values", dyn_size)
    call CFG_get_size(cfg, "ref_max_efield_at_delta", ref_size)

    if (dyn_size /= ref_size) then
       print *, "ref_efield_values and ref_delta_values have unequal size"
       stop
    end if

    allocate( E_dr_vs_efield(2, dyn_size) )
    call CFG_get(cfg, "ref_delta_values", E_dr_vs_efield(1, :))
    call CFG_get(cfg, "ref_max_efield_at_delta", E_dr_vs_efield(2, :))

    call CFG_get(cfg, "ref_min_elec_dens", E_ref_min_elec_dens)
    call CFG_get(cfg, "ref_max_levels", E_ref_max_lvl)
    call CFG_get(cfg, "ref_buffer_width", E_grid_buffer_width)
    call CFG_get(cfg, "grid_plasma_min_rel_pos", E_ref_min_xyz)
    call CFG_get(cfg, "grid_plasma_max_rel_pos", E_ref_max_xyz)
    E_ref_min_xyz = E_ref_min_xyz * PD_r_max
    E_ref_max_xyz = E_ref_max_xyz * PD_r_max

    allocate(root_grid)
    call set_grid(root_grid, 0, PD_size, (/0.0_dp, 0.0_dp, 0.0_dp/), PD_r_max)
  end subroutine E_initialize

  subroutine E_benchmark()
    integer :: ix, n_runs
    integer, parameter :: n_runs_max = 10, n_points_max = 1000*1000
    real(dp) :: time_start, time_mid, time_end, accum
    real(dp) :: coords(3, n_points_max), xyz(3)

    do ix = 1, n_points_max
       xyz = (ix-1.0_dp) / (n_points_max)
       xyz = xyz * (root_grid%r_max - root_grid%r_min) + root_grid%r_min
       coords(:, ix) = xyz
    end do
    accum = 0
    call E_set_vars((/E_i_tmp/), (/0.0_dp/))

    call cpu_time(time_start)
    do n_runs = 1, n_runs_max
       do ix = 1, n_points_max
          call E_add_to_var(E_i_tmp, coords(:, ix), 1.0_dp)
       end do
    end do
    call cpu_time(time_mid)
    do n_runs = 1, n_runs_max
       do ix = 1, n_points_max
          accum = accum + E_get_var(E_i_tmp, coords(:, ix))
       end do
    end do
    call cpu_time(time_end)

    print *, "sum", accum
    print *, "add", time_mid - time_start
    print *, "lookup", time_end - time_mid
    print *, "ns per point", (time_end - time_start) / (1.0e-9_dp * n_runs_max * n_points_max)
  end subroutine E_benchmark

  subroutine set_grid(amr_grid, lvl, Nr, r_min, r_max, children, parent)
    type(amr_grid_t), intent(inout)    :: amr_grid
    integer, intent(in)                :: lvl, Nr(3)
    real(dp), intent(in)               :: r_min(3), r_max(3)
    type(amr_grid_t), target, optional :: children(:), parent
    integer :: nc, i_min(3), i_max(3)

    amr_grid%lvl      = lvl
    amr_grid%Nr       = Nr
    amr_grid%r_min    = r_min
    amr_grid%r_max    = r_max
    amr_grid%dr       = (r_max - r_min) / (Nr - 1)
    amr_grid%inv_dr   = 1 / amr_grid%dr

    ! This sets all vars and child_ix to 0
    call allocate_grid_data(amr_grid)

    if (present(parent)) then
       amr_grid%parent   => parent
    else
       nullify(amr_grid%parent)
    end if

    if (present(children)) then
       amr_grid%children => children
       amr_grid%n_child  = size(children)

       ! Set child_ix
       do nc = 1, size(children)
          i_min = E_xyz_to_ix(amr_grid, children(nc)%r_min)
          i_max = E_xyz_to_ix(amr_grid, children(nc)%r_max)

          ! Mark child region
          amr_grid%child_ix(i_min(1):i_max(1), i_min(2):i_max(2), &
               i_min(3):i_max(3)) = nc
       end do
    else
       nullify(amr_grid%children)
       amr_grid%n_child  = 0
    end if
  end subroutine set_grid

  subroutine E_adjust_electrode(myrank, root, time)
    use m_electrode
    use m_units_constants
    integer, intent(in)   :: myrank, root
    real(dp), intent(in)  :: time
    integer               :: ix, nElecPoints
    real(dp)              :: xyz(3), max_diff_weight
    real(dp), allocatable :: diff_voltage(:)

    if (myrank == root) then
       print *, "Adjusting electrode"
       ! Store old potential, put it back at the end
       call copy_var_to_var_recursive(root_grid, E_i_pot, E_i_tmp)
       nElecPoints = EL_getNPoints()
       allocate(diff_voltage(nElecPoints))

       diff_voltage(:) = 0.0_dp

       ! Now compute the potential at electrode points when they have unit charge
       call E_set_vars((/E_i_src/), (/0.0_dp/))

       do ix = 1, nElecPoints
          call EL_getSurfacePoint(ix, xyz)
          call E_add_to_var(E_i_src, xyz, UC_elec_q_over_eps0)
       end do

       call compute_potential_recursive(root_grid, time, .true.)
       max_diff_weight = 0.0_dp

       do ix = 1, nElecPoints
          call EL_getSurfacePoint(ix, xyz)
          diff_voltage(ix) = E_get_var(E_i_pot, xyz) - diff_voltage(ix)
          if (abs(1 - EL_getWeightFactor(ix) * diff_voltage(ix)) > max_diff_weight) then
             max_diff_weight = abs(1 - EL_getWeightFactor(ix) * diff_voltage(ix))
          end if
          call EL_setWeightFactor(ix, 1.0d0 / diff_voltage(ix))
       end do
       call copy_var_to_var_recursive(root_grid, E_i_tmp, E_i_pot)
       ! print *, "Max change in weightfactor", max_diff_weight
    end if
  end subroutine E_adjust_electrode

  recursive subroutine sync_grid_structure_recursive(amr_grid, myrank, root)
    use mpi
    type(amr_grid_t), intent(inout), target :: amr_grid
    integer, intent(in)                     :: myrank, root
    integer                                 :: nc, ierr, n_points

    call MPI_BCAST(amr_grid%lvl, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(amr_grid%Nr, 3, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(amr_grid%n_child, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(amr_grid%dr, 3, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(amr_grid%inv_dr, 3, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(amr_grid%r_min, 3, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(amr_grid%r_max, 3, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)

    if (myrank /= root) then
       call allocate_grid_data(amr_grid)
       allocate(amr_grid%children(amr_grid%n_child))
    end if

    n_points = product(amr_grid%Nr)
    call MPI_BCAST(amr_grid%child_ix, n_points, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)

    do nc = 1, amr_grid%n_child
       amr_grid%children(nc)%parent => amr_grid
       call sync_grid_structure_recursive(amr_grid%children(nc), myrank, root)
    end do
  end subroutine sync_grid_structure_recursive

  subroutine E_readjust_elec_new_grid(myrank, root, time, cntr)
    use m_electrode
    use m_units_constants
    use mpi
    real(dp), intent(in) :: time
    integer, intent(in)  :: myrank, root
    integer, intent(out) :: cntr
    integer              :: ix, nElecPoints, ierror
    real(dp)             :: xyz(3), voltage, voltage_diff, newCharge, max_diff

    call mpi_collect_recursive(root_grid, (/E_i_elec, E_i_pion, E_i_nion, E_i_O2m/), myrank, root)
    call E_adjust_electrode(myrank, root, time)

    if (myrank == root) then
       nElecPoints = EL_getNPoints()
       call set_source_term_recursive(root_grid)

       do ix = 1, nElecPoints
          call EL_getSurfacePoint(ix, xyz)
          call E_add_to_var(E_i_src, xyz, EL_getCharge(ix))
       end do

       print *, "Readjusting electrode"
       voltage = EL_getVoltage(time)
       cntr    = 0

       do
          cntr = cntr + 1
          call compute_potential_recursive(root_grid, time)
          call check_elec_voltage(time, max_diff)
          if (max_diff < 1.05_dp * E_prev_max_diff) exit

          call set_source_term_recursive(root_grid)

          do ix = 1, nElecPoints
             call EL_getSurfacePoint(ix, xyz)
             voltage_diff = voltage - E_get_var(E_i_pot, xyz)
             newCharge    = EL_getCharge(ix) + voltage_diff * &
                  EL_getWeightFactor(ix) * UC_elec_q_over_eps0
             call EL_setCharge(ix, newCharge)

             call E_add_to_var(E_i_src, xyz, newCharge)
          end do
       end do
       print *, "Readjusting electrode done", max_diff, cntr
    end if
    call MPI_BCAST(cntr, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierror)
  end subroutine E_readjust_elec_new_grid

  subroutine check_elec_voltage(time, max_diff)
    use m_electrode
    real(dp), intent(in) :: time
    real(dp), intent(out) :: max_diff
    integer :: ix, nElecPoints
    real(dp) :: xyz(3), temp, voltage, elec_voltage

    voltage = 0.0d0
    max_diff = 0.0d0
    elec_voltage = EL_getVoltage(time)
    nElecPoints = EL_getNPoints()

    do ix = 1, nElecPoints
       call EL_getSurfacePoint(ix, xyz)
       temp = E_get_var(E_i_pot, xyz)
       voltage = voltage + temp
       if (abs(temp - elec_voltage) > max_diff) then
          max_diff = abs(temp - EL_getVoltage(time))
       end if
    end do
    print *, "Electrode: avg ", voltage / nElecPoints, "max diff", max_diff
  end subroutine check_elec_voltage

  subroutine E_compute_field(myrank, root, time)
    use m_electrode
    use m_units_constants
    use m_phys_domain
    real(dp), intent(in) :: time
    integer, intent(in)  :: myrank, root
    integer              :: ix, nElecPoints
    real(dp)             :: xyz(3), voltage_diff, newCharge, voltage

    ! Gather the charge densities
    call mpi_collect_recursive(root_grid, (/E_i_elec, E_i_pion, E_i_nion, E_i_O2m/), myrank, root)

    if (myrank == root) then
       call set_source_term_recursive(root_grid)

       ! Add the electrode charges based on the old potential
       if (PD_use_elec) then
          voltage = EL_getVoltage(time)
          nElecPoints = EL_getNPoints()

          do ix = 1, nElecPoints
             call EL_getSurfacePoint(ix, xyz)
             voltage_diff = voltage - E_get_var(E_i_pot, xyz)
             newCharge    = EL_getCharge(ix) + voltage_diff * EL_getWeightFactor(ix) * UC_elec_q_over_eps0
             call EL_setCharge(ix, newCharge)
             call E_add_to_var(E_i_src, xyz, newCharge)
          end do
       end if

       call compute_potential_recursive(root_grid, time)
       call compute_field_recursive(root_grid)

       if (PD_use_elec) call check_elec_voltage(time, E_prev_max_diff)
    end if

    call E_share_vars((/E_i_Ex, E_i_Ey, E_i_Ez/), root)
  end subroutine E_compute_field

  recursive subroutine compute_potential_recursive(amr_grid, time, zero_bc)
    type(amr_grid_t), intent(inout) :: amr_grid
    real(dp), intent(in)            :: time
    logical, intent(in), optional   :: zero_bc

    integer                         :: req_workspace
    integer                         :: ierror, Nx, Ny, Nz, nc
    real(dp)                        :: dummy(1)
    real(dp), save, allocatable     :: workspace(:)
    integer, save                   :: size_workspace = 0
    logical                         :: use_homogenuous_dirichlet

    interface
       subroutine hw3crt(xMin, xMax, xPanels, xBDCND, xBDCmin, xBDCmax, &
            yMin, yMax, yPanels, yBDCND, yBDCmin, yBDCmax, &
            zMin, zMax, zPanels, zBDCND, zBDCmin, zBDCmax, &
            lambda, sourceFirstDim, sourceSecondDim, source, &
            pertrb, ierror, workspace)
         import dp
         integer  :: xPanels, xBDCND, yPanels, yBDCND, zPanels, zBDCND
         integer  :: sourceFirstDim, sourceSecondDim, ierror
         real(dp) :: xMin, xMax, xBDCmin(*), xBDCmax(*)
         real(dp) :: yMin, yMax, yBDCmin(*), yBDCmax(*)
         real(dp) :: zMin, zMax, zBDCmin(*), zBDCmax(*)
         real(dp) :: lambda, source(sourceFirstDim, sourceSecondDim, *)
         real(dp) :: pertrb, workspace(*)
       end subroutine hw3crt
    end interface

    Nx = amr_grid%Nr(1)
    Ny = amr_grid%Nr(2)
    Nz = amr_grid%Nr(3)

    req_workspace = 30 + Nx + Ny + 5*Nz + max(Nx, Ny, Nz) + 7*(Nx/2 + Ny/2)

    if (req_workspace > size_workspace) then
       if (allocated(workspace)) deallocate(workspace)
       allocate(workspace(req_workspace))
       size_workspace = req_workspace
    end if

    use_homogenuous_dirichlet = .false.
    if (present(zero_bc)) use_homogenuous_dirichlet = zero_bc

    if (use_homogenuous_dirichlet) then
       amr_grid%vars(:,:,:, E_i_pot) = 0.0_dp
    else if (associated(amr_grid%parent)) then
       call set_bc_from_parent(amr_grid, (/E_i_pot/))
    else
       call set_bc_dirichlet(amr_grid, time)
    end if

    ! The input array holds the boundary condition on the sides and the source
    ! term in its interior
    amr_grid%vars(2:Nx-1, 2:Ny-1, 2:Nz-1, E_i_pot) = &
         amr_grid%vars(2:Nx-1, 2:Ny-1, 2:Nz-1, E_i_src)

    call hw3crt(amr_grid%r_min(1), amr_grid%r_max(1), Nx-1, 1, dummy, dummy, &
         amr_grid%r_min(2), amr_grid%r_max(2), Ny-1, 1, dummy, dummy, &
         amr_grid%r_min(3), amr_grid%r_max(3), Nz-1, 1, dummy, dummy, &
         0.0d0, Nx, Ny, amr_grid%vars(:,:,:, E_i_pot), dummy(1), &
         ierror, workspace)

    if (ierror /= 0) then
       print *, "HW3CRT: ierror = ", ierror
       stop
    end if

    ! Call this routine for children
    do nc = 1, amr_grid%n_child
       call compute_potential_recursive(amr_grid%children(nc), time)
    end do
  end subroutine compute_potential_recursive

  subroutine E_collect_mpi(v_ixs, myrank, root)
    integer, intent(in) :: v_ixs(:), myrank, root
    call mpi_collect_recursive(root_grid, v_ixs, myrank, root)
  end subroutine E_collect_mpi

  recursive subroutine mpi_collect_recursive(amr_grid, v_ixs, myrank, root)
    use mpi
    type(amr_grid_t), intent(inout) :: amr_grid
    integer, intent(in)             :: v_ixs(:), myrank, root
    integer                         :: ierr, n_points, ix, nc

    n_points = product(amr_grid%Nr)
    do ix = 1, size(v_ixs)
       if (myrank == root) then
          call MPI_REDUCE(MPI_IN_PLACE, amr_grid%vars(:, :, :, v_ixs(ix)), n_points, &
               MPI_DOUBLE_PRECISION, MPI_SUM, root, MPI_COMM_WORLD, ierr)
       else
          call MPI_REDUCE(amr_grid%vars(:, :, :, v_ixs(ix)), amr_grid%vars(:, :, :, v_ixs(ix)), n_points, &
               MPI_DOUBLE_PRECISION, MPI_SUM, root, MPI_COMM_WORLD, ierr)
          amr_grid%vars(:, :, :, v_ixs(ix)) = 0.0_dp
       end if
    end do

    do nc = 1, amr_grid%n_child
       call mpi_collect_recursive(amr_grid%children(nc), v_ixs, myrank, root)
    end do
  end subroutine mpi_collect_recursive

  subroutine E_share_vars(v_ixs, root)
    integer, intent(in) :: v_ixs(:), root
    call mpi_share_recursive(root_grid, v_ixs, root)
  end subroutine E_share_vars

  recursive subroutine mpi_share_recursive(amr_grid, v_ixs, root)
    use mpi
    type(amr_grid_t), intent(inout) :: amr_grid
    integer, intent(in)             :: v_ixs(:), root
    integer                         :: ierr, n_points, ix, nc

    n_points = product(amr_grid%Nr)
    do ix = 1, size(v_ixs)
       call MPI_BCAST(amr_grid%vars(:, :, :, v_ixs(ix)), n_points, &
            MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
    end do

    do nc = 1, amr_grid%n_child
       call mpi_share_recursive(amr_grid%children(nc), v_ixs, root)
    end do
  end subroutine mpi_share_recursive

  subroutine E_loop_over_grids(pc, rng, func_p)
    use m_particle_core
    use m_random
    class(PC_t), intent(inout) :: pc
    type(RNG_t), intent(inout) :: rng

    interface
       subroutine func_p(pc, rng, ag)
         import
         class(PC_t), intent(inout) :: pc
         type(RNG_t), intent(inout) :: rng
         type(amr_grid_t), intent(inout) :: ag
       end subroutine func_p
    end interface

    integer                       :: ix
    type(amr_grid_p), allocatable :: grid_list(:)

    call create_grid_list(root_grid, grid_list)
    do ix = 1, size(grid_list)
       call func_p(pc, rng, grid_list(ix)%ptr)
    end do
  end subroutine E_loop_over_grids

  ! After this routine the electric field has to be recomputed!
  subroutine E_update_grids(myrank, root, time)
    integer, intent(in)       :: myrank, root
    real(dp), intent(in)      :: time
    type(amr_grid_t), pointer :: new_grid

    ! print *, myrank, "updating grids"
    call mpi_collect_recursive(root_grid, &
         (/E_i_pion, E_i_elec, E_i_nion, E_i_O2m/), myrank, root)

    if (myrank == root) then
       allocate(new_grid)
       new_grid         = root_grid
       new_grid%n_child = 0
       nullify(new_grid%children)

       ! Boundary values are not correctly set in a PIC simulation, so get from parent
       call set_bc_from_parent_recursive(root_grid, &
            (/E_i_pion, E_i_elec, E_i_nion, E_i_O2m/))
       call set_new_grids_recursive(new_grid)
       call recursive_remove_amr_grid(root_grid)
       deallocate(root_grid)
       root_grid => new_grid
    else
       call recursive_remove_amr_grid(root_grid)
    end if

    call sync_grid_structure_recursive(root_grid, myrank, root)
  end subroutine E_update_grids

  logical function need_to_refine(amr_grid, ix)
    use m_electrode
    use m_phys_domain
    use m_lookup_table
    type(amr_grid_t), intent(in) :: amr_grid
    integer, intent(in)          :: ix(3)
    real(dp)                     :: xyz(3)
    real(dp)                     :: max_ef

    xyz = E_ix_to_xyz(amr_grid, ix)

    if (amr_grid%lvl == E_ref_max_lvl .or. any(xyz < E_ref_min_xyz .or. xyz > E_ref_max_xyz)) then
       need_to_refine = .false. ! Outside this region never refine
    else if (amr_grid%vars(ix(1), ix(2), ix(3), E_i_elec) < E_ref_min_elec_dens) then
       need_to_refine = .false.
    else
       ! Find maximum Efield for the delta at this level
       call LT_lin_interp_list(E_dr_vs_efield(1, :), E_dr_vs_efield(2, :),  maxval(amr_grid%dr), max_ef)
       if (norm2(E_get_field(xyz)) > max_ef) then
          need_to_refine = .true.
       else
          need_to_refine = .false.
       end if
    end if

    if (PD_use_elec) then
       if (EL_around_elec_tip(xyz)) then
          ! Refine around electrode tip
          need_to_refine = need_to_refine .or. amr_grid%lvl < E_ref_lvl_electrode
       end if
    end if
  end function need_to_refine

  recursive subroutine set_new_grids_recursive(amr_grid)
    type(amr_grid_t), intent(inout), target :: amr_grid

    integer                                 :: i, j, k, o_ix(3), ix_offset(3)
    integer                                 :: n_ref, n_child, nc, i_min(3), i_max(3)
    real(dp)                                :: xyz(3)
    type(amr_grid_t), pointer               :: old_grid
    logical                                 :: found_grid
    integer, allocatable                    :: ref_ixs(:, :)
    type(box_t), allocatable                :: box_list(:)
    integer, parameter                      :: v_ixs(*) = &
         (/E_i_pion, E_i_nion, E_i_O2m, E_i_elec/)

    if (associated(amr_grid%children) .or. amr_grid%n_child > 0) then
       print *, "set_new_grids_recursive error: children already present!"
       stop
    end if

    allocate(ref_ixs(3, product(amr_grid%Nr))) ! Allocate max amount needed
    n_ref      = 0
    found_grid = .false.

    do k = 1, amr_grid%Nr(3)
       do j = 1, amr_grid%Nr(2)
          do i = 1, amr_grid%Nr(1)
             xyz = E_ix_to_xyz(amr_grid, (/i, j, k/))

             ! Fill grid values. If possible, directly copy from old grid.
             if (found_grid) then
                o_ix = (/i, j, k/) + ix_offset
                if (any(o_ix < 1 .or. o_ix > old_grid%Nr)) then
                   found_grid = .false.
                else
                   amr_grid%vars(i,j,k, v_ixs) = old_grid%vars(o_ix(1), o_ix(2), o_ix(3), v_ixs)
                end if
             end if

             if (.not. found_grid) then
                old_grid => root_grid
                call get_grid_at_maxlvl(old_grid, xyz, amr_grid%lvl)
                if (old_grid%lvl == amr_grid%lvl) then
                   found_grid = .true.
                   ix_offset = nint((amr_grid%r_min - old_grid%r_min) * amr_grid%inv_dr)
                end if
                amr_grid%vars(i,j,k, v_ixs) = get_vars_at_grid(old_grid, v_ixs, xyz)
             end if

             ! Check whether we need to refine
             if (need_to_refine(amr_grid, (/i, j, k/))) then
                n_ref = n_ref + 1
                ref_ixs(:, n_ref) = (/i, j, k/)
             end if
          end do
       end do
    end do

    nullify(old_grid)

    ! By default set child region to 0 (which means no children)
    amr_grid%child_ix = 0

    if (n_ref > 0) then
       call get_boxes(ref_ixs(:, 1:n_ref), box_list, E_min_grid_separation**2, &
            (/2, 2, 2/), amr_grid%Nr-1, E_grid_buffer_width, E_min_grid_size)
       n_child          = size(box_list)
       amr_grid%n_child = n_child
       allocate(amr_grid%children(n_child))
       ! print *, "Refining, children", n_child, n_ref

       do nc = 1, n_child
          i_min = box_list(nc)%i_min
          i_max = box_list(nc)%i_max

          ! Mark child region
          amr_grid%child_ix(i_min(1):i_max(1), i_min(2):i_max(2), &
               i_min(3):i_max(3)) = nc

          ! Create child, has 2 * (i_max-i_min) + 1 points
          call set_grid(amr_grid%children(nc), amr_grid%lvl+1, &
               2 * (i_max-i_min) + 1, E_ix_to_xyz(amr_grid, i_min), &
               E_ix_to_xyz(amr_grid, i_max), parent=amr_grid)
          call set_new_grids_recursive(amr_grid%children(nc))
       end do
    end if
  end subroutine set_new_grids_recursive

  ! Given a list of points (point_list), returns a list of boxes (out_boxes)
  ! that enclose all the points. The squared distance between these boxes is at
  ! least min_dist2. The boxes furthermore are placed between i_min and i_max,
  ! have a buffer_width around the points they enclose and have minimal size
  ! min_size.
  subroutine get_boxes(point_list, out_boxes, min_dist2, i_min, i_max, &
       buffer_width, min_size)
    integer, intent(in)      :: point_list(:,:), min_dist2, i_min(3), i_max(3), buffer_width, min_size
    type(box_t), allocatable :: out_boxes(:)

    type(box_t), allocatable :: boxes(:), new_boxes(:)
    integer                  :: ix, nn, n_points, n_boxes, n_boxes_new, n_boxes_prev
    integer                  :: curr_size(3), add_size(3)
    logical                  :: add_box, already_adjusted

    n_points = size(point_list, 2)
    allocate(boxes(n_points))
    allocate(new_boxes(n_points))

    do ix = 1, n_points
       boxes(ix)%i_min = point_list(:, ix)
       boxes(ix)%i_max = point_list(:, ix)
       where (boxes(ix)%i_min < i_min) boxes(ix)%i_min = i_min
       where (boxes(ix)%i_max > i_max) boxes(ix)%i_max = i_max
    end do

    already_adjusted = .false.
    n_boxes          = n_points

    do
       n_boxes_prev = n_boxes
       n_boxes_new  = 0

       ! Replace boxes by a reduced list of boxes
       do ix = 1, n_boxes
          add_box = .true.

          do nn = n_boxes_new, 1, -1
             ! Reverse order, since latest new box is probably closer
             ! If close enough to other box, merge
             if (distance2_boxes(boxes(ix), new_boxes(nn)) < min_dist2) then
                new_boxes(nn) = merge_boxes(boxes(ix), new_boxes(nn))
                add_box = .false.
                exit
             end if
          end do

          ! If not merged with other box, add to new list of boxes
          if (add_box) then
             n_boxes_new = n_boxes_new + 1
             new_boxes(n_boxes_new) = boxes(ix)
          end if
       end do

       n_boxes          = n_boxes_new
       boxes(1:n_boxes) = new_boxes(1:n_boxes)

       if (n_boxes == n_boxes_prev) then
          if (already_adjusted) then
             exit
          else
             ! Make sure the boxes are large enough
             do ix = 1, n_boxes
                boxes(ix)%i_min = boxes(ix)%i_min - buffer_width
                boxes(ix)%i_max = boxes(ix)%i_max + buffer_width

                ! But they should not go past i_min or i_max
                do
                   where (boxes(ix)%i_min < i_min) boxes(ix)%i_min = i_min
                   where (boxes(ix)%i_max > i_max) boxes(ix)%i_max = i_max
                   curr_size = 1 + boxes(ix)%i_max - boxes(ix)%i_min
                   add_size = 1 + min_size - curr_size
                   if (all(add_size <= 0)) exit

                   where (add_size > 0)
                      boxes(ix)%i_min = boxes(ix)%i_min - 1
                      boxes(ix)%i_max = boxes(ix)%i_max + 1
                   end where
                end do
             end do
             already_adjusted = .true.
          end if
       end if

       ! If n_boxes /= n_boxes_prev, repeat the outer loop
    end do

    allocate(out_boxes(n_boxes))
    out_boxes = boxes(1:n_boxes)
  end subroutine get_boxes

  integer function distance2_boxes(box_a, box_b)
    type(box_t), intent(in) :: box_a, box_b
    integer                 :: dist(3)

    dist = 0
    where (box_a%i_max < box_b%i_min) dist = box_b%i_min - box_a%i_max
    where (box_a%i_min > box_b%i_max) dist = box_a%i_min - box_b%i_max
    distance2_boxes = sum(dist**2)
  end function distance2_boxes

  type(box_t) function merge_boxes(box_a, box_b)
    type(box_t), intent(in) :: box_a, box_b
    where (box_a%i_min < box_b%i_min)
       merge_boxes%i_min = box_a%i_min
    elsewhere
       merge_boxes%i_min = box_b%i_min
    end where

    where (box_a%i_max > box_b%i_max)
       merge_boxes%i_max = box_a%i_max
    elsewhere
       merge_boxes%i_max = box_b%i_max
    end where
  end function merge_boxes

  function E_ix_to_xyz(amr_grid, ix) result(xyz)
    type(amr_grid_t), intent(in) :: amr_grid
    integer, intent(in)          :: ix(3)
    real(dp)                     :: xyz(3)
    xyz = amr_grid%r_min + (ix-1) * amr_grid%dr
  end function E_ix_to_xyz

  function E_xyz_to_ix(amr_grid, xyz) result(ix)
    type(amr_grid_t), intent(in) :: amr_grid
    real(dp), intent(in)         :: xyz(3)
    integer                      :: ix(3)
    ix = nint((xyz - amr_grid%r_min) * amr_grid%inv_dr + 1) ! Nearest index
  end function E_xyz_to_ix

  subroutine set_bc_dirichlet(amr_grid, time)
    use m_electrode
    use m_phys_domain
    use m_lookup_table
    type(amr_grid_t), intent(inout) :: amr_grid
    real(dp), intent(in)            :: time

    integer                         :: i, j, k
    integer                         :: Nx, Ny, Nz
    real(dp)                        :: efield, voltage, temp, xyz(3)
    real(dp)                        :: elec_pos(3), decay_range, elec_radius

    Nx = amr_grid%Nr(1)
    Ny = amr_grid%Nr(2)
    Nz = amr_grid%Nr(3)

    if (PD_use_elec) then
       voltage = EL_getVoltage(time)
    else
       call LT_lin_interp_list(E_efield_times, E_efield_values, time, efield)
       voltage = -efield * (amr_grid%r_max(3) - amr_grid%r_min(3))
    end if

    if (E_bc_type == 1) then
       ! Potential linearly increasing along the sides, from voltage to zero
       do k = 1, Nz
          do j = 1, Ny
             do i = 1, Nx, Nx-1
                amr_grid%vars(i,j,k, E_i_pot) = (1.0d0 - (k-1)/dble(Nz-1)) * voltage
             end do
          end do
       end do

       do k = 1, Nz
          do j = 1, Ny, Ny-1
             do i = 1, Nx
                amr_grid%vars(i,j,k, E_i_pot) = (1.0d0 - (k-1)/dble(Nz-1)) * voltage
             end do
          end do
       end do

       amr_grid%vars(:,:,1, E_i_pot) = voltage
       amr_grid%vars(:,:,Nz, E_i_pot) = 0

    else if (E_bc_type == 2) then
       ! Potential zero at the sides, for use with an electrode
       if (.not. PD_use_elec) then
          print *, "This bnd. cnd. needs an electrode:", E_bc_type
          stop
       end if

       call EL_getBottomPos(elec_pos)
       elec_radius = EL_getRadius(amr_grid%r_min(3))
       decay_range = 3 * elec_radius

       amr_grid%vars(1,:,:, E_i_pot)  = 0.0d0
       amr_grid%vars(Nx,:,:, E_i_pot) = 0.0d0
       amr_grid%vars(:,1,:, E_i_pot)  = 0.0d0
       amr_grid%vars(:,Ny,:, E_i_pot) = 0.0d0
       amr_grid%vars(:,:,Nz, E_i_pot) = 0.0d0

       do j = 1, Ny
          do i = 1, Nx
             xyz = E_ix_to_xyz(amr_grid, (/i, j, 1/))
             temp = norm2(xyz(1:2) - elec_pos(1:2))

             if (temp <= elec_radius) then
                amr_grid%vars(i,j,1, E_i_pot) = voltage
             else
                amr_grid%vars(i,j,1, E_i_pot) = voltage * &
                     max(0.0_dp, 1 - (temp-elec_radius)/decay_range)
             end if
          end do
       end do

    else
       print *, "set_bc_dirichlet error, no correct boundary condition used"
       stop
    end if
  end subroutine set_bc_dirichlet

  recursive subroutine set_bc_from_parent_recursive(amr_grid, v_ixs)
    type(amr_grid_t), intent(inout) :: amr_grid
    integer, intent(in)             :: v_ixs(:)
    integer                         :: nc

    if (associated(amr_grid%parent)) call set_bc_from_parent(amr_grid, v_ixs)
    do nc = 1, amr_grid%n_child
       call set_bc_from_parent_recursive(amr_grid%children(nc), v_ixs)
    end do
  end subroutine set_bc_from_parent_recursive

  subroutine set_bc_from_parent(amr_grid, v_ixs)
    type(amr_grid_t), intent(inout) :: amr_grid
    integer, intent(in)             :: v_ixs(:)
    integer                         :: i, j, k, Nx, Ny, Nz
    real(dp)                        :: xyz(3)

    Nx       = amr_grid%Nr(1)
    Ny       = amr_grid%Nr(2)
    Nz       = amr_grid%Nr(3)

    do k = 1, Nz, Nz-1
       do j = 1, Ny
          do i = 1, Nx
             xyz = E_ix_to_xyz(amr_grid, (/i, j, k/))
             amr_grid%vars(i,j,k, v_ixs) = get_vars_at_grid(amr_grid%parent, &
                  v_ixs, xyz)
          end do
       end do
    end do

    do k = 1, Nz
       do j = 1, Ny, Ny-1
          do i = 1, Nx
             xyz = E_ix_to_xyz(amr_grid, (/i, j, k/))
             amr_grid%vars(i,j,k, v_ixs) = get_vars_at_grid(amr_grid%parent, &
                  v_ixs, xyz)
          end do
       end do
    end do

    do k = 1, Nz
       do j = 1, Ny
          do i = 1, Nx, Nx-1
             xyz = E_ix_to_xyz(amr_grid, (/i, j, k/))
             amr_grid%vars(i,j,k, v_ixs) = get_vars_at_grid(amr_grid%parent, &
                  v_ixs, xyz)
          end do
       end do
    end do
  end subroutine set_bc_from_parent

  subroutine get_grid_at(grid_ptr, xyz)
    real(dp), intent(in)                     :: xyz(3)
    type(amr_grid_t), intent(inout), pointer :: grid_ptr
    integer                                  :: child_ix
    do
       child_ix = get_child_ix_at_xyz(grid_ptr, xyz)
       if (child_ix == 0) exit
       grid_ptr => grid_ptr%children(child_ix)
    end do
  end subroutine get_grid_at

  subroutine get_grid_at_maxlvl(grid_ptr, xyz, maxlvl)
    real(dp), intent(in)                     :: xyz(3)
    type(amr_grid_t), intent(inout), pointer :: grid_ptr
    integer, intent(in)                      :: maxlvl
    integer                                  :: child_ix
    do
       child_ix = get_child_ix_at_xyz(grid_ptr, xyz)
       if (child_ix == 0 .or. grid_ptr%lvl == maxlvl) exit
       grid_ptr => grid_ptr%children(child_ix)
    end do
  end subroutine get_grid_at_maxlvl

  function get_vars(amr_grid, v_ixs, xyz) result(vars)
    type(amr_grid_t), intent(in), target :: amr_grid
    integer, intent(in)                  :: v_ixs(:)
    real(dp), intent(in)                 :: xyz(3)
    real(dp), dimension(size(v_ixs))     :: vars
    type(amr_grid_t), pointer            :: grid_ptr
    grid_ptr => amr_grid
    call get_grid_at(grid_ptr, xyz)
    vars = get_vars_at_grid(grid_ptr, v_ixs, xyz)
  end function get_vars

  function get_vars_x_dr(amr_grid, v_ixs, xyz) result(vars)
    type(amr_grid_t), intent(in), target :: amr_grid
    integer, intent(in)                  :: v_ixs(:)
    real(dp), intent(in)                 :: xyz(3)
    real(dp), dimension(size(v_ixs))     :: vars
    type(amr_grid_t), pointer            :: grid_ptr
    grid_ptr => amr_grid
    call get_grid_at(grid_ptr, xyz)
    vars = get_vars_at_grid(grid_ptr, v_ixs, xyz) * product(grid_ptr%dr)
  end function get_vars_x_dr

  integer function get_child_ix_at_xyz(amr_grid, xyz)
    type(amr_grid_t), intent(in) :: amr_grid
    real(dp), intent(in)         :: xyz(3)
    integer                      :: ixs(3)
    ixs = get_ixs(amr_grid, xyz)
    get_child_ix_at_xyz = minval(amr_grid%child_ix(ixs(1):ixs(1)+1, &
         ixs(2):ixs(2)+1, ixs(3):ixs(3)+1))
  end function get_child_ix_at_xyz

  integer function get_child_ix_at_ixs(amr_grid, ixs)
    type(amr_grid_t), intent(in) :: amr_grid
    integer, intent(in)          :: ixs(3)

    get_child_ix_at_ixs = minval(amr_grid%child_ix(ixs(1):ixs(1)+1, &
         ixs(2):ixs(2)+1, ixs(3):ixs(3)+1))
  end function get_child_ix_at_ixs

  function E_get_vars(v_ixs, xyz) result(vars)
    real(dp), intent(in)             :: xyz(3)
    integer, intent(in)              :: v_ixs(:)
    real(dp), dimension(size(v_ixs)) :: vars
    vars = get_vars(root_grid, v_ixs, xyz)
  end function E_get_vars

  subroutine E_get_max_of_vars(v_ixs, max_values)
    integer, intent(in)           :: v_ixs(:)
    real(dp), intent(out)         :: max_values(:)
    integer :: gg, vv
    type(amr_grid_p), allocatable :: grid_list(:)

    call create_grid_list(root_grid, grid_list)
    max_values = -huge(1.0_dp)

    do gg = 1, size(grid_list)
       do vv = 1, size(v_ixs)
          max_values(vv) = max(max_values(vv), &
               maxval(grid_list(gg)%ptr%vars(:,:,:,v_ixs(vv))))
       end do
    end do
  end subroutine E_get_max_of_vars

  real(dp) function E_get_var(v_ix, xyz)
    real(dp), intent(in) :: xyz(3)
    integer, intent(in)  :: v_ix
    real(dp)             :: temp(1)
    temp = get_vars(root_grid, (/v_ix/), xyz)
    E_get_var = temp(1)
  end function E_get_var

  function E_get_field(xyz) result(field)
    real(dp), intent(in) :: xyz(3)
    real(dp)             :: field(3)
    field = get_vars(root_grid, (/E_i_Ex, E_i_Ey, E_i_Ez/), xyz)
  end function E_get_field

  function E_get_accel_part(my_part) result(accel)
    use m_particle_core
    use m_units_constants
    type(PC_part_t), intent(in) :: my_part
    real(dp)                    :: field(3), accel(3)
    field = get_vars(root_grid, (/E_i_Ex, E_i_Ey, E_i_Ez/), my_part%x)
    accel = field * UC_elec_q_over_m
  end function E_get_accel_part

  real(dp) function E_get_var_x_dr(v_ix, xyz)
    real(dp), intent(in) :: xyz(3)
    integer, intent(in)  :: v_ix
    real(dp)             :: temp(1)
    temp = get_vars_x_dr(root_grid, (/v_ix/), xyz)
    E_get_var_x_dr = temp(1)
  end function E_get_var_x_dr

  function E_get_dr(xyz) result(dr)
    real(dp), intent(in) :: xyz(3)
    real(dp)             :: dr(3)
    dr = get_dr(root_grid, xyz)
  end function E_get_dr

  function E_get_smallest_dr() result(dr)
    real(dp)                      :: dr(3)
    type(amr_grid_p), allocatable :: grid_list(:)
    integer                       :: n

    dr(:) = huge(1.0_dp)
    call create_grid_list(root_grid, grid_list)

    do n = 1, size(grid_list)
       where (grid_list(n)%ptr%dr < dr) dr = grid_list(n)%ptr%dr
    end do
  end function E_get_smallest_dr

  function get_dr(amr_grid, xyz) result(dr)
    type(amr_grid_t), intent(in), target :: amr_grid
    real(dp), intent(in)                 :: xyz(3)
    real(dp)                             :: dr(3)
    type(amr_grid_t), pointer            :: grid_ptr
    grid_ptr => amr_grid
    call get_grid_at(grid_ptr, xyz)
    dr       = grid_ptr%dr
  end function get_dr

  function get_vars_at_grid(amr_grid, v_ixs, xyz) result(vars)
    type(amr_grid_t), intent(in)     :: amr_grid
    real(dp), intent(in)             :: xyz(3)
    integer, intent(in)              :: v_ixs(:)
    real(dp), dimension(size(v_ixs)) :: vars

    integer                          :: ixs(3), i
    real(dp)                         :: weights(2,2,2)

    call get_ixs_and_weights(amr_grid, xyz, ixs, weights)
    do i = 1, size(v_ixs)
       vars(i) = sum(amr_grid%vars(ixs(1):ixs(1)+1, ixs(2):ixs(2)+1, ixs(3):ixs(3)+1, v_ixs(i)) * weights)
    end do
  end function get_vars_at_grid

  recursive subroutine compute_field_recursive(amr_grid)
    type(amr_grid_t), intent(inout) :: amr_grid
    integer                         :: i, j, k, nc

    ! Do the interior part of the domain
    do k = 2, amr_grid%Nr(3)-1
       do j = 2, amr_grid%Nr(2)-1
          do i = 2, amr_grid%Nr(1)-1
             amr_grid%vars(i,j,k, (/E_i_Ex, E_i_Ey, E_i_Ez/)) = 0.5_dp * amr_grid%inv_dr * (/&
                  amr_grid%vars(i-1, j, k, E_i_pot) - amr_grid%vars(i+1, j, k, E_i_pot), &
                  amr_grid%vars(i, j-1, k, E_i_pot) - amr_grid%vars(i, j+1, k, E_i_pot), &
                  amr_grid%vars(i, j, k-1, E_i_pot) - amr_grid%vars(i, j, k+1, E_i_pot)/)
          end do
       end do
    end do

    if (associated(amr_grid%parent)) then
       call set_bc_from_parent(amr_grid, (/E_i_Ex, E_i_Ey, E_i_Ez/))
    else
       call set_bc_from_interior(amr_grid, (/E_i_Ex, E_i_Ey, E_i_Ez/))
    end if

    do nc = 1, amr_grid%n_child
       call compute_field_recursive(amr_grid%children(nc))
    end do
  end subroutine compute_field_recursive

  subroutine set_bc_from_interior(amr_grid, v_ixs)
    type(amr_grid_t), intent(inout) :: amr_grid
    integer, intent(in) :: v_ixs(:)

    integer :: i, j, k, Nx, Ny, Nz

    Nx = amr_grid%Nr(1)
    Ny = amr_grid%Nr(2)
    Nz = amr_grid%Nr(3)

    do j = 1, Ny
       do i = 1, Nx
          amr_grid%vars(i, j, 1, v_ixs)  = 2 * amr_grid%vars(i, j, 2, v_ixs) &
               - amr_grid%vars(i, j, 3, v_ixs)
          amr_grid%vars(i, j, Nz, v_ixs) = 2 * amr_grid%vars(i, j, Nz-1, v_ixs) &
               - amr_grid%vars(i, j, Nz-2, v_ixs)
       end do
    end do

    do k = 1, Nz
       do i = 1, Nx
          amr_grid%vars(i, 1, k, v_ixs)  = 2 * amr_grid%vars(i, 2, k, v_ixs) &
               - amr_grid%vars(i, 3, k, v_ixs)
          amr_grid%vars(i, Ny, k, v_ixs) = 2 * amr_grid%vars(i, Ny-1, k, v_ixs) &
               - amr_grid%vars(i, Ny-2, k, v_ixs)
       end do
    end do

    do k = 1, Nz
       do j = 1, Ny
          amr_grid%vars(1, j, k, v_ixs)  = 2 * amr_grid%vars(2, j, k, v_ixs) &
               - amr_grid%vars(3, j, k, v_ixs)
          amr_grid%vars(Nx, j, k, v_ixs) = 2 * amr_grid%vars(Nx-1, j, k, v_ixs) &
               - amr_grid%vars(Nx-2, j, k, v_ixs)
       end do
    end do
  end subroutine set_bc_from_interior

  subroutine allocate_grid_data(amr_grid)
    type(amr_grid_t), intent(inout) :: amr_grid
    integer                         :: Nr(3)
    Nr = amr_grid%Nr
    allocate(amr_grid%vars(Nr(1), Nr(2), Nr(3), E_n_vars))
    amr_grid%vars = 0.0_dp
    allocate(amr_grid%child_ix(Nr(1), Nr(2), Nr(3)))
    amr_grid%child_ix = 0
  end subroutine allocate_grid_data

  subroutine deallocate_grid_data(amr_grid)
    type(amr_grid_t), intent(inout)  :: amr_grid
    deallocate(amr_grid%vars)
    deallocate(amr_grid%child_ix)
  end subroutine deallocate_grid_data

  recursive subroutine recursive_remove_amr_grid(amr_grid)
    type(amr_grid_t), intent(inout) :: amr_grid
    integer                         :: nc

    do nc = 1, amr_grid%n_child
       call recursive_remove_amr_grid(amr_grid%children(nc))
    end do

    if (associated(amr_grid%children)) then
       deallocate(amr_grid%children)
    end if

    call deallocate_grid_data(amr_grid)
  end subroutine recursive_remove_amr_grid

  subroutine E_set_vars(v_ixs, amounts)
    integer, intent(in)  :: v_ixs(:)
    real(dp), intent(in) :: amounts(:)

    call set_vars_recursive(root_grid, v_ixs, amounts)
  end subroutine E_set_vars

  recursive subroutine copy_var_to_var_recursive(amr_grid, v_i, v_j)
    type(amr_grid_t), intent(inout) :: amr_grid
    integer, intent(in)             :: v_i, v_j
    integer                         :: nc

    amr_grid%vars(:,:,:, v_j) = amr_grid%vars(:,:,:, v_i)
    do nc = 1, amr_grid%n_child
       call copy_var_to_var_recursive(amr_grid%children(nc), v_i, v_j)
    end do
  end subroutine copy_var_to_var_recursive

  recursive subroutine set_vars_recursive(amr_grid, v_ixs, amounts)
    type(amr_grid_t), intent(inout) :: amr_grid
    integer, intent(in)             :: v_ixs(:)
    real(dp), intent(in)            :: amounts(:)
    integer                         :: ix, nc

    do ix = 1, size(v_ixs)
       amr_grid%vars(:,:,:, v_ixs(ix)) = amounts(ix)
    end do

    do nc = 1, amr_grid%n_child
       call set_vars_recursive(amr_grid%children(nc), v_ixs, amounts)
    end do
  end subroutine set_vars_recursive

  recursive subroutine set_source_term_recursive(amr_grid)
    use m_units_constants
    type(amr_grid_t), intent(inout) :: amr_grid
    integer                         :: nc

    amr_grid%vars(:,:,:, E_i_src) = (amr_grid%vars(:,:,:, E_i_elec) &
         + amr_grid%vars(:,:,:, E_i_O2m) + amr_grid%vars(:,:,:, E_i_nion) &
         - amr_grid%vars(:,:,:, E_i_pion)) * abs(UC_elec_q_over_eps0)

    do nc = 1, amr_grid%n_child
       call set_source_term_recursive(amr_grid%children(nc))
    end do
  end subroutine set_source_term_recursive

  subroutine E_add_to_var(v_ix, xyz, amount)
    real(dp), intent(in) :: xyz(3), amount
    integer, intent(in)  :: v_ix
    call add_to_var_recursive(root_grid, v_ix, xyz, amount)
  end subroutine E_add_to_var

  recursive subroutine add_to_var_recursive(amr_grid, v_ix, xyz, amount)
    type(amr_grid_t), intent(inout) :: amr_grid
    real(dp), intent(in)            :: xyz(3), amount
    integer, intent(in)             :: v_ix

    integer                         :: ixs(3), child_ix
    real(dp)                        :: weights(2,2,2)

    call get_ixs_and_weights(amr_grid, xyz, ixs, weights)
    weights = weights * amount * product(amr_grid%inv_dr)

    amr_grid%vars(ixs(1):ixs(1)+1, ixs(2):ixs(2)+1, ixs(3):ixs(3)+1, v_ix) = &
         amr_grid%vars(ixs(1):ixs(1)+1, ixs(2):ixs(2)+1, ixs(3):ixs(3)+1, v_ix) + weights

    child_ix = get_child_ix_at_ixs(amr_grid, ixs)
    if (child_ix > 0) call add_to_var_recursive(amr_grid%children(child_ix), v_ix, xyz, amount)
  end subroutine add_to_var_recursive

  function get_ixs(amr_grid, xyz) result(ixs)
    type(amr_grid_t), intent(in) :: amr_grid
    real(dp), intent(in)         :: xyz(3)
    integer                      :: ixs(3)

    ! Find the indices ixs so that r_rel lies in the cube given by [ixs, ixs+1],
    ! which corresponds to coordinates [(ixs-1)*dr, ixr*dr]
    ixs     = floor((xyz - amr_grid%r_min) * amr_grid%inv_dr) + 1
    where (ixs < 1) ixs = 1
    where (ixs >= amr_grid%Nr) ixs = amr_grid%Nr-1
  end function get_ixs

  subroutine get_ixs_and_weights(amr_grid, xyz, ixs, weights)
    type(amr_grid_t), intent(in) :: amr_grid
    real(dp), intent(in)         :: xyz(3)
    integer, intent(out)         :: ixs(:)
    real(dp), intent(out)        :: weights(:,:,:)
    real(dp)                     :: r_rel(3), l_coeff(3), u_coeff(3), temp(2)

    ! Find the indices ixs so that r_rel lies in the cube given by [ixs, ixs+1],
    ! which corresponds to coordinates [(ixs-1)*dr, ixr*dr]
    ixs     = get_ixs(amr_grid, xyz)
    r_rel   = xyz - amr_grid%r_min
    l_coeff = (ixs * amr_grid%dr - r_rel) * amr_grid%inv_dr
    u_coeff = 1 - l_coeff

    ! Now compute the coefficient of the charge on each of the 8 gridpoints at
    ! the corners of the cube, using linear interpolation (trilinear in this case)
    temp = (/l_coeff(1), u_coeff(1)/) * l_coeff(3)
    weights(1:2,1,1) = l_coeff(2) * temp
    weights(1:2,2,1) = u_coeff(2) * temp
    temp = (/l_coeff(1), u_coeff(1)/) * u_coeff(3)
    weights(1:2,1,2) = l_coeff(2) * temp
    weights(1:2,2,2) = u_coeff(2) * temp
  end subroutine get_ixs_and_weights

  recursive function count_grids(amr_grid) result(n_grids)
    type(amr_grid_t), intent(in) :: amr_grid
    integer                      :: n_grids, nc

    n_grids = 1 ! Self
    do nc = 1, amr_grid%n_child
       n_grids = n_grids + count_grids(amr_grid%children(nc))
    end do
  end function count_grids

  subroutine create_grid_list(amr_grid, grid_list)
    type(amr_grid_t), intent(in), target       :: amr_grid
    type(amr_grid_p), intent(out), allocatable :: grid_list(:)
    integer                                    :: n_grids, counter

    n_grids = count_grids(amr_grid)
    allocate(grid_list(n_grids))
    counter = 0
    call add_grid_to_list(amr_grid, grid_list, counter)
  end subroutine create_grid_list

  recursive subroutine add_grid_to_list(amr_grid, grid_list, counter)
    type(amr_grid_t), intent(in), target :: amr_grid
    type(amr_grid_p), intent(inout)      :: grid_list(:)
    integer, intent(inout)               :: counter
    integer                              :: nc

    counter = counter + 1
    grid_list(counter)%ptr => amr_grid

    do nc = 1, amr_grid%n_child
       call add_grid_to_list(amr_grid%children(nc), grid_list, counter)
    end do
  end subroutine add_grid_to_list

  subroutine E_write_grids(filename, myrank, root, n_cycle, time)
    use m_write_silo
    character(len=*), intent(in)    :: filename
    integer, intent(in)             :: myrank, root, n_cycle
    real(dp), intent(in)            :: time

    integer                         :: ix, v_ix, n_grids
    type(amr_grid_p), allocatable   :: grid_list(:)
    type(amr_grid_t), pointer       :: grid
    character(len=100), allocatable :: mesh_name_list(:), var_name_list(:, :)
    character(len=*), parameter     :: grid_name = "grid_", amr_name = "amr_grid"

    call mpi_collect_recursive(root_grid, &
         (/E_i_pion, E_i_elec, E_i_nion, E_i_O2m/), myrank, root)

    if (myrank == root) then
       call set_bc_from_parent_recursive(root_grid, &
            (/E_i_pion, E_i_elec, E_i_nion, E_i_O2m/))
       call create_grid_list(root_grid, grid_list)
       call SILO_create_file(filename)

       n_grids = size(grid_list)
       allocate(mesh_name_list(n_grids))
       allocate(var_name_list(E_n_vars, n_grids))

       do ix = 1, n_grids
          grid => grid_list(ix)%ptr
          write(mesh_name_list(ix), "(A,I0)") grid_name, ix
          call SILO_add_grid(filename, mesh_name_list(ix), 3, grid%Nr, grid%r_min, grid%dr)

          grid%vars(:,:,:, E_i_tmp) = sqrt(grid%vars(:,:,:, E_i_Ex)**2 + grid%vars(:,:,:, E_i_Ey)**2 &
               + grid%vars(:,:,:, E_i_Ez)**2)

          do v_ix = 1, E_n_vars
             write(var_name_list(v_ix, ix), "(A,I0)") trim(E_var_names(v_ix)) // "_", ix
             call SILO_add_var(filename, var_name_list(v_ix, ix), mesh_name_list(ix), &
                  pack(grid%vars(:,:,:, v_ix), .true.), grid%Nr, E_var_units(v_ix))
          end do
       end do

       call SILO_set_multimesh_grid(filename, amr_name, mesh_name_list, n_cycle, time)
       do v_ix = 1, E_n_vars
          call SILO_set_multimesh_var(filename, trim(E_var_names(v_ix)), amr_name, &
               var_name_list(v_ix, :), n_cycle, time)
       end do
    end if
  end subroutine E_write_grids

end module m_efield_amr
