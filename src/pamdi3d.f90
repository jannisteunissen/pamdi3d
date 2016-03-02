! Copyright 2005-2014, CWI
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

! This is pamdi3d version 0.02
program pamdi3d

  use m_units_constants
  use m_cross_sec
  use m_init_cond
  use m_particle
  use m_particle_core
  use m_particle_par
  use m_efield_amr
  use m_electrode
  use m_config
  use m_misc_process
  use m_photoi
  use m_gas
  use m_random
  use m_phys_domain
  use mpi

  implicit none
  integer, parameter             :: dp = kind(0.0d0)
  character(LEN=100)             :: tmp_name, sim_name, cfg_name
  character(LEN=100)             :: filename, prev_name

  integer                        :: n
  integer                        :: n_samples
  integer                        :: n_gas_comp
  integer                        :: step_cntr, steps_left_fld
  integer                        :: n_output, n_grid_adapt
  integer                        :: WCTime_start, WCTime_current
  integer                        :: n_part_rescale
  integer                        :: n_part_sum, n_part_sum_prev
  integer                        :: ierr, myrank, ntasks, root = 0, ix, n_its
  integer                        :: rng_seed(4)

  logical                        :: finished, flag_output, rescale, adapt_grid
  logical                        :: use_detach, use_photoi, limit_dens

  real(dp)                       :: max_ev

  ! Gas / cross sec parameters
  real(dp)                       :: pressure, temperature
  real(dp), allocatable          :: gas_fracs(:)
  character(len=80), allocatable :: gas_names(:), gas_files(:)
  type(CS_t), allocatable        :: cross_secs(:)

  ! Simulation variables
  real(dp)                       :: sim_time, end_time
  real(dp)                       :: time_per_output, time_per_grid_adapt
  real(dp)                       :: dt_fld, min_incr_rescale, prev_fld_time
  real(dp)                       :: fld_err
  real(dp)                       :: max_fld_err
  real(dp)                       :: n_elec_sum
  real(dp)                       :: dt_min, dt_max, dt, dt_next
  real(dp)                       :: mean_weight
  real(dp)                       :: cfl_num
  real(dp)                       :: max_density

  type(CFG_t)                    :: cfg
  type(PC_t)                     :: pc
  type(RNG_t)                    :: rng

  sim_time = 0.0D0

  call MPI_init(ierr)
  call MPI_comm_rank(MPI_COMM_WORLD, myrank, ierr)
  call MPI_comm_size(MPI_COMM_WORLD, ntasks, ierr)
  call create_config(cfg)

  sim_name = "sim"
  prev_name = ""
  do ix = 1, command_argument_count()
     call get_command_argument(ix, cfg_name)
     call CFG_read_file(cfg, trim(cfg_name))

     call CFG_get(cfg, "sim_name", tmp_name)
     if (tmp_name /= "" .and. tmp_name /= prev_name) &
          sim_name = trim(sim_name) // "_" // trim(tmp_name)
     prev_name = tmp_name
  end do

  if (myrank == root) then
     print *, "Starting simulation ** " // trim(sim_name) // " **"
     call CFG_write(cfg, "output/" // trim(sim_name) // "_config.txt")
  end if

  ! Set rng seed (note: particle model has its own rng)
  call CFG_get(cfg, "rng_seed", rng_seed)
  call rng%set_seed(rng_seed)

  ! Read in the crossection data from files
  call CFG_get_size(cfg, "gas_names", n_gas_comp)
  call CFG_get_size(cfg, "gas_fractions", n)
  if (n_gas_comp /= n) then
     print *, "gas_names and gas_fracs have unequal size"
     stop
  end if

  allocate(gas_names(n_gas_comp))
  allocate(gas_fracs(n_gas_comp))
  allocate(gas_files(n_gas_comp))
  call CFG_get(cfg, "gas_names", gas_names)
  call CFG_get(cfg, "gas_fractions", gas_fracs)
  call CFG_get(cfg, "gas_files", gas_files)
  call CFG_get(cfg, "gas_temperature", temperature)
  call CFG_get(cfg, "gas_pressure", pressure)
  call CFG_get(cfg, "part_max_ev", max_ev)
  call GAS_initialize(gas_names, gas_fracs, pressure, temperature)

  if (myrank == root) print *, "Reading crossection data"
  do n = 1, n_gas_comp
     call CS_add_from_file("input/" // gas_files(n), gas_names(n), &
          gas_fracs(n) * GAS_number_dens, max_ev, cross_secs)
  end do

  ! Initialize all modules that need initialization
  call PD_set(cfg)
  if (myrank == root) print *, " ~~~ initializing electric fld module"
  call E_initialize(cfg)
  if (myrank == root) print *, " ~~~ initializing electrode module"
  if (PD_use_elec) call EL_initialize(cfg, rng, PD_r_max, myrank, root)

  if (myrank == root) print *, " ~~~ initializing particle module"
  call PM_initialize(pc, cross_secs, cfg, myrank, ntasks)

  call CFG_get(cfg, "photoi_enabled", use_photoi)
  call CFG_get(cfg, "sim_use_o2m_detach", use_detach)

  if (myrank == root) print *, " ~~~ initializing misc module"
  call MISC_initialize(cfg)

  if (use_photoi) then
     if (myrank == root) print *, " ~~~ initializing photoi module"
     call pi_initialize(cfg)
  end if

  if (myrank == root) then
     call system_clock(COUNT=WCTime_start)
     print *, 'Setting up initial electron/ion positions and fld'
     call IC_set_init_cond(pc, cfg, rng)
  end if

  if (myrank == root) print *, "Sharing initial particles"
  call PP_share_mpi(pc, myrank, ntasks)

  if (myrank == root) print *, "Computing initial fld"
  call PM_particles_to_density(pc)
  call E_compute_field(myrank, root, sim_time)
  call E_update_grids(myrank, root, sim_time)
  call PM_particles_to_density(pc)
  if (PD_use_elec) call E_readjust_elec_new_grid(myrank, root, &
       sim_time, n_its)
  call E_compute_field(myrank, root, sim_time)

  if (PD_use_elec) then
     do ix = 1, 10
        call E_compute_field(myrank, root, sim_time)
     end do

     call E_update_grids(myrank, root, sim_time)
     call PM_particles_to_density(pc)
     call E_readjust_elec_new_grid(myrank, root, sim_time, n_its)

     call E_compute_field(myrank, root, sim_time)
     do ix = 1, 10
        call E_compute_field(myrank, root, sim_time)
     end do
  end if

  call E_update_grids(myrank, root, sim_time)
  call PM_particles_to_density(pc)
  if (PD_use_elec) call E_readjust_elec_new_grid(myrank, root, &
       sim_time, n_its)
  call E_compute_field(myrank, root, sim_time)
  call pc%set_accel(E_get_accel_part)

  ! Initial values
  finished        = .false.
  flag_output     = .false.
  rescale         = .false.
  adapt_grid      = .false.

  n_output        = 0
  n_grid_adapt    = 0
  step_cntr       = 0
  steps_left_fld  = 0
  n_part_sum_prev = PP_get_num_sim_part(pc)

  prev_fld_time   = 0.0D0
  max_fld_err     = 0.0D0
  fld_err         = 0.0D0

  call CFG_get(cfg, "part_n_rescale", n_part_rescale)
  call CFG_get(cfg, "part_n_samples", n_samples)
  call CFG_get(cfg, "sim_cfl_num", cfl_num)
  call CFG_get(cfg, "t_end", end_time)
  call CFG_get(cfg, "dt_output", time_per_output)
  call CFG_get(cfg, "dt_grid_adapt", time_per_grid_adapt)
  call CFG_get(cfg, "dt_min", dt_min)
  call CFG_get(cfg, "dt_max", dt_max)
  call CFG_get(cfg, "sim_max_fld_err", max_fld_err)
  call CFG_get(cfg, "sim_min_incr_rescale", min_incr_rescale)
  call CFG_get(cfg, "sim_max_density", max_density)

  dt_next    = dt_min
  dt_fld     = dt_min
  limit_dens = .true.

  if (myrank == root) then
     print *, "Starting simulation with name ", trim(sim_name), ","
     print *, "will simulate up to ", 1.0e9_dp * end_time, "ns"
  end if

  ! Write the output
  n_output = 1
  write(filename, "(A,I6.6,A)") "output/" // trim(sim_name) // "_", &
       n_output, ".silo"
  call E_write_grids(filename, myrank, root, n_output, sim_time)

  do while (.not. finished)
     step_cntr      = step_cntr + 1
     dt             = dt_next
     steps_left_fld = steps_left_fld - 1
     sim_time       = sim_time + dt

     call pc%advance(dt)
     call pc%correct_new_accel(dt, E_get_accel_part)

     if (steps_left_fld <= 0) then
        n_part_sum = PP_get_num_sim_part(pc)
        n_elec_sum = PP_get_num_real_part(pc)

        call MISC_gamma(sim_time - prev_fld_time)

        if (use_detach) call MISC_detachment(pc, rng, &
             sim_time - prev_fld_time, myrank, root)

        call PM_particles_to_density(pc)
        call PM_fld_error(pc, rng, n_samples, fld_err, store_samples=.true.)
        call E_compute_field(myrank, root, sim_time)
        call PM_fld_error(pc, rng, n_samples, fld_err, store_samples=.false.)
        if (myrank == root) print *, "Field error", fld_err
        call pc%set_accel(E_get_accel_part)

        ! Redistribute the particles
        call PP_share_mpi(pc, myrank, ntasks)

        dt_next = PM_get_max_dt(pc, rng, n_samples, cfl_num)
        dt_next = min(dt_max, dt_next)

        if (myrank == root) then
           dt_fld         = get_new_dt(dt_fld, fld_err, max_fld_err)
           steps_left_fld = ceiling(dt_fld / dt_next)
           dt_next        = dt_fld / steps_left_fld

           ! Centralize these flags because for heterogeneous systems
           finished    = sim_time > end_time
           flag_output = n_output <= sim_time / time_per_output .or. finished
           rescale     = n_part_sum > n_part_rescale .and. &
                n_part_sum > int(n_part_sum_prev * min_incr_rescale)
           rescale     = rescale
           adapt_grid  = n_grid_adapt <= sim_time / time_per_grid_adapt
        end if

        ! Communicate the new step sizes and flags
        call MPI_BCAST(dt_next, 1, MPI_DOUBLE_PRECISION, &
             root, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(steps_left_fld, 1, MPI_INTEGER, &
             root, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(flag_output, 1, MPI_LOGICAL, root, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(adapt_grid,  1, MPI_LOGICAL, root, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(finished,   1, MPI_LOGICAL, root, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(rescale,    1, MPI_LOGICAL, root, MPI_COMM_WORLD, ierr)

        ! Write the output
        if (flag_output) then
           n_output = n_output + 1
           write(filename, "(A,I6.6,A)") "output/" // trim(sim_name) // "_", &
                n_output, ".silo"
           call E_write_grids(filename, myrank, root, n_output, sim_time)
        end if

        if (adapt_grid) then
           n_grid_adapt = n_grid_adapt + 1
           call PM_particles_to_density(pc)
           call E_update_grids(myrank, root, sim_time)
           call PM_particles_to_density(pc)
           if (PD_use_elec) call E_readjust_elec_new_grid(myrank, root, &
                sim_time, n_its)
           call E_compute_field(myrank, root, sim_time)
        end if

        if (adapt_grid .or. rescale) then
           call E_share_vars((/E_i_elec/), root)

           if (limit_dens) then
              if (myrank == root) print *, "Limiting density"
              call PM_limit_dens(pc, rng, max_density)
           end if

           ! TODO: The above call modifies the electron density slightly,
           ! probably not very important but check it later.
           call PP_reorder_by_bins_mpi(pc, PM_binner, myrank, ntasks)

           if (myrank == root) then
              mean_weight = n_elec_sum / (n_part_sum + epsilon(1.0_dp))
              print *, step_cntr, "Rescaling particles, mean_weight:", &
                   mean_weight, "nParticles", n_part_sum
           end if

           call PM_fld_error(pc, rng, n_samples, fld_err, store_samples=.true.)
           call PM_adjust_weights(pc)

           n_part_sum = PP_get_num_sim_part(pc)
           n_part_sum_prev = n_part_sum

           call PM_particles_to_density(pc)
           call E_compute_field(myrank, root, sim_time)
           call PM_fld_error(pc, rng, n_samples, fld_err, store_samples=.false.)

           if (myrank == root) print *, "Fld error due to rescaling:", fld_err

           call pc%set_accel(E_get_accel_part)

           ! Redistribute the particles
           call PP_share_mpi(pc, myrank, ntasks)
        end if

        if (myrank == root) then
           call system_clock(COUNT=WCTime_current)
           write(*, FMT = "((I7),(A),(es10.4e2),(A),(es10.2e2),(A),(es10.2e2))") &
                n_output, "   t(s): ", sim_time, "   dt: ", dt_next, &
                "   dtPS", dt_fld
           write(*, FMT = "((A)(es10.2e2)(A),(I12),(A),(es10.2e2),(A),(es10.4e2))") &
                "       WC t(s): ", &
                & (WCTime_current - WCTime_start)/1.0D3, "   nPart:", n_part_sum, &
                & "  nElec: ", n_elec_sum
        end if
        prev_fld_time = sim_time
     end if

     ! The end of the simulation main loop
  end do

  if (myrank == root) print *, "The simulation has finished succesfully"
  call MPI_FINALIZE(ierr)

contains

  ! Given the old stepsize 'oldDt', the error 'err', the maximum allowed error
  ! 'maxErr', and the relative change in number of particles, return a new
  ! stepsize.
  real(dp) function get_new_dt(oldDt, err, maxErr)
    real(dp), intent(IN) :: oldDt
    real(dp), intent(IN) :: err, maxErr

    ! Adjust the timesteps to obtain an error close to the desired one. The
    ! interpolation error should be ~ interp_errFrac * max_fld_err
    if (err < 0.5D0 * maxErr) then
       get_new_dt = min(oldDt * (maxErr/(err+epsilon(err)))**0.1D0, 2.0D0 * oldDt)
    else if (err > maxErr) then
       get_new_dt = oldDt * (maxErr/err)
    else
       get_new_dt = oldDt
    end if

    get_new_dt = min(max(get_new_dt, dt_min), dt_max)

  end function get_new_dt

  subroutine create_config(cfg)
    use m_config
    use m_units_constants
    type(CFG_t), intent(inout) :: cfg

    ! General simulation parameters
    call CFG_add(cfg, "sim_name", "mcstr", &
         "The name of the simulation")
    call CFG_add(cfg, "sim_max_fld_err", 5.0D-2, &
         "The maximum allowed error in the electric field before recomputing it")
    call CFG_add(cfg, "sim_max_electrons", 1.0d99, &
         "Stop the simulation when this many electrons are in the simulation")
    CALL CFG_add(cfg, "sim_min_incr_rescale", 1.2_dp, &
         & "The minimum increase in number of particles before their weights are rescaled (superparticle creation)")
    CALL CFG_add(cfg, "sim_n_runs_max", 1, &
         & "The number of runs to perform for simulation types that require this")
    call CFG_add(cfg, "sim_use_o2m_detach", .false., &
         "Whether we use background O2- in the simulation")
    call CFG_add(cfg, "sim_cfl_num", 0.5_dp, &
         "CFL number used for the particles")
    call CFG_add(cfg, "sim_max_density", 2.5e21_dp, &
         "Maximum electron density; above this they are converted to ions")

    call CFG_add(cfg, "rng_seed", (/521288629, 362436069, 16163801, &
         1131199299/), "Seed for random number generator")

    call CFG_add(cfg, "ref_max_levels", 7, &
         "The maximum number of refinement levels in the electric field")
    call CFG_add(cfg, "ref_min_lvl_electrode", 3, &
         "The minimum refinement lvl around the electrode tip")
    call CFG_add(cfg, "ref_min_grid_size", 8, &
         "The minimum number of grid points in any direction of a refined grid")
    call CFG_add(cfg, "ref_buffer_width", 2, &
         "Increase refinement area by this size")
    call CFG_add(cfg, "ref_min_grid_separation", 2, &
         "The minimum distance between grids at the same level")
    CALL CFG_add(cfg, "ref_delta_values",  (/2.0d-6,   3.0d-6,  5.0d-6,  1.0d-5,  1.0d-4,  1.0d-3/), &
         & "List of spatial stepsizes for the electric fields in ref_EfieldValues", dyn_size = .true.)
    CALL CFG_add(cfg, "ref_max_efield_at_delta", (/3.0d7,    2.0d7,   1.5d7,   1.0d7,   3.0d6,   0.0d0/), &
         & "List of electric field values at which we specify the max. spatial stepsize", dyn_size = .true.)
    CALL CFG_add(cfg, "ref_min_elec_dens",  1.0d15, &
         & "Minimum electron density for considering refining in a region")

    !! Gas parameters
    call CFG_add(cfg, "gas_pressure", 1.0D0, &
         "The gas pressure (bar)")
    call CFG_add(cfg, "gas_temperature", 293.0D0, &
         "The gas temperature (Kelvin)")
    CALL CFG_add(cfg, "gas_names", (/"N2"/), &
         & "The names of the gases used in the simulation", dyn_size = .true.)
    CALL CFG_add(cfg, "gas_files", (/"cs_example.txt"/), &
         & "The files in which to find cross section data for each gas", dyn_size = .true.)
    CALL CFG_add(cfg, "gas_fractions", (/1.0_dp /), &
         & "The partial pressure of the gases (as if they were ideal gases)", dyn_size = .true.)

    ! Electric field parameters
    CALL CFG_add(cfg, "sim_efield_values", (/-7.0D6/), &
         "A list of values of the electric field applied in the z-direction", dyn_size = .true.)
    CALL CFG_add(cfg, "sim_efield_times", (/0.0D0/), &
         "The times (in increasing order) at which the values of the applied electric field are used", dyn_size = .true.)
    CALL CFG_add(cfg, "sim_constant_efield", .false., &
         & "Whether the electric field is constant in space and time (no space charge effect)")

    ! Initial conditions
    call CFG_add(cfg, "init_cond_type", "Gaussian", &
         "The type of initial condition")
    call CFG_add(cfg, "init_seed_pos_radius", 1.0D-4, &
         "The radius of the distribution of initial ion seeds")
    call CFG_add(cfg, "init_weight", 1.0_dp, &
         "The initial weight of electrons")
    call CFG_add(cfg, "init_n_ion_pairs", 1.0d2, &
         "The number of initial ion pairs")
    call CFG_add(cfg, "init_line_dens", 1.d10 , &
         "The density for some types of initial conditions")
    call CFG_add(cfg, "init_laser_direction", (/1.0d0, 0.0d0, 1.0d0/) , &
         "Direction of the initial laserbeam")
    call CFG_add(cfg, "init_laser_line_offset", (/0.5d0, 0.5d0, 0.6d0/) , &
         "Laser goes through this point")
    call CFG_add(cfg, "init_rel_pos", (/0.5D0, 0.5D0, 0.5D0/), &
         "The relative position of the ion pairs in the domain")
    call CFG_add(cfg, "init_bg_dens", 0.0D0, &
         "The background ion/electron density in the domain")
    call CFG_add(cfg, "init_o2m_bg_dens", 0.0D0, &
         "The background O2- density in the domain")

    ! Grid related parameters
    call CFG_add(cfg, "grid_size", (/33, 33, 33/), &
         "The number of grid cells in each direction")
    call CFG_add(cfg, "grid_delta", (/5.0e-4_dp, 5.0e-4_dp, 5.0e-4_dp/), &
         "The length of a grid cell in each direction")
    call CFG_add(cfg, "grid_plasma_min_rel_pos", (/0.1d0, 0.1d0, 0.1d0/), &
         "Plasma can only be above these relative coords")
    call CFG_add(cfg, "grid_plasma_max_rel_pos", (/0.9d0, 0.9d0, 0.9d0/), &
         "Plasma can only be below these relative coords")


    ! Particle model related parameters
    call CFG_add(cfg, "part_max_num", 5*1000*1000, &
         "The maximum total number of particles")
    call CFG_add(cfg, "part_n_samples", 1000, &
         "The number of samples to use (for estimates)")
    call CFG_add(cfg, "part_n_rescale", 1000*1000, &
         "The number of particles at which to start rescaling")
    call CFG_add(cfg, "part_max_ev", 1000.0D0, &
         "The maximum energy in ev for particles in the simulation")
    call CFG_add(cfg, "part_lkp_tbl_size", 10000, &
         "The table size to use in the particle model for lookup tables")
    CALL CFG_add(cfg, "part_per_cell", 64.0_dp, &
         & "The minimum number of particles in a cell to consider rescaling them")
    CALL CFG_add(cfg, "part_max_weight", 1.0e3_dp, &
         & "Only particles with a weight < maxWeight will be rescaled")
    call CFG_add(cfg, "part_v_rel_weight", 1.0e-12_dp, &
         & "Scale factor for v coords compared to x coords when merging")

    ! Parameters for the output of the time-integrated ionization
    call CFG_add(cfg, "gamma_decay_time", 1.0e-9_dp, &
         "Decay time for time-integrated ionization density")

    ! Photoionization parameters
    call CFG_add(cfg, "photoi_enabled", .false., &
         "Whether photoionization is used")
    CALL CFG_add(cfg, "photoi_efield_table", (/ 0.0D0, 0.25D7, 0.4D7, 0.75D7, 1.5D7, 3.75D7/), &
         & "The tabulated values of the electric field (for the photo-efficiency)")
    CALL CFG_add(cfg, "photoi_efficiency_table", (/0.0D0, 0.05D0, 0.12D0, 0.08D0, 0.06D0, 0.04D0/), &
         & "The tabulated values of the the photo-efficiency")
    CALL CFG_add(cfg, "photoi_frequencies", (/2.925D12, 3.059D12/), &
         & "The lower/upper bound for the frequency of the photons, not currently used")
    CALL CFG_add(cfg, "photoi_absorp_inv_lengths", (/3.5D0 / UC_torr_to_bar, 200D0 / UC_torr_to_bar/), &
         & "The inverse min/max absorption length, have to be scaled by oxygen partial pressure")
    call CFG_add(cfg, "photoi_mean_life_time_excited", 0.1d-10, &
         "The mean life time of excited species")

    ! Timesteps / end time
    call CFG_add(cfg, "t_end", 5.0D-9, &
         "The desired endtime in seconds of the simulation")
    call CFG_add(cfg, "dt_min", 1.0D-13, &
         "The minimal timestep in seconds")
    call CFG_add(cfg, "dt_max", 1.0D-11, &
         "The maximal timestep in seconds")
    call CFG_add(cfg, "dt_output", 1.0D-11, &
         "The timestep for writing output")
    call CFG_add(cfg, "dt_grid_adapt",  2.0d-11, &
         & "Time between consecutive adaptations of the grid structure")


    ! Electrode parameters
    call CFG_add(cfg, "elec_enabled", .false., &
         "Whether to use an electrode in the simulation")
    call CFG_add(cfg, "elec_bc_type", 1, &
         "Boundary cond. 1: linear increase along sides. 2: grounded sides.")
    call CFG_add(cfg, "elec_times", (/0.0d-9, 1.0d-9, 2.0d-9/), &
         "The times at which the electrode voltage is specified.", .true.)
    call CFG_add(cfg, "elec_voltages", (/0.0d3, 1.0d3, 5.0d3/), &
         "The voltages at the specified times.", .true.)
    call CFG_add(cfg, "elec_h_cone", 0.25D-3, &
         "Height of conical part of electrode")
    call CFG_add(cfg, "elec_r_cyl", 2.5D-4, &
         "Radius of cylindrical part of electrode")
    call CFG_add(cfg, "elec_rc_tip", 1.0D-4, &
         "Radius of curvature of the tip of the electrode")
    call CFG_add(cfg, "elec_voltage", 1.0d3, &
         "The voltage on the electrode")
    call CFG_add(cfg, "elec_rel_pos", (/0.5d0, 0.5d0, 0.25d0/), &
         "The relative position of the electrode in the domain")
    CALL CFG_add(cfg, "elec_rc_trans", 0.0D-5, &
         & "Radius of curvature of the transition from the cylindrical to conical region of the electrode")
    call CFG_add(cfg, "elec_spacing", 2.0d-6, &
         "The distance between points on the electrode surface")
    CALL CFG_add(cfg, "elec_n_surface_points", 100, &
         & "The number of points on the electrode surface")
    CALL CFG_add(cfg, "elec_surface_points_scale", 1.0D0, &
         & "The scale factor for the distribution of surface points on the electrode.")
    call CFG_sort(cfg)
  end subroutine create_config

end program pamdi3d
