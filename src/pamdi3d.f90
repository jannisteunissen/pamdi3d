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
  use m_gas
  use m_random
  use m_phys_domain
  use mpi

  implicit none
  integer, parameter :: dp = kind(0.0d0)
  character(LEN=80) :: sim_name, cfg_name, filename

  integer :: n
  integer :: n_fld_samples
  integer :: n_gas_comp
  integer  :: step_cntr, steps_left_fld
  integer  :: n_output, n_grid_adapt
  integer  :: WCTime_start, WCTime_current
  integer  :: n_part_rescale
  integer  :: n_part_sum, n_part_sum_prev
  integer  :: ierr, myrank, ntasks, root = 0, ix, n_its
  
  logical  :: finished, flag_output, rescale, adapt_grid
  logical :: use_detach, use_photoion, do_fake_attach, use_electrode
  
  ! Gas / cross sec parameters
  real(dp)                       :: pressure, temperature
  real(dp), allocatable          :: gas_fracs(:)
  character(len=80), allocatable :: gas_names(:), gas_files(:)
  type(CS_t), allocatable        :: cross_secs(:)

  ! Simulation variables
  real(dp) :: sim_time, end_time
  real(dp) :: time_per_output, time_per_grid_adapt
  real(dp) :: dt_fld, min_incr_rescale, prev_fld_time
  real(dp) :: interp_err, fld_err
  real(dp) :: max_fld_err, max_interp_err
  real(dp) :: n_elec_sum
  real(dp) :: dt_min, dt_max, dt, dt_next
  real(dp) :: meanWeight

  ! Particle binning
  real(dp) :: bin_args(2)
  integer :: n_bins

  type(PC_t) :: pc
  type(RNG_t) :: rng

  call MPI_init(ierr)

  call create_config()

  do ix = 1, command_argument_count()
     call get_command_argument(ix, cfg_name)
     call CFG_read_file(trim(cfg_name))
  end do

  call CFG_get("sim_name", sim_name)
  if (myrank == root) then
     call CFG_write("output/" // trim(sim_name) // "_config.txt")
  end if

  ! Read in the crossection data from files
  n_gas_comp = CFG_get_size("gas_component_names")
  if (n_gas_comp /= CFG_get_size("gas_component_fractions")) then
     print *, "gas_component_names and gas_component_fracs have unequal size"
     stop
  end if

  allocate(gas_names(n_gas_comp))
  allocate(gas_fracs(n_gas_comp))
  allocate(gas_files(n_gas_comp))
  call CFG_get("gas_names", gas_names)
  call CFG_get("gas_fractions", gas_fracs)
  call CFG_get("gas_files", gas_files)

  call GAS_initialize(gas_names, gas_fracs, pressure, temperature)
  print *, "Reading crossection data"
  do n = 1, n_gas_comp
     call CS_add_from_file(gas_files(n), gas_names(n), 1.0_dp, &
          gas_fracs(n) * GAS_get_number_dens(), &
          CFG_get_real("part_max_energy_eV"), cross_secs)
  end do
  call CFG_write( "output/" // trim(sim_name) // "_config.txt" )

  use_electrode = CFG_get_logic("sim_use_electrode")

  ! Initialize all modules that need initialization
  print *, " ~~~ initializing electric fld module"
  call E_initialize()
  print *, " ~~~ initializing electrode module"
  call EL_initialize(rng, myrank, root)
  
  print *, " ~~~ initializing particle module"
  call pc%initialize(UC_elec_mass, cross_secs, CFG_get_int("part_lkp_tbl_size"), &
       CFG_get_real("part_max_eV"), CFG_get_int("part_max_num"))

  use_photoion = CFG_get_logic("sim_use_photoionization")
  use_detach = CFG_get_logic("sim_use_detachment")
  if (use_photoion .or. use_detach) then
     print *, " ~~~ initializing photoionization module"
     call MISC_initialize()
  end if

  if (myrank == root) then
     call system_clock(COUNT=WCTime_start)
     print *, 'Setting up initial electron/ion positions and Efld'
     call IC_set_init_cond(rng)
  end if

  if (myrank == root) print *, "Sharing initial particles"
  call PP_share_mpi(pc, myrank, ntasks)

  if (myrank == root) print *, "Computing initial fld"
  call PM_particles_to_density(pc)
  call E_compute_field(myrank, root, sim_time)
  call E_update_grids(myrank, root, sim_time)
  call PM_particles_to_density(pc)
  if (use_electrode) call E_readjust_electrode_new_grid(myrank, root, &
       sim_time, n_its)
  call E_compute_field(myrank, root, sim_time)

  if (use_electrode) then
     do ix = 1, 10
        call E_compute_field(myrank, root, sim_time)
     end do

     call E_update_grids(myrank, root, sim_time)
     call PM_particles_to_density(pc)
     call E_readjust_electrode_new_grid(myrank, root, sim_time, n_its)

     call E_compute_field(myrank, root, sim_time)
     do ix = 1, 10
        call E_compute_field(myrank, root, sim_time)
     end do
  end if

  if (myrank == root) call E_benchmark()

  call E_update_grids(myrank, root, sim_time)
  call PM_particles_to_density(pc)
  if (use_electrode) call E_readjust_electrode_new_grid(myrank, root, &
       sim_time, n_its)
  call E_compute_field(myrank, root, sim_time)
  call pc%set_accel(E_get_accel_part)

  ! Initial values
  finished            = .false.
  flag_output          = .false.
  rescale             = .false.
  adapt_grid           = .false.

  n_output            = 0
  n_grid_adapt        = 0
  step_cntr           = 0
  steps_left_fld      = 0
  n_part_sum_prev     = PM_get_num_part_mpi(pc)

  sim_time            = 0.0D0
  prev_fld_time       = 0.0D0
  max_fld_err         = 0.0D0
  interp_err          = 0.0D0
  fld_err             = 0.0D0

  n_part_rescale      = CFG_get_int("part_n_rescale")
  end_time            = CFG_get_real("sim_end_time")
  time_per_output     = CFG_get_real("output_timestep")
  time_per_grid_adapt = CFG_get_real("ref_timePerAdaptation")
  dt_min              = CFG_get_real("sim_dt_min")
  dt_max              = CFG_get_real("sim_dt_max")
  max_fld_err         = CFG_get_real("sim_max_fld_err")
  max_interp_err      = CFG_get_real("sim_max_interp_err")
  min_incr_rescale    = CFG_get_real("sim_min_incr_rescale")
  use_detach          = CFG_get_logic("sim_use_O2Min_bg")

  bin_args(1) = 0
  bin_args(2) = 1/PD_dr(1)
  n_bins = PD_r_max(1) / PD_dr(1)

  dt_next             = dt_min
  dt_fld              = dt_min
  do_fake_attach      = .true.

  call E_benchmark()

  if (myrank == root) then
     print *, "Starting simulation with name ", trim(sim_name), ","
     print *, "will simulate up to ", 1.0e9_dp * end_time, "ns"
  end if

  ! Write the output
  n_output = 1
  write(filename, "(A,I6,A)") "output/" // trim(sim_name) // "_", &
       n_output, ".silo"
  print *, trim(filename)
  call E_write_grids(filename, myrank, root, n_output, sim_time)

  do while (.not. finished)
     step_cntr      = step_cntr + 1
     dt             = dt_next
     steps_left_fld = steps_left_fld - 1
     sim_time       = sim_time + dt

     call pc%advance(dt)

     if (steps_left_fld <= 0) then
        n_part_sum = PM_get_num_part_mpi(pc)

        if (use_photoion) call MISC_photoionization(pc, rng, sim_time - prev_fld_time, myrank, root)
        if (use_detach) call MISC_detachment(pc, rng, sim_time - prev_fld_time, myrank, root)

        call PM_particles_to_density(pc)
        call PM_fld_error(pc, rng, n_fld_samples, fld_err, only_store=.true.)
        call E_compute_field(myrank, root, sim_time)
        call PM_fld_error(pc, rng, n_fld_samples, fld_err, only_store=.false.)
        call pc%set_accel(E_get_accel_part)

        ! Redistribute the particles
        call PP_share_mpi(pc, myrank, ntasks)

        dt_next = PM_get_max_dt(pc, myrank, root)

        if (myrank == root) then
           dt_fld         = get_new_dt(dt_fld, fld_err, max_fld_err)
           dt_next        = min(dt_next, dt_fld)
           steps_left_fld = floor(dt_fld / dt_next)

           ! Centralize these flags because they might differ on heterogeneous systems
           finished          = sim_time > end_time
           flag_output        = n_output <= sim_time / time_per_output .or. finished
           rescale           = n_part_sum > n_part_rescale .and. &
                n_part_sum > int(n_part_sum_prev * min_incr_rescale)
           rescale           = rescale .or. step_cntr == 1
           adapt_grid        = n_grid_adapt <= sim_time / time_per_grid_adapt
        end if

        ! Communicate the new step sizes and flags
        call MPI_BCAST(dt_next, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(steps_left_fld, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(flag_output, 1, MPI_LOGICAL, root, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(adapt_grid,  1, MPI_LOGICAL, root, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(finished,   1, MPI_LOGICAL, root, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(rescale,    1, MPI_LOGICAL, root, MPI_COMM_WORLD, ierr)

        ! Write the output
        if (flag_output) then
           n_output = n_output + 1
           write(filename, "(A,I6,A)") "output/" // trim(sim_name) // "_", &
                n_output, ".silo"
           print *, trim(filename)
           call E_write_grids(filename, myrank, root, n_output, sim_time)
        end if

        if (adapt_grid) then
           n_grid_adapt = n_grid_adapt + 1
           call PM_particles_to_density(pc)
           call E_update_grids(myrank, root, sim_time)
           call PM_particles_to_density(pc)
           if (use_electrode) call E_readjust_electrode_new_grid(myrank, root, &
                sim_time, n_its)
           call E_compute_field(myrank, root, sim_time)
        end if

        if (adapt_grid .or. rescale) then
           call E_share_vars((/E_i_elec/), root)

           print *, "TODO"
           ! if (do_fake_attach) then
           !    if (myrank == root) print *, "Performing fake attachments!"
           !    call PM_convertElecForStability(dt_fld, 0.2d0) ! Second arg is max mobility estimate
           ! end if

           ! TODO: The above call modifies the electron density slightly,
           ! probably not very important but check it later.
           call PP_reorder_by_bins_mpi(pc, PM_bin_func, n_bins, bin_args, myrank, ntasks)

           if (myrank == root) then
              print *, step_cntr, "Rescaling particles, meanWeight:", meanWeight, "nParticles", n_part_sum
           end if

           call PM_fld_error(pc, rng, n_fld_samples, fld_err, only_store=.true.)
           call PM_mergeAndSplit(myrank)

           n_part_sum_prev  = PM_get_num_part_mpi(pc)
           n_part_sum             = n_part_sum_prev

           call PM_particles_to_density(pc)
           call E_compute_field(myrank, root, sim_time)
           call PM_fld_error(pc, rng, n_fld_samples, fld_err, only_store=.false.)

           if (myrank == root) then
              print *, "Fld error due to rescaling:", fld_err
           end if

           call pc%set_accel(E_get_accel_part)

           ! Redistribute the particles
           call PP_share_mpi(pc, myrank, ntasks)
        end if

        if (myrank == root) then
           call system_clock(COUNT=WCTime_current)
           write(*, FMT = "((I7),(A),(es10.4e2),(A),(es10.2e2),(A),(es10.2e2))") n_output, &
                & "   t(s): ", sim_time, "   dt: ", dt_next, "   dtPS", dt_fld
           write(*, FMT = "((A)(es10.2e2)(A),(I12),(A),(es10.2e2),(A),(es10.4e2))") "       WC t(s): ", &
                & (WCTime_current - WCTime_start)/1.0D3, "   nPart:", n_part_sum, &
                & "  nElec: ", n_elec_sum
           ! print *, "InterpErr:", interp_err, "fld_err:", fld_err
           ! print *, "maxPSerr:", maxPSerr
           ! write(*,*) ""
        end if
        prev_fld_time = sim_time
     end if

     ! The end of the simulation main loop
  end do

  call MPI_FINALIZE(ierr)
  print *, "The simulation has finished succesfully"

contains

  !> Given the old stepsize 'oldDt', the error 'err', the maximum allowed error 'maxErr',
  !! and the relative change in number of particles, return a new stepsize.
  real(dp) function get_new_dt(oldDt, err, maxErr)
    real(dp), intent(IN) :: oldDt       ! The old timestep
    real(dp), intent(IN) :: err, maxErr  ! The actual error and the maximum / prev. error

    print *, "err", err
    ! Adjust the timesteps to obtain an error close to the desired one.
    ! The interpolation error should be ~ interp_errFrac * max_fld_err
    if (err < 0.5D0 * maxErr) then
       get_new_dt = min(oldDt * (maxErr/(err+epsilon(err)))**0.1D0, 2.0D0 * oldDt)
    else if (err > maxErr) then
       get_new_dt = oldDt * (maxErr/err)
    else
       get_new_dt = oldDt
    end if

    get_new_dt = min(max(get_new_dt, dt_min), dt_max)

  end function get_new_dt

  subroutine create_config()
    use m_units_constants

    ! General simulation parameters
    CALL CFG_add("sim_end_time", 5.0D-9, "The desired endtime in seconds of the simulation")
    CALL CFG_add("sim_dt_min", 5.0D-15, "The minimal timestep in seconds")
    CALL CFG_add("sim_dt_max", 1.0D-11, "The maximal timestep in seconds")
    CALL CFG_add("sim_use_electrode", .false., "Whether to use an electrode in the simulation")
    CALL CFG_add("sim_name", "mcstr", "The name of the simulation")
    CALL CFG_add("sim_nSamplesEfield", 1000, "The number of samples to take to estimate the electric field error")
    CALL CFG_add("sim_maxEfieldErr", 1.0D-2, "The maximum allowed error in the electric field before recomputing it")
    CALL CFG_add("sim_maxElectrons", 1.0d99, "Stop the simulation when this many electrons are in the simulation")
    CALL CFG_add("sim_max_interp_err", 1.0D-2, &
         & "The maximum allowed error in the acceleration of particles, before interpolating a new value from the Efield")
    CALL CFG_add("sim_min_incr_rescale", 1.5_dp, &
         & "The minimum increase in number of particles before their weights are rescaled (superparticle creation)")
    CALL CFG_add("sim_nRunsMax", 1, &
         & "The number of runs to perform for simulation types that require this")
    CALL CFG_add("sim_limitIonDens", 5.0D20, "The ion density (1/m^3) from which to start limiting the Efield")
    CALL CFG_add("sim_EfieldErrThreshold", 1.0d-2, "Error threshold for electric field grid refinement")
    CALL CFG_add("sim_useTimeDepRNG", .false., "Whether the RNG is initialized using the current time")

    CALL CFG_add("ref_maxLevels", 4, "The maximum number of refinement levels in the electric field")
    CALL CFG_add("ref_min_lvl_electrode", 3, "The minimum refinement lvl around the electrode tip")
    CALL CFG_add("ref_outsideGridRatio", 4.0d0, &
         & "The ratio of the *boundary condition supplying grid* to the coarsest normal grid")
    CALL CFG_add("ref_min_grid_size", 8, "The minimum number of grid points in any direction of a refined grid")
    CALL CFG_add("ref_buffer_width", 2, "Increase refinement area by this size")
    CALL CFG_add("ref_min_grid_separation", 2, "The minimum distance between grids at the same level")
    CALL CFG_add("ref_deltaValues",  (/2.0d-6,   3.0d-6,  5.0d-6,  1.0d-5,  1.0d-4,  1.0d-3/), &
         & "List of spatial stepsizes for the electric fields in ref_EfieldValues", varSize = .true.)
    CALL CFG_add("ref_maxEfieldAtDelta", (/3.0d7,    2.0d7,   1.5d7,   1.0d7,   3.0d6,   0.0d0/), &
         & "List of electric field values at which we specify the max. spatial stepsize", varSize = .true.)
    CALL CFG_add("ref_timePerAdaptation",  2.0d-11, &
         & "Time between consecutive adaptations of the grid structure")
    CALL CFG_add("ref_MinElecDens",  1.0d15, &
         & "Minimum electron density for considering refining in a region")

    !! Ion transport parameters
    CALL CFG_add("sim_ionMobility", 2.5D-4, "The ion mobility constant")
    CALL CFG_add("sim_ionDiffusion", 5.0D-6, "The ion diffusion constant")

    !! Gas parameters
    CALL CFG_add("gas_pressure", 1.0D0, "The gas pressure (bar)")
    CALL CFG_add("gas_temperature", 293.0D0, "The gas temperature (Kelvin)")
    CALL CFG_add("gas_names", (/"N2", "O2", "Ar"/), &
         & "The names of the gases used in the simulation", varSize = .true.)
    CALL CFG_add("gas_files", (/"siglo.sec", "siglo.sec", "siglo.sec"/), &
         & "The files in which to find cross section data for each gas", varSize = .true.)
    CALL CFG_add("gas_fractions", (/0.7812d0, 0.20946d0, 0.009340d0 /), &
         & "The partial pressure of the gases (as if they were ideal gases)", varSize = .true.)
    CALL CFG_add("sim_crossSecScaleFactor", 1.0D0, "Scale factor for the cross sections in the input data")

    ! Electric field parameters
    CALL CFG_add("sim_efield_values", (/-1.0D7/), &
         "A list of values of the electric field applied in the z-direction", varSize = .true.)
    CALL CFG_add("sim_efield_times", (/0.0D0/), &
         "The times (in increasing order) at which the values of the applied electric field are used", varSize = .true.)
    CALL CFG_add("sim_constantEfield", .FALSE., &
         & "Whether the electric field is constant in space and time (no space charge effect)")

    ! Initial conditions
    CALL CFG_add("init_condType", "seed", "The type of initial condition")
    CALL CFG_add("init_seedPosRadius", 1.0D-4, "The radius of the distribution of initial ion seeds")
    CALL CFG_add("init_weightFactor", 1, "The initial weight of electrons")
    CALL CFG_add("init_nIonPairs", 1.d2 , "The number of initial ion pairs")
    CALL CFG_add("init_lineDens", 1.d10 , "The density for some types of initial conditions")
    CALL CFG_add("init_laserDirection", (/1.0d0, 0.0d0, 1.0d0/) , "Direction of the initial laserbeam")
    CALL CFG_add("init_laserLineOffset", (/0.5d0, 0.5d0, 0.6d0/) , "Laser goes through this point")
    CALL CFG_add("init_ionDistribution", 0, "The type of distribution of the initial ion pairs")
    CALL CFG_add("init_electronEnergy", 1.0D1 , "The initial energy of the electrons in eV")
    CALL CFG_add("init_relSeedPos", (/0.5D0, 0.5D0, 0.5D0/), "The relative position of the ion pairs in the domain")
    CALL CFG_add("init_backgroundDensity", 0.0D0, "The background ion/electron density in the domain")
    CALL CFG_add("init_O2MinBackgroundDens", 0.0D0, "The background O2- density in the domain")
    CALL CFG_add("sim_use_O2Min_bg", .FALSE., "Whether we use backgroundO2- in the simulation")
    CALL CFG_add("sim_ZenerInLiquidUsed", .FALSE., "Whether we use zener model in liquid!")
    CALL CFG_add("ref_epsilonRelative", 1.d0, "The relative permitivity in different media!")

    ! Poisson solver related parameters
    call CFG_add("grid_size", (/33, 33, 33/), "The number of grid cells in each direction")
    call CFG_add("grid_delta", (/5.0e-6_dp, 5.0e-6_dp, 5.0e-6_dp/), &
         "The length of a grid cell in each direction")
    call CFG_add("grid_plasma_min_rel_pos", (/0.0d0, 0.0d0, 0.0d0/), "Plasma can only be above these relative coords")
    call CFG_add("grid_plasma_max_rel_pos", (/1.0d0, 1.0d0, 1.0d0/), "Plasma can only be below these relative coords")


    ! Particle model related parameters
    CALL CFG_add("part_max_num", 5000*1000, "The maximum number of particles allowed per task")
    CALL CFG_add("part_n_rescale", 1000*1000, "The number of particles at which to start rescaling")
    CALL CFG_add("part_max_eV", 1000.0D0, "The maximum energy in eV for particles in the simulation")
    CALL CFG_add("part_lkp_tbl_size", 10000, "The table size to use in the particle model for lookup tables")
    CALL CFG_add("part_minMergeDens", 0.0d16, "The minimum electron density to consider merging particles")
    CALL CFG_add("part_per_cell", 64, &
         & "The minimum number of particles in a cell to consider rescaling them")
    CALL CFG_add("part_max_weight", 1000, &
         & "Only particles with a weight < maxWeight will be rescaled")
    CALL CFG_add("part_max_merge_dist", 1.0d0, &
         & "Only particles closer together than this var * delta can be merged")

    ! Photoionization parameters
    CALL CFG_add("PI_enabled", .FALSE., "Whether photoionization is used")
    CALL CFG_add("PI_EfieldTable", (/ 0.0D0, 0.25D7, 0.4D7, 0.75D7, 1.5D7, 3.75D7/), &
         & "The tabulated values of the electric field (for the photo-efficiency)")
    CALL CFG_add("PI_efficiencyTable", (/0.05D0, 0.05D0, 0.12D0, 0.08D0, 0.06D0, 0.4D0/), &
         & "The tabulated values of the the photo-efficiency")
    CALL CFG_add("PI_frequencies", (/2.925D12, 3.059D12/), &
         & "The lower/upper bound for the frequency of the photons, not currently used")
    CALL CFG_add("PI_absorpInvLengths", (/3.5D0 / UC_torr_to_bar, 200D0 / UC_torr_to_bar/), &
         & "The inverse min/max absorption length, have to be scaled by oxygen partial pressure")
    CALL CFG_add("PI_meanLifeTimeExcited", 0.1d-10, "The mean life time of excited species")

    ! Output parameters
    CALL CFG_add("output_timestep", 1.0D-11, "The timestep for writing output")
    CALL CFG_add("output_nIntervalsEEDF", 200, "The number of bins for writing the electron energy dist. function")
    CALL CFG_add("output_eVmaxEEDF", 9.9D2, "The maximum energy of the electron energy dist. function")

    ! Electrode parameters
    CALL CFG_add("elec_boundaryType", 1, "Boundary condition, 1: linear increase along sides. 2: grounded sides.")
    CALL CFG_add("elec_times", (/0.0d-9, 1.0d-9, 2.0d-9/), "The times at which the electrode voltage is specified.", .true.)
    CALL CFG_add("elec_voltages", (/0.0d3, 1.0d3, 5.0d3/), "The voltages at the specified times.", .true.)
    CALL CFG_add("electrode_Hcone", 0.25D-3, "Height of conical part of electrode")
    CALL CFG_add("electrode_Rcyl", 2.5D-4, "Radius of cylindrical part of electrode")
    CALL CFG_add("electrode_RcTip", 1.0D-4, "Radius of curvature of the tip of the electrode")
    CALL CFG_add("electrode_voltage", 1.0d3, "The voltage on the electrode")
    CALL CFG_add("electrode_xyzRelPos", (/0.5d0, 0.5d0, 0.25d0/), "The relative position of the electrode in the domain")
    CALL CFG_add("electrode_RcTrans", 0.0D-5, &
         & "Radius of curvature of the transition from the cylindrical to conical region of the electrode")
    CALL CFG_add("electrode_spacing", 2.0d-6, "The distance between points on the electrode surface")
    CALL CFG_add("electrode_nSurfacePoints", 100, &
         & "The number of points on the electrode surface that are used to create the CSM linear system of equations")
    CALL CFG_add("electrode_surfacePointsScale", 1.0D0, &
         & "The scale factor for the distribution of surface points on the electrode.")

    CALL CFG_add("sim_addBGElectrons", .FALSE., &
         & "Whether we need to add background electrons when they fly out)")

    CALL CFG_add("sim_peroidBundUsed", .FALSE., &
         & "Whether we use peroid boundary condition in z direction when electron flys out)")

    CALL CFG_add("sim_algorithmMfield", 1, "Choose the algorithm of updating particle with Magnetic field. &
         & 1: Boris rotation method; 2: Sasa's analytical solution")

    call CFG_sort()
  end subroutine create_config

end program pamdi3d
