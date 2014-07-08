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

!> Creates the parameters used in the program and gives them a default value
MODULE module_parameters
   USE module_config
   USE module_constants

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: createParameters

CONTAINS
   !> Create all the parameters used in the simulation and initialize them to their
   !! default values
   subroutine createParameters()
      INTEGER :: tempInt(3)
      DOUBLE PRECISION :: tempFloat(3)

      ! Set default values for the parameters
      ! General simulation parameters
      CALL CFG_addVar("sim_endTime", 5.0D-9, "The desired endtime in seconds of the simulation")
      CALL CFG_addVar("sim_minTimestep", 5.0D-15, "The minimal timestep in seconds")
      CALL CFG_addVar("sim_maxTimestep", 1.0D-11, "The maximal timestep in seconds")
      CALL CFG_addVar("sim_useElectrode", .false., "Whether to use an electrode in the simulation")
      CALL CFG_addVar("sim_name", "mcstr", "The name of the simulation")
      CALL CFG_addVar("sim_nSamplesEfield", 1000, "The number of samples to take to estimate the electric field error")
      CALL CFG_addVar("sim_maxEfieldErr", 1.0D-2, "The maximum allowed error in the electric field before recomputing it")
      CALL CFG_addVar("sim_maxElectrons", 1.0d99, "Stop the simulation when this many electrons are in the simulation")
      CALL CFG_addVar("sim_maxInterpErr", 1.0D-2, &
         & "The maximum allowed error in the acceleration of particles, before interpolating a new value from the Efield")
      CALL CFG_addVar("sim_minIncreaseRescale", 1.5, &
         & "The minimum increase in number of particles before their weights are rescaled (superparticle creation)")
      CALL CFG_addVar("sim_nRunsMax", 1, &
         & "The number of runs to perform for simulation types that require this")
      CALL CFG_addVar("sim_limitIonDens", 5.0D20, "The ion density (1/m^3) from which to start limiting the Efield")
      CALL CFG_addVar("sim_EfieldErrThreshold", 1.0d-2, "Error threshold for electric field grid refinement")
      CALL CFG_addVar("sim_useTimeDepRNG", .false., "Whether the RNG is initialized using the current time")

      CALL CFG_addVar("ref_maxLevels", 4, "The maximum number of refinement levels in the electric field")
      CALL CFG_addVar("ref_min_lvl_electrode", 3, "The minimum refinement lvl around the electrode tip")
      CALL CFG_addVar("ref_outsideGridRatio", 4.0d0, &
         & "The ratio of the *boundary condition supplying grid* to the coarsest normal grid")
      CALL CFG_addVar("ref_min_grid_size", 8, "The minimum number of grid points in any direction of a refined grid")
      CALL CFG_addVar("ref_buffer_width", 2, "Increase refinement area by this size")
      CALL CFG_addVar("ref_min_grid_separation", 2, "The minimum distance between grids at the same level")
      CALL CFG_addVar("ref_deltaValues",  (/2.0d-6,   3.0d-6,  5.0d-6,  1.0d-5,  1.0d-4,  1.0d-3/), &
         & "List of spatial stepsizes for the electric fields in ref_EfieldValues", varSize = .true.)
      CALL CFG_addVar("ref_maxEfieldAtDelta", (/3.0d7,    2.0d7,   1.5d7,   1.0d7,   3.0d6,   0.0d0/), &
         & "List of electric field values at which we specify the max. spatial stepsize", varSize = .true.)
      CALL CFG_addVar("ref_timePerAdaptation",  2.0d-11, &
         & "Time between consecutive adaptations of the grid structure")
      CALL CFG_addVar("ref_MinElecDens",  1.0d15, &
         & "Minimum electron density for considering refining in a region")

      !! Ion transport parameters
      CALL CFG_addVar("sim_ionMobility", 2.5D-4, "The ion mobility constant")
      CALL CFG_addVar("sim_ionDiffusion", 5.0D-6, "The ion diffusion constant")

      !! Gas parameters
      CALL CFG_addVar("sim_gasPressure", 1.0D0, "The gas pressure (bar)")
      CALL CFG_addVar("sim_gasTemperature", 293.0D0, "The gas temperature (Kelvin)")
      CALL CFG_addVar("sim_gasNames", (/"N2", "O2", "Ar"/), &
         & "The names of the gases used in the simulation", varSize = .true.)
      CALL CFG_addVar("sim_gasFiles", (/"siglo.sec", "siglo.sec", "siglo.sec"/), &
         & "The files in which to find cross section data for each gas", varSize = .true.)
      CALL CFG_addVar("sim_gasFractions", (/0.7812d0, 0.20946d0, 0.009340d0 /), &
         & "The partial pressure of the gases (as if they were ideal gases)", varSize = .true.)
      CALL CFG_addVar("sim_crossSecScaleFactor", 1.0D0, "Scale factor for the cross sections in the input data")

      ! Electric field parameters
      CALL CFG_addVar("sim_efield_values", (/-1.0D7/), &
         "A list of values of the electric field applied in the z-direction", varSize = .true.)
      CALL CFG_addVar("sim_efield_times", (/0.0D0/), &
         "The times (in increasing order) at which the values of the applied electric field are used", varSize = .true.)
      CALL CFG_addVar("sim_constantEfield", .FALSE., &
         & "Whether the electric field is constant in space and time (no space charge effect)")

      ! Initial conditions
      CALL CFG_addVar("init_condType", "seed", "The type of initial condition")
      CALL CFG_addVar("init_seedPosRadius", 1.0D-4, "The radius of the distribution of initial ion seeds")
      CALL CFG_addVar("init_weightFactor", 1, "The initial weight of electrons")
      CALL CFG_addVar("init_nIonPairs", 1.d2 , "The number of initial ion pairs")
      CALL CFG_addVar("init_lineDens", 1.d10 , "The density for some types of initial conditions")
      CALL CFG_addVar("init_laserDirection", (/1.0d0, 0.0d0, 1.0d0/) , "Direction of the initial laserbeam")
      CALL CFG_addVar("init_laserLineOffset", (/0.5d0, 0.5d0, 0.6d0/) , "Laser goes through this point")
      CALL CFG_addVar("init_ionDistribution", 0, "The type of distribution of the initial ion pairs")
      CALL CFG_addVar("init_electronEnergy", 1.0D1 , "The initial energy of the electrons in eV")
      CALL CFG_addVar("init_relSeedPos", (/0.5D0, 0.5D0, 0.5D0/), "The relative position of the ion pairs in the domain")
      CALL CFG_addVar("init_backgroundDensity", 0.0D0, "The background ion/electron density in the domain")
      CALL CFG_addVar("init_O2MinBackgroundDens", 0.0D0, "The background O2- density in the domain")
      CALL CFG_addVar("sim_use_O2Min_bg", .FALSE., "Whether we use backgroundO2- in the simulation")
      CALL CFG_addVar("sim_ZenerInLiquidUsed", .FALSE., "Whether we use zener model in liquid!")
      CALL CFG_addVar("ref_epsilonRelative", 1.d0, "The relative permitivity in different media!")

      ! Poisson solver related parameters

      tempInt = (/128, 128, 256/)
      tempFloat = electronMFP

      CALL CFG_addVar("grid_size", tempInt, "The number of grid cells in each direction")
      CALL CFG_addVar("grid_delta", tempFloat, "The length of a grid cell in each direction")
      call CFG_addVar("grid_plasma_min_rel_pos", (/0.0d0, 0.0d0, 0.0d0/), "Plasma can only be above these relative coords")
      call CFG_addVar("grid_plasma_max_rel_pos", (/1.0d0, 1.0d0, 1.0d0/), "Plasma can only be below these relative coords")


      ! Particle model related parameters
      CALL CFG_addVar("part_nParticlesMax", 5000*1000, "The maximum number of particles allowed per task")
      CALL CFG_addVar("part_nParticlesRescale", 1000*1000, "The number of particles at which to start rescaling")
      CALL CFG_addVar("part_weightingScheme", 1, "The order of the weighting scheme used to map particles to densities")
      CALL CFG_addVar("part_maxEnergyEv", 1000.0D0, "The maximum energy in eV for particles in the simulation")
      CALL CFG_addVar("part_tableSize", 10000, "The table size to use in the particle model for lookup tables")
      CALL CFG_addVar("part_minMergeDens", 0.0d16, "The minimum electron density to consider merging particles")
      CALL CFG_addVar("part_minPartForMeanEnergy", 10, "The minimum number of particles in a cell to output their mean energy")
      CALL CFG_addVar("part_minPartPerCell", 64, &
         & "The minimum number of particles in a cell to consider rescaling them")
      CALL CFG_addVar("part_rescaleAccelLimit", 10.0D0, &
         & "Only particles with an acceleration < part_rescaleAccelLimit * meanAccel will be rescaled")
      CALL CFG_addVar("part_rescaleVelDiff", 0.05D0, &
         & "Two particles can only be combined if their difference in velocity is smaller than rescaleVelDiff * meanVelocity")
      CALL CFG_addVar("part_maxWeight", 1000, &
         & "Only particles with a weight < maxWeight will be rescaled")
      CALL CFG_addVar("part_maxRelDistForMerge", 1.0d0, &
         & "Only particles closer together than this var * delta can be merged")

      ! KISS seeding
      CALL CFG_addVar("KISS_seed", (/123, 1234, 12345, 123456/), "Seeds for the KISS RNG")

      ! Photoionization parameters
      CALL CFG_addVar("PI_enabled", .FALSE., "Whether photoionization is used")
      CALL CFG_addVar("PI_EfieldTable", (/ 0.0D0, 0.25D7, 0.4D7, 0.75D7, 1.5D7, 3.75D7/), &
         & "The tabulated values of the electric field (for the photo-efficiency)")
      CALL CFG_addVar("PI_efficiencyTable", (/0.05D0, 0.05D0, 0.12D0, 0.08D0, 0.06D0, 0.4D0/), &
         & "The tabulated values of the the photo-efficiency")
      CALL CFG_addVar("PI_frequencies", (/2.925D12, 3.059D12/), &
         & "The lower/upper bound for the frequency of the photons, not currently used")
      CALL CFG_addVar("PI_absorpInvLengths", (/3.5D0 / TorrToBar, 200D0 / TorrToBar/), &
         & "The inverse min/max absorption length, have to be scaled by oxygen partial pressure")
      CALL CFG_addVar("PI_meanLifeTimeExcited", 0.1d-10, "The mean life time of excited species")

      ! Output parameters
      CALL CFG_addVar("output_timestep", 1.0D-11, "The timestep for writing output")
      CALL CFG_addVar("output_nIntervalsEEDF", 200, "The number of bins for writing the electron energy dist. function")
      CALL CFG_addVar("output_eVmaxEEDF", 9.9D2, "The maximum energy of the electron energy dist. function")

      ! Electrode parameters
      CALL CFG_addVar("elec_boundaryType", 1, "Boundary condition, 1: linear increase along sides. 2: grounded sides.")
      CALL CFG_addVar("elec_times", (/0.0d-9, 1.0d-9, 2.0d-9/), "The times at which the electrode voltage is specified.", .true.)
      CALL CFG_addVar("elec_voltages", (/0.0d3, 1.0d3, 5.0d3/), "The voltages at the specified times.", .true.)
      CALL CFG_addVar("electrode_Hcone", 0.25D-3, "Height of conical part of electrode")
      CALL CFG_addVar("electrode_Rcyl", 2.5D-4, "Radius of cylindrical part of electrode")
      CALL CFG_addVar("electrode_RcTip", 1.0D-4, "Radius of curvature of the tip of the electrode")
      CALL CFG_addVar("electrode_voltage", 1.0d3, "The voltage on the electrode")
      CALL CFG_addVar("electrode_xyzRelPos", (/0.5d0, 0.5d0, 0.25d0/), "The relative position of the electrode in the domain")
      CALL CFG_addVar("electrode_RcTrans", 0.0D-5, &
         & "Radius of curvature of the transition from the cylindrical to conical region of the electrode")
      CALL CFG_addVar("electrode_spacing", 2.0d-6, "The distance between points on the electrode surface")
      CALL CFG_addVar("electrode_nSurfacePoints", 100, &
         & "The number of points on the electrode surface that are used to create the CSM linear system of equations")
      CALL CFG_addVar("electrode_surfacePointsScale", 1.0D0, &
         & "The scale factor for the distribution of surface points on the electrode.")

      CALL CFG_addVar("sim_addBGElectrons", .FALSE., &
         & "Whether we need to add background electrons when they fly out)")

      CALL CFG_addVar("sim_peroidBundUsed", .FALSE., &
         & "Whether we use peroid boundary condition in z direction when electron flys out)")

      CALL CFG_addVar("sim_algorithmMfield", 1, "Choose the algorithm of updating particle with Magnetic field. &
         & 1: Boris rotation method; 2: Sasa's analytical solution")

    end subroutine createParameters

END MODULE module_parameters
