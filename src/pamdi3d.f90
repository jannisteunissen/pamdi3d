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

! This is pamdi3d version 0.01 (July 7th, 2014)
program pamdi3d

   use module_crossSections
   use module_initialCond
   use module_particle
   use m_efield_amr
   use module_electrode
   use module_parameters
   use module_config
   use module_globals
   use module_constants
   use module_photoionization
   use module_IO
   use module_electrode
   use module_kiss
   use module_gas
   use MPI

   implicit none
   character(LEN=80) :: simName, configName

   integer  :: timeStep, stepsLeftPoisson
   integer  :: nOutput, nGridAdapts
   integer  :: WCTime_start, WCTime_current
   integer  :: nParticlesRescale
   integer  :: nPartSum, nPartSumPrevRescale
   integer  :: KISSseed(4)
   integer  :: ierr, myrank, ntasks, root = 0, ix, n_its
   logical  :: finished, flagOutput, rescale, adaptGrid, &
        flagAddBG, useBGO2Minus, doFakeAttach, useElectrode

   double precision :: sim_time, endTime
   double precision :: timePerOutput, timePerGridAdapt
   double precision :: tauPoisson, minIncreaseRescale, prevPoissonTime
   double precision :: interpErr, PoissonErr, rescalePSERR
   double precision :: maxEfieldErr, maxInterpErr, maxPSErr
   double precision :: minTimestep, maxTimestep, tau, nextTau
   double precision :: nElecDiff, nElecSum, nElecSumPrev
   double precision :: meanWeight
   double precision :: maxElectrons
   double precision :: densRateAdded   ! desnity generate rate for background density

   finished             = .false.
   flagOutput           = .false.
   rescale              = .false.
   adaptGrid            = .false.

   nOutput              = 0
   nGridAdapts          = 0
   timeStep             = 0
   stepsLeftPoisson     = 0

   sim_time       = 0.0D0
   prevPoissonTime      = 0.0D0
   maxPSErr             = 0.0D0
   interpErr            = 0.0D0
   PoissonErr           = 0.0D0

   ! Set up MPI environment
   call MPI_INIT(ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)

   ! Create an empty configuration of initial size 100 (will increase automatically)
   call CFG_initialize(200)

   ! Root process reads in config file
   if (myrank == root) then
      print *, "Starting simulation with ", ntasks, "tasks"
      print *, " ~~~ initializing parameters and updating them from config file"
      call createParameters()

      do ix = 1, command_argument_count()
         call get_command_argument(ix, configName)
         call CFG_updateFromFile( trim(configName) )
      end do

      call CFG_sort()
   end if

   ! Now share the parameters with all other processes
   call CFG_shareMPI(myrank, root)

   ! Read in the crossection data from files
   call GAS_initialize()
   call CS_initialize()

   call CFG_getVar("sim_name", simName)

   ! The root process reads in cross section data
   if (myrank == root) then
      print *, "Reading crossection data"
      call CS_readFromFile()
      call CS_writeProcesses( "output/" // trim(simName) // "_collisions.txt" )
      call CFG_write( "output/" // trim(simName) // "_config.txt" )
   end if

   call CS_shareMPI(myrank, root)

   ! Set the global variables according to the parameters
   print *, " ~~~ initializing global variables module [to be removed]"
   call GL_initializeVariables(myrank, root)

   useElectrode = CFG_varLogic("sim_useElectrode")

   ! Initialize all modules that need initialization
   print *, " ~~~ initializing electric field module"
   call E_initialize()
   print *, " ~~~ initializing electrode module"
   call EL_initialize(myrank, root)
   print *, " ~~~ initializing particle module"
   call initializeParticle()

   if (GL_flagPhotoIonization) then
      print *, " ~~~ initializing photoionization module"
      call MISC_initialize()
   end if

   if (CFG_varLogic("sim_useTimeDepRNG")) then
      call itime(KISSseed)
      KISSseed(4) = 5346 + product(KISSseed(1:3))
      print *, "RNG Seed:", KISSseed
   else
      call CFG_getVar("KISS_seed", KISSseed)
   end if

   call initializeKISS(KISSseed, myrank)

   if (myrank == root) then
      call system_clock(COUNT=WCTime_start)
      print *, 'Setting up initial electron/ion positions and Efield'
      call setInitialConditions()
      call attachmentRateCounter(densRateAdded)
   end if

   if (myrank == root) print *, "Sharing initial particles"
   call shareParticles(myrank, ntasks)

   if (myrank == root) print *, "Computing initial field"
   call PM_particlesToDensity()
   call E_compute_field(myrank, root, sim_time)

   call E_update_grids(myrank, root, sim_time)
   call PM_particlesToDensity()
   if (useElectrode) call E_readjust_electrode_new_grid(myrank, root, &
        sim_time, n_its)
   call E_compute_field(myrank, root, sim_time)

   if (useElectrode) then
      do ix = 1, 10
         call E_compute_field(myrank, root, sim_time)
      end do

      call E_update_grids(myrank, root, sim_time)
      call PM_particlesToDensity()
      call E_readjust_electrode_new_grid(myrank, root, sim_time, n_its)

      call E_compute_field(myrank, root, sim_time)
      do ix = 1, 10
         call E_compute_field(myrank, root, sim_time)
      end do
   end if

   if (myrank == root) call E_benchmark()

   call E_update_grids(myrank, root, sim_time)
   call PM_particlesToDensity()
   if (useElectrode) call E_readjust_electrode_new_grid(myrank, root, &
        sim_time, n_its)

   call E_compute_field(myrank, root, sim_time)
   call updateParticleAccel()

   ! Initial values
   nParticlesRescale    = CFG_varInt("part_nParticlesRescale")
   endTime              = CFG_varDble("sim_endTime")
   timePerOutput        = CFG_varDble("output_timestep")
   timePerGridAdapt     = CFG_varDble("ref_timePerAdaptation")
   minTimestep          = CFG_varDble("sim_minTimestep")
   maxTimestep          = CFG_varDble("sim_maxTimestep")
   maxEfieldErr         = CFG_varDble("sim_maxEfieldErr")
   maxInterpErr         = CFG_varDble("sim_maxInterpErr")
   minIncreaseRescale   = CFG_varDble("sim_minIncreaseRescale")
   maxElectrons         = CFG_varDble("sim_maxElectrons")
   useBGO2Minus         = CFG_varLogic("sim_use_O2Min_bg")
   flagAddBG            = CFG_varLogic("sim_addBGElectrons")

   nextTau              = minTimestep
   tauPoisson           = minTimestep
   nElecSumPrev         = PM_get_num_elec_mpi()
   nPartSumPrevRescale  = PM_get_num_part_mpi()
   doFakeAttach         = .true.

   call E_benchmark()

   if (myrank == root) then
      print *, "Starting simulation with name ", trim(simName), ","
      print *, "will simulate up to ", endTime, " seconds"
   end if

   ! Write the output
   nOutput = 1
   call writeOutput(nOutput, 0.0d0, trim(simName), myrank, root)

   do while (.not. finished)

      timeStep          = timeStep + 1
      tau               = nextTau
      stepsLeftPoisson  = stepsLeftPoisson - 1
      sim_time    = sim_time + tau

      call updateParticlesMC(tau)

      nElecSum    = PM_get_num_elec_mpi()
      nElecDiff   = abs((nElecSum - nElecSumPrev) / (nElecSumPrev + smallNumber))

      if (myrank == root) then
         if (nElecDiff > 2.0d0) stepsLeftPoisson = 0
      end if
      call MPI_BCAST(stepsLeftPoisson, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)

      ! No need to update anything until we have particles, also don't write output
      if (nElecSum == 0) then
         do
            if (myrank == root) print *, "No particles, timestep: ", timestep, "time:", sim_time
            if(GL_flagPhotoIonization) call MISC_photoionization(sim_time - prevPoissonTime, myrank, root)
            if (useBGO2Minus) call MISC_detachment(sim_time - prevPoissonTime, myrank, root)
            call E_compute_field(myrank, root, sim_time)

            timeStep          = timeStep + 1
            tau               = maxTimeStep
            prevPoissonTime   = sim_time
            sim_time    = sim_time + tau
            nElecSum          = PM_get_num_elec_mpi()

            if (myrank == root) finished = sim_time > endTime
            call MPI_BCAST(finished, 1, MPI_LOGICAL, root, MPI_COMM_WORLD, ierr)
            if (nElecSum > 0 .or. finished) exit
         end do

         nOutput = int(sim_time / timePerOutput) ! Since the output was skipped, advance the counter.
         cycle
      end if


      if (stepsLeftPoisson <= 0) then
         ! Do an update of the electric field and the stepsizes,
         ! and possibly write output
         nPartSum    = PM_get_num_part_mpi()
         meanWeight  = nElecSum / (nPartSum + smallNumber)

         nElecSumPrev = nElecSum

         call checkinterpErrorMPI(interpErr, tau, nPartSum)
         call storeSampleOfEfield()

         ! Check for photoionization events
         if(GL_flagPhotoIonization) call MISC_photoionization(sim_time - prevPoissonTime, myrank, root)

         ! check for negative background ions events
         if (useBGO2Minus) call MISC_detachment(sim_time - prevPoissonTime, myrank, root)

         call PM_particlesToDensity()
         call E_compute_field(myrank, root, sim_time)

         call updateParticleAccel()
         call checkPoissonErrorMPI(PoissonErr, nPartSum)
         maxPSErr = max(maxPSErr, PoissonErr)

         ! Redistribute the particles
         call shareParticles(myrank, ntasks)

         if (myrank == root) then
            nextTau           = getNewTau(tau, interpErr, maxInterpErr, nElecDiff)
            tauPoisson        = getNewTau(tauPoisson, PoissonErr, maxEfieldErr, nElecDiff)
            nextTau           = min(nextTau, tauPoisson) ! Monte Carlo step can not be larger than Poisson step
            stepsLeftPoisson  = floor(tauPoisson / nextTau)

            ! Centralize these flags because they might differ on heterogeneous systems
            finished          = sim_time > endTime .or. nElecSum > maxElectrons
            flagOutput        = nOutput <= sim_time / timePerOutput .or. finished
            rescale           = nPartSum > nParticlesRescale .and. &
                 nPartSum > int(nPartSumPrevRescale * minIncreaseRescale)
            rescale           = rescale .or. timeStep == 1
            adaptGrid         = nGridAdapts <= sim_time / timePerGridAdapt
         end if

         ! Communicate the new step sizes and flags
         call MPI_BCAST(nextTau, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
         call MPI_BCAST(stepsLeftPoisson, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
         call MPI_BCAST(flagOutput, 1, MPI_LOGICAL, root, MPI_COMM_WORLD, ierr)
         call MPI_BCAST(adaptGrid,  1, MPI_LOGICAL, root, MPI_COMM_WORLD, ierr)
         call MPI_BCAST(finished,   1, MPI_LOGICAL, root, MPI_COMM_WORLD, ierr)
         call MPI_BCAST(rescale,    1, MPI_LOGICAL, root, MPI_COMM_WORLD, ierr)

         ! Write the output
         if (flagOutput) then
            nOutput = nOutput + 1
            call writeOutput(nOutput, sim_time, trim(simName), myrank, root)
         end if

         if (adaptGrid) then
            nGridAdapts = nGridAdapts + 1
            call PM_particlesToDensity()
            call E_update_grids(myrank, root, sim_time)
            call PM_particlesToDensity()
            if (useElectrode) call E_readjust_electrode_new_grid(myrank, root, &
                 sim_time, n_its)
            nextTau = nextTau / max(n_its-1, 1) ! For safety, new mesh can induce big currents
            call E_compute_field(myrank, root, sim_time)
         end if

         if (adaptGrid .or. rescale) then
            call E_share_vars((/E_i_elec/), root)

            if (doFakeAttach) then
              if (myrank == root) print *, "Performing fake attachments!"
              call PM_convertElecForStability(tauPoisson, 0.2d0) ! Second arg is max mobility estimate
           end if

            ! TODO: The above call modifies the electron density slightly, probably not very important but check it later.
            call PM_divideParticles(myrank, ntasks)

            if (myrank == root) then
               print *, timestep, "Rescaling particles, meanWeight:", meanWeight, "nParticles", nPartSum
            end if

            call storeSampleOfEfield()

            call PM_mergeAndSplit(myrank)

            nPartSumPrevRescale  = PM_get_num_part_mpi()
            nPartSum             = nPartSumPrevRescale

            call PM_particlesToDensity()
            call E_compute_field(myrank, root, sim_time)
            call checkPoissonErrorMPI(rescalePSERR, nPartSum)

            if (myrank == root) then
               print *, "Poisson error due to rescaling:", rescalePSERR
            end if

            call updateParticleAccel()

            ! Redistribute the particles
            call shareParticles(myrank, ntasks)
         end if

         if (myrank == root) then
            call system_clock(COUNT=WCTime_current)
            write(*, FMT = "((I7),(A),(es10.4e2),(A),(es10.2e2),(A),(es10.2e2))") nOutput, &
                 & "   t(s): ", sim_time, "   tau: ", nextTau, "   tauPS", tauPoisson
            write(*, FMT = "((A)(es10.2e2)(A),(I12),(A),(es10.2e2),(A),(es10.4e2))") "       WC t(s): ", &
                 & (WCTime_current - WCTime_start)/1.0D3, "   nPart:", nPartSum, &
                 & "  nElec: ", nElecSum
            ! print *, "InterpErr:", interpErr, "PoissonErr:", PoissonErr
            ! print *, "maxPSerr:", maxPSerr
            ! write(*,*) ""
         end if
         prevPoissonTime = sim_time
      end if

      ! The end of the simulation main loop
   end do

   call MPI_FINALIZE(ierr)

   ! The simulation is finished
   print *, "The simulation has finished succesfully"

contains

   !> Given the old stepsize 'oldTau', the error 'err', the maximum allowed error 'maxErr',
   !! and the relative change in number of particles, return a new stepsize.
   double precision function getNewTau(oldTau, err, maxErr, NRDiff)
      double precision, intent(IN) :: oldTau       ! The old timestep
      double precision, intent(IN) :: err, maxErr  ! The actual error and the maximum / prev. error
      double precision, intent(IN) :: NRDiff       ! The relative difference in number of particles between timesteps

      print *, "err", err
      ! Adjust the timesteps to obtain an error close to the desired one.
      ! The interpolation error should be ~ interpErrFrac * maxEfieldErr
      if (err < 0.5D0 * maxErr) then
         getNewTau = min(oldTau * (maxErr/(err+epsilon(err)))**0.1D0, 2.0D0 * oldTau)
      else if (err > maxErr) then
         getNewTau = oldTau * (maxErr/err)
      else
         getNewTau = oldTau
      end if

      getNewTau = min(max(getNewTau, minTimeStep), maxTimeStep)

   end function getNewTau

end program pamdi3d
