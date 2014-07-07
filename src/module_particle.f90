!> This module contains all of the Monte Carlo particle model
!!
!! In the simulation one will typically call updateParticlesMC(),
!! which starts
module module_particle

   use module_globals
   use module_constants
   use generalUtilities
   use module_kiss
   use module_config
   use m_efield_amr
   use module_crossSections
   use module_electrode
   use module_gas

   implicit none
   private

   integer, parameter :: dp = kind(0.0d0)

   !> The number of active particles in the simulation
   integer :: nParticles

   !> The electron type
   type Electron
      real(dp) :: T_left   ! The time until the next timestep
      real(dp) :: v_prev   ! The velocity at the last collision
      real(dp) :: x(3)     ! The position
      real(dp) :: v(3)     ! The velocity
      real(dp) :: accel(3) ! The electric field's acceleration
      integer  :: weight   ! The weight factor of the (super)particle
      logical  :: live     ! If .false. the particle will not be used (anymore)
   end type Electron

   !> MPI type for sending the 'electrons' between tasks.
   integer :: partTypeMPI

   !> How many samples to take from the electric field to estimate the error
   integer :: nSamplesEfield

   integer :: minPartForMeanEnergy

   !> Whether the simulation uses an electrode
   logical :: useElectrode

   !> The maximal total collision frequency
   real(dp) :: maxCollRate

   !> The inverse of maxCollRate, for efficiency
   real(dp) :: invMaxCollRate

   !> The minimum number of particles in a cell to be viable for rescaling (creating superparticles)
   integer :: PM_min_part_per_cell

   !> The maximum weight of particles to be candidates for rescaling
   integer :: PM_max_weight

   integer :: PM_maxNumberOfParticles

   real(dp) :: PM_minMergeDens

   real(dp) :: PM_maxRelDistForMerge
   real(dp) :: PM_maxPosIonDens

   real(dp) :: PM_part_min_xyz(3), PM_part_max_xyz(3)

   !> The list that contains all the particles in the simulation
   type(Electron), dimension(:), allocatable    :: pList

   !> A list that is used to share the cross sections at a given energy between routines
   real(dp), dimension(:), allocatable  :: colRateList

   !! To collect statistics during runtime
   integer(KIND=iKind18), dimension(:), allocatable   :: nCollisionSum
   real(dp), dimension(:,:), allocatable      :: samplePos, sampleEfield
   real(dp), dimension(:,:), allocatable      :: ElectronEDF

   ! Currently there are three mechanisms for loss of energy in the electrons:
   ! excitations, ionizations and attachment
   real(dp), dimension(3)                     :: energyLoss

   !! Lookup table variables
   integer                               :: tableSize
   real(dp)                              :: maxEnergyEv
   real(dp)                              :: maxVelocity, invMaxVelocity
   real(dp), dimension(:,:), allocatable :: colRateTable
   real(dp), dimension(:), allocatable   :: colRateTableVelocity
   real(dp), dimension(:), allocatable   :: colRateTableSum

   !initialize the maegetic field(uniform)
   real(dp), dimension(3) :: Mfield
   real(dp) :: cyclotronFreq

   double precision :: sinAngle, cosAngle, angle
   integer :: updateParticlesAlgotihmWithB

   !! Public functions
   public :: isOutOfGas
   public :: updateParticlesMC
   public :: correctParticleVelocities
   public :: removeDeadParticles
   public :: removeAllParticles

   public :: createIonPair, PM_createElectron
   public :: initializeParticle
   public :: PM_particlesToDensity
   public :: get_num_part, get_num_elec, getMeanEnergy, getMeanWeight
   public :: PM_get_num_elec_mpi, PM_get_num_part_mpi
   public :: resetNCollisions, getNCollisionsOfProcess, getNCollisionsOfType
   public :: storeSampleOfEfield
   public :: checkinterpErrorMPI, checkPoissonErrorMPI
   public :: updateParticleAccel
   public :: shareParticles
   public :: setEEDF_MPI, setEEDF
   public :: PM_divideParticles
   public :: PM_mergeAndSplit
   public :: PM_convertElecForStability

contains
   !> Initialization routine for the particle module, sets private members according
   !! to the configuration, and initializes arrays after allocation. Also sets up
   !! a lookup table for the cross section data, for fast lookup.
   subroutine initializeParticle()
      nSamplesEfield          = CFG_varInt("sim_nSamplesEfield")
      PM_maxNumberOfParticles = CFG_varInt("part_nParticlesMax")
      tableSize               = CFG_varInt("part_tableSize")
      minPartForMeanEnergy    = CFG_varInt("part_minPartForMeanEnergy")
      maxEnergyEv             = CFG_varDble("part_maxEnergyEv")
      PM_maxRelDistForMerge   = CFG_varDble("part_maxRelDistForMerge")
      maxVelocity             = en2v(maxEnergyEv*electronVolt)
      invMaxVelocity          = 1.0D0 / maxVelocity
      useElectrode            = CFG_varLogic("sim_useElectrode")
      PM_minMergeDens         = CFG_varDble("part_minMergeDens")
      PM_min_part_per_cell       = CFG_varInt("part_minPartPerCell")
      PM_max_weight            = CFG_varInt("part_maxWeight")

      call CFG_getVar("grid_plasma_min_rel_pos", PM_part_min_xyz)
      call CFG_getVar("grid_plasma_max_rel_pos", PM_part_max_xyz)
      PM_part_min_xyz = PM_part_min_xyz * GL_gridLength
      PM_part_max_xyz = PM_part_max_xyz * GL_gridLength

      allocate( pList(PM_maxNumberOfParticles) )
      allocate( samplePos(3, nSamplesEfield) )
      allocate( sampleEfield(3, nSamplesEfield) )
      allocate( colRateList(nProcesses) )
      allocate( nCollisionSum(0:nProcesses) ) ! Keep track of the total number of collisions in this array

      ! Create a lookup table for the cross sections of the collision processes
      call createColRateTable(maxVelocity, tableSize)
      call createPartTypeMPI(partTypeMPI)

      nParticles              = 0
      nCollisionSum(:)        = 0
      energyLoss(:)           = 0.0D0
      write(*, FMT="((A),(E10.4e2))") " Average collision time: ", invMaxCollRate

   end subroutine initializeParticle

   !> Loop over all the particles, and for each set the time left 'dt' for this step.
   !! Then they are feeded to the moveAndCollide() routine, which advances them in time
   !! and simulates collisions.
   subroutine updateParticlesMC(dt)
      real(dp), intent(IN)  :: dt
      integer                       :: ll
      real(dp)              :: timestep_copy

      pList(1:nParticles)%T_left = dt
      ll = 1
      do while (ll <= nParticles)
         timestep_copy = pList(ll)%T_left
         call moveAndCollide(ll, timestep_copy)
         if (pList(ll)%live) call correctParticleVelocity(ll, timestep_copy) ! Correct velocity with new acceleration
         ll = ll + 1
      end do

      call removeDeadParticles()
   end subroutine updateParticlesMC

   !> Perform a collision for an electron, either elastic, excitation, ionizationCollision,
   !! attachment or null. Also update GL_excDens according to newly created excited
   !! molecules if photoionizationCollision is enabled.
   subroutine moveAndCollide(ll, dt)
      integer, intent(IN)  :: ll
      real(dp), intent(IN) :: dt

      integer              :: cIx
      real(dp)             :: velocity, timeLeft, timeCollision

      ! Set the first collision time. Because the exponential distribution is memoryless,
      ! we don't have to keep track of the old collision time.
      timeCollision = collisionTime()
      timeLeft = dt

      pList(ll)%v_prev = norm2(pList(ll)%v)

      do while (timeCollision < timeLeft .and. pList(ll)%live)
         timeLeft = timeLeft - timeCollision

         ! Set x,v at the collision time
         call advanceParticle(ll, timeCollision)

         call remove_out_of_gas(ll)
         if (.not. pList(ll)%live) exit

         ! Take the mean of the current velocity and the velocity after the last collision
         ! as an estimate of the mean velocity
         velocity = norm2(pList(ll)%v)
         if (velocity >= maxVelocity) then
            print *, "Velocity too large"
            print *, "energy in eV: ", v2en(velocity)/electronVolt
            call E_add_to_var(E_i_nion, pList(ll)%x, dble(pList(ll)%weight))
            call killParticle(ll)
            stop
         end if

         cIx                = getCollisionIndex(velocity)
         nCollisionSum(cIx) = nCollisionSum(cIx) + pList(ll)%weight

         if (cIx > 0 ) then

            ! Perform the corresponding collision
            select case( crossSecTable(cIx)%colType )
            case (attachType)
               call attachmentCollision(ll)
            case (elasticType)
               call elasticCollision(ll, cIx, velocity)
            case (exciteType)
               call excitationCollision(ll, cIx, velocity)
            case (ionizeType)
               call ionizationCollision(ll, cIx, velocity)
            end select

         end if

         pList(ll)%v_prev = velocity
         timeCollision    = collisionTime()
      end do

      ! Update the particle position and velocity to the next timestep if there are no more
      ! collisions in the current timestep.
      call advanceParticle(ll, timeLeft)
      call remove_out_of_gas(ll)
   end subroutine moveAndCollide

   !> Checks whether pos is outside the computational domain or inside the electrode,
   !! and if so returns .TRUE., otherwise returns .FALSE.
   logical function isOutOfGas(pos)
      real(dp), dimension(3), intent(IN) :: pos

      ! Check whether any component of the position is outside of the domain
      isOutOfGas = any(pos <= PM_part_min_xyz .or. pos >= PM_part_max_xyz)

      ! Now if any component is outside the domain or inside the electrode, isOutOfGas = .TRUE.
      if (useElectrode .and. EL_insideElectrode(pos)) then
         isOutOfGas = .true.
      end if
   end function isOutOfGas

   subroutine remove_out_of_gas(ll)
      integer, intent(in) :: ll

      if (any(pList(ll)%x <= PM_part_min_xyz .or. pList(ll)%x >= PM_part_max_xyz)) then
         call killParticle(ll)
      else if (useElectrode) then
         if (EL_insideElectrode(pList(ll)%x)) then
            call killParticle(ll)
         end if
      end if
   end subroutine remove_out_of_gas

   !> Returns a sample from the exponential distribution of the collision times
   real(dp) function collisionTime()
      ! kiss_rand() is uniform on [0,1), but log(0) = nan, so we take 1 - kiss_rand()
      collisionTime = -log(1.0D0 - kiss_rand()) * invMaxCollRate
   end function collisionTime

   !> From the list crosssec(:) select the index of the process that will occur,
   !! or set colIndex = 0 if there is a null collision
   integer function getCollisionIndex(velocity)
      real(dp), intent(IN) :: velocity
      integer              :: j, tablePos
      real(dp)             :: freq, randFreq, sumCrossSec

      ! First check whether we need to do a lookup in the cross section table or that
      ! we already know that there will be a 'null' collision. Therefore, find the
      ! interval in tableSumCrossec that velocity lies in.

      tablePos = getTableIndex(velocity)

      ! Get an upper bound of the sum of the cross sections at velocity
      sumCrossSec = max(colRateTableSum(tablePos), colRateTableSum(tablePos + 1))

      ! Draw a random number
      randFreq = kiss_rand() * maxCollRate

      ! Check whether we surely have a null collision, and if so return
      if (randFreq > sumCrossSec) then
         getCollisionIndex = 0
         return
      end if

      ! Likely there will be no null collision, so randomly determine the collision process.
      ! We have to fill the colRateList array first with the appropriate values

      call lookupColRate(velocity, colRateList)

      getCollisionIndex = 0
      freq              = 0.0D0

      ! This is a very simple algorithm to convert a uniform [0,1) random number
      ! to an event. Here colRateList holds the event probabilities * maxCollRate
      do j = 1, nProcesses
         freq = freq + colRateList(j)
         if (randFreq <= freq) then
            getCollisionIndex = j
            exit
         end if
      end do

   end function getCollisionIndex

   !> Perform an elastic collision for particle 'll'
   subroutine elasticCollision(ll, collisionIx, velocity)
      integer, intent(IN)              :: ll, collisionIx
      real(dp), intent(INOUT)  :: velocity
      real(dp)                 :: mass_ratio, bg_vel(3), com_vel(3)

      mass_ratio = crossSecTable(collisionIx)%elecOverMoleculeMass
      bg_vel = 0.0_dp

      ! Compute center of mass velocity
      com_vel = (mass_ratio * pList(ll)%v + bg_vel) / (1 + mass_ratio)

      ! Scatter in center of mass coordinates
      pList(ll)%v = pList(ll)%v - com_vel
      call scatter_isotropic(ll, sqrt(sum(pList(ll)%v**2)))
      pList(ll)%v = pList(ll)%v + com_vel
   end subroutine elasticCollision

   subroutine scatter_isotropic(ll, vel_norm)
      integer, intent(in)  :: ll
      real(dp), intent(in) :: vel_norm
      real(dp), parameter  :: pi = acos(-1.0_dp)
      real(dp)             :: theta, phi

      theta             = acos(1.0_dp - 2.0_dp * kiss_rand())
      phi               = 2 * pi * kiss_rand()
      pList(ll)%v(1) = vel_norm * sin(theta) * cos(phi)
      pList(ll)%v(2) = vel_norm * sin(theta) * sin(phi)
      pList(ll)%v(3) = vel_norm * cos(theta)
   end subroutine scatter_isotropic

   !> Perform an excitation-collision for particle 'll'
   subroutine excitationCollision(ll, collisionIx, velocity)
      integer,intent(IN)               :: ll, collisionIx
      real(dp), intent(INOUT)  :: velocity
      real(dp)                 :: energy, oldEnergy

      oldEnergy = v2en(velocity)
      energy = oldEnergy - crossSecTable(collisionIx)%treshold_eV * electronVolt

      ! Because of the use of a lookup table we might get a tiny negative energy.
      if (energy < smallEnergy) energy = smallEnergy

      ! Set the new velocity
      velocity = en2v(energy)

      ! Set the new scattered velocity, which is lower than the old one.
      call scatter_isotropic(ll, velocity)

   end subroutine excitationCollision

   !> Perform an ionizing collision for particle 'll'
   subroutine ionizationCollision(ll, collisionIx, velocity)
      integer, intent(IN)              :: ll, collisionIx
      real(dp), intent(INOUT)  :: velocity

      real(dp) :: energyEv, oldEnergyEv, newVel
      real(dp) :: en_s1, en_s2

      nParticles   = nParticles + 1
      call checkNumParticles(nParticles)

      oldEnergyEv   = v2en(velocity) / electronVolt
      energyEv      = oldEnergyEv - crossSecTable(collisionIx)%treshold_eV

      ! Because of the use of a lookup table we might get a tiny negative energy
      if (energyEv < smallEnergyEv) energyEv = smallEnergyEv

      ! Equal sharing
      en_s1 = 0.5D0 * energyEv
      en_s2 = 0.5D0 * energyEv

      velocity = en2v(en_s1 * electronVolt)
      call scatter_isotropic(ll, velocity)

      newVel   = en2v(en_s2 * electronVolt)
      call scatter_isotropic(nParticles, newVel)

      pList(nParticles)%x       = pList(ll)%x
      pList(nParticles)%accel   = pList(ll)%accel
      pList(nParticles)%T_left  = pList(ll)%T_left
      pList(nParticles)%v_prev  = newVel
      pList(nParticles)%weight  = pList(ll)%weight
      pList(nParticles)%live    = .true.

      call E_add_to_var(E_i_pion, pList(nParticles)%x, dble(pList(nParticles)%weight))

      ! Generation of photons is proportional to the number of ionizations
      if (GL_flagPhotoionization) then
         call E_add_to_var(E_i_exc, pList(ll)%x, dble(pList(ll)%weight))
      end if

   end subroutine ionizationCollision

   !> Perform attachment of electron 'll'
   subroutine attachmentCollision(ll)
      integer, intent(IN) :: ll

      call killParticle(ll)

      ! Add an ion with negative charge, at this moment, all attached ions are set to O2-
      call E_add_to_var(E_i_O2m, pList(ll)%x, dble(pList(ll)%weight))

   end subroutine attachmentCollision

   !> Mark particle 'll' as inactive, so that it will be removed from the particle list
   !! at the next sweap.
   subroutine killParticle(ll)
      integer, intent(IN) :: ll

      pList(ll)%live = .false.

   end subroutine killParticle

   !> Give the particle with index ll a new velocity
   !!
   !! The original velocity is described by theta, phi with magnitude vel.
   !! The outgoing velocity is rotated by angles chi, psi relative to the incoming,
   !! The coordinate system where chi, psi are defined is obtained by first rotating
   !! around the z-axis over -phi, and then rotating around the y-axis over -theta.
   !! This way the original velocity is on the z-axis.
   subroutine scatterParticle(ll,theta,phi,chi,psi,vel)
      integer,intent(IN)          :: ll
      real(dp),intent(IN) :: theta,phi,chi,psi,vel

      real(dp)            :: costheta,cosphi,sintheta,sinphi,coschi
      real(dp)            :: cospsi,sinchi,sinpsi

      costheta = cos(theta);  sintheta = sin(theta)
      cosphi = cos(phi);      sinphi = sin(phi)
      coschi = cos(chi);      sinchi = sin(chi)
      cospsi = cos(psi);      sinpsi = sin(psi)

      ! The matrix on the right hand side is obtained by doing the following elemental rotations:
      ! Rotate around z-axis over -phi, rotate around y-axis over -theta,
      ! now vel lies on the z-axis, so only take z-component into next step
      ! rotate around y-axis over xhi, rotate around z-axis over psi
      ! Do the inverse of the first two rotations (y-axis over theta, z-axis over phi).

      pList(ll)%v(1) = vel * ( sintheta*cosphi*coschi + &
           sinchi*(costheta*cosphi*cospsi-sinphi*sinpsi) )
      pList(ll)%v(2) = vel * ( sintheta*sinphi*coschi + &
           sinchi*(costheta*sinphi*cospsi+cosphi*sinpsi) )
      pList(ll)%v(3) = vel * ( costheta*coschi - sintheta*sinchi*cospsi )

   end subroutine scatterParticle

   !> Advance the particle position and velocity over time tt
   subroutine advanceParticle(ll, tt)
      integer, intent(IN)           :: ll
      real(dp), intent(IN)  :: tt

      pList(ll)%x = pList(ll)%x + pList(ll)%v * tt + &
           0.5D0 * pList(ll)%accel * tt**2
      pList(ll)%v = pList(ll)%v + pList(ll)%accel * tt

      pList(ll)%T_left = pList(ll)%T_left - tt

   end subroutine advanceParticle

   !> Correct particle velocities for the previous timestep of 'dt'
   subroutine correctParticleVelocities(dt)
      real(dp), intent(IN) :: dt
      integer :: ll

      do ll = 1, nParticles
         call correctParticleVelocity(ll, dt)
      end do

   end subroutine correctParticleVelocities

   !> Correct particle velocities for the previous timestep of 'dt'
   !!
   !! During the timestep x,v have been advanced to:
   !! x(t+1) = x(t) + v(t)*dt + 0.5*a(t)*dt^2,
   !! v(t+1) = v(t) + a(t)*dt
   !! But the velocity at t+1 should be v(t+1) = v(t) + 0.5*(a(t) + a(t+1))*dt,
   !! to have a second order leapfrog scheme, so here we set it to that value.
   subroutine correctParticleVelocity(ll, dt)
      real(dp), intent(IN) :: dt
      integer, intent(IN) :: ll
      real(dp), dimension(3) :: newAccel

      newAccel = E_get_field(pList(ll)%x)
      newAccel          = newAccel * elecChargeOverMass
      pList(ll)%v       = pList(ll)%v + 0.5D0 * (newAccel - pList(ll)%accel) * dt
      pList(ll)%accel   = newAccel

   end subroutine correctParticleVelocity

   !> Remove all particles for which live == .FALSE. from the list, which can happen
   !! due moving outside the computational domain or attachment.
   subroutine removeDeadParticles()

      integer :: ll, ix_end

      ix_end = nParticles

      do ll = 1, nParticles

         if (.not. pList(ll)%live) then
            ! Find the first alive particle from the end of the list, and place it at index ll
            ! so that all the alive particles are placed at the beginning of pList
            do while (.not. pList(ix_end)%live .and. ix_end > ll)
               ix_end = ix_end - 1
            end do

            ! If there is no alive particle available in the range [ll+1 : nParticles] we have
            ! ll >= ix_end and we exit the routine. Else we put the alive particle at index ll
            if (ll >= ix_end) then
               nParticles = ll - 1
               exit
            else
               pList(ll)            = pList(ix_end)
               pList(ix_end)%live   = .false.
               ix_end               = ix_end - 1
            end if

         end if

      end do

   end subroutine removeDeadParticles

   !> Creates an electron at 'pos' with the given energy in 'eV', and a velocity in a
   !! random direction.
   subroutine PM_createElectron(pos, eV, weight)

      real(dp), dimension(3), intent(IN) :: pos
      real(dp), intent(IN)               :: eV
      integer, intent(IN)                        :: weight

      real(dp)     :: vel, theta, phi, accel(3)

      nParticles  = nParticles + 1
      call checkNumParticles(nParticles)
      vel         = en2v(eV * electronVolt)
      theta       = acos(1.0D0-2.0D0*kiss_rand())
      phi         = 2.0D0 * pi * kiss_rand()

      accel = E_get_field(pos)
      accel       = elecChargeOverMass * accel
      pList(nParticles)%weight   = weight

      if (weight < 1) then
         print *, "Error: trying to add electron/ion with negative weights!"
         stop
      end if

      pList(nParticles)%live     = .true.
      pList(nParticles)%T_left   = 0.0D0

      pList(nParticles)%x        = pos
      pList(nParticles)%accel    = accel

      pList(nParticles)%v(1)     = vel * sin(theta) * cos(phi)
      pList(nParticles)%v(2)     = vel * sin(theta) * sin(phi)
      pList(nParticles)%v(3)     = vel * cos(theta)
      pList(nParticles)%v_prev   = vel

   end subroutine PM_createElectron

   !> Creates an electron at 'pos' with the given energy in 'eV', and a velocity in a
   !! random direction. Also create a negative ion at this place.
   subroutine createIonPair(pos, eV, weight)

      real(dp), dimension(3), intent(IN) :: pos
      real(dp), intent(IN)               :: eV
      integer, intent(IN)                        :: weight

      call PM_createElectron(pos, eV, weight)
      call E_add_to_var(E_i_pion, pos, dble(weight))

   end subroutine createIonPair

   !> Sets the acceleration for the particles. Normally this happens automatically,
   ! in the routine that updates the velocity, but this is for initialization.
   subroutine updateParticleAccel()

      integer :: ll
      real(dp) :: Efield(3)

      do ll = 1, nParticles
         Efield = E_get_field(pList(ll)%x)
         pList(ll)%accel = elecChargeOverMass * Efield
      end do

   end subroutine updateParticleAccel

   !> Convert the particle list to a density, by mapping all the particles
   !! to densities on a grid.
   subroutine PM_particlesToDensity()
      integer :: n
      call E_set_vars((/E_i_elec/), (/0.0_dp/))
      do n = 1, nParticles
         call E_add_to_var(E_i_elec, pList(n)%x, dble(pList(n)%weight))
      end do
   end subroutine PM_particlesToDensity

   !> Compare the stores acceleration of particles with the acceleration obtained
   !! from the current electric field, and estimate the average error, which is returned
   !! through the argument 'interpErr'.
   subroutine checkinterpError(interpErr, dt)
      real(dp), intent(INOUT)  :: interpErr
      real(dp), intent(IN)     :: dt

      integer :: i, nTries, n
      real(dp), dimension(3)   :: Efield, oldEfield, newPos
      real(dp)                 :: diff, total

      if (nParticles == 0) then
         print *, "checkinterpError: No particles available"
         interpErr = 1.0D0
         return
      end if

      diff = 0.0D0
      total = 0.0D0

      i = 1
      nTries = 1
      do while (i <= nSamplesEfield .and. nTries < 10 * nSamplesEfield)
         nTries = nTries + 1

         ! Randomly select a particle
         n = int(kiss_rand() * nParticles) + 1

         ! Advance the position over dt
         newPos = pList(n)%x + pList(n)%v * dt + 0.5D0 * pList(n)%accel * dt**2

         if (isOutOfGas(newPos)) then
            cycle
         else
            i = i + 1

            Efield = E_get_field(newPos)
            oldEfield = pList(n)%accel * electronMass / electronCharge

            total = total + norm2(oldEfield)
            diff  = diff + norm2(Efield - oldEfield)
         end if

      end do
      interpErr = diff / (total + smallNumber)

   end subroutine checkinterpError

   subroutine checkinterpErrorMPI(interpErr, tau, nParticlesTotal)
      use mpi
      integer, intent(IN) :: nParticlesTotal
      real(dp), intent(INOUT)  :: interpErr
      real(dp), intent(IN)     :: tau
      real(dp) :: temp
      integer :: ierr

      call checkinterpError(temp, tau)
      temp = temp * nParticles / dble(nParticlesTotal)
      call MPI_ALLREDUCE(temp, interpErr, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

   end subroutine checkinterpErrorMPI

   !> Store samples of the electric field at particle positions, to estimate the error
   !! later on. The number of samples is controlled by the nSamplesEfield parameter.
   subroutine storeSampleOfEfield()
      integer :: i, n

      if (nParticles == 0) then
         return
      end if

      do i = 1, nSamplesEfield
         ! Randomly select a particle
         do
            n = int(kiss_rand() * nParticles) + 1
            if (.not. isOutOfGas(pList(n)%x)) exit
         end do
         samplePos(:,i) = pList(n)%x
         sampleEfield(:,i) = E_get_field(samplePos(:,i))
      end do

   end subroutine storeSampleOfEfield

   !> Compare old values of the Efield to the updated ones, to estimate the error
   !! in the electric field. The error is returned through the argument 'PoissonErr'.
   subroutine checkPoissonError(PoissonErr)

      real(dp), intent(INOUT) :: PoissonErr
      integer :: i
      real(dp), dimension(3) :: Efield
      real(dp) :: diff, total, maxDiff, maxRelDiff

      if (nParticles == 0) then
         PoissonErr = 0.0D0
         return
      end if

      maxDiff = 0.0D0
      maxRelDiff = 0.0D0
      !       nSamples = MAX(nSamplesEfield, nParticles / 100)

      !       DO i = 1, nSamplesEfield
      !          CALL PS_getEfield(samplePos(:,i), Efield)
      !          diff = norm2( Efield-sampleEfield(:,i) )
      !          IF (diff > maxDiff) THEN
      !             maxDiff = diff
      !             maxRelDiff = maxDiff / MAX(norm2(Efield), norm2(sampleEfield(:,i)))
      !          END IF
      !       END DO
      !
      !       PoissonErr = maxRelDiff
      total = 0.0D0
      diff = 0.0D0

      do i = 1, nSamplesEfield
         Efield = E_get_field(samplePos(:,i))
         total = total + norm2( sampleEfield(:,i) )
         diff  = diff + norm2( Efield-sampleEfield(:,i) )
      end do

      PoissonErr = diff / (total + smallNumber)

   end subroutine checkPoissonError

   subroutine checkPoissonErrorMPI(PoissonErr, nParticlesTotal)
      use mpi
      integer, intent(IN) :: nParticlesTotal
      real(dp), intent(INOUT) :: PoissonErr
      real(dp) :: temp
      integer :: ierr

      call checkPoissonError(temp)
      temp = temp * nParticles / dble(nParticlesTotal)

      call MPI_ALLREDUCE(temp, PoissonErr, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

   end subroutine checkPoissonErrorMPI

   subroutine selectParticlesForMerging(nPartMerge)
      integer, intent(inout) :: nPartMerge

      integer :: ll, ix_end
      real(dp) :: n_elec, max_weight
      real(dp), allocatable :: weightRatio(:)

      ix_end      = nParticles
      nPartMerge  = nParticles

      do ll = 1, nParticles
         if (.not. canParticleBeMerged(ll)) then
            ! Swap it with the first particle from the end of the list that is a candidate for merging

            do while (.not. canParticleBeMerged(ix_end) .and. ix_end > ll)
               ix_end = ix_end - 1
            end do

            ! If there is no merge-able particle available in the range [ll+1 : nParticles] we have
            ! ll >= ix_end and we exit the routine. Else we put the merge-able particle at index ll
            if (ll >= ix_end) then
               nPartMerge = ll - 1
               exit
            else
               call swapPart(pList(ix_end), pList(ll))
               ix_end               = ix_end - 1
            end if
         end if
      end do

      if (nPartMerge > 0) then
         ! Now sort the particles by their weight/max_weight ratio (low to high)
         allocate( weightRatio(nPartMerge) )

         do ll = 1, nPartMerge
            n_elec          = E_get_var_x_dr(E_i_elec, pList(ll)%x)
            max_weight      = n_elec / dble(PM_min_part_per_cell)
            weightRatio(ll) = dble(pList(ll)%weight) / max_weight
         end do

         call heapsortParticleByDbles(pList(1:nPartMerge), weightRatio)
         deallocate( weightRatio )
      end if

   end subroutine selectParticlesForMerging

   logical function canParticleBeMerged(ll)
      integer, intent(in) :: ll
      real(dp)            :: n_elec
      integer             :: max_weight

      n_elec = E_get_var_x_dr(E_i_elec, pList(ll)%x)
      max_weight = max(1, min(PM_max_weight, floor(n_elec / PM_min_part_per_cell)))

      canParticleBeMerged = (pList(ll)%weight <= 0.6667d0 * max_weight)
   end function canParticleBeMerged

   logical function shouldParticleBeSplit(ll)
      integer, intent(in)  :: ll

      real(dp)     :: n_elec
      integer              :: max_weight

      n_elec                   = E_get_var_x_dr(E_i_elec, pList(ll)%x)
      max_weight               = max(1, min(PM_max_weight, floor(n_elec / PM_min_part_per_cell)))

      shouldParticleBeSplit = (pList(ll)%weight > 1.5d0 * max_weight)
   end function shouldParticleBeSplit

   subroutine PM_mergeAndSplit(myrank, doMerge, doSplit)

      use kdtree2_module
      integer, intent(in)                 :: myrank
      logical, optional, intent(in)       :: doMerge, doSplit

      integer, parameter                  :: partDim = 4
      integer, parameter                  :: nNeighbours = 1
      type(kdtree2), pointer              :: kdTree
      type(kdtree2_result)                :: kdResults(nNeighbours)
      real(dp), allocatable       :: partCoords(:, :)
      logical, allocatable                :: partIsMerged(:)
      logical                             :: performMerge, performSplit

      real(dp), parameter         :: typicalLengthScale = 1.0d-6
      real(dp), parameter         :: typicalVelocityScale = 1.0d6
      real(dp), parameter         :: velocityToLength = typicalLengthScale / typicalVelocityScale

      integer                             :: pIx, nIx, nSuccess, nPartMerge
      integer                             :: tooFar, alreadyDone
      real(dp)                    :: w1, w2, temp, newPos(3), newVel(3), tempVec(3)
      real(dp)                    :: distance, maxDistance, delta_pos(3), separate_fac

      performMerge   = .true.
      performSplit   = .true.
      if (present(doMerge)) performMerge = doMerge
      if (present(doSplit)) performSplit = doSplit

      ! if (myrank == 0) print *, "MERGING", performMerge, "SPLITTING", performSplit

      call selectParticlesForMerging(nPartMerge)
      ! print *, myrank, "Number of candidates for merging:", nPartMerge, "/", nParticles

      if (performMerge .and. nPartMerge > 10) then ! Limit of the k-d tree code, nPartMerge should not be too small
         allocate( partCoords(partDim, nPartMerge) )
         allocate( partIsMerged(nPartMerge) )

         partIsMerged(:) = .false.

         do pIx = 1, nPartMerge
            partCoords(1:3, pIx) = pList(pIx)%x
            partCoords(4, pIx)   = velocityToLength * twoNorm(pList(pIx)%v)
         end do

         ! print *, myrank, nPartMerge, "Getting k-d tree"
         ! print *, any(partCoords /= partCoords), size(partCoords), maxval(partCoords), minval(partCoords)
         kdTree      => kdtree2_create(partCoords)
         ! print *, myrank, "Done"
         nSuccess    = 0
         alreadyDone = 0
         tooFar      = 0

         do pIx = 1, nPartMerge
            if (partIsMerged(pIx)) cycle

            call kdtree2_n_nearest_around_point(kdTree, idxin=pIx, nn=nNeighbours, correltime=1, results=kdResults)
            nIx = kdResults(1)%idx

            if (partIsMerged(nIx)) then
               alreadyDone = alreadyDone + 1
               cycle
            end if

            distance = twoNorm(partCoords(:, pIx) - partCoords(:, nIx))
            maxDistance = PM_maxRelDistForMerge * maxval(E_get_dr(pList(pIx)%x))

            if (distance > maxDistance) then
               tooFar = tooFar + 1
               cycle
            end if

            ! Merge particles pIx and nIx into a temporary particle
            w1                   = dble(pList(pIx)%weight)
            w2                   = dble(pList(nIx)%weight)

            newPos               = (w1 * pList(pIx)%x + w2 * pList(nIx)%x) / (w1 + w2)
            temp                 = sqrt( (w1*sum(pList(pIx)%v**2) + w2*sum(pList(nIx)%v**2)) / (w1 + w2) )
            if (kiss_rand() < w1 / (w1+w2)) then
               tempVec              = pList(pIx)%v
            else
               tempVec              = pList(nIx)%v
            end if
            newVel               = temp * tempVec / (sqrt(sum(tempVec**2)) + smallVelocity)

            pList(pIx)%weight    = pList(pIx)%weight + pList(nIx)%weight
            pList(pIx)%x         = newPos
            pList(pIx)%v         = newVel

            pList(nIx)%live      = .false.
            partIsMerged(pIx)    = .true.
            partIsMerged(nIx)    = .true.
            nSuccess             = nSuccess + 1

         end do

         ! print *, myrank, "destroy k-d tree"
         call kdtree2_destroy(kdTree)
         deallocate( partCoords )
         deallocate( partIsMerged )
         ! print *, myrank, "nSuccess", nSuccess, "tooFar", tooFar, "alreadyDone", alreadyDone
      end if

      ! Possibly split the particles that were not candidates to be merged
      separate_fac = (1.0_dp / PM_min_part_per_cell)**(1.0_dp/3)

      if (performSplit) then
         ! print *, "Splitting", nPartMerge, nParticles
         do pIx = nPartMerge + 1, nParticles
            if (shouldParticleBeSplit(pIx)) then
               call checkNumParticles(nParticles+1)
               nParticles               = nParticles + 1
               pList(nParticles)        = pList(pIx)
               pList(pIx)%weight        = (pList(pIx)%weight + 1) / 2
               pList(nParticles)%weight = pList(nParticles)%weight / 2

               ! Separate particles
               delta_pos = E_get_dr(pList(pIx)%x) * kiss_rand() * separate_fac
               pList(nParticles)%x = pList(nParticles)%x + delta_pos
               pList(pIx)%x = pList(pIx)%x - delta_pos
            end if
         end do
      end if

      ! print *, "Removing dead particles"
      call removeDeadParticles()
   end subroutine PM_mergeAndSplit

   subroutine checkNumParticles(nPart)
      integer, intent(in)         :: nPart
      integer                     :: new_size
      type(Electron), allocatable :: copy_particles(:)

      if (nPart > PM_maxNumberOfParticles) then
         print *, "Too many particles:", nPart

         allocate(copy_particles(PM_maxNumberOfParticles))
         copy_particles = pList
         deallocate(pList)

         new_size = max(ceiling(1.5d0 * nPart), nPart + 100*1000)
         allocate(pList(new_size))
         print *, "allocated more space:", new_size

         pList(1:PM_maxNumberOfParticles) = copy_particles
         PM_maxNumberOfParticles = new_size
      end if
   end subroutine checkNumParticles

   !> Sort list of particle indices based on particle velocity
   recursive subroutine QsortParticlesVelocity(list)

      integer, intent(INOUT) :: list(:)
      integer :: split

      if(size(list) > 1) then
         call partitionParticlesVelocity(list, split)
         call QsortParticlesVelocity( list(:split-1)  )
         call QsortParticlesVelocity( list(split:)    )
      end if

   end subroutine QsortParticlesVelocity

   !> Divide list of particle indices in two parts based on particle velocity
   subroutine partitionParticlesVelocity(list, marker)

      integer, intent(INOUT) :: list(:)
      integer, intent(OUT) :: marker

      ! Temporary variables
      integer :: i, left, right, iPivot, storeIndex
      real(dp) :: pivotValue

      left        = 1
      storeIndex  = 1
      right       = size(list)
      iPivot      = (left + right) / 2

      pivotValue = sum(pList(list(iPivot))%v**2)

      ! Swap pivot to the right
      call swapInt(list(iPivot), list(right))

      do i = left, right - 1
         if (sum(pList(list(i))%v**2) < pivotValue) then
            call swapInt(list(i), list(storeIndex))
            storeIndex = storeIndex + 1
         end if
      end do

      call swapInt(list(right), list(storeIndex))
      marker = min(storeIndex + 1, right)

   end subroutine partitionParticlesVelocity

   subroutine swapPart(a, b)
      type(Electron), intent(inout) :: a, b
      type(Electron)                :: c
      c = a
      a = b
      b = c
   end subroutine swapPart

   !> Return the number of particles in the simulation
   integer function get_num_part()
      get_num_part = nParticles
   end function get_num_part

   !> MPI: Share the particles evenly over the tasks, with no domain decomposition
   subroutine shareParticles(myrank, ntasks)
      include 'mpif.h'
      integer, intent(IN) :: myrank, ntasks

      integer :: iHasMore, iHasLess, nPartSend, nRecvs, nSends
      integer :: ierr, meanNP, iMin, tag
      integer :: mpi_status(MPI_STATUS_SIZE)
      integer, dimension(0:ntasks-1) :: nPartPerTask

      ! print *, myrank, "Before sharing: nParticles", nParticles
      ! Get the number of particles each task has
      call MPI_ALLGATHER(nParticles, 1, MPI_integer, nPartPerTask, 1, MPI_integer, MPI_COMM_WORLD, ierr)
      meanNP = (sum(nPartPerTask) + ntasks - 1) / ntasks

      if (maxval(abs(nPartPerTask - meanNP)) * 100 < meanNP) then
         ! Small difference, no need to correct for it
         return
      end if

      nSends = 0
      nRecvs = 0
      iHasLess = 0
      iHasMore = 0

      do while (iHasMore < ntasks)

         if (nPartPerTask(iHasMore) > meanNP) then

            ! Find first task with less than mean particles
            do while (nPartPerTask(iHasLess) >= meanNP)
               iHasLess = iHasLess + 1
            end do

            ! Determine amount of particles to send
            nPartSend = min(nPartPerTask(iHasMore) - meanNP, meanNP - nPartPerTask(iHasLess))

            if (myrank == iHasMore) then
               ! Send particles
               iMin     = nPartPerTask(iHasMore) - nPartSend + 1
               tag      = iHasMore
               nSends   = nSends + 1
               call MPI_SEND( pList(iMin), nPartSend, partTypeMPI, iHasLess, tag, &
                    & MPI_COMM_WORLD, ierr)
               ! print *, myrank, "sends", iMin, "--", iMin + nPartSend - 1, "to", iHasLess

            else if (myrank == iHasLess) then
               ! Receive particles
               iMin     = nPartPerTask(iHasLess) + 1
               call checkNumParticles(iMin + nPartSend - 1)
               tag      = iHasMore
               nRecvs   = nRecvs + 1
               call MPI_RECV( pList(iMin), nPartSend, partTypeMPI, iHasMore, tag, &
                    MPI_COMM_WORLD, mpi_status, ierr)
               ! print *, myrank, "recvs", iMin, "--", iMin + nPartSend - 1, "from", iHasMore
            end if

            ! Update the number of particles for each task
            nPartPerTask(iHasLess) = nPartPerTask(iHasLess) + nPartSend
            nPartPerTask(iHasMore) = nPartPerTask(iHasMore) - nPartSend
         end if

         if (nPartPerTask(iHasMore) <= meanNP) iHasMore = iHasMore + 1

      end do

      ! if (nRecvs > 0) call MPI_WAITALL(nRecvs, recvReqs, MPI_STATUSES_IGNORE, ierr)
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      nParticles = nPartPerTask(myrank)
      ! print *, myrank, "After sharing: nParticles", nParticles

   end subroutine shareParticles


   !> Divide the particles over the tasks based on their cell index
   subroutine PM_divideParticles(myrank, ntasks)
      include 'mpif.h'
      integer, intent(in) :: myrank, ntasks

      integer :: ll, ix, rank, ierr, iMin, nPartBehind, nPartInInterval, nPartSend
      integer :: tag, partCnt
      integer :: nSends, nRecvs, sender, recver, meanNP
      integer :: mpi_status(MPI_STATUS_SIZE)
      integer :: nPartPerTask(0:ntasks-1), nPartPerIntervalPerTask(0:ntasks-1, 0:ntasks-1)
      integer :: splitIx(-1:ntasks)
      integer, parameter :: n_cells = 10000
      integer, allocatable :: nPartPerCellIx(:)
      real(dp), allocatable :: cellIxs(:)
      real(dp) :: x_min, x_max, dx, inv_dx

      ! print *, "Before divide ", myrank, " has ", nParticles, " particles"

      ! Get the number of particles each task has
      call MPI_ALLGATHER(nParticles, 1, MPI_integer, nPartPerTask, 1, MPI_integer, MPI_COMM_WORLD, ierr)

      x_min = minval(pList(1:nParticles)%x(1))
      x_max = maxval(pList(1:nParticles)%x(1))
      call MPI_ALLREDUCE(MPI_IN_PLACE, x_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, x_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

      ! Add small extra space to avoid indexing problems at boundaries
      x_min = x_min * (1-1e-6_dp)
      x_max = x_max * (1+1e-6_dp)
      dx = (x_max - x_min) / n_cells
      inv_dx = 1 / dx

      meanNP      = (sum(nPartPerTask) + ntasks - 1) / ntasks
      allocate( nPartPerCellIx(n_cells) )
      allocate( cellIxs(nParticles) )

      nPartPerCellIx = 0

      ! Set the cell indexes for the particles
      do ll = 1, nParticles
         ix                   = int((pList(ll)%x(1) - x_min) * inv_dx) + 1
         cellIxs(ll)          = ix
         nPartPerCellIx(ix)   = nPartPerCellIx(ix) + 1
      end do

      call MPI_ALLREDUCE(MPI_IN_PLACE, nPartPerCellIx, n_cells, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

      call heapsortParticleByDbles(pList(1:nParticles), cellIxs)

      ! Find the splitIx(:) s.t. task n will hold particles with ix > splitIx(n-1) and ix <= splitIx(n)
      splitIx(-1)                   = 0
      splitIx(ntasks)               = n_cells
      rank                          = 0
      nPartBehind                   = 0
      nPartPerIntervalPerTask(:,:)  = 0
      nPartInInterval               = 0

      do ix = 1, n_cells
         nPartInInterval = nPartInInterval + nPartPerCellIx(ix)

         if (nPartInInterval >= meanNP .or. ix == n_cells) then
            splitIx(rank)     = ix
            nPartInInterval   = 0
            partCnt           = 0

            do ll = nPartBehind+1, nParticles
               !                print *, ll, cellIxs(ll)
               if (cellIxs(ll) <= ix) then
                  partCnt = partCnt + 1
               else
                  exit
               end if
            end do

            nPartPerIntervalPerTask(rank, myrank) = partCnt
            nPartBehind       = nPartBehind + partCnt
            rank              = rank + 1

            if (rank == ntasks) exit
         end if

      end do

      call MPI_ALLREDUCE(MPI_IN_PLACE, nPartPerIntervalPerTask, ntasks*ntasks, &
           MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

      ! Send particles to task that holds a region
      nSends = 0
      nRecvs = 0

      do recver = 0, ntasks - 1

         ! Receive the particles in tasks' region from all other tasks
         if (myrank == recver) then

            do sender = 0, ntasks - 1
               if (sender == myrank) cycle ! Don't have to send to ourselves

               nPartSend = nPartPerIntervalPerTask(recver, sender)

               if (nPartSend > 0) then
                  nRecvs = nRecvs + 1
                  iMin   = nParticles + 1
                  nParticles = nParticles + nPartSend
                  call checkNumParticles(nParticles)
                  tag    = sender

                  !                   print *, myrank, ": receives", nPartSend, "from", sender, "iMin", iMin
                  call MPI_RECV( pList(iMin), nPartSend, partTypeMPI, sender, tag, &
                       & MPI_COMM_WORLD, mpi_status, ierr)
               end if
            end do


         else ! Send particles to recver

            nPartSend = nPartPerIntervalPerTask(recver, myrank)

            if (nPartSend > 0) then
               nSends      = nSends + 1
               tag         = myrank

               ! Find index of first particle to send
               iMin        = sum( nPartPerIntervalPerTask(0:recver-1, myrank) ) + 1

               !                print *, myrank, ": sends", nPartSend, "to", recver, "iMin", iMin
               call MPI_SEND( pList(iMin), nPartSend, partTypeMPI, recver, tag, &
                    & MPI_COMM_WORLD, ierr)

               ! Mark the particles that we have send away as inactive
               do ll = iMin, iMin + nPartSend - 1
                  call killParticle(ll)
               end do

            end if
         end if

      end do

      ! if (nRecvs > 0) call MPI_WAITALL(nRecvs, recvReqs, MPI_STATUSES_IGNORE, ierr)
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      !       print *, "After divide ", myrank, " has ", nParticles, "nParticles (some dead)"
      call removeDeadParticles()
      ! print *, "After divide ", myrank, " has ", nParticles, " particles"

      deallocate( nPartPerCellIx )
      deallocate( cellIxs )

   end subroutine PM_divideParticles

   subroutine createPartTypeMPI(partTypeMPI)
      use mpi
      integer, intent(INOUT) :: partTypeMPI
      integer, parameter :: nf = 3
      integer :: blocklen(nf), types(nf), ierr, type_size, my_type_size
      integer(KIND=MPI_ADDRESS_KIND) :: displ(nf), lb, extent
      type(electron) :: electron_example

      blocklen = (/11, 1, 1/) ! 11 doubles, 1 integer, 1 logical

      displ(1) = 0
      call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION, lb, extent, ierr)
      displ(2) = displ(1) + blocklen(1) *  extent
      call MPI_TYPE_GET_EXTENT(MPI_integer, lb, extent, ierr)
      displ(3) = displ(2) + blocklen(2) * extent

      types = (/MPI_DOUBLE_PRECISION, MPI_integer, MPI_logical/)

      call MPI_TYPE_CREATE_STRUCT(nf, blocklen, displ, types, partTypeMPI, ierr)
      call MPI_TYPE_COMMIT(partTypeMPI, ierr)

      type_size = storage_size(electron_example)
      my_type_size = blocklen(1) * storage_size(0.0d0) + blocklen(2) * storage_size(1) + blocklen(3) * storage_size(.true.)

      if (type_size /= my_type_size .or. ierr /= 0) then
         print *, "createPartTypeMPI error!", type_size, my_type_size, ierr
         stop
      end if
   end subroutine createPartTypeMPI

   !> Return the number of electrons
   real(dp) function get_num_elec()
      ! Return the number of electrons simulated, different from the number of particles
      ! if rescaling has occurred.
      integer :: n
      integer(KIND=iKind18) :: sumWeight

      sumWeight = 0
      do n = 1, nParticles
         sumWeight = sumWeight + pList(n)%weight
      end do
      get_num_elec = dble(sumWeight)
   end function get_num_elec

   !> MPI: return the total number of electrons in the simulation
   real(dp) function PM_get_num_elec_mpi()
      use mpi
      integer :: ierr
      real(dp) :: partSum

      partSum = get_num_elec()
      call MPI_ALLREDUCE(partSum, PM_get_num_elec_mpi, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
           MPI_COMM_WORLD, ierr)

   end function PM_get_num_elec_mpi

   !> MPI: return the total number of particles in the simulation
   integer function PM_get_num_part_mpi()
      use mpi
      integer :: ierr
      integer :: partSum

      partSum = get_num_part()
      call MPI_ALLREDUCE(partSum, PM_get_num_part_mpi, 1, MPI_integer, MPI_SUM, &
           MPI_COMM_WORLD, ierr)

   end function PM_get_num_part_mpi

   !> Return the total energy in eV of the particles
   real(dp) function getTotalEnergy()
      integer :: n
      real(dp) :: sumE

      sumE = 0
      do n = 1, nParticles
         sumE = sumE + velocityToEnergy(pList(n)%v)*pList(n)%weight
      end do
      getTotalEnergy = sumE
   end function getTotalEnergy

   !> Return the mean energy in eV of the particles
   real(dp) function getMeanEnergy()
      getMeanEnergy = getTotalEnergy() / get_num_elec()
   end function getMeanEnergy

   !> Return the mean weight of the particles
   real(dp) function getMeanWeight()
      getMeanWeight = get_num_elec() / dble(nParticles)
   end function getMeanWeight

   subroutine heapsortParticleByDbles(partList, dbleList)
      type(Electron), intent(inout)    :: partList(0:)
      real(dp), intent(inout)  :: dbleList(0:)

      integer                          :: start, listSize, bottom

      listSize = size(partList)

      if (size(dbleList) /= listSize) then
         print *, "heapsortParticleByDbles error: arguments have different size"
         stop
      end if

      do start = (listSize - 2) / 2, 0, -1
         call siftdown(partList, dbleList, start, listSize);
      end do

      do bottom = listSize - 1, 1, -1
         call swapDble(dbleList(0), dbleList(bottom))
         call swapPart(partList(0), partList(bottom))
         call siftdown(partList, dbleList, 0, bottom)
      end do

   end subroutine heapsortParticleByDbles

   subroutine siftdown(partList, dbleList, start, bottom)
      type(Electron), intent(inout)    :: partList(0:)
      real(dp), intent(inout)  :: dbleList(0:)
      integer, intent(in)              :: start, bottom

      integer                          :: child, root

      root = start

      do while(root*2 + 1 < bottom)

         child = root * 2 + 1

         if (child + 1 < bottom) then
            if (dbleList(child) < dbleList(child+1)) then
               child = child + 1
            end if
         end if

         if (dbleList(root) < dbleList(child)) then
            call swapDble(dbleList(child), dbleList(root))
            call swapPart(partList(child), partList(root))
            root = child
         else
            return
         end if

      end do

   end subroutine siftdown

   !> Create a lookup table with cross sections for a number of energies
   subroutine createColRateTable(maxVel, nIntervals)
      real(dp), intent(IN)  :: maxVel
      integer, intent(IN)           :: nIntervals

      integer           :: i
      real(dp)  :: velocity, energy_eV, factor

      allocate( colRateTable(nProcesses, nIntervals) )
      allocate( colRateTableVelocity(nIntervals) )
      allocate( colRateTableSum(nIntervals) )

      do i = 1, nIntervals

         factor      = dble(i-1)/dble(nIntervals-1)
         velocity    = maxVel * (factor**2)
         energy_eV   = v2en(velocity) / electronVolt

         colRateTableVelocity(i) = velocity

         call getCrossSections(energy_eV, colRateList)

         ! To get a collision rate we multiply by the electron velocity.
         colRateTable(:, i) = colRateList(:) * velocity

         ! Keep track of the sum of the cross sections to later speed up some routines
         colRateTableSum(i) = sum( colRateTable(:, i) )
      end do

      maxCollRate = 1.01d0 * maxval(colRateTableSum)
      invMaxCollRate = 1.0d0 / maxCollRate
      print *, 'Maximum collision frequency = ', maxCollRate

   end subroutine createColRateTable

   !> Determine the index in a table for a given velocity. All tables are scaled
   !! in the same way so that we can quickly do this
   integer function getTableIndex(velocity)
      real(dp), intent(IN)  :: velocity
      real(dp)              :: tempDouble

      tempDouble = velocity * invMaxVelocity
      tempDouble = sqrt(tempDouble)
      getTableIndex = int(tempDouble * (tableSize-1)) + 1
      getTableIndex = min(getTableIndex, tableSize-1)

   end function getTableIndex

   !> Look up the collision rates for a given velocity
   subroutine lookupColRate(velocity, rateList)
      real(dp), intent(IN) :: velocity
      real(dp), dimension(:), intent(INOUT) :: rateList

      integer :: tablePos
      real(dp) :: tempDouble, lowerVelocity, upperVelocity

      ! Find the interval where we find the velocity, tablePos holds the lower index
      tablePos = getTableIndex(velocity)

      ! Find the velocity at the interval points
      lowerVelocity = colRateTableVelocity(tablePos)
      upperVelocity = colRateTableVelocity(tablePos + 1)

      ! Set the weight for the lower interval cross section data
      tempDouble = (upperVelocity - velocity) / (upperVelocity - lowerVelocity)

      ! Set the cross sections
      rateList = tempDouble * colRateTable(:, tablePos) &
           & + (1.0D0 - tempDouble) * colRateTable(:, tablePos + 1)

   end subroutine lookupColRate

   !> For a given energy, look up the cross sections of all the collision
   !! processes.
   !! This routine is only used at the start, later on we use a
   !! lookup table for much better performance.
   subroutine getCrossSections(energy, CSList)
      real(dp), intent(IN)                    :: energy
      real(dp), dimension(:), intent(INOUT)   :: CSList

      integer          :: n

      do n = 1, nProcesses

         call linearInterpolateList(   crossSecTable(n)%en_cs(1, :), &
              &  crossSecTable(n)%en_cs(2, :), &
              &  energy, CSList(n) )
      end do

   end subroutine getCrossSections

   !> Get the un-normalized electron energy distribution function at the argument
   !! 'EEDF', up to an energy of 'eVmax' electron Volt with 'nSteps' steps.
   subroutine getEEDF(EEDF, eVmax, nSteps)
      integer, intent(IN) :: nSteps
      real(dp), intent(IN) :: eVmax
      real(dp), dimension(2, nSteps), intent(INOUT) :: EEDF

      integer :: i, n, prevIndex
      integer, dimension(nParticles) :: ixList
      real(dp) :: maxEn, deltaEn

      EEDF(:,:) = 0.0D0

      if (nParticles < 2) then
         print *, "writeElectronEDF warning: no particles available"
         return
      end if

      ! Create a list of particle indices that is sorted by energy (from low to high)
      do n = 1, nParticles
         ixList(n) = n
      end do

      call QsortParticlesVelocity( ixList(1:nParticles) )

      ! Set the maximum energy and create a histogram of the EEDF
      maxEn = eVmax * electronVolt
      deltaEn = maxEn / (nSteps-1)

      prevIndex = 1

      do n = 1, nSteps
         EEDF(1,n) = n * deltaEn / electronVolt
      end do

      n = 1

      do while(n <= nSteps)
         i = prevIndex

         ! Find the number of electrons with energy less than n * deltaEn
         do while ( velocityToEnergy(pList(ixList(i))%v) < n * deltaEn)
            EEDF(2,n) = EEDF(2,n) + dble(pList(ixList(i))%weight)

            i = i + 1

            if (i == nParticles + 1) then
               ! Make sure to exit the outer loop by increasing n
               n = nSteps
               exit
            end if

         end do

         prevIndex = i
         n = n + 1
      end do

   end subroutine getEEDF

   !> Store the electron energy distribution function (EEDF) from all tasks
   subroutine setEEDF_MPI(eVmax, nSteps, myrank, rootMPI)
      use mpi
      real(dp), intent(IN) :: eVmax
      integer, intent(IN) :: nSteps, myrank, rootMPI

      integer :: ierr
      real(dp), dimension(:,:), allocatable :: EEDF

      if (allocated(electronEDF)) then
         if (size(electronEDF) /= 2 * nSteps) then
            deallocate(electronEDF)
            allocate( electronEDF(2, nSteps) )
         end if
      else
         allocate( electronEDF(2, nSteps) )
      end if

      allocate( EEDF(2, nSteps) )
      call getEEDF(EEDF, eVmax, nSteps)

      ! Call mpi_reduce on boths columns, then correct one later (seems easier)
      call MPI_REDUCE(EEDF, electronEDF, 2*nSteps, MPI_DOUBLE_PRECISION, MPI_SUM, &
           & rootMPI, MPI_COMM_WORLD, ierr)

      if (myrank == rootMPI) then
         ! Correct first column
         electronEDF(1,:) = EEDF(1,:)
         ! Normalize the EEDF
         electronEDF(2,:) = electronEDF(2,:) / sum(electronEDF(2,:))
      end if

      deallocate( EEDF )
   end subroutine setEEDF_MPI

   !> Store the electron energy distribution function (EEDF)
   subroutine setEEDF(eVmax, nSteps)
      real(dp), intent(IN) :: eVmax
      integer, intent(IN) :: nSteps

      if (allocated(electronEDF)) then
         if (size(electronEDF) /= 2 * nSteps) then
            deallocate(electronEDF)
            allocate( electronEDF(2, nSteps) )
         end if
      else
         allocate( electronEDF(2, nSteps) )
      end if

      call getEEDF(electronEDF, eVmax, nSteps)

   end subroutine setEEDF

   !> Get the number of collisions for process with index pIx, returns a double
   !! to not need non-standard integers in the calling routine.
   real(dp) function getNCollisionsOfProcess(pIx)
      integer, intent(IN) :: pIx
      getNCollisionsOfProcess = dble( nCollisionSum(pIx) )
   end function getNCollisionsOfProcess

   !> Get the number of collisions for processses of type colType, return a double
   !! to not need non-standard integers in the calling routine.
   real(dp) function getNCollisionsOfType(colType)
      integer, intent(IN) :: colType
      integer :: pIx
      real(dp) :: colsum

      colsum = 0.0D0
      do pIx = 1, nProcesses
         if (CrossSecTable(pIx)%colType == colType) then
            colsum = colsum + dble( nCollisionSum(pIx) )
         end if
      end do

      getNCollisionsOfType = colsum
   end function getNCollisionsOfType

   !> Reset the number of collisions
   subroutine resetNCollisions()
      nCollisionSum(:) = 0
   end subroutine resetNCollisions

   !> Remove all the particles
   subroutine removeAllParticles()
      nParticles = 0
   end subroutine removeAllParticles

   !> Perform attachment of electrons to stabilize the simulation in regions
   !! with a very high positive ion density.
   subroutine PM_convertElecForStability(dt, mu_max)
      use module_constants
      real(dp), intent(in) :: dt, mu_max

      integer             :: n, nAttach
      real(dp)            :: localDens, recombProb, att_threshold_dens, max_pion_dens(1)
      real(dp), parameter :: accelToEfield = 1.0d0 / elecChargeOverMass

      att_threshold_dens = epsilon0 / (mu_max * dt * elemCharge)
      call E_get_max_of_vars((/E_i_pion/), max_pion_dens)

      if (max_pion_dens(1) > att_threshold_dens) then
         print *, "Artificially attaching electrons for stability"
         print *, "Threshold density", att_threshold_dens
         nAttach = 0

         do n = 1, nParticles
            localDens   = E_get_var_x_dr(E_i_pion, pList(n)%x) ! Assume elec/ion dens are almost equal
            recombProb  = (localDens - att_threshold_dens) / PM_maxPosIonDens
            if (kiss_rand() < recombProb) then
               call E_add_to_var(E_i_nion, pList(n)%x, dble(pList(n)%weight))
               call killParticle(n)
               nAttach = nAttach + 1
            end if
         end do

         if (nAttach > 0) call removeDeadParticles()
      end if
   end subroutine PM_convertElecForStability

end module module_particle
