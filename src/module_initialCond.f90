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

!> Module to set up the initial particles
module module_initialCond

   use module_globals
   use generalUtilities
   use module_particle
   use module_kiss
   use module_config
   use module_constants
   use module_electrode
   use m_efield_amr

   implicit none
   private
   public :: setInitialConditions

contains

   !> Sets up the initial particles for a simulation
   subroutine setInitialConditions()
      character(LEN=40) :: initCond
      logical :: success
      integer :: ll, i, j, k, ix
      integer ::  initWeight, nBackgroundPairs
      double precision :: electronEnergy, radius, nIonPairs
      double precision :: temp
      double precision :: backgroundDensity, bgO2MinDensity
      double precision, dimension(3) :: pos, initSeedPos

      ! Laser line
      double precision :: orthVec1(3), orthVec2(3), lineLen, minFac, maxFac, tempVec(3), lineDirection(3), lineBase(3)


      nIonPairs         = CFG_varDble("init_nIonPairs")
      electronEnergy    = CFG_varDble("init_electronEnergy")
      initWeight        = CFG_varInt("init_weightFactor")
      backgroundDensity = CFG_varDble("init_backgroundDensity")
      bgO2MinDensity    = CFG_varDble("init_O2MinBackgroundDens")

      ! Set the background negative ions, for detachment
      call E_set_vars((/E_i_O2m, E_i_pion/), (/bgO2MinDensity, bgO2MinDensity/))

      ! Set the background ionization level
      nBackgroundPairs = int(backgroundDensity * GL_gridVolume * product(GL_gridSize))
      print *, "Background density creates:", nBackgroundPairs
      do ll = 1, nBackgroundPairs
         pos = (/kiss_rand(), kiss_rand(), kiss_rand()/)
         pos = pos * GL_gridLength
         if (.not. isOutOfGas(pos)) call createIonPair(pos, electronEnergy, 1)
      end do

      print *, "Setting up the initial condition"

      ! Set the initial conditions
      call CFG_getVar("init_condType", initCond)

      select case (initCond)
      case ('seed')
         ! Electron/ion pairs start at fixed position
         call CFG_getVar("init_relSeedPos", initSeedPos)
         do ll = 1, int (nIonPairs / initWeight)
            pos = initSeedPos * GL_gridLength
            call createIonPair(pos, electronEnergy, initWeight)
         end do

      case ('double_seed')
         ! Electron/ion pairs start at fixed position
         call CFG_getVar("init_relSeedPos", initSeedPos)
         radius = CFG_varDble("init_seedPosRadius")

         do ll = 1, int(nIonPairs / (2 * initWeight))
            pos = initSeedPos * GL_gridLength
            pos(1) = pos(1) - radius
            call createIonPair(pos, electronEnergy, initWeight)
         end do

         do ll = 1, int(nIonPairs / (2 * initWeight))
            pos = initSeedPos * GL_gridLength
            pos(1) = pos(1) + radius
            call createIonPair(pos, electronEnergy, initWeight)
         end do

      case ('ionCloud')

         call CFG_getVar("init_relSeedPos", initSeedPos)
         initSeedPos = initSeedPos * GL_gridLength
         radius = CFG_varDble("init_seedPosRadius")

         temp = nIonPairs / (((2.0D0*pi)**1.5D0) * radius**3) * product(GL_gridDelta)

         do k = 1, GL_gridsize(3)
            do j = 1, GL_gridsize(2)
               do i = 1, GL_gridsize(1)
                  pos = (/i-1, j-1, k-1/) * GL_gridDelta
                  call E_add_to_var(E_i_pion, pos, &
                       temp * exp( - sum((pos - initSeedPos)**2) / (2.0D0 * radius**2) ) )
               end do
            end do
         end do

         !             DO ll = 1, INT (nIonPairs *1.d-2 / initWeight)
         !                pos = (/kiss_rand(), kiss_rand(), kiss_rand()/)
         !                pos = initSeedPos + pos * radius
         !                CALL E_addElectron(pos, DBLE(initWeight))
         !             END DO

      case ('Gaussian')
         ! Electron/ion pairs are distributed as a 3D Gaussian distribution
         ! with mean initSeedPos and sigma init_seedPosRadius
         radius = CFG_varDble("init_seedPosRadius")
         call CFG_getVar("init_relSeedPos", initSeedPos)
         initSeedPos = initSeedPos * GL_gridLength

         do ll = 1, int (nIonPairs / initWeight)
            success = .false.
            do while (.not. success)
               pos = (/kiss_normal(), kiss_normal(), kiss_normal()/)
               pos = initSeedPos + pos * radius
               if (.not. isOutOfGas(pos)) then
                  call createIonPair(pos, electronEnergy, initWeight)
                  success = .true.
               end if
            end do
         end do

      case ('Gaussian2D')
         ! Electron/ion pairs are distributed as a 2D Gaussian distribution in x, y directions, while uniform distribution
         ! along z direction
         ! with mean initSeedPos and sigma init_seedPosRadius
         radius = CFG_varDble("init_seedPosRadius")
         call CFG_getVar("init_relSeedPos", initSeedPos)
         initSeedPos = initSeedPos * GL_gridLength

         do ll = 1, int (nIonPairs / initWeight)
            success = .false.
            do while (.not. success)
               ! X and y directions with gaussian distribution, z direction with uniform distribution with semispheric cap
               pos = (/kiss_normal(), kiss_normal(), (20.d0-2.d0) * (kiss_rand()-0.5d0) + kiss_normal() /)
               pos = initSeedPos + pos * radius
               if (.not. isOutOfGas(pos)) then
                  call createIonPair(pos, electronEnergy, initWeight)
                  success = .true.
               end if
            end do

         end do

      case ('laserLine')
         radius = CFG_varDble("init_seedPosRadius")
         call CFG_getVar("init_laserDirection", lineDirection)
         call CFG_getVar("init_laserLineOffset", lineBase)

         ! Find first orthogonal vector
         if (abs(lineDirection(2)) < epsilon(1.0d0)) then
            orthVec1 = (/0.0d0, 1.0d0, 0.0d0/)
         else
            orthVec1 = (/-lineDirection(3), 0.0d0, lineDirection(1)/)
         end if

         ! Find second orthogonal vector by taking cross product
         orthVec2 = (/lineDirection(2) * orthVec1(3) - lineDirection(3) * orthVec1(2), &
              & lineDirection(1) * orthVec1(3) - lineDirection(3) * orthVec1(1), &
              & lineDirection(1) * orthVec1(2) - lineDirection(2) * orthVec1(1) /)

         ! Normalize
         lineDirection = lineDirection / twoNorm(lineDirection)
         orthVec1 = orthVec1 / twoNorm(orthVec1)
         orthVec2 = orthVec2 / twoNorm(orthVec2)

         tempVec = lineDirection
         where (abs(tempVec) < epsilon(1.0d0)) tempVec = epsilon(1.0d0)

         minfac = -huge(1.0d0)
         maxfac = huge(1.0d0)

         ! Get closest intersection in the domain [0,0,0] to [1,1,1] of line_base + lambda * line_direction
         do ix = 1, 3
            temp = -lineBase(ix) / tempVec(ix)
            if (temp < 0 .and. temp > minfac) minfac = temp
            if (temp > 0 .and. temp < maxfac) maxfac = temp
            temp = (1 - lineBase(ix)) / tempVec(ix)
            if (temp < 0 .and. temp > minfac) minfac = temp
            if (temp > 0 .and. temp < maxfac) maxfac = temp
         end do

         lineLen = twoNorm((maxFac - minFac) * lineDirection * GL_gridLength)
         temp = CFG_varDble("init_lineDens")
         nIonPairs = temp * radius**2 * lineLen

         print *, "Laser nIonPairs", nIonPairs/initWeight, lineLen, minFac, maxFac
         print *, orthVec1, orthVec2

         do ll = 1, int(nIonPairs/initWeight)
            pos = (kiss_rand() * (maxFac-minFac) + minFac) * lineDirection + lineBase
            pos = pos * GL_gridLength
            pos = pos + (kiss_rand()-0.5d0) * radius * orthVec1
            pos = pos + (kiss_rand()-0.5d0) * radius * orthVec2
            !                print *, ll, pos
            if (.not. isOutOfGas(pos)) call createIonPair(pos, electronEnergy, initWeight)
         end do

      case ('Gaussian1D')
         ! Electron/ion pairs are distributed as a 1D Gaussian distribution
         ! with mean initSeedPos and sigma init_seedPosRadius
         radius = CFG_varDble("init_seedPosRadius")
         call CFG_getVar("init_relSeedPos", initSeedPos)
         initSeedPos = initSeedPos * GL_gridLength

         do ll = 1, int (nIonPairs / initWeight)
            success = .false.
            do while (.not. success)
               pos = (/0.0D0, 0.0D0, kiss_normal()/)
               pos = initSeedPos + pos * radius
               if (.not. isOutOfGas(pos)) then
                  call createIonPair(pos, electronEnergy, initWeight)
                  success = .true.
               end if
            end do

         end do


      case DEFAULT
         print *, "initial condition error: wrong condition name"
         stop
      end select

      print *, "Initial conditions created"

   end subroutine setInitialConditions

   !    > Anbang : Check the time for backgound electrons go cross half domain in efield direction
   !    subroutine timeCheckForAddBGE(time)
   !       double precision, intent(out) :: time
   !       double precision :: velocity, acceleration
   !       double precision,dimension(3) :: efield
   !       double precision :: electronEnergy
   !
   !       electronEnergy = CFG_varDble("init_electronEnergy")
   !       efield = CFG_varDble("sim_initialEfield")
   !
   !       velocity = sqrt(2.d0 * electronEnergy * electronVolt / electronMass)
   !       acceleration = abs(efield(3) * electronCharge / electronMass )
   !
   !       the root of equation V0*t + 1/2* a* t^2 = GL_gridLength(3) / 2
   !       time = (-2.d0 * velocity + sqrt ((2.d0 * velocity) **2 + 4.d0 * acceleration * GL_gridLength(3))) &
   !                & / (2.d0 * acceleration)
   !
   !
   !    end subroutine timeCheckForAddBGE

   !> Anbang: add electrons into domain, to make sure the background density is constant
   ! subroutine addBackGroundDens(time,densRate)

   !    double precision :: electronEnergy
   !    double precision,dimension(3) :: efield
   !    integer :: ll
   !    DOUBLE PRECISION, DIMENSION(3) :: pos
   !    double precision, intent(in)  :: time, densRate
   !    double precision, save :: accumulateElec = 0.0d0

   !    electronEnergy    = CFG_varDble("init_electronEnergy")
   !    efield = CFG_varDble("sim_initialEfield")


   !    accumulateElec = accumulateElec + time * densRate * GL_gridVolume * PRODUCT(GL_gridSize)

   !    print *, "New background electrons generated:", accumulateElec

   !    if (accumulateElec > 10.d0) then

   !       DO ll = 1, 10

   !          pos = (/kiss_rand(), kiss_rand(), kiss_rand()/)
   !          pos = pos * GL_gridLength

   !          ! Add the background from the upper or lower boundary
   !          if (efield(3) > 0.d0) then
   !             pos(3) = GL_gridLength(3) - GL_gridDelta(3)
   !          else
   !             pos(3) = GL_gridDelta(3)
   !          end if

   !          IF (.NOT. isOutOfGas(pos)) THEN
   !             CALL createIonPair(pos, electronEnergy, 1)
   !          END IF

   !       END DO

   !       accumulateElec = accumulateElec - 10.d0

   !    end if

   ! end subroutine addBackGroundDens

end module module_initialCond
