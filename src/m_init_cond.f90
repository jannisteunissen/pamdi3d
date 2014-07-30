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
module m_init_cond

   implicit none
   private

   integer, parameter :: dp = kind(0.0d0)
   
   public :: IC_set_init_cond

contains

   !> Sets up the initial particles for a simulation
  subroutine IC_set_init_cond(rng, grid_length)
    use m_efield_amr
    use m_random
    use m_config
    use m_phys_domain
    type(RNG_t), intent(inout) :: rng
    real(dp), intent(in) :: grid_length(3)
    character(LEN=40)      :: initCond
    logical                :: success
    integer                :: ll, i, j, k, ix
    integer                ::  initWeight, nBackgroundPairs
    real(dp)               :: electronEnergy, radius, nIonPairs
    real(dp)               :: temp
    real(dp)               :: backgroundDensity, bgO2MinDensity
    real(dp), dimension(3) :: pos, initSeedPos

    ! Laser line
    real(dp)               :: orthVec1(3), orthVec2(3), lineLen, minFac, &
         maxFac, tempVec(3), lineDirection(3), lineBase(3)


    nIonPairs = CFG_get_real("init_nIonPairs")
    electronEnergy    = CFG_get_real("init_electronEnergy")
    initWeight        = CFG_get_int("init_weightFactor")
    backgroundDensity = CFG_get_real("init_backgroundDensity")
    bgO2MinDensity    = CFG_get_real("init_O2MinBackgroundDens")

    ! Set the background negative ions, for detachment
    call E_set_vars((/E_i_O2m, E_i_pion/), (/bgO2MinDensity, bgO2MinDensity/))

    ! Set the background ionization level
    nBackgroundPairs = int(backgroundDensity * product(grid_length))
    print *, "Background density creates:", nBackgroundPairs
    do ll = 1, nBackgroundPairs
       pos = (/rng%uni_01(), rng%uni_01(), rng%uni_01()/)
       pos = pos * grid_length
       if (.not. PD_outside_domain(pos)) call createIonPair(pos, electronEnergy, 1)
    end do

    print *, "Setting up the initial condition"

    ! Set the initial conditions
    call CFG_getVar("init_condType", initCond)

    select case (initCond)
    case ('seed')
       ! Electron/ion pairs start at fixed position
       call CFG_getVar("init_relSeedPos", initSeedPos)
       do ll = 1, int (nIonPairs / initWeight)
          pos = initSeedPos * grid_length
          call createIonPair(pos, electronEnergy, initWeight)
       end do

    case ('double_seed')
       ! Electron/ion pairs start at fixed position
       call CFG_getVar("init_relSeedPos", initSeedPos)
       radius = CFG_get_real("init_seedPosRadius")

       do ll = 1, int(nIonPairs / (2 * initWeight))
          pos = initSeedPos * grid_length
          pos(1) = pos(1) - radius
          call createIonPair(pos, electronEnergy, initWeight)
       end do

       do ll = 1, int(nIonPairs / (2 * initWeight))
          pos = initSeedPos * grid_length
          pos(1) = pos(1) + radius
          call createIonPair(pos, electronEnergy, initWeight)
       end do

    case ('Gaussian')
       ! Electron/ion pairs are distributed as a 3D Gaussian distribution
       ! with mean initSeedPos and sigma init_seedPosRadius
       radius = CFG_get_real("init_seedPosRadius")
       call CFG_getVar("init_relSeedPos", initSeedPos)
       initSeedPos = initSeedPos * grid_length

       do ll = 1, int (nIonPairs / initWeight)
          success = .false.
          do while (.not. success)
             pos(1:2) = rng%two_normals()
             pos(2:3) = rng%two_normals()
             pos = initSeedPos + pos * radius
             if (.not. PD_outside_domain(pos)) then
                call createIonPair(pos, electronEnergy, initWeight)
                success = .true.
             end if
          end do
       end do

    case ('laserLine')
       radius = CFG_get_real("init_seedPosRadius")
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
       lineDirection = lineDirection / norm2(lineDirection)
       orthVec1 = orthVec1 / norm2(orthVec1)
       orthVec2 = orthVec2 / norm2(orthVec2)

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

       lineLen = norm2((maxFac - minFac) * lineDirection * grid_length)
       temp = CFG_get_real("init_lineDens")
       nIonPairs = temp * radius**2 * lineLen

       print *, "Laser nIonPairs", nIonPairs/initWeight, lineLen, minFac, maxFac
       print *, orthVec1, orthVec2

       do ll = 1, int(nIonPairs/initWeight)
          pos = (rng%uni_01() * (maxFac-minFac) + minFac) * lineDirection + lineBase
          pos = pos * grid_length
          pos = pos + (rng%uni_01()-0.5d0) * radius * orthVec1
          pos = pos + (rng%uni_01()-0.5d0) * radius * orthVec2
          !                print *, ll, pos
          if (.not. PD_outside_domain(pos)) call createIonPair(pos, electronEnergy, initWeight)
       end do

    case ('Gaussian1D')
       ! Electron/ion pairs are distributed as a 1D Gaussian distribution
       ! with mean initSeedPos and sigma init_seedPosRadius
       radius = CFG_get_real("init_seedPosRadius")
       call CFG_getVar("init_relSeedPos", initSeedPos)
       initSeedPos = initSeedPos * grid_length

       do ll = 1, int (nIonPairs / initWeight)
          success = .false.
          do while (.not. success)
             pos(2:3) = rng%two_normals()
             pos(1:2) = (/0.0D0, 0.0D0/)
             pos = initSeedPos + pos * radius
             if (.not. PD_outside_domain(pos)) then
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

  end subroutine IC_set_init_cond

end module m_init_cond
