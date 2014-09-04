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
  subroutine IC_set_init_cond(pc, cfg, rng)
    use m_efield_amr
    use m_random
    use m_config
    use m_phys_domain
    use m_particle_core
    use m_particle
    type(PC_t), intent(inout) :: pc
    type(CFG_t), intent(in) :: cfg
    type(RNG_t), intent(inout) :: rng
    character(LEN=40)      :: initCond
    logical                :: success
    integer                :: ll, ix
    integer                :: nBackgroundPairs
    real(dp) :: init_weight
    real(dp)               :: radius, nIonPairs
    real(dp)               :: temp
    real(dp)               :: backgroundDensity, bgO2MinDensity
    real(dp), dimension(3) :: pos, initSeedPos

    ! Laser line
    real(dp)               :: orthVec1(3), orthVec2(3), lineLen, minFac, &
         maxFac, tempVec(3), lineDirection(3), lineBase(3)


    call CFG_get(cfg, "init_n_ion_pairs", nIonPairs)
    call CFG_get(cfg, "init_weight", init_weight)
    call CFG_get(cfg, "init_bg_dens", backgroundDensity)
    call CFG_get(cfg, "init_o2m_bg_dens", bgO2MinDensity)

    ! Set the background negative ions, for detachment
    call E_set_vars((/E_i_O2m, E_i_pion/), (/bgO2MinDensity, bgO2MinDensity/))

    ! Set the background ionization level
    nBackgroundPairs = nint(backgroundDensity * product(PD_r_max))
    print *, "Background density creates:", nBackgroundPairs
    do ll = 1, nBackgroundPairs
       pos = (/rng%uni_01(), rng%uni_01(), rng%uni_01()/)
       pos = pos * PD_r_max
       if (.not. PD_outside_domain(pos)) call PM_create_ei_pair(pc, pos)
    end do

    print *, "Setting up the initial condition"

    ! Set the initial conditions
    call CFG_get(cfg, "init_cond_type", initCond)

    select case (initCond)
    case ('seed')
       ! Electron/ion pairs start at fixed position
       call CFG_get(cfg, "init_rel_pos", initSeedPos)
       do ll = 1, int (nIonPairs / init_weight)
          pos = initSeedPos * PD_r_max
          call PM_create_ei_pair(pc, pos, w=init_weight)
       end do

    case ('double_seed')
       ! Electron/ion pairs start at fixed position
       call CFG_get(cfg, "init_rel_pos", initSeedPos)
       call CFG_get(cfg, "init_seed_pos_radius", radius)

       do ll = 1, nint(nIonPairs / (2 * init_weight))
          pos = initSeedPos * PD_r_max
          pos(1) = pos(1) - radius
          call PM_create_ei_pair(pc, pos, w=init_weight)
       end do

       do ll = 1, nint(nIonPairs / (2 * init_weight))
          pos = initSeedPos * PD_r_max
          pos(1) = pos(1) + radius
          call PM_create_ei_pair(pc, pos, w=init_weight)
       end do

    case ('Gaussian')
       ! Electron/ion pairs are distributed as a 3D Gaussian distribution
       ! with mean initSeedPos and sigma init_seedPosRadius
       call CFG_get(cfg, "init_seed_pos_radius", radius)
       call CFG_get(cfg, "init_rel_pos", initSeedPos)
       initSeedPos = initSeedPos * PD_r_max

       do ll = 1, nint(nIonPairs / init_weight)
          pos(1:2) = rng%two_normals()
          pos(2:3) = rng%two_normals()
          pos = initSeedPos + pos * radius
          if (.not. PD_outside_domain(pos)) then
             call PM_create_ei_pair(pc, pos, w=init_weight)
          end if
       end do

    case ('laserLine')
       call CFG_get(cfg, "init_seed_pos_radius", radius)
       call CFG_get(cfg, "init_laser_direction", lineDirection)
       call CFG_get(cfg, "init_laser_line_offset", lineBase)

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

       lineLen = norm2((maxFac - minFac) * lineDirection * PD_r_max)
       call CFG_get(cfg, "init_line_dens", temp)
       nIonPairs = temp * radius**2 * lineLen

       print *, "Laser nIonPairs", nIonPairs/init_weight, lineLen, minFac, maxFac
       print *, orthVec1, orthVec2

       do ll = 1, nint(nIonPairs/init_weight)
          pos = (rng%uni_01() * (maxFac-minFac) + minFac) * lineDirection + lineBase
          pos = pos * PD_r_max
          pos = pos + (rng%uni_01()-0.5d0) * radius * orthVec1
          pos = pos + (rng%uni_01()-0.5d0) * radius * orthVec2
          !                print *, ll, pos
          if (.not. PD_outside_domain(pos)) call PM_create_ei_pair(pc, pos, w=init_weight)
       end do

    case ('Gaussian1D')
       ! Electron/ion pairs are distributed as a 1D Gaussian distribution
       ! with mean initSeedPos and sigma init_seedPosRadius
       call CFG_get(cfg, "init_seed_pos_radius", radius)
       call CFG_get(cfg, "init_rel_pos", initSeedPos)
       initSeedPos = initSeedPos * PD_r_max

       do ll = 1, int (nIonPairs / init_weight)
          success = .false.
          do while (.not. success)
             pos(2:3) = rng%two_normals()
             pos(1:2) = (/0.0D0, 0.0D0/)
             pos = initSeedPos + pos * radius
             if (.not. PD_outside_domain(pos)) then
                call PM_create_ei_pair(pc, pos, w=init_weight)
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
