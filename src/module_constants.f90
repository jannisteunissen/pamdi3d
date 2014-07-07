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

!> Module that contains physical and numerical constants
MODULE module_constants
   IMPLICIT NONE

   !! Datatype constants

   INTEGER, PARAMETER :: iKind18 = SELECTED_INT_KIND(18)

   !! Numerical constants

   DOUBLE PRECISION, PARAMETER :: pi         =  3.14159265358979324d0
   DOUBLE PRECISION, PARAMETER :: sqrt2      = 1.4142135623730951D0
   DOUBLE PRECISION, PARAMETER :: invSqrt2   = 1.0D0 / sqrt2

   !! Physical constants

   DOUBLE PRECISION, PARAMETER :: epsilon0         = 8.8541878176d-12 ! permitivity for air
   DOUBLE PRECISION, PARAMETER :: electronCharge   = -1.6022d-19     ! the electron charge in Coulombs
   DOUBLE PRECISION, PARAMETER :: elemCharge   = 1.6022d-19     ! the elementary charge in Coulombs
   DOUBLE PRECISION, PARAMETER :: electronVolt     = 1.6022d-19      ! the eV in joules
   DOUBLE PRECISION, PARAMETER :: electronMass     = 9.10938189d-31  ! the electron mass in kg
   DOUBLE PRECISION, PARAMETER :: atomicMass       = 1.66053886D-27  ! the atomic mass unit in kg
   DOUBLE PRECISION, PARAMETER :: N2_mass          = 28.0D0 * atomicMass ! The mass of a N2 molecule

   DOUBLE PRECISION, PARAMETER :: elecChargeOverMass  = electronCharge / electronMass
   DOUBLE PRECISION, PARAMETER :: elecChargeOverEps0  = electronCharge / epsilon0
   DOUBLE PRECISION, PARAMETER :: speedOfLight        = 299792458d0   ! the speed of light in m/s
   DOUBLE PRECISION, PARAMETER :: speedOfLight2       = speedOfLight ** 2   ! the speed of light squared
   DOUBLE PRECISION, PARAMETER :: BoltzmannConstant   = 1.3806503d-23   ! the Boltzmann constant
   DOUBLE PRECISION, PARAMETER :: BohrRadius          = 5.29d-11 !(m), the Bohr radius
   DOUBLE PRECISION, PARAMETER :: RydbergEnergy       = 13.6d0   !(eV), the Rydberg energy
   DOUBLE PRECISION, PARAMETER :: eleCS_Cons          = RydbergEnergy * pi * BohrRadius **2
   DOUBLE PRECISION, PARAMETER :: TorrToBar           = 133.322368 * 1.0D-5 ! one Torr in units of bar
   DOUBLE PRECISION, PARAMETER :: PlanckConstant      = 6.626068d-34   ! Planck constant in unit of m^2 kg /s


   !! Physical order-of-magnitude constants (not precise)

   DOUBLE PRECISION, PARAMETER :: electronMFP      = 2.3d-6  ! mean free path length in m
   DOUBLE PRECISION, PARAMETER :: splittingEnergy  = 27.21d0 * electronVolt    ! a constant in the energy spliting equation

   ! The smallest allowed energy and velocity for collisions
   DOUBLE PRECISION, PARAMETER :: smallEnergyEv =  1.0D-10
   DOUBLE PRECISION, PARAMETER :: smallEnergy   =  smallEnergyEv * electronVolt
   DOUBLE PRECISION, PARAMETER :: smallVelocity =  1.0D-10

   ! Small and large numbers
   DOUBLE PRECISION, PARAMETER :: smallNumber   = EPSILON(1.0D0)
   DOUBLE PRECISION, PARAMETER :: hugeNumber    = HUGE(1.0D0)

END MODULE module_constants
