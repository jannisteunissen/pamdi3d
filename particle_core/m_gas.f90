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

!> Module for gas properties (type, fraction, pressure, temperature etc)
module m_gas

  implicit none
  public

  integer, parameter, private               :: dp            = kind(0.0d0)
  integer, protected                        :: GAS_num_gases = 0
  real(dp), allocatable, protected          :: GAS_comp_fracs(:)
  character(len=20), allocatable, protected :: GAS_comp_names(:)
  real(dp), protected                       :: GAS_pressure
  real(dp), protected                       :: GAS_temperature
  real(dp), protected                       :: GAS_number_dens

  public :: GAS_initialize
  public :: GAS_get_fraction

contains

  subroutine GAS_initialize(comp_names, comp_fracs, pressure, temperature)
    use m_units_constants
    character(len=*), intent(in) :: comp_names(:)
    real(dp), intent(in)         :: comp_fracs(:), pressure, temperature
    real(dp), parameter :: err_threshold = 1.0e-4_dp

    GAS_num_gases   = size(comp_names)
    GAS_pressure    = pressure
    GAS_temperature = temperature
    allocate(GAS_comp_fracs(GAS_num_gases))
    allocate(GAS_comp_names(GAS_num_gases))
    GAS_comp_names  = comp_names
    GAS_comp_fracs  = comp_fracs

    if (abs(sum(GAS_comp_fracs) - 1) > err_threshold) then
       print *, "Gas comps are not normalized: ", comp_fracs
       stop
    end if

    ! Ideal gas law, pressure is in bar
    GAS_number_dens = 1.0D5 * GAS_pressure / &
         (UC_boltzmann_const * GAS_temperature)
  end subroutine GAS_initialize

  real(dp) function GAS_get_fraction(comp_name)
    character(LEN=*), intent(IN) :: comp_name
    integer                      :: ix

    do ix = 1, GAS_num_gases
       if (comp_name == trim(GAS_comp_names(ix))) then
          GAS_get_fraction = GAS_comp_fracs(ix)
          return
       end if
    end do

    print *, "GAS_get_fraction: " // comp_name // " not found"
    GAS_get_fraction = 0.0_dp
    stop
  end function GAS_get_fraction

end module m_gas
