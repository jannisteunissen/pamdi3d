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
  private

  integer, parameter             :: dp            = kind(0.0d0)
  integer                        :: GAS_num_gases = 0
  real(dp), allocatable          :: GAS_comp_fracs(:)
  character(len=20), allocatable :: GAS_comp_names(:)

  real(dp) :: GAS_pressure, GAS_temperature, GAS_number_dens

  public :: GAS_initialize
  public :: GAS_get_fraction
  public :: GAS_get_num_gases
  public :: GAS_get_pressure
  public :: GAS_get_temperature
  public :: GAS_get_number_dens

contains

  subroutine GAS_initialize(comp_names, comp_fracs, pressure, temperature)
    use m_units_constants

    character(len=*), intent(in) :: comp_names(:)
    real(dp), intent(in) :: comp_fracs(size(gas_comp_names)), pressure, temperature

    GAS_num_gases = size(gas_comp_names)
    allocate( GAS_comp_fracs(GAS_num_gases) )
    allocate( GAS_comp_names(GAS_num_gases) )

    GAS_comp_names = comp_names
    GAS_comp_fracs = comp_fracs
    GAS_pressure = pressure
    GAS_temperature = temperature

    ! Ideal gas law, pressure is in bar
    GAS_number_dens = 1.0D5 * GAS_pressure / (UC_boltzmann_const * GAS_temperature)

  end subroutine GAS_initialize

  real(dp) function GAS_get_fraction(comp_name)
    character(LEN=*), intent(IN) :: comp_name
    integer :: ix

    do ix = 1, GAS_num_gases
       if (comp_name == trim(GAS_comp_names(ix))) exit
    end do

    if (ix == GAS_num_gases + 1) then
       print *, "GAS_get_fraction: " // comp_name // " not found"
       stop
    else
       GAS_get_fraction = GAS_comp_fracs(ix)
    end if
  end function GAS_get_fraction

  integer function GAS_get_num_gases()
    GAS_get_num_gases = GAS_num_gases
  end function GAS_get_num_gases

  real(dp) function GAS_get_pressure()
    GAS_get_pressure = GAS_pressure
  end function GAS_get_pressure

  real(dp) function GAS_get_temperature()
    GAS_get_temperature = GAS_temperature
  end function GAS_get_temperature

  real(dp) function GAS_get_number_dens()
    GAS_get_number_dens = GAS_number_dens
  end function GAS_get_number_dens

end module m_gas
