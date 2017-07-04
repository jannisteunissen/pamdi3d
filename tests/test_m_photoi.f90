!  Casper Rutjes
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

!> Program to test the photo ionization module

!> Set-up:
!> Electrons move to negative z and are delected at z=0.
!> Photo ionizations create new electrons at positive z locations
!> which will move again to negative z and disappear at z=0.0.

program test_m_photoi
  use m_pc_all

  implicit none

  integer, parameter :: dp = kind(0.0d0)

  integer                        :: init_num_part = 1000
  integer                        :: max_num_steps = 10000
  real(dp)                       :: delta_t = 1.0d-11
  integer                        :: lookup_table_size = 10000
  real(dp)                       :: max_energy_eV = 1000.0D0
  integer                        :: max_n_particles = 100000
  real(dp)                       :: pressure = 1.0d0
  real(dp)                       :: temperature = 293.0d0
  real(dp), allocatable          :: gas_fracs(:)
  character(len=80), allocatable :: molecule_names(:), cross_files(:)
  integer                        :: n_gas_comp = 2
  character(len=*), parameter     :: cs_file = "test_m_photoi_cs.txt"


  type(CS_t), allocatable        :: cross_secs(:)
  type(PC_t)                     :: pc

  real(dp)                    :: efield
  integer                    :: ll, step, n

  real(dp)                               :: quench_fac
  real(dp)                               :: min_inv_abs_len_resc, max_inv_abs_len_resc
  integer                                :: size_photo_eff_table
  real(dp), dimension(:, :), allocatable :: photo_eff_table


  print *, "Testing m_photoi.f90 implementation"

  allocate(molecule_names(n_gas_comp))
  molecule_names = ["N2", "O2"]
  allocate(gas_fracs(n_gas_comp))
  gas_fracs = [0.8d0, 0.2d0]
  allocate(cross_files(n_gas_comp))
  cross_files = ["test_m_photoi_cs.txt", "test_m_photoi_cs.txt"]

  ! Load GAS and cross sections, needed for photoi
  call GAS_initialize(molecule_names, gas_fracs, pressure, temperature)
  do n = 1, n_gas_comp
     call CS_add_from_file(cs_file, molecule_names(n), &
          gas_fracs(n) * GAS_number_dens, max_energy_eV, cross_secs)
  end do

  min_inv_abs_len_resc = 2.0D3 / GAS_get_fraction("O2") / GAS_pressure
  max_inv_abs_len_resc = 2.0D3 / GAS_get_fraction("O2") / GAS_pressure

  size_photo_eff_table = 2
  allocate(photo_eff_table(2,size_photo_eff_table))
  photo_eff_table(1,:) = [0.0d0,   1d8] ! for any field
  photo_eff_table(2,:) = [1.0d0, 1.0d0] ! an efficiency 1 per ionization
  quench_fac = 1.0d0

  call PI_initialize(quench_fac,           &
                     min_inv_abs_len_resc, &
                     max_inv_abs_len_resc, &
                     size_photo_eff_table, &
                     photo_eff_table)

  call pc%add_ionization_callback(PI_do_photoi)

  print *, "Initializing particle module"
  call pc%initialize(    UC_elec_mass, & ! mass of particle
  &                        cross_secs, & ! cross section list
  &                 lookup_table_size, & ! lookup table size
  &                     max_energy_eV, & ! max energy (eV)
  &                   max_n_particles)  ! max number of particles
  pc%outside_check => part_outside_check


  efield = 3.5d6
  do ll = 1, init_num_part
    call pc%create_part( [ 0.0D0, 0.0D0, 1.0D-4], & ! position (m)
    &                    [ 0.0D0, 0.0D0, 0.0D0], & ! velocity (m/s)
    &                    [ 0.0_dp, 0.0_dp, efield*UC_elec_q_over_m], & ! accelaration (m/s2)
    &                                       1.0D0, & ! weight
    &                                       0.0D0)   ! time left for move_and_collide (s)
  end do

  if (pc%n_part>0.0d0) write(*,'(A14,A14,A14,A14,A14)') &
  "Time (ns)", "Number","zmin (mm)", "z_av (mm)", "zmax (mm)"
  do step = 1, max_num_steps
     if (pc%n_part>1 .and. mod(step,10).eq.0) then

       call pc%advance(delta_t)

       write(*,'(F14.1,I14,F14.3,F14.3,F14.3)') &
       step*delta_t*1.0d9, pc%n_part, &
       minval(pc%particles(1:pc%n_part)%x(3))*1.0d3, &
       sum(pc%particles(1:pc%n_part)%x(3))/pc%n_part,&
       maxval(pc%particles(1:pc%n_part)%x(3))*1.0d3
     end if
  end do


contains
  logical function part_outside_check(my_part)
    implicit none
    type(PC_part_t), intent(in) :: my_part

    part_outside_check = .false.
    if ( my_part%x(3)<0.0d0) part_outside_check=.true.
  end function part_outside_check
end program test_m_photoi
