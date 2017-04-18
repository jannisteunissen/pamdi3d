
!> electrons move in negative z direction
!> we compute the positive velocity due to photo ionization

program test_m_photoi
  use m_pc_all

  implicit none

  integer, parameter :: dp = kind(0.0d0)
  
  integer                        :: init_num_part = 20000
  integer                        :: max_num_steps = 100
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
  
  real(dp)                    :: efield, wall, av_pos1,av_pos2
  integer                    :: ll, step, n
  

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
  
  call PI_initialize(pc) ! default values are taken defined inside PI_initialize

  print *, "Initializing particle module"
  call pc%initialize(    UC_elec_mass, & ! mass of particle
  &                        cross_secs, & ! cross section list
  &                 lookup_table_size, & ! lookup table size
  &                     max_energy_eV, & ! max energy (eV)
  &                   max_n_particles)  ! max number of particles
  pc%outside_check => part_outside_check

  
  efield = 4.0d6
  wall = 0.0d0
  do ll = 1, init_num_part
    call pc%create_part( [ 0.0D0, 0.0D0, 2.0D-5], & ! position (m)
    &                    [ 0.0D0, 0.0D0, 0.0D0], & ! velocity (m/s)
    &                    [ 0.0_dp, 0.0_dp, efield*UC_elec_q_over_m], & ! accelaration (m/s2)
    &                                       1.0D0, & ! weight
    &                                       0.0D0)   ! time left for move_and_collide (s)
  end do

  if (pc%n_part>0.0d0) write(*,'(A10,A10,A10,A10,A10)') &
  "Time (ns)", "Number","v (mm/ns)", "zmin (mm)"
  do step = 1, max_num_steps
     if (pc%n_part>1) then
       av_pos1 = sum(pc%particles(1:pc%n_part)%x(3))/pc%n_part
       call pc%advance(delta_t)
       av_pos2 = sum(pc%particles(1:pc%n_part)%x(3))/pc%n_part
       write(*,'(F10.2,I10,F10.2,F10.2)') &
       step*delta_t*1.0d9, pc%n_part, (av_pos2-av_pos1)/delta_t*1.0d-6, &
       minval(pc%particles(1:pc%n_part)%x(3))*1.0d3
     end if
  end do


contains
  logical function part_outside_check(my_part)
    implicit none
    type(PC_part_t), intent(in) :: my_part

    part_outside_check = .false.
    if ( my_part%x(3)<wall) part_outside_check=.true.
  end function part_outside_check
end program test_m_photoi