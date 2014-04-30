program test_m_particle_core
   use m_particle_core
   use m_cross_sec
   use m_units_constants

   implicit none
   integer, parameter :: dp = kind(0.0d0)
   character(len=*), parameter :: cs_file = "test_m_particle_core_cs.txt"
   character(len=*), parameter :: gas_name = "druyv_gas"

   integer, parameter :: max_num_part = 5*1000*1000
   integer, parameter :: init_num_part = 1000
   integer, parameter :: max_num_steps = 10
   integer, parameter :: lkp_tbl_size = 10*1000
   real(dp), parameter :: delta_t = 1.0e-8_dp
   real(dp), parameter :: max_en_eV = 1.0e-1_dp
   real(dp), parameter :: neutral_dens = 2.5e25_dp
   real(dp), parameter :: part_mass = UC_elec_mass
   real(dp), parameter :: init_accel(3) = (/0.0_dp, 0.0_dp, 1.0e12_dp/)
   real(dp) :: norm_cross_sec, mass_ratio
   real(dp) :: pos(3), vel(3), accel(3), weight
   integer :: ll, step, num_colls
   type(CS_type), allocatable :: cross_secs(:)

   print *, "Testing m_particle_core.f90 implementation"

   ! Use constant momentum transfer cross section, so that we get a Druyvesteyn distribution
   print *, "Reading in cross sections from ", trim(cs_file)
   call CS_read_file(cs_file, gas_name, 1.0_dp, neutral_dens, 1.0e3_dp)
   call CS_get_cross_secs(cross_secs)
   call CS_write_summary("test_m_particle_core_cs_summ.txt")
   norm_cross_sec = cross_secs(1)%en_cs(2,1) ! First entry of first (and only) cross sec
   mass_ratio = cross_secs(1)%spec_value

   print *, "Initializing particle module"
   print *, part_mass
   call PC_initialize(part_mass, cross_secs, lkp_tbl_size, max_en_eV)
   num_colls = PC_get_num_colls()
   deallocate(cross_secs)

   print *, "Creating initial particles"
   do ll = 1, init_num_part
      pos = 0.0_dp
      vel = 0.0_dp
      accel = init_accel
      weight = 1
      call PC_create_part(pos, vel, accel, weight, 0.0_dp)
   end do

   do step = 1, max_num_steps
      print *, ""
      print *, "at step", step, " and time ", (step-1) * delta_t
      call print_stats()
      call PC_advance(delta_t)
   end do

   call print_stats()
   call PC_merge_and_split((/0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/), &
        1.0e-6_dp, get_weight_2, merge_particles, split_particles)
   call print_stats()
   call PC_merge_and_split((/0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/), &
        1.0e-6_dp, get_weight_1, merge_particles, split_particles)
   call print_stats()

contains

   subroutine print_stats()
      integer :: nn, n_part
      type(PC_part_t) :: my_part
      real(dp) :: sum_x(3), sum_v(3), sum_a(3), sum_en
      real(dp) :: coll_count(num_colls), sum_weight

      sum_x = 0.0_dp
      sum_v = 0.0_dp
      sum_a = 0.0_dp
      sum_en = 0.0_dp
      sum_weight = 0.0_dp

      n_part = PC_get_num_sim_part()
      do nn = 1, n_part
         call PC_get_part(nn, my_part)
         sum_x = sum_x + my_part%x * my_part%weight
         sum_v = sum_v + my_part%v * my_part%weight
         sum_a = sum_a + my_part%a * my_part%weight
         sum_en = sum_en + sum(my_part%v**2) * my_part%weight
         sum_weight = sum_weight + my_part%weight
      end do

      call PC_get_coll_count(coll_count)

      print *, "mean position", sum_x / sum_weight
      print *, "mean velocity", sum_v / sum_weight
      print *, "mean energy (eV)              ", 0.5_dp * part_mass * sum_en / (sum_weight * UC_elec_volt)
      print *, "Druyvesteyn mean energy: (eV) ", 0.5_dp * UC_elec_mass * 0.739669_dp / &
           (sqrt(3*norm_cross_sec**2 * (2 / (1 + 1/mass_ratio)) / (8 * init_accel(3)**2)) *  UC_elec_volt)
      print *, "Number of particles           ", n_part
      print *, "Total weight                  ", sum_weight
      print *, "collision count               ", coll_count
   end subroutine print_stats

   real(dp) function get_weight_1(my_part)
      type(PC_part_t), intent(in) :: my_part
      get_weight_1 = 1
   end function get_weight_1

   real(dp) function get_weight_2(my_part)
      type(PC_part_t), intent(in) :: my_part
      get_weight_2 = 2
   end function get_weight_2

   subroutine merge_particles(part_a, part_b)
      type(PC_part_t), intent(inout) :: part_a, part_b
      part_a%weight = part_a%weight + part_b%weight
      part_b%live = .false.
   end subroutine merge_particles

   subroutine split_particles(part_a, part_b)
      type(PC_part_t), intent(inout) :: part_a, part_b
      real(dp) :: old_weight
      old_weight = part_a%weight
      part_a%weight = (old_weight + 1)/2
      part_b = part_a
      part_b%weight = old_weight/2
   end subroutine split_particles

end program test_m_particle_core
