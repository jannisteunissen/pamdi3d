program test_m_particle_core
   use m_particle_core
   use m_cross_sec
   use m_units_constants

   implicit none

   integer, parameter :: dp = kind(0.0d0)
   character(len=*), parameter :: cs_file = "test_m_particle_core_cs.txt"
   character(len=*), parameter :: gas_name = "druyv_gas"

   integer, parameter         :: max_num_part  = 1000*1000
   integer, parameter         :: init_num_part = 10*1000
   integer, parameter         :: max_num_steps = 100
   integer, parameter         :: lkp_tbl_size  = 1000
   integer, parameter         :: num_lists = 12
   integer, parameter :: n_bins = 100
   real(dp), parameter        :: delta_t       = 5.0e-9_dp
   real(dp), parameter        :: max_en_eV     = 1.0e-1_dp
   real(dp), parameter        :: neutral_dens  = 2.5e25_dp
   real(dp), parameter        :: part_mass     = UC_elec_mass
   real(dp), parameter        :: init_accel(3) = (/0.0_dp, 0.0_dp, 1.0e12_dp/)
   real(dp)                   :: norm_cross_sec, mass_ratio
   real(dp)                   :: pos(3), vel(3), accel(3), weight
   real(dp) :: bin_args(1)
   integer                    :: ll, step, num_colls, i
   type(CS_t), allocatable :: cross_secs(:)
   integer, parameter :: n_pc = 4
   type(PC_t) :: pmodels(n_pc)

   print *, "Testing m_particle_core.f90 implementation"

   ! Use constant momentum transfer cross section, so that we get a Druyvesteyn
   ! distribution
   print *, "Reading in cross sections from ", trim(cs_file)
   call CS_add_from_file(cs_file, gas_name, 1.0_dp, neutral_dens, &
        1.0e3_dp, cross_secs)

   ! All cross sections should be constant, simply take the first value
   norm_cross_sec = 0
   do ll = 1, size(cross_secs)
      norm_cross_sec = norm_cross_sec + cross_secs(ll)%en_cs(2,1)
   end do
   print *, "Total cross sec", norm_cross_sec
   mass_ratio = cross_secs(1)%coll%rel_mass

   print *, "Initializing particle module"
   print *, part_mass
   do i = 1, n_pc
      call pmodels(i)%initialize(part_mass, cross_secs, lkp_tbl_size, &
           max_en_eV, max_num_part)
   end do
   
   num_colls = pmodels(1)%get_num_colls()
   deallocate(cross_secs)

   print *, "Creating initial particles"
   do ll = 1, init_num_part
      pos = 0.0_dp
      vel = 0.0_dp
      accel = init_accel
      weight = 1
      ! i = 1 + mod(ll, n_pc)
      call pmodels(1)%create_part(pos, vel, accel, weight, 0.0_dp)
   end do

   call PC_share_particles(pmodels)
   ! call PC_reorder_by_bins(pmodels, bin_func, n_bins, bin_args)

   do step = 1, max_num_steps
      print *, ""
      print *, "at step", step, " and time ", (step-1) * delta_t
      call print_stats()

      !$omp parallel do
      do i = 1, n_pc
         call pmodels(i)%advance(delta_t)
      end do
      !$omp end parallel do
   end do

   call print_stats()
   ! call pmodel%merge_and_split((/0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/), &
   !      1.0e-6_dp, get_weight_2, PC_merge_part_rxv, PC_split_part)
   ! call print_stats()
   ! call pmodel%merge_and_split((/0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/), &
   !      1.0e-6_dp, get_weight_2, PC_merge_part_rxv, PC_split_part)
   ! call print_stats()

contains

   subroutine part_stats(part, vec)
      type(PC_part_t), intent(in) :: part
      real(dp), intent(out) :: vec(:)
      vec(1:3) = part%w * part%x
      vec(4:6) = part%w * part%v
      vec(7:9) = part%w * part%a
      vec(10) = part%w * PC_v_to_en(part%v, part_mass)
      vec(11) = part%w
   end subroutine part_stats

   subroutine print_stats()
      integer :: n_part
      real(dp) :: sum_x(3), sum_v(3), sum_a(3), sum_en
      real(dp) :: sum_weight
      real(dp) :: sum_vec(11), tmp_vec(11)

      n_part = 0
      sum_vec = 0
      
      do i = 1, n_pc
         n_part = n_part + pmodels(i)%get_num_sim_part()
         call pmodels(i)%compute_vector_sum(part_stats, tmp_vec)
         sum_vec = sum_vec + tmp_vec
      end do
      
      sum_weight = sum_vec(11)
      sum_x = sum_vec(1:3)
      sum_v = sum_vec(4:6)
      sum_a = sum_vec(7:9)
      sum_en =  sum_vec(10)

      ! print *, "mean position", sum_x / sum_weight
      ! print *, "mean velocity", sum_v / sum_weight
      print *, "mean energy (eV)              ", sum_en / (sum_weight * UC_elec_volt)
      print *, "Druyvesteyn mean energy: (eV) ", 0.5_dp * UC_elec_mass * 0.739669_dp / &
           (sqrt(3*norm_cross_sec**2 * (2 / (1 + 1/mass_ratio)) / (8 * init_accel(3)**2)) *  UC_elec_volt)
      print *, "Number of particles           ", n_part, sum_weight
   end subroutine print_stats

   real(dp) function get_weight_1(my_part)
      type(PC_part_t), intent(in) :: my_part
      get_weight_1 = 1
   end function get_weight_1

   real(dp) function get_weight_2(my_part)
      type(PC_part_t), intent(in) :: my_part
      get_weight_2 = 2
   end function get_weight_2

   integer function bin_func(my_part, bin_args)
     type(PC_part_t), intent(in) :: my_part
     real(dp), intent(in) :: bin_args(:)
     real(dp) :: tmp
     call random_number(tmp)
     bin_func = 1 + tmp * n_bins
     ! print *, bin_func
   end function bin_func

 end program test_m_particle_core
