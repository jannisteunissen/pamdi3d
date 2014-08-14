program test_m_cross_sec
   use m_cross_sec

   integer, parameter :: dp = kind(0.0d0)
   type(CS_t), allocatable :: cross_secs(:)

   character(len=*), parameter :: first_in = "test_m_cross_sec_input.txt"
   character(len=*), parameter :: first_out = "test_m_cross_sec_all_output.txt"
   character(len=*), parameter :: first_summ = "test_m_cross_sec_summary_output.txt"
   character(len=*), parameter :: second_out = "test_m_cross_sec_all_2_output.txt"
   character(len=*), parameter :: gas_1 = "Ar", gas_2 = "N2", gas_3 = "O2"

   call CS_add_from_file(first_in, gas_1, 1.0_dp, 1.0e3_dp, cross_secs)
   call CS_add_from_file(first_in, gas_2, 1.0_dp, 1.0e3_dp, cross_secs)
   call CS_add_from_file(first_in, gas_3, 1.0_dp, 1.0e3_dp, cross_secs)

   call CS_write_summary(cross_secs, first_summ)
   print *, "First pass: read in ", size(cross_secs), " cross sections"

   deallocate(cross_secs)
end program test_m_cross_sec
