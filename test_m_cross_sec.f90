program test_m_cross_sec
   use m_cross_sec

   integer, parameter :: dp = kind(0.0d0)
   type(CS_type), allocatable :: cross_secs(:)

   character(len=*), parameter :: first_in = "test_m_cross_sec_input.txt"
   character(len=*), parameter :: first_out = "test_m_cross_sec_all_output.txt"
   character(len=*), parameter :: first_summ = "test_m_cross_sec_summary_output.txt"
   character(len=*), parameter :: second_out = "test_m_cross_sec_all_2_output.txt"
   character(len=*), parameter :: gas_1 = "Ar", gas_2 = "N2", gas_3 = "O2"

   call CS_read_file(first_in, gas_1, 1.0_dp, 1.0_dp, 1.0e3_dp)
   call CS_read_file(first_in, gas_2, 1.0_dp, 1.0_dp, 1.0e3_dp)
   call CS_read_file(first_in, gas_3, 1.0_dp, 1.0_dp, 1.0e3_dp)
   call CS_write_summary(first_summ)
   call CS_write_all(first_out)
   call CS_get_cross_secs(cross_secs)
   print *, "First pass: read in ", size(cross_secs), " cross sections"

   call CS_reset()
   deallocate(cross_secs)

   call CS_read_file(first_out, gas_1, 1.0_dp, 1.0_dp, 1.0e3_dp)
   call CS_read_file(first_out, gas_2, 1.0_dp, 1.0_dp, 1.0e3_dp)
   call CS_read_file(first_out, gas_3, 1.0_dp, 1.0_dp, 1.0e3_dp)
   call CS_write_all(second_out)

   call CS_get_cross_secs(cross_secs)
   print *, "Second pass: read in ", size(cross_secs), " cross sections"
   print *, "The following two files should be the same:"
   print *, first_out
   print *, second_out
   deallocate(cross_secs)

end program test_m_cross_sec
