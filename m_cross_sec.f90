! Copyright 2005-2012, Chao Li, Anbang Sun, Jannis Teunissen
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

!> Module that contains routines to read in cross section data from textfiles.
module m_cross_sec

   implicit none
   private

   integer, parameter :: dp = kind(0.0d0)
   integer, parameter :: tiny_len = 20
   integer, parameter :: name_len = 40
   integer, parameter :: line_len = 200

   !> The type of cross section table
   type CS_type
      real(dp), allocatable   :: en_cs(:,:)    ! Stores the energy vs cross sec table
      integer                 :: n_rows        ! Number of rows in the table
      integer                 :: col_type      ! Type of collision, defined below
      real(dp)                :: spec_value    ! Specific value that can be used
      real(dp)                :: min_energy    ! Minimum energy for the collision
      real(dp)                :: max_energy    ! Maximum energy for the collision
      character(LEN=tiny_len) :: gas_name      ! Name of the colliding neutral molecule
      character(LEN=line_len) :: description   ! Description of the collision
      character(LEN=line_len) :: comment       ! Additional comments
   end type CS_type

   type(CS_type), allocatable :: CS_table(:)

   integer, parameter :: CS_attach_t = 1, CS_elastic_t = 2, CS_excite_t = 3, CS_ionize_t = 4
   integer, parameter :: max_num_cols_per_gas = 50
   integer, parameter :: max_num_rows = 400
   integer            :: CS_num_colls_found = 0

   ! Public variables
   public :: CS_type
   public :: CS_attach_t, CS_elastic_t, CS_excite_t, CS_ionize_t

   ! Methods
   public :: CS_reset
   public :: CS_get_cross_secs
   public :: CS_read_file
   public :: CS_write_summary
   public :: CS_write_all

contains

   ! Reset all the cross sections that have been found
   subroutine CS_reset()
      deallocate(CS_table)
      CS_num_colls_found = 0
   end subroutine CS_reset

   ! Return a copy of the cross sections that have been found
   subroutine CS_get_cross_secs(cross_secs)
      type(CS_type), allocatable :: cross_secs(:)
      allocate(cross_secs(CS_num_colls_found))
      cross_secs = CS_table(1:CS_num_colls_found)
   end subroutine CS_get_cross_secs

   ! Search 'filename' for cross section data concerning 'gas_name'
   subroutine CS_read_file(filename, gas_name, x_normalization, y_normalization, req_energy)
      character(LEN=*), intent(IN) :: gas_name, filename
      real(dp), intent(in)         :: x_normalization, y_normalization, req_energy
      integer                      :: n, cIx, nL, n_rows, col_type
      integer                      :: my_unit, io_state, len_gas_name
      character(LEN=name_len)      :: lineFMT
      character(LEN=line_len)      :: line, prev_line
      real(dp)                     :: tempArray(2, max_num_rows)
      real(dp)                     :: x_scaling, y_scaling

      my_unit      = 333
      nL           = 0 ! Set the number of lines to 0
      len_gas_name = len(trim(gas_name))

      ! Set the line format to read, only depends on line_len currently
      write(lineFMT, FMT = "(A,I0,A)") "(A", line_len, ")"

      ! Open 'filename' (with error checking)
      open(my_unit, FILE = filename, STATUS = "OLD", ACTION = "READ", ERR = 999, IOSTAT = io_state)

      ! Look for collision processes with the correct gas name in the file,
      ! which should look for example like:

      !     ATTACHMENT                    [description of the type of process, always in CAPS]
      !     H2O -> H2O^-                  [the gas name possibly followed by the result of the process]
      !     COMMENT: total attachment     [possibly comments]
      !     UPDATED: 2010-06-24 15:04:36
      !     SCALING: 1.0 1.0              [optionally scale factors for the columns]
      !     ------------------            [at least 5 dashes]
      !     xxx   xxx                     [cross section data in two column format]
      !     ...   ...
      !     xxx   xxx
      !     ------------------

      ! So if we find the gas name the previous line holds the type of collision, then
      ! there is possibly some extra information and between the dashes the actual cross
      ! sections are found.

      ! The outer DO loop, running until the end of the file is reached
      do
         ! Search for 'gas_name' in the file
         line = ' '
         do
            prev_line = line
            read(my_unit, FMT = lineFMT, ERR = 999, end = 666) line; nL = nL+1
            line = adjustl(line)
            if (line(1:len_gas_name) == gas_name) exit
         end do

         ! Check prev_line for the type of collision
         select case (prev_line)
         case ("ATTACHMENT")
            col_type = CS_attach_t
         case ("ELASTIC", "MOMENTUM")
            col_type = CS_elastic_t
         case ("EFFECTIVE")
            col_type = CS_elastic_t
            print *, "CS_read_file warning: using EFFECTIVE elastic cross section for ", &
                 gas_name, ", this should not be used in particle simulations."
         case ("EXCITATION")
            col_type = CS_excite_t
         case ("IONIZATION")
            col_type = CS_ionize_t
         case ("COMMENT")
            cycle
         case DEFAULT
            print *, "CS_read_file warning: ignoring unknown process type for ", &
                 gas_name, " in ", filename, " at line ", nL
            cycle
         end select

         ! Update the number of processes and set the gas name and collision type
         CS_num_colls_found = CS_num_colls_found + 1
         cIx = CS_num_colls_found
         call ensure_free_storage(cIx)

         ! Add the reaction description to the table
         CS_table(cIx)%description = adjustl(line)

         ! For all collisions except attachment, there is a specific value on the next line
         if (col_type /= CS_attach_t) then
            read(my_unit, FMT = lineFMT, ERR = 999, end = 555) line; nL = nL+1
            read(line, FMT = *, ERR = 999, end = 555) CS_table(cIx)%spec_value
         else
            CS_table(cIx)%spec_value = 0.0_dp
         end if

         CS_table(cIx)%gas_name = gas_name
         CS_table(cIx)%col_type = col_type
         CS_table(cIx)%comment = "COMMENT: (empty)"
         x_scaling = 1.0_dp
         y_scaling = 1.0_dp

         ! Now we can check whether there is a ZDPLASKIN statement and a reaction description,
         ! while scanning lines until dashes are found, which indicate the start of the cross section data
         do
            read(my_unit, FMT = lineFMT, ERR = 999, end = 555) line; nL = nL+1
            line = adjustl(line)
            if ( line(1:9) == "ZDPLASKIN" ) then
               CS_table(cIx)%description = trim(gas_name) // " [" // trim(adjustl(line(11:))) // "]"
            else if ( line(1:7) == "COMMENT") then
               CS_table(cIx)%comment = line
            else if ( line(1:7) == "SCALING") then
               read(line(9:), *) x_scaling, y_scaling
            else if ( line(1:5) == "-----" ) then
               exit
            end if
         end do

         ! Read the cross section data into a temporary array
         n_rows = 0
         do
            read(my_unit, FMT = lineFMT, ERR = 999, end = 555) line; nL = nL+1
            line = adjustl(line)
            if ( line(1:5) == "-----" ) then
               exit  ! Dashes mark the end of the data
            else if (trim(line) == "" .or. line(1:1) == "#") then
               cycle ! Ignore whitespace or comments
            else if (n_rows < max_num_rows) then
               n_rows = n_rows + 1
               read(line, FMT = *, ERR = 999, end = 555) tempArray(:, n_rows)
            else
               print *, "CS_read_file error: too many rows in ", filename, " at line ", nL
               stop
            end if
         end do

         if (n_rows < 2) then
            print *, "CS_read_file error: need at least two values in ", &
              filename, " at line number ", nL
            stop
         end if

         ! Store the data in the actual table
         allocate( CS_table(cIx)%en_cs(2, n_rows) )
         CS_table(cIx)%n_rows = n_rows
         CS_table(cIx)%en_cs(1,:) = tempArray(1, 1:n_rows) * x_normalization * x_scaling
         CS_table(cIx)%en_cs(2,:) = tempArray(2, 1:n_rows) * y_normalization * y_scaling
         CS_table(cIx)%max_energy = CS_table(cIx)%en_cs(1, n_rows)

         ! Locate minimum energy (first value followed by non-zero cross sec)
         CS_table(cIx)%min_energy = 0.0_dp
         do n = 1, n_rows-1
            if (CS_table(cIx)%en_cs(2, n+1) > 0.0_dp) then
               CS_table(cIx)%min_energy = CS_table(cIx)%en_cs(1, n)
               exit
            end if
         end do

         ! Locate maximum energy (last value preceded by non-zero)
         CS_table(cIx)%max_energy = 0.0_dp
         do n = n_rows, 2, -1
            if (CS_table(cIx)%en_cs(2, n-1) > 0.0_dp) then
               CS_table(cIx)%max_energy = CS_table(cIx)%en_cs(1, n)
               exit
            end if
         end do

         ! Check whether the tables that have been read in go up to high enough energies for our
         ! simulation. They can also have 0.0 as their highest listed cross section, in which
         ! case we assume the cross section is 0.0 for all higher energies

         if ( CS_table(cIx)%en_cs(1, n_rows) < req_energy .and. &
              & CS_table(cIx)%en_cs(2, n_rows) > 0.0D0 ) then
            print *, "CS_read_file error: cross section data at line ", nL, &
                 " does not go up to high enough x-values (energy)"
            stop
         end if

      end do

555   continue ! Routine ends here if the end of "filename" is reached erroneously
      close(my_unit, ERR = 999, IOSTAT = io_state)
      print *, "CS_read_file error, reached end of file while searching. ", &
           "io_state = ", io_state, " while reading from [", filename, "] at line ", nL
      stop
      return

666   continue ! Routine ends here if the end of "filename" is reached correctly
      close(my_unit, ERR = 999, IOSTAT = io_state)
      return

999   continue ! If there was an error, the routine will end here
      print *, "CS_read_file error at line ", nL, &
           " io_state = ", io_state, " while searching [", gas_name, "] in [", filename, "]"
      stop

   end subroutine CS_read_file

   subroutine ensure_free_storage(req_size)
      integer, intent(in) :: req_size
      type(CS_type), allocatable :: tbl_copy(:)
      integer :: curr_size

      if (allocated(CS_table)) then
         curr_size = size(CS_table)
         if (curr_size < req_size) then
            allocate(tbl_copy(curr_size))
            tbl_copy = CS_table
            deallocate(CS_table)
            allocate(CS_table(max(req_size, 2 * curr_size)))
            CS_table(1:curr_size) = tbl_copy
         end if
      else
         allocate(CS_table(req_size))
      end if
   end subroutine ensure_free_storage

   !       if (CS_table(n)%description == "O2+O2->O2^-+O2" .or. &
   !            & CS_table(n)%description == "O2 -> products" .or. &
   !            & CS_table(n)%description == "O2 -> O2^-") then

   !          ! 3-body attachment, should be multiplied by gas density / cm^3
   !          CS_table(n)%en_cs(2,:) = 1.0D-6 * GAS_numberDensity * CS_table(n)%en_cs(2,:)
   !       CS_table(n)%en_cs(2,:) = CS_table(n)%en_cs(2,:) * GAS_fractions(CS_table(n)%gas_index) * GAS_numberDensity

   !> Write a list of all the collision processes in the simulation to a file 'filename',
   !! together with their type (elastic, excitation, etc) and a short description.
   subroutine CS_write_all(filename)
      character(LEN=*), intent(in) :: filename
      integer                      :: ics, n, io_state, my_unit

      my_unit = 333
      open(my_unit, FILE = filename, ACTION = "WRITE", ERR = 999, IOSTAT = io_state)

      write(my_unit, *, ERR = 999) "# A list of all the cross sections that have been read in by the"
      write(my_unit, *, ERR = 999) "# m_cross_sec.f90 module. You can use this file as input again."
      write(my_unit, *, ERR = 999) ""

      do ics = 1, CS_num_colls_found
         select case (CS_table(ics)%col_type)
         case (CS_elastic_t)
            write(my_unit, *, ERR = 999) "ELASTIC"
            write(my_unit, *, ERR = 999) trim(CS_table(ics)%description)
            write(my_unit, *, ERR = 999) CS_table(ics)%spec_value, " / mass ratio"
         case (CS_excite_t)
            write(my_unit, *, ERR = 999) "EXCITATION"
            write(my_unit, *, ERR = 999) trim(CS_table(ics)%description)
            write(my_unit, *, ERR = 999) CS_table(ics)%spec_value, " / threshold energy"
         case (CS_attach_t)
            write(my_unit, *, ERR = 999) "ATTACHMENT"
            write(my_unit, *, ERR = 999) trim(CS_table(ics)%description)
         case (CS_ionize_t)
            write(my_unit, *, ERR = 999) "IONIZATION"
            write(my_unit, *, ERR = 999) trim(CS_table(ics)%description)
            write(my_unit, *, ERR = 999) CS_table(ics)%spec_value, " / threshold energy"
         end select

         write(my_unit, *, ERR = 999) trim(CS_table(ics)%comment)
         write(my_unit, *, ERR = 999) "------------------------"
         do n = 1, CS_table(ics)%n_rows
            write(my_unit, *, ERR = 999) CS_table(ics)%en_cs(:, n)
         end do
         write(my_unit, *, ERR = 999) "------------------------"
         write(my_unit, *, ERR = 999) ""
      end do

      close(my_unit, ERR = 999, IOSTAT = io_state)
      return

999   continue ! If there was an error, the routine will end here
      print *, "CS_write_all error, io_state = ", io_state, " while writing to ", filename
      stop

   end subroutine CS_write_all

   subroutine CS_write_summary(filename)
      character(LEN=*), intent(in) :: filename
      character(LEN=name_len)      :: col_name
      integer                      :: n, io_state, my_unit
      my_unit = 333

      open(my_unit, FILE = filename, ACTION = "WRITE", ERR = 999, IOSTAT = io_state)

      write(my_unit, ERR = 999, FMT = "(A)") "# List of collision processes"
      write(my_unit, ERR = 999, FMT = "(A)") "Index      Gasname            Coltype     Description"
      write(my_unit, ERR = 999, FMT = "(A)") "-----------------------------------------------------"

      do n = 1, CS_num_colls_found
         select case (CS_table(n)%col_type)
         case (CS_elastic_t)
            col_name = "Elastic"
         case (CS_excite_t)
            col_name = "Excitation"
         case (CS_attach_t)
            col_name = "Attachment"
         case (CS_ionize_t)
            col_name = "Ionization"
         end select

         write(my_unit, ERR = 999, FMT = "((I4),(A),(A12),(A),(A15),(A),(A25))") &
              n, "    ", trim(CS_table(n)%gas_name), "  ", trim(col_name), "     ", CS_table(n)%description
      end do

      close(my_unit, STATUS = "KEEP", ERR = 999, IOSTAT = io_state)
      return

999   continue ! If there was an error, the routine will end here
      print *, "CS_write_summary error, io_state = ", io_state, " while writing to ", filename
      stop

   end subroutine CS_write_summary

end module m_cross_sec
