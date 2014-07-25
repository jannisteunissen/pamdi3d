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

  type CS_coll_t
     integer  :: type       = -1
     real(dp) :: part_mass  = 0
     real(dp) :: rel_mass   = 0
     real(dp) :: en_loss    = 0
  end type CS_coll_t

  !> The type of cross section table
  type CS_t
     type(CS_coll_t) :: coll
     real(dp), allocatable   :: en_cs(:,:)    ! Stores the energy vs cross sec table
     integer                 :: n_rows        ! Number of rows in the table
     real(dp)                :: min_energy    ! Minimum energy for the collision
     real(dp)                :: max_energy    ! Maximum energy for the collision
     character(LEN=tiny_len) :: gas_name      ! Name of the colliding neutral molecule
     character(LEN=line_len) :: description   ! Description of the collision
     character(LEN=line_len) :: comment       ! Additional comments
  end type CS_t

  integer, parameter :: CS_attach_t = 1, &
       CS_elastic_t = 2, &
       CS_excite_t = 3, &
       CS_ionize_t = 4, &
       CS_num_types = 4
  integer, parameter :: max_num_cols_per_gas = 50
  integer, parameter :: max_num_rows = 400

  ! Public variables
  public :: CS_t
  public :: CS_coll_t
  public :: CS_attach_t, CS_elastic_t, CS_excite_t, CS_ionize_t

  ! Methods
  public :: CS_add_from_file
  public :: CS_write_summary

contains

  ! Search 'filename' for cross section data concerning 'gas_name'
  subroutine CS_add_from_file(filename, gas_name, x_normalization, &
       y_normalization, req_energy, cross_secs)
    use m_units_constants
    character(len=*), intent(IN) :: gas_name, filename
    real(dp), intent(in)         :: x_normalization, y_normalization
    real(dp), intent(in)         :: req_energy
    type(CS_t), intent(inout), allocatable :: cross_secs(:)
    type(CS_t), allocatable :: cs_cpy(:)
    type(CS_t) :: cs_buf(max_num_cols_per_gas)
    integer                      :: n, cIx, nL, n_rows, col_type
    integer                      :: my_unit, io_state, len_gas_name
    character(LEN=name_len)      :: lineFMT
    character(LEN=line_len)      :: line, prev_line
    real(dp)                     :: tempArray(2, max_num_rows)
    real(dp)                     :: x_scaling, y_scaling, tmp_value

    my_unit      = 333
    nL           = 0 ! Set the number of lines to 0
    cIx = 0
    len_gas_name = len(trim(gas_name))

    ! Set the line format to read, only depends on line_len currently
    write(lineFMT, FMT = "(A,I0,A)") "(A", line_len, ")"

    ! Open 'filename' (with error checking)
    open(my_unit, FILE = trim(filename), ACTION = "READ", &
         ERR = 999, IOSTAT = io_state)

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
               trim(gas_name), ", this should not be used in particle simulations."
       case ("EXCITATION")
          col_type = CS_excite_t
       case ("IONIZATION")
          col_type = CS_ionize_t
       case ("COMMENT")
          cycle
       case DEFAULT
          print *, "CS_read_file warning: ignoring unknown process type for ", &
               trim(gas_name), " in ", filename, " at line ", nL
          cycle
       end select

       ! Update the number of processes and set the gas name and collision type
       cIx = cIx + 1

       ! Add the reaction description to the table
       cs_buf(cIx)%description = adjustl(line)

       ! For all collisions except attachment, there is a value on the next line
       if (col_type /= CS_attach_t) then
          read(my_unit, FMT = lineFMT, ERR = 999, end = 555) line; nL = nL+1
          read(line, FMT = *, ERR = 999, end = 555) tmp_value
       else
          tmp_value = 0
       end if

       cs_buf(cIx)%coll%type = col_type
       cs_buf(cIx)%coll%part_mass = UC_elec_mass

       select case(col_type)
       case (CS_elastic_t)
          cs_buf(cIx)%coll%rel_mass = tmp_value ! Store relative mass e/M
       case (CS_excite_t)
          cs_buf(cIx)%coll%en_loss  = tmp_value ! Energy loss
       case (CS_ionize_t)
          cs_buf(cIx)%coll%en_loss  = tmp_value ! Energy loss
       end select

       cs_buf(cIx)%gas_name = gas_name
       cs_buf(cIx)%comment = "COMMENT: (empty)"
       x_scaling = 1.0_dp
       y_scaling = 1.0_dp

       ! Now we can check whether there is a ZDPLASKIN statement and a
       ! reaction description, while scanning lines until dashes are
       ! found, which indicate the start of the cross section data
       do
          read(my_unit, FMT = lineFMT, ERR = 999, end = 555) line; nL = nL+1
          line = adjustl(line)
          if ( line(1:9) == "ZDPLASKIN" ) then
             cs_buf(cIx)%description = trim(gas_name) // " [" // trim(adjustl(line(11:))) // "]"
          else if ( line(1:7) == "COMMENT") then
             cs_buf(cIx)%comment = line
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
       allocate(cs_buf(cIx)%en_cs(2, n_rows))
       cs_buf(cIx)%n_rows = n_rows
       cs_buf(cIx)%en_cs(1,:) = tempArray(1, 1:n_rows) * x_normalization * x_scaling
       cs_buf(cIx)%en_cs(2,:) = tempArray(2, 1:n_rows) * y_normalization * y_scaling
       cs_buf(cIx)%max_energy = cs_buf(cIx)%en_cs(1, n_rows)

       ! Locate minimum energy (first value followed by non-zero cross sec)
       cs_buf(cIx)%min_energy = 0.0_dp
       do n = 1, n_rows-1
          if (cs_buf(cIx)%en_cs(2, n+1) > 0.0_dp) then
             cs_buf(cIx)%min_energy = cs_buf(cIx)%en_cs(1, n)
             exit
          end if
       end do

       ! Locate maximum energy (last value preceded by non-zero)
       cs_buf(cIx)%max_energy = 0.0_dp
       do n = n_rows, 2, -1
          if (cs_buf(cIx)%en_cs(2, n-1) > 0.0_dp) then
             cs_buf(cIx)%max_energy = cs_buf(cIx)%en_cs(1, n)
             exit
          end if
       end do

       ! Check whether the tables that have been read in go up to high enough energies for our
       ! simulation. They can also have 0.0 as their highest listed cross section, in which
       ! case we assume the cross section is 0.0 for all higher energies

       if ( cs_buf(cIx)%en_cs(1, n_rows) < req_energy .and. &
            & cs_buf(cIx)%en_cs(2, n_rows) > 0.0D0 ) then
          print *, "CS_read_file error: cross section data at line ", nL, &
               " does not go up to high enough x-values (energy)"
          stop
       end if

    end do

555 continue ! Routine ends here if the end of "filename" is reached erroneously
    close(my_unit, ERR = 999, IOSTAT = io_state)
    print *, "CS_read_file error, reached end of file while searching. ", &
         "io_state = ", io_state, " while reading from [", filename, "] at line ", nL
    stop
    return

666 continue ! Routine ends here if the end of "filename" is reached correctly
    close(my_unit, ERR = 999, IOSTAT = io_state)

    ! Set the output data
    if (allocated(cross_secs)) then
       ! Resize array
       n = size(cross_secs)
       allocate(cs_cpy(n))
       cs_cpy = cross_secs
       deallocate(cross_secs)
       allocate(cross_secs(n + cIx))
       cross_secs(1:n) = cs_cpy
       cross_secs(n+1:) = cs_buf(1:cIx)
    else
       allocate(cross_secs(cIx))
       cross_secs(:) = cs_buf(1:cIx)
    end if

    return

999 continue ! If there was an error, the routine will end here
    print *, "CS_read_file error at line ", nL, " io_state = ", io_state, &
         " while searching [", gas_name, "] in [", filename, "]"
    stop

  end subroutine CS_add_from_file

  subroutine CS_write_summary(cross_secs, filename)
    type(CS_t), intent(in)       :: cross_secs(:)
    character(LEN=*), intent(in) :: filename
    character(LEN=name_len)      :: col_name
    integer                      :: n, io_state, my_unit
    my_unit = 333

    open(my_unit, FILE = filename, ACTION = "WRITE", &
         ERR = 999, IOSTAT = io_state)

    write(my_unit, ERR = 999, FMT = "(A)") "# List of collision processes"
    write(my_unit, ERR = 999, FMT = "(A)") &
         "Index      Gasname            Coltype     Description"
    write(my_unit, ERR = 999, FMT = "(A)") &
         "-----------------------------------------------------"

    do n = 1, size(cross_secs)
       select case (cross_secs(n)%coll%type)
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
            n, "    ", trim(cross_secs(n)%gas_name), "  ", trim(col_name), &
            "     ", cross_secs(n)%description
    end do

    close(my_unit, STATUS = "KEEP", ERR = 999, IOSTAT = io_state)
    return

999 continue ! If there was an error, the routine will end here
    print *, "CS_write_summary error, io_state = ", io_state, &
         " while writing to ", filename
    stop

  end subroutine CS_write_summary

end module m_cross_sec
