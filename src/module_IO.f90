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

!> Very basic version of a module that provides basic file input and output
module module_IO
   use module_globals
   use module_config
   use module_particle
   use m_efield_amr
   use m_write_silo
   use MPI

   implicit none
   private

   !> Maximum length of a file/group/dataset name
   integer, parameter :: MAXLENGTH = 80

   public :: writeOutput
   public :: writeDoubleToText

contains
   !> Writes the program output to HDF files
   subroutine writeOutput(timeStep, simulationTime, simname, myrank, root)
      integer, intent(IN)           :: timeStep, myrank, root
      double precision, intent(IN)  :: simulationTime
      character(LEN=*), intent(IN)  :: simname

      character(LEN=MAXLENGTH) :: filename

      ! Only root writes output currently, so the other tasks can return after
      ! they have shared their data

      ! Create the filename
      write(filename, FMT = "(I6.6)") timeStep
      filename = "output/" // trim(simname) // "_" // trim(filename) // ".silo"
      call PM_particlesToDensity()
      call E_write_grids(filename, myrank, root, timeStep, simulationTime)
   end subroutine writeOutput

   !> Write a
   subroutine writeDoubleToText(array1D, array2D, header, filename)
      double precision, dimension(:), optional, intent(IN)   :: array1D
      double precision, dimension(:,:), optional, intent(IN) :: array2D
      character(LEN=*), intent(IN), optional                 :: header
      character(LEN=*), intent(IN)                           :: filename

      integer :: i, ios

      open(UNIT=1, FILE=filename, IOSTAT=ios)

      if (ios /= 0) then
         print *, "writeDoubleToText error, iostat = ", ios, " while writing ", filename
         return
      end if

      if (present(header)) then
         write(1, FMT = *) header
      end if

      if (present(array1D)) then
         write(1, FMT = '(E16.4e4)') array1D
      else if (present(array2D)) then
         do i = 1, size(array2D, 1)
            write(1, *) array2D(i, :)
         end do
      end if

      close(UNIT=1)

   end subroutine writeDoubleToText

end module module_IO
