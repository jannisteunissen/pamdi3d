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

module m_write_silo

   implicit none
   private

   include 'silo.inc'

   integer, parameter :: dp      = kind(0.0d0)
   integer, parameter :: lineLen = 200
   integer, parameter :: DB_TYPE = DB_PDB

   public :: create_file
   public :: add_grid
   public :: add_data
   public :: set_multimesh

contains

   subroutine create_file(filename)
      character(len=*), intent(in)     :: filename
      integer, intent(inout), optional :: fileIx
      integer                          :: ierr, dbix
      character(len=lineLen)           :: fileinfo

      fileinfo = "A silo file"
      ierr = dbcreate(trim(filename), len_trim(filename), DB_CLOBBER, DB_LOCAL, &
           fileinfo, len_trim(fileinfo), DB_TYPE, dbix)
      if (ierr /= 0) print *, "Error creating file", trim(filename)
      call close_file(dbix)
   end subroutine create_file

   function open_file(filename) result(dbix)
      character(len=*), intent(in) :: filename
      integer :: dbix, ierr

      ierr = dbopen(trim(filename), len_trim(filename), DB_TYPE, DB_APPEND, dbix)
      if (ierr /= 0) print *, "Error opening file", trim(filename)
   end function open_file

   subroutine close_file(dbix)
      integer, intent(in) :: dbix
      integer :: ierr

      ierr = dbclose(dbix)
      if (ierr /= 0) print *, "Error closing file with index", dbix
   end subroutine close_file

   subroutine set_dir_create(dbix, dirname)
      character(len=*), intent(in) :: dirname
      integer, intent(in) :: dbix
      integer :: ierr

      ierr = dbmkdir(fileTable(ix)%dbix, trim(dirname), len_trim(dirname), iostat)
      if (ierr /= 0) then
         ierr = dbmkdir(dbix, trim(dirname), len_trim(dirname), iostat)
         if (ierr /= 0) print *, "Error creating directory", dirname
         ierr = dbsetdir(dbix, trim(dirname), len_trim(dirname))
         if (ierr /= 0) print *, "Could not go to directory", dirname
      end if
   end subroutine set_dir_create

   subroutine add_grid(filename, dirname, gridname, n_dim, N_r, r_min, dr)
      character(len=*), intent(in) :: filename, dirname, gridname
      integer, intent(in)          :: n_dim, N_r(:)
      real(dp), intent(in)         :: r_min(:), dr(:)

      real(dp), allocatable        :: x_coords(:), y_coords(:), z_coords(:)
      integer                      :: i, ierr, iostat, dboptix, dbix

      interface
         function dbputqm(dbid, name, lname, xname, lxname, yname, &
              lyname, zname, lzname, x, y, z, dims, ndims, &
              datatype, coordtype, optlist_id, status)
            use, intrinsic :: iso_c_binding
            integer(c_int) :: dbid, lname, lxname, lyname, lzname, dims(*), ndims
            integer(c_int) :: datatype, coordtype, optlist_id, status, dbputqm
            real(c_double) :: x(*), y(*), z(*)
            character(kind=c_char) :: name(*), xname(*), yname(*), zname(*)
         end function dbputqm
      end interface

      if (n_dim < 1 .or. n_dim > 3) then
         print *, "Cannot add grid for which n_dim < 1 or n_dim > 3"
         return
      end if

      allocate(x_coords(N_r(1)))
      do i = 1, N_r(1)
         x_coords(i) = r_min(1) + (i-1) * dr(1)
      end do

      if (n_dim > 1) then
         allocate(y_coords(N_r(2)))
         do i = 1, N_r(2)
            y_coords(i) = r_min(2) + (i-1) * dr(2)
         end do
      end if

      if (n_dim > 2) then
         allocate(z_coords(N_r(2)))
         do i = 1, N_r(3)
            z_coords(i) = r_min(3) + (i-1) * dr(3)
         end do
      end if

      ! Make option list
      ierr = dbmkoptlist(20, dboptix)

      ! Set string options
      ierr = dbaddcopt(dboptix, DBOPT_XLABEL, 'X'//char(0), 2)
      ierr = dbaddcopt(dboptix, DBOPT_YLABEL, 'Y'//char(0), 2)
      ierr = dbaddcopt(dboptix, DBOPT_ZLABEL, 'Z'//char(0), 2)
      ierr = dbaddcopt(dboptix, DBOPT_XUNITS, 'm'//char(0), 2)
      ierr = dbaddcopt(dboptix, DBOPT_YUNITS, 'm'//char(0), 2)
      ierr = dbaddcopt(dboptix, DBOPT_ZUNITS, 'm'//char(0), 2)

      ! Set integer options
      ierr = dbaddiopt(dboptix, DBOPT_NSPACE, n_dim)
      ierr = dbaddiopt(dboptix, DBOPT_HIDE_FROM_GUI, 1)
      !       ierr = dbaddiopt(dboptix, DBOPT_MAJORORDER, 1)

      dbix = open_file(filename)
      call set_dir_create(dirname)

      ! Write the grid structure
      ierr = dbputqm(dbix, trim(gridname), len_trim(gridname), 'X', 1, 'Y', 1, 'Z', 1, x_coords, y_coords, z_coords, &
           & N_r, n_dim, DB_DOUBLE, DB_COLLINEAR, dboptix, iostat)

      call close_file(dbix)
      ierr = dbfreeoptlist(dboptix)
   end subroutine add_grid

   subroutine add_data(filename, dirname, dataname, gridname, data_packed, data_shape, data_unit)
      character(len=*), intent(in) :: filename, dirname, gridname, dataname, data_unit
      real(dp), intent(in)         :: data_packed(:)
      integer, intent(in)          :: data_shape(:)

      integer                      :: dboptix, ierr, iostat, dbix
      real(dp)                     :: dummy(1)

      interface
         function dbputqv1(dbid, name, lname, meshname, lmeshname, &
              var, dims, ndims, mixvar, mixlen, datatype, &
              centering, optlist_id, status)
            use, intrinsic :: iso_c_binding
            integer(c_int) :: dbid, lname, lmeshname, dims(*), ndims, mixlen
            integer(c_int) :: centering, optlist_id, status, datatype, dbputqv1
            real(c_double) :: var(*), mixvar(*)
            character(kind=c_char) :: name(*), meshname(*)
         end function dbputqv1
      end interface

      if (size(data_packed) /= product(data_shape)) then
         print *, "Error: data_packed does not correspond to data_shape"
         return
      end if

      if (size(data_shape) < 1 .or. size(data_shape) > 3) then
         print *, "Error: size(data_shape) < 1 or size(data_shape) > 3"
         return
      end if

      ierr = dbmkoptlist(10, dboptix)
      ierr = dbaddcopt(dboptix, DBOPT_UNITS, trim(data_unit)//char(0), len_trim(data_unit))
      ierr = dbaddiopt(dboptix, DBOPT_HIDE_FROM_GUI, 1)
      !       ierr = dbaddiopt(dboptix, DBOPT_MAJORORDER, 1)

      ! Get ix corresponding to filename
      dbix = open_file(filename)
      call set_dir_create(dirname, dbix)

      ! Write the data to the grid
      ierr = dbputqv1(dbix, trim(dataname), len_trim(dataname), trim(gridname), len_trim(gridname), &
           & packed_data, data_shape, size(data_shape), dummy, 0, DB_DOUBLE, DB_NODECENT, dboptix, iostat)

      ! Close the file and remove option list
      call close_file(dbix)
      ierr = dbfreeoptlist(dboptix)
   end subroutine add_data

   subroutine set_multimesh(filename, mmdir, mmname, gridnames, datanames, n_cycle, time)
      character(len=*), intent(in) :: filename, mmdir, mmname, gridnames(:), datanames(:)
      integer, intent(in), optional :: n_cycle
      real(dp), intent(in), optional :: time

      integer :: int_dummy(1), ierr, dbix, dboptix, iostat, old_str_len
      integer, allocatable :: mtypes(:), name_lengths(:)
      character, allocatable :: mnames(:), dnames(:)

      interface
         function dbputmmesh(dbid, name, lname, nmesh, meshnames, lmeshnames, &
              meshtypes, optlist_id, status)
            use, intrinsic :: iso_c_binding
            integer(c_int) :: dbputmmesh, lmeshnames(*)
            integer(c_int) :: dbid, lname, nmesh, meshtypes(*), optlist_id, status
            character(kind=c_char) :: name(*), meshnames(*)
         end function dbputmmesh

         function dbputmvar(dbid, name, lname, nlevels, meshnames, lmnames, meshtypes, optlist_id, status)
            use, intrinsic :: iso_c_binding
            integer(c_int) :: dbputmvar, lmnames(*)
            integer(c_int) :: dbid, lname, nlevels, meshtypes(*), optlist_id, status
            character(kind=c_char) :: name(*), meshnames(*)
         end function dbputmvar

         function dbmkmrgtree(mesh_type, info_bits, max_children, optlist_id, tree_id)
            use, intrinsic :: iso_c_binding
            integer(c_int) :: mesh_type, info_bits, max_children, optlist_id, tree_id, dbmkmrgtree
         end function dbmkmrgtree

         function dbwrite(dbid, varname, lvarname, var, dims, ndims, datatype)
            use, intrinsic :: iso_c_binding
            integer(c_int) :: dbid, lvarname, var(*), dims(*), ndims, datatype, dbwrite
            character(kind=c_char) :: varname(*)
         end function dbwrite

      end interface

      n_grids = size(gridnames)
      if (size(datanames) /= n_grids) then
         print *, "Error, size(datanames) /= size(gridnames)"
         return
      end if

      if (n_grids < 1) then
         print *, "Error, too few grids (<1)"
         return
      end if

      ! Write grid multimesh
      name_len = len(gridnames(1))
      allocate(mnames(name_len * n_grids))
      allocate(name_lenghts(n_grids))
      allocate(mtypes(n_grids))

      mnames = pack(gridnames, .true.)
      old_str_len = dbset2dstrlen(name_len)
      m_types = DB_QUADMESH
      name_lengths = name_len

      dbix = open_file(filename)
      ierr = dbmkoptlist(10, dboptix)
      if (present(n_cycle)) ierr = dbaddiopt(dboptix, DBOPT_CYCLE, n_cycle)
      if (present(time)) ierr = dbaddiopt(dboptix, DBOPT_DTIME, time)

      ierr = dbputmmesh(dbix, trim(mmname), len_trim(mmname), n_grids, &
           mnames, name_lengths, mtypes, dboptix, iostat)

      ! Write the data multimesh(es)
      name_len = len(datanames(1))
      allocate(dnames(name_len * n_grids))

      dnames = pack(datanames, .true.)
      ierr = dbset2dstrlen(name_len)
      name_lengths = name_len

      ierr = dbaddcopt(dboptix, DBOPT_MMESH_NAME, trim(mmname)//char(0), len_trim(mmname))
      ierr = dbputmvar(dbix, trim(mvname), len_trim(mvname), nLevels, &
              dnames, name_lengths, mtypes(1:nLevels), dboptix, iostat)

      call close_file(dbix)
      ierr = dbfreeoptlist(dboptix)
      ierr = dbset2dstrlen(old_str_len)
   end subroutine set_multimesh

end module module_writeSilo
