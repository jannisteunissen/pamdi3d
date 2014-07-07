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

!> Module that contains routines to read in cross section data from textfiles.
module module_crossSections

   use module_constants
   use generalUtilities
   use module_config
   use MPI
   use module_gas

   implicit none
   private

   !> The type of cross section table
   type CStable
      double precision, pointer :: en_cs(:,:)
      integer              :: tSize
      integer              :: colType
      integer              :: gasIndex
      double precision     :: treshold_eV
      double precision     :: elecOverMoleculeMass
      double precision     :: OpalSplitEnergy
      character(LEN=40)    :: gasName
      character(LEN=100)   :: comment
   end type CStable

   type(CStable), dimension(:), allocatable :: crossSecTable

   double precision :: maxEnergyEv
   double precision :: crossSecScaleFactor

   integer, parameter :: attachType = 1, elasticType = 2, exciteType = 3, ionizeType = 4, nColTypes = 4

   integer :: nProcesses
   integer, parameter :: maxNProcessPerGas = 50
   integer, parameter :: tSizeMax = 200

   ! Public variables
   public :: CStable, crossSecTable, attachType, elasticType, exciteType, ionizeType, &
        nColTypes, nProcesses
   ! Methods
   public :: CS_initialize
   public :: CS_readFromFile
   public :: CS_writeProcesses
   public :: CS_shareMPI


contains

   !> Set up the reading of textfiles for cross section data.
   !!
   !! Get the filenames, number of gases etc. from the configuration,
   !! and allocate space for the cross sections table.
   subroutine CS_initialize()

      nProcesses = 0

      allocate( crossSecTable(GAS_nGases * maxNProcessPerGas) )

      maxEnergyEv          = CFG_varDble("part_maxEnergyEv")
      crossSecScaleFactor  = CFG_varDble("sim_crossSecScaleFactor")

   end subroutine CS_initialize

   subroutine CS_shareMPI(myrank, root)
      integer, intent(IN) :: myrank, root

      integer :: n, ierr, tSize

      ! Share the number of processes that have been read in
      call MPI_BCAST(nProcesses, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)

      do n = 1, nProcesses
         call MPI_BCAST(crossSecTable(n)%tSize, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
         call MPI_BCAST(crossSecTable(n)%colType, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
         call MPI_BCAST(crossSecTable(n)%gasIndex, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
         call MPI_BCAST(crossSecTable(n)%treshold_eV, 1, &
              & MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
         call MPI_BCAST(crossSecTable(n)%elecOverMoleculeMass, 1, &
              & MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
         call MPI_BCAST(crossSecTable(n)%OpalSplitEnergy, 1, &
              & MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
         call MPI_BCAST(crossSecTable(n)%gasName, 40, MPI_CHARACTER, root, MPI_COMM_WORLD, ierr)
         call MPI_BCAST(crossSecTable(n)%comment, 100, MPI_CHARACTER, root, MPI_COMM_WORLD, ierr)

         tSize = crossSecTable(n)%tSize
         if (myrank /= root) then
            allocate( crossSecTable(n)%en_cs(2, tSize) )
         end if
         call MPI_BCAST(crossSecTable(n)%en_cs, 2*tSize, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
      end do

   end subroutine CS_shareMPI

   !> Read in the cross section data for the gases in 'gasnames' from their respective files,
   !! and then adjust the cross sections to fit our simulation
   subroutine  CS_readFromFile()
      integer :: n

      do n = 1, GAS_nGases
         call readCrossSecFromFile( trim(GAS_names(n)), "input/" // trim(GAS_files(n)), n)
      end do

      call adjustCSTable()

   end subroutine  CS_readFromFile

   subroutine readCrossSecFromFile(gasName, filename, gasIndex)
      ! Search 'filename' for cross section data concerning 'gasName', and if found
      ! add the cross sections to oldCSTable
      character(LEN=*), intent(IN)  :: gasName, filename
      integer, intent(IN)           :: gasIndex

      ! Temporary variables
      integer, parameter         :: lineLength = 100
      integer                    :: gasNameLength
      integer                    :: cIx
      integer                    :: ioState, nL
      integer                    :: tSize, collisionType
      character(LEN=40)          :: lineFMT
      character(LEN=lineLength)  :: line, prevLine
      double precision, dimension(2,tSizeMax) :: tempArray

      ! Set the number of lines to 0
      nL = 0
      gasNameLength = len(gasName)

      ! Set the line format to read, only depends on lineLength currently
      write(lineFMT, FMT = "(I6)") lineLength
      lineFMT = "(A" // trim(adjustl(lineFMT)) // ")"

      ! Open 'filename' (with error checking)
      open(UNIT = 1, FILE = filename, STATUS = "OLD", ACTION = "READ", ERR = 999, IOSTAT = ioState)

      ! Look for collision processes with the correct gas name in the file,
      ! which should look for example like:

      !     ATTACHMENT                    [description of the type of process, always in CAPS]
      !     H2O -> H2O^-                  [the gas name possibly followed by the result of the process]
      !     COMMENT: total attachment     [possibly comments]
      !     UPDATED: 2010-06-24 15:04:36
      !     ------------------            [at least 5 dashes]
      !     xxx [eV]  xxx [m2]            [cross section data in two column format]
      !     ...       ...
      !     xxx [eV]  xxx [m2]
      !     ------------------

      ! So if we find the gas name the previous line holds the type of collision, then
      ! there is possibly some extra information and between the dashes the actual cross
      ! sections are found. The first column holds the energy in eV and the second the
      ! cross section at that energy in m^2.

      ! The outer DO loop, running until the end of the file is reached
      do
         ! Search for 'gasName' in the file
         line = ' '
         do
            prevLine = line
            read(UNIT = 1, FMT = lineFMT, ERR = 999, end = 666) line; nL = nL+1
            if ( line(1:gasNameLength) == gasName ) exit
         end do

         ! Check prevLine for the type of collision
         select case ( prevLine(1:15) )
         case ("ATTACHMENT")
            collisionType = attachType
         case ("ELASTIC", "MOMENTUM")
            collisionType = elasticType
         case ("EXCITATION")
            collisionType = exciteType
         case ("IONIZATION")
            collisionType = ionizeType
         case ("COMMENT")
            cycle
         case DEFAULT
            print *, "readCrossSec warning: Can not determine process type for "
            print *, " ", gasName, " in ", filename, " at line ", nL
            cycle
         end select

         ! Update the number of processes and set the gas name and collision type
         nProcesses = nProcesses + 1
         cIx = nProcesses

         crossSecTable(cIx)%gasName = gasName
         crossSecTable(cIx)%colType = collisionType
         crossSecTable(cIx)%gasIndex = gasIndex

         ! Add the reaction description to the table
         crossSecTable(cIx)%comment = line

         ! Set the energy constant for Opal's formula
         select case (gasName)
         case ("O2")
            crossSecTable(cIx)%OpalSplitEnergy = 17.4D0 ! (eV)
            crossSecTable(cIx)%elecOverMoleculeMass = electronMass / (32 * atomicMass)
         case ("N2")
            crossSecTable(cIx)%OpalSplitEnergy = 13.0D0 ! (eV)
            crossSecTable(cIx)%elecOverMoleculeMass = electronMass / (28 * atomicMass)
         case ("Ar")
            crossSecTable(cIx)%OpalSplitEnergy = 10.0D0 ! (eV)
            crossSecTable(cIx)%elecOverMoleculeMass = electronMass / (40 * atomicMass)
         case DEFAULT
            print *, "readCrossSec warning: I don't have Opal splitting energy"
            print *, "and molecule mass for ", gasName, ", using default values (inaccurate)"

            crossSecTable(cIx)%OpalSplitEnergy = 10.0D0 ! (eV)
            crossSecTable(cIx)%elecOverMoleculeMass = electronMass / (20 * atomicMass)
         end select

         ! In case of excitations and ionizations, store the threshold energy, which
         ! is on the next line.
         if (collisionType == exciteType .or. collisionType == ionizeType) then
            read(UNIT = 1, FMT = lineFMT, ERR = 999, end = 555) line; nL = nL+1
            read(line, FMT = *, ERR = 999, end = 555) crossSecTable(cIx)%treshold_eV
         else
            crossSecTable(cIx)%treshold_eV = 0.0D0
         end if

         ! Now we can check whether there is a ZDPLASKIN statement and a reaction description,
         ! while scanning lines until dashes are found, which indicate the start of the cross section data
         do
            read(UNIT = 1, FMT = lineFMT, ERR = 999, end = 555) line; nL = nL+1
            if ( line(1:9) == "ZDPLASKIN" ) then
               crossSecTable(cIx)%comment = adjustl( line(12:) )
            end if
            if ( line(1:5) == "-----" ) exit
         end do

         ! Read the cross section data into a temporary array
         tSize = 0
         do
            read(UNIT = 1, FMT = lineFMT, ERR = 999, end = 555) line; nL = nL+1
            if ( line(1:5) == "-----" ) then    ! Dashes mark the end of the data
               exit
            else if (tSize < tSizeMax) then
               tSize = tSize + 1
               read(line, FMT = *, ERR = 999, end = 555) tempArray(:, tSize)
            end if
         end do

         ! Store the data in the actual table
         allocate( crossSecTable(cIx)%en_cs(2, tSize) )
         crossSecTable(cIx)%tSize = tSize
         crossSecTable(cIx)%en_cs = tempArray(:, 1:tSize)

         ! Check whether the tables that have been read in go up to high enough energies for our
         ! simulation. They can also have 0.0 as their highest listed cross section, in which
         ! case we assume the cross section is 0.0 for all higher energies

         if ( crossSecTable(cIx)%en_cs(1, tSize) < maxEnergyEv .and. &
              & crossSecTable(cIx)%en_cs(2, tSize) > 0.0D0 ) then
            print *, "readCrossSecFromFile warning: cross section data at line ", nL
            print *, "does not go up to high enough energies"
         end if

      end do

555   continue ! Routine ends here if the end of "filename" is reached erroneously
      print *, "readCrossSec error, reached end of file while searching:"
      print *, "ioState = ", ioState, " while reading from [", filename, "] at line ", nL
      close(UNIT = 1, ERR = 999, IOSTAT = ioState)
      return


666   continue ! Routine ends here if the end of "filename" is reached correctly
      close(UNIT = 1, ERR = 999, IOSTAT = ioState)
      return


999   continue ! If there was an error, the routine will end here
      print *, "readCrossSec error at line", nL
      print *, "ioState = ", ioState, " while searching [", gasName, "] in [", filename, "]"
      stop

   end subroutine readCrossSecFromFile

   subroutine adjustCSTable()
      ! Adjusts the cross sections that have been read in to a form that is suitable for
      ! the Monte Carlo simulations
      integer :: n

      ! Transform elastic momentum transfer cross sections to total elastic cross sections.
      ! Usually the total momentum transfer cross sections are given, which are defined as
      ! the expectation value of [1 - Cos(Xhi)] over the differential cross section,
      ! with Xhi the scattering angle. If one makes assumptions about the differential cross
      ! section, one can deduce the total cross section from the momentum transfer cross section,
      ! see for example: 'Electron anisotropic scattering in gases: A formula for Monte Carlo
      ! simulations', by Okhrimovskyy et al. - DOI 10.1103/PhysRevE.65.037402 and
      ! 'Self-consistent model of a direct-current glow discharge: Treatment of fast electrons'
      ! by Surendra et al. - DOI 10.1103/PhysRevA.41.1112

      ! Check for 3-body cross sections that have to be multiplied by the gas density
      do n = 1, nProcesses
         if (crossSecTable(n)%comment == "O2+O2->O2^-+O2") then

            ! 3-body attachment, should be multiplied by gas density / cm^3
            crossSecTable(n)%en_cs(2,:) = 1.0D-6 * GAS_numberDensity * crossSecTable(n)%en_cs(2,:)

         end if
      end do

      ! Multiply the total cross section (number of collision of one scattering center per unit time
      ! per unit incoming flux) by the number of scattering centers (gas number density * gas fraction),
      ! and possibly by a scale factor different from 1

      do n = 1, nProcesses
         crossSecTable(n)%en_cs(2,:) = crossSecTable(n)%en_cs(2,:) * crossSecScaleFactor * &
              & GAS_fractions(crossSecTable(n)%gasIndex) * GAS_numberDensity
      end do

   end subroutine adjustCSTable

   !> Write a list of all the collision processes in the simulation to a file 'filename',
   !! together with their type (elastic, excitation, etc) and a short description.
   subroutine CS_writeProcesses(filename)
      character(LEN=*) :: filename
      character(LEN=15) :: colName
      integer :: n, ioState

      open(UNIT = 1, FILE = filename, ACTION = "WRITE", ERR = 999, IOSTAT = ioState)

      write(UNIT = 1, ERR = 999, FMT = "(A)") "# List of collision processes"
      write(UNIT = 1, ERR = 999, FMT = "(A)") "Index Gasname   Coltype     Description"
      write(UNIT = 1, ERR = 999, FMT = "(A)") "---------------------------------------"

      write(UNIT = 1, ERR = 999, FMT = "((I4),(A),(A6),(A),(A12),(A),(A25))") &
           0, "  ", "None  ", "  ", "None      ", "  ", "Null collision               "
      do n = 1, nProcesses
         select case (crossSecTable(n)%colType)
         case (elasticType)
            colName = "Elastic"
         case (exciteType)
            colName = "Excitation"
         case (attachType)
            colName = "Attachment"
         case (ionizeType)
            colName = "Ionization"
         end select

         write(UNIT = 1, ERR = 999, FMT = "((I4),(A),(A6),(A),(A17),(A),(A25))") &
              n, "  ", crossSecTable(n)%gasName, "  ",colName, "  ", crossSecTable(n)%comment
      end do

      close(UNIT = 1, STATUS = "KEEP", ERR = 999, IOSTAT = ioState)

      return

999   continue ! If there was an error, the routine will end here
      print *, "CS_writeProcesses error:"
      print *, "ioState = ", ioState, " while writing to ", filename

   end subroutine CS_writeProcesses

   !> Anbang: This routine is for verifying the code when magnetic field is included,
   ! We use Reid ramp model, while the elastic and inelastic cross section are defined specificly
   !    subroutine CS_reid_ramp_model()
   !
   !       Integer :: n
   !
   !       DO n = 1, nProcesses
   !          IF (crossSecTable(n)%colType == elasticType) THEN
   !             crossSecTable(n)%en_cs(2,:) = 6.D-20    !elastic corss section
   !             crossSecTable(n)%en_cs(2,:) = crossSecTable(n)%en_cs(2,:) * crossSecScaleFactor * GAS_numberDensity * &
   !             & GAS_fractions(crossSecTable(n)%gasIndex)
   !          ELSE
   !             crossSecTable(n)%en_cs(2,:) = 10 * (crossSecTable(n)%en_cs(1,:)-0.2d0)*1.d-20   !inelastic cross section
   !             crossSecTable(n)%en_cs(2,:) = crossSecTable(n)%en_cs(2,:) * crossSecScaleFactor * GAS_numberDensity * &
   !             & GAS_fractions(crossSecTable(n)%gasIndex)
   !          END IF
   !       END DO
   !    end subroutine CS_reid_ramp_model

end module module_crossSections
