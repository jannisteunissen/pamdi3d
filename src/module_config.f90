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

!> Module that allows working with a configuration file.
module module_config

  implicit none

  ! By default everything is private in this module, public interfaces are defined at the bottom
  private

  integer, parameter :: tinyLen = 20
  integer, parameter :: nameLen = 40
  integer, parameter :: lineLen = 1000
  integer, parameter :: maxVarSize = 50
  integer, parameter :: intType = 0, floatType = 1, charType = 2, logicType = 3

  type CFG_varType
     character(len=nameLen)           :: pName
     character(len=lineLen)           :: comment
     integer                          :: pType
     integer                          :: pSize
     logical                          :: varSize
     double precision, pointer        :: floatValue(:)
     integer, pointer                 :: intValue(:)
     character(len=nameLen), pointer  :: charValue(:)
     logical, pointer                 :: logicValue(:)
  end type CFG_varType

  logical  :: CFG_initialized    = .false.
  logical  :: CFG_sorted         = .false.

  integer  :: CFG_nVars
  integer  :: CFG_maxVars

  type(CFG_varType), allocatable :: CFG_varList(:)

  !> Overloaded interface to add config variables
  interface CFG_addVar
     module procedure CFG_addVarDouble, CFG_addVarDoubleArray, CFG_addVarInt, CFG_addVarIntArray, &
          CFG_addVarChar, CFG_addVarCharArray, CFG_addVarlogical, CFG_addVarlogicalArray
  end interface CFG_addVar

  !> Overloaded interface to get variables
  interface CFG_getVar
     module procedure CFG_getVarDouble, CFG_getVarDoubleArray, CFG_getVarInt, CFG_getVarIntArray, &
          CFG_getVarChar, CFG_getVarCharArray, CFG_getVarlogical, CFG_getVarlogicalArray
  end interface CFG_getVar

  public :: CFG_initialize, CFG_destruct
  public :: CFG_addVar, CFG_updateFromFile
  public :: CFG_varInt, CFG_varDble, CFG_varLogic
  public :: CFG_getSize, CFG_getType, CFG_getVar
  public :: CFG_sort, CFG_write
  public :: CFG_shareMPI

contains

  !> Prepare the configuration module for usage, with initially space
  !! for 'maxVars' variables
  subroutine CFG_initialize(maxVars)
    integer, intent(in), optional :: maxVars

    if (CFG_initialized) then
       print *, "CFG_initialize: Already CFG_initialized"
       return
    end if

    CFG_maxVars  = 100
    if (present(maxVars)) CFG_maxVars = maxVars

    allocate(CFG_varList(CFG_maxVars))

    CFG_nVars         = 0
    CFG_initialized   = .true.

  end subroutine CFG_initialize

  !> Destroy the configuration and deallocate the space used
  subroutine CFG_destruct

    integer :: i

    if (.not. CFG_initialized) then
       print *, "CFG_destruct: Not CFG_initialized yet"
       return
    end if

    do i = 1, CFG_nVars

       select case (CFG_varList(i)%pType)
       case (floatType)
          deallocate(CFG_varList(i)%floatValue)
       case (intType)
          deallocate(CFG_varList(i)%intValue)
       case (charType)
          deallocate(CFG_varList(i)%charValue)
       case (logicType)
          deallocate(CFG_varList(i)%logicValue)
       end select

    end do

    deallocate(CFG_varList)
    CFG_initialized = .false.

  end subroutine CFG_destruct

  !> MPI: broadcast the configuration at the root task to the other tasks. This routine is not very
  !! efficient, but that is usually not necessary, as it is executed only once.
  subroutine CFG_shareMPI(myrank, root)
    use mpi
    integer, intent(in)  :: myrank, root

    integer              :: n, i, ierr, pSize, newSize

    if (myrank == root) newSize = CFG_nVars

    ! Share the number of variables that have been read in
    call MPI_BCAST(newSize, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)

    if (myrank /= root) then
       call resizeVarList(newSize)
       CFG_nVars = newSize
    end if

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    do n = 1, CFG_nVars

       call MPI_BCAST(CFG_varList(n)%pType, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
       call MPI_BCAST(CFG_varList(n)%pSize, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
       call MPI_BCAST(CFG_varList(n)%pName, nameLen, MPI_CHARACTER, root, MPI_COMM_WORLD, ierr)
       call MPI_BCAST(CFG_varList(n)%comment, lineLen, MPI_CHARACTER, root, MPI_COMM_WORLD, ierr)
       call MPI_BCAST(CFG_varList(n)%varSize, 1, MPI_LOGICAL, root, MPI_COMM_WORLD, ierr)


       pSize = CFG_varList(n)%pSize

       if (myrank /= root) then
          select case (CFG_varList(n)%pType)
          case (intType)
             allocate( CFG_varList(n)%intValue(pSize) )
          case (floatType)
             allocate( CFG_varList(n)%floatValue(pSize) )
          case (logicType)
             allocate( CFG_varList(n)%logicValue(pSize) )
          case (charType)
             allocate( CFG_varList(n)%charValue(pSize) )
          end select
       end if

       select case (CFG_varList(n)%pType)
       case (intType)
          call MPI_BCAST(CFG_varList(n)%intValue, pSize, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
       case (floatType)
          call MPI_BCAST(CFG_varList(n)%floatValue, pSize, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
       case (logicType)
          call MPI_BCAST(CFG_varList(n)%logicValue, pSize, MPI_LOGICAL, root, MPI_COMM_WORLD, ierr)
       case (charType)
          do i = 1, pSize
             call MPI_BCAST(CFG_varList(n)%charValue(i), nameLen, MPI_CHARACTER, root, MPI_COMM_WORLD, ierr)
          end do
       end select

    end do

  end subroutine CFG_shareMPI

  !> Sort the variables for faster lookup
  subroutine CFG_sort()
    call QsortCFG(CFG_varList(1:CFG_nVars))
    CFG_sorted = .true.
  end subroutine CFG_sort

  !> Simple implementation of quicksort algorithm to sort the variable list alphabetically.
  recursive subroutine QsortCFG(list)
    type(CFG_varType), intent(inout) :: list(:)

    integer                          :: splitPos

    if (size(list) > 1) then
       call partitionVarList(list, splitPos)
       call QsortCFG( list(:splitPos-1) )
       call QsortCFG( list(splitPos:) )
    end if

  end subroutine QsortCFG

  !> Helper routine for quicksort, to perform partitioning
  subroutine partitionVarList(list, marker)
    type(CFG_varType), intent(inout) :: list(:)
    integer, intent(out)             :: marker

    integer                 :: left, right, pivot
    type(CFG_varType)       :: temp
    character(len=nameLen)  :: pivotValue

    left = 0
    right = size(list) + 1

    ! Take the middle element as pivot
    pivot = size(list) / 2
    pivotValue = list(pivot)%pName

    do while (left < right)

       right = right - 1
       do while (lgt(list(right)%pName, pivotValue))
          right = right - 1
       end do

       left = left + 1
       do while (lgt(pivotValue, list(left)%pName))
          left = left + 1
       end do

       if (left < right) then
          temp = list(left)
          list(left) = list(right)
          list(right) = temp
       end if
    end do

    if (left == right) then
       marker = left + 1
    else
       marker = left
    end if

  end subroutine partitionVarList

  !> Update the variables in the configartion with the values found in 'filename'
  subroutine CFG_updateFromFile(filename)
    character(len=*), intent(in)  :: filename

    integer                       :: ioState, equalSignPos, lineEnd
    integer                       :: n, ix, nL
    integer                       :: startIxs(maxVarSize), endIxs(maxVarSize), nEntries
    character(len=tinyLen)        :: configFormat
    character(len=nameLen)        :: pName
    character(len=lineLen)        :: line, tempString

    if (.not. CFG_initialized) then
       print *, "CFG_updateFromFile: Not CFG_initialized yet"
       return
    end if

    open(UNIT = 1, FILE = filename, STATUS = "OLD", ACTION = "READ", ERR = 999, IOSTAT = ioState)
    nL = 0

    ! Set the line format to read, only depends on lineLen currently
    write(configFormat, FMT = "(I6)") lineLen
    configFormat = "(A" // trim(adjustl(configFormat)) // ")"

    do
       read( UNIT = 1, FMT = configFormat, ERR = 999, end = 666) line; nL = nL + 1
       ! Line format should be "variable = value(s) [#Comment]"

       lineEnd = scan(line, '#') - 1       ! Don't use the part of the line after '#'
       if (lineEnd == -1) lineEnd = lineLen
       line = line(1:lineEnd)

       equalSignPos = scan(line, '=')      ! Locate the '=' sign, if there is no such sign then skip the line
       if (equalSignPos == 0) cycle

       pName = line(1 : equalSignPos - 1)  ! Set variable name
       pName = adjustl(pName)              ! Remove leading blanks

       line = line(equalSignPos + 1 : lineEnd)   ! Set line to the values behind the '=' sign
       line = adjustl(line)                      ! Remove leading blanks

       ! Find variable corresponding to name in file
       ix = getVarIndex(pName)

       if (ix > 0) then

          if (CFG_varList(ix)%varSize) then

             ! Get the start and end positions of the line content, and the number of entries
             call splitLineByDelims(line, " ,'"""//char(9), maxVarSize, nEntries, startIxs, endIxs)

             select case (CFG_varList(ix)%pType)

             case (intType)
                if (CFG_varList(ix)%pSize /= nEntries) then
                   deallocate( CFG_varList(ix)%intValue )
                   allocate( CFG_varList(ix)%intValue(nEntries) )
                   CFG_varList(ix)%pSize = nEntries
                end if

                do n = 1, nEntries
                   read(line(startIxs(n):endIxs(n)), *, ERR = 999) CFG_varList(ix)%intValue(n)
                end do

             case (floatType)
                if (CFG_varList(ix)%pSize /= nEntries) then
                   deallocate( CFG_varList(ix)%floatValue )
                   allocate( CFG_varList(ix)%floatValue(nEntries) )
                   CFG_varList(ix)%pSize = nEntries
                end if

                do n = 1, nEntries
                   read(line(startIxs(n):endIxs(n)), *, ERR = 999) CFG_varList(ix)%floatValue(n)
                end do

             case (charType)
                if (CFG_varList(ix)%pSize /= nEntries) then
                   deallocate( CFG_varList(ix)%charValue )
                   allocate( CFG_varList(ix)%charValue(nEntries) )
                   CFG_varList(ix)%pSize = nEntries
                end if

                do n = 1, nEntries
                   read(line(startIxs(n):endIxs(n)), *, ERR = 999) CFG_varList(ix)%charValue(n)
                end do

             case (logicType)
                if (CFG_varList(ix)%pSize /= nEntries) then
                   deallocate( CFG_varList(ix)%logicValue )
                   allocate( CFG_varList(ix)%logicValue(nEntries) )
                   CFG_varList(ix)%pSize = nEntries
                end if

                do n = 1, nEntries
                   read(line(startIxs(n):endIxs(n)), *, ERR = 999) CFG_varList(ix)%logicValue(n)
                end do

             end select

          else
             select case (CFG_varList(ix)%pType)
             case (intType)
                read(line, *, ERR = 999, end = 998) CFG_varList(ix)%intValue
             case (floatType)
                read(line, *, ERR = 999, end = 998) CFG_varList(ix)%floatValue
             case (charType)
                read(line, *, ERR = 999, end = 998) CFG_varList(ix)%charValue
             case (logicType)
                read(line, *, ERR = 999, end = 998) CFG_varList(ix)%logicValue
             end select
          end if

       else

          write(tempString, *) "variable [", trim(pName), "] from ", filename, " could not be updated"
          call printError("CFG_updateFromFile", tempString)

       end if

       ! Go to the next iteration
       cycle

    end do

666 continue ! Routine ends here if the end of "filename" is reached
    close(UNIT = 1, STATUS = "KEEP", ERR = 999, IOSTAT = ioState)
    return

    ! The routine only ends up here through an error
998 continue
    write(tempString, *) "Not enough values for variable [", trim(pName), "] in ", filename
    call printError("CFG_updateFromFile", tempString)
    return

999 continue
    write(tempString, *) "ioState = ", ioState, " while reading from ", filename, " at line ", nL
    call printError("CFG_updateFromFile", tempString)
    return

  end subroutine CFG_updateFromFile

  !> Return the index of the variable with name 'pName', or -1 if not found.
  integer function getVarIndex(pName)
    character(len=*), intent(in)  :: pName

    integer                       :: i

    if (.not. CFG_sorted) then

       getVarIndex = -1
       do i = 1, CFG_nVars
          if (CFG_varList(i)%pName == pName) then
             getVarIndex = i
             exit
          end if
       end do

    else if (CFG_sorted) then

       getVarIndex = bSearchVar(pName)

    end if

  end function getVarIndex

  !> Performa a binary search for the variable 'pName'
  integer function bSearchVar(pName)
    character(len=*), intent(in)  :: pName

    integer                       :: iMin, iMax, iMiddle
    logical                       :: success

    iMin     = 1
    iMax     = CFG_nVars
    success  = .false.

    do while (iMin <= iMax)
       iMiddle = iMin + (iMax - iMin) / 2

       if ( lgt(CFG_varList(iMiddle)%pName, pName) ) then
          iMax = iMiddle - 1
       else if ( llt(CFG_varList(iMiddle)%pName, pName) ) then
          iMin = iMiddle + 1
       else
          bSearchVar = iMiddle
          success = .true.
          exit
       end if

    end do

    if (.not. success) then
       bSearchVar = - 1
       print *, "bSearchVar warning: ", pName, " not found"
    end if

  end function bSearchVar

  !> Helper routine to create variables. This is useful because a lot of the same
  !! code is executed for the different types of variables.
  subroutine prepareStoreVar(pName, pType, pSize, comment, varSize)
    character(len=*), intent(in)  :: pName, comment
    integer, intent(in)           :: pType, pSize
    logical, intent(in), optional :: varSize

    call checkInitialized()
    CFG_sorted = .false.

    if (CFG_nVars >= CFG_maxVars ) call resizeVarList(2 * CFG_nVars)
    CFG_nVars = CFG_nVars + 1

    CFG_varList(CFG_nVars)%pName       = pName
    CFG_varList(CFG_nVars)%comment     = comment
    CFG_varList(CFG_nVars)%pType       = pType
    CFG_varList(CFG_nVars)%varSize   = .false.
    if (present(varSize)) CFG_varList(CFG_nVars)%varSize = varSize
    CFG_varList(CFG_nVars)%pSize       = pSize

    select case (CFG_varList(CFG_nVars)%pType)

    case (intType)
       allocate( CFG_varList(CFG_nVars)%intValue(pSize) )

    case (floatType)
       allocate( CFG_varList(CFG_nVars)%floatValue(pSize) )

    case (charType)
       allocate( CFG_varList(CFG_nVars)%charValue(pSize) )

    case (logicType)
       allocate( CFG_varList(CFG_nVars)%logicValue(pSize) )

    end select

  end subroutine prepareStoreVar

  integer function prepareGetVar(pName, pType, pSize)
    character(len=*), intent(in)  :: pName
    integer, intent(in)           :: pType, pSize

    integer                       :: ix
    logical                       :: correctType, correctSize

    ix                   = getVarIndex(pName)
    prepareGetVar  = ix
    correctType          = .false.
    correctSize          = .false.

    if (ix /= -1) then
       if (CFG_varList(ix)%pType == pType) correctType = .true.
       if (CFG_varList(ix)%pSize == pSize) correctSize = .true.
    end if

    if (ix == -1) then
       call printError("CFG_getVar", "prepareGetVar: ["//trim(pName)//"] not found")
    else if (.not. correctType) then
       call printError("CFG_getVar", "prepareGetVar: ["//trim(pName)//"] has different type")
    else if (.not. correctSize) then
       call printError("CFG_getVar", "prepareGetVar: ["//trim(pName)//"] has different size")
    end if

  end function prepareGetVar

  subroutine CFG_addVarDouble(pName, floatValue, comment)
    character(len=*), intent(in)  :: pName, comment
    double precision, intent(in)  :: floatValue
    call prepareStoreVar(pName, floatType, 1, comment)
    CFG_varList(CFG_nVars)%floatValue(1) = floatValue
  end subroutine CFG_addVarDouble

  subroutine CFG_addVarDoubleArray(pName, floatValue, comment, varSize)
    character(len=*), intent(in)  :: pName, comment
    double precision, intent(in)  :: floatValue(:)
    logical, intent(in), optional :: varSize
    call prepareStoreVar(pName, floatType, size(floatValue), comment, varSize)
    CFG_varList(CFG_nVars)%floatValue = floatValue
  end subroutine CFG_addVarDoubleArray

  subroutine CFG_addVarInt(pName, intValue, comment)
    character(len=*), intent(in)  :: pName, comment
    integer, intent(in)  :: intValue
    call prepareStoreVar(pName, intType, 1, comment)
    CFG_varList(CFG_nVars)%intValue(1) = intValue
  end subroutine CFG_addVarInt

  subroutine CFG_addVarIntArray(pName, intValue, comment, varSize)
    character(len=*), intent(in)  :: pName, comment
    integer, intent(in)           :: intValue(:)
    logical, intent(in), optional :: varSize
    call prepareStoreVar(pName, intType, size(intValue), comment, varSize)
    CFG_varList(CFG_nVars)%intValue = intValue
  end subroutine CFG_addVarIntArray

  subroutine CFG_addVarChar(pName, charValue, comment)
    character(len=*), intent(in)  :: pName, comment, charValue
    call prepareStoreVar(pName, charType, 1, comment)
    CFG_varList(CFG_nVars)%charValue(1) = charValue
  end subroutine CFG_addVarChar

  subroutine CFG_addVarCharArray(pName, charValue, comment, varSize)
    character(len=*), intent(in)  :: pName, comment, charValue(:)
    logical, intent(in), optional :: varSize
    call prepareStoreVar(pName, charType, size(charValue), comment, varSize)
    CFG_varList(CFG_nVars)%charValue = charValue
  end subroutine CFG_addVarCharArray

  subroutine CFG_addVarlogical(pName, logicValue, comment)
    character(len=*), intent(in)  :: pName, comment
    logical, intent(in)           :: logicValue
    call prepareStoreVar(pName, logicType, 1, comment)
    CFG_varList(CFG_nVars)%logicValue(1) = logicValue
  end subroutine CFG_addVarlogical

  subroutine CFG_addVarlogicalArray(pName, logicValue, comment, varSize)
    character(len=*), intent(in)  :: pName, comment
    logical, intent(in)           :: logicValue(:)
    logical, intent(in), optional :: varSize
    call prepareStoreVar(pName, logicType, size(logicValue), comment, varSize)
    CFG_varList(CFG_nVars)%logicValue = logicValue
  end subroutine CFG_addVarlogicalArray



  subroutine CFG_getVarDouble(pName, floatValue)
    character(len=*), intent(in)     :: pName
    double precision, intent(inout)  :: floatValue
    integer :: ix
    ix          = prepareGetVar(pName, floatType, 1)
    floatValue  = CFG_varList(ix)%floatValue(1)
  end subroutine CFG_getVarDouble

  subroutine CFG_getVarDoubleArray(pName, floatValue)
    character(len=*), intent(in)     :: pName
    double precision, intent(inout)  :: floatValue(:)
    integer :: ix
    ix          = prepareGetVar(pName, floatType, size(floatValue))
    floatValue  = CFG_varList(ix)%floatValue
  end subroutine CFG_getVarDoubleArray

  subroutine CFG_getVarInt(pName, intValue)
    character(len=*), intent(in)     :: pName
    integer, intent(inout)           :: intValue
    integer :: ix
    ix          = prepareGetVar(pName, intType, 1)
    intValue    = CFG_varList(ix)%intValue(1)
  end subroutine CFG_getVarInt

  subroutine CFG_getVarIntArray(pName, intValue)
    character(len=*), intent(in)     :: pName
    integer, intent(inout)           :: intValue(:)
    integer :: ix
    ix          = prepareGetVar(pName, intType, size(intValue))
    intValue    = CFG_varList(ix)%intValue
  end subroutine CFG_getVarIntArray

  subroutine CFG_getVarChar(pName, charValue)
    character(len=*), intent(in)     :: pName
    character(len=*), intent(inout)  :: charValue
    integer :: ix
    ix          = prepareGetVar(pName, charType, 1)
    charValue   = CFG_varList(ix)%charValue(1)
  end subroutine CFG_getVarChar

  subroutine CFG_getVarCharArray(pName, charValue)
    character(len=*), intent(in)     :: pName
    character(len=*), intent(inout)  :: charValue(:)
    integer :: ix
    ix          = prepareGetVar(pName, charType, size(charValue))
    charValue   = CFG_varList(ix)%charValue
  end subroutine CFG_getVarCharArray

  subroutine CFG_getVarlogical(pName, logicValue)
    character(len=*), intent(in)     :: pName
    logical, intent(inout)           :: logicValue
    integer :: ix
    ix          = prepareGetVar(pName, logicType, 1)
    logicValue   = CFG_varList(ix)%logicValue(1)
  end subroutine CFG_getVarlogical

  subroutine CFG_getVarlogicalArray(pName, logicValue)
    character(len=*), intent(in)     :: pName
    logical, intent(inout)           :: logicValue(:)
    integer :: ix
    ix          = prepareGetVar(pName, logicType, size(logicValue))
    logicValue   = CFG_varList(ix)%logicValue
  end subroutine CFG_getVarlogicalArray

  subroutine prepareVarTo(pName, pType, argIndex, ix, pIndex)
    character(len=*), intent(in)  :: pName
    integer, intent(in)           :: pType
    integer, intent(inout)        :: ix, pIndex
    integer, intent(in), optional :: argIndex

    logical :: correctSize, correctType
    correctSize = .false.
    correctType = .false.

    pIndex = 1
    if (present(argIndex)) pIndex = argIndex

    ix = getVarIndex(pName)

    if (ix /= -1) then
       if (CFG_varList(ix)%pType == pType) correctType = .true.
       if (pIndex <= CFG_varList(ix)%pSize .and. (present(argIndex) .or. CFG_varList(ix)%pSize == 1)) &
            & correctSize = .true.
    end if

    if (ix == -1) then
       call printError("CFG_varXXX", "prepareVarTo: ["//trim(pName)//"] not found")
    else if (.not. correctType) then
       call printError("CFG_varXXX", "prepareVarTo: ["//trim(pName)//"] has different type")
    else if (.not. correctSize) then
       call printError("CFG_varXXX", "prepareVarTo: ["//trim(pName)//"] has smaller size")
    end if

  end subroutine prepareVarTo

  double precision function CFG_varDble(pName, argIndex)
    character(len=*), intent(in)  :: pName
    integer, intent(in), optional :: argIndex
    integer :: ix, pIndex
    call prepareVarTo(pName, floatType, argIndex, ix, pIndex)
    CFG_varDble = CFG_varList(ix)%floatValue(pIndex)
  end function CFG_varDble

  integer function CFG_varInt(pName, argIndex)
    character(len=*), intent(in)  :: pName
    integer, intent(in), optional :: argIndex
    integer :: ix, pIndex
    call prepareVarTo(pName, intType, argIndex, ix, pIndex)
    CFG_varInt = CFG_varList(ix)%intValue(pIndex)
  end function CFG_varInt

  logical function CFG_varLogic(pName, argIndex)
    character(len=*), intent(in)  :: pName
    integer, intent(in), optional :: argIndex
    integer :: ix, pIndex
    call prepareVarTo(pName, logicType, argIndex, ix, pIndex)
    CFG_varLogic = CFG_varList(ix)%logicValue(pIndex)
  end function CFG_varLogic

  integer function CFG_getSize(pName)
    character(len=*), intent(in)  :: pName
    integer                       :: ix
    ix = getVarIndex(pName)
    if (ix /= -1) then
       CFG_getSize = CFG_varList(ix)%pSize
    else
       call printError("CFG_getSize", "CFG_getSize: variable ["//trim(pName)//"] not found")
    end if
  end function CFG_getSize

  integer function CFG_getType(pName)
    character(len=*), intent(in)  :: pName
    integer                       :: ix
    ix = getVarIndex(pName)
    if (ix /= -1) then
       CFG_getType = CFG_varList(ix)%pType
    else
       call printError("CFG_getType", "CFG_getSize: variable ["//trim(pName)//"] not found")
    end if
  end function CFG_getType

  subroutine CFG_write(filename)
    character(len=*), intent(in)  :: filename

    integer                       :: i, j, ioState
    character(len=tinyLen)        :: nameFormat
    integer, parameter            :: myUnit = 333

    write(nameFormat, FMT = "(I6)") nameLen
    nameFormat = "(A,A" // trim(adjustl(nameFormat)) // ",A)"

    if (.not. CFG_initialized) then
       print *, "CFG_updateFromFile: Not CFG_initialized yet"
       return
    end if

    open(myUnit, FILE = filename, ACTION = "WRITE", ERR = 999, IOSTAT = ioState)

    write(myUnit, ERR = 999, FMT = "(A)") " ##############################################"
    write(myUnit, ERR = 999, FMT = "(A)") " ###          Configuration file            ###"
    write(myUnit, ERR = 999, FMT = "(A)") " ##############################################"
    write(myUnit, ERR = 999, FMT = "(A)") ""

    do i = 1, CFG_nVars
       write(myUnit, ERR = 999, FMT = "(A,A,A)") " # ", trim(CFG_varList(i)%comment), ":"
       write(myUnit, ADVANCE = "NO", ERR = 999, FMT = nameFormat) " ", CFG_varList(i)%pName, " = "

       select case(CFG_varList(i)%pType)
       case (intType)
          do j = 1, CFG_varList(i)%pSize
             write(myUnit, ADVANCE = "NO", ERR = 999, FMT = "(I10, A) ") CFG_varList(i)%intValue(j), "  "
          end do
       case (floatType)
          do j = 1, CFG_varList(i)%pSize
             write(myUnit, ADVANCE = "NO", ERR = 999, FMT = "(E10.4, A) ") CFG_varList(i)%floatValue(j), "  "
          end do
       case (charType)
          do j = 1, CFG_varList(i)%pSize
             write(myUnit, ADVANCE = "NO", ERR = 999, FMT = "(A, A)") &
                  & '"' // trim(CFG_varList(i)%charValue(j)) // '"', "  "
          end do
       case (logicType)
          do j = 1, CFG_varList(i)%pSize
             write(myUnit, ADVANCE = "NO", ERR = 999, FMT = "(L10, A) ") CFG_varList(i)%logicValue(j), "  "
          end do
       end select
       write(myUnit, ERR = 999, FMT = "(A)") ""
       write(myUnit, ERR = 999, FMT = "(A)") ""
    end do
    write(myUnit, ERR = 999, FMT = "(A)") ""

    close(myUnit, STATUS = "KEEP", ERR = 999, IOSTAT = ioState)

    return

999 continue ! If there was an error, the routine will end here
    print *, "CFG_write error:"
    print *, "ioState = ", ioState, " while writing to ", filename

  end subroutine CFG_write

  subroutine resizeVarList(newSize)
    integer, intent(in)              :: newSize
    type(CFG_varType), allocatable   :: CFG_varListCopy(:)

    allocate(CFG_varListCopy(CFG_nVars))
    CFG_varListCopy = CFG_varList(1:CFG_nVars)
    deallocate(CFG_varList)

    allocate(CFG_varList(newSize))
    CFG_varList(1:CFG_nVars) = CFG_varListCopy
    deallocate(CFG_varListCopy)

    CFG_maxVars = newSize
  end subroutine resizeVarList

  subroutine checkInitialized()
    if( .not. CFG_initialized ) call CFG_initialize(100)
  end subroutine checkInitialized

  subroutine printError(fName, message)
    character(len=*) :: fName, message

    print *, "In function <", fname, "> an error occurred, with the following message:"
    print *, message
    print *, "Will halt now"
    stop

  end subroutine printError

  !> Routine to help read a variable number of entries from a line.
  ! In arguments:
  !  line           The line from which we want to read
  !  delims         A string with delimiters. For example delims = " ,'"""//char(9)
  !                 " ,'"""//char(9) = space, comma, single/double quotation marks, tab
  !  nEntriesMax    Maximum number of entries to read in
  !
  ! Out arguments:
  !  nEntries       Number of entries found
  !  startIxs       On return, startIxs(i) holds the starting point of entry i
  !  endIxs         On return, endIxs(i) holds the end point of entry i
  subroutine splitLineByDelims(line, delims, nEntriesMax, nEntries, startIxs, endIxs)
    character(len=*), intent(in)  :: line, delims
    integer, intent(in)           :: nEntriesMax
    integer, intent(inout)        :: nEntries, startIxs(nEntriesMax), endIxs(nEntriesMax)

    integer                       :: ix, prevIx

    prevIx   = 0
    nEntries = 0

    do while (nEntries < nEntriesMax)

       ! Find the starting point of the next entry (a non-delimiter value)
       ix                   = verify(line(prevIx+1:), delims)
       if (ix == 0) exit                      ! No more entries

       nEntries             = nEntries + 1
       startIxs(nEntries)   = prevIx + ix     ! This is the absolute position in 'line'

       ! Get the end point of the current entry (next delimiter index minus one)
       ix = scan(line(startIxs(nEntries)+1:), delims) - 1

       if (ix == -1) then                     ! If there is no last delimiter,
          endIxs(nEntries)  = len(line)       ! the end of the line is the endpoint
       else
          endIxs(nEntries)  = startIxs(nEntries) + ix
       end if

       prevIx = endIxs(nEntries)              ! We continue to search from here
    end do

  end subroutine splitLineByDelims

end module module_config
