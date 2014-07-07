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


module module_EfieldAMR

   use m_write_silo
   use module_constants
   use module_config
   use generalUtilities
   use module_gas
   use mpi

   implicit none
   private

   type amrGrid
      logical :: hasChild

      integer :: nPoints(3)
      integer :: child_ixMin(3)
      integer :: child_ixMax(3)
      integer :: cellIxOffset

      double precision :: delta(3), invDelta(3)
      double precision :: posMin(3), posMax(3)
      double precision, allocatable :: potential(:,:,:)
      double precision, allocatable :: elecDens(:,:,:)
      double precision, allocatable :: posIonDens(:,:,:)
      double precision, allocatable :: negIonDens(:,:,:)
      double precision, allocatable :: O2MinDens(:,:,:)
      double precision, allocatable :: excDens(:,:,:)
      double precision, allocatable :: sourceTerm(:,:,:)
      double precision, allocatable :: Efield(:,:,:,:)
      double precision, allocatable :: temp(:,:,:)
   end type amrGrid

   !> The size of the workspace array needed by the external Poisson solver
   integer           :: nWorkspace_allocated = 0
   integer           :: EFA_minGridSize
   integer           :: EFA_maxAllowedLvl
   integer           :: EFA_maxLvl
   integer           :: EFA_noDeRefLvl
   integer           :: EFA_init_RefUpToLvl
   integer           :: EFA_bcType

   double precision  :: EFA_minElecDensForRef
   double precision  :: EFA_maxPosGridVel(3)
   double precision  :: EFA_maxNegGridVel(3)
   double precision  :: EFA_outsideGridRatio

   logical           :: EFA_doCorrect
   logical           :: EFA_useElectrode

   !> Workspace for the external Poisson solver
   double precision, allocatable :: workspace(:)

   !> Electric field values
   double precision, allocatable :: EFA_efield_values(:)

   !> Electric field times (at which the values given above are used)
   double precision, allocatable :: EFA_efield_times(:)

   !> List of the max. delta vs. electric field
   double precision, allocatable :: EFA_deltaVsEfield(:,:)

   !> The list of grids
   type(amrGrid), allocatable    :: grids(:)

   !> Anbang: Relative permitivity
   double precision :: epsilonRelative

   public :: EFA_initialize, EFA_compute_Efield, EFA_get_potential
   public :: EFA_writeGridsToFile, EFA_getEfield, EFA_printGridStructure

   public :: EFA_addElectron, EFA_addPosIon, EFA_addNegIon, EFA_addO2MinIon
   public :: EFA_resetElecDens, EFA_resetIonDens
   public :: EFA_addGrid, EFA_removeGrid, EFA_reshapeGrid
   public :: EFA_adaptGridStructure
   public :: EFA_getCellIndexAt, EFA_getMaxCellIx
   public :: EFA_getNumElecAt
   public :: EFA_getElecDensAt
   public :: EFA_getMaxDeltaAt
   public :: EFA_getMaxLvl
   public :: EFA_allCollectTemp
   public :: EFA_resetTemp
   public :: EFA_addParticleCount
   public :: EFA_setNDecays
   public :: EFA_getDelta
   public :: EFA_getGridSize
   public :: EFA_getPosMin
   public :: EFA_addExcited
   public :: EFA_EfieldMax
   public :: EFA_setNegIonLoss
   public :: EFA_calibrateElectrode
   public :: EFA_shareElecDens
   public :: EFA_getPosIonDensAt
   public :: EFA_addO2MinDens
   public :: EFA_addPosIonDens
   public :: EFA_collectDensitiesAtRoot
   public :: EFA_setNeutralMoleLoss

contains

   subroutine EFA_initialize()

      integer              :: gridSize(3), varSize, lvl
      double precision     :: gridLength(3), relPosAndSize(6), posMin(3), posMax(3)
      double precision     :: parentLength(3), childLength(3), newCenter(3)
      double precision     :: coarseMin(3), coarseMax(3), coarseDelta(3)
      character(len=100)   :: cfgString

      call CFG_getVar("grid_size", gridSize)
      call CFG_getVar("grid_delta", gridLength)
      gridLength = gridLength * (gridSize - 1)
      call CFG_getVar("ref_minGridSize", EFA_minGridSize)

      varSize = CFG_getSize("sim_efield_times")
      allocate(EFA_efield_times(varSize))
      allocate(EFA_efield_values(varSize))
      CALL CFG_getVar("sim_efield_values", EFA_efield_values)
      CALL CFG_getVar("sim_efield_times", EFA_efield_times)

      EFA_bcType           = CFG_varInt("elec_boundaryType")
      EFA_maxAllowedLvl    = CFG_varInt("ref_maxLevels")
      EFA_minElecDensForRef= CFG_varDble("ref_MinElecDens")
      EFA_maxLvl           = -2
      EFA_doCorrect        = CFG_varLogic("ref_doCorrectionStep")
      EFA_useElectrode     = CFG_varLogic("sim_useElectrode")
      EFA_noDeRefLvl       = CFG_varInt("ref_NoDerefineUpToLvl")
      EFA_init_RefUpToLvl  = CFG_varInt("ref_init_UpToLvl")
      EFA_outsideGridRatio = CFG_varDble("ref_outsideGridRatio")
      varSize              = CFG_getSize("ref_deltaValues")
      epsilonRelative      = CFG_varDble("ref_epsilonRelative")

      if (EFA_init_refUpToLvl > EFA_noDeRefLvl .or. EFA_init_refUpToLvl > 4) then
         print *, "EFA_init_refUpToLvl should not be larger than EFA_noDeRefLvl or larger than 4"
         stop
      end if

      if (varSize /= CFG_getSize("ref_maxEfieldAtDelta")) then
         print *, "Make sure ref_EfieldValues and ref_deltaValues have the same size"
         stop
      end if

      allocate( EFA_deltaVsEfield(2, varSize) )
      call CFG_getVar("ref_deltaValues", EFA_deltaVsEfield(1, :))
      call CFG_getVar("ref_maxEfieldAtDelta", EFA_deltaVsEfield(2, :))

      call CFG_getVar("ref_maxNegativeGridVelocity", EFA_maxNegGridVel)
      call CFG_getVar("ref_maxPositiveGridVelocity", EFA_maxPosGridVel)

      if (allocated(grids)) deallocate(grids)
      allocate( grids(-1:EFA_maxAllowedLvl) )

      if (EFA_init_RefUpToLvl > EFA_maxAllowedLvl) then
         print *, "Error: EFA_init_RefUpToLvl is larger than EFA_maxAllowedLvl"
         stop
      end if

      ! Level -1 is there to provide boundary conditions to lvl 0
      coarseDelta = EFA_outsideGridRatio * gridLength / (gridSize - 1)
      coarseMin   = nint((-(EFA_outsideGridRatio - 1.0d0)/(2.0d0 * EFA_outsideGridRatio)) * gridSize) * coarseDelta
      coarseMax   = coarseMin + (gridSize-1) * coarseDelta
      call setGridlevel(-1, gridSize, coarseMin, coarseMax, .true.)

      ! Level 0 is the default unrefined grid
      call setGridlevel(0, gridSize, (/0.0d0, 0.0d0, 0.0d0/), gridLength, .false.)

      ! Create the initially refined grids
      do lvl = 1, EFA_init_refUpToLvl

         ! Create the name of the configuration variable that controls the size of the refined grid
         write(cfgString, "(I1)") lvl
         cfgString = "ref_init_relPosAndSize_" // trim(adjustl(cfgString))

         ! Determine the center of the new grid relative to its parent, and its length.
         call CFG_getVar(trim(cfgString), relPosAndSize)

         parentLength   = grids(lvl-1)%posMax - grids(lvl-1)%posMin
         childLength    = relPosAndSize(4:6) * parentLength
         newCenter      = grids(lvl-1)%posMin + relPosAndSize(1:3) * parentLength

         posMin         = newCenter - 0.5d0 * childLength
         posMax         = newCenter + 0.5d0 * childLength

         call EFA_addGrid(lvl, posMin, posMax)
      end do

   end subroutine EFA_initialize

   subroutine EFA_calibrateElectrode(myrank, root)
      use module_electrode
      integer, intent(in) :: myrank, root

      integer :: ix, lvl, maxlvl, nElecPoints
      double precision :: xyz(3)
      double precision, allocatable :: delta_voltage(:)

      nElecPoints = EL_getNPoints()

      ! Gather the electron and ion densities
      do lvl = -1, EFA_maxLvl
         call EFA_collectDensitiesAtRoot(lvl, myrank, root)

         if (myrank == root) then
            call setSourceTermFromDens(lvl)
         end if
      end do

      ! Make a copy of the potential, restore it later
      do lvl = -1, EFA_maxLvl
         grids(lvl)%temp = grids(lvl)%potential
      end do

      ! Only the root task is responsible for the electrode charges and executes the code below
      if (myrank /= root) return

      allocate( delta_voltage(nElecPoints) )

      ! Compute potential at electrode points without any electrode charges
      call compute_potential_allLevels(0.0d0)

      do ix = 1, nElecPoints
         call EL_getSurfacePoint(ix, xyz)
         maxlvl = getMaxLevelFrom(xyz, -1)

         delta_voltage(ix) = get_value_at(grids(maxlvl), grids(maxlvl)%potential, xyz)
         !          print *, delta_voltage(ix)
      end do

      ! Now compute the potential at electrode points when they have unit charge
      do ix = 1, nElecPoints
         call EL_getSurfacePoint(ix, xyz)
         maxlvl = getMaxLevelFrom(xyz, -1)

         ! We can also add electrode outside the 'normal' domain
         do lvl = -1, maxlvl
            call addParticleTo(grids(lvl), grids(lvl)%sourceTerm, xyz, elecChargeOverEps0 / epsilonRelative)
         end do
      end do

      call compute_potential_allLevels(0.0d0)

      do ix = 1, nElecPoints
         call EL_getSurfacePoint(ix, xyz)
         maxlvl = getMaxLevelFrom(xyz, -1)

         delta_voltage(ix) = get_value_at(grids(maxlvl), grids(maxlvl)%potential, xyz) - delta_voltage(ix)
         call EL_setWeightFactor(ix, 1.0d0 / delta_voltage(ix))
      end do

      ! Copy back old potential
      do lvl = -1, EFA_maxLvl
         grids(lvl)%potential = grids(lvl)%temp
      end do
      deallocate( delta_voltage )

   end subroutine EFA_calibrateElectrode

   subroutine syncGridStructure(myrank, root)
      integer, intent(in)           :: myrank, root

      integer                       :: lvl, ix, dataSize, ierr
      integer, allocatable          :: nPointsArray(:), childIxMinArray(:), childIxMaxArray(:)
      double precision, allocatable :: posMinArray(:), posMaxArray(:), deltaArray(:)

      CALL MPI_BCAST(EFA_maxLvl, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)

      ! Actually there is no need to send data about lvl 0, but it keeps the code simpler
      dataSize = 3 * (EFA_maxLvl + 1)

      allocate( nPointsArray(dataSize) )
      allocate( childIxMinArray(dataSize) )
      allocate( childIxMaxArray(dataSize) )
      allocate( posMinArray(dataSize) )
      allocate( posMaxArray(dataSize) )
      allocate( deltaArray(dataSize) )

      if (myrank == root) then
         do lvl = 0, EFA_maxLvl
            ix = lvl * 3 + 1
            nPointsArray(ix:ix+2)      = grids(lvl)%nPoints
            childIxMinArray(ix:ix+2)   = grids(lvl)%child_ixMin
            childIxMaxArray(ix:ix+2)   = grids(lvl)%child_ixMax
            posMinArray(ix:ix+2)       = grids(lvl)%posMin
            posMaxArray(ix:ix+2)       = grids(lvl)%posMax
            deltaArray(ix:ix+2)        = grids(lvl)%delta
         end do
      end if

      CALL MPI_BCAST(nPointsArray, dataSize, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(childIxMinArray, dataSize, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(childIxMaxArray, dataSize, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)

      CALL MPI_BCAST(posMinArray, dataSize, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(posMaxArray, dataSize, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(deltaArray, dataSize, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)

      if (myrank /= root) then

         grids(0:EFA_maxLvl-1)%hasChild = .true.
         grids(EFA_maxLvl:)%hasChild   = .false.

         do lvl = 0, EFA_maxLvl
            ix = lvl * 3 + 1
            grids(lvl)%nPoints         = nPointsArray(ix:ix+2)
            grids(lvl)%child_ixMin     = childIxMinArray(ix:ix+2)
            grids(lvl)%child_ixMax     = childIxMaxArray(ix:ix+2)
            grids(lvl)%posMin          = posMinArray(ix:ix+2)
            grids(lvl)%posMax          = posMaxArray(ix:ix+2)
            grids(lvl)%delta           = deltaArray(ix:ix+2)
            grids(lvl)%invDelta        = 1.0d0 / deltaArray(ix:ix+2)
            call ensureGridIsAllocated(grids(lvl))
         end do

      end if

      deallocate( nPointsArray )
      deallocate( childIxMinArray)
      deallocate( childIxMaxArray )
      deallocate( posMinArray )
      deallocate( posMaxArray )
      deallocate( deltaArray )

   end subroutine syncGridStructure

   subroutine EFA_adaptGridStructure(dt, myrank, root)
      use module_electrode
      double precision, intent(in)  :: dt
      integer, intent(in)           :: myrank, root

      logical              :: refine, accurate
      integer              :: i, j, k, lvl
      integer              :: nIx(3), ix(3)
      integer              :: iRefMin(3), iRefMax(3), diffSize(3)
      integer              :: maxIncrPerStep(3), maxDecrPerStep(3)
      integer              :: nElecPoints, maxlvl
      double precision     :: xyz(3), delta(3), posMin(3)
      double precision     :: Ef1(3), Ef2(3), maxValue
      double precision     :: truncError, maxTruncError, maxEfieldAtLvl

      ! Gather the electron and ion densities
      do lvl = 0, EFA_maxLvl
         call EFA_collectDensitiesAtRoot(lvl, myrank, root)
         call collectExcDensAtRoot(lvl, myrank, root)
      end do

      if (myrank /= root) then
         call syncGridStructure(myrank, root)
         call setCellIxOffsets()
         return
      end if

      ! Estimate the Efield error in the current maximum level
      lvl      = 0
      accurate = .false.
      maxValue = 0.0d0

      do while (lvl <= EFA_init_refUpToLvl .or. (.not. accurate .and. lvl <= EFA_maxAllowedLvl))
         nIx            = grids(lvl)%nPoints
         delta          = grids(lvl)%delta
         posMin         = grids(lvl)%posMin

         call setSourceTermFromDens(lvl)

         if (EFA_useElectrode) then
            nElecPoints = EL_getNPoints()

            do i = 1, nElecPoints
               call EL_getSurfacePoint(i, xyz)
               maxlvl = getMaxLevelFrom(xyz, -1)

               ! Add electrode charge if it is inside the current level
               if (maxlvl >= lvl) then
                  call addParticleTo(grids(lvl), grids(lvl)%sourceTerm, xyz, EL_getCharge(i))
               end if
            end do
         end if

         if (lvl > 0) call setBndCndFromParent(lvl)
         call computePotentialAtLevel(lvl)
         call computeEfieldAtLvl(lvl)

         ! Find maximum Efield for the delta at this level
         call linearInterpolateList(EFA_deltaVsEfield(1, :), EFA_deltaVsEfield(2, :), &
              maxval(delta), maxEfieldAtLvl)

         iRefMin        = nIx+1
         iRefMax        = 0
         accurate       = .true.
         maxTruncError  = 0.0d0

         do k = 2, nIx(3)-1
            do j = 2, nIx(2)-1
               do i = 2, nIx(1)-1
                  ix = (/i, j, k/)

                  Ef1 = grids(lvl)%Efield(:, i, j, k)

                  if (grids(lvl)%elecDens(i,j,k) >= EFA_minElecDensForRef .and. &
                       sqrt(sum(Ef1**2)) > maxEfieldAtLvl) then
                     where (ix < iRefMin)
                        iRefMin = ix
                     end where

                     where (ix > iRefMax)
                        iRefMax = ix
                     end where

                     accurate = .false.
                  end if

               end do
            end do
         end do

         iRefMin = max(2, iRefMin - EFA_minGridSize/2)
         iRefMax = min(nIx-1, iRefMax + EFA_minGridSize/2)

         ! Up to EFA_noDeRefLvl grids do not de-refine
         if (accurate .and. lvl >= EFA_noDeRefLvl) then

            if (grids(lvl)%hasChild) then
               call EFA_removeGrid(lvl+1)
            end if

         else if (.not. accurate .and. lvl < EFA_maxAllowedLvl) then
            ! Create or move the grid's child
            if (.not. grids(lvl)%hasChild) then
               ! Make sure the refined region is large enough
               do i = 1, 3
                  do
                     diffSize(i) = EFA_minGridSize - (iRefMax(i) - iRefMin(i)) + 2
                     if (diffSize(i) <= 0) exit
                     iRefMin(i) = max(2, iRefMin(i)-1)
                     iRefMax(i) = min(nIx(i)-1, iRefMax(i)+1)
                  end do
               end do

               call createChildForLvl(lvl, iRefMin, iRefMax)
            else
               maxIncrPerStep = ceiling(EFA_maxPosGridVel * dt / delta)
               maxDecrPerStep = ceiling(EFA_maxNegGridVel * dt / delta)

               ! Limit the movement of the child
               where (iRefMin < grids(lvl)%child_ixMin - maxDecrPerStep)
                  iRefMin = grids(lvl)%child_ixMin - maxDecrPerStep
               end where

               where (iRefMin > grids(lvl)%child_ixMin + maxIncrPerStep)
                  iRefMin = grids(lvl)%child_ixMin + maxIncrPerStep
               end where

               where (iRefMax < grids(lvl)%child_ixMax - maxDecrPerStep)
                  iRefMax = grids(lvl)%child_ixMax - maxDecrPerStep
               end where

               where (iRefMax > grids(lvl)%child_ixMax + maxIncrPerStep)
                  iRefMax = grids(lvl)%child_ixMax + maxIncrPerStep
               end where

               ! Up to EFA_noDeRefLvl grids do not de-refine
               if (lvl < EFA_noDeRefLvl) then
                  where (iRefMin > grids(lvl)%child_ixMin)
                     iRefMin = grids(lvl)%child_ixMin
                  end where

                  where (iRefMax < grids(lvl)%child_ixMax)
                     iRefMax = grids(lvl)%child_ixMax
                  end where
               end if

               ! Make sure the refined region is large enough
               do i = 1, 3
                  do
                     diffSize(i) = EFA_minGridSize - (iRefMax(i) - iRefMin(i)) + 2
                     if (diffSize(i) <= 0) exit
                     iRefMin(i) = max(2, iRefMin(i)-1)
                     iRefMax(i) = min(nIx(i)-1, iRefMax(i)+1)
                  end do
               end do

               call EFA_reshapeGrid(lvl+1, (iRefMin-1) * delta + posMin, &
                    (iRefMax-1) * delta + posMin)
            end if

         end if

         lvl = lvl + 1

      end do

      call syncGridStructure(myrank, root)
      call setCellIxOffsets()

      ! do lvl = 0, EFA_maxLvl
      !    print *, "At lvl", lvl, "size:", grids(lvl)%nPoints
      ! end do

   end subroutine EFA_adaptGridStructure

   subroutine EFA_addGrid(lvl, posMin, posMax)
      integer, intent(in) :: lvl
      double precision, intent(in) :: posMin(3), posMax(3)

      integer :: pIxMin(3), pIxMax(3)

      ! The grid should be added to the previously highest level
      if (lvl /= EFA_maxLvl + 1) then
         print *, "EFA_addGrid error: not adding to highest current level"
         stop
      end if

      ! Determine the indices of the new grid in the parent grid at lvl-1
      pIxMin = nint((posMin - grids(lvl-1)%posMin) / grids(lvl-1)%delta) + 1
      pIxMax = nint((posMax - grids(lvl-1)%posMin) / grids(lvl-1)%delta) + 1

      if (any(pIxMin < 2) .or. any(pIxMax > grids(lvl-1)%nPoints-1)) then
         print *, "EFA_addGrid error: cannot add bigger grid than parent"
         stop
      end if

      call createChildForLvl(lvl-1, pIxMin, pIxMax)

   end subroutine EFA_addGrid

   recursive subroutine EFA_removeGrid(lvl)
      integer, intent(in) :: lvl

      if (lvl <= 0 .or. lvl > EFA_maxLvl) then
         print *, "removeGrid error: cannot remove that grid"
         print *, "Lvl:", lvl, "EFA_maxLvl:", EFA_maxLvl
         stop
      end if

      if (grids(lvl)%hasChild) then
         call EFA_removeGrid(lvl+1)
      end if

      call removeGridData(grids(lvl))
      grids(lvl-1)%hasChild = .false.

      EFA_maxLvl = lvl - 1

   end subroutine EFA_removeGrid

   recursive subroutine EFA_reshapeGrid(lvl, wantedPosMin, wantedPosMax)
      integer, intent(in) :: lvl
      double precision, intent(in) :: wantedPosMin(3), wantedPosMax(3)

      integer :: pIxMin(3), pIxMax(3), nIx(3), ixOffset(3)
      integer :: oldIxMin(3), oldIxMax(3), newChildMinIx(3), newChildMaxIx(3)
      integer :: i, j, k, i0, j0, k0
      double precision :: xyz(3), newPosMin(3)
      double precision, allocatable :: tempPIonDens(:,:,:), tempNIonDens(:,:,:), tempO2MinDens(:, :, :)
      double precision, allocatable :: tempExcDens(:,:,:), tempElecDens(:,:,:)

      if (lvl <= 0 .or. lvl > EFA_maxLvl) then
         print *, "EFA_reshapeGrid error: cannot reshape that grid"
         print *, "Lvl:", lvl, "EFA_maxLvl:", EFA_maxLvl
         stop
      end if

      ! Determine the new indices of the grid in the parent grid at lvl-1
      pIxMin = nint((wantedPosMin - grids(lvl-1)%posMin) / grids(lvl-1)%delta) + 1
      pIxMax = nint((wantedPosMax - grids(lvl-1)%posMin) / grids(lvl-1)%delta) + 1
      newPosMin = (pIxMin-1) * grids(lvl-1)%delta + grids(lvl-1)%posMin

      if (any(pIxMin < 2) .or. any(pIxMax > grids(lvl-1)%nPoints-1)) then
         print *, "EFA_reshapeGrid error: cannot cannot move grid outside parent"
         print *, "lvl", lvl, "pIxMin:", pIxMin, "pIxMax:", pIxMax
         stop
      end if

      ! Determine the new range of indices of the grid
      nIx = 2 * (pIxMax - pIxMin) + 1

      if (any(nIx < EFA_minGridSize)) then
         print *, "EFA_reshapeGrid error: new grid too small"
         print *, nIx, EFA_minGridSize
         stop
      end if

      ! Determine the range of indices of the moved grid that already existed in the unmoved grid
      oldIxMin = 1   - 2 * (pIxMin - grids(lvl-1)%child_ixMin)
      oldIxMax = nIx - 2 * (pIxMax - grids(lvl-1)%child_ixMax)

      ! Set the parent's child coordinates to the new values
      grids(lvl-1)%child_ixMin = pIxMin
      grids(lvl-1)%child_ixMax = pIxMax

      ! Determine the offset in indices
      i0 = oldIxMin(1) - 1
      j0 = oldIxMin(2) - 1
      k0 = oldIxMin(3) - 1

      allocate( tempPIonDens(nIx(1), nIx(2), nIx(3)) )
      allocate( tempNIonDens(nIx(1), nIx(2), nIx(3)) )
      allocate( tempElecDens(nIx(1), nIx(2), nIx(3)) )
      allocate( tempExcDens(nIx(1), nIx(2), nIx(3)) )
      allocate( tempO2MinDens(nIx(1), nIx(2), nIx(3)) )

      tempPIonDens   = 0.0d0
      tempNIonDens   = 0.0d0
      tempO2MinDens  = 0.0d0
      tempElecDens   = 0.0d0
      tempExcDens    = 0.0d0

      do k = 1, nIx(3)
         do j = 1, nIx(2)
            do i = 1, nIx(1)

               ! Use this part of the 'old' grid for the moved grid.
               if (all((/i,j,k/) > oldIxMin) .and. all((/i,j,k/) < oldIxMax)) then
                  tempPIonDens(i,j,k)  = grids(lvl)%posIonDens( i-i0, j-j0, k-k0)
                  tempNIonDens(i,j,k)  = grids(lvl)%negIonDens( i-i0, j-j0, k-k0)
                  tempO2MinDens(i,j,k) = grids(lvl)%O2MinDens( i-i0, j-j0, k-k0)
                  tempElecDens(i,j,k)  = grids(lvl)%elecDens(i-i0, j-j0, k-k0)
                  tempExcDens(i,j,k)   = grids(lvl)%excDens(i-i0, j-j0, k-k0)
               else  ! Interpolate the ion density from the parent grid in newly refined parts
                  xyz = (/i-1, j-1, k-1/) * grids(lvl)%delta + newPosMin
                  tempPIonDens(i,j,k)  = get_value_at(grids(lvl-1), grids(lvl-1)%posIonDens, xyz)
                  tempNIonDens(i,j,k)  = get_value_at(grids(lvl-1), grids(lvl-1)%negIonDens, xyz)
                  tempO2MinDens(i,j,k) = get_value_at(grids(lvl-1), grids(lvl-1)%O2MinDens, xyz)
                  tempElecDens(i,j,k)  = get_value_at(grids(lvl-1), grids(lvl-1)%elecDens, xyz)
                  tempExcDens(i,j,k)   = get_value_at(grids(lvl-1), grids(lvl-1)%excDens, xyz)
               end if

            end do
         end do
      end do

      grids(lvl)%nPoints      = nIx
      grids(lvl)%posMin       = newPosMin
      grids(lvl)%posMax       = (pIxMax-1) * grids(lvl-1)%delta + grids(lvl-1)%posMin

      call ensureGridIsAllocated(grids(lvl))
      grids(lvl)%posIonDens   = tempPIonDens
      grids(lvl)%negIonDens   = tempNIonDens
      grids(lvl)%O2MinDens    = tempO2MinDens
      grids(lvl)%elecDens     = tempElecDens
      grids(lvl)%excDens      = tempExcDens

      deallocate(tempPIonDens)
      deallocate(tempNIonDens)
      deallocate(tempElecDens)
      deallocate(tempExcDens)
      deallocate(tempO2MinDens)

      if (grids(lvl)%hasChild) then

         ! The children of this grid have a new relative position, differing by oldIxMin-1
         grids(lvl)%child_ixMin = grids(lvl)%child_ixMin + oldIxMin - 1
         grids(lvl)%child_ixMax = grids(lvl)%child_ixMax + oldIxMin - 1

         if (any(grids(lvl)%child_ixMin < 2) .or. any(grids(lvl)%child_ixMax > nIx-1)) then
            print *, "EFA_reshapeGrid warning: reshaping children outside parent region"

            newChildMinIx = grids(lvl)%child_ixMin
            newChildMaxIx = grids(lvl)%child_ixMax

            where (newChildMinIx < 2)
               newChildMinIx = 2
            end where

            where (newChildMaxIx > nIx-1)
               newChildMaxIx = nIx-1
            end where

            ! Put the child in the right position
            call EFA_reshapeGrid(lvl+1, (newChildMinIx-1) * grids(lvl)%delta + newPosMin, &
                 (newChildMaxIx-1) * grids(lvl)%delta + newPosMin)
         end if

      end if

   end subroutine EFA_reshapeGrid

   subroutine EFA_compute_Efield(time, myrank, root)
      use module_electrode
      double precision, intent(in) :: time
      integer, intent(in)  :: myrank, root
      integer              :: ix, lvl, maxlvl, nElecPoints, pointsUsed
      integer              :: dataSize, ierr
      double precision     :: xyz(3), voltageDiff, oldCharge, newCharge, voltage
      double precision     :: voltageRMS, temp, maxErrPos(3)

      ! Gather the electron and ion densities
      do lvl = -1, EFA_maxLvl
         call EFA_collectDensitiesAtRoot(lvl, myrank, root)

         if (myrank == root) then
            call setSourceTermFromDens(lvl)
         end if
      end do

      if (EFA_useElectrode .and. myrank == root) then ! Add the electrode charges
         voltage = EL_getVoltage(time)
         !          call compute_potential_allLevels()

         nElecPoints = EL_getNPoints()

         do lvl = -1, EFA_maxLvl
            grids(lvl)%temp = 0.0d0
         end do

         do ix = 1, nElecPoints
            call EL_getSurfacePoint(ix, xyz)
            maxlvl = getMaxLevelFrom(xyz, -1)

            voltageDiff = voltage - get_value_at(grids(maxlvl), grids(maxlvl)%potential, xyz)
            oldCharge   = EL_getCharge(ix)
            newCharge   = oldCharge + voltageDiff * EL_getWeightFactor(ix) * elecChargeOverEps0 / epsilonRelative
            call EL_setCharge(ix, newCharge)

            do lvl = -1, maxlvl
               call addParticleTo(grids(lvl), grids(lvl)%temp, xyz, newCharge)
            end do
         end do

         do lvl = -1, EFA_maxLvl
            grids(lvl)%sourceTerm = grids(lvl)%sourceTerm + grids(lvl)%temp
         end do
      end if

      ! The root task computes the potential
      if (myrank == root) then
         call compute_potential_allLevels(time)

         do lvl = -1, EFA_maxLvl
            call computeEfieldAtLvl(lvl)
         end do

         if (EFA_useElectrode) then
            voltage = 0.0d0
            voltageDiff = 0.0d0
            temp = 0.0d0
            voltageRMS = 0.0d0
            pointsUsed = 0
            do ix = 1, nElecPoints
               call EL_getSurfacePoint(ix, xyz)
               if (xyz(3) > 0.0d0) then
                  pointsUsed = pointsUsed + 1
                  maxlvl = getMaxLevelFrom(xyz, -1)
                  temp = get_value_at(grids(maxlvl), grids(maxlvl)%potential, xyz)
                  voltage = voltage + temp
                  if (abs(temp - EL_getVoltage(time)) > voltageDiff) then
                     voltageDiff = abs(temp - EL_getVoltage(time))
                     maxErrPos = xyz
                  end if
               end if
            end do
            ! print *, "Average voltage", voltage / pointsUsed, "max diff", voltageDiff
         end if

      end if

      do lvl = 0, EFA_maxLvl
         dataSize = size(grids(lvl)%Efield)
         call MPI_BCAST(grids(lvl)%Efield, dataSize, MPI_DOUBLE_PRECISION, &
              root, MPI_COMM_WORLD, ierr)
      end do

   end subroutine EFA_compute_Efield

   subroutine compute_potential_allLevels(time)
      double precision, intent(in) :: time
      integer :: lvl

      call setDirichletBndCnd(time, lvl = -1)
      call computePotentialAtLevel(lvl = -1)

      call setBndCndFromParentInterp(lvl = 0)
      call computePotentialAtLevel(lvl = 0)

      !       print *, "First two levels are computed"

      do lvl = 1, EFA_maxLvl
         call setBndCndFromParent(lvl)
         call computePotentialAtLevel(lvl)
      end do


      if (EFA_doCorrect) then

         ! Upward sweep that propagates correction charge layers down
         do lvl = EFA_maxLvl-1, 0, -1
            call addSingleLayerFromChild(lvl)
            call computePotentialAtLevel(lvl)
         end do

         ! Downward sweep that uses the corrected boundary conditions
         do lvl = 1, EFA_maxLvl
            call setBndCndFromParent(lvl)
            call computePotentialAtLevel(lvl)
         end do

      end if

   end subroutine compute_potential_allLevels

   subroutine addSingleLayerFromChild(lvl)
      integer, intent(in) :: lvl

      integer :: i, j, k
      integer :: Nx, Ny, Nz
      integer :: cIx, cIy, cIz
      integer :: ciMin, ciMax, cjMin, cjMax, ckMin, ckMax

      double precision :: invDx2, invDy2, invDz2
      double precision :: temp, tempSum

      ciMin = grids(lvl)%child_ixMin(1)
      cjMin = grids(lvl)%child_ixMin(2)
      ckMin = grids(lvl)%child_ixMin(3)

      ciMax = grids(lvl)%child_ixMax(1)
      cjMax = grids(lvl)%child_ixMax(2)
      ckMax = grids(lvl)%child_ixMax(3)

      invDx2 = grids(lvl)%invDelta(1)**2
      invDy2 = grids(lvl)%invDelta(2)**2
      invDz2 = grids(lvl)%invDelta(3)**2

      Nx = grids(lvl)%nPoints(1)
      Ny = grids(lvl)%nPoints(2)
      Nz = grids(lvl)%nPoints(3)

      tempSum = 0.0d0

      ! Compute Laplacian (u_coarse - u_composite) on the boundary which is the source term
      ! for the correction. Only the component normal to the refinement boundary is non-zero.

      ! Do faces in z-direction
      do j = cjMin+1, cjMax-1
         do i = ciMin+1, ciMax-1

            cIx = 2 * (i - ciMin) + 1
            cIy = 2 * (j - cjMin) + 1
            cIz = 3

            temp = (grids(lvl)%potential(i, j, ckMin+1) - grids(lvl+1)%potential(cIx, cIy, cIz)) * invDz2
            grids(lvl)%sourceTerm(i, j, ckMin+1) = grids(lvl)%sourceTerm(i, j, ckMin+1) + temp
            tempSum = tempSum + temp

            cIz = grids(lvl+1)%nPoints(3)-2

            temp = (grids(lvl)%potential(i, j, ckMax-1) - grids(lvl+1)%potential(cIx, cIy, cIz)) * invDz2
            grids(lvl)%sourceTerm(i, j, ckMax-1) = grids(lvl)%sourceTerm(i, j, ckMax-1) + temp
            tempSum = tempSum + temp
         end do
      end do

      ! Do faces in y-direction
      do k = ckMin+1, ckMax-1
         do i = ciMin+1, ciMax-1

            cIx = 2 * (i - ciMin) + 1
            cIy = 3
            cIz = 2 * (k - ckMin) + 1

            temp = (grids(lvl)%potential(i, cjMin+1, k) - grids(lvl+1)%potential(cIx, cIy, cIz)) * invDy2
            grids(lvl)%sourceTerm(i, cjMin+1, k) = grids(lvl)%sourceTerm(i, cjMin+1, k) + temp
            tempSum = tempSum + temp

            cIy = grids(lvl+1)%nPoints(2)-2

            temp = (grids(lvl)%potential(i, cjMax-1, k) - grids(lvl+1)%potential(cIx, cIy, cIz)) * invDy2
            grids(lvl)%sourceTerm(i, cjMax-1, k) = grids(lvl)%sourceTerm(i, cjMax-1, k) + temp
            tempSum = tempSum + temp
         end do
      end do

      ! Do faces in x-direction
      do k = ckMin+1, ckMax-1
         do j = cjMin+1, cjMax-1

            cIx = 3
            cIy = 2 * (j - cjMin) + 1
            cIz = 2 * (k - ckMin) + 1

            temp = (grids(lvl)%potential(ciMin+1, j, k) - grids(lvl+1)%potential(cIx, cIy, cIz)) * invDx2
            grids(lvl)%sourceTerm(ciMin+1, j, k) = grids(lvl)%sourceTerm(ciMin+1, j, k) + temp
            tempSum = tempSum + temp

            cIx = grids(lvl+1)%nPoints(1)-2

            temp = (grids(lvl)%potential(ciMax-1, j, k) - grids(lvl+1)%potential(cIx, cIy, cIz)) * invDx2
            grids(lvl)%sourceTerm(ciMax-1, j, k) = grids(lvl)%sourceTerm(ciMax-1, j, k) + temp
            tempSum = tempSum + temp
         end do
      end do

      tempSum = tempSum / product(dble(grids(lvl)%nPoints-2))
      grids(lvl)%sourceTerm(2:Nx-1, 2:Ny-1, 2:Nz-1) = grids(lvl)%sourceTerm(2:Nx-1, 2:Ny-1, 2:Nz-1) - tempSum

   end subroutine addSingleLayerFromChild

   subroutine computePotentialAtLevel(lvl)
      integer, intent(in) :: lvl

      interface
         subroutine hw3crt(xMin, xMax, xPanels, xBDCND, xBDCmin, xBDCmax, &
              yMin, yMax, yPanels, yBDCND, yBDCmin, yBDCmax, &
              zMin, zMax, zPanels, zBDCND, zBDCmin, zBDCmax, &
              lambda, sourceFirstDim, sourceSecondDim, source, &
              pertrb, ierror, workspace)

            integer          :: xPanels, xBDCND, yPanels, yBDCND, zPanels, zBDCND
            integer          :: sourceFirstDim, sourceSecondDim, ierror
            double precision :: xMin, xMax, xBDCmin(*), xBDCmax(*)
            double precision :: yMin, yMax, yBDCmin(*), yBDCmax(*)
            double precision :: zMin, zMax, zBDCmin(*), zBDCmax(*)
            double precision :: lambda, source(sourceFirstDim, sourceSecondDim, *)
            double precision :: pertrb, workspace(*)
         end subroutine hw3crt
      end interface

      integer :: ierror, Nx, Ny, Nz
      double precision :: dummy(1)

      call ensureWorkspace(lvl)

      Nx = grids(lvl)%nPoints(1)
      Ny = grids(lvl)%nPoints(2)
      Nz = grids(lvl)%nPoints(3)

      ! The input array holds the boundary condition on the sides and the source term in its interior
      ! Note: source term is minus the charge density. electronCharge has a sign!
      grids(lvl)%potential(2:Nx-1, 2:Ny-1, 2:Nz-1) = grids(lvl)%sourceTerm(2:Nx-1, 2:Ny-1, 2:Nz-1)

      call hw3crt(grids(lvl)%posMin(1), grids(lvl)%posMax(1), Nx-1, 1, dummy, dummy, &
           grids(lvl)%posMin(2), grids(lvl)%posMax(2), Ny-1, 1, dummy, dummy, &
           grids(lvl)%posMin(3), grids(lvl)%posMax(3), Nz-1, 1, dummy, dummy, &
           0.0d0, Nx, Ny, grids(lvl)%potential, dummy(1), ierror, workspace)

      if (ierror /= 0) then
         print *, "computePotentialAtLevel: ierror = ", ierror
         stop
      end if

   end subroutine computePotentialAtLevel

   !> Collect electron and ion densities.
   !! Important: for the particle code, elecDens needs to be available on every task.
   !! Use the function EFA_shareElecDens() for that after this routine is called
   subroutine EFA_collectDensitiesAtRoot(lvl, myrank, root)
      integer, intent(in)           :: lvl, myrank, root

      integer                       :: ierr, dataSize
      dataSize = product(grids(lvl)%nPoints)

      if (myrank == root) then
         call MPI_REDUCE(MPI_IN_PLACE, grids(lvl)%posIonDens, dataSize, &
              MPI_DOUBLE_PRECISION, MPI_SUM, root, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(MPI_IN_PLACE, grids(lvl)%negIonDens, dataSize, &
              MPI_DOUBLE_PRECISION, MPI_SUM, root, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(MPI_IN_PLACE, grids(lvl)%O2MinDens, dataSize, &
              MPI_DOUBLE_PRECISION, MPI_SUM, root, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(MPI_IN_PLACE, grids(lvl)%elecDens, dataSize, &
              MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      else
         call MPI_REDUCE(grids(lvl)%posIonDens, grids(lvl)%posIonDens, dataSize, &
              MPI_DOUBLE_PRECISION, MPI_SUM, root, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(grids(lvl)%negIonDens, grids(lvl)%negIonDens, dataSize, &
              MPI_DOUBLE_PRECISION, MPI_SUM, root, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(grids(lvl)%O2MinDens, grids(lvl)%O2MinDens, dataSize, &
              MPI_DOUBLE_PRECISION, MPI_SUM, root, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(grids(lvl)%elecDens, grids(lvl)%elecDens, dataSize, &
              MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

         ! We collect ions at the root task, remove them locally
         grids(lvl)%posIonDens = 0.0d0
         grids(lvl)%negIonDens = 0.0d0
         grids(lvl)%O2MinDens  = 0.0d0
         grids(lvl)%elecDens   = 0.0d0
      end if

   end subroutine EFA_collectDensitiesAtRoot

   subroutine EFA_shareElecDens(myrank, root)
      integer, intent(in)  :: myrank, root
      integer              :: lvl, dataSize, ierr

      do lvl = 0, EFA_maxLvl
         dataSize = product(grids(lvl)%nPoints)
         call MPI_BCAST(grids(lvl)%elecDens, dataSize, MPI_DOUBLE_PRECISION, &
              root, MPI_COMM_WORLD, ierr)
      end do
   end subroutine EFA_shareElecDens

   subroutine collectExcDensAtRoot(lvl, myrank, root)
      integer, intent(in)           :: lvl, myrank, root

      integer                       :: ierr, dataSize

      dataSize = product(grids(lvl)%nPoints)

      if (myrank == root) then
         call MPI_REDUCE(MPI_IN_PLACE, grids(lvl)%excDens, dataSize, &
              MPI_DOUBLE_PRECISION, MPI_SUM, root, MPI_COMM_WORLD, ierr)

      else
         call MPI_REDUCE(grids(lvl)%excDens, grids(lvl)%excDens, dataSize, &
              MPI_DOUBLE_PRECISION, MPI_SUM, root, MPI_COMM_WORLD, ierr)

         ! We collect excited species at the root task, remove them locally
         grids(lvl)%excDens = 0.0d0
      end if

   end subroutine collectExcDensAtRoot

   subroutine EFA_allCollectTemp()
      integer              :: ierr, dataSize, lvl

      do lvl = 0, EFA_maxLvl
         dataSize = product(grids(lvl)%nPoints)

         call MPI_ALLREDUCE(MPI_IN_PLACE, grids(lvl)%temp, dataSize, &
              MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      end do

   end subroutine EFA_allCollectTemp

   subroutine EFA_addO2MinDens(value)
      double precision, intent(in) :: value
      integer :: lvl
      do lvl = -1, EFA_maxLvl
         grids(lvl)%O2MinDens = grids(lvl)%O2MinDens + value
      end do
   end subroutine EFA_addO2MinDens

   subroutine EFA_addPosIonDens(value)
      double precision, intent(in) :: value
      integer :: lvl
      do lvl = -1, EFA_maxLvl
         grids(lvl)%posIonDens = grids(lvl)%posIonDens + value
      end do
   end subroutine EFA_addPosIonDens

   subroutine EFA_setNDecays(lvl, lossRatio, nDecays, myrank, root)
      integer, intent(in)             :: lvl, myrank, root
      double precision, intent(in)    :: lossRatio
      double precision, intent(out)   :: nDecays(:,:,:)

      integer :: i, j, k, nIx(3), ix(3)
      double precision :: cellVolume

      nIx = grids(lvl)%nPoints

      if (myrank == root ) then
         if(any(nIx /= shape(nDecays))) then
            print *, "EFA_setNDecays: argument nDecays has wrong shape"
            stop
         end if
      end if

      cellVolume = product(grids(lvl)%delta)

      call collectExcDensAtRoot(lvl, myrank, root)

      if (myrank == root) then

         do k = 1, nIx(3)
            do j = 1, nIx(2)
               do i = 1, nIx(1)
                  ix = (/i, j, k/)

                  ! If not in child region, set nDecays. In child region it is zero.
                  if (any(ix < grids(lvl)%child_ixMin) .or. any(ix > grids(lvl)%child_ixMax)) then
                     nDecays(i,j,k) = lossRatio * grids(lvl)%excDens(i,j,k) * cellVolume
                  else
                     nDecays(i,j,k) = 0.0d0
                  end if

                  grids(lvl)%excDens(i,j,k) = grids(lvl)%excDens(i,j,k) * (1.0d0 - lossRatio)

               end do
            end do
         end do

      end if

   end subroutine EFA_setNDecays

   ! Before this routine is called, the densities should be collected at the root.
   ! It should only be called by the root.
   subroutine EFA_setNegIonLoss(lvl, dt, negIonLoss)
      integer, intent(in)             :: lvl
      double precision, intent(in)    :: dt
      double precision, intent(out)   :: negIonLoss(:,:,:)

      integer :: i, j, k, nIx(3), ix(3)
      double precision :: cellVolume, O2_bgdens, N2_bgdens

      nIx = grids(lvl)%nPoints

      if (any(nIx /= shape(negIonLoss))) then
         print *, "EFA_setNegIonLoss: argument negIonLoss has wrong shape"
         stop
      end if

      cellVolume = product(grids(lvl)%delta)
      O2_bgdens = GAS_numberDensity * getGasFraction("O2")
      N2_bgdens = GAS_numberDensity * getGasFraction("N2")
      do k = 1, nIx(3)
         do j = 1, nIx(2)
            do i = 1, nIx(1)
               ix = (/i, j, k/)
               negIonLoss(i,j,k) = getNegIonLoss(grids(lvl)%O2MinDens(i,j,k), &
                    O2_bgdens, N2_bgdens, grids(lvl)%Efield(:,i,j,k), dt)
               grids(lvl)%O2MinDens(i,j,k) = grids(lvl)%O2MinDens(i,j,k) - negIonLoss(i,j,k)

               ! In child region we should use values from the child, so set to zero.
               if (grids(lvl)%hasChild .and. &
                    all(ix >= grids(lvl)%child_ixMin) .and. all(ix <= grids(lvl)%child_ixMax)) then
                  negIonLoss(i,j,k) = 0.0d0
               end if
            end do
         end do
      end do

      negIonLoss = negIonLoss * cellVolume
   end subroutine EFA_setNegIonLoss

   double precision function getNegIonLoss(O2min_dens, O2_bgdens, N2_bgdens, Efield, dt)
      double precision, intent(in) :: O2min_dens, O2_bgdens, N2_bgdens, Efield(3), dt

      double precision :: T_eff, EfieldNorm

      EfieldNorm = sqrt(sum(Efield**2))
      T_eff = GAS_temperature + N2_mass / (3.0d0 * BoltzmannConstant) * (EfieldNorm * 2.0d-4)**2

      getNegIonLoss = 1.9d-12 * sqrt(T_eff/3.0d2) * exp(-4.990e3/T_eff) * N2_bgdens * 1.0d-6 + &
           2.7d-10 * sqrt(T_eff/3.0d2) * exp(-5.590e3/T_eff) * O2_bgdens * 1.0d-6
      getNegIonLoss = getNegIonLoss * O2min_dens * dt
   end function getNegIonLoss

   subroutine setSourceTermFromParent(lvl)
      integer, intent(in)           :: lvl

      integer :: i, j, k, nIx(3)
      double precision :: xyz(3)

      if (lvl < 0 .or. lvl > EFA_maxLvl) then
         print *, "setSourceTermFromParent error: level does not exist", lvl
         stop
      end if

      nIx = grids(lvl)%nPoints

      do k = 1, nIx(3)
         do j = 1, nIx(2)
            do i = 1, nIx(1)
               xyz = (/i-1, j-1, k-1/) * grids(lvl)%delta + grids(lvl)%posMin
               grids(lvl)%sourceTerm(i, j, k) = get_value_at(grids(lvl-1), grids(lvl-1)%sourceTerm, xyz)
            end do
         end do
      end do

   end subroutine setSourceTermFromParent

   subroutine setSourceTermFromDens(lvl)
      integer, intent(in) :: lvl
      grids(lvl)%sourceTerm = elecChargeOverEps0 / epsilonRelative * &
           (grids(lvl)%posIonDens - grids(lvl)%negIonDens - grids(lvl)%O2MinDens - grids(lvl)%elecDens)
   end subroutine setSourceTermFromDens

   subroutine setDirichletBndCnd(time, lvl)
      use module_electrode
      double precision, intent(in) :: time
      integer, intent(in) :: lvl

      integer :: i, j, k
      integer :: Nx, Ny, Nz
      double precision :: xyz(3), posMin(3), delta(3), efield, voltage, temp, temp2
      double precision :: elec_pos(3), decay_range, elec_radius, elec_length

      delta    = grids(lvl)%delta
      posMin   = grids(lvl)%posMin

      Nx = grids(lvl)%nPoints(1)
      Ny = grids(lvl)%nPoints(2)
      Nz = grids(lvl)%nPoints(3)


      if (EFA_useElectrode) then
         voltage = EL_getVoltage(time)
      else
         call linearInterpolateList(EFA_efield_times, EFA_efield_values, time, efield)
         voltage = -efield * (grids(lvl)%posMax(3) - grids(lvl)%posMin(3))
         print *, "time/efield", time, efield
      end if

      if (.not. EFA_useElectrode .or. EFA_bcType == 1) then
         ! Potential linearly de/increasing along the sides, from voltage to zero

         do k = 1, Nz
            do j = 1, Ny
               do i = 1, Nx, Nx-1
                  grids(lvl)%potential(i,j,k) = (1.0d0 - (k-1)/dble(Nz-1)) * voltage
               end do
            end do
         end do

         do k = 1, Nz
            do j = 1, Ny, Ny-1
               do i = 1, Nx
                  grids(lvl)%potential(i,j,k) = (1.0d0 - (k-1)/dble(Nz-1)) * voltage
               end do
            end do
         end do

         do k = 1, Nz, Nz-1
            do j = 1, Ny
               do i = 1, Nx
                  grids(lvl)%potential(i,j,k) = (1.0d0 - (k-1)/dble(Nz-1)) * voltage
               end do
            end do
         end do

      else if (EFA_useElectrode .and. EFA_bcType == 2) then
         ! Potential zero at the sides.
         elec_radius = EL_getRadius(posMin(3))
         elec_length = 0.5d0 * delta(3) * (Nz-1)
         decay_range = 0.25d0 * min(delta(1) * (Nx-1), delta(2) * (Ny-1))
         call EL_getBottomPos(elec_pos)

         do k = 1, Nz
            do j = 1, Ny
               do i = 1, Nx, Nx-1
                  grids(lvl)%potential(i,j,k) = 0.0d0
               end do
            end do
         end do

         do k = 1, Nz
            do j = 1, Ny, Ny-1
               do i = 1, Nx
                  grids(lvl)%potential(i,j,k) = 0.0d0
               end do
            end do
         end do

         do j = 1, Ny
            do i = 1, Nx
               temp = twoNorm(posMin(1:2) +  (/i-1, j-1/) * delta(1:2) - elec_pos(1:2) )

               if (temp <= elec_radius) then
                  grids(lvl)%potential(i,j,1) = voltage
               else if (temp <= decay_range) then
                  temp2 = voltage / log( (elec_length + sqrt(elec_length**2 + elec_radius**2)) / &
                       & (-elec_length + sqrt(elec_length**2 + elec_radius**2)) )
                  grids(lvl)%potential(i,j,1) = temp2 * log( (elec_length + sqrt(elec_length**2 + temp**2)) / &
                       & (-elec_length + sqrt(elec_length**2 + temp**2)) )
               else
                  grids(lvl)%potential(i,j,1) = 0.0d0
               end if

               grids(lvl)%potential(i,j,Nz) = 0.0d0
            end do
         end do

      else
         print *, "setDirichletBndCnd error, no correct boundary condition used"
         stop
      end if

   end subroutine setDirichletBndCnd

   subroutine setBndCndFromParentInterp(lvl)
      integer, intent(in) :: lvl

      integer :: i, j, k
      integer :: Nx, Ny, Nz
      double precision :: xyz(3), delta(3), posmin(3)

      Nx       = grids(lvl)%nPoints(1)
      Ny       = grids(lvl)%nPoints(2)
      Nz       = grids(lvl)%nPoints(3)

      delta    = grids(lvl)%delta
      posmin   = grids(lvl)%posMin

      do k = 1, Nz, Nz-1
         do j = 1, Ny
            do i = 1, Nx
               xyz = (/i-1, j-1, k-1/) * delta + posMin
               grids(lvl)%potential(i,j,k) = get_value_at(grids(lvl-1), grids(lvl-1)%potential, xyz)
            end do
         end do
      end do

      do k = 1, Nz
         do j = 1, Ny, Ny-1
            do i = 1, Nx
               xyz = (/i-1, j-1, k-1/) * delta + posMin
               grids(lvl)%potential(i,j,k) = get_value_at(grids(lvl-1), grids(lvl-1)%potential, xyz)
            end do
         end do
      end do

      do k = 1, Nz
         do j = 1, Ny
            do i = 1, Nx, Nx-1
               xyz = (/i-1, j-1, k-1/) * delta + posMin
               grids(lvl)%potential(i,j,k) = get_value_at(grids(lvl-1), grids(lvl-1)%potential, xyz)
            end do
         end do
      end do

   end subroutine setBndCndFromParentInterp

   subroutine setBndCndFromParent(lvl)
      integer, intent(in) :: lvl

      integer :: i, j, k
      integer :: Nx, Ny, Nz
      integer :: pIx, pIy, pIz
      integer :: pIx0, pIy0, pIz0
      integer :: pIx1, pIy1, pIz1

      double precision, allocatable :: xplane0(:,:), yplane0(:,:), zplane0(:,:)
      double precision, allocatable :: xplane1(:,:), yplane1(:,:), zplane1(:,:)

      double precision :: temp

      if (lvl <= 0) then
         print *, "setBndCndFromParent: lvl should be > 0"
         stop
      end if

      pIx0 = grids(lvl-1)%child_ixMin(1)
      pIy0 = grids(lvl-1)%child_ixMin(2)
      pIz0 = grids(lvl-1)%child_ixMin(3)

      pIx1 = grids(lvl-1)%child_ixMax(1)
      pIy1 = grids(lvl-1)%child_ixMax(2)
      pIz1 = grids(lvl-1)%child_ixMax(3)

      Nx = grids(lvl)%nPoints(1)
      Ny = grids(lvl)%nPoints(2)
      Nz = grids(lvl)%nPoints(3)

      ! These will hold a temporary copy of the boundary planes. They are somewhat larger,
      ! to simplify the interpolation code.
      allocate( xplane0(-1:Ny+2, -1:Nz+2) )
      allocate( yplane0(-1:Nx+2, -1:Nz+2) )
      allocate( zplane0(-1:Nx+2, -1:Ny+2) )

      allocate( xplane1(-1:Ny+2, -1:Nz+2) )
      allocate( yplane1(-1:Nx+2, -1:Nz+2) )
      allocate( zplane1(-1:Nx+2, -1:Ny+2) )

      ! We can get the boundary conditions from the parent. Because we go from N -> 2*N-1 points
      ! and have dx -> 0.5*dx we can get half the points directly from the parent (they overlap).
      ! For the other we interpolate, using 1/16 * (f_-3/2 + 9 * f_-1/2 + 9 * f_1/2 - f_3/2) to
      ! get f_0, which is 4th order accurate.

      ! Set points that are overlapping on x-planes
      do k = -1, Nz+2, 2
         do j = -1, Ny+2, 2
            pIy = pIy0 + (j-1)/2
            pIz = pIz0 + (k-1)/2
            xplane0(j,k) = grids(lvl-1)%potential(pIx0, pIy, pIz)
            xplane1(j,k) = grids(lvl-1)%potential(pIx1, pIy, pIz)
         end do
      end do

      ! Set points that are overlapping on y-planes and z-planes
      do i = -1, Nx+2, 2
         pIx = pIx0 + (i-1)/2

         do k = -1, Nz+2, 2
            pIz = pIz0 + (k-1)/2
            yplane0(i,k) = grids(lvl-1)%potential(pIx, pIy0, pIz)
            yplane1(i,k) = grids(lvl-1)%potential(pIx, pIy1, pIz)
         end do

         do j = -1, Ny+2, 2
            pIy = pIy0 + (j-1)/2
            zplane0(i,j) = grids(lvl-1)%potential(pIx, pIy, pIz0)
            zplane1(i,j) = grids(lvl-1)%potential(pIx, pIy, pIz1)
         end do
      end do


      ! Set non-overlapping points on the x-planes
      do k = -1, Nz+2, 2
         ! Set points that are right between two parents points using a 4th order method
         do j = 2, Ny-1, 2
            xplane0(j,k) = 0.5625d0 * (xplane0(j-1,k) + xplane0(j+1,k)) &
                 - 0.0625d0 * (xplane0(j-3,k) + xplane0(j+3,k))
            xplane1(j,k) = 0.5625d0 * (xplane1(j-1,k) + xplane1(j+1,k)) &
                 - 0.0625d0 * (xplane1(j-3,k) + xplane1(j+3,k))
         end do
      end do

      do k = 2, Nz-1, 2
         ! Use the values computed in the loop above to set values in between 4 parent points
         do j = 2, Ny-1, 2
            xplane0(j,k) = 0.5625d0 * (xplane0(j,k-1) + xplane0(j,k+1)) &
                 - 0.0625d0 * (xplane0(j,k-3) + xplane0(j,k+3))
            xplane1(j,k) = 0.5625d0 * (xplane1(j,k-1) + xplane1(j,k+1)) &
                 - 0.0625d0 * (xplane1(j,k-3) + xplane1(j,k+3))
         end do
      end do

      ! Set points that are right between two parents points using a 4th order method
      do k = 2, Nz-1, 2
         do j = 1, Ny, 2
            xplane0(j,k) = 0.5625d0 * (xplane0(j,k-1) + xplane0(j,k+1)) &
                 - 0.0625d0 * (xplane0(j,k-3) + xplane0(j,k+3))
            xplane1(j,k) = 0.5625d0 * (xplane1(j,k-1) + xplane1(j,k+1)) &
                 - 0.0625d0 * (xplane1(j,k-3) + xplane1(j,k+3))
         end do
      end do


      do i = 2, Nx-1, 2

         ! Set non-overlapping points on the y-planes
         ! Set points that are right between two parents points using a 4th order method
         do k = -1, Nz+2, 2
            yplane0(i,k) = 0.5625d0 * (yplane0(i-1,k) + yplane0(i+1,k)) &
                 - 0.0625d0 * (yplane0(i-3,k) + yplane0(i+3,k))
            yplane1(i,k) = 0.5625d0 * (yplane1(i-1,k) + yplane1(i+1,k)) &
                 - 0.0625d0 * (yplane1(i-3,k) + yplane1(i+3,k))
         end do

         ! Use the values computed in the loop above to set values in between 4 parent points
         do k = 2, Nz-1, 2
            yplane0(i,k) = 0.5625d0 * (yplane0(i,k-1) + yplane0(i,k+1)) &
                 - 0.0625d0 * (yplane0(i,k-3) + yplane0(i,k+3))
            yplane1(i,k) = 0.5625d0 * (yplane1(i,k-1) + yplane1(i,k+1)) &
                 - 0.0625d0 * (yplane1(i,k-3) + yplane1(i,k+3))
         end do

         ! Set non-overlapping points on the z-planes
         ! Set points that are right between two parents points using a 4th order method
         do j = -1, Ny+2, 2
            zplane0(i,j) = 0.5625d0 * (zplane0(i-1,j) + zplane0(i+1,j)) &
                 - 0.0625d0 * (zplane0(i-3,j) + zplane0(i+3,j))
            zplane1(i,j) = 0.5625d0 * (zplane1(i-1,j) + zplane1(i+1,j)) &
                 - 0.0625d0 * (zplane1(i-3,j) + zplane1(i+3,j))
         end do

         ! Use the values computed in the loop above to set values in between 4 parent points
         do j = 2, Ny-1, 2
            zplane0(i,j) = 0.5625d0 * (zplane0(i,j-1) + zplane0(i,j+1)) &
                 - 0.0625d0 * (zplane0(i,j-3) + zplane0(i,j+3))
            zplane1(i,j) = 0.5625d0 * (zplane1(i,j-1) + zplane1(i,j+1)) &
                 - 0.0625d0 * (zplane1(i,j-3) + zplane1(i,j+3))
         end do
      end do

      do k = 2, Nz-1, 2
         do i = 1, Nx, 2
            yplane0(i,k) = 0.5625d0 * (yplane0(i,k-1) + yplane0(i,k+1)) &
                 - 0.0625d0 * (yplane0(i,k-3) + yplane0(i,k+3))
            yplane1(i,k) = 0.5625d0 * (yplane1(i,k-1) + yplane1(i,k+1)) &
                 - 0.0625d0 * (yplane1(i,k-3) + yplane1(i,k+3))
         end do
      end do

      do j = 2, Ny-1, 2
         do i = 1, Nx, 2
            zplane0(i,j) = 0.5625d0 * (zplane0(i,j-1) + zplane0(i,j+1)) &
                 - 0.0625d0 * (zplane0(i,j-3) + zplane0(i,j+3))
            zplane1(i,j) = 0.5625d0 * (zplane1(i,j-1) + zplane1(i,j+1)) &
                 - 0.0625d0 * (zplane1(i,j-3) + zplane1(i,j+3))
         end do
      end do


      ! Now we have all the needed boundary values
      grids(lvl)%potential(1,:,:)   = xplane0(1:Ny, 1:Nz)
      grids(lvl)%potential(Nx,:,:)  = xplane1(1:Ny, 1:Nz)

      grids(lvl)%potential(:,1,:)   = yplane0(1:Nx, 1:Nz)
      grids(lvl)%potential(:,Ny,:)  = yplane1(1:Nx, 1:Nz)

      grids(lvl)%potential(:,:,1)   = zplane0(1:Nx, 1:Ny)
      grids(lvl)%potential(:,:,Nz)  = zplane1(1:Nx, 1:Ny)

      deallocate( xplane0 )
      deallocate( yplane0 )
      deallocate( zplane0 )

      deallocate( xplane1 )
      deallocate( yplane1 )
      deallocate( zplane1 )

   end subroutine setBndCndFromParent

   double precision function get_value_at(myGrid, array3D, xyz)
      type(amrGrid), intent(inout)     :: myGrid
      double precision, intent(in)     :: xyz(3)
      double precision, intent(inout)  :: array3D(:,:,:)

      integer                          :: i, ixs(3)
      double precision                 :: delta(3), invDelta(3)
      double precision                 :: rRel(3), rmin(3), rmax(3), tmp, tmp2

      delta       = myGrid%delta
      invDelta    = myGrid%invDelta
      rRel        = xyz - myGrid%posMin
      ixs         = floor(rRel * invDelta) + 1

      ! Bounds checking
      do i = 1, 3
         ixs(i)   = max(1, min(ixs(i), myGrid%nPoints(i)-1))
      end do

      rmax        = ixs * delta
      rmin        = (ixs-1) * delta

      ! Perform trilinear interpolation of the sourceTerm at position xyz from
      ! the 8 corners of the cube in which rRel lies. This is done by first doing
      ! two bilinear interpolations on the y,z-planes and then linearly interpolating
      ! in the x-direction.
      call bilinearInterp(tmp, rRel(2), rRel(3), rmin(2), rmin(3), rmax(2), rmax(3), &
           array3D(ixs(1), ixs(2),   ixs(3)), array3D(ixs(1), ixs(2),   ixs(3)+1), &
           array3D(ixs(1), ixs(2)+1, ixs(3)), array3D(ixs(1), ixs(2)+1, ixs(3)+1) )

      call bilinearInterp(tmp2, rRel(2), rRel(3), rmin(2), rmin(3), rmax(2), rmax(3), &
           array3D(ixs(1)+1, ixs(2),   ixs(3)), array3D(ixs(1)+1, ixs(2),   ixs(3)+1), &
           array3D(ixs(1)+1, ixs(2)+1, ixs(3)), array3D(ixs(1)+1, ixs(2)+1, ixs(3)+1) )

      call linearInterp(get_value_at, rRel(1), rmin(1), rmax(1), tmp, tmp2)

   end function get_value_at

   subroutine computeEfieldAtLvl(lvl)
      integer, intent(in) :: lvl

      integer :: i, j, k
      integer :: Nx, Ny, Nz
      double precision :: delta(3), halfInvDelta(3), Efield(3), xyz(3), posMin(3)

      Nx = grids(lvl)%nPoints(1)
      Ny = grids(lvl)%nPoints(2)
      Nz = grids(lvl)%nPoints(3)

      posMin         = grids(lvl)%posMin
      delta          = grids(lvl)%delta
      halfInvDelta   = 0.5d0 * grids(lvl)%invDelta

      ! Do the interior part of the domain
      do k = 2, Nz-1
         do j = 2, Ny-1
            do i = 2, Nx-1
               Efield(1) = grids(lvl)%potential(i-1, j, k) - grids(lvl)%potential(i+1, j, k)
               Efield(2) = grids(lvl)%potential(i, j-1, k) - grids(lvl)%potential(i, j+1, k)
               Efield(3) = grids(lvl)%potential(i, j, k-1) - grids(lvl)%potential(i, j, k+1)

               grids(lvl)%Efield(:,i,j,k) = halfInvDelta * Efield
            end do
         end do
      end do

      ! Do the boundaries. Either interpolate from a coarser level or interpolate from the inside of
      ! the domain.

      if (lvl >= 0) then

         ! Interpolate electric field from parent grid

         do k = 1, Nz, Nz-1
            do j = 1, Ny
               do i = 1, Nx
                  xyz = (/i-1, j-1, k-1/) * delta + posMin
                  call get_Efield_at_lvl(lvl-1, xyz, grids(lvl)%Efield(:,i, j, k))
               end do
            end do
         end do

         do k = 1, Nz
            do j = 1, Ny, Ny-1
               do i = 1, Nx
                  xyz = (/i-1, j-1, k-1/) * delta + posMin
                  call get_Efield_at_lvl(lvl-1, xyz, grids(lvl)%Efield(:,i, j, k))
               end do
            end do
         end do

         do k = 1, Nz
            do j = 1, Ny
               do i = 1, Nx, Nx-1
                  xyz = (/i-1, j-1, k-1/) * delta + posMin
                  call get_Efield_at_lvl(lvl-1, xyz, grids(lvl)%Efield(:,i, j, k))
               end do
            end do
         end do

      else ! lvl < 0

         ! Interpolate from the interior region
         do j = 1, Ny
            do i = 1, Nx
               grids(lvl)%Efield(:,i, j, 1)  = 2.0d0 * grids(lvl)%Efield(:,i, j, 2) &
                    - grids(lvl)%Efield(:,i, j, 3)
               grids(lvl)%Efield(:,i, j, Nz) = 2.0d0 * grids(lvl)%Efield(:,i, j, Nz-1) &
                    - grids(lvl)%Efield(:,i, j, Nz-2)
            end do
         end do

         do k = 1, Nz
            do i = 1, Nx
               grids(lvl)%Efield(:,i, 1, k)  = 2.0d0 * grids(lvl)%Efield(:,i, 2, k) &
                    - grids(lvl)%Efield(:,i, 3, k)
               grids(lvl)%Efield(:,i, Ny, k) = 2.0d0 * grids(lvl)%Efield(:,i, Ny-1, k) &
                    - grids(lvl)%Efield(:,i, Ny-2, k)
            end do
         end do

         do k = 1, Nz
            do j = 1, Ny
               grids(lvl)%Efield(:,1, j, k)  = 2.0d0 * grids(lvl)%Efield(:,2, j, k) &
                    - grids(lvl)%Efield(:,3, j, k)
               grids(lvl)%Efield(:,Nx, j, k) = 2.0d0 * grids(lvl)%Efield(:,Nx-1, j, k) &
                    - grids(lvl)%Efield(:,Nx-2, j, k)
            end do
         end do

      end if

   end subroutine computeEfieldAtLvl

   subroutine get_Efield_at_lvl(lvl, xyz, Efield)
      integer, intent(in)              :: lvl
      double precision, intent(in)     :: xyz(3)
      double precision, intent(inout)  :: Efield(3)

      integer                          :: i, ixs(3)
      double precision                 :: delta(3), invDelta(3)
      double precision                 :: rRel(3), rmin(3), rmax(3), Efield1(3), Efield2(3)

      delta       = grids(lvl)%delta
      invDelta    = grids(lvl)%invDelta
      rRel        = xyz - grids(lvl)%posMin

      ixs         = floor(rRel * invDelta) + 1

      ! Bounds checkings
      do i = 1, 3
         ixs(i) = max(1, min(ixs(i), grids(lvl)%nPoints(i)-1))
      end do

      rmax  = ixs * delta
      rmin  = (ixs-1) * delta

      ! Perform trilinear interpolation of the potential at position xyz from
      ! the 8 corners of the cube in which rRel lies. This is done by first doing
      ! two bilinear interpolations on the y,z-planes and then linearly interpolating
      ! in the x-direction.
      call bilinearInterp3(Efield1, rRel(2), rRel(3), rmin(2), rmin(3), rmax(2), rmax(3), &
           grids(lvl)%Efield(:,ixs(1), ixs(2),   ixs(3)), grids(lvl)%Efield(:,ixs(1), ixs(2),   ixs(3)+1), &
           grids(lvl)%Efield(:,ixs(1), ixs(2)+1, ixs(3)), grids(lvl)%Efield(:,ixs(1), ixs(2)+1, ixs(3)+1) )

      call bilinearInterp3(Efield2, rRel(2), rRel(3), rmin(2), rmin(3), rmax(2), rmax(3), &
           grids(lvl)%Efield(:,ixs(1)+1, ixs(2),   ixs(3)), grids(lvl)%Efield(:,ixs(1)+1, ixs(2),   ixs(3)+1), &
           grids(lvl)%Efield(:,ixs(1)+1, ixs(2)+1, ixs(3)), grids(lvl)%Efield(:,ixs(1)+1, ixs(2)+1, ixs(3)+1) )

      call linearInterp3(Efield, rRel(1), rmin(1), rmax(1), Efield1, Efield2)

   end subroutine get_Efield_at_lvl

   double precision function EFA_get_potential(xyz)
      double precision, intent(in) :: xyz(3)
      integer :: lvl

      lvl               = getMaxLevelFromZero(xyz)
      EFA_get_potential = get_value_at(grids(lvl), grids(lvl)%potential, xyz)

   end function EFA_get_potential

   subroutine EFA_getEfield(xyz, Efield)
      double precision, intent(in)     :: xyz(3)
      double precision, intent(inout)  :: Efield(3)
      integer                          :: lvl

      lvl = getMaxLevelFromZero(xyz)
      call get_Efield_at_lvl(lvl, xyz, Efield)

   end subroutine EFA_getEfield

   integer function getMaxLevelFromZero(xyz)
      double precision, intent(in)  :: xyz(3)
      getMaxLevelFromZero = getMaxLevelFrom(xyz, 0)
   end function getMaxLevelFromZero

   integer function getMaxLevelFrom(xyz, startLvl)
      double precision, intent(in)  :: xyz(3)
      integer, intent(in)           :: startLvl
      integer                       :: lvl

      lvl = startLvl
      do while (grids(lvl)%hasChild)
         if (any(xyz <= grids(lvl+1)%posMin) .or. any(xyz >= grids(lvl+1)%posMax)) exit
         lvl = lvl + 1
      end do

      getMaxLevelFrom = lvl
   end function getMaxLevelFrom

   integer function getMaxLevelInsideAt(xyz)
      double precision, intent(in)  :: xyz(3)
      integer                       :: lvl

      lvl = 0
      do while (grids(lvl)%hasChild)
         if (any(xyz <= grids(lvl+1)%posMin + grids(lvl+1)%delta) .or. &
              & any(xyz >= grids(lvl+1)%posMax - grids(lvl+1)%delta)) exit
         lvl = lvl + 1
      end do

      getMaxLevelInsideAt = lvl
   end function getMaxLevelInsideAt

   subroutine setGridlevel(lvl, nPoints, posMin, posMax, hasChild, ixOffset)
      integer, intent(in)           :: lvl, nPoints(3)
      double precision, intent(in)  :: posMin(3), posMax(3)
      logical, intent(in), optional :: hasChild
      integer, intent(in), optional :: ixOffset

      grids(lvl)%hasChild     = .false.
      if (present(hasChild)) grids(lvl)%hasChild = hasChild

      grids(lvl)%cellIxOffset = 0
      if (present(ixOffset)) grids(lvl)%cellIxOffset = ixOffset

      grids(lvl)%nPoints      = nPoints
      grids(lvl)%posMin       = posMin
      grids(lvl)%posMax       = posMax
      grids(lvl)%delta        = (posMax - posMin) / dble(nPoints - 1)
      grids(lvl)%invDelta     = 1.0d0 / grids(lvl)%delta

      call ensureGridIsAllocated(grids(lvl))

      EFA_maxLvl = max(EFA_maxLvl, lvl)

   end subroutine setGridlevel

   subroutine createChildForLvl(lvl, ixMin, ixMax)
      integer, intent(in) :: lvl, ixMin(3), ixMax(3)

      integer :: cnPoints(3), i, j, k
      double precision :: cPosMin(3), cPosMax(3), xyz(3), cDelta(3)

      grids(lvl)%hasChild = .true.
      grids(lvl)%child_ixMin = ixMin
      grids(lvl)%child_ixMax = ixMax

      ! Create a child grid at lvl+1
      cPosMin  = (ixMin-1) * grids(lvl)%delta + grids(lvl)%posMin
      cPosMax  = (ixMax-1) * grids(lvl)%delta + grids(lvl)%posMin
      cnPoints = 2 * (ixMax - ixMin) + 1
      cDelta   = grids(lvl)%delta * 0.5d0

      call setGridlevel(lvl+1, cnPoints, cPosMin, cPosMax)

      do k = 1, cnPoints(3)
         do j = 1, cnPoints(2)
            do i = 1, cnPoints(1)

               ! Interpolate the ion density from the parent grid
               xyz = (/i-1, j-1, k-1/) * cDelta + cPosMin
               grids(lvl+1)%posIonDens(i,j,k)   = get_value_at(grids(lvl), grids(lvl)%posIonDens, xyz)
               grids(lvl+1)%negIonDens(i,j,k)   = get_value_at(grids(lvl), grids(lvl)%negIonDens, xyz)
               grids(lvl+1)%O2MinDens(i,j,k)    = get_value_at(grids(lvl), grids(lvl)%O2MinDens, xyz)
               grids(lvl+1)%elecDens(i,j,k)     = get_value_at(grids(lvl), grids(lvl)%elecDens, xyz)
               grids(lvl+1)%excDens(i,j,k)      = get_value_at(grids(lvl), grids(lvl)%excDens, xyz)
            end do
         end do
      end do

      EFA_maxLvl = max(EFA_maxLvl, lvl+1)

   end subroutine createChildForLvl

   ! Ensure there is enough workspace for a grid of size Nx * Ny
   subroutine ensureWorkspace(lvl)
      integer, intent(in) :: lvl
      integer :: Nx, Ny, Nz, nWorkspace_required

      Nx = grids(lvl)%nPoints(1)
      Ny = grids(lvl)%nPoints(2)
      Nz = grids(lvl)%nPoints(3)

      nWorkspace_required = 30 + Nx + Ny + 5*Nz + MAX(Nx, Ny, Nz) + 7*(Nx/2 + Ny/2)

      if (nWorkspace_required > nWorkspace_allocated) then
         if (allocated(workspace)) deallocate(workspace)
         allocate( workspace(nWorkspace_required) )
         nWorkspace_allocated = nWorkspace_required
      end if

   end subroutine ensureWorkspace

   subroutine ensureGridIsAllocated(myGrid)
      type(amrGrid), intent(inout)  :: myGrid
      integer :: Nx, Ny, Nz

      Nx = myGrid%nPoints(1)
      Ny = myGrid%nPoints(2)
      Nz = myGrid%nPoints(3)

      if ( allocated(myGrid%Efield) ) then
         if ( any(shape(myGrid%Efield) /= (/3, Nx, Ny, Nz/)) ) then
            deallocate(myGrid%Efield)
         end if
      end if

      if ( .not. allocated(myGrid%Efield) ) allocate( myGrid%Efield(3, Nx, Ny, Nz) )
      myGrid%Efield = 0.0d0

      if ( allocated(myGrid%potential) ) then
         if ( any(shape(myGrid%potential) /= (/Nx, Ny, Nz/)) ) then
            deallocate(myGrid%potential)
         end if
      end if

      if ( .not. allocated(myGrid%potential) ) allocate( myGrid%potential(Nx, Ny, Nz) )
      myGrid%potential = 0.0d0

      if ( allocated(myGrid%elecDens) ) then
         if ( any(shape(myGrid%elecDens) /= (/Nx, Ny, Nz/)) ) then
            deallocate(myGrid%elecDens)
         end if
      end if

      if ( .not. allocated(myGrid%elecDens) ) allocate( myGrid%elecDens(Nx, Ny, Nz) )
      myGrid%elecDens = 0.0d0

      if ( allocated(myGrid%excDens) ) then
         if ( any(shape(myGrid%excDens) /= (/Nx, Ny, Nz/)) ) then
            deallocate(myGrid%excDens)
         end if
      end if

      if ( .not. allocated(myGrid%excDens) ) allocate( myGrid%excDens(Nx, Ny, Nz) )
      myGrid%excDens = 0.0d0

      if ( allocated(myGrid%posIonDens) ) then
         if ( any(shape(myGrid%posIonDens) /= (/Nx, Ny, Nz/)) ) then
            deallocate(myGrid%posIonDens)
         end if
      end if

      if ( .not. allocated(myGrid%posIonDens) ) allocate( myGrid%posIonDens(Nx, Ny, Nz) )
      myGrid%posIonDens = 0.0d0

      if ( allocated(myGrid%negIonDens) ) then
         if ( any(shape(myGrid%negIonDens) /= (/Nx, Ny, Nz/)) ) then
            deallocate(myGrid%negIonDens)
         end if
      end if

      if ( .not. allocated(myGrid%negIonDens) ) allocate( myGrid%negIonDens(Nx, Ny, Nz) )
      myGrid%negIonDens = 0.0d0

      if ( allocated(myGrid%O2MinDens) ) then
         if ( any(shape(myGrid%O2MinDens) /= (/Nx, Ny, Nz/)) ) then
            deallocate(myGrid%O2MinDens)
         end if
      end if

      if ( .not. allocated(myGrid%O2MinDens) ) allocate( myGrid%O2MinDens(Nx, Ny, Nz) )
      myGrid%O2MinDens = 0.0d0

      if ( allocated(myGrid%sourceTerm) ) then
         if ( any(shape(myGrid%sourceTerm) /= (/Nx, Ny, Nz/)) ) then
            deallocate(myGrid%sourceTerm)
         end if
      end if

      if ( .not. allocated(myGrid%sourceTerm) ) allocate( myGrid%sourceTerm(Nx, Ny, Nz) )
      myGrid%sourceTerm  = 0.0d0

      if ( allocated(myGrid%temp) ) then
         if ( any(shape(myGrid%temp) /= (/Nx, Ny, Nz/)) ) then
            deallocate(myGrid%temp)
         end if
      end if

      if ( .not. allocated(myGrid%temp) ) allocate( myGrid%temp(Nx, Ny, Nz) )
      myGrid%temp = 0.0d0

   end subroutine ensureGridIsAllocated

   subroutine removeGridData(myGrid)
      type(amrGrid), intent(inout) :: myGrid

      if ( allocated(myGrid%potential) ) then
         deallocate(myGrid%potential)
      end if

      if ( allocated(myGrid%elecDens) ) then
         deallocate(myGrid%elecDens)
      end if

      if ( allocated(myGrid%posIonDens) ) then
         deallocate(myGrid%posIonDens)
      end if

      if ( allocated(myGrid%negIonDens) ) then
         deallocate(myGrid%negIonDens)
      end if

      if ( allocated(myGrid%O2MinDens) ) then
         deallocate(myGrid%O2MinDens)
      end if

      if ( allocated(myGrid%temp) ) then
         deallocate(myGrid%temp)
      end if

      if ( allocated(myGrid%excDens) ) then
         deallocate(myGrid%excDens)
      end if

      if ( allocated(myGrid%sourceTerm) ) then
         deallocate(myGrid%sourceTerm)
      end if

      if ( allocated(myGrid%Efield) ) then
         deallocate(myGrid%Efield)
      end if

   end subroutine removeGridData

   !> Linear interpolation of f(x) between  x0 and x1, with f0 = f(x0) and f1 = f(x1)
   subroutine linearInterp(interp, x, x0, x1, f0, f1)
      DOUBLE PRECISION, INTENT(IN)     :: x, x0, x1
      DOUBLE PRECISION, INTENT(IN)     :: f0, f1
      DOUBLE PRECISION, INTENT(INOUT)  :: interp

      interp = f0 + (f1-f0) * ((x-x0)/(x1-x0))

   END subroutine linearInterp

   subroutine bilinearInterp(interp, x, y, x1, y1, x2, y2, f11, f12, f21, f22)
      double precision, intent(in)     :: x, y, x1, y1, x2, y2
      double precision, intent(in)     :: f11, f21, f12, f22
      double precision, intent(inout)  :: interp

      interp =    ( f11*(x2-x) + f21*(x-x1) ) * (y2-y) + &
           & ( f12*(x2-x) + f22*(x-x1) ) * (y-y1)
      interp = interp / ((x2-x1)*(y2-y1))

   end subroutine bilinearInterp

   !> Linear interpolation of vector function f(x) between x0 and x1, with f0 = f(x0) and f1 = f(x1)
   SUBROUTINE linearInterp3(interp, x, x0, x1, f0, f1)
      DOUBLE PRECISION, INTENT(IN)                 :: x, x0, x1
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN)   :: f0, f1
      DOUBLE PRECISION, DIMENSION(3), INTENT(INOUT):: interp

      interp = f0 + (f1-f0) * ((x-x0)/(x1-x0))

   END SUBROUTINE linearInterp3

   !> Bilinear interpolation of vector function f(x,y) on the domain  x1 <= x <= x2, y1 <= y <= y2
   SUBROUTINE bilinearInterp3(interp, x, y, x1, y1, x2, y2, f11, f12, f21, f22)
      DOUBLE PRECISION, INTENT(IN)                 :: x, y, x1, y1, x2, y2
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN)   :: f11, f21, f12, f22
      DOUBLE PRECISION, DIMENSION(3), INTENT(INOUT):: interp

      interp =    ( f11*(x2-x) + f21*(x-x1) ) * (y2-y) + &
           & ( f12*(x2-x) + f22*(x-x1) ) * (y-y1)
      interp = interp / ((x2-x1)*(y2-y1))

   END SUBROUTINE bilinearInterp3

   subroutine EFA_resetElecDens()
      integer :: lvl

      do lvl = -1, EFA_maxLvl
         grids(lvl)%elecDens = 0.0d0
      end do

   end subroutine EFA_resetElecDens

   subroutine EFA_resetIonDens()
      integer :: lvl

      do lvl = -1, EFA_maxLvl
         grids(lvl)%posIonDens = 0.0d0
         grids(lvl)%negIonDens = 0.0d0
         grids(lvl)%O2MinDens = 0.0d0
      end do

   end subroutine EFA_resetIonDens

   subroutine EFA_resetTemp()
      integer :: lvl

      do lvl = -1, EFA_maxLvl
         grids(lvl)%temp = 0.0d0
      end do

   end subroutine EFA_resetTemp

   subroutine EFA_addElectron(xyz, weight)
      double precision, intent(in)  :: xyz(3), weight
      integer                       :: lvl, maxlvl

      maxlvl = getMaxLevelFromZero(xyz)

      do lvl = -1, maxlvl
         call addParticleTo(grids(lvl), grids(lvl)%elecDens, xyz, weight)
      end do

   end subroutine EFA_addElectron

   subroutine EFA_addExcited(xyz, weight)
      double precision, intent(in)  :: xyz(3), weight
      integer                       :: lvl, maxlvl

      maxlvl = getMaxLevelFromZero(xyz)

      do lvl = -1, maxlvl
         call addParticleTo(grids(lvl), grids(lvl)%excDens, xyz, weight)
      end do

   end subroutine EFA_addExcited

   subroutine EFA_addParticleCount(xyz)
      double precision, intent(in)  :: xyz(3)
      integer                       :: lvl, maxlvl

      maxlvl = getMaxLevelFromZero(xyz)

      do lvl = -1, maxlvl
         call addParticleTo(grids(lvl), grids(lvl)%temp, xyz, 1.0d0)
      end do

   end subroutine EFA_addParticleCount

   subroutine EFA_addPosIon(xyz, weight)
      double precision, intent(in)  :: xyz(3), weight
      integer                       :: lvl, maxlvl

      maxlvl = getMaxLevelFromZero(xyz)

      do lvl = -1, maxlvl
         call addParticleTo(grids(lvl), grids(lvl)%posIonDens, xyz, weight)
      end do

   end subroutine EFA_addPosIon

   subroutine EFA_addNegIon(xyz, weight)
      double precision, intent(in)  :: xyz(3), weight
      integer                       :: lvl, maxlvl

      maxlvl = getMaxLevelFromZero(xyz)

      do lvl = -1, maxlvl
         call addParticleTo(grids(lvl), grids(lvl)%negIonDens, xyz, weight)
      end do

   end subroutine EFA_addNegIon

   subroutine EFA_addO2MinIon(xyz, weight)
      double precision, intent(in)  :: xyz(3), weight
      integer                       :: lvl, maxlvl

      maxlvl = getMaxLevelFromZero(xyz)

      do lvl = -1, maxlvl
         call addParticleTo(grids(lvl), grids(lvl)%O2MinDens, xyz, weight)
      end do

   end subroutine EFA_addO2MinIon

   subroutine addParticleTo(myGrid, array3D, xyz, weight)
      type(amrGrid), intent(inout)     :: myGrid
      double precision, intent(in)     :: xyz(3), weight
      double precision, intent(inout)  :: array3D(:,:,:)

      integer                       :: low(3), i
      double precision              :: rRel(3), delta(3), invDelta(3)
      double precision              :: rmin(3), rmax(3)
      double precision              :: rho(2,2,2), lowerCoeff(3), upperCoeff(3)
      double precision              :: temp(4), temp2(2)

      rRel     = xyz - myGrid%posMin
      delta    = myGrid%delta
      invDelta = myGrid%invDelta

      ! Find the indices low so that rRel lies in the cube given by [low, low+1]
      low   = floor(rRel * invDelta) + 1

      ! Bounds checking
      if (any(low < 1 .or. low > myGrid%nPoints-1)) return

      ! Set the min en max coordinates of the cube in which rRel lies
      rmax        = low * delta
      rmin        = rmax - delta
      lowerCoeff  = rmax - rRel
      upperCoeff  = rRel - rmin

      ! Now compute the coefficient of the charge on each of the 8 gridpoints at
      ! the corners of the cube, using linear interpolation (trilinear in this case)
      temp(1) = lowerCoeff(1) * lowerCoeff(2)
      temp(2) = lowerCoeff(1) * upperCoeff(2)
      temp(3) = upperCoeff(1) * lowerCoeff(2)
      temp(4) = upperCoeff(1) * upperCoeff(2)
      temp2   = (/lowerCoeff(3), upperCoeff(3)/)

      rho(1,1,1:2) = temp(1) * temp2
      rho(1,2,1:2) = temp(2) * temp2
      rho(2,1,1:2) = temp(3) * temp2
      rho(2,2,1:2) = temp(4) * temp2

      ! Scale rho to the right units and add it to the density
      rho = rho * weight / product(delta)**2

      array3D(low(1):low(1)+1, low(2):low(2)+1, low(3):low(3)+1) = &
           array3D(low(1):low(1)+1, low(2):low(2)+1, low(3):low(3)+1) + rho

   end subroutine addParticleTo

   subroutine setCellIxOffsets()
      integer :: lvl, cIxOffset

      cIxOffset = 0

      do lvl = 0, EFA_maxLvl
         grids(lvl)%cellIxOffset = cIxOffset
         cIxOffset               = cIxOffset + product(grids(lvl)%nPoints)
         !          print *, "At lvl", lvl, "offset:", grids(lvl)%cellIxOffset
      end do

   end subroutine setCellIxOffsets

   !> Return the unique index of the finest grid cell containing point xyz
   integer function EFA_getCellIndexAt(xyz)
      double precision, intent(in) :: xyz(3)

      integer :: lvl, nPoints(3), cellIx(3)

      lvl      = getMaxLevelFromZero(xyz)
      nPoints  = grids(lvl)%nPoints

      ! CellIx goes from 0 - nPoints-1
      cellIx   = floor( (xyz - grids(lvl)%posMin) * grids(lvl)%invDelta )

      ! Reverse order because of Fortran array order.
      EFA_getCellIndexAt = grids(lvl)%cellIxOffset + &
           cellIx(3) * nPoints(1) * nPoints(2) + &
           cellIx(2) * nPoints(1) + cellIx(1) + 1

   end function EFA_getCellIndexAt

   double precision function EFA_getNumElecAt(xyz)
      double precision, intent(in) :: xyz(3)
      integer                      :: lvl
      double precision             :: density

      lvl               = getMaxLevelFromZero(xyz)
      density           = EFA_getElecDensAt(xyz)
      EFA_getNumElecAt  = density * product(grids(lvl)%delta)
   end function EFA_getNumElecAt

   subroutine EFA_getDelta(lvl, delta)
      integer, intent(in)              :: lvl
      double precision, intent(inout)  :: delta(3)

      delta = grids(lvl)%delta
   end subroutine EFA_getDelta

   subroutine EFA_getGridSize(lvl, nIx)
      integer, intent(in)     :: lvl
      integer, intent(inout)  :: nIx(3)

      nIx = grids(lvl)%nPoints
   end subroutine EFA_getGridSize

   subroutine EFA_getPosMin(lvl, posMin)
      integer, intent(in)              :: lvl
      double precision, intent(inout)  :: posMin(3)

      posMin = grids(lvl)%posMin
   end subroutine EFA_getPosMin

   double precision function EFA_getMaxDeltaAt(xyz)
      double precision, intent(in)  :: xyz(3)
      integer                       :: lvl

      lvl               = getMaxLevelFromZero(xyz)
      EFA_getMaxDeltaAt = maxval(grids(lvl)%delta)
   end function EFA_getMaxDeltaAt

   integer function EFA_getMaxLvl()
      EFA_getMaxLvl = EFA_maxLvl
   end function EFA_getMaxLvl

   double precision function EFA_getElecDensAt(xyz)
      double precision, intent(in) :: xyz(3)
      integer                      :: lvl

      lvl               = getMaxLevelInsideAt(xyz)
      EFA_getElecDensAt = get_value_at(grids(lvl), grids(lvl)%elecDens, xyz)
   end function EFA_getElecDensAt

   double precision function EFA_getPosIonDensAt(xyz)
      double precision, intent(in) :: xyz(3)
      integer                      :: lvl

      lvl                  = getMaxLevelInsideAt(xyz)
      EFA_getPosIonDensAt  = get_value_at(grids(lvl), grids(lvl)%posIonDens, xyz)
   end function EFA_getPosIonDensAt

   integer function EFA_getMaxCellIx()
      EFA_getMaxCellIx = grids(EFA_maxLvl)%cellIxOffset + &
           product(grids(EFA_maxLvl)%nPoints)
   end function EFA_getMaxCellIx

   subroutine EFA_writeGridsToFile(filename, dirname, cycleNumber, time)
      character(len=*), intent(in) :: filename, dirname
      double precision, intent(in) :: time
      integer, intent(in)          :: cycleNumber

      integer :: lvl, i, j, k, Nx, Ny, Nz
      double precision :: dummy(2,2,2), xyz(3), delta(3), posmin(3)

      dummy = 0.0d0

      call computeEfieldAtLvl(lvl = -1)

      ! Interpolate electron and ion densities from parent grids
      do lvl = 1, EFA_maxLvl

         Nx       = grids(lvl)%nPoints(1)
         Ny       = grids(lvl)%nPoints(2)
         Nz       = grids(lvl)%nPoints(3)

         delta    = grids(lvl)%delta
         posmin   = grids(lvl)%posMin
         do k = 1, Nz, Nz-1
            do j = 1, Ny
               do i = 1, Nx
                  xyz = (/i-1, j-1, k-1/) * delta + posMin
                  grids(lvl)%elecDens(i,j,k) = get_value_at(grids(lvl-1), grids(lvl-1)%elecDens, xyz)
                  grids(lvl)%posIonDens(i,j,k) = get_value_at(grids(lvl-1), grids(lvl-1)%posIonDens, xyz)
                  grids(lvl)%negIonDens(i,j,k) = get_value_at(grids(lvl-1), grids(lvl-1)%negIonDens, xyz)
                  grids(lvl)%O2MinDens(i,j,k) = get_value_at(grids(lvl-1), grids(lvl-1)%O2MinDens, xyz)
                  grids(lvl)%temp(i,j,k) = get_value_at(grids(lvl-1), grids(lvl-1)%temp, xyz)
               end do
            end do
         end do

         do k = 1, Nz
            do j = 1, Ny, Ny-1
               do i = 1, Nx
                  xyz = (/i-1, j-1, k-1/) * delta + posMin
                  grids(lvl)%elecDens(i,j,k) = get_value_at(grids(lvl-1), grids(lvl-1)%elecDens, xyz)
                  grids(lvl)%posIonDens(i,j,k) = get_value_at(grids(lvl-1), grids(lvl-1)%posIonDens, xyz)
                  grids(lvl)%negIonDens(i,j,k) = get_value_at(grids(lvl-1), grids(lvl-1)%negIonDens, xyz)
                  grids(lvl)%O2MinDens(i,j,k) = get_value_at(grids(lvl-1), grids(lvl-1)%O2MinDens, xyz)
                  grids(lvl)%temp(i,j,k) = get_value_at(grids(lvl-1), grids(lvl-1)%temp, xyz)
               end do
            end do
         end do

         do k = 1, Nz
            do j = 1, Ny
               do i = 1, Nx, Nx-1
                  xyz = (/i-1, j-1, k-1/) * delta + posMin
                  grids(lvl)%elecDens(i,j,k) = get_value_at(grids(lvl-1), grids(lvl-1)%elecDens, xyz)
                  grids(lvl)%posIonDens(i,j,k) = get_value_at(grids(lvl-1), grids(lvl-1)%posIonDens, xyz)
                  grids(lvl)%negIonDens(i,j,k) = get_value_at(grids(lvl-1), grids(lvl-1)%negIonDens, xyz)
                  grids(lvl)%O2MinDens(i,j,k) = get_value_at(grids(lvl-1), grids(lvl-1)%O2MinDens, xyz)
                  grids(lvl)%temp(i,j,k) = get_value_at(grids(lvl-1), grids(lvl-1)%temp, xyz)
               end do
            end do
         end do

      end do

      do lvl = -1, EFA_maxLvl

         call addGridStructure3D(filename, dirname, "grid", grids(lvl)%potential, &
              grids(lvl)%posMin, grids(lvl)%posMax, cycleNumber, time, lvl+1)
         call addDataToGrid3D(filename, dirname, "potential", "grid", grids(lvl)%potential, "V", lvl+1)
         call addDataToGrid3D(filename, dirname, "elecDensx", "grid", grids(lvl)%elecDens, "1/m^3", lvl+1)
         call addDataToGrid3D(filename, dirname, "posIonDns", "grid", grids(lvl)%posIonDens, "1/m^3", lvl+1)
         call addDataToGrid3D(filename, dirname, "negIonDns", "grid", grids(lvl)%negIonDens, "1/m^3", lvl+1)
         call addDataToGrid3D(filename, dirname, "O2MinDens", "grid", grids(lvl)%O2MinDens, "1/m^3", lvl+1)
         call addDataToGrid3D(filename, dirname, "excDensxx", "grid", grids(lvl)%excDens, "1/m^3", lvl+1)
         call addDataToGrid3D(filename, dirname, "partDensx", "grid", grids(lvl)%temp, "[X]", lvl+1)
         call addDataToGrid3D(filename, dirname, "chargeDen", "grid", &
              grids(lvl)%posIonDens - grids(lvl)%negIonDens - grids(lvl)%elecDens, "1/m^3", lvl+1)
         call addDataToGrid3D(filename, dirname, "Elecfield", "grid", &
              sqrt(grids(lvl)%Efield(1,:,:,:)**2 &
              + grids(lvl)%Efield(2,:,:,:)**2 &
              + grids(lvl)%Efield(3,:,:,:)**2), "V/m", lvl+1)
      end do

      call setAmrStructure2D(filename, dirname, "amr", "grid", &
           (/"potential", "Elecfield", "elecDensx", "posIonDns", "negIonDns", &
           "O2MinDens", "excDensxx", "partDensx", "chargeDen"/), cycleNumber, time)

      call setAmrStructure2D(filename, dirname, "amr2", "grid", &
           (/"potential", "Elecfield", "elecDensx", "posIonDns", "negIonDns", &
           "O2MinDens", "excDensxx", "partDensx", "chargeDen"/), cycleNumber, time, forceMinLvl = 1)

   end subroutine EFA_writeGridsToFile

   subroutine EFA_printGridStructure(myrank)
      integer, intent(in) :: myrank

      integer :: lvl, rank, ntasks, ierr

      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)

      do rank = 0, ntasks-1

         if (rank == myrank) then
            do lvl = -1, EFA_maxLvl
               print *, "Rank:", myrank, "level:", lvl
               print *, "Size:", grids(lvl)%nPoints
               write(*, '(A,3E12.4)') "posMin:", grids(lvl)%posMin
               write(*, '(A,3E12.4)') "posMax:", grids(lvl)%posMax
               if (grids(lvl)%hasChild) then
                  print *, 'Child min:', grids(lvl)%child_ixMin
                  print *, 'Child max:', grids(lvl)%child_ixMax
               end if
               print *, ""
            end do
         end if

         call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      end do

   end subroutine EFA_printGridStructure

   !> To get the maximum electric field and position in the whole domain
   subroutine EFA_EfieldMax(maxWholeDomain, location)
      double precision, intent(out) :: maxWholeDomain
      double precision, intent(out) :: location

      integer :: subscript(3)

      subscript = maxloc(abs(grids(EFA_maxlvl)%Efield(3,:,:,:)))
      maxWholeDomain = maxval(abs(grids(EFA_maxlvl)%Efield(3,:,:,:)))

      location = grids(EFA_maxlvl)%posmin(3) + (subscript(3) - 1) * grids(EFA_maxlvl)%delta(3)

   end subroutine EFA_EfieldMax

   !> Anbang: Zener model in liquid, similar as detachment in air
   ! Before this routine is called, the densities should be collected at the root.
   ! It should only be called by the root.
   subroutine EFA_setNeutralMoleLoss(lvl, dt, neutralMoleLoss)
      integer, intent(in)             :: lvl
      double precision, intent(in)    :: dt
      double precision, intent(out)   :: neutralMoleLoss(:,:,:)

      integer :: i, j, k, nIx(3), ix(3)
      double precision :: cellVolume

      nIx = grids(lvl)%nPoints

      if (any(nIx /= shape(neutralMoleLoss))) then
         print *, "EFA_setNegIonLoss: argument negIonLoss has wrong shape"
         stop
      end if

      cellVolume = product(grids(lvl)%delta)
      do k = 1, nIx(3)
         do j = 1, nIx(2)
            do i = 1, nIx(1)
               ix = (/i, j, k/)
               neutralMoleLoss(i,j,k) = getNeutralMoleLoss(grids(lvl)%Efield(:,i,j,k), dt)

               ! In child region we should use values from the child, so set to zero.
               if (grids(lvl)%hasChild .and. &
                    all(ix >= grids(lvl)%child_ixMin) .and. all(ix <= grids(lvl)%child_ixMax)) then
                  neutralMoleLoss(i,j,k) = 0.0d0
               end if
            end do
         end do
      end do

      neutralMoleLoss = neutralMoleLoss * cellVolume
   end subroutine EFA_setNeutralMoleLoss

   double precision function getNeutralMoleLoss(Efield, dt)
      double precision, intent(in) :: Efield(3), dt

      double precision  :: densityIonizable   ! number of density of ionizable species, the number has to be checked
      double precision  :: distanceMolecular   ! molecular separation distance
      double precision  :: massEffective      ! effective electron mass, 0.1* me
      double precision  :: ionizationPotential  ! energy in Joule
      double precision  :: EfieldNorm


      EfieldNorm = sqrt(sum(Efield**2))
      densityIonizable = 1.d23
      distanceMolecular = 3.d-10
      massEffective = 0.1d0 * electronMass
      ionizationPotential = 2.2545d-18

      getNeutralMoleLoss = abs(electronCharge) * densityIonizable * distanceMolecular * EfieldNorm / PlanckConstant * &
            exp(-pi**2 * massEffective * distanceMolecular * (ionizationPotential)**2 / &
            (-electronCharge * PlanckConstant**2 * EfieldNorm))
      getNeutralMoleLoss = getNeutralMoleLoss * dt
   end function getNeutralMoleLoss

end module module_EfieldAMR
