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

!> Module that contains variables that are shared between modules ('globals')
!!
!! This module contains the variables that are shared between modules, together
!! with some function that operate on them.

MODULE module_globals

   USE module_config
   USE module_constants
   use generalUtilities

   IMPLICIT NONE
   PUBLIC

   !! Global variables

   !> Whether photoionization is included
   LOGICAL     :: GL_flagPhotoIonization

   !> Vector of the sizes of the grid (in the x,y,z direction)
   INTEGER, DIMENSION(3) :: GL_gridSize

   !> Volume of a grid cell
   DOUBLE PRECISION :: GL_gridVolume

   !> Inverse volume of a grid cell
   DOUBLE PRECISION :: GL_invGridVolume

   !> The total number of grid cells for the potential, electric field, density etc.
   INTEGER          :: GL_nGridCells

   !> Vector of the spatial stepsize on the grid (in x,y,z direction)
   DOUBLE PRECISION, DIMENSION(3) :: GL_gridDelta
   DOUBLE PRECISION, DIMENSION(3) :: GL_invGridDelta

   !> Vector of the total length of the grid (in x,y,z direction)
   DOUBLE PRECISION, DIMENSION(3) :: GL_gridLength

   !> Controls the weighting scheme used, supported are zeroth-order (or nearest grid point, NGP)
   !! and first order (or cloud in cell, CIC).
   INTEGER :: weightingScheme

   PRIVATE :: CICaddToDens, NGPaddToDens, weightingScheme

CONTAINS

   !> Initialization routine for the global variables
   !!
   !! Allocates space for the arrays and initializes values based on the configuration
   !! that has been read in.
   SUBROUTINE GL_initializeVariables(myrank, rootMPI)
      use mpi
      ! Allocates all the non-fixed sized arrays and do initialization
      INTEGER, INTENT(IN) :: myrank, rootMPI

      ! Set the global parameters
      GL_flagPhotoIonization  = CFG_varLogic("PI_enabled")

      CALL CFG_getVar("grid_size", GL_gridSize)
      CALL CFG_getVar("grid_delta", GL_gridDelta)
      GL_invGridDelta = 1.0D0 / GL_gridDelta
      GL_nGridCells = PRODUCT(GL_gridSize + 1)

      GL_gridLength     = (GL_gridSize - 1) * GL_gridDelta
      GL_gridVolume     = PRODUCT(GL_gridDelta)
      GL_invGridVolume  = 1.0D0 / GL_gridVolume

      weightingScheme = CFG_varInt("part_weightingScheme")
   END SUBROUTINE GL_initializeVariables

   !> Store the nearest grid cell at position 'pos' in 'loc'
   SUBROUTINE GL_findIndex(pos, loc)
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN)   :: pos
      INTEGER, DIMENSION(3), INTENT(OUT)           :: loc

      loc = NINT(pos*GL_invGridDelta) + 1
   END SUBROUTINE GL_findIndex

   !> Return true if 'pos' is inside the simulation domain (including any electrode)
   LOGICAL FUNCTION GL_inDomain(pos)
      DOUBLE PRECISION, INTENT(IN), DIMENSION(3) :: pos
      GL_inDomain = ALL( pos >= 0.0D0 ) .AND. ALL( pos <= GL_gridLength)
   END FUNCTION GL_inDomain

   SUBROUTINE NGPaddToDens(dens, pos, weightFactor)
      ! NGP (nearest grid point) weighting scheme to add particles to a density
      ! This scheme is zeroth order in the weighting
      DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(INOUT)  :: dens
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN)         :: pos
      DOUBLE PRECISION, INTENT(IN)                       :: weightFactor

      INTEGER, DIMENSION(3) :: ix

      CALL GL_findIndex(pos, ix)

      dens(ix(1), ix(2), ix(3)) = &
         & dens(ix(1), ix(2), ix(3)) + weightFactor * GL_invGridVolume

   END SUBROUTINE NGPaddToDens

   SUBROUTINE CICaddToDens(dens, pos, weightFactor)
      ! CIC (cloud in cell) weighting scheme to add particles to a density
      ! This scheme is first order in the weighting
      DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(INOUT)  :: dens
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN)         :: pos
      DOUBLE PRECISION, INTENT(IN)                       :: weightFactor

      ! Temporary variables
      INTEGER, DIMENSION(3)               :: low, upp
      DOUBLE PRECISION, DIMENSION(3)      :: rmin, rmax
      DOUBLE PRECISION, DIMENSION(2,2,2)  :: rho

      ! Find the indices low, upp so that pos lies in the cube given by [low, upp]
      low = FLOOR(pos/GL_gridDelta) + 1
      upp = low + 1
      ! Set the min en max coordinates of the cube in which pos lies
      rmin = (low-1) * GL_gridDelta
      rmax = (upp-1) * GL_gridDelta

      ! Now compute the coefficient of the charge on each of the 8 gridpoints at
      ! the corners of the cube, using linear interpolation (trilinear in this case)
      rho(1,1,1) = (rmax(1) - pos(1)) * (rmax(2) - pos(2)) * (rmax(3) - pos(3))
      rho(1,1,2) = (rmax(1) - pos(1)) * (rmax(2) - pos(2)) * (pos(3) - rmin(3))
      rho(1,2,1) = (rmax(1) - pos(1)) * (pos(2) - rmin(2)) * (rmax(3) - pos(3))
      rho(1,2,2) = (rmax(1) - pos(1)) * (pos(2) - rmin(2)) * (pos(3) - rmin(3))

      rho(2,1,1) = (pos(1) - rmin(1)) * (rmax(2) - pos(2)) * (rmax(3) - pos(3))
      rho(2,1,2) = (pos(1) - rmin(1)) * (rmax(2) - pos(2)) * (pos(3) - rmin(3))
      rho(2,2,1) = (pos(1) - rmin(1)) * (pos(2) - rmin(2)) * (rmax(3) - pos(3))
      rho(2,2,2) = (pos(1) - rmin(1)) * (pos(2) - rmin(2)) * (pos(3) - rmin(3))

      ! Scale rho to the right units and add it to the density
      rho = rho * weightFactor * (GL_invGridVolume**2)

      dens(low(1), low(2), low(3)) = dens(low(1), low(2), low(3)) + rho(1,1,1)
      dens(low(1), low(2), upp(3)) = dens(low(1), low(2), upp(3)) + rho(1,1,2)
      dens(low(1), upp(2), low(3)) = dens(low(1), upp(2), low(3)) + rho(1,2,1)
      dens(low(1), upp(2), upp(3)) = dens(low(1), upp(2), upp(3)) + rho(1,2,2)

      dens(upp(1), low(2), low(3)) = dens(upp(1), low(2), low(3)) + rho(2,1,1)
      dens(upp(1), low(2), upp(3)) = dens(upp(1), low(2), upp(3)) + rho(2,1,2)
      dens(upp(1), upp(2), low(3)) = dens(upp(1), upp(2), low(3)) + rho(2,2,1)
      dens(upp(1), upp(2), upp(3)) = dens(upp(1), upp(2), upp(3)) + rho(2,2,2)

   END SUBROUTINE CICaddToDens

   SUBROUTINE CICgetCoeff(coeff, pos, weightFactor)
      DOUBLE PRECISION, DIMENSION(8), INTENT(INOUT)      :: coeff
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN)         :: pos
      DOUBLE PRECISION, INTENT(IN)                       :: weightFactor

      ! Temporary variables
      INTEGER, DIMENSION(3)               :: low, upp
      DOUBLE PRECISION, DIMENSION(3)      :: rmin, rmax

      ! Find the indices low, upp so that pos lies in the cube given by [low, upp]
      low = FLOOR(pos/GL_gridDelta) + 1
      upp = low + 1
      ! Set the min en max coordinates of the cube in which pos lies
      rmin = (low-1) * GL_gridDelta
      rmax = (upp-1) * GL_gridDelta

      ! Now compute the coefficient of the charge on each of the 8 gridpoints at
      ! the corners of the cube, using linear interpolation (trilinear in this case)
      coeff(1) = (rmax(1) - pos(1)) * (rmax(2) - pos(2)) * (rmax(3) - pos(3))
      coeff(2) = (rmax(1) - pos(1)) * (rmax(2) - pos(2)) * (pos(3) - rmin(3))
      coeff(3) = (rmax(1) - pos(1)) * (pos(2) - rmin(2)) * (rmax(3) - pos(3))
      coeff(4) = (rmax(1) - pos(1)) * (pos(2) - rmin(2)) * (pos(3) - rmin(3))

      coeff(5) = (pos(1) - rmin(1)) * (rmax(2) - pos(2)) * (rmax(3) - pos(3))
      coeff(6) = (pos(1) - rmin(1)) * (rmax(2) - pos(2)) * (pos(3) - rmin(3))
      coeff(7) = (pos(1) - rmin(1)) * (pos(2) - rmin(2)) * (rmax(3) - pos(3))
      coeff(8) = (pos(1) - rmin(1)) * (pos(2) - rmin(2)) * (pos(3) - rmin(3))

      ! Scale rho to the right units and add it to the density
      coeff = coeff * weightFactor * (GL_invGridVolume**2)

   END SUBROUTINE CICgetCoeff

   !> Add a particle with weigth 'weightFactor' to the density 'dens' at position 'pos'.
   !!
   !! The weighting scheme used depends on the parameter 'weightingScheme', currently
   !! zeroth-order and first order are supported.
   SUBROUTINE GL_addToDensity(dens, pos, weightFactor)
      DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(INOUT)  :: dens
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN)         :: pos
      DOUBLE PRECISION, INTENT(IN)                       :: weightFactor

      ! Select the weighting method
      SELECT CASE (weightingScheme)
         CASE (0)
            CALL NGPaddToDens(dens, pos, weightFactor)
         CASE (1)
            CALL CICaddToDens(dens, pos, weightFactor)
      END SELECT

   END SUBROUTINE GL_addToDensity

   !> Set the density 'dens' to zero
   SUBROUTINE GL_resetDensity(dens)
      DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(INOUT) :: dens

      dens(:,:,:) = 0.0D0

   END SUBROUTINE GL_resetDensity


END MODULE module_globals
