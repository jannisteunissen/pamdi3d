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


MODULE generalUtilities

   use module_constants

   IMPLICIT NONE
   PRIVATE

   interface binarySearch
      module procedure binarySearchDble, binarySearchInt
   end interface binarySearch

   PUBLIC :: twoNorm
   PUBLIC :: v2en
   PUBLIC :: velocityToEnergy
   PUBLIC :: en2v
   PUBLIC :: binarySearch
   PUBLIC :: linearInterpolateList
   PUBLIC :: intInRange
   PUBLIC :: swapInt
   PUBLIC :: swapDble
   PUBLIC :: smoothArray3D
   PUBLIC :: crossProduct

CONTAINS

   SUBROUTINE swapInt(a, b)
      INTEGER, INTENT(INOUT) :: a, b
      INTEGER                :: temp
      temp = a
      a = b
      b = temp
   END SUBROUTINE swapInt

   SUBROUTINE swapDble(a, b)
      DOUBLE PRECISION, INTENT(INOUT) :: a, b
      DOUBLE PRECISION                :: temp
      temp = a
      a = b
      b = temp
   END SUBROUTINE swapDble

   !> Smooth the array in 3D
   SUBROUTINE smoothArray3D(array3d)
      DOUBLE PRECISION, INTENT(INOUT) :: array3d(:,:,:)

      DOUBLE PRECISION, ALLOCATABLE :: tempX(:), tempY(:), tempZ(:)

      integer :: i, j, k, nx, ny, nz
      double precision :: filter(3) = (/0.25D0, 0.5D0, 0.25D0/)

      nx = size(array3d, 1);
      ny = size(array3d, 2);
      nz = size(array3d, 3);

      allocate( tempX(nx) )
      allocate( tempY(ny) )
      allocate( tempZ(nz) )

      ! Apply filter in x-direction
      DO k = 1, nz
         DO j = 1, ny
            CALL smoothArray1D(array3d(:, j, k), tempX, filter)
            array3d(:, j, k) = tempX
         END DO
      END DO

      ! Apply filter in y-direction
      DO k = 1, nz
         DO i = 1, nx
            CALL smoothArray1D(array3d(i, :, k), tempY, filter)
            array3d(i, :, k) = tempY
         END DO
      END DO

      ! Apply filter in z-direction
      DO j = 1, ny
         DO i = 1, nx
            CALL smoothArray1D(array3d(i, j, :), tempZ, filter)
            array3d(i, j, :) = tempZ
         END DO
      END DO

      deallocate( tempX )
      deallocate( tempY )
      deallocate( tempZ )

   END SUBROUTINE smoothArray3D

   !> Smooth the array in 1D
   SUBROUTINE smoothArray1D(array_in, array_out, filter)
      DOUBLE PRECISION, INTENT(IN) :: array_in(:), filter(:)
      DOUBLE PRECISION, INTENT(INOUT) :: array_out(:)

      integer :: i, k, nx, nf, mid

      nx = size(array_in)
      nf = size(filter)
      mid = (nf-1)/2

      array_out(:) = 0.0D0

      ! Apply filter
      DO i = 1+mid, nx-mid
         DO k = -mid, mid
            array_out(i) = array_out(i) + array_in(i+k) * filter(k+mid+1)
         END DO
      END DO
   END SUBROUTINE smoothArray1D

   INTEGER FUNCTION binarySearchDble(list, value)
      ! Searches list for the interval containing value, such that
      ! list(i) <= value < list(i+1), and returns i
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: list
      DOUBLE PRECISION, INTENT(IN) :: value

      INTEGER :: iMin, iMax, iMiddle

      iMin = 1
      iMax = size(list)

      ! Check case where list has length 0
      if (iMax == 0) then
         binarySearchDble = 0
         return
      end if

      DO WHILE (iMax - iMin > 1)
         iMiddle = iMin + (iMax - iMin) / 2

         IF ( value < list(iMiddle) ) THEN
            iMax = iMiddle
         ELSE
            iMin = iMiddle
         END IF
      END DO

      binarySearchDble = iMin
   END FUNCTION binarySearchDble

   INTEGER FUNCTION binarySearchInt(list, value)
      ! Searches list for the interval containing value, such that
      ! list(i) <= value < list(i+1), and returns i
      integer, DIMENSION(:), INTENT(IN) :: list
      integer, INTENT(IN) :: value

      INTEGER :: iMin, iMax, iMiddle

      iMin = 1
      iMax = size(list)

      ! Check case where list has length 0
      if (iMax == 0) then
         binarySearchInt = 0
         return
      end if

      DO WHILE (iMax - iMin > 1)
         iMiddle = iMin + (iMax - iMin) / 2

         IF ( value < list(iMiddle) ) THEN
            iMax = iMiddle
         ELSE
            iMin = iMiddle
         END IF
      END DO

      binarySearchInt = iMin
   END FUNCTION binarySearchInt

   LOGICAL FUNCTION intInRange(i, iMin, iMax)
      INTEGER, INTENT(IN), DIMENSION(:) :: i, iMin, iMax
      intInRange = ALL( i >= iMin ) .AND. ALL( i <= iMax)
   END FUNCTION intInRange

   SUBROUTINE linearInterpolateList(xList, yList, xValue, yValue)
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN)   :: xList, yList
      DOUBLE PRECISION, INTENT(IN)                 :: xValue
      DOUBLE PRECISION, INTENT(INOUT)              :: yValue

      INTEGER                                      :: i, iMin, iMax
      DOUBLE PRECISION                             :: temp

      iMin = LBOUND(xList, 1)
      iMax = UBOUND(xList, 1)

      IF (xValue <= xList(iMin)) THEN
         yValue   = yList(iMin)
      ELSE IF (xValue >= xList(iMax)) THEN
         yValue   = yList(iMax)
      ELSE
         i        = binarySearch(xList, xValue)
         temp     = (xValue - xList(i)) / (xList(i+1) - xList(i))
         yValue   = (1.0D0 - temp) * yList(i) + temp * yList(i+1)
      END IF

   END SUBROUTINE linearInterpolateList

   DOUBLE PRECISION FUNCTION twoNorm(vec)
      ! Computes the two-norm of a vector
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: vec

      twoNorm = SQRT( SUM(vec**2) )
   END FUNCTION twoNorm

   DOUBLE PRECISION FUNCTION v2en(v)
      ! Return the kinetic energy corresponding to an electron velocity v
      DOUBLE PRECISION,intent(IN) :: v
      v2en=0.5D0*electronMass*v**2
   END FUNCTION v2en

   DOUBLE PRECISION FUNCTION velocityToEnergy(velocity)
      ! Return the kinetic energy corresponding to an electron velocity
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: velocity
      velocityToEnergy = 0.5D0 * electronMass * SUM(velocity**2)
   END FUNCTION velocityToEnergy

   DOUBLE PRECISION FUNCTION en2v(en)
      ! Returns the electron velocity corresponding to a kinetic energy of en
      DOUBLE PRECISION, INTENT(IN) :: en
      en2v=sqrt(2.0D0*en/electronMass)
   END FUNCTION en2v

   subroutine crossProduct(vect1,vect2,vect)
      double precision, dimension(3),intent(in) :: vect1, vect2
      double precision, dimension(3),intent(out) :: vect

      vect(1) = vect1(2) * vect2(3) - vect1(3) * vect2(2)
      vect(2) = vect1(3) * vect2(1) - vect1(1) * vect2(3)
      vect(3) = vect1(1) * vect2(2) - vect1(2) * vect2(1)


   end subroutine crossProduct

END MODULE generalUtilities
