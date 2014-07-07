! TODO: check the license on this code... Came somewhere from the web.

!> Module for the random number generation

MODULE module_kiss
   ! Jannis 17 feb 2011: made the kiss_rand() function cleaner and faster
   ! and introduced the kisset() function.
   IMPLICIT NONE
   PRIVATE
   PUBLIC  :: kiss 
   PUBLIC  :: kiss_rand
   PUBLIC  :: kiss_normal
   PUBLIC  :: kiss_Poisson
   PUBLIC  :: initializeKISS

   INTEGER :: xx=123456789, yy=362436069, zz=521288629, ww=916191069

   DOUBLE PRECISION, PARAMETER :: maxIntPlusOne = DBLE(HUGE(1)) + 1.0D0
   DOUBLE PRECISION, PARAMETER :: toUniformFactor = 0.5D0 / maxIntPlusOne

CONTAINS

   !> Initialize the seeds for the KISS random number generator
   SUBROUTINE initializeKISS(seeds, myrank)
      INTEGER, INTENT(IN) :: seeds(4), myrank
      INTEGER :: temp(4), n
      
      CALL seedKiss(seeds)
      
      ! Make sure different tasks are initialized differently
      DO n = 1, myrank
         temp = (/kiss(), kiss(), kiss(), kiss()/)
         CALL seedKiss(temp)
      END DO
      
   END SUBROUTINE initializeKISS
   
   !> Sets the seeds for the KISS RNG
   SUBROUTINE seedKiss(seeds)
      INTEGER, INTENT(IN) :: seeds(4)
      xx = seeds(1)
      yy = seeds(2)
      zz = seeds(3)
      ww = seeds(4)
   END SUBROUTINE seedKiss

   !> The  KISS (Keep It Simple Stupid) random number generator. Combines:
   !! (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
   !! (2) A 3-shift shift-register generator, period 2^32-1,
   !! (3) Two 16-bit multiply-with-carry generators, period 597273182964842497>2^59
   !!  Overall period>2^123;  Default seeds x,y,z,w.
   !!  Set your own seeds with statement i=kisset(ix,iy,iz,iw).
   INTEGER FUNCTION kiss()
      xx = 69069 * xx + 1327217885
      yy = mm (mm (mm (yy, 13), - 17), 5)
      zz = 18000 * iand (zz, 65535) + ishft (zz, - 16)
      ww = 30903 * iand (ww, 65535) + ishft (ww, - 16)
      kiss = xx + yy + ishft (zz, 16) + ww
   CONTAINS
      !> Small helper function for kiss() RNG.
      INTEGER FUNCTION mm(k, n)
         INTEGER :: k, n
         mm = ieor (k, ishft (k, n) )
      END FUNCTION mm
   END FUNCTION kiss

   !> Generate a uniform deviate between 0.0 (inclusive) and 1.0 (exclusive)
   !! RINT / (INT_MAX + 1) lies in [-1, 1), so 0.5 * RINT / (INT_MAX + 1) + 0.5
   !! should be in the interval [0,1)
   DOUBLE PRECISION FUNCTION kiss_rand()
      kiss_rand = DBLE(kiss()) * toUniformFactor + 0.5D0
   END FUNCTION kiss_rand
   
   !> Return normal random variate with mean 0 and variance 1
   DOUBLE PRECISION FUNCTION kiss_normal()
      DOUBLE PRECISION, SAVE :: prevNumber
      LOGICAL, SAVE :: haveNumber = .FALSE.
      
      DOUBLE PRECISION :: z0, z1, u0, u1
      
      IF (haveNumber) THEN
         haveNumber = .FALSE.
         kiss_normal = prevNumber
         RETURN
      ELSE
         u0 = kiss_rand()
         u1 = kiss_rand()
         z0 = sqrt(-2.0D0 * log(1.0D0 - u0)) * cos(2.0D0 * 3.14159265358979324D0 * u1)
         z1 = sqrt(-2.0D0 * log(1.0D0 - u0)) * sin(2.0D0 * 3.14159265358979324D0 * u1)
         
         haveNumber = .TRUE.
         prevNumber = z1
         kiss_normal = z0
         RETURN
      END IF
      
   END FUNCTION kiss_normal
   
   !> Return Poisson random variate with rate labda. Works well for labda < 30 or so.
   !! For labda >> 1 it can produce wrong results due to roundoff error.
   INTEGER FUNCTION kiss_Poisson(labda)
      DOUBLE PRECISION, INTENT(IN) :: labda
      INTEGER :: k
      DOUBLE PRECISION :: expL, p
      
      expL = exp(-labda)
      k = 0
      p = kiss_rand()
      
      DO WHILE (p > expL)
         k = k + 1
         p = p * kiss_rand()
      END DO
      kiss_Poisson = k
      
   END FUNCTION kiss_Poisson

END MODULE module_kiss
