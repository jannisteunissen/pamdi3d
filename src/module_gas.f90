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


MODULE module_gas
   ! Module for gas properties (type, fraction, pressure, temperature etc)
   USE module_config
   USE module_constants
   
   IMPLICIT NONE
   
   PRIVATE
   
   INTEGER :: GAS_nGases = 0
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: GAS_fractions
   CHARACTER(LEN=40), DIMENSION(:), ALLOCATABLE :: GAS_names
   CHARACTER(LEN=40), DIMENSION(:), ALLOCATABLE :: GAS_files
   
   DOUBLE PRECISION :: GAS_pressure, GAS_temperature, GAS_numberDensity
   
   ! Methods
   PUBLIC :: getGasFraction, GAS_initialize
   ! Variables
   PUBLIC :: GAS_nGases, GAS_pressure, GAS_temperature
   PUBLIC :: GAS_fractions, GAS_names, GAS_files, GAS_numberDensity
   PUBLIC :: attachmentRateCounter
   
CONTAINS
   SUBROUTINE GAS_initialize()
      
      ! Find the number of gases used
      GAS_nGases = CFG_getSize("sim_gasFractions")
      
      ! Check whether the config entries have the same length
      if (CFG_getSize("sim_gasNames") /= GAS_nGases .or. &
          CFG_getSize("sim_gasFiles") /= GAS_nGases) then
         print *, "GAS_initialize error: make sure that sim_gasNames, ", &
                  "sim_gasFiles and sim_gasFractions have the same length"
         stop
      end if
      
      ALLOCATE( GAS_fractions(GAS_nGases) )
      ALLOCATE( GAS_names(GAS_nGases) )
      ALLOCATE( GAS_files(GAS_nGases) )
      
      CALL CFG_getVar("sim_gasFractions", GAS_fractions)
      CALL CFG_getVar("sim_gasNames", GAS_names)
      CALL CFG_getVar("sim_gasFiles", GAS_files)
      
      GAS_pressure = CFG_varDble("sim_gasPressure")
      GAS_temperature = CFG_varDble("sim_gasTemperature")
      
      ! Ideal gas law, pressure is in bar
      GAS_numberDensity = 1.0D5 * GAS_pressure / (BoltzmannConstant * GAS_temperature)
      print *, 'GAS_numberDensity = ', GAS_numberDensity
      
   END SUBROUTINE GAS_initialize
   
   DOUBLE PRECISION FUNCTION getGasFraction(gasName)
      CHARACTER(LEN=*), INTENT(IN) :: gasName
      INTEGER :: n
      DO n = 1, GAS_nGases
         IF ( gasName == TRIM(GAS_names(n)) ) THEN
            getGasFraction = GAS_fractions(n)
            RETURN
         END IF
      END DO
      
      print *, "getGasFraction warning: ", gasName, " is not found"
      getGasFraction = 0.0D0
      stop

   END FUNCTION getGasFraction
   
  !> the background electron generate rate: due to the loss of electrons of attachment
  !> paper "S. Pancheshnyi. Plasma Sources Sci. T. 14:645, 2005"
  subroutine attachmentRateCounter(attachmentRate)
   
      double precision :: fractionO2,fractionN2,backgroundDensity
      integer :: ll
      double precision,dimension(2) :: attachCoeff
      double precision, intent(out) :: attachmentRate
      
      backgroundDensity = CFG_varDble("init_backgroundDensity")
      
      attachCoeff(1) =  2.d-42      !unit: m^6 s^-1; e+O2+O2--O2- + O2
      attachCoeff(2) =  8.d-44      !unit: m^6 s^-1; e+O2+N2--O2- + N2
      
      do ll = 1, GAS_nGases
         if (GAS_names(ll) == "O2") then
            fractionO2 = GAS_fractions(ll)
         else if (GAS_names(ll) == "N2") then
            fractionN2 = GAS_fractions(ll)
         end if
      end do
           
      attachmentRate = attachCoeff(1) * backgroundDensity * (fractionO2 * GAS_numberDensity) **2 &
               & + attachCoeff(2) * backgroundDensity * (fractionN2 *fractionO2) * GAS_numberDensity **2
            
   end subroutine attachmentRateCounter
   
END MODULE module_gas