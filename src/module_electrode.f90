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

!> Module that enables the use of an electrode
module module_electrode

  use module_config
  use module_globals
  use module_constants
  use generalUtilities
  use module_kiss

  implicit none
  private

  integer :: EL_nPoints
  double precision :: EL_voltage, EL_xyzPos(3)
  double precision :: EL_topAngle, EL_topZ, EL_TransCurveStartZ, EL_TransCurveEndZ, EL_spacing
  double precision :: EL_tipCurveBeginZ, EL_Rcyl, EL_RcTrans, EL_Hcone, EL_Hcyl, EL_RcTip

  double precision, allocatable :: EL_surfacePoints(:,:)
  double precision, allocatable :: EL_charges(:)
  double precision, allocatable :: EL_weightFactors(:)
  double precision, allocatable :: elec_time_voltage(:,:)

  public :: EL_getVoltage, EL_getCharge, EL_getNPoints, EL_getWeightFactor, EL_getSurfacePoint
  public :: EL_setCharge, EL_setWeightFactor, EL_initialize, EL_insideElectrode, EL_getRadius
  public :: EL_getBottomPos
  public :: EL_getTopPos
  public :: EL_around_electrode_tip

contains

  !> Generates a list of points that lie on the electrode surface, where we want to have the electrode potential.
  !!
  !! Currently this list of point remains constant during the simulation.
  subroutine EL_initialize(myrank, root)
    integer, intent(in) :: myrank, root
    double precision :: H_cone_eff
    integer :: n_time_points

    ! Set up a basic electrode, with conical part and cylindrical part
    ! This should be moved to a subroutine, so that we can select different electrodes
    EL_voltage     = CFG_varDble("electrode_voltage")
    EL_Hcone     = CFG_varDble("electrode_Hcone")
    !       EL_Hcyl      = CFG_varDble("electrode_Hcyl")

    EL_Rcyl      = CFG_varDble("electrode_Rcyl")
    EL_RcTip     = CFG_varDble("electrode_RcTip")
    EL_RcTrans   = CFG_varDble("electrode_RcTrans")
    EL_spacing   = CFG_varDble("electrode_spacing")

    call CFG_getVar("electrode_xyzRelPos", EL_xyzPos)
    EL_topAngle  = atan(EL_Rcyl / EL_Hcone)
    H_cone_eff = EL_Hcone - EL_Rctip * (1.0D0/sin(EL_topAngle) - 1.0D0)
    EL_Hcyl      = EL_xyzPos(3) * GL_gridLength(3) - H_cone_eff - EL_spacing
    EL_topZ      = EL_Hcyl + H_cone_eff
    EL_xyzPos    = EL_xyzPos * GL_gridLength - (/0.0d0, 0.0d0, EL_topZ/)

    ! Determine the transitions between rounded and straight parts
    EL_TransCurveStartZ = EL_Hcyl - EL_RcTrans * tan(EL_topAngle/2.0D0)
    EL_TransCurveEndZ   = EL_TransCurveStartZ + EL_RcTrans * sin(EL_topAngle)
    EL_tipCurveBeginZ   = EL_Hcyl + EL_Hcone - EL_Rctip * (1.0D0/sin(EL_topAngle) - sin(EL_topAngle))
    n_time_points = CFG_getSize("elec_times")
    allocate(elec_time_voltage(n_time_points, 2))
    call CFG_getVar("elec_times", elec_time_voltage(:, 1))
    call CFG_getVar("elec_voltages", elec_time_voltage(:,2))

    if (myrank == root) call generate_electrode_surface()

  end subroutine EL_initialize

  subroutine EL_getBottomPos(coordinates)
    double precision, intent(out) :: coordinates(3)
    coordinates = EL_xyzPos
  end subroutine EL_getBottomPos

  subroutine EL_getTopPos(coordinates)
    double precision, intent(out) :: coordinates(3)
    coordinates = EL_xyzPos
    coordinates(3) = coordinates(3) + EL_topZ
 end subroutine EL_getTopPos

  double precision function EL_getRadius(zCoord)
    double precision, intent(in) :: zCoord

    if (zCoord < EL_TransCurveStartZ) then ! Cylindrical region
       EL_getRadius = EL_Rcyl
    else if (zCoord < EL_TransCurveEndZ) then ! Rounded part between conical and cylindrical region
       EL_getRadius = EL_RcTrans * sqrt(1.0D0 - (zCoord - EL_TransCurveStartZ)**2 / EL_RcTrans**2)
       EL_getRadius = EL_getRadius - EL_RcTrans + EL_Rcyl
    else if (zCoord < EL_tipCurveBeginZ) then ! Conical part
       EL_getRadius = (1.0D0 - (zCoord - EL_Hcyl)/EL_Hcone) * EL_Rcyl
    else ! Rounded tip
       EL_getRadius = sqrt(abs(EL_RcTip**2 - (zCoord - EL_topZ + EL_RcTip)**2))
    end if
  end function EL_getRadius

  subroutine generate_electrode_surface()
    integer :: ix, iy, nPoints, nVertical
    integer, allocatable :: nHorizontal(:)
    double precision :: height, radius, xPos, yPos
    double precision :: randAngle, usedSpacing

    double precision, parameter :: min_rel_radius = 0.01d0

    ! Determine the number of points
    height         = 0.0d0
    nVertical      = 0
    do while (height < EL_topZ)
       radius      = EL_getRadius(height)
       nVertical   = nVertical + 1
       if (height < EL_TransCurveStartZ) then
          usedSpacing = EL_spacing * (1.0d0 + 7.0d0 * (EL_TransCurveStartZ - height) / EL_TransCurveStartZ)
       else
          usedSpacing = EL_spacing * max(min_rel_radius, radius / EL_Rcyl)
       end if
       height = height + usedSpacing
    end do

    allocate( nHorizontal(nVertical) )

    height         = 0.0d0
    ix             = 0
    do while (height < EL_topZ)
       radius      = EL_getRadius(height)
       ix          = ix + 1
       if (height < EL_TransCurveStartZ) then
          usedSpacing = EL_spacing * (1.0d0 + 7.0d0 * (EL_TransCurveStartZ - height) / EL_TransCurveStartZ)
       else
          usedSpacing = EL_spacing * max(min_rel_radius, radius / EL_Rcyl)
       end if
       nHorizontal(ix) = nint(2.0d0 * pi * max(radius, min_rel_radius*EL_Rcyl) / usedSpacing) + 1
       height = height + usedSpacing
    end do

    ! Now store the surface points
    EL_nPoints = sum(nHorizontal)
    print *, "EL_nPoints: ", EL_nPoints
    allocate( EL_surfacePoints(3, EL_nPoints) )
    allocate( EL_charges(EL_nPoints) )
    allocate( EL_weightFactors(EL_nPoints) )

    EL_charges        = 0.0d0
    EL_weightFactors  = 0.0d0
    nPoints           = 0

    height         = 0.0d0
    ix             = 0
    do while (height < EL_topZ)
       radius      = EL_getRadius(height)
       ix          = ix + 1
       randAngle   = kiss_rand() * 2.0d0 * pi
       do iy = 1, nHorizontal(ix)
          xPos = radius * cos(randAngle + dble(iy-1) * 2.0D0 * pi / nHorizontal(ix))
          yPos = radius * sin(randAngle + dble(iy-1) * 2.0D0 * pi / nHorizontal(ix))
          nPoints = nPoints + 1
          EL_surfacePoints(:, nPoints) = EL_xyzPos + (/xPos, yPos, height/)
       end do

       if (height < EL_TransCurveStartZ) then
          usedSpacing = EL_spacing * (1.0d0 + 7.0d0 * (EL_TransCurveStartZ - height) / EL_TransCurveStartZ)
       else
          usedSpacing = EL_spacing * max(min_rel_radius, radius / EL_Rcyl)
       end if
       height = height + usedSpacing
    end do

    deallocate( nHorizontal )
  end subroutine generate_electrode_surface

  !> Return .TRUE. if the location (x,y,z) is inside the simulated electrode
  logical function EL_insideElectrode(pos)
    double precision, intent(IN) :: pos(3)
    double precision :: dist, radius, relPos(3)

    EL_insideElectrode = .false.
    relPos = pos - EL_xyzPos

    if (relPos(3) < EL_topZ) then
       dist = twoNorm(relPos(1:2))
       if (dist < EL_Rcyl) then
          if (relPos(3) < EL_TransCurveStartZ) then ! Cylindrical region
             radius = EL_Rcyl
          else if (relPos(3) < EL_TransCurveEndZ) then ! Rounded part between conical and cylindrical region
             radius = EL_RcTrans * sqrt(1.0D0 - (relPos(3) - EL_TransCurveStartZ)**2 / EL_RcTrans**2)
             radius = radius - EL_RcTrans + EL_Rcyl
          else if (relPos(3) < EL_tipCurveBeginZ) then ! Conical part
             radius = (1.0D0 - (relPos(3) - EL_Hcyl)/EL_Hcone) * EL_Rcyl
          else ! Rounded tip
             radius = sqrt(abs(EL_RcTip**2 - (relPos(3) - EL_topZ + EL_RcTip)**2))
          end if
          if (dist < radius) EL_insideElectrode = .true.
       end if
    end if

  end function EL_insideElectrode

  logical function EL_around_electrode_tip(pos)
     double precision, intent(IN) :: pos(3)
     double precision             :: relPos(3)
     ! Compute distance to tip
     relPos = pos - EL_xyzPos
     relPos(3) = relPos(3) - EL_topZ
     EL_around_electrode_tip = (norm2(relPos) < 0.25d0 * EL_Hcone)
  end function EL_around_electrode_tip

  subroutine EL_setWeightFactor(ix, factor)
    integer, intent(in) :: ix
    double precision, intent(in) :: factor

    EL_weightFactors(ix) = factor
  end subroutine EL_setWeightFactor

  double precision function EL_getWeightFactor(ix)
    integer, intent(in) :: ix
    EL_getWeightFactor = EL_weightFactors(ix)
  end function EL_getWeightFactor

  subroutine EL_setCharge(ix, charge)
    integer, intent(in)          :: ix
    double precision, intent(in) :: charge

    EL_charges(ix) = charge
  end subroutine EL_setCharge

  double precision function EL_getCharge(ix)
    integer, intent(in) :: ix
    EL_getCharge = EL_charges(ix)
  end function EL_getCharge

  integer function EL_getNPoints()
    EL_getNPoints = EL_nPoints
  end function EL_getNPoints

  double precision function EL_getVoltage(time)
    double precision, intent(in) :: time

    call linearInterpolateList(elec_time_voltage(:,1), elec_time_voltage(:,2), time, EL_getVoltage)
  end function EL_getVoltage

  subroutine EL_getSurfacePoint(ix, xyz)
    integer, intent(in) :: ix
    double precision, intent(out) :: xyz(3)

    xyz = EL_surfacePoints(:, ix)
  end subroutine EL_getSurfacePoint

end module module_electrode
