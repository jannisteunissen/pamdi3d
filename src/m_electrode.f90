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
module m_electrode

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)
  
  integer :: EL_nPoints
  real(dp) :: EL_voltage, EL_xyzPos(3)
  real(dp) :: EL_topAngle, EL_topZ, EL_TransCurveStartZ, EL_TransCurveEndZ, EL_spacing
  real(dp) :: EL_tipCurveBeginZ, EL_Rcyl, EL_RcTrans, EL_Hcone, EL_Hcyl, EL_RcTip

  real(dp), allocatable :: EL_surfacePoints(:,:)
  real(dp), allocatable :: EL_charges(:)
  real(dp), allocatable :: EL_weightFactors(:)
  real(dp), allocatable :: elec_time_voltage(:,:)

  public :: EL_getVoltage, EL_getCharge, EL_getNPoints, EL_getWeightFactor, EL_getSurfacePoint
  public :: EL_setCharge, EL_setWeightFactor, EL_initialize, EL_inside_elec, EL_getRadius
  public :: EL_getBottomPos
  public :: EL_getTopPos
  public :: EL_around_elec_tip

contains

  ! Generates a list of points that lie on the electrode surface, where we want
  ! to have the electrode potential.
  subroutine EL_initialize(cfg, rng, r_max, myrank, root)
    use m_config
    use m_random
    type(CFG_t), intent(inout) :: cfg
    type(RNG_t), intent(inout) :: rng
    integer, intent(in)        :: myrank, root
    real(dp), intent(in)       :: r_max(3)
    real(dp)                   :: H_cone_eff
    real(dp), parameter        :: bottom_rel_sep = 1.0e-3_dp
    integer                    :: n_time_points

    ! Set up a basic electrode, with conical part and cylindrical part
    ! This should be moved to a subroutine, so that we can select different electrodes
    call CFG_get(cfg, "elec_voltage", EL_voltage)
    call CFG_get(cfg, "elec_h_cone", EL_Hcone)
    call CFG_get(cfg, "elec_r_cyl", EL_Rcyl)
    call CFG_get(cfg, "elec_rc_tip", EL_RcTip)
    call CFG_get(cfg, "elec_rc_trans", EL_RcTrans)
    call CFG_get(cfg, "elec_spacing", EL_spacing)

    call CFG_get(cfg, "elec_rel_pos", EL_xyzPos)
    EL_topAngle  = atan(EL_Rcyl / EL_Hcone)
    H_cone_eff = EL_Hcone - EL_Rctip * (1.0D0/sin(EL_topAngle) - 1.0D0)
    EL_Hcyl      = EL_xyzPos(3) * r_max(3) - H_cone_eff - &
         bottom_rel_sep * r_max(3)
    EL_topZ      = EL_Hcyl + H_cone_eff
    EL_xyzPos    = EL_xyzPos * r_max - (/0.0d0, 0.0d0, EL_topZ/)

    ! Determine the transitions between rounded and straight parts
    EL_TransCurveStartZ = EL_Hcyl - EL_RcTrans * tan(EL_topAngle/2.0D0)
    EL_TransCurveEndZ   = EL_TransCurveStartZ + EL_RcTrans * sin(EL_topAngle)
    EL_tipCurveBeginZ   = EL_Hcyl + EL_Hcone - EL_Rctip * (1.0D0/sin(EL_topAngle) - sin(EL_topAngle))
    call CFG_get_size(cfg, "elec_times", n_time_points)
    allocate(elec_time_voltage(n_time_points, 2))
    call CFG_get(cfg, "elec_times", elec_time_voltage(:, 1))
    call CFG_get(cfg, "elec_voltages", elec_time_voltage(:,2))

    if (myrank == root) call generate_elec_surface(rng)

  end subroutine EL_initialize

  subroutine EL_getBottomPos(coordinates)
    real(dp), intent(out) :: coordinates(3)
    coordinates = EL_xyzPos
  end subroutine EL_getBottomPos

  subroutine EL_getTopPos(coordinates)
    real(dp), intent(out) :: coordinates(3)
    coordinates = EL_xyzPos
    coordinates(3) = coordinates(3) + EL_topZ
 end subroutine EL_getTopPos

  real(dp) function EL_getRadius(zCoord)
    real(dp), intent(in) :: zCoord

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

  subroutine generate_elec_surface(rng)
    use m_units_constants
    use m_random
    type(RNG_t), intent(inout) :: rng
    integer :: ix, iy, nPoints, nVertical
    integer, allocatable :: nHorizontal(:)
    real(dp) :: height, radius, xPos, yPos
    real(dp) :: randAngle, usedSpacing

    real(dp), parameter :: min_rel_radius = 0.01d0

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
       nHorizontal(ix) = nint(2.0d0 * UC_pi * max(radius, min_rel_radius*EL_Rcyl) / usedSpacing) + 1
       height = height + usedSpacing
    end do

    ! Now store the surface points
    EL_nPoints = sum(nHorizontal)
    print *, "Electrode #points: ", EL_nPoints
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
       randAngle   = rng%unif_01() * 2.0d0 * UC_pi
       do iy = 1, nHorizontal(ix)
          xPos = radius * cos(randAngle + dble(iy-1) * 2.0D0 * UC_pi / nHorizontal(ix))
          yPos = radius * sin(randAngle + dble(iy-1) * 2.0D0 * UC_pi / nHorizontal(ix))
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
  end subroutine generate_elec_surface

  !> Return .TRUE. if the location (x,y,z) is inside the simulated electrode
  logical function EL_inside_elec(pos)
    real(dp), intent(IN) :: pos(3)
    real(dp) :: dist, radius, relPos(3)

    EL_inside_elec = .false.
    relPos = pos - EL_xyzPos

    if (relPos(3) < EL_topZ) then
       dist = norm2(relPos(1:2))
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
          if (dist < radius) EL_inside_elec = .true.
       end if
    end if

  end function EL_inside_elec

  logical function EL_around_elec_tip(pos)
     real(dp), intent(IN) :: pos(3)
     real(dp)             :: relPos(3)
     ! Compute distance to tip
     relPos = pos - EL_xyzPos
     relPos(3) = relPos(3) - EL_topZ
     EL_around_elec_tip = (norm2(relPos) < 0.25d0 * EL_Hcone)
  end function EL_around_elec_tip

  subroutine EL_setWeightFactor(ix, factor)
    integer, intent(in) :: ix
    real(dp), intent(in) :: factor

    EL_weightFactors(ix) = factor
  end subroutine EL_setWeightFactor

  real(dp) function EL_getWeightFactor(ix)
    integer, intent(in) :: ix
    EL_getWeightFactor = EL_weightFactors(ix)
  end function EL_getWeightFactor

  subroutine EL_setCharge(ix, charge)
    integer, intent(in)          :: ix
    real(dp), intent(in) :: charge

    EL_charges(ix) = charge
  end subroutine EL_setCharge

  real(dp) function EL_getCharge(ix)
    integer, intent(in) :: ix
    EL_getCharge = EL_charges(ix)
  end function EL_getCharge

  integer function EL_getNPoints()
    EL_getNPoints = EL_nPoints
  end function EL_getNPoints

  real(dp) function EL_getVoltage(time)
    use m_lookup_table
    real(dp), intent(in) :: time

    call LT_lin_interp_list(elec_time_voltage(:,1), elec_time_voltage(:,2), time, EL_getVoltage)
  end function EL_getVoltage

  subroutine EL_getSurfacePoint(ix, xyz)
    integer, intent(in) :: ix
    real(dp), intent(out) :: xyz(3)

    xyz = EL_surfacePoints(:, ix)
  end subroutine EL_getSurfacePoint

end module m_electrode
