!/*****************************************************************************/
! *
! *  Elmer/Ice, a glaciological add-on to Elmer
! *  http://elmerice.elmerfem.org
! *
! * 
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! * 
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the 
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/
! ******************************************************************************
! *
! *  Authors: Denis Cohen, Thomas Zwinger, and Peter Raback
! *  Email:   denis.cohen@gmail.com
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 
! *   2016/02/16. Denis Cohen
! *****************************************************************************
! Contains three functions to model permafrost in rock beneath ice sheet:
! - GetEnthalpy
! - GetDensity
! - GetConductivity
! Permafrost is modeled as a three-component mixture (rock + ice + water) with
! an effective heat capacity that depends on water content to mimick phase
! change near the freezing point.
! Two models for water content as a function of temperature are available:
! - power-law (e.g. Cutler et al, 2000, A numerical injvestigation of 
!   ice-lobe-permafrost interaction around the southern Laurentide ice sheet, 
!   J. Glaciol., 46, 311-325)
! - exponential (e.g. Willeit and Ganopolski, 2015, Coupled Northern Hemisphere 
!   permafrost-ice-sheet evolution over the last glacial cycle, Clim. Past, 11,
!   1165-1180)
! *****************************************************************************

!==============================================================================
FUNCTION PermafrostEnthalpy(Model, Node, Temp) RESULT(enthalpy)
!==============================================================================

  USE DefUtils

  IMPLICIT None

  TYPE(Model_t) :: Model
  TYPE(ValueList_t), POINTER :: Material

  INTEGER :: Node
  REAL(KIND=dp) :: Temp
  REAL(KIND=dp) :: enthalpy

  ! Local variables

  TYPE(Variable_t), POINTER :: DepthVar, DepthVar2
  REAL(KIND=dp) :: Depth, Depth2

  REAL(KIND=dp) :: por                 ! Porosity
  REAL(KIND=dp) :: porscale            ! Porosity length scale
  REAL(KIND=dp) :: pordepth            ! Porosity as a function of depth
  REAL(KIND=dp) :: phir, phiw, phii    ! Volume fractions
  ! r = rock, w = water, i = ice
  REAL(KIND=dp) :: rhor, rhow, rhoi    ! Densities 
  REAL(KIND=dp) :: Cr, Cw, Ci          ! Heat capacities
  REAL(KIND=dp) :: L                   ! Latent heat
  REAL(KIND=dp) :: Tpmp                ! Melting point temperature of ice
  REAL(KIND=dp) :: a, b, Tstar, dT     ! Params for powerlaw/exponential model
  REAL(KIND=dp) :: fw                  ! Params for exponential model
  REAL(KIND=dp) :: iceDepth, pice, prock, press ! For computing pressures

  CHARACTER(LEN=MAX_NAME_LEN) :: PermafrostModel

  LOGICAL :: FirstTime = .TRUE.
  LOGICAL :: Found

  SAVE por, porscale
  SAVE rhor, rhow, rhoi
  SAVE Cr, Cw, Ci, L
  SAVE a, b, dT
  SAVE FirstTime, PermafrostModel

  !-----------------------------------------------
  ! Read parameters from sif file Material section
  !-----------------------------------------------
  IF (FirstTime) THEN
    Material => GetMaterial(Model % CurrentElement)

    IF (.NOT.ASSOCIATED(Material)) THEN
      CALL FATAL('Permafrost', 'No Material found')
    END IF

    PermafrostModel = GetString( Material, 'Permafrost Model', Found )
    IF (.NOT. Found) THEN
      CALL FATAL('Permafrost', 'Cound not find Permafrost Model')
    ENDIF

    por = GetCReal( Material, 'Permafrost Porosity', Found )
    IF (.NOT. Found) THEN
      CALL FATAL('Permafrost', 'Cound not find Porosity')
    ENDIF

    porscale = GetCReal( Material, 'Permafrost Porosity Depth Scale', Found )
    IF (.NOT. Found) THEN
      porscale = -1.0 ! Negative value means no depth dependence for por
    ENDIF

    L = GetCReal( Material, 'Latent Heat', Found )
    IF (.NOT. Found) THEN
      CALL FATAL('Permafrost', 'Cound not find Latent Heat')
    ENDIF

    !--- Rock parameters ---
    rhor = GetCReal( Material, 'Permafrost Density Rock', Found )
    IF (.NOT. Found) THEN
      CALL FATAL('Permafrost', 'Cound not find Permafrost Density Rock')
    ENDIF

    Cr = GetCReal( Material, 'Permafrost Heat Capacity Rock', Found )
    IF (.NOT. Found) THEN
      CALL FATAL('Permafrost', 'Cound not find Permafrost Heat Capacity Rock')
    ENDIF

    !--- Water parameters ---
    rhow = GetCReal( Material, 'Permafrost Density Water', Found )
    IF (.NOT. Found) THEN
      CALL FATAL('Permafrost', 'Cound not find Permafrost Density Water')
    ENDIF
    write (*,*) 'rhow = ', rhow

    Cw = GetCReal( Material, 'Permafrost Heat Capacity Water', Found )
    IF (.NOT. Found) THEN
      CALL FATAL('Permafrost', 'Cound not find Permafrost Heat Capacity Water')
    ENDIF
    write (*,*) 'Cw = ', Cw

    !--- Ice parameters ---
    rhoi = GetCReal( Material, 'Permafrost Density Ice', Found )
    IF (.NOT. Found) THEN
      CALL FATAL('Permafrost', 'Cound not find Permafrost Density Ice')
    ENDIF

    Ci = GetCReal( Material, 'Permafrost Heat Capacity Ice', Found )
    IF (.NOT. Found) THEN
      CALL FATAL('Permafrost', 'Cound not find Permafrost Heat Capacity Ice')
    ENDIF

    !--- Power Law model parameters ---
    IF (TRIM(PermafrostModel) .EQ. "power law") THEN
      a = GetCReal( Material, 'Permafrost Power Law Factor', Found )
      IF (.NOT. Found) THEN
        CALL FATAL('Permafrost', 'Cound not find Permafrost Power Law Factor')
      ENDIF
      b = GetCReal( Material, 'Permafrost Power Law Exponent', Found )
      IF (.NOT. Found) THEN
        CALL FATAL('Permafrost', 'Cound not find Permafrost Power Law Exponent')
      ENDIF
      dT = GetCReal( Material, 'Permafrost Power law Temperature Offset', Found )
      IF (.NOT. Found) THEN
        CALL FATAL('Permafrost', 'Cound not find Permafrost Power law Temperature Offset')
      ENDIF

      !--- Exponential model parameters ---
    ELSE IF (TRIM(PermafrostModel) .EQ. "exponential") THEN
      a = GetCReal( Material, 'Permafrost Exponential Temperature Interval', Found )
      IF (.NOT. Found) THEN
        CALL FATAL('Permafrost', 'Cound not find Permafrost Exponential Temperature Interval')
      ENDIF
    ELSE
      CALL FATAL('Permafrost', 'Unknown Permafrost Model')
    ENDIF

    FirstTime = .FALSE.
  ENDIF

  !-----------------------------------------------
  ! Get the depth of lower layer
  !-----------------------------------------------
  DepthVar => VariableGet(Model % Mesh % Variables, "lower depth")
  IF ( ASSOCIATED(DepthVar) ) THEN
    Depth = DepthVar % Values ( DepthVar % Perm(Node) )
  ELSE 
    Depth = 0.0
  END IF

  !-----------------------------------------------
  ! Get the total depth below all layers
  !-----------------------------------------------
  DepthVar2 => VariableGet(Model % Mesh % Variables, "depth")
  IF ( ASSOCIATED(DepthVar2) ) THEN
    Depth2 = DepthVar2 % Values ( DepthVar2 % Perm(Node) )
  ELSE
    Depth2 = 0.0
  END IF

  ! In case there is no lower layer, use the only layer
  if (Depth == 0.0) Depth = Depth2

  WRITE(MESSAGE,*) Node, Temp, Depth, Depth2
  CALL Info('Permafrost', MESSAGE, level=11)

  !-------------------
  ! Start calculations
  !-------------------
  IF (porscale .LE. 0.0) THEN
    pordepth = por
  ELSE
    pordepth = por * EXP(-Depth/porscale)
  ENDIF

  phir = 1 - pordepth

  !----------------------
  ! Compute Tpmp at depth
  !----------------------
  iceDepth = Depth2 - Depth
  pice = iceDepth * rhoi * 9.81
  prock = Depth * rhor * 9.81
  press = pice + prock
 !Tpmp = 273.15 ! Should be changed to depend on pressure
  Tpmp = 273.15 - 9.8E-08*press

  IF (TRIM(PermafrostModel) .EQ. "power law") THEN
    Tstar = Tpmp - Temp
    IF (Tstar <= dT) THEN
      phiw =  pordepth
    ELSE 
      phiw = (rhor * (1 - pordepth)/ rhow) * a * (Tstar)**b
      !phiw = (rhor / rhow) * a * (Tstar)**b
    ENDIF

  ELSE IF (TRIM(PermafrostModel) .EQ. "exponential") THEN
    IF (Temp > Tpmp) then 
      fw = 1.0
    ELSE
      fw = EXP(-((Temp - Tpmp)/a)**2)
    ENDIF
    phiw =  pordepth * fw
  ELSE
    CALL FATAL('Permafrost', 'Unknown Permafrost Model')

  ENDIF

  ! Check that phiw does not exceed porosity
  IF (phiw > pordepth) THEN
    phiw = pordepth
  ENDIF

  enthalpy = (phir*rhor*Cr + (pordepth-phiw)*rhoi*Ci + phiw*rhow*Cw)*(Temp) + phiw*rhow*L
  !WRITE (*,'(F12.10,1X,F10.5,1X,F10.5,1X,F12.10,1X,E17.10)') pordepth, Temp, Tstar, phiw, enthalpy

END FUNCTION PermafrostEnthalpy

!==============================================================================
FUNCTION PermafrostDensity(Model, Node, Temp) RESULT(Dens)
!==============================================================================

  USE DefUtils

  IMPLICIT None

  TYPE(Model_t) :: Model
  TYPE(ValueList_t), POINTER :: Material

  INTEGER :: Node
  REAL(KIND=dp) :: Temp
  REAL(KIND=dp) :: Dens

  ! Local variables

  TYPE(Variable_t), POINTER :: DepthVar, DepthVar2
  REAL(KIND=dp) :: Depth, Depth2

  REAL(KIND=dp) :: por                 ! Porosity
  REAL(KIND=dp) :: porscale            ! Porosity length scale
  REAL(KIND=dp) :: pordepth            ! Porosity as a function of depth
  REAL(KIND=dp) :: phir, phiw, phii    ! Volume fractions
  ! r = rock, w = water, i = ice
  REAL(KIND=dp) :: rhor, rhow, rhoi    ! Densities 
  REAL(KIND=dp) :: Tpmp                ! Melting point temperature of ice
  REAL(KIND=dp) :: a, b, Tstar, dT     ! Params for powerlaw/exponential model
  REAL(KIND=dp) :: fw                  ! Params for exponential model
  REAL(KIND=dp) :: iceDepth, pice, prock, press

  CHARACTER(LEN=MAX_NAME_LEN) :: PermafrostModel

  LOGICAL :: FirstTime = .TRUE.
  LOGICAL :: Found

  SAVE por, porscale
  SAVE rhor, rhow, rhoi
  SAVE a, b, dT
  SAVE FirstTime, PermafrostModel

  !-----------------------------------------------
  ! Read parameters from sif file Material section
  !-----------------------------------------------
  IF (FirstTime) THEN
    Material => GetMaterial(Model % CurrentElement)

    IF (.NOT.ASSOCIATED(Material)) THEN
      CALL FATAL('Permafrost',"No Material found")
    END IF

    PermafrostModel = GetString( Material, 'Permafrost Model', Found )
    IF (.NOT. Found) THEN
      CALL FATAL('Permafrost', 'Cound not find Permafrost Model')
    ENDIF

    por = GetCReal( Material, 'Permafrost Porosity', Found )
    IF (.NOT. Found) THEN
      CALL FATAL('Permafrost', 'Cound not find Porosity')
    ENDIF

    porscale = GetCReal( Material, 'Permafrost Porosity Depth Scale', Found )
    IF (.NOT. Found) THEN
      porscale = -1.0 ! Negative value means no depth depedence for por
    ENDIF

    rhor = GetCReal( Material, 'Permafrost Density Rock', Found )
    IF (.NOT. Found) THEN
      CALL FATAL('Permafrost', 'Cound not find Permafrost Density Rock')
    ENDIF

    rhow = GetCReal( Material, 'Permafrost Density Water', Found )
    IF (.NOT. Found) THEN
      CALL FATAL('Permafrost', 'Cound not find Permafrost Density Water')
    ENDIF

    rhoi = GetCReal( Material, 'Permafrost Density Ice', Found )
    IF (.NOT. Found) THEN
      CALL FATAL('Permafrost', 'Cound not find Permafrost Density Ice')
    ENDIF

    !--- Power Law model parameters ---
    IF (TRIM(PermafrostModel) .EQ. "power law") THEN
      a = GetCReal( Material, 'Permafrost Power Law Factor', Found )
      IF (.NOT. Found) THEN
        CALL FATAL('Permafrost', 'Cound not find Permafrost Power Law Factor')
      ENDIF
      b = GetCReal( Material, 'Permafrost Power Law Exponent', Found )
      IF (.NOT. Found) THEN
        CALL FATAL('Permafrost', 'Cound not find Permafrost Power Law Exponent')
      ENDIF
      dT = GetCReal( Material, 'Permafrost Power law Temperature Offset', Found )
      IF (.NOT. Found) THEN
        CALL FATAL('Permafrost', 'Cound not find Permafrost Power law Temperature Offset')
      ENDIF

      !--- Exponential model parameters ---
    ELSE IF (TRIM(PermafrostModel) .EQ. "exponential") THEN
      a = GetCReal( Material, 'Permafrost Exponential Temperature Interval', Found )
      IF (.NOT. Found) THEN
        CALL FATAL('Permafrost', 'Cound not find Permafrost Exponential Temperature Interval')
      ENDIF
    ELSE
      CALL FATAL('Permafrost', 'Unknown Permafrost Model')
    ENDIF

    FirstTime = .FALSE.
  ENDIF

  !-----------------------------------------------
  ! Get the depth of lower layer
  !-----------------------------------------------
  DepthVar => VariableGet(Model % Mesh % Variables, "lower depth")
  IF ( ASSOCIATED(DepthVar) ) THEN
    Depth = DepthVar % Values ( DepthVar % Perm(Node) )
  ELSE 
    Depth = 0.0
  END IF

  !-----------------------------------------------
  ! Get the total depth below all layers
  !-----------------------------------------------
  DepthVar2 => VariableGet(Model % Mesh % Variables, "depth")
  IF ( ASSOCIATED(DepthVar2) ) THEN
    Depth2 = DepthVar2 % Values ( DepthVar2 % Perm(Node) )
  ELSE
    Depth2 = 0.0
  END IF

  ! In case there is no lower layer, use the only layer
  if (Depth == 0.0) Depth = Depth2

  WRITE(MESSAGE,*) Node, Temp, Depth, Depth2
  CALL Info('Permafrost', MESSAGE, level=11)

  !-------------------
  ! Start calculations
  !-------------------
  IF (porscale .LE. 0.0) THEN
    pordepth = por
  ELSE
    pordepth = por * EXP(-Depth/porscale)
  ENDIF

  phir = 1 - pordepth

  !----------------------
  ! Compute Tpmp at depth
  !----------------------
  iceDepth = Depth2 - Depth
  pice = iceDepth * rhoi * 9.81
  prock = Depth * rhor * 9.81
  press = pice + prock
 !Tpmp = 273.15 ! Should be changed to depend on pressure
  Tpmp = 273.15 - 9.8E-08*press

  IF (TRIM(PermafrostModel) .EQ. "power law") THEN
    Tstar = Tpmp - Temp
    IF (Tstar <= dT) THEN
      phiw =  pordepth
    ELSE 
      phiw = (rhor * (1 - pordepth)/ rhow) * a * (Tstar)**b
      !phiw = (rhor / rhow) * a * (Tstar)**b
    ENDIF

  ELSE IF (TRIM(PermafrostModel) .EQ. "exponential") THEN
    IF (Temp > Tpmp) then
      fw = 1.0
    ELSE
      fw = EXP(-((Temp - Tpmp)/a)**2)
    ENDIF
    phiw =  pordepth * fw
  ELSE
    CALL FATAL('Permafrost', 'Unknown Permafrost Model')

  ENDIF

  ! Check that phiw does not exceed porosity
  IF (phiw > pordepth) THEN
    phiw = pordepth
  ENDIF

  Dens = phir*rhor + (pordepth-phiw)*rhoi + phiw*rhow

END FUNCTION PermafrostDensity

!==============================================================================
FUNCTION PermafrostConductivity(Model, Node, Temp) RESULT(Cond)
!==============================================================================

  USE DefUtils

  IMPLICIT None

  TYPE(Model_t) :: Model
  TYPE(ValueList_t), POINTER :: Material

  INTEGER :: Node
  REAL(KIND=dp) :: Temp
  REAL(KIND=dp) :: Cond

  ! Local variables

  TYPE(Variable_t), POINTER :: DepthVar, DepthVar2
  REAL(KIND=dp) :: Depth, Depth2

  REAL(KIND=dp) :: por                 ! Porosity
  REAL(KIND=dp) :: porscale            ! Porosity length scale
  REAL(KIND=dp) :: pordepth            ! Porosity as a function of depth
  REAL(KIND=dp) :: phir, phiw, phii    ! Volume fractions
  ! r = rock, w = water, i = ice
  REAL(KIND=dp) :: rhoi, rhor, rhow          ! Densities 
  REAL(KIND=dp) :: Kr, Kw, Ki          ! Heat conductivities
  REAL(KIND=dp) :: L                   ! Latent heat
  REAL(KIND=dp) :: Tpmp                ! Melting point temperature of ice
  REAL(KIND=dp) :: a, b, Tstar, dT     ! Params for powerlaw/exponential model
  REAL(KIND=dp) :: fw                  ! Params for exponential model
  REAL(KIND=dp) :: iceDepth, pice, prock, press

  CHARACTER(LEN=MAX_NAME_LEN) :: PermafrostModel

  LOGICAL :: FirstTime = .TRUE.
  LOGICAL :: Found

  SAVE por, porscale
  SAVE rhoi, rhor, rhow
  SAVE Kr, Kw, Ki, L
  SAVE a, b, dT
  SAVE FirstTime, PermafrostModel

  !-----------------------------------------------
  ! Read parameters from sif file Material section
  !-----------------------------------------------
  IF (FirstTime) THEN
    Material => GetMaterial(Model % CurrentElement)

    IF (.NOT.ASSOCIATED(Material)) THEN
      CALL FATAL('Permafrost',"No Material found")
    END IF

    PermafrostModel = GetString( Material, 'Permafrost Model', Found )
    IF (.NOT. Found) THEN
      CALL FATAL('Permafrost', 'Cound not find Permafrost Model')
    ENDIF

    por = GetCReal( Material, 'Permafrost Porosity', Found )
    IF (.NOT. Found) THEN
      CALL FATAL('Permafrost', 'Cound not find Porosity')
    ENDIF

    porscale = GetCReal( Material, 'Permafrost Porosity Depth Scale', Found )
    IF (.NOT. Found) THEN
      porscale = -1.0 ! Negative value means no depth depedence for por
    ENDIF

    rhor = GetCReal( Material, 'Permafrost Density Rock', Found )
    IF (.NOT. Found) THEN
      CALL FATAL('Permafrost', 'Cound not find Permafrost Density Rock')
    ENDIF

    rhow = GetCReal( Material, 'Permafrost Density Water', Found )
    IF (.NOT. Found) THEN
      CALL FATAL('Permafrost', 'Cound not find Permafrost Density Water')
    ENDIF

    rhoi = GetCReal( Material, 'Permafrost Density Ice', Found )
    IF (.NOT. Found) THEN
      CALL FATAL('Permafrost', 'Cound not find Permafrost Density Ice')
    ENDIF

    Kr = GetCReal( Material, 'Permafrost Heat Conductivity Rock', Found )
    IF (.NOT. Found) THEN
      CALL FATAL('Permafrost', 'Cound not find Permafrost Heat Conductivity Rock')
    ENDIF

    Kw = GetCReal( Material, 'Permafrost Heat Conductivity Water', Found )
    IF (.NOT. Found) THEN
      CALL FATAL('Permafrost', 'Cound not find Permafrost Heat Conductivity Water')
    ENDIF
    write (*,*) 'Kw = ', Kw

    Ki = GetCReal( Material, 'Permafrost Heat Conductivity Ice', Found )
    IF (.NOT. Found) THEN
      CALL FATAL('Permafrost', 'Cound not find Permafrost Heat Conductivity Ice')
    ENDIF

    !--- Power Law model parameters ---
    IF (TRIM(PermafrostModel) .EQ. "power law") THEN
      a = GetCReal( Material, 'Permafrost Power Law Factor', Found )
      IF (.NOT. Found) THEN
        CALL FATAL('Permafrost', 'Cound not find Permafrost Power Law Factor')
      ENDIF
      b = GetCReal( Material, 'Permafrost Power Law Exponent', Found )
      IF (.NOT. Found) THEN
        CALL FATAL('Permafrost', 'Cound not find Permafrost Power Law Exponent')
      ENDIF
      dT = GetCReal( Material, 'Permafrost Power law Temperature Offset', Found )
      IF (.NOT. Found) THEN
        CALL FATAL('Permafrost', 'Cound not find Permafrost Power law Temperature Offset')
      ENDIF

      !--- Exponential model parameters ---
    ELSE IF (TRIM(PermafrostModel) .EQ. "exponential") THEN
      a = GetCReal( Material, 'Permafrost Exponential Temperature Interval', Found )
      IF (.NOT. Found) THEN
        CALL FATAL('Permafrost', 'Cound not find Permafrost Exponential Temperature Interval')
      ENDIF
    ELSE
      CALL FATAL('Permafrost', 'Unknown Permafrost Model')
    ENDIF

    FirstTime = .FALSE.
  ENDIF

  !-----------------------------------------------
  ! Get the depth of lower layer
  !-----------------------------------------------
  DepthVar => VariableGet(Model % Mesh % Variables, "lower depth")
  IF ( ASSOCIATED(DepthVar) ) THEN
    Depth = DepthVar % Values ( DepthVar % Perm(Node) )
  ELSE 
    Depth = 0.0
  END IF

  !-----------------------------------------------
  ! Get the total depth below all layers
  !-----------------------------------------------
  DepthVar2 => VariableGet(Model % Mesh % Variables, "depth")
  IF ( ASSOCIATED(DepthVar2) ) THEN
    Depth2 = DepthVar2 % Values ( DepthVar2 % Perm(Node) )
  ELSE
    Depth2 = 0.0
  END IF

  ! In case there is no lower layer, use the only layer
  if (Depth == 0.0) Depth = Depth2

  WRITE(MESSAGE,*) Node, Temp, Depth, Depth2
  CALL Info('Permafrost', MESSAGE, level=11)

  !-------------------
  ! Start calculations
  !-------------------
  IF (porscale .LE. 0.0) THEN
    pordepth = por
  ELSE
    pordepth = por * EXP(-Depth/porscale)
  ENDIF

  phir = 1 - pordepth

  !----------------------
  ! Compute Tpmp at depth
  !----------------------
  iceDepth = Depth2 - Depth
  pice = iceDepth * rhoi * 9.81
  prock = Depth * rhor * 9.81
  press = pice + prock
 !Tpmp = 273.15 ! Should be changed to depend on pressure
  Tpmp = 273.15 - 9.8E-08*press

  IF (TRIM(PermafrostModel) .EQ. "power law") THEN
    Tstar = Tpmp - Temp
    IF (Tstar <= dT) THEN
      phiw =  pordepth
    ELSE 
      phiw = (rhor * (1 - pordepth)/ rhow) * a * (Tstar)**b
      !phiw = (rhor / rhow) * a * (Tstar)**b
    ENDIF

  ELSE IF (TRIM(PermafrostModel) .EQ. "exponential") THEN
    IF (Temp > Tpmp) then
      fw = 1.0
    ELSE
      fw = EXP(-((Temp - Tpmp)/a)**2)
    ENDIF
    phiw =  pordepth * fw
  ELSE
    CALL FATAL('Permafrost', 'Unknown Permafrost Model')

  ENDIF

  ! Check that phiw does not exceed porosity
  IF (phiw > pordepth) THEN
    phiw = pordepth
  ENDIF

  Cond = Kr**phir * Ki**(pordepth-phiw) * Kw**phiw

END FUNCTION PermafrostConductivity
