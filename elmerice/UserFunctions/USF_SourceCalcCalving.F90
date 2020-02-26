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
! *  Authors: Samuel Cook
! *  Email:   
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 
! *   2018/05/30. Samuel Cook
! *****************************************************************************
!------------------------------------------------------------------------------
  !> Works out a total source (meltwater input) term for GlaDS calculations.
  !> Uses surface runoff (from RACMO, etc.) and temp residual values from Elmer.
  !> Additional keywords needed in Body Force section of SIF:
  !>
  !>   Surface Melt = Logical True or False
  !>     If True, will look for a surface melt variable to add to source
  !>
  !>   Internal Melt = Logical True or False
  !>     If True, will calculate melt due to friction heating using residuals
  !>     from the TemperateIceSolver
  !>
  !>   Surface Melt Variable Name = String "[name of variable]"
  !>     If Surface Melt = Logical True, the name of the relevant variable to
  !>     search for surface melt values in needs to be provided. Make sure this 
  !>     variable is defined throughout the vertical extent of the mesh, not
  !>     just on the surface, or it won't be found
  !>
  !>   Internal Melt Variable Name = String "[name of variable]"
  !>     The name of the variable in which the residual temperature values are
  !>     stored, if Internal Melt = Logical True
  !> 
!------------------------------------------------------------------------------
  FUNCTION SourceCalc (Model, NodeNumber, SomeVariable) RESULT(Source)
!------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE

    TYPE(Model_t) :: Model
    INTEGER :: NodeNumber
    REAL(KIND=dp) :: SomeVariable

    !----------------------------------------------------------
    TYPE(Element_t), POINTER :: CurrentElement
    TYPE(Variable_t), POINTER :: IMVar, SMVar, Weights
    TYPE(ValueList_t), POINTER :: BodyForce
    INTEGER ::  i, HPSolver
    REAL(KIND=dp) :: Source, InternalMelt, SurfaceMelt
    LOGICAL :: Found, SMLogical, IMLogical, FirstTime=.TRUE., Calving
    CHARACTER(LEN=20) :: IMVarName, SMVarName
  
    SAVE FirstTime, HPSolver, Calving
  !------------------------------------------------------------------------
    Found = .FALSE.
    SMLogical = .FALSE.
    IMLogical = .FALSE.

    CurrentElement => Model % CurrentElement
    BodyForce => GetBodyForce(CurrentElement)

    !Construct table of element areas and get some mesh parameters.
    !Only needed first time function called; subsequently saved
    IF (FirstTime) THEN
      Calving = GetLogical(Model % Simulation, 'Calving', Found)
      IF(.NOT. Found) Calving = .FALSE.
      DO i=1, Model % NumberOfSolvers
        IF(Model % Solvers(i) % Variable % Name == 'hydraulic potential') THEN
          HPSolver = i
          EXIT
        END IF
      END DO

      !Calculate nodal weights for scaling purposes
      CALL CalculateNodalWeights(Model % Solvers(HPSolver), .FALSE., VarName='Weights')
      FirstTime = .FALSE.
    END IF !FirstTime

    !Get internal melt variable - in Elmer, will be the temperature residual
    !Name and logical defined in Body Force section of sif
    IMLogical = GetLogical(BodyForce,'Internal Melt',Found)
    IF(.NOT. Found) CALL Fatal('SourceCalc', "You've not specified if you want internal melt, silly")
    IF (IMLogical) THEN
      IMVarName = GetString(BodyForce,'Internal Melt Variable Name',Found)
      IF(Calving) THEN
        IMVar => VariableGet(Model % Solvers(HPSolver) % Mesh % Variables,&
               IMVarName, ThisOnly = .TRUE., UnfoundFatal = .TRUE.)
        Weights => VariableGet(Model % Solvers(HPSolver) % Mesh % Variables,&
               'Weights', ThisOnly = .TRUE., UnfoundFatal = .TRUE.)
      ELSE
        IMVar => VariableGet(Model % Variables, IMVarName, ThisOnly = .TRUE., UnfoundFatal = .TRUE.)
        Weights => VariableGet(Model % Variables, 'Weights', ThisOnly = .TRUE., UnfoundFatal = .TRUE.)
      END IF
      Found = .FALSE.
    END IF

    !Get Surface Melt - Logical and name set in Body Force in sif
    !Needs to be surface runoff values in m/year
    !Make sure is defined on basal nodes/throughout glacier, rather than just
    !surface
    SMVarName = GetString(BodyForce,'Surface Melt Variable Name',Found)
    SMLogical = GetLogical(BodyForce,'Surface Melt',Found)
    IF(.NOT. Found) CALL Fatal('SourceCalc', "You've not defined surface melt, you muppet")
    IF (SMLogical) THEN
      IF(Calving) THEN
        SMVar => VariableGet(Model % Solvers(HPSolver) % Mesh % Variables,&
                 SMVarName, ThisOnly = .TRUE., UnfoundFatal = .TRUE.)
      ELSE
        SMVar => VariableGet(Model % Variables, SMVarName, ThisOnly = .TRUE., UnfoundFatal = .TRUE.)
      END IF
      Found = .FALSE.
    END IF

    IF(IMLogical) THEN
      !If Temp Residual overall positive, no internal melting
      !Temp Res is essentially difference between temperature and upper or lower
      !limit
      !If temp exceeded upper limit, then is negative; if below lower limit,
      !is positive
      !Here, should therefore all be negative, unless glacier dropped below
      !absolute zero....
      IF(IMVar % Values(IMVar % Perm(NodeNumber))>=0.0) THEN
        InternalMelt = 0.0
      ELSE
        !Latent heat of fusion of water is 333.55 J/g, so dividing by that gives
        ! g of ice melted.
        !TempRes in MJ, though (probably), so dividing by 333.55 gives Mg of ice
        ! melted
        !1 Mg is 1 t, which is 1000 kg, so 1000 l, so 1 m3 (all per year), so
        !that's it
        !Also need to divide by element area to get m
        InternalMelt = (ABS(IMVar % Values(IMVar % Perm(NodeNumber)))/Weights % Values(Weights % Perm(NodeNumber)))/333.55
      END IF
    ELSE
      InternalMelt = 0.0
    END IF

    IF(SMLogical) THEN
      IF(Calving) THEN
        IF(SMVar % Values(SMVar % Perm(NodeNumber))>0.0) THEN
          SurfaceMelt = SMVar % Values(SMVar % Perm(NodeNumber))
        ELSE
          SurfaceMelt = 0.0
        END IF
      ELSE
        IF(SMVar % Values(SMVar % Perm(NodeNumber))>0.0) THEN
          SurfaceMelt = SMVar % Values(SMVar % Perm(NodeNumber))
        ELSE
          SurfaceMelt = 0.0
        END IF
      END IF
    ELSE
      SurfaceMelt = 0.0
    END IF

    !Work out total
    Source = InternalMelt + SurfaceMelt

!------------------------------------------------------------------------------
  END FUNCTION SourceCalc
!------------------------------------------------------------------------------
