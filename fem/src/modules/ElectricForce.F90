!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
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
!
!/******************************************************************************
! *
! *  Authors: Juha Ruokolainen, Antti Pursula
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 08 Jun 1997
! *
! *****************************************************************************/
 

!------------------------------------------------------------------------------
!>  Calculates the force due to static electric field by integrating 
!>  Maxwell stress tensor over specified boundaries. Nodal forces added later.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE StatElecForce( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE Types
  USE Lists
  USE MeshUtils
  USE Integration
  USE ElementDescription
  USE SolverUtils
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver  !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model            !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt               !< Timestep size for time dependent simulations
  LOGICAL :: TransientSimulation    !< Steady state or transient simulation

!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------

  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Variable_t), POINTER :: FieldVar, NodalVar, PotVar, DVar
  TYPE(Nodes_t) :: ElementNodes, ParentNodes
  TYPE(Element_t), POINTER   :: CurrentElement, Parent
  TYPE(ValueList_t), POINTER :: Material
  REAL(KIND=dp), ALLOCATABLE :: LocalPotential(:), Permittivity(:,:,:) 
  REAL(KIND=dp), POINTER :: Potential(:), Pwrk(:,:,:), ForceDensity(:), NodalForce(:)
  REAL(KIND=dp) :: Force(3), MomentAbout(3), Moment(3), &
    Area, PermittivityOfVacuum, sf
  INTEGER, POINTER :: NodeIndexes(:), PotentialPerm(:)
  REAL(KIND=dp), ALLOCATABLE :: NodalWeight(:)
  INTEGER :: DIM, t, pn, n, k, s, i, j, istat
  LOGICAL :: stat, CalculateMoment, ActiveBoundary, DoDisplacedBoundaries, &
      FirstTime = .TRUE., CalculateNodal, CalculateField
  CHARACTER(LEN=MAX_NAME_LEN) :: PotName, VariableName

!------------------------------------------------------------------------------

  CALL Info('StatElecForce','Computing electric force on boundaries')
  
  DIM = CoordinateSystemDimension()
  Mesh => Solver % Mesh
  NULLIFY( Pwrk )

!------------------------------------------------------------------------------
!    Figure out the mesh that potential solver is using
!------------------------------------------------------------------------------
  
  PotName = ListGetString( Solver % Values, 'Potential Field Name', stat )
  IF(.NOT. Stat) PotName = 'Potential'
  IF ( .NOT. stat ) THEN
    PotVar => VariableGet( Mesh % Variables, PotName ) 
    IF(ASSOCIATED(PotVar)) THEN
      Potential => PotVar % Values
      PotentialPerm => PotVar % Perm
    ELSE
      CALL Fatal('ElectricForce','Potential field does not exist:'//TRIM(PotName))
    END IF
  END IF

  
  CalculateNodal = ListGetLogical( Solver % Values,'Calculate Nodal Force',stat )
  NodalVar => VariableGet( Mesh % Variables, 'Electric Nodal Force' )
  IF(.NOT. ASSOCIATED(NodalVar) .AND. CalculateNodal ) THEN
    CALL VariableAddVector( Mesh % Variables, Mesh, Solver, &
        'Electric Nodal Force', Dim, Perm = PotentialPerm )
    NodalVar => VariableGet( Mesh % Variables, 'Electric Nodal Force' )
  ELSE IF( ASSOCIATED( NodalVar ) ) THEN
    CalculateNodal = .TRUE.
  END IF
  IF( CalculateNodal ) NodalForce => NodalVar % Values    
  
  CalculateField = ListGetLogical( Solver % Values,'Calculate Force Density',stat )
  IF(.NOT. CalculateNodal ) CalculateField = .TRUE.
  FieldVar => VariableGet( Mesh % Variables, 'Electric Force Density' )
  IF(.NOT. ASSOCIATED(FieldVar) .AND. CalculateField ) THEN
    CALL VariableAddVector( Mesh % Variables, Mesh, Solver, &
        'Electric Force Density', Dim, Perm = PotentialPerm )
    FieldVar => VariableGet( Mesh % Variables, 'Electric Force Density' )
  ELSE IF( ASSOCIATED( FieldVar ) ) THEN
    CalculateField = .TRUE.
  END IF
  IF( CalculateField ) ForceDensity => FieldVar % Values
      
  n = Mesh % MaxElementNodes
  ALLOCATE( ElementNodes % x( n ), &
	ElementNodes % y(n), ElementNodes % z(n) )
  ALLOCATE( ParentNodes % x( n ), &
	ParentNodes % y(n), ParentNodes % z(n) )
  ALLOCATE( LocalPotential( n ), Permittivity( 3, 3, n ) ) 
  ALLOCATE( NodalWeight( Mesh % NumberOfNodes ) )

!------------------------------------------------------------------------------

  DoDisplacedBoundaries = ListGetLogical(Solver % Values,'All Displaced Boundaries',stat)
  IF(DoDisplacedBoundaries) THEN
    VariableName = GetString( Solver % Values, 'Displacement Field Name', stat )
    IF ( .NOT. stat )  THEN
      DVar => VariableGet( Mesh % Variables, 'Displacement', .TRUE. )
    ELSE
      DVar => VariableGet( Mesh % Variables, VariableName, .TRUE. )
    END IF
  END IF

  PermittivityOfVacuum = ListGetConstReal( Model % Constants, &
      'Permittivity Of Vacuum', stat )
  IF ( .NOT. stat ) THEN
    CALL Warn( 'StatElecForce', 'Permittivity of Vacuum not given, using 1' )
    PermittivityOfVacuum = 1.0d0
  END IF

  MomentAbout(1) = ListGetConstReal( Model % Simulation,'Moment About 1', stat )
  CalculateMoment = stat 
  MomentAbout(2) = ListGetConstReal( Model % Simulation,'Moment About 2', stat )
  CalculateMoment = stat .OR. CalculateMoment
  MomentAbout(3) = ListGetConstReal( Model % Simulation,'Moment About 3', stat )
  CalculateMoment = stat .OR. CalculateMoment

  ElementNodes % x = 0.0d0
  ElementNodes % y = 0.0d0
  ElementNodes % z = 0.0d0

  ParentNodes % x = 0.0d0
  ParentNodes % y = 0.0d0
  ParentNodes % z = 0.0d0

  Area = 0.0d0
  Force  = 0.0d0
  Moment = 0.0d0
  ForceDensity = 0.0d0
  NodalWeight = 0.0_dp

  DO t = Mesh % NumberOfBulkElements + 1, &
      Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

!------------------------------------------------------------------------------
    CurrentElement => Mesh % Elements(t)
!------------------------------------------------------------------------------
!    Set the current element pointer in the model structure to 
!    reflect the element being processed
!------------------------------------------------------------------------------
    Model % CurrentElement => Mesh % Elements(t)
!------------------------------------------------------------------------------
    n = CurrentElement % TYPE % NumberOfNodes
    NodeIndexes => CurrentElement % NodeIndexes

    IF ( CurrentElement % TYPE % ElementCode == 101 ) CYCLE

    ActiveBoundary = .FALSE.
    
    IF(DoDisplacedBoundaries) THEN
      ActiveBoundary = ALL( DVar % Perm( NodeIndexes ) > 0 )
    ELSE 
      DO k=1, Model % NumberOfBCs
        IF ( Model % BCs(k) % Tag /= CurrentElement % BoundaryInfo & 
            % Constraint ) CYCLE
        IF ( ListGetLogical(Model % BCs(k) % Values, &
            'Calculate Electric Force', stat ) ) ActiveBoundary = .TRUE.
      END DO
    END IF
    
    IF(.NOT. ActiveBoundary) CYCLE
      
    ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
    ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
    ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)

!------------------------------------------------------------------------------
!     Need parent element to determine material number
!------------------------------------------------------------------------------
    Parent => CurrentElement % BoundaryInfo % Left
    
    stat = ASSOCIATED( Parent )
    IF ( stat ) THEN
      i = Parent % BodyId
      j = ListGetInteger( Model % Bodies(i) % Values, 'Material', &
          minv=1, maxv=Model % NumberOfMaterials )
      Material => Model % Materials(j) % Values
      
      CALL ListGetRealArray( Material, 'Relative Permittivity', Pwrk, n, &
          NodeIndexes, stat )
      IF ( .NOT. stat )  CALL ListGetRealArray( Material, &
          'Permittivity', Pwrk, n, NodeIndexes, stat )
    END IF
    IF ( .NOT. stat ) THEN
      Parent => CurrentElement % BoundaryInfo % Right
      IF ( .NOT. ASSOCIATED( Parent ) ) THEN
        WRITE( Message, * ) 'No permittivity found on specified boundary'
        CALL Fatal( 'StatElecForce', Message )
      END IF
      i = Parent % BodyId
      j = ListGetInteger( Model % Bodies(i) % Values, 'Material', &
          minv=1, maxv=Model % NumberOfMaterials )
      Material => Model % Materials(j) % Values
      
      CALL ListGetRealArray( Material, 'Relative Permittivity', Pwrk, n, &
          NodeIndexes, stat )
      IF ( .NOT. stat )  CALL ListGetRealArray( Material, &
          'Permittivity', Pwrk, n, NodeIndexes, stat )
      
      IF ( .NOT. stat ) THEN
        WRITE( Message, *) 'No permittivity found on specified boundary'
        CALL Fatal( 'StatElecForce', Message )
      END IF
    END IF
      
!------------------------------------------------------------------------------

    stat = ALL( PotVar % Perm( NodeIndexes ) > 0 )
    
    IF ( .NOT. stat ) THEN
      WRITE( Message, *) 'No potential available for specified boundary'
      CALL Fatal( 'StatElecForce', Message )
    END IF
      
!------------------------------------------------------------------------------
    pn = Parent % TYPE % NumberOfNodes
      
    ParentNodes % x(1:pn) = Mesh % Nodes % x(Parent % NodeIndexes)
    ParentNodes % y(1:pn) = Mesh % Nodes % y(Parent % NodeIndexes)
    ParentNodes % z(1:pn) = Mesh % Nodes % z(Parent % NodeIndexes)

!------------------------------------------------------------------------------

    Permittivity = 0.0d0
    IF ( SIZE(Pwrk,1) == 1 ) THEN
      DO i=1,3
        Permittivity( i,i,1:n ) = Pwrk( 1,1,1:n )
      END DO
    ELSEIF ( SIZE(Pwrk,2) == 1 ) THEN
      DO i=1,MIN(3,SIZE(Pwrk,1))
        Permittivity(i,i,1:n) = Pwrk(i,1,1:n)
      END DO
    ELSE
      DO i=1,MIN(3,SIZE(Pwrk,1))
        DO j=1,MIN(3,SIZE(Pwrk,2))
          Permittivity( i,j,1:n ) = Pwrk(i,j,1:n)
        END DO
      END DO
    END IF

!------------------------------------------------------------------------------
      
    LocalPotential = 0.0d0
    LocalPotential(1:pn) = Potential( PotentialPerm( Parent % NodeIndexes ) )
    
    CALL MaxwellStressTensorIntegrate( Force, Moment, Area )      
  END DO

  
  DO i= 1, Mesh % NumberOfNodes
    IF ( NodalWeight(i) > 0 ) THEN
      DO j = 1, Dim
        ForceDensity(Dim*(i-1)+j) = ForceDensity(Dim*i-Dim+j) / NodalWeight(i)
      END DO
    END IF
  END DO

!------------------------------------------------------------------------------

  CALL Info( 'StatElecForce', ' ', Level=4 )
  IF(dim == 2 ) THEN
    WRITE( Message, '("Net electric force  : ", 3ES15.6E2 )' ) Force(1:2)
  ELSE
    WRITE( Message, '("Net electric force  : ", 3ES15.6E2 )' ) Force
  END IF
  CALL Info( 'StatElecForce', Message, Level=4 )
  
  sf = SQRT(SUM( Force**2 ) )
  WRITE( Message, '("Resultant force  : ", ES15.6 )' ) sf
  CALL Info( 'StatElecForce', Message, Level=4 )

  IF(CalculateMoment) THEN
    CALL ListAddConstReal( Model % Simulation, &
        'res: Electric Moment 3', Moment(3)  )    
    CALL ListAddConstReal( Model % Simulation, &
        'res: Electric Moment 2', Moment(2)  )
    CALL ListAddConstReal( Model % Simulation, &
        'res: Electric Moment 1', Moment(1)  )
  END IF
 
  IF(DIM > 2) CALL ListAddConstReal( Model % Simulation, &
      'res: Electric Force 3', Force(3)  )
  CALL ListAddConstReal( Model % Simulation, &
      'res: Electric Force 2', Force(2)  )
  CALL ListAddConstReal( Model % Simulation, &
      'res: Electric Force 1', Force(1)  )
   
  DEALLOCATE( ElementNodes % x, ElementNodes % y, ElementNodes % z )
  DEALLOCATE( ParentNodes % x, ParentNodes % y, ParentNodes % z )
  DEALLOCATE( LocalPotential, Permittivity )
  DEALLOCATE( NodalWeight )

  ! These are some obsolete stuff but let's do this anyways
  IF( CalculateField ) THEN
    FieldVar % PrimaryMesh => Mesh
    CALL InvalidateVariable( Model % Meshes,  Mesh, 'Electric Force Density' )
    FieldVar % Valid = .TRUE.
  END IF
  IF( CalculateNodal ) THEN
    NodalVar % PrimaryMesh => Mesh
    CALL InvalidateVariable( Model % Meshes,  Mesh, 'Electric Force Density' )
    NodalVar % Valid = .TRUE.
  END IF
    
  CALL SetCurrentMesh( Model, Solver % Mesh )

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE MaxwellStressTensorIntegrate( Force, Moment, Area )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Force(3), Moment(3), Area
!------------------------------------------------------------------------------

    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    REAL(KIND=dp), POINTER :: U_Integ(:), V_Integ(:), W_Integ(:), S_Integ(:)
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),ddBasisddx(n,3,3)
    REAL(KIND=dp) :: ParentBasis(pn),ParentdBasisdx(pn,3)
    REAL(KIND=dp) :: SqrtMetric,Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3)
    REAL(KIND=dp) :: u,v,w,s, detJ, x(n), y(n), z(n)
    REAL(KIND=dp) :: Tensor(3,3), Normal(3), EField(3), DFlux(3), Lforce(3), &
      LMoment(3), Radius(3)
    REAL(KIND=dp) :: ElementWeight(n), ElementForce(3,n), xpos, ypos, zpos
    INTEGER :: N_Integ
    INTEGER :: i,j,l
    LOGICAL :: stat

    ElementWeight = 0.0d0
    ElementForce = 0.0d0

!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
    IntegStuff = GaussPoints( CurrentElement )

    U_Integ => IntegStuff % u
    V_Integ => IntegStuff % v
    W_Integ => IntegStuff % w
    S_Integ => IntegStuff % s
    N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------
!     Over integration points
!------------------------------------------------------------------------------
    DO l=1,N_Integ
!------------------------------------------------------------------------------
      u = U_Integ(l)
      v = V_Integ(l)
      w = W_Integ(l)
!------------------------------------------------------------------------------
!      Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
      stat = ElementInfo( CurrentElement, ElementNodes, u, v, w, &
          detJ, Basis, dBasisdx, ddBasisddx, .FALSE., .FALSE. )

!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
      s = 1.0d0
      IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
        xpos = SUM( ElementNodes % x(1:n) * Basis(1:n) )
        ypos = SUM( ElementNodes % y(1:n) * Basis(1:n) )
        zpos = SUM( ElementNodes % z(1:n) * Basis(1:n) )
        s = 2*PI
      END IF
         
      CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,xpos,ypos,zpos )
 
      s = s * SqrtMetric * detJ * S_Integ(l)

      Normal = Normalvector( CurrentElement,ElementNodes, u,v, .TRUE. )
!------------------------------------------------------------------------------
! Need parent element basis etc., for computing normal derivatives
! on boundary.
!------------------------------------------------------------------------------
      DO i = 1,n
        DO j = 1,pn
          IF ( CurrentElement % NodeIndexes(i) == &
              Parent % NodeIndexes(j) ) THEN
            x(i) = Parent % TYPE % NodeU(j)
            y(i) = Parent % TYPE % NodeV(j)
            z(i) = Parent % TYPE % NodeW(j)
            EXIT
          END IF
        END DO
      END DO

      u = SUM( Basis(1:n) * x(1:n) )
      v = SUM( Basis(1:n) * y(1:n) )
      w = SUM( Basis(1:n) * z(1:n) )

      stat = ElementInfo( Parent, ParentNodes, u, v, w, detJ, ParentBasis, &
          ParentdBasisdx, ddBasisddx, .FALSE., .FALSE. )

!------------------------------------------------------------------------------
      EField = 0.0d0
      DFlux = 0.0d0
      DO i=1,DIM
        EField(i) = SUM( ParentdBasisdx(1:pn,i) * LocalPotential(1:pn) )
        DO j=1,DIM
          DFlux(i) = DFlux(i)+ SUM( ParentdBasisdx(1:pn,j) * LocalPotential(1:pn) ) &
              * SUM( Permittivity(i,j,1:n) * Basis(1:n) )
        END DO
      END DO
      DFlux = PermittivityOfVacuum * DFlux

      Tensor = 0.0d0
      DO i=1, DIM
        DO j=1, DIM
          Tensor(i,j) = - DFlux(i) * EField(j)
        END DO
      END DO
      DO i=1, DIM
        Tensor(i,i) = Tensor(i,i) + SUM( DFlux * EField ) / 2.0d0
      END DO

      LForce =  MATMUL( Tensor, Normal )
      Force  = Force  + s * LForce
      Area = Area + s * SUM( Basis(1:n ) )

      IF(CalculateMoment) THEN
       Radius(1) = SUM( ElementNodes % x(1:n) * Basis(1:n) ) - MomentAbout(1)
       Radius(2) = SUM( ElementNodes % y(1:n) * Basis(1:n) ) - MomentAbout(2)
       Radius(3) = SUM( ElementNodes % z(1:n) * Basis(1:n) ) - MomentAbout(3)

        LMoment(1) = Radius(2) * LForce(3) - Radius(3) * LForce(2)
        LMoment(2) = Radius(3) * LForce(1) - Radius(1) * LForce(3)
        LMoment(3) = Radius(1) * LForce(2) - Radius(2) * LForce(1)
        Moment = Moment + s * Lmoment
      END IF

      DO i=1,n      
        ElementForce(:,i) = ElementForce(:,i) + s * LForce * Basis(i)
        ElementWeight(i) = ElementWeight(i) + s * Basis(i)
      END DO

!------------------------------------------------------------------------------
    END DO
    
    DO i=1, n
      k = PotentialPerm( CurrentElement % NodeIndexes(i))
      DO j = 1, Dim
        IF( CalculateNodal ) THEN
          NodalForce( Dim*(k-1)+j ) = NodalForce(Dim*(k-1)+j) &
              + ElementForce(j,i)
        END IF
        IF( CalculateField ) THEN
          ForceDensity( Dim*(k-1)+j ) = ForceDensity(Dim*(k-1)+j) &
              + ElementForce(j,i) 
        END IF
      END DO
      NodalWeight(k) = NodalWeight(k) + ElementWeight(i)
    END DO
  
!------------------------------------------------------------------------------
  END SUBROUTINE MaxwellStressTensorIntegrate
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE StatElecForce
!------------------------------------------------------------------------------
