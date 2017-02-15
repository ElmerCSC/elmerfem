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
! *  Authors: Olivier Gagliardini
! *  Email:   
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 
! *  Date Modifications:
! *          2009/03/17 GAG generalised for 3D problems and non-linear elements
! *****************************************************************************
!> 
!> 
!> 
!>  
! *****************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>  
!>    Buoyancy user functions
!>  
!>   SeaPressure returns the pressure -rho_w. g . (Hsl - S - b.Ns.dt)
!>    where Ns = sqrt(1 + (dS/dx)^2 + (dS/dy)^2)
!>          b = melting below the shelf (normal flux)
!>   External Pressure = ....
!>  
!>   SeaSpring returns the viscous spring rho_w.g.Ns.dt 
!>   Normal-Tangential set to true
!>   Slip Coefficient 1 = ....
!>  
!>  Adapted from SeaPressureOld. Allows accretion melting below the shelf.
!>   2009/03/17 GAG
!>   - generalised for 3D problems and non-linear elements
FUNCTION SeaPressure ( Model, nodenumber, y) RESULT(pw)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   TYPE(Solver_t):: Solver
   TYPE(Nodes_t), SAVE :: Nodes
   TYPE(variable_t), POINTER :: Timevar
   TYPE(Element_t), POINTER ::  BoundaryElement, BCElement, CurElement, ParentElement
   TYPE(ValueList_t), POINTER :: BC, material, ParentMaterial, BodyForce
   INTEGER :: NBoundary, NParent, BoundaryElementNode, ParentElementNode, body_id, other_body_id, material_id
   INTEGER :: nodenumber, NumberOfNodesOnBoundary 
   INTEGER, ALLOCATABLE :: NodeOnBoundary(:)
   INTEGER :: Nn, i, j, p, n, Nmax, bf_id, DIM, bf_id_FS 
   REAL(KIND=dp) :: y, pw, t, told, dt, Bu, Bv
   REAL(KIND=dp) :: Zsl, rhow, gravity
   REAL(KIND=dp), ALLOCATABLE :: S(:), auxReal(:), Ns(:),  a_perp(:), SourceFunc(:), normal(:,:)
   LOGICAL :: FirstTime = .TRUE., NewTime, GotIt, ComputeS,  NormalFlux = .TRUE., UnFoundFatal=.TRUE.
   CHARACTER(LEN=MAX_NAME_LEN)  :: BottomSurfaceName
       
   SAVE told, FirstTime, NewTime, Nn, dt, Ns, Bodyforce, DIM
   SAVE S, rhow, gravity, Zsl, auxReal, NormalFlux, a_perp, SourceFunc
   SAVE NumberOfNodesOnBoundary, NodeOnBoundary, normal
   SAVE BottomSurfaceName, bf_id_FS 
   


   Timevar => VariableGet( Model % Variables,'Time')
   t = TimeVar % Values(1)
   dt = Model % Solver % dt 


   IF (FirstTime) THEN
      FirstTime = .FALSE.
      NewTime = .TRUE.
      told = t
      ALLOCATE( NodeOnBoundary( Model % Mesh % NumberOfNodes ))
      n = Model % MaxElementNodes 
      ALLOCATE( auxReal(n), SourceFunc(n) )
      DIM = CoordinateSystemDimension()

      rhow = GetConstReal( Model % Constants, 'Water Density', GotIt )
      IF (.NOT.GotIt) THEN
         WRITE(Message,'(A)') 'Variable Water Density not found. &
              &Setting to 1.03225e-18'
         CALL INFO('SeaPressure', Message, level=20)
         rhow = 1.03225e-18_dp
      ELSE
         WRITE(Message,'(A,F10.4)') 'Water Density = ', rhow
         CALL INFO('SeaPressure', Message , level = 20)
      END IF
      
!-----------------------------------------------------------------
! Is there basal melt to account for          
!-----------------------------------------------------------------
      NormalFlux = GetLogical( Model % Constants, 'Buoyancy Use Basal Melt', GotIt )    
      IF (.NOT.GotIt) THEN
         WRITE(Message,'(A)') 'Variable Buoyancy Use Basal Melt not found. &
              &Setting to FALSE'
         CALL INFO('SeaPressure', Message, level=20)
         NormalFlux = .FALSE.
      ELSE
         WRITE(Message,*) 'Buoyancy Use Basal Melt = ', NormalFlux
         CALL INFO('SeaPressure', Message , level = 20)
      END IF

!-----------------------------------------------------------------
! Find the bf_id of the FS solvers if NormalFlux is TRUE
! and read the name of the Bottom Free surface solver
!-----------------------------------------------------------------
      IF (NormalFlux) THEN
        BottomSurfaceName = GetString( Model % Constants, 'Bottom Surface Name', GotIt )    
        IF (.NOT.GotIt) THEN
           WRITE(Message,'(A)') 'Variable Bottom SUrface Name not found. &
                &Setting to DyBottom'
           CALL INFO('SeaPressure', Message, level=20)
           BottomSurfaceName = 'DyBottom'
        ELSE
           WRITE(Message,*)'Bottom Surface Name = ',TRIM(BottomSurfaceName)
           CALL INFO('SeaPressure', Message , level = 20)
        END IF

        bf_id_FS = -1
        DO bf_id=1,Model % NumberOFBodyForces
           IF( ListCheckPresent( Model % BodyForces(bf_id) % Values, &
                TRIM(BottomSurfaceName)//' Accumulation') ) bf_id_FS = bf_id 
        END DO

        IF (bf_id_FS<0) THEN 
           CALL FATAL('Sea Pressure','No Basal Melt found in any body Force')
        END IF 
      END IF

!-----------------------------------------------------------------
! Read the gravity 
!-----------------------------------------------------------------
      bf_id = ListGetInteger( Model % Bodies(1) % Values,&
           'Body Force',  minv=1, maxv=Model % NumberOFMaterials)

      IF (DIM==2) THEN    !NOT the same keyword if the shape factor is used
        gravity = -GetConstReal( Model % BodyForces(bf_id) % Values, 'Lateral Friction Gravity 2',Gotit)
        IF (.NOT. Gotit) THEN
          gravity = -GetConstReal( Model % BodyForces(bf_id) % Values, 'Flow BodyForce 2',Gotit)
        END IF
      ELSE
        gravity = -GetConstReal( Model % BodyForces(bf_id) % Values, 'Flow BodyForce 3',Gotit)
      END IF
 
!-----------------------------------------------------------------
! Number of nodes concerned (needed to compute Ns from the normal)
!-----------------------------------------------------------------

      CurElement => Model % CurrentElement
      NodeOnBoundary = 0
      NumberOfNodesOnBoundary = 0
      DO p = 1, Model % NumberOfBoundaryElements
         BCElement => GetBoundaryElement(p) 
         BC => GetBC( BCElement ) 
         IF( GetElementFamily(BCElement) == 1 ) CYCLE
         ComputeS = GetLogical( BC, 'Compute Sea Pressure', GotIt)
         IF (.Not.GotIt) ComputeS = .FALSE.
         IF (ComputeS) THEN
            n = BCElement % Type % NumberOfNodes
            DO i = 1, n
               j = BCElement % NodeIndexes (i)
               IF (NodeOnBoundary(j)==0) THEN
                 NumberOfNodesOnBoundary = NumberOfNodesOnBoundary + 1
                 NodeOnBoundary(j) = NumberOfNodesOnBoundary
               END IF
            END DO
         END IF
      END DO
      Model % CurrentElement => CurElement 

      ALLOCATE( S( NumberOfNodesOnBoundary ), Ns( NumberOfNodesOnBoundary ))
      IF (NormalFlux) ALLOCATE( a_perp ( NumberOfNodesOnBoundary ))
      
   ELSE
      IF (t > told) THEN
         NewTime = .TRUE.
         told = t
      END IF
   ENDIF  ! FirstTime
   
   IF (NewTime) THEN
      NewTime = .FALSE.

!-----------------------------------------------------------------
! get some information upon active boundary element and its parent
!-----------------------------------------------------------------
      BoundaryElement => Model % CurrentElement
      IF ( .NOT. ASSOCIATED(BoundaryElement) ) THEN
         CALL FATAL('Sea Pressure','No boundary element found')
      END IF
      other_body_id = BoundaryElement % BoundaryInfo % outbody
      IF (other_body_id < 1) THEN ! only one body in calculation
         ParentElement => BoundaryElement % BoundaryInfo % Right
         IF ( .NOT. ASSOCIATED(ParentElement) ) ParentElement => BoundaryElement % BoundaryInfo % Left
      ELSE ! we are dealing with a body-body boundary and asume that the normal is pointing outwards
         ParentElement => BoundaryElement % BoundaryInfo % Right
         IF (ParentElement % BodyId == other_body_id) ParentElement => BoundaryElement % BoundaryInfo % Left
      END IF
      
      body_id = ParentElement % BodyId
      material_id = ListGetInteger(Model % Bodies(body_id) % Values, 'Material', GotIt, UnFoundFatal=UnFoundFatal)
      ParentMaterial => Model % Materials(material_id) % Values
      IF ((.NOT. ASSOCIATED(ParentMaterial))) THEN
         WRITE(Message,'(A,I10,A,I10)')&
              'No material values found for body no ', body_id,&
              ' under material id ', material_id
         CALL FATAL('Sea Pressure',Message)
      END IF
      NParent = ParentElement % Type % NumberOfNodes
      
!-------------------------
! Get material parameters
! Sea level for that time
!-------------------------
     
      Zsl = GetCReal( ParentMaterial, 'Sea level', GotIt )
      IF (.NOT.GotIt) THEN
         WRITE(Message,'(A)') 'Variable Sea level not found. &
              &Setting to 0.0'
         CALL INFO('SeaPressure', Message, level=2)
         Zsl = 0.0_dp
      ELSE
         WRITE(Message,'(A,F10.4)') 'Sea level = ', Zsl 
         CALL INFO('SeaPressure', Message , level = 2)
      END IF
                     
!---------------------------     
! Compute S for new t
!---------------------------     
      CurElement => Model % CurrentElement
      S = 0.0
      DO p = 1, Model % NumberOfBoundaryElements
         BCElement => GetBoundaryElement(p) 
         BC => GetBC( BCElement ) 
         IF( GetElementFamily(BCElement) == 1 ) CYCLE
         ComputeS = GetLogical( BC, 'Compute Sea Pressure', GotIt)
         IF (.Not.GotIt) ComputeS = .FALSE.
         IF (ComputeS) THEN
            n = BCElement % Type % NumberOfNodes
            DO i = 1, n
               j = BCElement % NodeIndexes (i)
               IF (DIM==2) THEN
                 S(NodeOnBoundary(j)) = Model % Nodes % y (j)
               ELSE
                 S(NodeOnBoundary(j)) = Model % Nodes % z (j)
               END IF
            END DO
         END IF
      END DO
     Model % CurrentElement => CurElement 
!----------------------------------------------------------------------------
! Compute Ns and Flux (melting) for new t
! Restricted to the surface below the shelf 
! The melting is not taken into account where dS/dx is infinite (calving front)
!----------------------------------------------------------------------------
     Ns = 0.0_dp
     IF (NormalFlux) THEN
        a_perp = 0.0_dp
        ALLOCATE( Normal(3, NumberOfNodesOnBoundary))
        Normal = 0.0_dp
        CurElement => Model % CurrentElement
        DO p = 1, Model % NumberOfBoundaryElements
           BCElement => GetBoundaryElement(p) 
           BC => GetBC( BCElement ) 
           IF( GetElementFamily(BCElement) == 1 ) CYCLE
           ComputeS = GetLogical( BC, 'Compute Sea Spring', GotIt)
           IF (.Not.GotIt) ComputeS = .FALSE.
         
     ! we only compute the term Ns.b if the Spring is on         
           SourceFunc = 0.0_dp
           IF (ComputeS) THEN
               n = BCElement % Type % NumberOfNodes
               SourceFunc(1:n) = GetReal( Model % BodyForces(bf_id_FS) % Values, &
                     TRIM(BottomSurfaceName)//' Accumulation', GotIt)

             CALL GetElementNodes( Nodes , BCElement )
             DO i = 1,n
                 j = BCElement % NodeIndexes( i )
                 Bu = BCElement % Type % NodeU(i)
                 IF ( BCElement % Type % Dimension > 1 ) THEN
                    Bv = BCElement % Type % NodeV(i)
                 ELSE
                    Bv = 0.0D0
                 END IF
                 Normal(:,NodeOnBoundary(j))  = Normal(:,NodeOnBoundary(j)) + &
                                NormalVector(BCElement, Nodes, Bu, Bv, .TRUE.)
                 a_perp( NodeOnBoundary(j) ) = SourceFunc (i)
             END DO
           END IF
        END DO

        DO i=1, NumberOfNodesOnBoundary
          IF (ABS(Normal(DIM,i)) > 1.0e-20_dp) THEN
            Ns(i) = 1.0_dp + (Normal(1,i)/Normal(DIM,i))**2.0_dp
            IF (DIM>2) THEN
                 Ns(i) = Ns(i) + (Normal(2,i)/Normal(3,i))**2.0_dp
            END IF
            Ns(i) = SQRT(Ns(i))
          ELSE
            Ns(i) = -999.0
          END IF
        END DO
        DEALLOCATE (Normal)
        Model % CurrentElement => CurElement 
     END IF
   ENDIF  ! new dt
    
   j = NodeOnBoundary( nodenumber )
   IF ( NormalFlux .AND. (Ns(j) > 0.0_dp)) THEN 
      pw = -gravity * rhow * (Zsl - S(j) - a_perp(j) * dt * Ns(j) )  
   ELSE
      pw = -gravity * rhow * (Zsl - S(j))  
   END IF
   IF (pw > 0.0) pw = 0.0
END FUNCTION SeaPressure



FUNCTION SeaSpring ( Model, nodenumber, y) RESULT(C)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   TYPE(Solver_t):: Solver
   TYPE(Nodes_t), SAVE :: Nodes
   TYPE(variable_t), POINTER :: Timevar
   TYPE(Element_t), POINTER ::  BoundaryElement, BCElement, CurElement, ParentElement
   TYPE(ValueList_t), POINTER :: BC, material, ParentMaterial, BodyForce
   INTEGER :: NBoundary, NParent, BoundaryElementNode, ParentElementNode, body_id, other_body_id, material_id
   INTEGER :: nodenumber, NumberOfNodesOnBoundary 
   INTEGER, ALLOCATABLE :: NodeOnBoundary(:)
   INTEGER :: Nn, i, j, p, n, Nmax, bf_id, DIM 
   REAL(KIND=dp) :: y, C, t, told, dt, Bu, Bv
   REAL(KIND=dp) :: rhow, gravity
   REAL(KIND=dp), ALLOCATABLE :: auxReal(:), Ns(:), normal(:,:)
   LOGICAL :: FirstTime = .TRUE., NewTime, GotIt, ComputeS   
       
   SAVE told, FirstTime, NewTime, Nn, dt, Ns, Bodyforce, DIM
   SAVE rhow, gravity, auxReal
   SAVE NumberOfNodesOnBoundary, NodeOnBoundary, normal 

   Timevar => VariableGet( Model % Variables,'Time')
   t = TimeVar % Values(1)
   dt = Model % Solver % dt 

   IF (FirstTime) THEN
      FirstTime = .FALSE.
      NewTime = .TRUE.
      told = t
      ALLOCATE( NodeOnBoundary( Model % Mesh % NumberOfNodes ))
      n = Model % MaxElementNodes 
      ALLOCATE( auxReal(n) )
      DIM = CoordinateSystemDimension()

      rhow = GetConstReal( Model % Constants, 'Water Density', GotIt )
      IF (.NOT.GotIt) THEN
         WRITE(Message,'(A)') 'Variable Water Density not found. &
              &Setting to 1.03225e-18'
         CALL INFO('SeaSpring', Message, level=20)
         rhow = 1.03225e-18_dp
      ELSE
         WRITE(Message,'(A,F10.4)') 'Water Density = ', rhow
         CALL INFO('SeaSpring', Message , level = 20)
      END IF
      
      !-----------------------------------------------------------------
      ! get some information upon active boundary element and its parent
      !-----------------------------------------------------------------
      BoundaryElement => Model % CurrentElement
      IF ( .NOT. ASSOCIATED(BoundaryElement) ) THEN
         CALL FATAL('Sea Pressure','No boundary element found')
      END IF
      other_body_id = BoundaryElement % BoundaryInfo % outbody
      IF (other_body_id < 1) THEN ! only one body in calculation
         ParentElement => BoundaryElement % BoundaryInfo % Right
         IF ( .NOT. ASSOCIATED(ParentElement) ) ParentElement => BoundaryElement % BoundaryInfo % Left
      ELSE ! we are dealing with a body-body boundary and asume that the normal is pointing outwards
         ParentElement => BoundaryElement % BoundaryInfo % Right
         IF (ParentElement % BodyId == other_body_id) ParentElement => BoundaryElement % BoundaryInfo % Left
      END IF

      BodyForce => GetBodyForce(ParentElement)  
      bf_id = ListGetInteger( Model % Bodies(1) % Values,&
           'Body Force',  minv=1, maxv=Model % NumberOFMaterials)
      IF (DIM==2) THEN !NOT the same keyword if the shape factor is used
           gravity = -GetConstReal( Model % BodyForces(bf_id) % Values, 'Lateral Friction Gravity 2',Gotit)
           IF (.NOT. Gotit) THEN
              gravity = -GetConstReal( Model % BodyForces(bf_id) % Values, 'Flow BodyForce 2',Gotit)
           END IF
      ELSE
         gravity = -GetConstReal( Model % BodyForces(bf_id) % Values, 'Flow BodyForce 3',Gotit)
      END IF

! Number of node concerned 

      CurElement => Model % CurrentElement
      NodeOnBoundary = 0
      NumberOfNodesOnBoundary = 0
      DO p = 1, Model % NumberOfBoundaryElements
         BCElement => GetBoundaryElement(p) 
         BC => GetBC( BCElement ) 
         IF( GetElementFamily(BCElement) == 1 ) CYCLE
         ComputeS = GetLogical( BC, 'Compute Sea Spring', GotIt)
         IF (.Not.GotIt) ComputeS = .FALSE.
         IF (ComputeS) THEN
            n = BCElement % Type % NumberOfNodes
            DO i = 1, n
               j = BCElement % NodeIndexes (i)
               IF (NodeOnBoundary(j)==0) THEN
                 NumberOfNodesOnBoundary = NumberOfNodesOnBoundary + 1
                 NodeOnBoundary(j) = NumberOfNodesOnBoundary
               END IF
            END DO
         END IF
      END DO
      Model % CurrentElement => CurElement 

      ALLOCATE( Ns( NumberOfNodesOnBoundary ))
      
   ELSE
      IF (t > told) THEN
         NewTime = .TRUE.
         told = t
      END IF
   ENDIF
   
   IF (NewTime) THEN
      NewTime = .FALSE.
      
      ! Compute Ns for new t
      ALLOCATE( Normal(3, NumberOfNodesOnBoundary))
      Ns = 0.0_dp
      Normal = 0.0_dp
      CurElement => Model % CurrentElement
      DO p = 1, Model % NumberOfBoundaryElements
         BCElement => GetBoundaryElement(p) 
         BC => GetBC( BCElement ) 
         IF( GetElementFamily(BCElement) == 1 ) CYCLE
         ComputeS = GetLogical( BC, 'Compute Sea Spring', GotIt)
         IF (.Not.GotIt) ComputeS = .FALSE.
         IF (ComputeS) THEN
           CALL GetElementNodes( Nodes , BCElement )
           n = BCElement % Type % NumberOfNodes
           DO i = 1,n
             j = BCElement % NodeIndexes( i )
             Bu = BCElement % Type % NodeU(i)
             IF ( BCElement % Type % Dimension > 1 ) THEN
                 Bv = BCElement % Type % NodeV(i)
             ELSE
                Bv = 0.0D0
             END IF
             Normal(:,NodeOnBoundary(j))  = Normal(:,NodeOnBoundary(j)) + &
                                NormalVector(BCElement, Nodes, Bu, Bv, .TRUE.)
           END DO
         END IF
      END DO

      DO i=1, NumberOfNodesOnBoundary
      IF (ABS(Normal(DIM,i)) > 1.0e-20_dp) THEN
         Ns(i) = 1.0_dp + (Normal(1,i)/Normal(DIM,i))**2.0_dp
         IF (DIM>2) THEN
               Ns(i) = Ns(i) + (Normal(2,i)/Normal(3,i))**2.0_dp
         END IF
         Ns(i) = SQRT(Ns(i))
      ELSE
         Ns(i) = -999.0
      END IF
      END DO
      DEALLOCATE (Normal)
      Model % CurrentElement => CurElement 

   ENDIF  ! new dt
    
   j = NodeOnBoundary( nodenumber )
   IF (Ns(j) > 0.0_dp) THEN 
      C = gravity * rhow * dt * Ns(j)                                   
   ELSE
      C = 1.0e20_dp
   END IF
END FUNCTION SeaSpring



FUNCTION BedPressure ( Model, nodenumber, y) RESULT(pw)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   TYPE(Solver_t):: Solver
   TYPE(Nodes_t), SAVE :: Nodes
   TYPE(variable_t), POINTER :: Timevar
   TYPE(Element_t), POINTER ::  BoundaryElement, BCElement, CurElement, ParentElement
   TYPE(ValueList_t), POINTER :: BC, material, ParentMaterial, BodyForce
   INTEGER :: NBoundary, NParent, BoundaryElementNode, ParentElementNode, body_id, other_body_id, material_id
   INTEGER :: nodenumber, NumberOfNodesOnBoundary 
   INTEGER, ALLOCATABLE :: NodeOnBoundary(:)
   INTEGER :: Nn, i, j, p, n, Nmax, bf_id, DIM, bf_id_FS 
   REAL(KIND=dp) :: y, pw, t, told, dt, Bu, Bv
   REAL(KIND=dp) :: Zsl, rhow, gravity
   REAL(KIND=dp), ALLOCATABLE :: S(:), auxReal(:), Ns(:),  a_perp(:), SourceFunc(:), normal(:,:)
   LOGICAL :: FirstTime = .TRUE., NewTime, GotIt, ComputeS,  NormalFlux = .TRUE. 
   CHARACTER(LEN=MAX_NAME_LEN)  :: BottomSurfaceName
       
!------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Variable_t), POINTER :: bedrockVar
  LOGICAL :: AllocationsDone = .FALSE.

  INTEGER :: bedrockSource
  INTEGER, POINTER ::  bedrockPerm(:)
  CHARACTER(LEN=MAX_NAME_LEN) ::  bedrockName
  INTEGER,PARAMETER :: MATERIAL_DEFAULT = 1, MATERIAL_NAMED = 2, VARIABLE = 3

  SAVE AllocationsDone

!------------------------------------------------------------


   SAVE told, FirstTime, NewTime, Nn, dt, Ns, Bodyforce, DIM
   SAVE S, rhow, gravity, Zsl, auxReal, NormalFlux, a_perp, SourceFunc
   SAVE NumberOfNodesOnBoundary, NodeOnBoundary, normal
   SAVE BottomSurfaceName, bf_id_FS 
   
  !--------------------------------------------------------------
  ! Allocate some permanent storage:
  !--------------------------------------------------------------
   Timevar => VariableGet( Model % Variables,'Time')
   t = TimeVar % Values(1)
   dt = Model % Solver % dt 

   IF (FirstTime) THEN
      FirstTime = .FALSE.
      NewTime = .TRUE.
      told = t
      ALLOCATE( NodeOnBoundary( Model % Mesh % NumberOfNodes ))
      n = Model % MaxElementNodes 
      ALLOCATE( auxReal(n), SourceFunc(n) )
      DIM = CoordinateSystemDimension()

      rhow = GetConstReal( Model % Constants, 'Water Density', GotIt )
      IF (.NOT.GotIt) THEN
         WRITE(Message,'(A)') 'Variable Water Density not found. &
              &Setting to 1.03225e-18'
         CALL INFO('BedPressure', Message, level=20)
         rhow = 1.03225e-18_dp
      ELSE
         WRITE(Message,'(A,F10.4)') 'Water Density = ', rhow
         CALL INFO('BedPressure', Message , level = 20)
      END IF
      
!-----------------------------------------------------------------
! Is there basal melt to account for          
!-----------------------------------------------------------------
      NormalFlux = GetLogical( Model % Constants, 'Buoyancy Use Basal Melt', GotIt )    
      IF (.NOT.GotIt) THEN
         WRITE(Message,'(A)') 'Variable Buoyancy Use Basal Melt not found. &
              &Setting to FALSE'
         CALL INFO('BedPressure', Message, level=20)
         NormalFlux = .FALSE.
      ELSE
         WRITE(Message,*) 'Buoyancy Use Basal Melt = ', NormalFlux
         CALL INFO('BedPressure', Message , level = 20)
      END IF

!-----------------------------------------------------------------
! Find the bf_id of the FS solvers if NormalFlux is TRUE
! and read the name of the Bottom Free surface solver
!-----------------------------------------------------------------
      IF (NormalFlux) THEN
        BottomSurfaceName = GetString( Model % Constants, 'Bottom Surface Name', GotIt )    
        IF (.NOT.GotIt) THEN
           WRITE(Message,'(A)') 'Variable Bottom SUrface Name not found. &
                &Setting to DyBottom'
           CALL INFO('BedPressure', Message, level=20)
           BottomSurfaceName = 'DyBottom'
        ELSE
           WRITE(Message,*)'Bottom Surface Name = ',TRIM(BottomSurfaceName)
           CALL INFO('BedPressure', Message , level = 20)
        END IF

        bf_id_FS = -1
        DO bf_id=1,Model % NumberOFBodyForces
           IF( ListCheckPresent( Model % BodyForces(bf_id) % Values, &
                TRIM(BottomSurfaceName)//' Accumulation') ) bf_id_FS = bf_id 
        END DO

        IF (bf_id_FS<0) THEN 
           CALL FATAL('BedPressure','No Basal Melt found in any body Force')
        END IF 
      END IF

!-----------------------------------------------------------------
! Read the gravity 
!-----------------------------------------------------------------
      bf_id = ListGetInteger( Model % Bodies(1) % Values,&
           'Body Force',  minv=1, maxv=Model % NumberOFMaterials)

      IF (DIM==2) THEN    !NOT the same keyword if the shape factor is used
        gravity = -GetConstReal( Model % BodyForces(bf_id) % Values, 'Lateral Friction Gravity 2',Gotit)
        IF (.NOT. Gotit) THEN
          gravity = -GetConstReal( Model % BodyForces(bf_id) % Values, 'Flow BodyForce 2',Gotit)
        END IF
      ELSE
        gravity = -GetConstReal( Model % BodyForces(bf_id) % Values, 'Flow BodyForce 3',Gotit)
      END IF
 
!-----------------------------------------------------------------
! Number of nodes concerned (needed to compute Ns from the normal)
!-----------------------------------------------------------------

      CurElement => Model % CurrentElement
      NodeOnBoundary = 0
      NumberOfNodesOnBoundary = 0
      DO p = 1, Model % NumberOfBoundaryElements
         BCElement => GetBoundaryElement(p) 
         BC => GetBC( BCElement ) 
         IF( GetElementFamily(BCElement) == 1 ) CYCLE
         ComputeS = GetLogical( BC, 'Compute Bed Pressure', GotIt)
         IF (.Not.GotIt) ComputeS = .FALSE.
         IF (ComputeS) THEN
            n = BCElement % Type % NumberOfNodes
            DO i = 1, n
               j = BCElement % NodeIndexes (i)
               IF (NodeOnBoundary(j)==0) THEN
                 NumberOfNodesOnBoundary = NumberOfNodesOnBoundary + 1
                 NodeOnBoundary(j) = NumberOfNodesOnBoundary
               END IF
            END DO
         END IF
      END DO
      Model % CurrentElement => CurElement 

      ALLOCATE( S( NumberOfNodesOnBoundary ), Ns( NumberOfNodesOnBoundary ))
      IF (NormalFlux) ALLOCATE( a_perp ( NumberOfNodesOnBoundary ))
      
   ELSE
      IF (t > told) THEN
         NewTime = .TRUE.
         told = t
      END IF
   ENDIF  ! FirstTime
   
   IF (NewTime) THEN
      NewTime = .FALSE.

!-----------------------------------------------------------------
! get some information upon active boundary element and its parent
!-----------------------------------------------------------------
      BoundaryElement => Model % CurrentElement
      IF ( .NOT. ASSOCIATED(BoundaryElement) ) THEN
         CALL FATAL('BedPressure','No boundary element found')
      END IF
      other_body_id = BoundaryElement % BoundaryInfo % outbody
      IF (other_body_id < 1) THEN ! only one body in calculation
         ParentElement => BoundaryElement % BoundaryInfo % Right
         IF ( .NOT. ASSOCIATED(ParentElement) ) ParentElement => BoundaryElement % BoundaryInfo % Left
      ELSE ! we are dealing with a body-body boundary and asume that the normal is pointing outwards
         ParentElement => BoundaryElement % BoundaryInfo % Right
         IF (ParentElement % BodyId == other_body_id) ParentElement => BoundaryElement % BoundaryInfo % Left
      END IF
      
      body_id = ParentElement % BodyId
      material_id = ListGetInteger(Model % Bodies(body_id) % Values, 'Material', GotIt)
      ParentMaterial => Model % Materials(material_id) % Values
      IF ((.NOT. ASSOCIATED(ParentMaterial)) .OR. (.NOT. GotIt)) THEN
         WRITE(Message,'(A,I10,A,I10)')&
              'No material values found for body no ', body_id,&
              ' under material id ', material_id
         CALL FATAL('BedPressure',Message)
      END IF
      NParent = ParentElement % Type % NumberOfNodes
      
!-------------------------
! Get material parameters
! Sea level for that time
!-------------------------
     
      Zsl = GetConstReal( ParentMaterial, 'Sea level', GotIt )

      IF (.NOT.GotIt) THEN
         WRITE(Message,'(A)') 'Variable Sea level not found. &
              &Setting to 0.0'
         CALL INFO('BedPressure', Message, level=2)
         Zsl = 0.0_dp
      ELSE
         WRITE(Message,'(A,F10.4)') 'Sea level = ', Zsl 
         CALL INFO('BedPressure', Message , level = 2)
      END IF
           

!-------------------------
!  Get bedrock topography
!-------------------------
      Solver = Model % Solver
      IF ( (.NOT. AllocationsDone)) THEN
         CALL INFO( 'BedPressure', 'Memory allocation done.',Level=1 )
         AllocationsDone = .TRUE.
         SolverParams => GetSolverParams()
        IF (.NOT.GotIt) THEN
           CALL FATAL('BedPressure', 'No tolerance given for the Grounded Mask.')
        END IF

        bedrockName = GetString(SolverParams, 'Bedrock Variable', GotIt)
        IF (GotIt) THEN
           bedrockSource = VARIABLE
           CALL info('BedPressure', 'Bedrock Variable name found', level=8)
        ELSE
           bedrockName = GetString(SolverParams, 'Bedrock Material', GotIt)
           IF (GotIt) THEN
              bedrockSource = MATERIAL_NAMED
              CALL info('BedPressure', 'Bedrock Material name found', level=8)
           ELSE
              bedrockSource = MATERIAL_DEFAULT     
              CALL info('BedPressure', 'No Bedrock Variable or Material; searching for material \"Min Zs Bottom\".', level=8)
           END IF
        END IF
          CurElement => Model % CurrentElement

         SELECT CASE(bedrockSource)
         CASE (VARIABLE)
            bedrockVar => VariableGet(Model % Mesh % Variables, bedrockName )
            IF (.NOT. ASSOCIATED(bedrockVar)) CALL FATAL('BedPressure',"Could not find bedrock variable")
            bedrockPerm => bedrockVar % Perm
            DO j = 1, NumberOfNodesOnBoundary
              S(j) =  bedrockVar % values(bedrockPerm(j))
            END DO
            NULLIFY(bedrockPerm)
            NULLIFY(bedrockVar)
         CASE (MATERIAL_NAMED)
!             Material => GetMaterial( BCElement )
!             S = ListGetReal( Material,bedrockName, NumberOfNodesOnBoundary , & 
!                  BCElement % NodeIndexes, GotIt )
!             IF (.NOT. GotIt) CALL FATAL('BedPressure',"Could not find bedrock material")
         CASE (MATERIAL_DEFAULT)
!             Material => GetMaterial( BCElement )
!             S = ListGetReal( Material,'Min Zs Bottom',NumberOfNodesOnBoundary , & 
!                  BCElement % NodeIndexes, GotIt )
!             IF (.NOT. GotIt) CALL FATAL('BedPressure',"Could not find bedrock material")
         END SELECT

      END IF  
   ENDIF  ! new dt
    
   j = NodeOnBoundary( nodenumber )
   pw = -gravity * rhow * (Zsl - S(j))  
   IF (pw > 0.0) pw = 0.0
END FUNCTION BedPressure
