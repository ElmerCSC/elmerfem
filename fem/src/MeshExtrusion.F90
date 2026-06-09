!/******************************************************************************
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! *
! *  This library is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU Lesser General Public
! *  License as published by the Free Software Foundation; either
! *  version 2.1 of the License, or (at your option) any later version.
! *
! *  This library is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! *  Lesser General Public License for more details.
! *
! *  You should have received a copy of the GNU Lesser General Public
! *  License along with this library (in file ../LGPL-2.1); if not, write
! *  to the Free Software Foundation, Inc., 51 Franklin Street,
! *  Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Authors: Juha Ruokolainen, Peter Raback
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland
! *
! *****************************************************************************/

!> \ingroup ElmerLib
!> \{

!------------------------------------------------------------------------------
!>  Mesh extrusion: helpers and 2D-to-3D extrusion functions.
!>  Extracted from MeshUtils.
!------------------------------------------------------------------------------

MODULE MeshExtrusion

    USE MeshBasics
    USE MeshLoad, ONLY : PrepareMesh
    IMPLICIT NONE

CONTAINS



  
  FUNCTION MapExtrudedMaterial(Vlist,mat0,ilayer,EndLayer) RESULT ( mat )
    TYPE(Valuelist_t), POINTER :: Vlist
    INTEGER :: mat0, mat
    INTEGER, OPTIONAL :: ilayer
    LOGICAL, OPTIONAL :: EndLayer

    TYPE(ValueList_t), POINTER, SAVE :: PrevList
    LOGICAL, SAVE :: EndMat, ExtMat
    INTEGER, POINTER, SAVE :: ExtrudedElements(:)
    INTEGER, ALLOCATABLE, SAVE :: InvExtrudedElements(:)
    INTEGER, SAVE :: nDiv, nElems
    INTEGER :: i,j
    LOGICAL :: SetMat
    
    IF(.NOT. ASSOCIATED(PrevList,Vlist)) THEN
      IF(ALLOCATED(InvExtrudedElements))  DEALLOCATE(InvExtrudedElements)
           
      PrevList => Vlist
      EndMat = ListCheckPresent( Vlist,'Extruded Mesh End Map')
      IF(EndMat) THEN
        CALL Info('MapExtrudedMaterial','Extruded Mesh will be mapped at the ends!')
      END IF        

      ExtrudedElements => ListGetIntegerArray(Vlist,'Extruded Elements',ExtMat)                
      IF(ExtMat) THEN
        nDiv = SIZE(ExtrudedElements)
        nElems = SUM(ExtrudedElements)
        ALLOCATE(InvExtrudedElements(nElems))
        InvExtrudedElements = 0
        j = 0
        DO i=1,nDiv
          IF( ExtrudedElements(i) == 0) CYCLE
          InvExtrudedElements(j+1:j+ExtrudedElements(i)) = i
          j = j+ExtrudedElements(i)
        END DO
      ELSE
        nElems = ListGetInteger(Vlist,'Extruded Mesh Layers',UnfoundFatal=.TRUE.)
      END IF
    END IF

    mat = mat0
    IF( EndMat ) THEN
      IF(ilayer < 1 .OR. ilayer > nElems ) THEN
        CALL Fatal('MapExtrudedMaterial','Invalid body id: '//I2S(ilayer))
      END IF

      SetMat = .FALSE.
      IF( ExtMat ) THEN
        j = InvExtrudedElements(ilayer)
        SetMat = (j==1 .OR. j==nDiv)
      ELSE
        SetMat = (ilayer==1 .OR. ilayer==nElems)
      END IF             
      IF(SetMat) mat = NINT(ListGetFun( Vlist,'Extruded Mesh End Map',1.0_dp * mat0 ) )
    END IF

  END FUNCTION MapExtrudedMaterial

  
  
  SUBROUTINE CheckPointElementParents(Mesh)
    TYPE(Mesh_t) :: Mesh
    LOGICAL :: Found
    INTEGER :: Misses(3)
    TYPE(Element_t), POINTER :: Element, Parent
    INTEGER :: i,j,t,t1,t2
    INTEGER, ALLOCATABLE :: OneOwner(:)
    
    t1 = Mesh % NumberOfBulkElements
    t2 = Mesh % NumberOfBoundaryElements

    Misses = 0
    DO t=t1+1,t1+t2
      Element => Mesh % Elements(t)
      IF(Element % TYPE % NumberOfNodes > 1) CYCLE
      IF(.NOT. ASSOCIATED(Element % BoundaryInfo)) THEN
        Misses(1) = Misses(1) + 1
        CYCLE
      END IF
      Parent => Element % BoundaryInfo % Left
      IF(ASSOCIATED(Parent)) THEN
        i = Element % NodeIndexes(1)
        IF(ALL(Parent % NodeIndexes /= i)) Misses(2) = Misses(2) + 1
      END IF
      Parent => Element % BoundaryInfo % Right
      IF(ASSOCIATED(Parent)) THEN
        i = Element % NodeIndexes(1)
        IF(ALL(Parent % NodeIndexes /= i)) Misses(3) = Misses(3) + 1
      END IF      
    END DO

    i = SUM(Misses)
    IF(i == 0) RETURN

    IF( i > 0 ) THEN
      CALL Info('CheckPointElementParents',&
          'We have point '//I2S(i)//' elements with faulty parents!')
    END IF
    
    ALLOCATE(OneOwner(Mesh % NumberOfNodes))
    OneOwner = 0
    DO t=1,t1      
      Element => Mesh % Elements(t)
      DO i=1,Element % TYPE % NumberOfNodes
        j = Element % NodeIndexes(i)
        IF(OneOwner(j)==0) OneOwner(j) = t
      END DO
    END DO

    DO t=t1+1,t1+t2
      Element => Mesh % Elements(t)
      IF(Element % TYPE % NumberOfNodes > 1) CYCLE
      Parent => Element % BoundaryInfo % Left
      IF(ASSOCIATED(Parent)) THEN
        i = Element % NodeIndexes(1)
        IF(ALL(Parent % NodeIndexes /= i)) THEN
          Element % BoundaryInfo % Left => Mesh % Elements(OneOwner(i))
        END IF
      END IF
      Element % ElementIndex = t
    END DO
    
  END SUBROUTINE CheckPointElementParents

  
  ! Collect here the routines that defines the division in the exruded direction.
  !-----------------------------------------------------------------------------  
  SUBROUTINE ExtrudedDivision(Vlist, nlevels, Wtable)
    TYPE(ValueList_t), POINTER :: Vlist
    INTEGER :: nlevels
    REAL(KIND=dp), ALLOCATABLE :: Wtable(:)

    LOGICAL :: Found, GotLimits    
    REAL(KIND=dp) :: q,zmin,zmax,z
    INTEGER :: i,j,k,nDiv
    REAL(KIND=dp), POINTER :: ExtrudedLimits(:,:), ExtrudedSizes(:,:), ExtrudedRatios(:,:)
    INTEGER, POINTER :: ExtrudedElements(:)
    REAL(KIND=dp), ALLOCATABLE :: Wtmp(:)
    
    nlevels = ListGetInteger(Vlist,'Extruded Mesh Layers',Found)
    IF( .NOT. Found ) THEN
      nlevels = ListGetInteger(Vlist,'Extruded Mesh Levels',Found)-1 
      IF(Found) THEN
        CALL ListAddNewInteger(Vlist,'Extruded Mesh Layers',nlevels)     
      END IF
    END IF
    IF(Found ) THEN
      q = ListGetCReal(Vlist,'Extruded Mesh Ratio',Found )
      IF(.NOT. Found) q = 1.0_dp
      ALLOCATE(Wtable(0:nlevels))
      CALL UnitSegmentDivision(Wtable,nlevels,Vlist)
      zmin = ListGetCReal(Vlist,'Extruded Min Coordinate',Found )
      zmax = ListGetCReal(Vlist,'Extruded Max Coordinate',Found )
      IF(.NOT. Found) zmax = 1.0_dp

      Wtable = zmin + (zmax-zmin) * Wtable      
    ELSE
      ExtrudedElements => ListGetIntegerArray(Vlist,'Extruded Elements',Found)                
      IF(.NOT. Found ) CALL Fatal('ExtrudedDivision','We should not even be here!')      
      nDiv = SIZE(ExtrudedElements) 
      
      ExtrudedLimits => ListGetConstRealArray(Vlist,'Extruded Limits',GotLimits) 
      IF(GotLimits) THEN
        IF(SIZE(ExtrudedLimits,1) /= nDiv+1 .OR. SIZE(ExtrudedLimits,2) /= 1) THEN
          CALL Fatal('ExtrudedDivision','Incompatible size for "Extruded Limits"')
        END IF
      ELSE
        ExtrudedSizes => ListGetConstRealArray(Vlist,'Extruded Sizes',Found ) 
        IF(.NOT. Found) THEN
          CALL Fatal('ExtrudedDivision','Give either "Extruded Limits" or "Extruded Sizes"!')
        END IF
        IF(SIZE(ExtrudedSizes,1) /= nDiv .OR. SIZE(ExtrudedSizes,2) /= 1) THEN
          CALL Fatal('ExtrudedDivision','Incompatible size for "Extruded Sizes"')
        END IF
      END IF

      ExtrudedRatios => ListGetConstRealArray(Vlist,'Extruded Ratios',Found)                
      IF(Found) THEN
        IF(SIZE(ExtrudedRatios,1) /= nDiv .OR. SIZE(ExtrudedRatios,2) /= 1) THEN
          CALL Fatal('ExtrudedDivision','Incompatible size for "Extruded Elements"')
        END IF
      END IF

      i = MAXVAL(ExtrudedElements)
      nlevels = SUM(ExtrudedElements)
      ALLOCATE(Wtable(0:nlevels),Wtmp(0:i))
      j = 0
      q = 1.0_dp
      DO i=1,nDiv
        IF(ASSOCIATED(ExtrudedRatios)) q = ExtrudedRatios(i,1)        

        k = ExtrudedElements(i)
        CALL GeometricUnitDivision(Wtmp,k,q)

        IF(GotLimits) THEN          
          Wtable(j:j+k) = ExtrudedLimits(i,1) + &
              Wtmp(0:k)*(ExtrudedLimits(i+1,1)-ExtrudedLimits(i,1))
        ELSE
          Wtable(j:j+k) = z + ExtrudedSizes(i,1)*Wtmp(0:k) 
          z = z + ExtrudedSizes(i,1)
        END IF
        j = j + k
      END DO
      DO i=0,nlevels
        WRITE( Message, '(A,I0,A,ES12.4)') 'w(',i,') : ',wTable(i)
        CALL Info('ExtrudedDivision', Message )
      END DO

      !CALL ListAddNewConstReal(Vlist,'Extruded Min Coordinate',Wtable(0) )
      !CALL ListAddNewConstReal(Vlist,'Extruded Max Coordinate',Wtable(nlevels) )
    END IF
    
    IF(nlevels < 2) THEN
      CALL Fatal('ExtrudedDivision','There must be at least two "Extruded Mesh Layers"!')
    END IF    
      
  END SUBROUTINE ExtrudedDivision


  ! Enable skew for extruded or initially 3D mesh, mainly intended for electrical
  ! machines. This is a library routine since we may want to perform skew right
  ! after the extrusion, if the mesh is further to be split into other elements. 
  !-----------------------------------------------------------------------------  
  SUBROUTINE SetMeshSkew(Mesh, Vlist )
    TYPE(ValueList_t), POINTER :: Vlist
    TYPE(Mesh_t) :: Mesh
    REAL(KIND=dp) :: RotorRad, AngleCoeff, RotorSkew, StatorSkew
    REAL(KIND=dp) :: zmin, zmax, Coord(3), zloc, alpha, minskew, maxskew
    LOGICAL :: Found, GotSkewFun, GotSkew, IsRotor
    LOGICAL, ALLOCATABLE :: NodeDone(:)
    INTEGER :: NoNodes, elem, n, i, j, NodeIndex(1)
    LOGICAL :: SkewDone = .FALSE.        
    TYPE(Element_t), POINTER :: Element
    INTEGER, POINTER :: RotorBodies(:)
    CHARACTER(*), PARAMETER :: Caller="SetMeshSkew"   

    SAVE SkewDone

    IF(SkewDone) THEN
      CALL Info(Caller,'Skew already done!',Level=10)
      RETURN
    END IF

    RotorBodies => ListGetIntegerArray( Vlist,'Rotor Bodies',Found )
    IF(.NOT. ASSOCIATED(RotorBodies) ) THEN
      RotorRad = ListGetCReal(Vlist,'Rotor Radius',Found )
      IF(.NOT. Found) THEN
        CALL Info(Caller,'Neither "Rotor Radius" or "Rotor Bodies" given!',Level=10)
        RETURN
      END IF
    END IF
    
    IF( ListGetLogical( Vlist,'Rotate in Radians',Found ) ) THEN
      CALL Info(Caller,'Using radians for skew!',Level=10)
      AngleCoeff = 1.0_dp
    ELSE
      CALL Info(Caller,'Using degrees for skew!',Level=10)
      AngleCoeff = PI / 180.0_dp
    END IF

    RotorSkew = AngleCoeff * ListGetCReal(Vlist,'Rotor Skew',GotSkew )
    GotSkewFun = ListCheckPresent( Vlist,'Rotor Skew Function')
    StatorSkew = AngleCoeff * ListGetCReal(Vlist,'Stator Skew',Found )
    GotSkew = GotSkew .OR. GotSkewFun .OR. Found
    IF(.NOT. GotSkew) THEN
      CALL Info(Caller,'No settings for skew given!',Level=10)
      RETURN
    END IF

    NoNodes = Mesh % NumberOfNodes

    zmin = ListGetCReal( Vlist,'Rotor Skew Min Coordinate',Found ) 
    IF(.NOT. Found) THEN
      zmin = ListGetCReal( Vlist,'Extruded Min Coordinate',Found ) 
    END IF
    IF(.NOT. Found) THEN
      zmin = MINVAL(Mesh % Nodes % z(1:NoNodes))
      zmin = ParallelReduction(zmin,1)
    END IF

    zmax = ListGetCReal( Vlist,'Rotor Skew Max Coordinate',Found ) 
    IF(.NOT. Found) THEN
      zmax = ListGetCReal( Vlist,'Extruded Max Coordinate',Found ) 
    END IF
    IF(.NOT. Found) THEN
      zmax = MAXVAL(Mesh % Nodes % z(1:NoNodes))
      zmax = ParallelReduction(zmax,2)
    END IF
    
    WRITE(Message,'(A,2ES12.3)') 'Coordinate range for extrusion:',zmin,zmax    
    CALL Info(Caller,Message)
    
    ALLOCATE(NodeDone(NoNodes))
    NodeDone = .FALSE.

    maxskew = -HUGE(maxskew)
    minskew = HUGE(minskew)
    
    DO elem = 1,Mesh % NumberOfBulkElements      
      Element => Mesh % Elements(elem)
      n = Element % TYPE % NumberOfNodes

      Coord(1) = SUM(Mesh % Nodes % x(Element % NodeIndexes)) / n
      Coord(2) = SUM(Mesh % Nodes % y(Element % NodeIndexes)) / n
      Coord(3) = SUM(Mesh % Nodes % z(Element % NodeIndexes)) / n

      IF(ASSOCIATED(RotorBodies)) THEN
        IsRotor = ANY( RotorBodies == Element % BodyId ) 
      ELSE
        IsRotor = (Coord(1)**2+Coord(2)**2 < RotorRad**2) 
      END IF

      DO i=1,n
        j = Element % NodeIndexes(i)
        NodeIndex(1) = j
        IF(.NOT. NodeDone(j)) THEN
          Coord(1) = Mesh % Nodes % x(j)
          Coord(2) = Mesh % Nodes % y(j)
          Coord(3) = Mesh % Nodes % z(j)

          ! Skew is not constant, perform it for each node 1st if requested. 
          zloc = (coord(3)-zmin)/(zmax-zmin)

          ! By construction this must be in [0,1]
          zloc = MAX(0.0_dp,MIN(1.0_dp,zloc))

          IF( IsRotor ) THEN
            IF(GotSkewFun) THEN
              alpha = AngleCoeff * ListGetFun( Vlist,'Rotor Skew Function',zloc)                
            ELSE
              alpha = (zloc-0.5_dp) * RotorSkew
            END IF
          ELSE
            alpha = (zloc-0.5_dp) * StatorSkew 
          END IF

          maxskew = MAX(alpha, maxskew)
          minskew = MIN(alpha, minskew)
          
          Mesh % Nodes % x(j) = Coord(1)*COS(alpha) - Coord(2)*SIN(alpha)
          Mesh % Nodes % y(j) = Coord(1)*SIN(alpha) + Coord(2)*COS(alpha)        
          NodeDone(j) = .TRUE.
        END IF
      END DO
    END DO

    SkewDone = .TRUE.

    IF(InfoActive(10)) THEN
      IF(GotSkewFun) THEN
        minskew = (180.0/PI) * ParallelReduction(minskew,1)
        maxskew = (180.0/PI) * ParallelReduction(maxskew,2)
        WRITE(Message,'(A,2ES12.3)') 'Rotor skew done with range (degrees): ',minskew,maxskew
      ELSE
        WRITE(Message,'(A,2ES12.3)') 'Rotor skew done with total angle: ',(180.0/PI)*RotorSkew
      END IF
      CALL Info(Caller,Message)
    END IF

  END SUBROUTINE SetMeshSkew
    
  
!------------------------------------------------------------------------------
!> Given a 2D mesh extrude it to be 3D. The 3rd coordinate will always
!> be at the interval [0,1]. Therefore the adaptation for different shapes
!> must be done with StructuredMeshMapper, or some similar utility. 
!> The top and bottom surface will be assigned Boundary Condition tags
!> with indexes one larger than the maximum used on by the 2D mesh. 
!> NOTE: This function handles NDOFs of the element structure in a way
!>       which is not consistent with "Element = n:N ...", with N>1 
!------------------------------------------------------------------------------
  FUNCTION MeshExtrude(Mesh_in, Vlist) RESULT(Mesh_out)
!------------------------------------------------------------------------------
    TYPE(Mesh_t), TARGET :: Mesh_in
    TYPE(Mesh_t), POINTER :: Mesh_out
    TYPE(ValueList_t), POINTER :: Vlist
!------------------------------------------------------------------------------
    CHARACTER(:), ALLOCATABLE :: ExtrudedMeshName
    INTEGER :: i,j,k,l,n,cnt,ind(8),max_baseline_bid,max_bid,l_n,max_body,&
        ExtrudedCoord,dg_n,totalnumberofelements
    TYPE(Element_t), POINTER :: Elem_in, Elem_out
    TYPE(ParallelInfo_t), POINTER :: PI_in, PI_out
    INTEGER :: nnodes,gnodes,gelements,ierr,bcignored,cnt101
    LOGICAL :: isParallel, Found, PreserveBaseline, Rotational, Rotate2Pi, CollectExtrudedBCs
    REAL(KIND=dp)::CurrCoord 
    REAL(KIND=dp), POINTER :: ActiveCoord(:)
    REAL(KIND=dp), ALLOCATABLE :: Wtable(:)
    INTEGER, POINTER :: BCLayers(:), TmpLayers(:)
    INTEGER :: NoBCLayers, bcoffset, baseline0, bclevel, BaseLineLayer, bcind, &
        m, max_bid0, in_levels, nlev
    INTEGER :: BcCounter(100)    
    LOGICAL :: GotBCLayers, DoCount
    CHARACTER(*), PARAMETER :: Caller="MeshExtrude"   

!------------------------------------------------------------------------------

    Mesh_out => AllocateMesh()
    isParallel = ( ParEnv % PEs > 1 )

    ! Create the division for the 1D mesh
    !--------------------------------------------
    CALL ExtrudedDivision(Vlist,nlev,Wtable)    
    CALL Info(Caller,'Extruding '//I2S(nlev)//' element layers on: '//TRIM(Mesh_in % Name),Level=10)
    in_levels = nlev-1
        
    ! Generate volume nodal points:
    ! -----------------------------
    n = Mesh_in % NumberOfNodes
    nnodes = (in_levels+2)*n
    gnodes = nnodes

    ALLOCATE( Mesh_out % Nodes % x(nnodes) )
    ALLOCATE( Mesh_out % Nodes % y(nnodes) )
    ALLOCATE( Mesh_out % Nodes % z(nnodes) )

    gelements = Mesh_in % NumberOfBulkElements

    ! There are some meshes with corrupted owners for 101 elements!
    ! This checks these nodes. 
    CALL CheckPointElementParents(Mesh_in)
    
    IF (isParallel) THEN
      PI_in  => Mesh_in % ParallelInfo
      PI_out => Mesh_out % ParallelInfo
      
      IF(.NOT. ASSOCIATED( PI_in ) ) CALL Fatal(Caller,'PI_in not associated!')
      IF(.NOT. ASSOCIATED( PI_out ) ) CALL Fatal(Caller,'PI_out not associated!')
            
      ALLOCATE(PI_out % NeighbourList(nnodes))
      ALLOCATE(PI_out % GInterface(nnodes))
      ALLOCATE(PI_out % GlobalDOFs(nnodes))

      IF(.NOT. ASSOCIATED( PI_in % NeighbourList ) ) THEN
        CALL Fatal(Caller,'Neighnours not associated!')
      END IF

      ! For unset neighbours just set the this partition to be the only owner
      DO i=1,Mesh_in % NumberOfNodes
        IF (.NOT.ASSOCIATED(PI_in % NeighbourList(i) % Neighbours)) THEN
          CALL AllocateVector(PI_in % NeighbourList(i) % Neighbours,1)
          PI_in % NeighbourList(i) % Neighbours(1) = ParEnv % Mype
        END IF
      END DO
          
      j=0
      DO i=1,Mesh_in % NumberOfNodes
        IF (PI_in % NeighbourList(i) % &
            Neighbours(1) == ParEnv % MyPE ) j=j+1
      END DO

      CALL MPI_ALLREDUCE(j,gnodes,1, &
           MPI_INTEGER,MPI_SUM,ELMER_COMM_WORLD,ierr)
      
      j=0
      DO i=1,Mesh_in % NumberOfBulkElements
        IF (Mesh_in % Elements(i) % PartIndex == ParEnv % MyPE) j=j+1
      END DO
      
      CALL MPI_ALLREDUCE(j,gelements,1, &
           MPI_INTEGER,MPI_SUM,ELMER_COMM_WORLD,ierr)
    END IF

    CALL Info(Caller,'Global count of original elements: '//I2S(gelements),Level=12)
    CALL Info(Caller,'Number of nodes for extruded mesh: '//I2S(nnodes),Level=12)

    ExtrudedCoord = ListGetInteger( CurrentModel % Simulation,'Extruded Coordinate Index', &
        Found, minv=1,maxv=3 )
    IF(.NOT. Found) ExtrudedCoord = MIN(3,Mesh_in % MeshDim + 1)
    CALL Info(Caller,'Extrusion in direction of dimension: '//I2S(ExtrudedCoord),Level=12)
    
    IF( ExtrudedCoord == 1 ) THEN
      ActiveCoord => Mesh_out % Nodes % x
    ELSE IF( ExtrudedCoord == 2 ) THEN
      ActiveCoord => Mesh_out % Nodes % y
    ELSE IF( ExtrudedCoord == 3 ) THEN
      ActiveCoord => Mesh_out % Nodes % z
    END IF

    PreserveBaseline = ListGetLogical( CurrentModel % Simulation,'Preserve Baseline',Found )

    CollectExtrudedBCs = ListGetLogical( CurrentModel % Simulation,'Extruded BCs Collect',Found )
    
    Rotate2Pi = .FALSE.
    Rotational = ListGetLogical( CurrentModel % Simulation,'Extruded Mesh Rotational',Found )    
    IF( Rotational ) THEN
      Rotate2Pi = ( ABS(MAXVAL(Wtable)-MINVAL(Wtable) - 2*PI) < 1.0d-3*PI )
      IF( Rotate2Pi ) CALL Info(Caller,'Perfoming full 2Pi rotation',Level=6)
    END IF

    ! This sets the BC layers.
    ! We honor the old way of assuming just bottom and top layer so the internal BCs are
    ! set as additional layers between.
    TmpLayers => ListGetIntegerArray( CurrentModel % Simulation,'Extruded BC Layers', GotBCLayers ) 
    IF( GotBCLayers ) THEN
      NoBCLayers = 2 + SIZE( TmpLayers )      
    ELSE
      NoBCLayers = 2
    END IF
    ALLOCATE(BCLayers(NoBCLayers))
    BCLayers(1) = 0
    BCLayers(NoBCLayers) = in_levels+1    
    IF( GotBCLayers ) THEN
      CALL Info(Caller,'There will be total of '//I2S(NoBCLayers)//' layers with BCs',Level=8)
      BCLayers(2:NoBCLayers-1) = TmpLayers
      DO i=1,NoBCLayers-1
        IF(BCLayers(i) >= BCLayers(i+1)) THEN
          CALL Fatal(Caller,'BC layers should be in increasing order')
        END IF
      END DO
    END IF    

    DoCount = .FALSE.
    BCCounter = 0 
    BaseLineLayer = 0
    baseline0 = 0
    IF( PreserveBaseline ) THEN
      BaseLineLayer = ListGetInteger( CurrentModel % Simulation,'Extruded Baseline Layer', Found )
      IF(.NOT. Found) BaseLineLayer = 1
      IF( BaseLineLayer > NoBCLayers ) THEN
        CALL Fatal(Caller,"'Extruded Baseline Layer' cannot exceed: "//I2S(NoBCLayers)) 
      END IF
      CALL Info(Caller,'Baseline will be set to layer '//I2S(BaselineLayer),Level=8)
    END IF
    
    max_body=0
    DO i=1,Mesh_in % NumberOfBulkElements
      max_body = MAX(max_body,Mesh_in % Elements(i) % Bodyid)
    END DO
    IF(isParallel) THEN
      j=max_body
      CALL MPI_ALLREDUCE(j,max_body,1,MPI_INTEGER,MPI_MAX,ELMER_COMM_WORLD,ierr)
    END IF
    CALL Info(Caller,'Maximum body index in original mesh: '//I2S(max_body),Level=6)

    max_bid0 = 0
    DO j=1,Mesh_in % NumberOfBoundaryElements
      k = j + Mesh_in % NumberOfBulkElements
      Elem_in => Mesh_in % Elements(k)
      IF(.NOT. ASSOCIATED(Elem_in % BoundaryInfo)) CYCLE
      bcind = Elem_in % BoundaryInfo % constraint 
      max_bid0 = MAX(max_bid0,bcind)
    END DO        
    IF(isParallel) THEN
      j = max_bid0
      CALL MPI_ALLREDUCE(j,max_bid0,1,MPI_INTEGER,MPI_MAX,ELMER_COMM_WORLD,ierr)
    END IF
    CALL Info(Caller,'Maximum boundary index in original mesh: '//I2S(max_bid0),Level=6)
    

    ! Create the nodes (and in parallel their global indexes).
    ! This assumes exacyly same distribution for each extruded node. 
    cnt=0
    DO i=0,in_levels+1

      ! If we rotate full 2Pi then we have natural closure!
      IF( Rotate2Pi ) THEN
        IF( i == in_levels+1) EXIT
      END IF      
      CurrCoord = Wtable( i ) 
      
      DO j=1,Mesh_in % NumberOfNodes

        cnt = cnt + 1

        Mesh_out % Nodes % x(cnt) = Mesh_in % Nodes % x(j) 
        Mesh_out % Nodes % y(cnt) = Mesh_in % Nodes % y(j) 
        Mesh_out % Nodes % z(cnt) = Mesh_in % Nodes % z(j) 

        ! Override the coordinate in the extruded direction by the value on the layer.
        ActiveCoord(cnt) = CurrCoord

        IF (isParallel) THEN
          PI_out % GInterface(cnt) = PI_in % GInterface(j)

          ALLOCATE(PI_out % NeighbourList(cnt) % Neighbours(&
               SIZE(PI_in % NeighbourList(j) % Neighbours)))
          PI_out % NeighbourList(cnt) % Neighbours = &
            PI_in % NeighbourList(j) % Neighbours
          PI_out % GlobalDOFs(cnt) = PI_in % GlobalDOFs(j)+i*gnodes
        END IF

      END DO
    END DO
    Mesh_out % NumberOfNodes = cnt
    Mesh_out % Nodes % NumberOfNodes = cnt

    ! For rotational geometry map the coordinates. 
    IF( Rotational ) THEN
      BLOCK
        REAL(KIND=DP) :: x,y,z,r        
        DO i=1,cnt          
          x = Mesh_out % Nodes % x(i)
          y = Mesh_out % Nodes % y(i)
          z = Mesh_out % Nodes % z(i)

          Mesh_out % Nodes % x(i) = COS(z) * x
          Mesh_out % Nodes % y(i) = SIN(z) * x
          Mesh_out % Nodes % z(i) = y
        END DO
      END BLOCK
    END IF
        
    ! Warn about 101 elements:
    ! -------------------------
    cnt101 = 0
    bcignored = 0
    DO i=Mesh_in % NumberOfBulkElements+1, &
        Mesh_in % NumberOfBulkElements+Mesh_in % NumberOfBoundaryElements
      Elem_in => Mesh_in % Elements(i)
      IF(Elem_in % TYPE % ElementCode == 101) cnt101 = cnt101 + 1
      IF(Elem_in % BoundaryInfo % Constraint == 0 ) bcignored = bcignored + 1
    END DO

    IF(isParallel) THEN
      j=cnt101
      CALL MPI_ALLREDUCE(j,cnt101,1,MPI_INTEGER,MPI_SUM,ELMER_COMM_WORLD,ierr)
      j=bcignored
      CALL MPI_ALLREDUCE(j,bcignored,1,MPI_INTEGER,MPI_SUM,ELMER_COMM_WORLD,ierr)
    END IF
       
    IF( bcignored > 0 ) THEN
      CALL Info(Caller,"WARNING: We are skipping '//I2S(bcignored)//&
          ' non-defined BC elements in extrusion!",Level=3)
    END IF       
    IF( cnt101 > 0 ) THEN
      CALL Info(Caller,"WARNING: Historically 101's were extruded as is, now they become 202's!",Level=3)
    END IF
    
    ! Compute total number of elements needed
    ! extruded bulk + extruded bc elements
    n = Mesh_in % NumberOfBulkElements + Mesh_in % NumberOfBoundaryElements
    totalnumberofelements = n*(in_levels+1) 
    IF(.NOT. Rotate2Pi ) THEN
      ! new layer bc's
      totalnumberofelements = totalnumberofelements + NoBCLayers * Mesh_in % NumberOfBulkElements
    END IF
    IF (PreserveBaseline) THEN
      ! additional baseline elements, if requested
      totalnumberofelements = totalnumberofelements + Mesh_in % NumberOfBoundaryElements
    END IF
    ALLOCATE(Mesh_out % Elements(totalnumberofelements))
    
    ! Initialize all elements to zero
    DO i = 1, totalnumberofelements
      Elem_out => Mesh_out % Elements(i)      
      Elem_out % DGDOFs = 0
      Elem_out % NDOFs = 0
      Elem_out % BodyId = 0
      Elem_out % DGIndexes => NULL()
      Elem_out % PDefs => NULL()
      Elem_out % EdgeIndexes => NULL()
      Elem_out % FaceIndexes => NULL()
      Elem_out % BubbleIndexes => NULL()
    END DO
    Mesh_out % MaxElementNodes = 0

    
    ! Generate volume bulk elements:
    ! ------------------------------
    n = Mesh_in % NumberOfNodes
    cnt=0
    DO i=0,in_levels
      DO j=1,Mesh_in % NumberOfBulkElements
        cnt = cnt+1
        Elem_in => Mesh_in % Elements(j)
        Elem_out => Mesh_out % Elements(cnt)

        !Elem_out % BodyId = Elem_in % BodyId
        Elem_out % BodyId = MapExtrudedMaterial(Vlist,Elem_in % BodyId,i+1)        

        Elem_out % PartIndex = Elem_in % PartIndex
               
        ! If we have internal BC layers then find the correct index for the body
        IF( NoBCLayers > 2 ) THEN
          DO k=1,NoBCLayers-1
            IF(i < BCLayers(k+1) ) EXIT
          END DO
          Elem_out % BodyId = Elem_out % BodyId + max_body*(k-1)
        END IF
        
        m = Elem_in % TYPE % NumberOfNodes
        ind(1:m) = Elem_in % NodeIndexes(1:m) + i*n

        IF( Rotate2Pi .AND. i==in_levels ) THEN
          ind(m+1:2*m) = Elem_in % NodeIndexes(1:m)
        ELSE
          ind(m+1:2*m) = Elem_in % NodeIndexes(1:m)+(i+1)*n
        END IF
        m = 2*m
                
        Elem_out % NDOFs = m
        Mesh_out % MaxElementNodes = MAX(Mesh_out % MaxElementNodes,m)

        SELECT CASE(m)
        CASE(4)
          Elem_out % TYPE => GetElementType(404)
          ! We need to reorder for the quad element!
          k = ind(3); ind(3)=ind(4); ind(4) = k
        CASE(6)
          Elem_out % TYPE => GetElementType(706)
        CASE(8)
          Elem_out % TYPE => GetElementType(808)
        END SELECT

        Elem_out % GElementIndex = Elem_in % GelementIndex + gelements*i

        Elem_out % ElementIndex = cnt
        ALLOCATE(Elem_out % NodeIndexes(m)) 
        Elem_out % NodeIndexes = ind(1:m)
      END DO
    END DO
    Mesh_out % NumberOfBulkElements = cnt
    CALL Info(Caller,'Number of extruded bulk elements: '//I2S(cnt),Level=8)
      
    ! Add side boundaries with the bottom mesh boundary id's:
    ! (or shift ids if preserving the baseline boundary)
    ! -------------------------------------------------------
    max_bid = 0
    bcoffset = 0
    IF( PreserveBaseline ) THEN
      CALL Info(Caller,'Preserving original '//I2S(max_bid0)//' BCs',Level=8)
      bcoffset = max_bid0
    END IF

    CALL Info(Caller,'First extruded boundary element index: '//I2S(cnt+1),Level=20)
    
    DO i=0,in_levels
      DO j=1,Mesh_in % NumberOfBoundaryElements
        k = j + Mesh_in % NumberOfBulkElements

        Elem_in => Mesh_in % Elements(k)
        bcind = Elem_in % BoundaryInfo % constraint

        ! Do not include BCs that are originally not activated
        IF(bcind==0) CYCLE

        cnt = cnt+1
        Elem_out => Mesh_out % Elements(cnt)  
        
        Elem_out = Elem_in

        Elem_out % ElementIndex = cnt        
        ALLOCATE(Elem_out % BoundaryInfo)
        Elem_out % BoundaryInfo = Elem_in % BoundaryInfo
        Elem_out % PartIndex = Elem_in % PartIndex
        
        ! Offset from possible baseline         
        bcind = bcind + bcoffset

        ! If we have internal BC layers then find the correct index for the body
        IF( NoBCLayers > 2 ) THEN
          DO k=1,NoBCLayers-1
            IF(i < BCLayers(k+1) ) EXIT
          END DO
          bcind = bcind + max_bid0*(k-1)
        END IF
        IF(DoCount .AND. bcind <= 100) BcCounter(bcind) = BcCounter(bcind) + 1
        
        Elem_out % BoundaryInfo % constraint = bcind
        max_bid = MAX(max_bid,bcind )

        m = Elem_in % TYPE % ElementCode / 100
        IF(m == 2) THEN
          Elem_out % NDOFs = 4
          ALLOCATE(Mesh_out % Elements(cnt) % NodeIndexes(4)) 

          ind(1) = Elem_in % NodeIndexes(1)+i*n
          ind(2) = Elem_in % NodeIndexes(2)+i*n
          IF( Rotate2Pi .AND. i==in_levels ) THEN
            ind(3) = Elem_in % NodeIndexes(2)
            ind(4) = Elem_in % NodeIndexes(1)
          ELSE
            ind(3) = Elem_in % NodeIndexes(2)+(i+1)*n
            ind(4) = Elem_in % NodeIndexes(1)+(i+1)*n
          END IF
          
          Elem_out % NodeIndexes = ind(1:4)
          Elem_out % TYPE => GetElementType(404)
        ELSE IF(m == 1) THEN
          Elem_out % NDOFs = 2
          ALLOCATE(Elem_out % NodeIndexes(2))
          
          ind(1) = Elem_in % NodeIndexes(1)+i*n
          ind(2) = Elem_in % NodeIndexes(1)+(i+1)*n

          Elem_out % NodeIndexes = ind(1:2)
          Elem_out % TYPE => GetElementType(202)
        ELSE
          CALL Fatal(Caller,'Invalid number of nodes: '//I2S(m))
        END IF

        IF( bcind <= CurrentModel % NumberOfBCs) THEN
          k = ListGetInteger(CurrentModel % BCs(bcind) % Values,'Body Id',Found)
          IF(Found) Elem_out % BodyId = k
        END IF

        IF(ASSOCIATED(Elem_in % BoundaryInfo % Left)) THEN
          l = Elem_in % BoundaryInfo % Left % ElementIndex
          Elem_out % BoundaryInfo % Left => &
             Mesh_out % Elements(Mesh_in % NumberOfBulkElements*i+l)
        END IF
        IF(ASSOCIATED(Elem_in % BoundaryInfo % Right)) THEN
          l = Elem_in % BoundaryInfo % Right % ElementIndex
          Elem_out % BoundaryInfo % Right => &
             Mesh_out % Elements(Mesh_in % NumberOfBulkElements*i+l)
        END IF

        ! Just check that we have correct parents. We had some issues here with
        ! corrupted initial meshes.
        BLOCK
          INTEGER :: ii,jj
          IF(ASSOCIATED(Elem_in % BoundaryInfo % Left)) THEN
            DO ii = 1, Elem_out % TYPE % NumberOfNodes
              jj = Elem_out % NodeIndexes(ii)
              IF( ALL( Elem_out % BoundaryInfo % Left % NodeIndexes /= jj ) ) THEN
                CALL Warn(Caller,'Node not available in left parent!')
              END IF
            END DO
          END IF
          IF(ASSOCIATED(Elem_in % BoundaryInfo % Right)) THEN
            DO ii = 1, Elem_out % TYPE % NumberOfNodes
              jj = Elem_out % NodeIndexes(ii)
              IF( ALL( Elem_out % BoundaryInfo % Right % NodeIndexes /= jj ) ) THEN
                CALL Warn(Caller,'Node not available in right parent!')
              END IF
            END DO
          END IF
        END BLOCK

        
      END DO
    END DO
        
    IF(isParallel) THEN
      j=max_bid
      CALL MPI_ALLREDUCE(j,max_bid,1,MPI_INTEGER,MPI_MAX,ELMER_COMM_WORLD,ierr)
    END IF
    CALL Info(Caller,'Largest bc index after extruded BCs: '//I2S(max_bid),Level=8)
    IF(DoCount) PRINT *,'BCInd1:',BcCounter(1:20)
    

    ! Add start and finish planes except if we have a full rotational symmetry
    IF(Rotate2Pi ) GOTO 100 
    
    ! Add bottom, top, and possible mid boundaries:
    ! ---------------------------------------------
    CALL Info(Caller,'First plane boundary element index: '//I2S(cnt+1),Level=20)
    bcoffset = max_bid
    DO k=1,NoBCLayers

      bclevel = BCLayers(k)

      IF( PreserveBaseline ) THEN
        ! Register the starting point for parents of baseline elements
        IF(k == BaselineLayer ) THEN
          baseline0 = cnt
          CALL Info(Caller,'Baseline elements parents start from element index: '//I2S(cnt),Level=8)
        END IF
      END IF

      DO i=1,Mesh_in % NumberOfBulkElements
        cnt=cnt+1
        
        Elem_in => Mesh_in % Elements(i)
        Elem_out => Mesh_out % Elements(cnt)

        Elem_out = Elem_in
        Elem_out % PartIndex = Elem_in % PartIndex
        
        m = Elem_in % TYPE % NumberOfNodes
        Elem_out % NDOFs = m
        ALLOCATE(Elem_out % NodeIndexes(m))        
        ALLOCATE(Elem_out % BoundaryInfo)
        
        Elem_out % BoundaryInfo % Right => NULL()
        IF( bclevel == in_levels+1 ) THEN
          Elem_out % BoundaryInfo % Left => &
              Mesh_out % Elements((bclevel-1) * Mesh_in % NumberOfBulkElements+i)          
        ELSE
          Elem_out % BoundaryInfo % Left => &
              Mesh_out % Elements(bclevel * Mesh_in % NumberOfBulkElements+i)
          IF(bclevel > 0 ) THEN
            ! for internal BCs add the 2nd parent also!
            Elem_out % BoundaryInfo % Right => &
                Mesh_out % Elements((bclevel-1) * Mesh_in % NumberOfBulkElements+i)
          END IF
        END IF

        IF( CollectExtrudedBCs ) THEN
          bcind = bcoffset + k
        ELSE
          bcind = bcoffset + (k-1)*max_body + Elem_in % BodyId
        END IF
        IF(DoCount .AND. bcind <= 100) BcCounter(bcind) = BcCounter(bcind) + 1
        
        max_bid = MAX(max_bid,bcind )
        
        Elem_out % BoundaryInfo % Constraint = bcind        
        Elem_out % BodyId = 0

        IF( bcind <= CurrentModel % NumberOfBCs) THEN
          j = ListGetInteger(CurrentModel % BCs(bcind) % Values,'Body Id',Found)
          IF(Found) Elem_out % BodyId = j
        END IF

        Elem_out % NodeIndexes = Elem_in % NodeIndexes + bclevel * n  

        Elem_out % ElementIndex = cnt
        Elem_out % TYPE => Elem_in % TYPE
      END DO
    END DO

    IF(isParallel) THEN
      j=max_bid
      CALL MPI_ALLREDUCE(j,max_bid,1,MPI_INTEGER,MPI_MAX,ELMER_COMM_WORLD,ierr)
    END IF
    CALL Info(Caller,'Largest bc index after layer BCs: '//I2S(max_bid),Level=8)
    IF(DoCount) PRINT *,'BCInd2:',BcCounter(1:20)

    
    ! If baseline preservation is requested, these will be
    ! available in the given layer with original bc tags.
    ! We do this at the end but still use the smallest (original)
    ! bc constraint indeces here.
    ! -------------------------------------------------------
    CALL Info(Caller,'First plane boundary element index: '//I2S(cnt+1),Level=20)
    IF (PreserveBaseline ) THEN
      DO j=1,Mesh_in % NumberOfBoundaryElements
        k = j + Mesh_in % NumberOfBulkElements
        
        Elem_In => Mesh_in % Elements(k)
        bcind = Elem_In % BoundaryInfo % Constraint 
        IF(bcind==0) CYCLE
        IF(DoCount .AND. bcind <= 100) BcCounter(bcind) = BcCounter(bcind) + 1
        
        cnt = cnt+1
        Elem_out => Mesh_out % Elements(cnt) 
        
        ALLOCATE(Elem_out % BoundaryInfo)
        Elem_out % BoundaryInfo = Elem_In % BoundaryInfo        
        Elem_out % BoundaryInfo % Constraint = bcind
        Elem_out % PartIndex = Elem_in % PartIndex
        
        Elem_out % TYPE => Elem_In % TYPE
        m = Elem_out % TYPE % ElementCode / 100
        Elem_out % NDOFs = m 

        k = BCLayers(BaselineLayer) * Mesh_in % NumberOfNodes
        ind(1:m) = Elem_in % NodeIndexes(1:m) + k
        
        ALLOCATE(Elem_out % NodeIndexes(m)) 
        Elem_out % NodeIndexes(1:m) = ind(1:m)

        Elem_out % ElementIndex = cnt
        
        IF(ASSOCIATED(Elem_In % BoundaryInfo % Left)) THEN
          l = Elem_in % BoundaryInfo % Left % ElementIndex + baseline0
          Elem_out % BoundaryInfo % Left => Mesh_out % Elements(l)
        END IF
        IF(ASSOCIATED(Elem_In % BoundaryInfo % Right)) THEN
          l = Elem_in % BoundaryInfo % Right % ElementIndex + baseline0
          Elem_out % BoundaryInfo % Right => Mesh_out % Elements(l)
        END IF
      END DO

      CALL Info(Caller,'Original baseline given by BCs: '//I2S(max_bid0))
    END IF
    IF(DoCount) PRINT *,'BCInd3:',BcCounter(1:20)
    CALL Info(Caller,'Last boundary element index: '//I2S(cnt),Level=20)

    DO i=1,cnt
      Elem_out => Mesh_out % Elements(i)
      IF( Elem_out % ElementIndex /= i) PRINT *,'mismatch: ',i,Elem_out % ElementIndex
    END DO

    
    
100 Mesh_out % NumberOfBoundaryElements = cnt-Mesh_out % NumberOfBulkElements
    
    Mesh_out % Name = Mesh_in % Name
    Mesh_out % DiscontMesh = Mesh_in % DiscontMesh
    Mesh_out % MaxElementDOFs  = Mesh_out % MaxElementNodes
    Mesh_out % Stabilize = Mesh_in % Stabilize

    Mesh_out % MeshDim = MIN(3, Mesh_in % MeshDim + 1)
    CurrentModel % Dimension = MIN( CurrentModel % Dimension+1, 3 )
   
    DEALLOCATE( BCLayers ) 

    ! Check whether the *.sif file has included enough BCs.
    ! If not then add some for convenience.
    j = 0
    DO i=Mesh_out % NumberOfBulkElements+1, &
        Mesh_out % NumberOfBulkElements+Mesh_out % NumberOfBoundaryElements
      Elem_out => Mesh_out % Elements(i)      
      bcind = Elem_out % BoundaryInfo % Constraint
      IF(bcind==0) CYCLE
      j = MAX(bcind,j)
    END DO
    CALL Info(Caller,'Maximum bc constraint in extruded mesh: '//I2S(j))
    IF( j > CurrentModel % NumberOfBCs ) THEN
      CALL AppendMissingBCs(CurrentModel,j)
    END IF

    CALL SetMeshSkew(Mesh_out, CurrentModel % Simulation )
    
    CALL PrepareMesh( CurrentModel, Mesh_out, isParallel )
    
    ExtrudedMeshName = ListGetString(CurrentModel % Simulation,'Extruded Mesh Name',Found)
    IF(Found) THEN
      IF( ParEnv % PEs > 1 ) THEN
        ! Or WriteMeshToDiskPartitioned ? 
        CALL WriteMeshToDisk2( CurrentModel, Mesh_out, ExtrudedMeshName, ParEnv % MyPe )
      ELSE        
        CALL WriteMeshToDisk(Mesh_out, ExtrudedMeshName)
      END IF
    END IF
    
  CONTAINS

    SUBROUTINE AppendMissingBCs(Model,maxbc)
       TYPE(Model_t) :: Model
       INTEGER :: maxbc
       
       INTEGER :: i, NoBCs, tag 
       TYPE(BoundaryConditionArray_t), POINTER :: OldBCs(:) => NULL()
       
       NoBcs = Model % NumberOfBCs
       IF(NoBCs >= maxbc ) RETURN

       CALL Info(Caller,'Generating '//I2S(maxbc-NoBCs)//' dummy list BCs for convenience!',Level=5)
       
       OldBCs => Model % BCs(:)

       NULLIFY( Model % BCs )
       ALLOCATE( Model % BCs(maxbc) )
       
       DO i=1,NoBCs
         Model % BCs(i) % Values => OldBCs(i) % Values        
         tag = OldBCs(i) % Tag
         IF(tag == 0) tag = i 
         Model % BCs(i) % Tag = tag
       END DO
       IF (ASSOCIATED(OldBCs) .AND. NoBCs > 0) DEALLOCATE( OldBCs ) 
       DO i=NoBCs+1,maxbc
         Model % BCs(i) % Tag = i
       END DO
       DO i=1,maxbc
         IF(.NOT.ASSOCIATED(Model % BCs(i) % Values) .OR. i > NoBCs) THEN
           Model % BCs(i) % Values => ListAllocate()
           CALL ListAddString( Model % BCs(i) % Values,'Name','BC'//I2S(i))
         END IF
       END DO
       Model % NumberOfBCs = maxbc
       
     END SUBROUTINE AppendMissingBCs
    
!------------------------------------------------------------------------------
  END FUNCTION MeshExtrude
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> As the previous one except the extrusion is done in parallel for single meshes
!> that each take an internal in the extruded direction. This affects the coordinates
!> but also the communication pattern. A separate routine was made in order to avoid
!> introducing of bugs as the internal extrusion is a widely used feature. 
!------------------------------------------------------------------------------
  FUNCTION MeshExtrudeSlices(Mesh_in, Vlist) RESULT(Mesh_out)
!------------------------------------------------------------------------------
    TYPE(Mesh_t), TARGET :: Mesh_in
    TYPE(Mesh_t), POINTER :: Mesh_out
    TYPE(ValueList_t), POINTER :: Vlist
!------------------------------------------------------------------------------
    CHARACTER(:), ALLOCATABLE :: ExtrudedMeshName
    INTEGER :: i,j,k,l,n,m,cnt,ind(8),bid,max_bid,l_n,max_body,bcid,&
        ExtrudedCoord,dg_n,totalnumberofelements,lastbc
    INTEGER, POINTER :: pInds(:)
    TYPE(ParallelInfo_t), POINTER :: PI_in, PI_out
    TYPE(Element_t), POINTER :: Element
    INTEGER :: nnodes,gnodes,gelements,ierr,nlev,ilev,&
        nParMesh,nParExt,OrigPart,ElemCode,bodyid, newbcs
    LOGICAL :: isParallel, SingleIn, Found, TopBC, BotBC, &
        CollectExtrudedBCs, SeparateSlices, CreateInternalBCs, InternalBC
    INTEGER,ALLOCATABLE :: ChildBCs(:)
    REAL(KIND=dp)::CurrCoord 
    REAL(KIND=dp), POINTER :: ActiveCoord(:)
    REAL(KIND=dp), ALLOCATABLE :: Wtable(:)
    CHARACTER(*), PARAMETER :: Caller="MeshExtrudeSlices"   
!------------------------------------------------------------------------------

    ! The historical choice in_levels in annoying when we want to split the divisions.
    
    IF( ListGetLogical( CurrentModel % Simulation,'Preserve Baseline',Found ) ) &
        CALL Fatal(Caller,'The slice version cannot handle "Preserve Baseline"!')
    
    IF( ListGetLogical( CurrentModel % Simulation,'Extruded Mesh Rotational',Found ) ) &
        CALL Fatal(Caller,'The slice version cannot handle "Extruded Mesh Rotational"!')    
    
    isParallel = ( ParEnv % PEs > 1 )
    SingleIn = Mesh_in % SingleMesh
    
    ! Create the division for the 1D mesh
    !--------------------------------------------
    CALL ExtrudedDivision(Vlist,nlev,Wtable)    
    CALL Info(Caller,'Extruding '//I2S(nlev)//' element layers on: '//TRIM(Mesh_in % Name),Level=10)

    
    ! In parallel let us pick only our own share of the
    ! division. This logic makes it possible to have nonuniform divisions easily.
    ! The number of element layers is evenly distributed among partitions. 
    !-------------------------------------------------------------------------------
    IF( isParallel ) THEN
      nParExt = ParEnv % PEs 
      nParMesh = ListGetInteger( CurrentModel % Simulation,'Parallel Mesh Modulo',Found)
      IF(.NOT. Found) THEN
        nParMesh = 1
        IF(.NOT. SingleIn ) THEN
          CALL Fatal(Caller,'This routine expects either Mesh Modulo or Single Mesh!')
        END IF
      END IF
      
      nParExt = nParExt / nParMesh                    
      IF( MODULO(nlev,nParExt) /= 0 ) THEN
        CALL Fatal(Caller,'Number of element layers '//I2S(nlev)//&
            ' not divisible by '//I2S(ParEnv % PEs))
      END IF
      nlev = nlev / nParExt
      IF(nlev < 2) THEN
        CALL Fatal(Caller,'At least two element layers needed in each partition!')
      END IF
      ilev = ( ParEnv % MyPe / nParMesh ) * nlev
      Wtable(0:nlev) = Wtable(ilev:nlev+ilev) 
    ELSE
      nParExt = 1
      nParMesh = 1 
      ilev = 0
    END IF
        
    ! Allocate extruded mesh:
    ! We do this only after splitting the division.
    ! ---------------------------------------------
    n = Mesh_in % NumberOfNodes
    nnodes = (nlev+1)*n

    Mesh_out => AllocateMesh()
    ALLOCATE( Mesh_out % Nodes % x(nnodes) )
    ALLOCATE( Mesh_out % Nodes % y(nnodes) )
    ALLOCATE( Mesh_out % Nodes % z(nnodes) )
    
    gnodes = Mesh_in % NumberOfNodes
    gelements = Mesh_in % NumberOfBulkElements

    Mesh_out % SingleMesh = .FALSE.
    SeparateSlices = .FALSE.
    CreateInternalBCs = .FALSE.
    
    ExtrudedCoord = ListGetInteger( CurrentModel % Simulation,'Extruded Coordinate Index', &
        Found, minv=1,maxv=3 )
    IF(.NOT. Found) ExtrudedCoord = 3 
    
    IF( ExtrudedCoord == 1 ) THEN
      ActiveCoord => Mesh_out % Nodes % x
    ELSE IF( ExtrudedCoord == 2 ) THEN
      ActiveCoord => Mesh_out % Nodes % y
    ELSE IF( ExtrudedCoord == 3 ) THEN
      ActiveCoord => Mesh_out % Nodes % z
    END IF

    CollectExtrudedBCs = ListGetLogical( CurrentModel % Simulation,'Extruded BCs Collect',Found )
    
    IF (isParallel) THEN
      PI_in  => Mesh_in % ParallelInfo
      PI_out => Mesh_out % ParallelInfo

      IF(.NOT. ASSOCIATED( PI_in ) ) CALL Fatal(Caller,'PI_in not associated!')
      IF(.NOT. ASSOCIATED( PI_out ) ) CALL Fatal(Caller,'PI_out not associated!')
            
      SeparateSlices = ListGetLogical( CurrentModel % Simulation,'Extruded Mesh Slices Separate',Found )
      CreateInternalBCs = ListGetLogical( CurrentModel % Simulation,'Extruded BCs Internal',Found ) 
      IF(.NOT. Found) CreateInternalBCs = SeparateSlices 
      
      ALLOCATE(PI_out % NeighbourList(nnodes))
      ALLOCATE(PI_out % GInterface(nnodes))
      ALLOCATE(PI_out % GlobalDOFs(nnodes))

      IF(.NOT. SingleIn ) THEN
        IF(.NOT. ASSOCIATED( PI_in % NeighbourList ) ) THEN
          CALL Fatal(Caller,'Neighnours not associated in initial mesh!')
        END IF

        ! Count own nodes
        j=0
        DO i=1,Mesh_in % NumberOfNodes
          IF(.NOT. ASSOCIATED(PI_in % NeighbourList(i) % Neighbours ) ) THEN
            j = j + 1
          ELSE IF (PI_in % NeighbourList(i) % Neighbours(1) == ParEnv % MyPE ) THEN
            j=j+1
          END IF
        END DO
        CALL MPI_ALLREDUCE(j,gnodes,1, &
            MPI_INTEGER,MPI_SUM,ELMER_COMM_WORLD,ierr)
        gnodes = gnodes / nParExt
        
        j=0
        DO i=1,Mesh_in % NumberOfBulkElements
          IF (Mesh_in % Elements(i) % PartIndex == ParEnv % MyPE) j=j+1
        END DO
        CALL MPI_ALLREDUCE(j,gelements,1, &
            MPI_INTEGER,MPI_SUM,ELMER_COMM_WORLD,ierr)
        gelements = gelements / nParExt

        !PRINT *,'nParExt:',ParEnv % Mype, nParExt, nParMesh, gnodes,gelements        
      END IF

      Mesh_out % ParallelInfo % NothingShared = ( SingleIn .AND. SeparateSlices )
    END IF

    CALL Info(Caller,'Number of nodes in layer: '//I2S(gnodes),Level=12)
    CALL Info(Caller,'Number of elements in layer: '//I2S(gelements),Level=12)
    
    cnt=0
    DO i=0,nlev

      CurrCoord = Wtable(i)  
      
      DO j=1,Mesh_in % NumberOfNodes

        cnt = cnt + 1

        Mesh_out % Nodes % x(cnt) = Mesh_in % Nodes % x(j) 
        Mesh_out % Nodes % y(cnt) = Mesh_in % Nodes % y(j) 
        Mesh_out % Nodes % z(cnt) = Mesh_in % Nodes % z(j) 

        ! Override the coordinate in the extruded direction by the value on the layer.
        ActiveCoord(cnt) = CurrCoord

        IF (isParallel) THEN
          m = 1
          IF( nParMesh > 1 ) THEN
            IF( ASSOCIATED( PI_in % NeighbourList(j) % Neighbours ) ) THEN
              m = SIZE(PI_in % NeighbourList(j) % Neighbours)
            END IF
          END IF
          IF( SeparateSlices ) THEN
            k = m
          ELSE IF(i==0 .AND. ParEnv % MyPe > (nParMesh-1) ) THEN            
            k = 2*m
          ELSE IF(i==nlev .AND. ParEnv % MyPe < ParEnv % PEs- nParMesh ) THEN
            k = 2*m
          ELSE
            k = m
          END IF

          ALLOCATE(PI_out % NeighbourList(cnt) % Neighbours(k))
          PI_out % GInterface(cnt) = (k>1)
                    
          DO k=1,m
            IF(m>1) THEN
              OrigPart = PI_in % NeighbourList(j) % Neighbours(k)
            ELSE
              OrigPart = ParEnv % MyPe
            END IF                       

            IF( SeparateSlices ) THEN
              l = ( ParEnv % MyPe / nParMesh ) * ( nlev + 1 ) * gnodes
            ELSE
              l = ilev * gnodes
            END IF
                            
            IF(SingleIn) THEN
              l = l + j + i * gnodes
            ELSE
              l = l + MODULO(PI_in % GlobalDOFs(j)-1,gnodes)+1 + i * gnodes 
            END IF
            PI_out % GlobalDOFs(cnt) = l

            IF( SeparateSlices ) THEN
              PI_out % NeighbourList(cnt) % Neighbours(k) = OrigPart               
            ELSE IF(i==0 .AND. ParEnv % MyPe > nParMesh-1 ) THEN
              PI_out % NeighbourList(cnt) % Neighbours(2*k-1) = OrigPart
              PI_out % NeighbourList(cnt) % Neighbours(2*k) = OrigPart-1            
            ELSE IF(i==nlev .AND. ParEnv % MyPe < ParEnv % PEs-nParMesh ) THEN
              PI_out % NeighbourList(cnt) % Neighbours(2*k-1) = OrigPart+1
              PI_out % NeighbourList(cnt) % Neighbours(2*k) = OrigPart
            ELSE
              PI_out % NeighbourList(cnt) % Neighbours(k) = OrigPart 
            END IF                       
          END DO
          
        END IF
      END DO
    END DO
    
    Mesh_out % NumberOfNodes = cnt
    Mesh_out % Nodes % NumberOfNodes = cnt

    ! Calculate exactly and allocate the number of extruded elements
    n = Mesh_in % NumberOfBulkElements + Mesh_in % NumberOfBoundaryElements
    totalnumberofelements = n*nlev

    IF( CreateInternalBCs ) THEN
      totalnumberofelements = &
          totalnumberofelements + 2 * Mesh_in % NumberOfBulkElements 
    ELSE
      IF( ParEnv % MyPe < nParMesh ) totalnumberofelements = &
          totalnumberofelements + Mesh_in % NumberOfBulkElements 
      IF( ParEnv % MyPe >= ParEnv % PEs-nParMesh ) totalnumberofelements = &
          totalnumberofelements + Mesh_in % NumberOfBulkElements 
    END IF
      
    ALLOCATE(Mesh_out % Elements(totalnumberofelements))

    
    ! Generate volume bulk elements:
    ! ------------------------------
    Mesh_out % MaxElementNodes = 0
    n = Mesh_in % NumberOfNodes
    cnt=0; dg_n  = 0

    DO i=0,nlev-1
      DO j=1,Mesh_in % NumberOfBulkElements

        cnt = cnt+1
        Element => Mesh_out % Elements(cnt)
        Element = Mesh_in % Elements(j)

        bodyid = Element % BodyId
        Element % BodyId = MapExtrudedMaterial(Vlist,bodyid,ilev+i+1)        
        
        l_n = Mesh_in % Elements(j) % TYPE % NumberOfNodes
        ind(1:l_n) = Mesh_in % Elements(j) % NodeIndexes(1:l_n)+i*n
        ind(l_n+1:2*l_n) = Mesh_in % Elements(j) % NodeIndexes(1:l_n)+(i+1)*n
        l_n = 2*l_n
        Element % NDOFs = l_n
        Mesh_out % MaxElementNodes = MAX(Mesh_out % MaxElementNodes,l_n)

        SELECT CASE(l_n)
        CASE(6)
          Element % TYPE => GetElementType(706)
        CASE(8)
          Element % TYPE => GetElementType(808)
        END SELECT

        IF( isParallel ) THEN
          IF(SingleIn) THEN
            l = j + (ilev+i) * gelements
          ELSE
            l = MODULO(Mesh_in % Elements(j) % GElementIndex-1,gelements)+1 + (ilev+i) * gelements 
          END IF
          Element % GElementIndex = l
        ELSE
          Element % GElementIndex = cnt
        END IF
          
        Element % ElementIndex = cnt
        ALLOCATE(Element % NodeIndexes(l_n)) 
        Element % NodeIndexes = ind(1:l_n)
      END DO
    END DO
    Mesh_out % NumberOfBulkElements = cnt

    
    ! Add side boundaries with the bottom mesh boundary id's:
    ! -------------------------------------------------------
    max_bid = 0
    DO i=0,nlev-1
      DO j=1,Mesh_in % NumberOfBoundaryElements
        k = j + Mesh_in % NumberOfBulkElements

        cnt=cnt+1

        Element => Mesh_out % Elements(cnt)
        Element = Mesh_in % Elements(k)        
        ALLOCATE(Element % BoundaryInfo)

        Element % BoundaryInfo = Mesh_in % Elements(k) % BoundaryInfo

        bid = Mesh_in % Elements(k) % BoundaryInfo % Constraint
        max_bid = MAX(max_bid, bid )

        IF(ASSOCIATED(Mesh_in % Elements(k) % BoundaryInfo % Left)) THEN
          l = Mesh_in % Elements(k) % BoundaryInfo % Left % ElementIndex
          Element % BoundaryInfo % Left => &
              Mesh_out % Elements(Mesh_in % NumberOfBulkElements*i+l)
        END IF
        IF(ASSOCIATED(Mesh_in % Elements(k) % BoundaryInfo % Right)) THEN
          l = Mesh_in % Elements(k) % BoundaryInfo % Right % ElementIndex
          Element % BoundaryInfo % Right => &
             Mesh_out % Elements(Mesh_in % NumberOfBulkElements*i+l)
        END IF

        ElemCode = Mesh_in % Elements(k) % TYPE % ElementCode        
        m = 2*MODULO(ElemCode,100)        
        Element % NDOFs = m
        ALLOCATE(Element % NodeIndexes(m))
        pInds => Element % NodeIndexes
               
        IF(ElemCode == 202) THEN
          pInds(1) = Mesh_in % Elements(k) % NodeIndexes(1)+i*n
          pInds(2) = Mesh_in % Elements(k) % NodeIndexes(2)+i*n
          pInds(3) = Mesh_in % Elements(k) % NodeIndexes(2)+(i+1)*n
          pInds(4) = Mesh_in % Elements(k) % NodeIndexes(1)+(i+1)*n
          Mesh_out % Elements(cnt) % TYPE => GetElementType(404)
        ELSE IF(ElemCode == 101 ) THEN
          pInds(1) = Mesh_in % Elements(k) % NodeIndexes(1) +i*n
          pInds(2) = Mesh_in % Elements(k) % NodeIndexes(1) +(i+1)*n
        ELSE
          CALL Fatal(Caller,'Cannot extrude boundary element: '//I2S(ElemCode))
        END IF
        Element % ElementIndex = cnt
      END DO
    END DO

    IF(.NOT. SingleIn .AND. isParallel) THEN
      j=max_bid
      CALL MPI_ALLREDUCE(j,max_bid,1, &
          MPI_INTEGER,MPI_MAX,ELMER_COMM_WORLD,ierr)
    END IF
   
    CALL Info(Caller,'First Extruded BC set to: '//I2S(max_bid+1),Level=6)
    lastbc = max_bid+1

    max_body=0
    DO i=1,Mesh_in % NumberOfBulkElements
      max_body = MAX(max_body,Mesh_in % Elements(i) % Bodyid)
    END DO
    IF(.NOT. SingleIn .AND. isParallel) THEN
      j=max_body
      CALL MPI_ALLREDUCE(j,max_body,1, &
          MPI_INTEGER,MPI_MAX,ELMER_COMM_WORLD,ierr)
    END IF

    IF( CollectExtrudedBCs ) THEN
      CALL Info(Caller,'Number of new BCs for each layer: 1',Level=6)
    ELSE
      CALL Info(Caller,'Number of new BCs for each layer: '//I2S(max_body),Level=6)
    END IF
    
    IF( CollectExtrudedBCs ) THEN
      newbcs = 2
    ELSE
      newbcs = 2 * max_body
    END IF

    IF( CreateInternalBCs ) THEN
      CALL Info(Caller,'Internal bottom boundary: '//I2S(max_bid+newbcs+1),Level=6)
      CALL Info(Caller,'Internal top boundary: '//I2S(max_bid+newbcs+2),Level=6)
    END IF
    
    ALLOCATE(ChildBCs(2*max_body))
    ChildBCs = -1
           
    ! Add bottom boundary:
    ! --------------------
    IF( ParEnv % PEs == 1 .OR. ParEnv % MyPe < nParMesh .OR. CreateInternalBCs ) THEN  
      InternalBC = (ParEnv % PEs > 1 .AND. ParEnv % MyPe >= nParMesh )
      DO i=1,Mesh_in % NumberOfBulkElements
        cnt=cnt+1
        Element => Mesh_out % Elements(cnt) 
        
        Element = Mesh_in % Elements(i)

        l_n = Mesh_in % Elements(i) % TYPE % NumberOfNodes
        Element % NDOFs = l_n

        ALLOCATE(Element % BoundaryInfo)
        Element % BoundaryInfo % Left => Mesh_out % Elements(i)
        Element % BoundaryInfo % Right => NULL()

        bodyid = Mesh_in % Elements(i) % BodyId                
        IF( InternalBC ) THEN
          bcid = max_bid + newbcs + 1
        ELSE IF( CollectExtrudedBCs ) THEN
          bcid = max_bid + 1
        ELSE
          bcid = max_bid + bodyid
        END IF
        Element % BoundaryInfo % Constraint = bcid

        IF(.NOT. InternalBC) ChildBCs(2*bodyid-1) = bcid
        lastbc = MAX(lastbc,bcid)

        Element % BodyId = 0
        IF( bcid <= CurrentModel % NumberOfBCs) THEN
          j = ListGetInteger(CurrentModel % BCs(bcid) % Values,'Body Id',Found)
          IF(Found) Element % BodyId = j
        END IF

        ALLOCATE(Element % NodeIndexes(l_n))
        Element % NodeIndexes = Mesh_in % Elements(i) % NodeIndexes
        Element % ElementIndex = cnt
        Element % TYPE => Mesh_in % Elements(i) % TYPE
      END DO
    END IF

    
    ! Add top boundary:
    ! -----------------
    IF( ParEnv % PEs == 1 .OR. ParEnv % MyPe >= ParEnv % PEs - nParMesh .OR. CreateInternalBCs ) THEN
      InternalBC = (ParEnv % PEs > 1 .AND. ParEnv % MyPe < ParEnv % PEs - nParMesh )
      DO i=1,Mesh_in % NumberOfBulkElements
        cnt=cnt+1
        Element => Mesh_out % Elements(cnt) 
        
        Element = Mesh_in % Elements(i)

        l_n = Mesh_in % Elements(i) % TYPE % NumberOfNodes
        Element % NDOFs = l_n

        ALLOCATE(Element % BoundaryInfo)
        Element % BoundaryInfo % Left => &
            Mesh_out % Elements((nlev-1)*Mesh_in % NumberOfBulkElements+i)
        Element % BoundaryInfo % Right => NULL()
        
        bodyid = Mesh_in % Elements(i) % BodyId                
        IF( InternalBC ) THEN
          bcid = max_bid + newbcs + 2
        ELSE IF( CollectExtrudedBCs ) THEN
          bcid = max_bid + 2
        ELSE
          bcid = max_bid + bodyid + max_body
        END IF
        Element % BoundaryInfo % Constraint = bcid

        IF(.NOT. InternalBC) ChildBCs(2*bodyid) = bcid 
        lastbc = MAX(lastbc,bcid)
        
        Element % BodyId = 0
        IF( bcid<=CurrentModel % NumberOfBCs) THEN
          j = ListGetInteger(CurrentModel % BCs(bcid) % Values,'Body Id',Found)
          IF(Found) Element % BodyId = j
        END IF

        ALLOCATE(Element % NodeIndexes(l_n))
        Element % NodeIndexes = Mesh_in % Elements(i) % NodeIndexes+nlev*n
        Element % ElementIndex = cnt
        Element % TYPE => Mesh_in % Elements(i) % TYPE
      END DO
    END IF

    IF(.NOT. SingleIn .AND. isParallel) THEN
      j=lastbc
      CALL MPI_ALLREDUCE(j,lastbc,1, &
          MPI_INTEGER,MPI_MAX,ELMER_COMM_WORLD,ierr)
    END IF
    CALL Info(Caller,'Last Extruded BC set to: '//I2S(lastbc),Level=6)
    
    IF( cnt /= totalnumberofelements ) THEN
      CALL Fatal(Caller,'Mismatch between allocated and set elements: '//&
          I2S(totalnumberofelements)//' vs. '//I2S(cnt))
    END IF

    ! Set some unset stuff to be on the safe side
    DO i=1,cnt
      Element => Mesh_out % Elements(i)
      Element % DGDOFs = 0
      Element % DGIndexes => NULL()
      Element % PDefs => NULL()
      Element % EdgeIndexes => NULL()
      Element % FaceIndexes => NULL()
      Element % BubbleIndexes => NULL()
    END DO
         
    Mesh_out % NumberOfBoundaryElements = cnt - Mesh_out % NumberOfBulkElements
    
    Mesh_out % Name = Mesh_in % Name
    Mesh_out % DiscontMesh = Mesh_in % DiscontMesh
    Mesh_out % MaxElementDOFs = Mesh_out % MaxElementNodes
    Mesh_out % Stabilize = Mesh_in % Stabilize
    Mesh_out % MeshDim = 3
    CurrentModel % DIMENSION = 3


    ! Let us mark the child BCs to the bodies that they originate from.
    BLOCK
      INTEGER, POINTER :: TmpPair(:), TmpBCs(:) 
      TYPE(ValueList_t), POINTER :: vList

      ALLOCATE(TmpBCs(2*max_body))
      TmpBCs = ChildBCs

      IF( ParEnv % PEs > 1 ) THEN
        CALL MPI_ALLREDUCE(TmpBCs,ChildBCs,2*max_body, &
            MPI_INTEGER,MPI_MAX,ELMER_COMM_WORLD,ierr)
      END IF

      DO i=1,CurrentModel % NumberOfBodies
        vList => CurrentModel % Bodies(i) % Values
        IF( ASSOCIATED(vList) ) THEN
          NULLIFY(TmpPair)
          ALLOCATE(TmpPair(2))
          TmpPair(1) = ChildBCs(2*i-1)
          TmpPair(2) = ChildBCs(2*i)
          CALL ListAddIntegerArray(vList,'Extruded Child BCs',2,TmpPair)

          IF( InfoActive(10) ) THEN
            CALL Info(Caller,'Setting Body '//I2S(i)//' "Extruded Child BCs" to '&
                //I2S(TmpPair(1))//' '//I2S(TmpPair(2)))
          END IF
          NULLIFY(TmpPair)
        END IF
      END DO

      DEALLOCATE(TmpBCs)
    END BLOCK
      
    CALL SetMeshSkew(Mesh_out, CurrentModel % Simulation )
    
    ExtrudedMeshName = ListGetString(CurrentModel % Simulation,'Extruded Mesh Name',Found)
    IF(Found) THEN
      IF( ParEnv % PEs == 1 ) THEN
        CALL WriteMeshToDisk(Mesh_out, ExtrudedMeshName)
      ELSE
        CALL WriteMeshToDisk2(CurrentModel, Mesh_out, ExtrudedMeshName, ParEnv % MyPe )
      END IF
    END IF

    CALL PrepareMesh( CurrentModel, Mesh_out, isParallel )
    
!------------------------------------------------------------------------------
  END FUNCTION MeshExtrudeSlices
!------------------------------------------------------------------------------


END MODULE MeshExtrusion
