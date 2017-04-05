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

!------------------------------------------------------------------------------
!> Module for saving result in the old DX format 
!------------------------------------------------------------------------------
   MODULE DXFile

      USE MeshUtils
      USE ElementDescription
      
      IMPLICIT NONE
      !    PRIVATE
      SAVE
      
      PUBLIC :: WriteDXFiles
      
      INTEGER, PARAMETER :: MAX_VERTIX = 4, MAX_PART_ELEM = 21
      
    CONTAINS
      
      SUBROUTINE WriteDXFiles( Prefix, Model, SubtractDisp, nTime )
        CHARACTER(LEN=*), INTENT(IN) :: Prefix
        TYPE(Model_t) :: Model 
        LOGICAL, INTENT(IN) :: SubtractDisp
        INTEGER, INTENT(IN) :: nTime
        TYPE(Variable_t), POINTER :: Var,Var1
        CHARACTER(LEN=512) :: str
        INTEGER :: i
        INTEGER, PARAMETER :: MasterUnit = 58
        
        IF ( nTime == 1) THEN
          CALL WriteGrid( Prefix, Model, SubtractDisp )
          
          OPEN( MasterUnit, FILE = Prefix // "Master.dx", STATUS="unknown" )
          WRITE( MasterUnit, '(A)') 'object "group" class group'
        END IF
        
        Var => Model % Variables
        DO WHILE( ASSOCIATED( Var ) )
          IF ( .NOT.Var % Output ) THEN
            Var => Var % Next
            CYCLE
          END IF
          
          IF ( SIZE( Var % Values ) == Var % DOFs ) THEN
            Var => Var % Next
            CYCLE
          END IF
          
          SELECT CASE( Var % Name )
            
          CASE( 'mesh update' )
            Var1 => Model % Variables
            
            DO WHILE( ASSOCIATED( Var1 ) )
              IF ( TRIM( Var1 % Name ) == 'displacement' ) EXIT
              Var1 => Var1 % Next
            END DO
            
            IF ( .NOT.ASSOCIATED( Var1 ) ) THEN
              CALL WriteVariable( "MeshUpdate", Var, &
                  Model % NumberOfNodes, Var % DOFs, 0,  &
                  nTime, MasterUnit, Prefix )
            END IF
            
          CASE( 'mesh update 1','mesh update 2', 'mesh update 3' )
            
          CASE( 'displacement' )
            CALL WriteDisplacement( Var, Model, nTime, MasterUnit, Prefix )
          CASE( 'displacement 1','displacement 2','displacement 3' )
            
          CASE( 'flow solution' )
            CALL WriteVariable( "Velocity", Var, Model % NumberOfNodes, &
                Var % DOFs-1, 0, nTime, MasterUnit, Prefix )
            CALL WriteVariable( "Pressure", Var, Model % NumberOfNodes, 1, &
                Var % DOFs-1, nTime, MasterUnit, Prefix )
          CASE( 'velocity 1','velocity 2','velocity 3','pressure' )
            
          CASE( 'magnetic field' )
            CALL WriteVariable( "MagField", Var, Model % NumberOfNodes, &
                Var % DOFs, 0, nTime, MasterUnit, Prefix )
          CASE( 'magnetic field 1','magnetic field 2', 'magnetic field 3' )
            
          CASE( 'electric current' )
            CALL WriteVariable( "Current", Var, Model % NumberOfNodes, &
                Var % DOFs, 0, nTime, MasterUnit, Prefix )
          CASE('electric current 1','electric current 2','electric current 3')
            
          CASE( 'coordinate 1','coordinate 2','coordinate 3' )
            
          CASE( 'magnetic flux density' )
            CALL WriteVariable( "MagneticFlux", Var, Model % NumberOfNodes,&
                Var % DOFs, 0, nTime, MasterUnit, Prefix )
          CASE( 'magnetic flux density 1','magnetic flux density 2', &
              'magnetic flux density 3' )
          CASE DEFAULT
            DO i=1,Var % NameLen
              str(i:i) = Var % Name(i:i)
              IF( str(i:i) == ' ' ) str(i:i) = '_'
            END DO
            str(1:1) = CHAR(ICHAR(str(1:1))-ICHAR('a')+ICHAR('A'))
            
            CALL WriteVariable( TRIM(str), Var, Model % NumberOfNodes, &
                Var % DOFs, 0,  nTime, MasterUnit, Prefix )
          END SELECT
          Var => Var % Next
        END DO
        
        IF( nTime == 1) THEN
          CLOSE( MasterUnit )
        END IF
      END SUBROUTINE WriteDXFiles
      

      SUBROUTINE WriteVariable( VarName, Var, nNodes, SelfDOF, Offset, nTime, &
          MasterUnit, Prefix )
        CHARACTER(*), INTENT(IN) :: VarName
        TYPE(Variable_t), INTENT(IN) :: Var
        INTEGER, INTENT(IN) :: nNodes, SelfDOF, Offset, nTime, MasterUnit
        CHARACTER(*), INTENT(IN) :: Prefix
        INTEGER :: i, j, k
        INTEGER :: FUnit
        CHARACTER(MAX_NAME_LEN) :: FName, RelativeFName, MeshFile
        CHARACTER(7) :: VType
        
        FUnit = MasterUnit + 1
        FName = Prefix // VarName // ".dx"
        
        i = INDEX(FName, '/', BACK=.TRUE.)
        RelativeFName = FName(i+1:)
        
        i = INDEX(Prefix, '/', BACK=.TRUE.)
        MeshFile = Prefix(i+1:) // "Mesh.dx"
        
        IF( nTime == 1 ) THEN
          
            IF ( SelfDOF == 1 ) THEN
              VType = "-scalar"
            ELSE
              VType = "-vector"
            END IF
            WRITE( MasterUnit, '(A)' ) 'member "' // VarName // VType     &
                // '" value file "' // TRIM(RelativeFName) // '","'  &
                // VarName // 'series' // '"'
            
            OPEN( FUnit, FILE=FName, STATUS="unknown" )
          ELSE
            OPEN( FUnit, FILE=FName, STATUS="old", POSITION="append" )
            DO i = 1, 2+(nTime-1)*7
              BACKSPACE FUnit   ! Yuck!
            END DO
          END IF
          
          IF( SelfDOF == 1 )THEN
            WRITE( FUnit,'("object ",I0," class array type double rank 0 &
                &items ",I0," data follows")') nTime, nNodes
          ELSE
            WRITE( FUnit,'("object ",I0," class array type double rank 1 shape &
                &",I0," items ",I0," data follows")') nTime, & 
                SelfDOF, nNodes
          END IF
          
          DO i = 1, nNodes
            k = i
            IF( ASSOCIATED( Var % Perm ) ) k = Var % Perm(k)
            IF( k > 0 )THEN
              DO j=1, SelfDOF
                WRITE( FUnit,'(ES16.7E3)',ADVANCE='NO' ) &
                    Var % Values(Var % DOFs*(k-1)+j+Offset)
              END DO
              WRITE( FUnit, * ) 
            ELSE
              WRITE( FUnit, '(9F4.1)' ) (/ (0.0, k=1,SelfDOF) /)
            END IF
          END DO
          
          WRITE( FUnit, '(A)' ) 'attribute "dep" string "positions"'
          WRITE( FUnit, * ) 
          
          DO i = nTime + 1, 2*nTime
            WRITE( FUnit,'("object ",I0," class field")' ) i
            WRITE( FUnit,'(A,I0)' ) 'component "data" value ', i - nTime
            WRITE( FUnit,'(A,A,A)' ) 'component "positions" value file "',&
                TRIM( MeshFile ), '",1'
            WRITE( FUnit,'(A,A,A)' ) 'component "connections" value file "',&
                TRIM( MeshFile ), '",2'
            WRITE( FUnit, '(A,A,A)' ) 'attribute "name" string "', VarName,'"'
            WRITE( FUnit, * ) 
          END DO
          
          WRITE( FUnit,'(A)' ) 'object "'//VarName//'series'//'" class series'
          DO i = 1, nTime
            WRITE( FUnit, '("member ",I0," value ",I0," position ",I0)' ) &
                i-1, i+nTime, i
          END DO
          WRITE( FUnit, '("end")' ) 
          
          CLOSE( FUnit )
          
        END SUBROUTINE WriteVariable
        
    
        ! WriteDisplacement is like WriteVariable, but specialized for
        ! Displacements; displacements need special treatment.
        
        SUBROUTINE WriteDisplacement( Var, Model, nTime, MasterUnit, Prefix )
          TYPE(Variable_t), INTENT(IN) :: Var
          TYPE(Model_t), INTENT(IN) :: Model
          INTEGER, INTENT(IN) :: nTime, MasterUnit
          CHARACTER(*), INTENT(IN) :: Prefix
          INTEGER :: i, j, k
          INTEGER :: FUnit
          CHARACTER(MAX_NAME_LEN) :: FName, RelativeFName, MeshFile
          TYPE(Variable_t), POINTER :: Var1
          
          FUnit = MasterUnit + 1
          FName = Prefix // "Displacement.dx"
          
          i = INDEX(FName, '/', BACK=.TRUE.)
          RelativeFName = FName(i+1:)
          
          i = INDEX(Prefix, '/', BACK=.TRUE.)
          MeshFile = Prefix(i+1:) // "Mesh.dx"
          
          IF( nTime == 1 ) THEN
            WRITE( MasterUnit,'(A)' ) &
                'member "Displacement-vector" value file "' &
                // TRIM(RelativeFName) // '", "Displacementseries"'
            
            OPEN( FUnit, FILE=FName, STATUS="unknown" )
          ELSE
            OPEN( FUnit, FILE=FName, STATUS="old", POSITION="append" )
            DO i = 1, 2+(nTime-1)*7
              BACKSPACE FUnit   ! Yuck!
            END DO
          END IF
          
          WRITE( FUnit,'("object ",I0," class array type double rank 1 shape ",  &
              &I0," items ",I0," data follows")') nTime, Var % DOFs, &
              Model % NumberOfNodes
          DO i = 1, Model % NumberOfNodes
            k = i
            IF( ASSOCIATED( Var % Perm ) ) k = Var % Perm(k)
            IF( k > 0 )THEN
              DO j=1, Var % DOFs
                WRITE( FUnit,'(ES16.7E3)',ADVANCE='NO' ) &
                    Var % Values(Var % DOFs*(k-1)+j)
              END DO
              WRITE( FUnit, * ) 
            ELSE
              Var1 => Model % Variables
              DO WHILE( ASSOCIATED( Var1 ) )
                IF ( TRIM( Var1 % Name ) == 'mesh update' ) EXIT
                Var1 => Var1 % Next
              END DO
              IF( ASSOCIATED( Var1 ) )THEN
                k = i
                IF( ASSOCIATED(Var1 % Perm ) ) k = Var1 % Perm(k)
                IF( k > 0 )THEN
                  DO j=1,Var1 % DOFs
                    WRITE( FUnit,'(ES16.7E3)',ADVANCE='NO' ) &
                        Var1 % Values(Var1 % DOFs*(k-1)+j)
                  END DO
                  WRITE( FUnit, * ) 
                ELSE
                  WRITE( FUnit, '(9F4.1)' ) (/ (0.0, k=1,Var % DOFs) /)
                END IF
              ELSE
                WRITE( FUnit, '(9F4.1)' ) (/ (0.0, k=1,Var % DOFs) /)
              END IF
            END IF
          END DO
          
          WRITE( FUnit, '(A)' ) 'attribute "dep" string "positions"'
          WRITE( FUnit, * ) 
          
          DO i = nTime + 1, 2*nTime
            WRITE( FUnit,'("object ",I0," class field")' ) i
            WRITE( FUnit,'(A,I0)' ) 'component "data" value ', i - nTime
            WRITE( FUnit,'(A,A,A)' ) 'component "positions" value file "',&
                TRIM( MeshFile ), '",1'
            WRITE( FUnit,'(A,A,A)' ) 'component "connections" value file "',&
                TRIM( MeshFile ), '",2'
            WRITE( FUnit,'(A,A,A)' ) 'attribute "name" string "Displacement"'
            WRITE( FUnit,* ) 
          END DO
          
          WRITE( FUnit,'(A)' ) 'object "Displacementseries" class series'
          DO i = 1, nTime
            WRITE( FUnit,'("member ",I0," value ",I0," position ",I0)' ) &
                i-1, i+nTime, i
          END DO
          WRITE( FUnit, '("end")' ) 
          
          CLOSE( FUnit )
          
        END SUBROUTINE WriteDisplacement
        
        
        SUBROUTINE WriteGrid ( PRefix, Model, SubtractDisp )
          CHARACTER(*), INTENT(IN) :: Prefix
          TYPE(Model_t), INTENT(IN) :: Model
          LOGICAL, INTENT(IN) :: SubtractDisp ! Subtract Displacement from Coords.
          TYPE(Variable_t), POINTER :: Displacement, MeshUpdate, Var1, Var2
          INTEGER :: i, k, l, nElem, nVertix, dim
          CHARACTER(MAX_NAME_LEN) :: FName
          INTEGER, PARAMETER :: FUnit = 58
          INTEGER :: NodeIndex(MAX_PART_ELEM, MAX_VERTIX)
          REAL(KIND=dp) :: Coord(3)
          
          FName = Prefix // "Mesh.dx"
          OPEN( UNIT=FUnit, FILE=FName, STATUS="unknown", ACTION="write" )
          
          WRITE ( FUnit,'("# ElmerSolver output; started at ",A)' ) &
              TRIM( FormatDate() )
          !
          ! Coordinates: 
          !
          WRITE ( FUnit, '("# Node Coordinates")' )
          WRITE ( FUnit, '("object 1 class array type float rank 1 shape 3 &
              &items ",I0," data follows")' ) Model % NumberOfNodes
          
          ! First, look for displacements
          dim = Model % Mesh % MeshDim

          Displacement => NULL()
          MeshUpdate => NULL()
          IF( SubtractDisp )THEN
            Var1 => Model % Variables
            DO WHILE( ASSOCIATED( Var1 ) )
              IF( .NOT.Var1 % Output )THEN
                Var1 => Var1 % Next
                CYCLE
              END IF
              
              SELECT CASE( Var1 % Name )
              CASE( 'mesh update' )
                Var2 => Model % Variables
                DO WHILE( ASSOCIATED( Var2 ) )
                  IF ( TRIM( Var2 % Name ) == 'displacement' ) EXIT
                  Var2 => Var2 % Next
                END DO
                
                IF( .NOT. ASSOCIATED( Var2 ) )THEN
                  Displacement => Var1
                ELSE
                  MeshUpdate   => Var1
                END IF
                
              CASE( 'displacement' )
                Displacement => Var1
              END SELECT
              
              Var1 => Var1 % Next
            END DO
          END IF
          
          DO i = 1, Model % NumberOfNodes

            Coord(1) =  Model % Nodes % x(i)
            Coord(2) =  Model % Nodes % y(i)
            Coord(3) =  Model % Nodes % z(i)

            IF( SubtractDisp ) THEN
              k = 0
              IF( ASSOCIATED( Displacement ) ) k = Displacement % Perm(i)
              l = 0
              IF ( ASSOCIATED( MeshUpdate ) ) l = MeshUpdate % Perm(i)
              
              IF( k > 0 ) THEN
                k = Displacement % DOFs * (k-1)
                Coord(1) = Displacement % Values(k+1)
                IF( Displacement % DOFs >= 2) Coord(2) = Displacement % Values(k+2)
                IF( Displacement % DOFs == 2) Coord(3) = Displacement % Values(k+3)
              END IF
              IF( k == 0 .AND. l > 0 ) THEN
                l = MeshUpdate % DOFs * (l-1)
                Coord(1) = MeshUpdate % Values(l+1)
                IF( MeshUpdate % DOFs >= 2) Coord(2) = MeshUpdate % Values(l+2)
                IF( MeshUpdate % DOFs == 2) Coord(3) = MeshUpdate % Values(l+3)                
              END IF
            END IF

            IF( dim == 3 ) THEN
              WRITE( FUnit,'(3ES16.7E3)' ) Coord
            ELSE
              WRITE( FUnit,'(2ES16.7E3,A)' ) Coord(1:2),' 0.0'
            END IF
          END DO
          
          WRITE ( FUnit, * )
          
          ! Elements.
          !
          ! The only DX elements we use are triangles (for 2D) and tetrahedra
          ! (3D).  All other element types are translated into compounds of either
          ! of these.
          !
          ! Note that, at the moment, either all elements have to be 2D, or they
          ! all have to be 3D; therefore, boundary elements are not written.
          
          WRITE ( FUnit, '("# Element definitions")' )
          
          CALL GetNElem( Model, nElem, nVertix )
          WRITE( FUnit,'("object 2 class array type int rank 1 shape ", &
              &I0, " items ", I0, " data follows")' ) nVertix, nElem
          
          DO i = 1, Model % NumberOfBulkElements
            CALL TranslateElem( Model % Elements(i), NodeIndex, nElem  )
            DO k = 1, nElem
              WRITE( FUnit, '(4(" ",I0))') NodeIndex(k,1:nVertix)
            END DO
          END DO
          
          IF ( nVertix == 3 ) THEN
            WRITE( FUnit,'(A)') 'attribute "element type" string "triangles"'
          ELSE
            WRITE( FUnit,'(A)') 'attribute "element type" string "tetrahedra"'
          END IF
          WRITE( FUnit, '("end")' )
          
          
        CONTAINS
          
          ! Here's a bunch of routines to deal with ELMER -> DX element
          ! translations.
          
          ! nElem is the number of DX elements, and nVertix is the numbe of
          ! vertices per element (either 3 for triangle or 4 for tetrahedra).
          
          SUBROUTINE GetNElem( Model, nElem, nVertix )
            TYPE(Model_t), INTENT(IN) :: Model
            INTEGER, INTENT(OUT) :: nElem, nVertix
            
            IF (Model % Elements(1) % TYPE % ElementCode < 500) THEN
              nVertix = 3
            ELSE
              nVertix = 4
            END IF
            
            nElem = 0
            DO i = 1, Model%NumberOfBulkElements
              SELECT CASE ( Model % Elements(i) % TYPE % ElementCode )
              CASE( 303 )
                nElem = nElem + 1
              CASE( 404 )
                nElem = nElem + 2
              CASE( 409 )
                nElem = nElem + 8
              CASE( 504 )
                nElem = nElem + 1
              CASE( 605 )
                nElem = nElem + 2
              CASE( 613 )
                nElem = nElem + 14
              CASE( 706 )
                nElem = nElem + 3
              CASE( 715 )
                nElem = nElem + 21
              CASE( 808 )
                nElem = nElem + 4
              CASE DEFAULT
                WRITE (0, *) 'Sorry! Element type ',                    &
                    Model % Elements(i) % TYPE % ElementCode, &
                    ' not supported! '
                STOP
              END SELECT
            END DO
          END SUBROUTINE GetNElem
          
          
          ! Translate an ELMER element to a (compound of) DX element(s).
          
          SUBROUTINE TranslateElem( Elem, NodeIndex, nElem )
            TYPE(Element_t), INTENT(IN) :: Elem
            INTEGER, INTENT(OUT) :: NodeIndex(MAX_PART_ELEM, MAX_VERTIX)
            INTEGER, INTENT(OUT) :: nElem
            
            SELECT CASE ( Elem % TYPE % ElementCode )
            CASE( 303 ) ! Triangle
              nElem = 1
              NodeIndex(1,:3) = Elem % NodeIndexes(1:3)
            CASE( 404 ) ! Quad; translate into 2 triangles
              nElem = 2
              NodeIndex(1,:3) = Elem % NodeIndexes(1:3)
              NodeIndex(2,:3) = Elem % NodeIndexes((/ 1,3,4 /))
            CASE( 409 )
              nElem = 8
              NodeIndex(1,:3) = Elem % NodeIndexes((/ 1,5,9 /))
              NodeIndex(2,:3) = Elem % NodeIndexes((/ 5,2,9 /))
              NodeIndex(3,:3) = Elem % NodeIndexes((/ 2,6,9 /))
              NodeIndex(4,:3) = Elem % NodeIndexes((/ 6,3,9 /))
              NodeIndex(5,:3) = Elem % NodeIndexes((/ 3,7,9 /))
              NodeIndex(6,:3) = Elem % NodeIndexes((/ 7,4,9 /))
              NodeIndex(7,:3) = Elem % NodeIndexes((/ 4,8,9 /))
              NodeIndex(8,:3) = Elem % NodeIndexes((/ 8,1,9 /))
            CASE( 504 ) ! Tetrahedron
              nElem = 1
              NodeIndex(1,:) = Elem % NodeIndexes(1:4)
            CASE( 605 ) ! Pyramid
              nElem = 2
              NodeIndex(1,:) = Elem % NodeIndexes((/ 3,5,4,1 /))
              NodeIndex(2,:) = Elem % NodeIndexes((/ 3,5,2,1 /))
            CASE( 613 )
              nElem = 14
              NodeIndex(1,:) = Elem % NodeIndexes((/ 7,8,3,12 /))
              NodeIndex(2,:) = Elem % NodeIndexes((/ 10,11,12,5 /))
              NodeIndex(3,:) = Elem % NodeIndexes((/ 10,13,12,5 /))
              NodeIndex(4,:) = Elem % NodeIndexes((/ 9,10,13,11 /))
              NodeIndex(5,:) = Elem % NodeIndexes((/ 9,13,11,12 /))
              NodeIndex(6,:) = Elem % NodeIndexes((/ 9,10,11,1 /))
              NodeIndex(7,:) = Elem % NodeIndexes((/ 9,6,11,1 /))
              NodeIndex(8,:) = Elem % NodeIndexes((/ 9,6,11,12 /))
              NodeIndex(9,:) = Elem % NodeIndexes((/ 9,8,12,4 /))
              NodeIndex(10,:) = Elem % NodeIndexes((/ 9,13,12,4 /))
              NodeIndex(11,:) = Elem % NodeIndexes((/ 7,9,8,12 /))
              NodeIndex(12,:) = Elem % NodeIndexes((/ 7,9,6,12 /))
              NodeIndex(13,:) = Elem % NodeIndexes((/ 7,6,11,12 /))
              NodeIndex(14,:) = Elem % NodeIndexes((/ 7,6,11,2 /))
            CASE( 706 ) ! Wedge
              nElem = 3
              NodeIndex(1,:) = Elem % NodeIndexes((/ 5,4,3,1 /))
              NodeIndex(2,:) = Elem % NodeIndexes((/ 5,3,2,1 /))
              NodeIndex(3,:) = Elem % NodeIndexes((/ 5,6,4,3 /))
            CASE( 715 )
              nElem = 21
              NodeIndex(1,:) = Elem % NodeIndexes((/ 10,11,5,2 /))
              NodeIndex(2,:) = Elem % NodeIndexes((/ 12,11,6,3 /))
              NodeIndex(3,:) = Elem % NodeIndexes((/ 12,10,4,1 /))
              NodeIndex(4,:) = Elem % NodeIndexes((/ 7,8,11,2 /))
              NodeIndex(5,:) = Elem % NodeIndexes((/ 7,10,11,2 /))
              NodeIndex(6,:) = Elem % NodeIndexes((/ 13,14,11,5 /))
              NodeIndex(7,:) = Elem % NodeIndexes((/ 13,10,11,5 /))
              NodeIndex(8,:) = Elem % NodeIndexes((/ 9,10,8,11 /))
              NodeIndex(9,:) = Elem % NodeIndexes((/ 9,7,10,8 /))
              NodeIndex(10,:) = Elem % NodeIndexes((/ 9,12,10,11 /))
              NodeIndex(11,:) = Elem % NodeIndexes((/ 9,12,10,1 /))
              NodeIndex(12,:) = Elem % NodeIndexes((/ 9,7,10,1 /))
              NodeIndex(13,:) = Elem % NodeIndexes((/ 9,8,11,3 /))
              NodeIndex(14,:) = Elem % NodeIndexes((/ 9,12,11,3 /))
              NodeIndex(15,:) = Elem % NodeIndexes((/ 15,12,10,11 /))
              NodeIndex(16,:) = Elem % NodeIndexes((/ 15,12,10,4 /))
              NodeIndex(17,:) = Elem % NodeIndexes((/ 15,10,14,11 /))
              NodeIndex(18,:) = Elem % NodeIndexes((/ 15,13,10,14 /))
              NodeIndex(19,:) = Elem % NodeIndexes((/ 15,13,10,4 /))
              NodeIndex(20,:) = Elem % NodeIndexes((/ 15,14,11,6 /))
              NodeIndex(21,:) = Elem % NodeIndexes((/ 15,12,11,6 /))
            CASE( 808 ) ! Hexahedron; translate into 4 tetrahedra
              nElem = 4
              NodeIndex(1,:) = Elem % NodeIndexes((/ 1,2,4,5 /))
              NodeIndex(2,:) = Elem % NodeIndexes((/ 2,3,4,7 /))
              NodeIndex(3,:) = Elem % NodeIndexes((/ 2,7,5,6 /))
              NodeIndex(4,:) = Elem % NodeIndexes((/ 4,5,7,8 /))
            END SELECT
            
            NodeIndex = NodeIndex - 1
            
          END SUBROUTINE TranslateElem
          
        END SUBROUTINE WriteGrid
        
      END MODULE DXFile
      
      
!------------------------------------------------------------------------------
!> Module for the DX result writer.
!------------------------------------------------------------------------------
      SUBROUTINE DXOutputSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------

        USE DefUtils 
        USE DXFile
        
        IMPLICIT NONE
        TYPE(Solver_t) :: Solver
        TYPE(Model_t) :: Model
        REAL(dp) :: dt
        LOGICAL :: TransientSimulation
        
        INTEGER, SAVE :: nTime = 0
        LOGICAL :: GotIt
        CHARACTER(MAX_NAME_LEN), SAVE :: FilePrefix
        
        ! Avoid compiler warings about unused variables
        IF ( TransientSimulation ) THEN; ENDIF
          IF ( dt > 0.0 ) THEN; ENDIF
            
            IF ( nTime == 0 ) THEN
              FilePrefix = GetString( Solver % Values,'Output File Name',GotIt )
              IF ( .NOT.GotIt ) FilePrefix = "Output"
            END IF
            nTime = nTime + 1
            
            CALL WriteData( TRIM(FilePrefix), Model, nTime )
            

CONTAINS

  SUBROUTINE WriteData( Prefix, Model, nTime )
    CHARACTER(*), INTENT(IN) :: Prefix
    TYPE(Model_t) :: Model
    INTEGER, INTENT(IN) :: nTime
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Variable_t), POINTER :: Var
    INTEGER :: i, j, k
    LOGICAL :: EigAnal
    REAL(dp), POINTER :: OldValues(:)
    CHARACTER(MAX_NAME_LEN) :: Dir
    
    Mesh => Model % Mesh
      
    IF (LEN_TRIM(Mesh % Name) > 0 ) THEN
      Dir = TRIM(Mesh % Name) // "/"
    ELSE
      Dir = "./"
    END IF
      
    EigAnal = .FALSE.
    
    Solvers: DO i = 1, Model % NumberOfSolvers
      EigAnal = ListGetLogical( Model % Solvers(i) % Values, &
          "Eigen Analysis", GotIt )
      Var => Model % Solvers(i) % Variable
      IF ( EigAnal .AND. ASSOCIATED(Var % EigenValues) ) THEN
        DO j = 1, Model % Solvers(i) % NOfEigenValues
          OldValues => Var % Values
          
          IF ( Model % Solvers(i) % Matrix % COMPLEX ) THEN
            ALLOCATE( Var % Values(2*SIZE(Var%EigenVectors,2)) )
            FORALL ( k = 1:SIZE(Var % Values)/2 )
              Var%Values(2*k-1) = REAL(Var%EigenVectors(j,k))
              Var%Values(2*k) = AIMAG(Var%EigenVectors(j,k))
            END FORALL
          ELSE
            ALLOCATE( Var % Values(SIZE(Var % EigenVectors,2)) )
            Var % Values = Var % EigenVectors(j,:)
          END IF
          
          CALL WriteDXFiles( TRIM(Dir)//Prefix,Model,.FALSE.,j )
          
          DEALLOCATE( Var % Values )
          Var % Values => OldValues
        END DO
        EXIT Solvers
      END IF
    END DO Solvers
    
    IF ( .NOT.EigAnal ) THEN
      CALL WriteDXFiles( TRIM(Dir)//Prefix, Model, .TRUE., nTime )
    END IF

  END SUBROUTINE WriteData
  
END SUBROUTINE DXOutputSolver


