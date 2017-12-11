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
! *  Authors: Peter Råback
! *  Email:   Peter.Raback@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 29.09.2016
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
!> Module for reading in OpenFOAM cell centers and writing data interpolated on them.
!------------------------------------------------------------------------------
SUBROUTINE Elmer2OpenFoamWrite( Model,Solver,dt,TransientSimulation )
  
  USE DefUtils
  USE Interpolation
  USE MeshUtils
  USE ElementUtils
  USE ParticleUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
  
! local variables
!------------------------------------------------------------------------------  
  TYPE(ValueList_t), POINTER :: Params
  TYPE(Variable_t), POINTER :: Var, OFVar
  TYPE(Mesh_t), POINTER :: Mesh, OFMesh
  CHARACTER(LEN=MAX_NAME_LEN) :: VarName, FileName, DirName, BaseDir, OFfile
  INTEGER :: i, NoDir, IOStatus, PassiveCoord
  LOGICAL :: Found, Visited = .FALSE., UseProjFound, UseProjSave
  REAL(KIND=dp) :: MinF, MaxF, MeanF
  
  SAVE OFMesh, OFVar, NoDir, Visited
  

  CALL Info('Elmer2OpenFoamWrite','-----------------------------------------', Level=4 )
  CALL Info('Elmer2OpenFoamWrite','Projecting field to OpenFOAM cell centers',Level=4) 

  
  ! The variable containing the field contributions
  !------------------------------------------------------------------------
  Params => GetSolverParams()
  VarName = GetString( Params,'Target Variable',Found)
  IF(.NOT. Found ) THEN
    CALL Fatal('Elmer2OpenFoamWrite','> Target Variable < must exist for the solver!')
  END IF

  ! Save the status of the "Use Mesh Projector" keyword.
  UseProjSave = ListGetLogical( Model % Simulation,'Use Mesh Projector',UseProjFound )
  CALL Info('Elmer2OpenFOAMWrite','Enforcing mapping without projector matrix!',Level=6)
  CALL ListAddLogical( Model % Simulation,'Use Mesh Projector',.FALSE.)
    
  
  Mesh => GetMesh()
  ! Test that the variable exists in the primary mesh
  Var => VariableGet(Mesh % Variables, VarName )
  IF(.NOT. ASSOCIATED( Var ) ) THEN
    CALL Fatal('Elmer2OpenFoamWrite','Variable does not exist in Elmer mesh!')
  END IF

  ! If we visit this the second time, then destroy the structures that were saved last time.
  IF( Visited ) THEN
    CALL ReleaseMesh( OFMesh )
    DEALLOCATE( OFMesh ) 
  END IF
  OFMesh => AllocateMesh()
  
  
  ! If the Elmer mesh has different dimension we may make a simple
  ! dimensional reduction for the OpenFOAM mesh.
  !-------------------------------------------------------------------------
  PassiveCoord = ListGetInteger( Solver % Values,'Passive OpenFOAM Coordinate',Found )
  IF( .NOT. Found .AND. Mesh % MeshDim < 3 ) THEN
    CALL Warn('Elmer2OpenFOAM','Dimension of Elmer mesh is reduced, and OpenFOAM not?!')
  END IF
  
  
  ! This is just for helping to write Elmer cell centers for testing purposes
  FileName = GetString( Params,'Elmer Center Filename',Found)
  IF(Found ) THEN
    CALL WriteMeshCenters( Mesh, Filename )
    RETURN
  END IF
  
  BaseDir = GetString( Params,'OpenFOAM Directory',Found)
  IF( Found ) THEN
    CALL Info('Elmer2OpenFoamWrite','Using given > OpenFOAM Directory < : '//TRIM(BaseDir),Level=6)
  ELSE
    CALL Fatal('Elmer2OpenFoamWrite','> OpenFOAM Directory < must exist for the solver!')
  END IF

  
  ! If the blocks do not exist find them 
  ! When they are stored as keywords the user may give them also manually
  IF( .NOT. ListCheckPresent(Params,'OpenFOAM Mesh 1') ) THEN
    CALL OpenFOAMBlocks()
  END IF
  NoDir = NINT( ParallelReduction(1.0_dp * NoDir ) )
  CALL Info('Elmer2OpenFOAMWrite','Number of active OpenFOAM blocks: '//TRIM(I2S(NoDir)),Level=5)

  
  DO i = 1, NoDir
    
    IF( ParEnv % MyPe == 0 ) THEN
      IF( NoDir > 1 ) THEN
        CALL Info('Elmer2OpenFOAMWrite','Treating OpenFOAM block: '//TRIM(I2S(i)),Level=5)
      END IF
      DirName = ListGetString(Params,'OpenFOAM Mesh '//TRIM(I2S(i)),Found)
      IF(.NOT. Found ) CALL Fatal('Elmer2OpenFoamWrite','Could not find keyword: '//TRIM(DirName))
      
      FileName = TRIM(DirName)//'C'    
      CALL Info('Elmer2OpenFoamWrite','Projecting to OpenFOAM nodes in file: '//TRIM(FileName),Level=5)
    END IF

    
    CALL CreateFOAMMesh(FileName,OFMesh)

    CurrentModel % Mesh => OFMesh

    IF( PassiveCoord == 1 ) THEN
      OFMesh % Nodes % x = 0.0_dp
    ELSE IF( PassiveCoord == 2 ) THEN
      OFMesh % Nodes % y = 0.0_dp
    ELSE IF( PassiveCoord == 3 ) THEN
      OFMesh % Nodes % z = 0.0_dp
    END IF

    CALL Info('Elmer2OpenFoamWrite','Mapping data to the temporal mesh using libary routines',Level=6)
    
    OFVar => VariableGet(OFMesh % Variables, VarName )
    
    ! Put the solver variable so that we can study norms etc.
    Solver % Variable => OFVar
    
    IF( ParEnv % MyPe == 0 ) THEN
      MinF = MINVAL( OFVar % Values )
      MaxF = MAXVAL( OFVar % Values )
      MeanF = SUM( OFVar % Values ) / SIZE( OFVar % Values ) 

      WRITE( Message,'(A,ES12.5)') 'Minimum field value: ',MinF
      CALL Info('Elmer2OpenFoamWrite',Message,Level=6)
      WRITE( Message,'(A,ES12.5)') 'Maximum field value: ',MaxF
      CALL Info('Elmer2OpenFoamWrite',Message,Level=6)
      WRITE( Message,'(A,ES12.5)') 'Average field value: ',MeanF
      CALL Info('Elmer2OpenFoamWrite',Message,Level=6)
      
      OFfile = GetString( Params,'OpenFOAM file',Found)
      IF(.NOT. Found ) Offile = 'fieldSolidHS.dat'
      FileName = TRIM(DirName)//TRIM(Offile)
      CALL Info('Elmer2OpenFoamWrite','Writing projected field to file: '//TRIM(FileName),Level=5)
      CALL WriteOFField(Filename, OFMesh, OFVar)
    END IF

    ! We can only have one OpenFOAM mesh at a time, hence release the structures if we have a second mesh.
    IF( i < NoDir ) CALL ReleaseMesh( OFMesh )
  END DO
    
  ! Restore the pointer to the initial Elmer mesh
  CurrentModel % Mesh => Mesh

  ! Restore the status of the "Use Mesh Projector" keyword.
  IF( UseProjFound ) THEN
    CALL ListAddLogical( Model % Simulation,'Use Mesh Projector',UseProjSave )
  ELSE
    CALL ListRemove( Model % Simulation,'Use Mesh Projector')
  END IF

  Visited = .TRUE.
    
  CALL Info('Elmer2OpenFoamWrite','All done', Level=4 )
  CALL Info('Elmer2OpenFoamWrite','-----------------------------------------', Level=4 )  
  
  
CONTAINS 
  
  SUBROUTINE OpenFOAMBlocks( )
    
    CHARACTER(LEN=MAX_NAME_LEN) :: DirCommand
    LOGICAL :: FileExists
    INTEGER, PARAMETER :: InFileUnit = 28
    
    
    NoDir = 0
    IF( ParEnv % MyPe /= 0 ) RETURN

#ifdef __INTEL_COMPILER
    ! Fortran standard states that inquiry for a file returns true if the queried entity is a file
    INQUIRE( Directory = TRIM(BaseDir), Exist = FileExists )
#else
    INQUIRE( File = TRIM(BaseDir), Exist = FileExists )
#endif
    IF(.NOT. FileExists ) THEN
      CALL Fatal('Elmer2OpenFoamWrite','OpenFOAM directory does not exist: '//TRIM(BaseDir))
    END IF

    DirName = TRIM(BaseDir)//'/0/'
#ifdef __INTEL_COMPILER
    INQUIRE( Directory = DirName, Exist = FileExists )
#else
    INQUIRE( File = DirName, Exist = FileExists )
#endif
    IF(.NOT. FileExists ) THEN
      CALL Fatal('Elmer2OpenFoamWrite','OpenFOAM mesh does not exist: '//TRIM(DirName))
    END IF
    
    FileName = TRIM(DirName)//'C'
    CALL Info('Elmer2OpenFoamWrite','Inquire file: '//TRIM(FileName),Level=12)
    INQUIRE( File = FileName, Exist = FileExists )
   
    IF( FileExists ) THEN
      CALL Info('Elmer2OpenFoamWrite','Using OpenFOAM centers in: '//TRIM(FileName),Level=10)
      CALL ListAddString( Params, 'OpenFOAM Mesh 1', DirName, .FALSE.)
      NoDir = 1
      RETURN
    END IF

    DirCommand = 'ls -d '//TRIM(DirName)//'*/ > OpenFOAMBlocks.txt' 
    CALL Info('Elmer2OpenFoamWrite','Performing command: '//TRIM(DirCommand),Level=12)
    CALL SystemCommand( DirCommand )

    OPEN(InFileUnit,File='OpenFOAMBlocks.txt',IOStat=IOstatus)
    IF(IOStatus /= 0 ) THEN
      CALL Fatal('Elmer2OpenFoamWrite','Could not open file: OpenFOAMBlocks.txt')
    END IF
     
    DO
      READ(InFileUnit,'(A)',IOStat = IOStatus) DirName
      IF( IOStatus /= 0 ) EXIT
      FileName = TRIM(DirName)//'C'
      CALL Info('Elmer2OpenFoamWrite','Inquire file: '//TRIM(FileName),Level=12)
      INQUIRE( File = FileName, Exist = FileExists )
      IF( FileExists ) THEN
        NoDir = NoDir + 1
        CALL Info('Elmer2OpenFoamWrite','Using OpenFOAM centers in: '//TRIM(FileName),Level=10)
        CALL ListAddString( Params, 'OpenFOAM Mesh '//TRIM(I2S(NoDir)), DirName, .FALSE.)
      ELSE
        CALL Info('Elmer2OpenFoamWrite','No OpenFOAM center in: '//TRIM(DirName),Level=12)
      END IF
    END DO
    CLOSE(InFileUnit)
    
    IF( NoDir == 0 ) THEN
      CALL Fatal('Elmer2OpenFoamWrite','No OpenFOAM mesh blocks found!')
    ELSE
      CALL Info('Elmer2OpenFoamWrite','Number of OpenFOAM blocks: '//TRIM(I2S(NoDir)),Level=10)
    END IF
    
  END SUBROUTINE OpenFOAMBlocks
    

  SUBROUTINE WriteMeshCenters( Mesh, Filename )

    CHARACTER(LEN=MAX_NAME_LEN) :: FileName
    TYPE(Mesh_t), POINTER :: Mesh
    
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: x,y,z
    INTEGER, PARAMETER :: OutFileUnit = 29
    INTEGER :: i,n,IOstatus


    OPEN(OutFileUnit,FILE = Filename, IOSTAT=IOStatus)
    IF( IOStatus /= 0 ) THEN
      CALL Fatal('Elmer2OpenFoamWrite','Could not open file for writing: '//TRIM(FileName))
    END IF
    
    CALL Info('Elmer2OpenFoamWrite','Writing Elmer element centers to file: '//TRIM(FileName),Level=6)

    DO i=1, Mesh % NumberOfBulkElements
      Element => Mesh % Elements(i)
      n = Element % TYPE % NumberOfNodes
      x = SUM( Mesh % Nodes % x( Element % NodeIndexes ) ) / n
      y = SUM( Mesh % Nodes % y( Element % NodeIndexes ) ) / n
      z = SUM( Mesh % Nodes % z( Element % NodeIndexes ) ) / n
      WRITE( OutFileUnit, * ) '(',x,y,z,')'
    END DO

    CLOSE( OutFileUnit ) 
    
  END SUBROUTINE WriteMeshCenters


  
  !------------------------------------------------------------------------
  !> Open file in OpenFOAM format and read the cell centers from there.
  !-------------------------------------------------------------------------
  SUBROUTINE CreateFOAMMesh( Filename, Mesh ) 
    
    CHARACTER(LEN=MAX_NAME_LEN) :: FileName
    TYPE(Mesh_t), POINTER :: Mesh
   
    INTEGER :: line,i,j,k,n
    REAL(KIND=dp) :: x,y,z
    INTEGER :: NumberOfNodes, IOStatus
    INTEGER, PARAMETER :: InFileUnit = 28
    CHARACTER(LEN=:), ALLOCATABLE :: ReadStr

    ALLOCATE( Mesh % Nodes )
    ALLOCATE( Mesh % Variables )
    Mesh % NumberOfBulkElements = 0
    Mesh % NumberOfBoundaryElements = 0
    
    ! Partition zero does all the work!
    IF( ParEnv % MyPe /= 0 ) THEN
      Mesh % NumberOfNodes = 0
      GOTO 100      
    END IF
    
    
    ALLOCATE(CHARACTER(MAX_STRING_LEN)::ReadStr)
                
    OPEN(InFileUnit,FILE = Filename, STATUS='old', IOSTAT=IOstatus)
    IF( IOStatus /= 0 ) THEN
      CALL Fatal('Elmer2OpenFoamWrite','Could not open file for reading: '//TRIM(FileName))
    END IF
    
    CALL Info('Elmer2OpenFoamWrite','Reading data points from file: '//TRIM(FileName),Level=6)
    
    j = 0
    DO Line = 1, 100
      READ( InFileUnit,'(A)',IOSTAT=IOStatus ) ReadStr
      IF( IOStatus /= 0 ) THEN
        CALL Warn('Elmer2OpenFoamWrite','End of file after '//TRIM(I2S(Line))//' lines')
        EXIT
      END IF

      j =  INDEX( ReadStr,'internalField',.TRUE.) 
      IF( j > 0 ) EXIT
    END DO

    IF( j == 0 ) THEN
      CALL Fatal('Elmer2OpenFoamWrite','Could not find > internalField < in header!')
    ELSE
      CALL Info('Elmer2OpenFoamWrite','internalField found at line: '//TRIM(I2S(Line)),Level=7)    
    END IF
      
    READ(InFileUnit,*,IOSTAT=IOStatus) NumberOfNodes    
    IF( IOStatus /= 0 ) THEN
      CALL Fatal('Elmer2OpenFoamWrite','Could not read number of nodes!')
    END IF
    CALL Info('Elmer2OpenFoamWrite','Number of OpenFOAM nodes: '&
        //TRIM(I2S(NumberOfNodes)))

    i = ListGetInteger(Params,'Number of cells',Found)
    IF( i > 0 .AND. i < NumberOfNodes ) THEN
      NumberOfNodes = i
      CALL Info('Elmer2OpenFoamWrite','Limiting number of OpenFOAM nodes: '&
          //TRIM(I2S(NumberOfNodes)))
    END IF


    n = NumberOfNodes    
    ALLOCATE( Mesh % Nodes % x(n), &          
        Mesh % Nodes % y(n), &
        Mesh % Nodes % z(n) )

    Mesh % Nodes % x(1:n) = 0.0_dp
    Mesh % Nodes % y(1:n) = 0.0_dp
    Mesh % Nodes % z(1:n) = 0.0_dp
    
    Mesh % NumberOfNodes = n
         

    ! This is just empty left paranthesis
    READ( InFileUnit,'(A)',IOSTAT=IOStatus ) ReadStr
    !PRINT *,'EmptyLine:',TRIM(ReadStr)
   
    DO i=1,n
      READ( InFileUnit,'(A)',IOSTAT=IOStatus ) ReadStr
      IF( IOStatus /= 0 ) THEN
        CALL Fatal('Elmer2OpenFoamWrite','Could not read coordinate line: '//TRIM(I2S(i)))
      END IF
      
      j =  INDEX( ReadStr,'(',.TRUE.) 
      IF( j == 0 ) THEN
        CALL Fatal('Elmer2OpenFoamWrite',&
            'Expecting a paranthesis at the start of OpenFOAM line: '//TRIM(I2S(i)))
      END IF
      k =  INDEX( ReadStr,')',.TRUE.) 
      IF( k == 0 ) THEN
        CALL Fatal('Elmer2OpenFoamWrite',&
            'Expecting a paranthesis at the end of OpenFOAM line: '//TRIM(I2S(i)))
      END IF
      
      READ( ReadStr(j+1:k-1),*,IOSTAT=IOStatus ) x,y,z
      IF( IOStatus /= 0 ) THEN
        CALL Fatal('Elmer2OpenFoamWrite','Could not read coordinate values: '//TRIM(I2S(i)))
      END IF
      Mesh % Nodes % x(i) = x
      Mesh % Nodes % y(i) = y
      Mesh % Nodes % z(i) = z
    END DO
    CLOSE( InFileUnit ) 

    
    CALL Info('Elmer2OpenFoamWrite','Creating coordinates for temporal mesh',Level=7)

100 CALL VariableAdd( Mesh % Variables, Mesh, Solver, &
        'Coordinate 1',1, Mesh % Nodes % x )
    
    CALL VariableAdd( Mesh % Variables, Mesh, Solver, &
        'Coordinate 2',1, Mesh % Nodes % y )
    
    CALL VariableAdd( Mesh % Variables, Mesh, Solver, &
        'Coordinate 3',1,Mesh % Nodes % z )
    
    CALL Info('Elmer2OpenFoamWrite','Created temporal OpenFOAM mesh just for nodes',Level=8)
        
  END SUBROUTINE CreateFOAMMesh


  !------------------------------------------------------------------------
  !> Write the data to the OpenFOAM cell centers. This assumes same order
  !> for the cells as was used in reading in the temporal mesh.
  !-------------------------------------------------------------------------
  SUBROUTINE WriteOFField( Filename, Mesh, Var )
    CHARACTER(LEN=MAX_NAME_LEN) :: FileName
    TYPE(Mesh_t) :: Mesh
    TYPE(Variable_t) :: Var
    INTEGER :: i,j,IOStatus
    
    INTEGER, PARAMETER :: OutFileUnit = 29
    CHARACTER :: NL


    NL = NEW_LINE('A')
    
    CALL Info('Elmer2OpenFoamWrite','Writing a file with field value for OpenFOAM: '//TRIM(FileName),Level=6)
    
    OPEN(OutFileUnit,FILE = Filename, IOSTAT=IOStatus)
    
    WRITE(OutFileUnit,'(A)') &        
        '/*--------------------------------*- C++ -*----------------------------------*\ '//NL//&
        '| =========                 |                                                 | '//NL//&
        '| \\      /  F ield         | OpenFOAM: The OPEN Source CFD Toolbox           | '//NL//&
        '|  \\    /   O peration     | Version:  dev                                   | '//NL//&
        '|   \\  /    A nd           | Web:      www.OpenFOAM.org                      | '//NL//&
        '|    \\/     M anipulation  |                                                 | '//NL//&
        '|                                                                             | '//NL//&
        '|  This file written by Elmer2OpenFOAM coupling module in Elmer               | '//NL//&
        '\*---------------------------------------------------------------------------*/ '
    WRITE(OutFileUnit,'(A)') &        
        'FoamFile                                '//NL//&
        '{                                       '//NL//&
        '    version     2.0;                    '//NL//& 
        '    format      ascii;                  '//NL//&
        '    class       volScalarField;         '//NL//&   
        '    object      fieldSolidHS;           '//NL//&
        '}'


      WRITE(OutFileUnit,'(A)') &        
        'dimensions      [1 -1 -3 0 0 0 0];      '        
      
      WRITE(OutFileUnit,'(A)') &        
          'boundaryField                           '//NL//&
          '{                                       '//NL//&    
          '    ".*"                                '//NL//&
          '    {                                   '//NL//& 
          '        TYPE zeroGradient;              '//NL//&
          '    }                                   '//NL//&
          '}'
          
      WRITE(OutFileUnit,'(A)') &        
          'internalField   nonuniform List<scalar> '//NL//&
          TRIM(I2S(Mesh % NumberOfNodes))           //NL//&
          '(' 
    
      DO i=1,Mesh % NumberOfNodes
        IF( ASSOCIATED( Var % Perm ) ) THEN
          j = Var % Perm(i)
        ELSE
          j = i
        END IF
        IF( j < 0 .OR. j > Mesh % NumberOfNodes ) THEN
          CALL Fatal('Elmer2OpenFoamWrite','We have a troubling ENTRY: '//TRIM(I2S(i))//','//TRIM(I2S(j)))
        END IF
        
        WRITE( OutFileUnit,*) Var % Values(j)
      END DO

      WRITE(OutFileUnit,'(A)') ');'             
      CLOSE( OutFileUnit ) 
      CALL Info('Elmer2OpenFoamWrite','Created a file with OpenFOAM field values',Level=8)

    END SUBROUTINE WriteOFField


END SUBROUTINE Elmer2OpenFoamWrite
