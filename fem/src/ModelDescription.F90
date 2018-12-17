!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation; either
! * version 2.1 of the License, or (at your option) any later version.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! * 
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library (in file ../LGPL-2.1); if not, write 
! * to the Free Software Foundation, Inc., 51 Franklin Street, 
! * Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 01 Oct 1996
! *
! *****************************************************************************/

!> \ingroup ElmerLib
!> \{

!-----------------------------------------------------------------------------
!>  Module that defined the Model: reads in the command file, mesh and results etc.
!-----------------------------------------------------------------------------
#include "../config.h"

MODULE ModelDescription

#ifdef USE_ISO_C_BINDINGS
    USE LoadMod
#endif
    USE MeshUtils
    USE ElementDescription
    USE BinIO
    USE Messages
 
    IMPLICIT NONE
#ifndef USE_ISO_C_BINDINGS
    INTERFACE
      FUNCTION LoadFunction( Quiet,Abort_not_found,Libname,Procname ) RESULT(Proc)
        USE, INTRINSIC :: ISO_C_BINDING
        USE Types
        CHARACTER(C_CHAR) :: Libname(*),Procname(*)
        INTEGER(KIND=AddrInt) :: Proc
        INTEGER(C_INT) :: Quiet, Abort_not_found
      END FUNCTION LoadFunction
    END INTERFACE
#endif

    CHARACTER(LEN=1024) :: IncludePath = ' ', OutputPath = ' ', SimulationId = ' '

    INTEGER, PARAMETER :: PosUnit = 32, OutputUnit = 31, RestartUnit = 30,&
                          PostFileUnit = 29, InFileUnit = 28

    INTEGER, PARAMETER, PRIVATE :: MAX_OUTPUT_VARS = 1000

CONTAINS

!------------------------------------------------------------------------------
!> Loads a dynamic object (e.g. solver or user defined function) and returns
!> its address in order to be able to call it. 
!------------------------------------------------------------------------------
  FUNCTION GetProcAddr( str, Quiet, Abort ) RESULT( Proc )
!------------------------------------------------------------------------------
    CHARACTER(LEN=*) :: str
    LOGICAL, OPTIONAL :: Quiet, Abort

    INTEGER(KIND=AddrInt) :: Proc
    INTEGER   :: i,j,slen,q,a
    CHARACTER :: Libname(MAX_NAME_LEN),Procname(MAX_NAME_LEN)
!------------------------------------------------------------------------------

    DO slen=LEN(str),1,-1
      IF ( str(slen:slen) /= ' ' ) EXIT
    END DO

    i = 1
    DO WHILE( i <= slen )
      IF ( str(i:i) == ' ' ) EXIT
      Libname(i) = str(i:i)
      i = i + 1
    END DO
    Libname(i) = CHAR(0)

    DO WHILE( i <= slen )
       IF ( str(i:i) /= ' ' ) EXIT
       i = i + 1
    END DO

    j = 1
    DO WHILE( i <= slen )
      IF (  str(i:i) == ' ' ) EXIT
      Procname(j) = str(i:i)
      i = i + 1
      j = j + 1
    END DO
    ProcName(j) = CHAR(0)

    q = 0
    IF ( OutputPE < 0 ) THEN
      q=1
    ELSE IF( PRESENT(Quiet) ) THEN
      IF ( Quiet ) q = 1
    ELSE
      IF ( .NOT. OutputLevelMask(4) ) q = 1
    END IF

    a = 1
    IF ( PRESENT(abort) ) THEN
       IF ( .NOT. abort ) a=0
    END IF

    Proc = LoadFunction( q,a,Libname,Procname )
  END FUNCTION GetProcAddr
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
  SUBROUTINE LoadIncludeFile( Model,InFileUnit,FileName,MeshDir,MeshName,ScanOnly )
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     INTEGER :: InFileUnit
     LOGICAL :: ScanOnly
     CHARACTER(LEN=*) :: FileName
     CHARACTER(LEN=*) :: MeshDir,MeshName
!------------------------------------------------------------------------------
     CHARACTER(LEN=MAX_STRING_LEN) :: FName
 
     INTEGER :: k,k0,k1,l
!------------------------------------------------------------------------------

     CALL Info('LoadIncludeFile','Loading include file: '//TRIM(FileName),Level=8)
     
     IF ( .NOT. FileNameQualified(FileName) ) THEN
       k0 = 1
       k1 = INDEX( IncludePath, ';' )
       DO WHILE( k1 >= k0 )
         DO k = k1-1,k0,-1
           IF ( IncludePath(k:k) /= ' ' ) EXIT
         END DO 

         IF ( k >= k0 ) THEN
           WRITE( FName, '(a,a,a)' ) IncludePath(k0:k), '/', &
              TRIM( FileName )
           OPEN( InFileUnit, FILE=TRIM(FName), STATUS='OLD',ERR=10 )
           CALL LoadInputFile( Model, InFileUnit, FName, &
                 MeshDir, MeshName, .FALSE., ScanOnly )
           CLOSE( InFileUnit )
           RETURN
         END IF

10       CONTINUE

         k0 = k1+1
         k1 = INDEX( IncludePath(k0:), ';' ) + k0 - 1
       END DO

       IF ( LEN_TRIM(IncludePath) > 0 ) THEN
         WRITE( FName, '(a,a,a)' ) TRIM(IncludePath(k0:)), '/', &
            TRIM( FileName )

         OPEN( InFileUnit, FILE=TRIM(FName), STATUS='OLD',ERR=20 )
         CALL LoadInputFile( Model, InFileUnit, FName, &
                MeshDir, MeshName, .FALSE., ScanOnly )
         CLOSE( InFileUnit )
         RETURN
       END IF

20     CONTINUE

       OPEN( InFileUnit, FILE=TRIM(FileName), STATUS='OLD' )
       CALL LoadInputFile( Model, InFileUnit, FileName, &
              MeshDir, MeshName, .FALSE., ScanOnly )
       CLOSE( InFileUnit )

     ELSE
       OPEN( InFileUnit, FILE=TRIM(FileName), STATUS='OLD' )
       CALL LoadInputFile( Model, InFileUnit, FileName, &
            MeshDir, MeshName, .FALSE., ScanOnly )
       CLOSE( InFileUnit )
     END IF

!------------------------------------------------------------------------------
  END SUBROUTINE LoadIncludeFile
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> This subroutine is used to reload input from the file to allow
!> multiple parameter sets within the same simulation.
!------------------------------------------------------------------------------
  FUNCTION ReloadInputFile( Model ) RESULT(got)
!------------------------------------------------------------------------------
    LOGICAL :: got
    TYPE(Model_t) :: Model
    INTEGER :: pos, posn
    CHARACTER(LEN=MAX_NAME_LEN) :: MeshDir, MeshName

    CALL Info('ReloadInputFile','Realoading input file',Level=7)
    MeshDir  = ' '
    Meshname = ' '
    CALL LoadInputFile( Model, InFileUnit, ' ', &
        MeshDir, MeshName, .FALSE., .FALSE., got )
!------------------------------------------------------------------------------
  END FUNCTION ReloadInputFile
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Subroutine for loading the fields from the input file. This may be 
!> used multiple times. On first calling the ScanOnly should be true 
!> to allow only the reading of number of fields. 
!> Performs also some simple sanity tests for the lists. 
!------------------------------------------------------------------------------
  RECURSIVE SUBROUTINE LoadInputFile( Model, InFileUnit, FileName, &
         MeshDir, MeshName, BaseLoad, ScanOnly, Runc )
!------------------------------------------------------------------------------

    CHARACTER(LEN=*) :: FileName
    TYPE(Model_t), TARGET :: Model
    INTEGER :: InFileUnit, iostat
    LOGICAL :: BaseLoad
    LOGICAL :: ScanOnly
    LOGICAL, OPTIONAL :: runc
    CHARACTER(LEN=*) :: MeshDir,MeshName
!------------------------------------------------------------------------------

    TYPE( ValueList_t ), POINTER :: List

    INTEGER :: i,j,k,n,Arrayn,TYPE,Sect,N1,N2,BoundaryIndex

    INTEGER(KIND=AddrInt) :: Proc

    CHARACTER(LEN=:), ALLOCATABLE :: section, name, str

    LOGICAL :: Found, SizeGiven, FoundName
    LOGICAL :: FreeNames=.FALSE., Echo = .FALSE., Numbering = .TRUE.
    INTEGER :: CheckAbort = 0

    TYPE(Solver_t), POINTER :: ASolvers(:)

    TYPE(ComponentArray_t), POINTER  :: AComponent(:)
    TYPE(MaterialArray_t), POINTER  :: AMaterial(:)
    TYPE(EquationArray_t), POINTER  :: AEquation(:)
    TYPE(BodyArray_t), POINTER      :: ABody(:)
    TYPE(BodyForceArray_t), POINTER :: ABF(:)
    TYPE(InitialConditionArray_t), POINTER  :: AIC(:)
    TYPE(BoundaryConditionArray_t), POINTER :: ABC(:)
    LOGICAL, ALLOCATABLE :: EntryUsed(:)

    LOGICAL :: FirstTime = .TRUE.

    INTEGER :: nlen, BCcount, BodyCount, EqCount, MatCount, BfCount, &
        IcCount, SolverCount, LineCount, ComponentCount
    REAL(KIND=dp) :: Val

!------------------------------------------------------------------------------

    ALLOCATE(CHARACTER(MAX_STRING_LEN)::section)
    ALLOCATE(CHARACTER(MAX_STRING_LEN)::str)
    ALLOCATE(CHARACTER(MAX_STRING_LEN)::name)

    IF( ScanOnly ) THEN
      CALL Info('LoadInputFile','Scanning input file: '//TRIM(FileName),Level=7)
    ELSE
      CALL Info('LoadInputFile','Loading input file: '//TRIM(FileName),Level=7)
    END IF
      
!------------------------------------------------------------------------------
!   Read model header first
!------------------------------------------------------------------------------
    CheckAbort = 1
    IF ( FirstTime ) THEN
       CALL CheckKeyword( 'coordinate system', 'string', &
           CheckAbort,FreeNames,'simulation' )
    END IF
 
    IF ( BaseLoad ) THEN
      DO WHILE( ReadAndTrim( InFileUnit, Name, Echo ) )
        IF ( Name=='' .OR. Name==' ' )   CYCLE
        IF (Name == 'end' ) EXIT

        IF ( SEQL(Name,'check keywords') ) THEN
           k = 16
           IF ( Name(k:k) == '"' ) k = k + 1
           SELECT CASE(Name(k:))
           CASE('ignore')
             CheckAbort = 0
           CASE('warn')
             CheckAbort = 1
           CASE('silent')
             CheckAbort = 2
           CASE('abort')
             CheckAbort = 3
           END SELECT
        ELSE IF ( Name == 'echo on' ) THEN
           Echo = .TRUE.
        ELSE IF ( Name == 'echo off' ) THEN
           Echo = .FALSE.
        ELSE IF ( Name == 'numbering on' ) THEN
           Numbering = .TRUE.
        ELSE IF ( Name == 'numbering off' ) THEN
           Numbering = .FALSE.
        ELSE IF ( SEQL(Name, 'bodies') ) THEN
        ELSE IF ( SEQL(Name, 'initial conditions') ) THEN
        ELSE IF ( SEQL(Name, 'boundaries') ) THEN
        ELSE IF ( SEQL(Name, 'boundary conditions') ) THEN
        ELSE IF ( SEQL(Name, 'components') ) THEN
        ELSE IF ( SEQL(Name, 'equations') ) THEN
        ELSE IF ( SEQL(Name, 'solvers') ) THEN
        ELSE IF ( SEQL(Name, 'materials') ) THEN
        ELSE IF ( SEQL(Name, 'body forces') ) THEN
        ELSE IF ( SEQL(Name, 'mesh db') ) THEN
          k = 9
          i = 1
          nlen = LEN_TRIM(Name)
          DO WHILE( Name(k:k) /= ' ' )
            MeshDir(i:i)  = Name(k:k)
            Meshname(i:i) = Name(k:k)
            k = k + 1
            i = i + 1
          END DO
          MeshDir(i:i) = CHAR(0)

          DO WHILE( k<=nlen .AND. Name(k:k) == ' ' )
            k = k + 1
          END DO

          IF ( k<=nlen ) THEN
             MeshName(i:i) = '/'
             i = i + 1
             DO WHILE( Name(k:k) /= ' ' )
               MeshName(i:i) = Name(k:k)
               k = k + 1
               i = i + 1
             END DO
          ELSE
             MeshDir = "." // CHAR(0)
          END IF
          MeshName(i:i) = CHAR(0)
        ELSE IF ( SEQL(Name,'header') ) THEN
        ELSE IF ( SEQL(Name,'include path') ) THEN
           IncludePath = Name(14:)
        ELSE IF ( SEQL(Name,'results directory') ) THEN
           OutputPath = Name(19:)
        ELSE IF ( SEQL(Name,'simulation id') ) THEN
           SimulationId = Name(15:)
        ELSE
          WRITE( Message, * ) 'Unknown input field in header section: ' // TRIM(Name)
          CALL Fatal( 'Model Input',  Message )
        END IF
      END DO

      Model % BCs => NULL()
      Model % ICs => NULL()
      Model % Bodies => NULL()
      Model % Solvers => NULL()
      Model % Components => NULL()
      Model % Equations => NULL()
      Model % Materials => NULL()
      Model % BodyForces => NULL()
      Model % Boundaries => NULL()
      Model % Constants => NULL()
      Model % Simulation => NULL()
    END IF

!------------------------------------------------------------------------------
    IF ( .NOT. ScanOnly ) THEN
       IF ( .NOT.ASSOCIATED( Model % Boundaries ) ) THEN
         ALLOCATE( Model % Boundaries(Model % NumberOfBoundaries) )
         ALLOCATE( Model % BoundaryId(Model % NumberOfBoundaries) )
         DO i=1,Model % NumberOfBoundaries
           Model % BoundaryId(i) = 0
         END DO
         BoundaryIndex = 0
       END IF
    END IF

    BCcount = 0
    BodyCount = 0
    BfCount = 0
    EqCount = 0
    MatCount = 0
    IcCount = 0
    SolverCount = 0
    ComponentCount = 0
    LineCount = 0
    
    IF ( PRESENT(runc) ) runc = .FALSE.
!------------------------------------------------------------------------------
    DO WHILE(ReadAndTrim(InFileUnit,Section,Echo))
!------------------------------------------------------------------------------
      IF ( Section == '' .OR. Section == ' ' ) CYCLE

      IF ( SEQL(Section,'include') ) THEN
        CALL LoadIncludeFile( Model, InFileUnit-1, Section(9:), &
                    MeshDir, MeshName, ScanOnly )
        CYCLE
      END IF

      IF ( SEQL(Section, 'header') ) THEN
         DO WHILE( ReadAndTrim( InFileUnit, Section, Echo ) )
            IF ( Section == 'end' ) EXIT
         END DO
         CYCLE
      ELSE IF ( SEQL(Section, 'echo ') .OR. SEQL(Section, 'check ')) THEN
         CYCLE
      ELSE IF ( Section == 'run' ) THEN
         IF ( PRESENT(runc) ) runc=.TRUE.
         EXIT
      END IF

      FreeNames = ( CheckAbort <= 0 )
      ArrayN = 0
      LineCount = LineCount + 1
      
      IF ( SEQL(Section, 'constants') ) THEN
        IF ( .NOT. ScanOnly ) THEN
          ArrayN = 1
          IF(.NOT.ASSOCIATED(Model % Constants)) &
              Model % Constants => ListAllocate()
          List => Model % Constants
        END IF
        
      ELSE IF ( SEQL(Section, 'simulation') ) THEN        
        IF ( .NOT. ScanOnly ) THEN
          ArrayN = 1
          IF(.NOT.ASSOCIATED(Model % Simulation)) &
              Model % Simulation=>ListAllocate()
          List => Model % Simulation
        END IF
        
      ELSE IF ( SEQL(Section, 'boundary condition') ) THEN
        
        READ( Section(19:),*,iostat=iostat ) Arrayn
        IF( iostat /= 0 ) THEN
          IF( Numbering ) THEN
            CALL Fatal('LoadInputFile','Problem reading section '&
                //TRIM(I2S(LineCount))//': '//TRIM(Section))
          END IF
          BCcount = BCcount + 1
          ArrayN = BCcount 
          IF( ScanOnly ) THEN
            CALL Info('LoadInputFile','Giving an empty > Boundary Condition < index next value: &
                '//TRIM(I2S(ArrayN)),Level=4)
          END IF
        ELSE
          BcCount = MAX( BcCount, ArrayN )
        END IF
          
        IF ( ScanOnly ) THEN
          Model % NumberOFBCs = MAX( Model % NumberOfBCs, Arrayn )
        ELSE          
          IF ( .NOT.ASSOCIATED( Model % BCs ) ) THEN
            ALLOCATE( Model % BCs(Model % NumberOfBCs) )
          ELSE             
            Model % NumberOfBCs = MAX( Arrayn, Model % NumberOfBCs )
            
            IF ( SIZE( Model % BCs ) < Model % NumberOfBCs ) THEN
              ALLOCATE( ABC(Model % NumberOfBCs) )
              DO i=1,SIZE(Model % BCs)
                ABC(i) % Values => Model % BCs(i) % Values
              END DO
              DEALLOCATE( Model % BCs )
              Model % BCs => ABC
            END IF
          END IF
          
          DO i=1,Model % NUmberOfBCs
            IF(.NOT.ASSOCIATED(Model % BCs(i) % Values)) &
                Model % BCs(i) % Values => ListAllocate()
          END DO
          
          IF ( Arrayn <= 0 .OR. Arrayn > Model % NumberOfBCs ) THEN
            WRITE( Message, * ) 'Boundary Condition section number ('//TRIM(I2S(Arrayn))// &
                ') exceeds number of BCs ('//TRIM(I2S(Model % NumberOfBCs))//')'
            CALL Fatal( 'Model Input', Message )
          END IF
          Model % BCs(ArrayN) % Tag = ArrayN
          List => Model % BCs(Arrayn) % Values
        END IF
        
        FreeNames = .TRUE.

      ELSE IF ( SEQL(Section, 'boundary') ) THEN

        IF ( ScanOnly ) THEN
          Model % NumberOfBoundaries = Model % NumberOfBoundaries + 1
        ELSE
          IF ( .NOT.ASSOCIATED( Model % Boundaries ) ) THEN
            ALLOCATE( Model % Boundaries(Model % NumberOfBoundaries) )
            ALLOCATE( Model % BoundaryId(Model % NumberOfBoundaries) )
            DO i=1,Model % NumberOfBoundaries
              Model % BoundaryId(i) = 0
            END DO
          END IF

          READ( Section(9:),*,iostat=iostat ) Arrayn
          IF( iostat /= 0 ) THEN
            CALL Fatal('LoadInputFile','Problem reading section '&
                //TRIM(I2S(LineCount))//': '//TRIM(Section))
          END IF

          BoundaryIndex = BoundaryIndex + 1
          IF ( BoundaryIndex <= 0 .OR. BoundaryIndex >  &
              Model % NumberOfBoundaries ) THEN
            WRITE( Message, * ) 'Boundary section number: ',BoundaryIndex, &
                ' exeeds header value.'
            CALL Fatal( 'Model Input', Message )
          END IF
          Model % BoundaryId(BoundaryIndex) = Arrayn
          List => Model % Boundaries(BoundaryIndex) % Values
        END IF

      ELSE IF ( SEQL(Section, 'initial condition') ) THEN

        READ( Section(18:),*,iostat=iostat ) Arrayn
        IF( iostat /= 0 ) THEN
          IF( Numbering ) THEN               
            CALL Fatal('LoadInputFile','Problem reading section '&
                //TRIM(I2S(LineCount))//': '//TRIM(Section))
          END IF
          IcCount = IcCount + 1
          ArrayN = IcCount 
          IF( ScanOnly ) THEN
            CALL Info('LoadInputFile','Giving an empty > Initial Condition < index next value: &
                '//TRIM(I2S(ArrayN)),Level=4)
          END IF
        ELSE
          IcCount = MAX( IcCount, ArrayN ) 
        END IF
        
        IF ( ScanOnly ) THEN
          Model % NumberOFICs = MAX( Model % NumberOfICs, ArrayN )
        ELSE
          IF ( .NOT.ASSOCIATED( Model % ICs ) ) THEN
            ALLOCATE( Model % ICs(Model % NumberOfICs) )
          ELSE
            Model % NumberOfICs = MAX( Model % NumberOfICs, Arrayn )
            IF ( SIZE( Model % ICs ) < Model % NumberOfICs ) THEN
              ALLOCATE( AIC(Model % NumberOfICs) )
              DO i=1,SIZE(Model % ICs)
                AIC(i) % Values => Model % ICs(i) % Values
              END DO
              DEALLOCATE( Model % ICs )
              Model % ICs => AIC
            END IF
          END IF
          
          DO i=1,Model % NUmberOfICs
            IF(.NOT.ASSOCIATED(Model % ICs(i) % Values)) &
                Model % ICs(i) % Values => ListAllocate()
          END DO
          
          IF ( Arrayn <= 0 .OR. Arrayn > Model % NumberOfICs ) THEN
            WRITE( Message, * ) 'Initial Condition section number: ',Arrayn, &
                ' exeeds header value.'
            CALL Fatal( 'Model Input', Message )
          END IF
          Model % ICs(ArrayN) % Tag = ArrayN
          List => Model % ICs(Arrayn) % Values
        END IF

        FreeNames = .TRUE.

      ELSE IF ( SEQL(Section, 'material') ) THEN

        READ( Section(9:),*,iostat=iostat ) Arrayn
        
        IF( iostat /= 0 ) THEN
          IF( Numbering ) THEN               
            CALL Fatal('LoadInputFile','Problem reading section '&
                //TRIM(I2S(LineCount))//': '//TRIM(Section))
          END IF
          MatCount = MatCount + 1
          ArrayN = MatCount 
          IF( ScanOnly ) THEN
            CALL Info('LoadInputFile','Giving an empty > Material < index next value: &
                '//TRIM(I2S(ArrayN)),Level=4)
          END IF
        ELSE
          MatCount = MAX( MatCount, ArrayN ) 
        END IF
        
        IF ( ScanOnly ) THEN
          Model % NumberOFMaterials = MAX( Model % NumberOfMaterials, ArrayN )
        ELSE
          IF ( .NOT.ASSOCIATED( Model % Materials ) ) THEN
            ALLOCATE( Model % Materials(Model % NumberOfMaterials) )
          ELSE
            Model % NumberOfMaterials = MAX( Arrayn, Model % NumberOFMaterials ) 
            IF ( SIZE( Model % Materials ) < Model % NumberOfMaterials ) THEN
              ALLOCATE( AMaterial(Model % NumberOfMaterials) )
              DO i=1,SIZE(Model % Materials)
                AMaterial(i) % Values => Model % Materials(i) % Values
              END DO
              DEALLOCATE( Model % Materials )
              Model % Materials => AMaterial
            END IF
          END IF
          
          DO i=1,Model % NumberOfMaterials
            IF(.NOT.ASSOCIATED(Model % Materials(i) % Values)) &
                Model % Materials(i) % Values => ListAllocate()
          END DO
          
          IF ( Arrayn <= 0 .OR. Arrayn > Model % NumberOfMaterials ) THEN
            WRITE( Message, * ) 'Material section number: ',Arrayn, &
                ' exeeds header value.'
            CALL Fatal( 'Model Input', Message )
          END IF
          List => Model % Materials(Arrayn) % Values
        END IF

      ELSE IF ( SEQL(Section, 'body force') ) THEN

        READ( Section(12:),*,iostat=iostat ) Arrayn
        IF( iostat /= 0 ) THEN
          IF( Numbering ) THEN               
            CALL Fatal('LoadInputFile','Problem reading section '&
                //TRIM(I2S(LineCount))//': '//TRIM(Section))
          END IF
          BfCount = BfCount + 1
          ArrayN = BfCount 
          IF( ScanOnly ) THEN
            CALL Info('LoadInputFile','Giving an empty > Body Force < index next value: &
                '//TRIM(I2S(ArrayN)),Level=4)
          END IF
        ELSE
          BfCount = MAX( ArrayN, BfCount ) 
        END IF
        
        IF ( ScanOnly ) THEN
          Model % NumberOFBodyForces = BfCount
        ELSE
          IF ( .NOT.ASSOCIATED( Model % BodyForces ) ) THEN
            ALLOCATE( Model % BodyForces(Model % NumberOfBodyForces) )
          ELSE
            Model % NumberOFBodyForces = MAX( Arrayn, Model % NumberOfBodyForces )
            IF ( SIZE( Model % BodyForces ) < Model % NumberOfBodyForces ) THEN
              ALLOCATE( ABF(Model % NumberOfBodyForces) )
              DO i=1,SIZE(Model % BodyForces)
                ABF(i) % Values => Model % BodyForces(i) % Values
              END DO
              DEALLOCATE( Model % BodyForces )
              Model % BodyForces => ABF
            END IF
          END IF
          
          DO i=1,Model % NumberOfBodyForces
            IF(.NOT.ASSOCIATED(Model % BodyForces(i) % Values)) &
                Model % BodyForces(i) % Values => ListAllocate()
          END DO
          
          IF ( Arrayn <= 0 .OR. Arrayn > Model % NumberOfBodyForces ) THEN
            WRITE( Message, * ) 'Body Force section number: ',Arrayn, &
                ' exeeds header value.'
            CALL Fatal( 'Model Input', Message )
          END IF
          List => Model % BodyForces(Arrayn) % Values
        END IF
        
      ELSE IF ( SEQL(Section, 'equation') ) THEN

        READ( Section(9:),*,iostat=iostat ) Arrayn
        IF( iostat /= 0 ) THEN
          IF( Numbering ) THEN               
            CALL Fatal('LoadInputFile','Problem reading section '&
                //TRIM(I2S(LineCount))//': '//TRIM(Section))
          END IF
          EqCount = EqCount + 1
          ArrayN = EqCount 
          IF( ScanOnly ) THEN
            CALL Info('LoadInputFile','Giving an empty > Equation < index next value: &
                '//TRIM(I2S(ArrayN)),Level=4)
          END IF
        ELSE
          EqCount = MAX( EqCount, ArrayN )
        END IF
          
        IF ( ScanOnly ) THEN
          Model % NUmberOfEquations = MAX( Model % NumberOFEquations, ArrayN )
        ELSE
          IF ( .NOT.ASSOCIATED( Model % Equations ) ) THEN
            ALLOCATE( Model % Equations(Model % NumberOfEquations) )
          ELSE
            Model % NumberOFEquations = MAX( Arrayn, Model % NumberOFEquations )
            IF ( SIZE( Model % Equations ) < Model % NumberOfEquations ) THEN
              ALLOCATE( AEquation(Model % NumberOfEquations) )
              DO i=1,SIZE(Model % Equations)
                AEquation(i) % Values => Model % Equations(i) % Values
              END DO
              DEALLOCATE( Model % Equations )
              Model % Equations => AEquation
            END IF
          END IF
          
          DO i=1,Model % NumberOfEquations
            IF(.NOT.ASSOCIATED(Model % Equations(i) % Values)) &
                Model % Equations(i) % Values => ListAllocate()
          END DO
          
          IF ( Arrayn <= 0 .OR. Arrayn > Model % NumberOfEquations ) THEN
            WRITE( Message, * ) 'Equation section number: ',Arrayn, &
                ' exeeds header value.'
            CALL Fatal( 'Model Input', Message )
          END IF
          List => Model % Equations(ArrayN) % Values
        END IF

        FreeNames = .TRUE.

      ELSE IF ( SEQL(Section, 'body') ) THEN

        READ( Section(5:),*,iostat=iostat ) Arrayn
        IF( iostat /= 0 ) THEN
          IF( Numbering ) THEN               
            CALL Fatal('LoadInputFile','Problem reading section '&
                //TRIM(I2S(LineCount))//': '//TRIM(Section))
          END IF
          BodyCount = BodyCount + 1
          ArrayN = BodyCount 
          IF( ScanOnly ) THEN
            CALL Info('LoadInputFile','Giving an empty > Body < index next value: &
                '//TRIM(I2S(ArrayN)),Level=4)
          END IF
        ELSE
          BodyCount = MAX( BodyCount, ArrayN ) 
        END IF
        
        IF ( ScanOnly ) THEN
          Model % NumberOFBodies = BodyCount 
        ELSE
          IF ( .NOT.ASSOCIATED( Model % Bodies ) ) THEN
            ALLOCATE( Model % Bodies(Model % NumberOfBodies) )
          ELSE
            Model % NumberOFBodies = MAX( Arrayn, Model % NumberOFBodies )
            IF ( SIZE( Model % Bodies ) < Model % NumberOfBodies ) THEN
              ALLOCATE( ABody(Model % NumberOfBodies) )
              DO i=1,SIZE(Model % Bodies)
                ABody(i) % Values => Model % Bodies(i) % Values
              END DO
              DEALLOCATE( Model % Bodies )
              Model % Bodies => ABody
            END IF
          END IF
          
          DO i=1,Model % NumberOfBodies
            IF(.NOT.ASSOCIATED(Model % Bodies(i) % Values)) &
                Model % Bodies(i) % Values => ListAllocate()
          END DO
          
          IF ( Arrayn <= 0 .OR. Arrayn > Model % NumberOfBodies ) THEN
            WRITE( Message, * ) 'Body section number: ',Arrayn, &
                ' exeeds header value. Aborting. '
            CALL Fatal( 'Model Input', Message )
          END IF
          List => Model % Bodies(Arrayn) % Values
        END IF

      ELSE IF ( SEQL(Section, 'component') ) THEN

        READ( Section(10:),*,iostat=iostat ) Arrayn
        
        IF( iostat /= 0 ) THEN
          IF( Numbering ) THEN
            CALL Fatal('LoadInputFile','Problem reading section '&
                //TRIM(I2S(LineCount))//': '//TRIM(Section))
          END IF
          ComponentCount = ComponentCount + 1
          ArrayN = ComponentCount 
          IF( ScanOnly ) THEN
            CALL Info('LoadInputFile','Giving an empty > Component < index next value: &
                '//TRIM(I2S(ArrayN)),Level=4)
          END IF
        ELSE
          ComponentCount = MAX( ComponentCount, ArrayN )
        END IF
        
        IF ( ScanOnly ) THEN
          Model % NumberOFComponents = ComponentCount  
        ELSE
          IF ( .NOT.ASSOCIATED( Model % Components ) ) THEN
            ALLOCATE( Model % Components(Model % NumberOfComponents) )
          ELSE
            Model % NumberOFComponents = MAX( Arrayn, Model % NumberOFComponents )
            IF ( SIZE( Model % Components ) < Model % NumberOfComponents ) THEN
              ALLOCATE( ABody(Model % NumberOfComponents) )
              DO i=1,SIZE(Model % Components)
                AComponent(i) % Values % head => Model % Components(i) % Values % head
              END DO
              DEALLOCATE( Model % Components )
              Model % Components => AComponent
            END IF
          END IF
          
          DO i=1,Model % NumberOfComponents
            IF(.NOT.ASSOCIATED(Model % Components(i) % Values)) &
                Model % Components(i) % Values => ListAllocate()
          END DO
          
          IF ( Arrayn <= 0 .OR. Arrayn > Model % NumberOfComponents ) THEN
            WRITE( Message, * ) 'Component section number: ',Arrayn, &
                ' exeeds header value. Aborting. '
            CALL Fatal( 'Model Input', Message )
          END IF
          List => Model % Components(Arrayn) % Values
        END IF

      ELSE IF ( SEQL(Section, 'solver') ) THEN

        READ( Section(7:),*,iostat=iostat ) Arrayn
        IF( iostat /= 0 ) THEN
          IF( Numbering ) THEN
            CALL Fatal('LoadInputFile','Problem reading section '&
                //TRIM(I2S(LineCount))//': '//TRIM(Section))
          END IF
          SolverCount = SolverCount + 1
          ArrayN = SolverCount
          IF( ScanOnly ) THEN
            CALL Info('LoadInputFile','Giving an empty > Solver < index next value: &
                '//TRIM(I2S(ArrayN)),Level=4)
          END IF
        ELSE
          SolverCount = MAX( SolverCount, ArrayN )
        END IF
        
        IF ( ScanOnly ) THEN
          Model % NumberOfSolvers = SolverCount
        ELSE
          IF ( .NOT.ASSOCIATED( Model % Solvers ) ) THEN
            ALLOCATE( Model % Solvers(Model % NumberOfSolvers) )
            DO i=1,Model % NumberOfSolvers
              Model % Solvers(i) % PROCEDURE = 0
              NULLIFY( Model % Solvers(i) % Matrix )
              NULLIFY( Model % Solvers(i) % Variable )
              NULLIFY( Model % Solvers(i) % ActiveElements )
              Model % Solvers(i) % NumberOfActiveElements = 0
              Model % Solvers(i) % SolverId = i
            END DO
          ELSE
            Model % NumberOfSolvers = MAX( Arrayn, Model % NumberOfSolvers )
            IF ( SIZE(Model % Solvers) < Model % NumberOfSolvers ) THEN
              ALLOCATE( ASolvers(Model % NumberOfSolvers) )
              DO i=1,SIZE(Model % Solvers)
                ASolvers(i) = Model % Solvers(i)
              END DO
              DO i=SIZE(Model % Solvers)+1,Model % NumberOfSolvers
                ASolvers(i) % PROCEDURE = 0
                NULLIFY( ASolvers(i) % Matrix )
                NULLIFY( ASolvers(i) % Mesh )
                NULLIFY( ASolvers(i) % Variable )
                NULLIFY( ASolvers(i) % ActiveElements )
                ASolvers(i) % NumberOfActiveElements = 0
              END DO
              DEALLOCATE( Model % Solvers )
              Model % Solvers => ASolvers
            END IF
          END IF
          
          DO i=1,Model % NumberOfSolvers
            IF(.NOT.ASSOCIATED(Model % Solvers(i) % Values)) &
                Model % Solvers(i) % Values => ListAllocate()
          END DO
          
          IF ( Arrayn <= 0 .OR. Arrayn > Model % NumberOfSolvers ) THEN
            WRITE( Message, * ) 'Solver section number: ',Arrayn, &
                ' exeeds header value. Aborting. '
            CALL Fatal( 'Model Input', Message )
          END IF
          List => Model % Solvers(Arrayn) % Values
        END IF
      ELSE
        WRITE( Message, * ) 'Unknown input section name: ',TRIM(Section)
        CALL Fatal( 'Model Input', Message )
      END IF
!------------------------------------------------------------------------------
      IF ( .NOT. ScanOnly .AND. ArrayN == 0 ) CYCLE

      CALL SectionContents( Model, List, CheckAbort, FreeNames, &
                Section, InFileUnit, ScanOnly, Echo )
!------------------------------------------------------------------------------
    END DO
!------------------------------------------------------------------------------

    IF ( BaseLoad .AND. .NOT. ScanOnly )  THEN

      ! Make some sanity checks that all the entries have been defined
      ! The might be missing entries due to duplicate numbering etc.
      ! For some fields this is not detrimental, thus just a warning.
      !--------------------------------------------------------------------
      DO i = 1, Model % NumberOFBCs
        IF( ListEmpty(Model % BCs(i) % Values) ) THEN
          WRITE( Message,'(A,I0)') 'Entry missing for: Boundary Condition ',i
          CALL Warn('LoadInputFile',Message)
        END IF
      END DO
      
      DO i = 1, Model % NumberOfBodyForces
        IF( ListEmpty(Model % BodyForces(i) % Values) ) THEN
          WRITE( Message,'(A,I0)') 'Entry missing for: Body Force ',i
          CALL Warn('LoadInputFile',Message)
        END IF
      END DO

      DO i = 1, Model % NumberOfICs
        IF( ListEmpty(Model % ICs(i) % Values) ) THEN
          WRITE( Message,'(A,I0)') 'Entry missing for: Initial Condition ',i
          CALL Warn('LoadInputFile',Message)
        END IF
      END DO
      
      DO i = 1, Model % NumberOfMaterials         
        IF( ListEmpty(Model % Materials(i) % Values) ) THEN
          WRITE( Message,'(A,I0)') 'Entry missing for: Material ',i
          CALL Warn('LoadInputFile',Message)
        END IF
      END DO

      DO i = 1, Model % NumberOFEquations 
        IF( ListEmpty(Model % Equations(i) % Values) ) THEN
          WRITE( Message,'(A,I0)') 'Entry missing for: Equation ',i
          CALL Fatal('LoadInputFile',Message)
        END IF
      END DO
    
      DO i = 1, Model % NumberOfSolvers         
        IF( ListEmpty(Model % Solvers(i) % Values) ) THEN
          WRITE( Message,'(A,I0)') 'Entry missing for: Solver ',i
          CALL Fatal('LoadInputFile',Message)
        END IF
      END DO

      DO i = 1, Model % NumberOfBodies
        IF( ListEmpty(Model % Bodies(i) % Values) ) THEN
          WRITE( Message,'(A,I0)') 'Entry missing for: Body ',i
          CALL Fatal('LoadInputFile',Message)
        END IF
      END DO

      ! Check that the same name is not used twice,
      ! each solver should be uniquely defined by its name.
      ! If the same equation name is used twice it may lead to difficult-to-find 
      ! problems later on. 
      !-----------------------------------------------------------------------------------
      DO i = 1, Model % NumberOfSolvers
        str = ListGetString( Model % Solvers(i) % Values,'Equation',Found )
        IF(.NOT. Found ) CYCLE
        DO j = i+1, Model % NumberOfSolvers
          IF( TRIM(str) == TRIM( ListGetString( Model % Solvers(j) % Values,'Equation',Found ))) THEN
            CALL Fatal('LoadInputFile','Solvers '//TRIM(I2S(i))//' and '//TRIM(I2S(j))//&
                ' have the same Equation name!')
          END IF
        END DO
      END DO
            
      ! If automatic numbering is used map the names to numbers
      !------------------------------------------------------------
      IF( .NOT. Numbering ) THEN
        DO i = 1, Model % NumberOfBodies
          IF( .NOT. ListCheckPresent( Model % Bodies(i) % Values,'Material') ) THEN
            name = ListGetString( Model % Bodies(i) % Values,'Material Name', Found )
            IF(.NOT. Found ) CYCLE
            FoundName = .FALSE.
            DO j = 1,Model % NumberOfMaterials
              str = ListGetString( Model % Materials(j) % Values,'Name',Found )
              IF(.NOT. Found ) CYCLE
              IF( str == name ) THEN
                CALL ListAddInteger( Model % Bodies(i) % Values,'Material',j)
                CALL Info('LoadInputFile','Giving material > '//TRIM(Name)//' < index: '//TRIM(I2S(j)),Level=5)
                FoundName = .TRUE.
                EXIT
              END IF
            END DO
            IF(.NOT. FoundName ) THEN
              CALL Fatal('LoadInputFile','> Material Name = '//TRIM(name)//&
                  ' < given but no such material exists!')
            END IF            
          END IF
            
          IF( .NOT. ListCheckPresent( Model % Bodies(i) % Values,'Equation') ) THEN
            name = ListGetString( Model % Bodies(i) % Values,'Equation Name', Found )
            IF(.NOT. Found ) CYCLE
            FoundName = .FALSE.
            DO j = 1,Model % NumberOfEquations
              str = ListGetString( Model % Equations(j) % Values,'Name',Found )
              IF(.NOT. Found ) CYCLE
              IF( str == name ) THEN
                CALL ListAddInteger( Model % Bodies(i) % Values,'Equation',j)
                CALL Info('LoadInputFile','Giving equation > '//TRIM(Name)//' < index: '//TRIM(I2S(j)),Level=5)
                FoundName = .TRUE.
                EXIT
              END IF
            END DO
            IF(.NOT. FoundName ) THEN
              CALL Fatal('LoadInputFile','> Equation Name = '//TRIM(name)//&
                  ' < given but no such equation exists!')
            END IF            
          END IF

          IF( .NOT. ListCheckPresent( Model % Bodies(i) % Values,'Body Force') ) THEN
            name = ListGetString( Model % Bodies(i) % Values,'Body Force Name', Found )
            IF(.NOT. Found ) CYCLE
            FoundName = .FALSE.
            DO j = 1,Model % NumberOfBodyForces
              str = ListGetString( Model % BodyForces(j) % Values,'Name',Found )
              IF(.NOT. Found ) CYCLE
              IF( str == name ) THEN
                CALL ListAddInteger( Model % Bodies(i) % Values,'Body Force',j)
                CALL Info('LoadInputFile','Giving body force > '//TRIM(Name)//' < index: '//TRIM(I2S(j)),Level=5)
                FoundName = .TRUE.
                EXIT
              END IF
            END DO
            IF(.NOT. FoundName ) THEN
              CALL Fatal('LoadInputFile','> Body Force Name = '//TRIM(name)//&
                  ' < given but no such body force exists!')
            END IF
          END IF

          IF( .NOT. ListCheckPresent( Model % Bodies(i) % Values,'Initial Condition') ) THEN
            name = ListGetString( Model % Bodies(i) % Values,'Initial Condition Name', Found )
            IF(.NOT. Found ) CYCLE
            FoundName = .FALSE.
            DO j = 1,Model % NumberOfICs
              str = ListGetString( Model % ICs(j) % Values,'Name',Found )
              IF(.NOT. Found ) CYCLE
              IF( str == name ) THEN
                CALL ListAddInteger( Model % Bodies(i) % Values,'Initial Condition',j)
                CALL Info('LoadInputFile','Giving initial condition > '//TRIM(Name)//' < index: '//TRIM(I2S(j)),Level=5)
                FoundName = .TRUE.
                EXIT
              END IF
            END DO
            IF(.NOT. FoundName ) THEN
              CALL Fatal('LoadInputFile','> Initial Condition Name = '//TRIM(name)//&
                  ' < given but no such initial condition exists!')
            END IF            
          END IF
               
        END DO ! number of bodies
      END IF
      

      ! Make sanity check that all Material, Body Force and Equation is associated to some
      ! body. This is not detrimental so a warning suffices.
      !-----------------------------------------------------------------------------------
      n = MAX( Model % NumberOfMaterials, &
          Model % NumberOfBodyForces, &
          Model % NumberOfEquations )
      ALLOCATE( EntryUsed( n ) )

      IF( Model % NumberOfMaterials > 0 ) THEN
        EntryUsed = .FALSE.
        DO i = 1, Model % NumberOfBodies
          j = ListGetInteger( Model % Bodies(i) % Values,'Material',Found )
          IF( Found ) EntryUsed(j) = .TRUE.
        END DO
        DO i = 1,Model % NumberOfMaterials
          IF( .NOT. EntryUsed(i) ) THEN
            CALL Warn('LoadInputFile','> Material '// TRIM(I2S(i)) //' < not used in any Body!')
          END IF
        END DO
      END IF

      IF( Model % NumberOfBodyForces > 0 ) THEN
        EntryUsed = .FALSE.
        DO i = 1, Model % NumberOfBodies
          j = ListGetInteger( Model % Bodies(i) % Values,'Body Force',Found )
          IF( Found ) EntryUsed(j) = .TRUE.
        END DO
        DO i = 1,Model % NumberOfBodyForces
          IF( .NOT. EntryUsed(i) ) THEN
            CALL Warn('LoadInputFile','> Body Force '// TRIM(I2S(i)) //' < not used in any Body!')
          END IF
        END DO
      END IF

      IF( Model % NumberOfEquations > 0 ) THEN
        EntryUsed = .FALSE.
        DO i = 1, Model % NumberOfBodies
          j = ListGetInteger( Model % Bodies(i) % Values,'Equation',Found )
          IF( Found ) EntryUsed(j) = .TRUE.
        END DO
        DO i = 1,Model % NumberOfEquations
          IF( .NOT. EntryUsed(i) ) THEN
            CALL Warn('LoadInputFile','> Equation '// TRIM(I2S(i)) //' < not used in any Body!')
          END IF
        END DO
      END IF

      DEALLOCATE( EntryUsed ) 
      !--- sanity checks done

      ! Add default equation, material, ic, bodyforce, and body if not given:
      ! ---------------------------------------------------------------------
      IF ( Model % NumberOFEquations <= 0 ) THEN
        Model % NumberOfEquations = 1
        ALLOCATE( Model % Equations(1) )
        Model % Equations(1) % Values => ListAllocate()
        CALL ListAddIntegerArray( Model % Equations(1) % Values, 'Active Solvers', &
            Model % NumberOFSolvers, (/ (i,i=1,Model % NumberOfSolvers) /) )
        CALL ListAddString ( Model % Equations(1) % Values, 'Name', 'Default Equation 1' )
      END IF

      IF ( Model % NumberOfMaterials <= 0 ) THEN
        Model % NumberOfMaterials = 1
        ALLOCATE( Model % Materials(1) )
        Model % Materials(1) % Values => ListAllocate()
        CALL ListAddString ( Model % Materials(1) % Values, 'Name', 'Default Material 1' )
      END IF

      IF ( Model % NumberOfBodyForces <= 0 ) THEN
        Model % NumberOfBodyForces = 1
        ALLOCATE( Model % BodyForces(1) )
        Model % BodyForces(1) % Values => ListAllocate()
        CALL ListAddString ( Model % BodyForces(1) % Values, 'Name','Default Body Force 1' )
      END IF

      IF ( Model % NumberOfICs <= 0 ) THEN
        Model % NumberOfICs = 1
        ALLOCATE( Model % ICs(1) )
        Model % ICs(1) % Values => ListAllocate()
        CALL ListAddString ( Model % ICs(1) % Values, 'Name','Default IC 1' )
      END IF

      IF ( Model % NumberOfBodies <= 0 ) THEN
        Model % NumberOfBodies = 1
        ALLOCATE( Model % Bodies(1) )
        Model % Bodies(1) % Values => ListAllocate()
        CALL ListAddString(  Model % Bodies(1) % Values, 'Name', 'Default Body 1' )
        CALL ListAddInteger( Model % Bodies(1) % Values, 'Equation',   1 )
        CALL ListAddInteger( Model % Bodies(1) % Values, 'Material',   1 )
        CALL ListAddInteger( Model % Bodies(1) % Values, 'Body Force', 1 )
        CALL ListAddInteger( Model % Bodies(1) % Values, 'Initial Condition', 1 )
      END IF
      ! -- done adding default fields
               
    END IF
    !--------------------------------------------------------------------

    FirstTime = .FALSE.

    RETURN

10  CONTINUE

    WRITE( Message, * ) 'Cannot find input file: ', TRIM(FileName)
    CALL Warn( 'Model Input', Message )

CONTAINS

    SUBROUTINE CheckKeyWord( Name,TYPE,CheckAbort,FreeNames,Section,ReturnType )
       USE HashTable

       CHARACTER(LEN=*) :: Name,TYPE,Section
       INTEGER :: CheckAbort
       LOGICAL, OPTIONAL :: ReturnType,FreeNames

       INTEGER :: i,j, k,n, istat
       TYPE(HashTable_t), POINTER, SAVE :: Hash =>NULL()
       TYPE(HashValue_t), POINTER :: Val
       LOGICAL :: FirstTime = .TRUE.,lstat, fexist
       CHARACTER(LEN=:), ALLOCATABLE :: str
       CHARACTER(LEN=MAX_STRING_LEN) :: str1
!       EXTERNAL ENVIR

       IF ( PRESENT( ReturnType ) ) ReturnType = .FALSE.

       IF ( CheckAbort <= 0 ) RETURN

       ALLOCATE(CHARACTER(MAX_STRING_LEN) :: str)

       IF ( FirstTime ) THEN
!
!         First time in, read the SOLVER.KEYWORDS database, and
!         build up a local hash table for it:
! 
!         Priority is in ELMER_LIB, ELMER_HOME, and finally, if all else fails
!         use the compilation time prefix.
!         ------------------------------------------------------

#ifdef USE_ISO_C_BINDINGS
          str = 'ELMER_LIB'
#else
          str = 'ELMER_LIB'//CHAR(0)
#endif
          CALL envir( str,str1,k ) 

	  fexist = .FALSE.
          IF ( k > 0  ) THEN
             str1 = str1(1:k) // '/SOLVER.KEYWORDS'
             INQUIRE(FILE=TRIM(str1), EXIST=fexist)
          END IF
          IF (.NOT. fexist) THEN
#ifdef USE_ISO_C_BINDINGS
             str = 'ELMER_HOME'
#else
             str = 'ELMER_HOME'//CHAR(0)
#endif
             CALL envir( str,str1,k ) 
             IF ( k > 0 ) THEN
                str1 = str1(1:k) // '/share/elmersolver/lib/' // 'SOLVER.KEYWORDS'
                INQUIRE(FILE=TRIM(str1), EXIST=fexist)
             END IF
             IF ((.NOT. fexist) .AND. k>0) THEN
                str1 = str1(1:k) // '/SOLVER.KEYWORDS'
                INQUIRE(FILE=TRIM(str1), EXIST=fexist)
             END IF
          END IF
          IF (.NOT. fexist) THEN
             CALL GetSolverHome(str1, n)
             str1 = str1(1:n) // '/lib/SOLVER.KEYWORDS'
             INQUIRE(FILE=TRIM(str1), EXIST=fexist)
          END IF
          IF (.NOT. fexist) THEN
             CALL Fatal('CheckKeyWord', 'SOLVER.KEYWORDS not found')
          END IF

          OPEN( 1, FILE=TRIM(str1), STATUS='OLD', ERR=10 )

!
!         Initially 50 buckets, on avarage MAX 4 entries / bucket:
!         --------------------------------------------------------
          hash => HashCreate( 50,4 )
          IF ( .NOT. ASSOCIATED( hash ) ) THEN
             IF ( CheckAbort <= 2 ) THEN
               CALL Warn( 'Model Input', 'Can not create the hash table for SOLVER.KEYWORDS.' )
               CALL Warn( 'Model Input', 'keyword checking disabled.' )
               CheckAbort = 0
               RETURN
             ELSE
               CALL Fatal( 'Model Input','Can not create the hash table for SOLVER.KEYWORDS.' )
             END IF
          END IF

5         CONTINUE

!
!         Read the keywords file row by row and add to the hash table:
!         ------------------------------------------------------------
          DO WHILE( ReadAndTrim( 1, str ) )

             i = INDEX( str, ':' )
             j = INDEX( str, "'" )
             IF ( i <= 0 .OR. j<= 0 ) CYCLE
             str1 = str(1:i-1) // ':' //  str(j+1:LEN_TRIM(str)-1)

             ALLOCATE( Val, STAT=istat )

             IF ( istat /= 0 ) THEN
                IF ( CheckAbort <= 2 ) THEN
                  CALL Warn( 'Model Input', 'Can not allocate the hash table entry for SOLVER.KEYWORDS.' )
                  CALL Warn( 'Model Input', ' keyword checking disabled.' )
                  CheckAbort = 0
                  RETURN
                ELSE
                  CALL Fatal( 'Model Input', 'Can not allocate the hash table entry for SOLVER.KEYWORDS.' )
                END IF
             END IF

             Val % TYPE = str(i+1:j-3)

             lstat = HashAdd( hash, str1, Val )
             IF ( .NOT. lstat ) THEN
                IF ( CheckAbort <= 2 ) THEN
                   CALL Warn( 'Model Input', 'Hash table build error. Keyword checking disabled.' )
                   CheckAbort = 0
                   RETURN
                ELSE
                   CALL Fatal( 'Model Input', 'Hash table build error.' )
                END IF
             END IF
          END DO
          CLOSE(1)

          IF ( FirstTime ) THEN
             FirstTime = .FALSE.
             OPEN( 1, FILE='SOLVER.KEYWORDS', STATUS='OLD', ERR=6 )
             CALL Info( 'Model Input', 'Found local SOLVER.KEYWORDS file, ' // &
                        'adding keywords to runtime database.' )
             GOTO 5
6            CONTINUE
          END IF
       END IF

!------------------------------------------------------------------------------

        IF ( SEQL(Section, 'constants') ) THEN
           str =  'constants: '
        ELSE IF ( SEQL(Section, 'simulation') ) THEN
           str =  'simulation: '
        ELSE IF ( SEQL(Section,  'boundary condition') ) THEN
           str =  'bc: '
        ELSE IF ( SEQL(Section, 'boundary') ) THEN
          str =  'boundary: '
        ELSE IF ( SEQL(Section, 'initial condition') ) THEN
          str =  'ic: '
        ELSE IF ( SEQL(Section, 'material') ) THEN
          str =  'material: '
        ELSE IF ( SEQL(Section, 'body force') ) THEN
          str =  'bodyforce: '
        ELSE IF ( SEQL(Section, 'equation') ) THEN
          str =  'equation: '
        ELSE IF ( SEQL(Section, 'body') ) THEN
          str =  'body: '
        ELSE IF ( SEQL(Section, 'solver') ) THEN
          str =  'solver: '
        ELSE IF ( SEQL(Section, 'component') ) THEN
          str =  'component: '
        END IF

        i= INDEX(Name,':') + 1
        IF (Name(i:i)==' ') i=i+1
        str = TRIM(str) // Name(i:LEN_TRIM(Name))

!------------------------------------------------------------------------------

       Val => HashValue( Hash, str )
       IF ( ASSOCIATED( Val ) ) THEN
          IF ( PRESENT( ReturnType ) ) THEN
             ReturnType = .TRUE.
             TYPE = Val % TYPE
          END IF
          IF  ( HashEqualKeys( Val % TYPE, TYPE ) ) RETURN
       END IF

       IF ( PRESENT( ReturnType ) ) ReturnType = .FALSE.

       IF ( .NOT.ASSOCIATED(Val) .AND. (CheckAbort <= 2 .OR. FreeNames) ) THEN
          IF ( .NOT. ( ScanOnly .OR. CheckAbort == 2) ) THEN
            WRITE( Message, * ) 'Unlisted keyword: [', TRIM(name), &
                      '] in section: [', TRIM(Section), ']'
            CALL Info( 'Model Input', Message )

            ! This is intended to be activated when new keywords are checked 
            ! Generally it can be set false
            !---------------------------------------------------------------
#ifdef DEVEL_KEYWORDMISSES
            OPEN( 10,File='../SOLVER.KEYWORDS.byname',&
                STATUS='UNKNOWN',POSITION='APPEND' )
            WRITE( 10,'(A,T40,A)') TRIM(Name),TRIM(str)
            CLOSE(10)

            i = INDEX( str,':' )
            OPEN( 10,File='../SOLVER.KEYWORDS.bysection',&
                STATUS='UNKNOWN',POSITION='APPEND' )
            WRITE( 10,'(A,T22,A)') str(1:i)//TRIM(TYPE)//':',"'"//TRIM(Name)//"'"
            CLOSE(10 )
#endif
          END IF
       ELSE IF ( ASSOCIATED( Val ) ) THEN
         ! Difference between types 'string' and 'file' is just that 
         ! file is case sensitive while string is not. Hence both are ok. 
         !----------------------------------------------------------------
         IF( TRIM(Val % TYPE) == 'string' .AND. TRIM(TYPE) == 'file') THEN
           RETURN
         ELSE
           WRITE( Message, * ) 'Keyword: [', TRIM(name), &
               '] in section: [', TRIM(Section), ']',  &
               ' is given wrong type: [', TRIM(TYPE),  &
               '], should be of type: [', TRIM(Val % TYPE),']'
          CALL Fatal( 'Model Input', Message )
        END IF
       ELSE
         WRITE( Message, * ) 'Unlisted keyword: [', TRIM(name), &
             '] in section: [', TRIM(Section), '].'
         CALL Fatal( 'Model Input', Message )
       END IF

       RETURN

10     CONTINUE

       IF ( CheckAbort <= 2 ) THEN
          CALL Warn( 'Model Input', 'Keyword check requested, but SOLVER.KEYWORDS' // &
                 ' database not available.' )
       ELSE
          CALL Fatal( 'Model Input', 'Keyword check requested, but SOLVER.KEYWORDS' // &
                 ' database not available.' )
       END IF
!------------------------------------------------------------------------------
    END SUBROUTINE CheckKeyword
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
    RECURSIVE SUBROUTINE SectionContents( Model,List, CheckAbort,FreeNames, &
              Section, InFileUnit, ScanOnly, Echo )
!------------------------------------------------------------------------------
      TYPE(ValueList_t), POINTER :: List,ll
      INTEGER :: InFileUnit,CheckAbort
      TYPE(Model_t) :: Model
      LOGICAL :: FreeNames,Echo
      CHARACTER(LEN=*)  :: Section

      INTEGER, ALLOCATABLE  :: IValues(:)
      REAL(KIND=dp), ALLOCATABLE :: Atx(:,:,:), ATt(:)

      CHARACTER(LEN=MAX_NAME_LEN) :: TypeString,Keyword
      CHARACTER(LEN=:), ALLOCATABLE :: Name,str, Depname
      LOGICAL :: ReturnType, ScanOnly, String_literal,  SizeGiven, Cubic, AllInt, Monotone

      INTEGER(KIND=AddrInt) :: Proc
      INTEGER :: i,j,k,l,n,slen, str_beg, str_end, n1,n2, TYPE, &
          abuflen=0, maxbuflen=0, iostat

      ALLOCATE( ATt(1), ATx(1,1,1), IValues(1) )
      ALLOCATE(CHARACTER(MAX_STRING_LEN)::Name)
      ALLOCATE(CHARACTER(MAX_STRING_LEN)::str)
      ALLOCATE(CHARACTER(MAX_STRING_LEN)::Depname)


      Name = ''
      DO WHILE( ReadAndTrim( InFileUnit,Name,Echo ) )

        IF ( Name == '' .OR. Name == ' ')  CYCLE
        IF ( SEQL(Name,'end') ) THEN
          EXIT
        END IF

        IF ( SEQL(Name,'include') ) THEN
          OPEN( InFileUnit-1,FILE=TRIM(Name(9:)),STATUS='OLD',ERR=10 )
          CALL SectionContents( Model,List,CheckAbort,FreeNames, &
                  Section,InFileUnit-1,ScanOnly, Echo )
          CLOSE( InFileUnit-1 )
          CYCLE

10        CONTINUE

          WRITE( Message, * ) 'Cannot find include file: ', Name(9:)
          CALL Warn( 'Model Input', Message )
          CYCLE
        END IF

        TYPE = LIST_TYPE_CONSTANT_SCALAR
        N1   = 1
        N2   = 1
        SizeGiven = .FALSE.
        DO WHILE( ReadAndTrim(InFileUnit,str,echo,string_literal) ) 

          IF ( string_literal ) THEN
            ReturnType = .TRUE.
            CALL CheckKeyWord( Name, TypeString, CheckAbort,FreeNames,Section,ReturnType )
            IF ( .NOT. ReturnType ) THEN
              CALL SyntaxError( Section, Name,str )
            ELSE IF (TypeString/='string' .AND. TypeString/='file') THEN 
              CALL SyntaxError( Section, Name,str )
            ELSE
              str = TRIM(TypeString) // ' ' // TRIM(str)
            END IF
          END IF

20        CONTINUE

          slen = LEN_TRIM(str)
          j = slen
          DO i=1,slen
            IF ( str(i:i)==' ') EXIT
            j = i
          END DO
          Keyword=str(1:j)
          str_beg = j+2

          SELECT CASE(Keyword)
          CASE('real')
             CALL CheckKeyWord( Name,'real',CheckAbort,FreeNames,Section )

             Proc = 0
             IF ( SEQL(str(str_beg:),'procedure ') ) THEN

               IF ( .NOT. ScanOnly ) THEN
                  Proc = GetProcAddr( str(str_beg+10:) )

                  SELECT CASE( TYPE )
                  CASE( LIST_TYPE_CONSTANT_SCALAR )

                    IF ( SizeGiven ) THEN
                      CALL ListAddConstRealArray( List,Name,N1,N2, &
                             ATx(1:N1,1:N2,1),Proc )
                    ELSE
                      CALL ListAddConstReal( List,Name,Val,Proc )
                    END IF
   
                  CASE( LIST_TYPE_VARIABLE_SCALAR )

                    IF ( SizeGiven ) THEN
                      CALL ListAddDepRealArray( List,Name,Depname,1,ATt, &
                            N1,N2,ATx(1:N1,1:N2,1:1),Proc )
                    ELSE
                      CALL ListAddDepReal( List,Name,Depname,1,ATt,ATx,Proc )
                    END IF
                  END SELECT
               END IF

             ELSE IF ( SEQL(str(str_beg:),'matc ') ) THEN

               IF ( .NOT. ScanOnly ) THEN

                  SELECT CASE( TYPE )
                  CASE( LIST_TYPE_CONSTANT_SCALAR )

                    IF ( SizeGiven ) THEN
                       CALL ListAddConstRealArray( List,Name,N1,N2, &
                         ATx(1:N1,1:N2,1), Proc, str(str_beg+5:) )
                    ELSE
                      CALL ListAddConstReal( List,Name,Val,Proc, &
                                  str(str_beg+5:) )
                    END IF
   
                  CASE( LIST_TYPE_VARIABLE_SCALAR )

                    IF ( SizeGiven ) THEN
                      CALL ListAddDepRealArray( List,Name,Depname,1,ATt, &
                            N1,N2,ATx(1:N1,1:N2,1:n),Proc, str(str_beg+5:) )
                    ELSE
                      CALL ListAddDepReal( List,Name,Depname,1,ATt,ATx, &
                                  Proc, str(str_beg+5:) )
                    END IF
                  END SELECT
               END IF

#ifdef HAVE_LUA
               ! TODO: Here comes the Lua part. Actually create the lua functions here by calling
               !       some routine that transforms str(str_beg+4:) to lua function. But that function needs a name.
             ELSE IF( SEQL(str(str_beg:), 'lua ') ) THEN
               
               IF ( .NOT. ScanOnly ) THEN 
                 SELECT CASE ( TYPE )
                 CASE (LIST_TYPE_CONSTANT_SCALAR )
                   call Fatal('ModelDescription', 'Constant expressions are not supported with Lua. &
                       Please provide at least a dummy argument.')

                   IF ( SizeGiven ) THEN
                     CALL ListAddConstRealArray( List, Name, N1, N2, &
                         ATx(1:N1,1:N2,1), Proc, str(str_beg+4:) )
                   ELSE
                     CALL ListAddConstReal(List, Name, Val, Proc, &
                         str(str_beg+4:))
                   END IF

                 CASE( LIST_TYPE_VARIABLE_SCALAR )
                   block
                     TYPE(ValueListEntry_t), POINTER :: v_ptr
                     CHARACTER(len=:, kind=c_char), pointer :: lua_fname
                     integer :: fname_len, lstat
                     !$OMP PARALLEL default(shared)
                     !$OMP CRITICAL
                     lstat = lua_dostring(LuaState, &
                         'return create_new_fun("'//trim(name)//'", "' // &
                         trim(str(str_beg+4:)) // '")'// c_null_char, 1)
                     lua_fname => lua_popstring(LuaState, fname_len)
                     !$OMP END CRITICAL
                     !$OMP END PARALLEL
                     IF ( SizeGiven ) THEN 
                       CALL ListAddDepRealArray( List, Name, Depname, 1, Att, &
                           N1, N2, Atx(1:N1, 1:N2, 1:n), proc, lua_fname(1:fname_len) // c_null_char)
                     ELSE
                       CALL ListAddDepReal( List, Name, Depname, 1, ATt, ATx, &
                           Proc, lua_fname(1:fname_len) // c_null_char)
                     END IF
                     v_ptr => ListFind(list, name)
                     v_ptr % LuaFun = .true.
                   end block
                   END SELECT

               END IF
#endif
             ELSE

               SELECT CASE( TYPE )
               CASE( LIST_TYPE_CONSTANT_SCALAR )

                  k = 0
                  DO i=1,N1
                     DO j=1,N2
                        DO WHILE( k <= slen )
                           k = k + 1
                           IF ( str(k:k) == ' ' ) EXIT
                        END DO

                        DO WHILE( k <= slen )
                          k = k + 1
                             IF ( str(k:k) /= ' ' ) EXIT
                        END DO

                        IF ( k > slen ) THEN
                          IF ( .NOT. ReadAndTrim( InFileUnit,str,Echo) ) &
                            CALL SyntaxError( Section,Name,str )
                          k = 1
                          slen = LEN_TRIM(str)
                        END IF

                        IF (.NOT.ScanOnly ) THEN
                          READ( str(k:),*,iostat=iostat ) ATx(i,j,1)
                          IF( iostat /= 0 ) THEN
                            CALL Fatal('SectionContents','Problem reading real keyword: '//TRIM(Name)//': '//str(k:)) 
                          END IF
                        END IF
                     END DO
                  END DO
 
                  IF ( .NOT. ScanOnly ) THEN
                     IF ( SizeGiven ) THEN
                       CALL ListAddConstRealArray( List,Name,N1,N2, &
                              ATx(1:N1,1:N2,1) )
                     ELSE
                       CALL ListAddConstReal( List,Name,ATx(1,1,1) )
                     END IF
                  END IF
  
               CASE( LIST_TYPE_VARIABLE_SCALAR )

                 IF (ScanOnly) THEN
                   AbufLen=0
                 ELSE
                   IF (ALLOCATED(ATt) ) DEALLOCATE(ATt,ATx)
                   ALLOCATE( ATt(MaxBufLen), ATx(n1,n2,MaxBufLen) )
                 END IF
                 
                 ! Enable both "cubic monotone" and "monotone cubic"
                 Cubic = SEQL(str(str_beg:),'cubic')
                 IF(Cubic) THEN
                   monotone = SEQL(str(str_beg+6:),'monotone')
                 ELSE
                   monotone = SEQL(str(str_beg:),'monotone')
                   IF( Monotone ) THEN
                     Cubic = SEQL(str(str_beg+9:),'cubic')
                     IF( .NOT. Cubic ) CALL Warn('SectionContents','Monotone curves only applicable to cubic splines!')
                   END IF
                 END IF

                 n = 0
                 DO WHILE( ReadAndTrim(InFileUnit,str,Echo) )

                   IF ( str == '' .OR. str==' '  ) CYCLE
                   IF ( SEQL(str,'end') ) EXIT
 
                   slen = LEN_TRIM(str)
                   IF ( .NOT. ScanOnly ) THEN
                     n = n + 1

                     READ( str,*,iostat=iostat ) ATt(n)
                     IF( iostat /= 0 ) THEN
                       CALL Fatal('SectionContents','Problem reading real keyword: '//TRIM(Name)//': '//str) 
                     END IF

                   ELSE
                     AbufLen = AbufLen+1
                   END IF

                   k = 0
                   DO i=1,N1
                     DO j=1,N2
                       DO WHILE( k <= slen )
                         k = k + 1
                         IF ( str(k:k) == ' ' ) EXIT
                       END DO

                       DO WHILE( k <= slen )
                         k = k + 1
                         IF ( str(k:k) /= ' ' ) EXIT
                       END DO

                       IF ( k > slen ) THEN
                          IF ( .NOT. ReadAndTrim( InFileUnit,str,Echo) ) &
                            CALL SyntaxError( Section,Name,str )

                          k = 1
                          slen = LEN_TRIM(str)
                       END IF

                       IF ( .NOT. ScanOnly ) THEN
                         READ( str(k:),*,iostat=iostat ) ATx(i,j,n)
                         IF( iostat /= 0 ) THEN
                           CALL Fatal('SectionContents','Problem reading real keyword: '//TRIM(Name)//': '//str(k:)) 
                         END IF
                       END IF

                     END DO
                   END DO
                 END DO
 
                 IF (.NOT. ScanOnly ) THEN
                   IF ( SizeGiven ) THEN
                     CALL ListAddDepRealArray( List,Name,Depname,n,ATt(1:n), &
                              N1,N2,ATx(1:N1,1:N2,1:n) )
                   ELSE
                     CALL ListAddDepReal( List,Name,Depname,n,ATt(1:n), &
                         ATx(1,1,1:n),CubicTable=Cubic, Monotone=monotone )
                   END IF
                 END IF
                 MaxBufLen = MAX(MaxBuflen, Abuflen)
               END SELECT
             END IF
             EXIT

          CASE('logical')

            CALL CheckKeyWord( Name,'logical',CheckAbort,FreeNames,Section )

            IF ( .NOT. ScanOnly ) THEN
               IF ( SEQL(str(str_beg:),'true') .OR. &
                 str(str_beg:str_beg) == '1' ) THEN
                 CALL ListAddLogical( List,Name,.TRUE. )
               ELSE IF ( SEQL(str(str_beg:),'false') .OR. &
                 str(str_beg:str_beg) == '0' ) THEN
                 CALL ListAddLogical( List,Name,.FALSE. )
               ELSE 
                 CALL Fatal('SectionContents','Problem reading logical keyword: '//TRIM(Name)//': '//TRIM(str(str_beg:)))
               END IF
            END IF
            EXIT

          CASE('integer')

            CALL CheckKeyWord( Name,'integer',CheckAbort,FreeNames,Section )

             Proc = 0
             IF ( SEQL(str(str_beg:),'procedure ') ) THEN
               IF ( .NOT. ScanOnly ) THEN
                 Proc = GetProcAddr( str(str_beg+10:) )
                 IF ( SizeGiven ) THEN
                   CALL ListAddIntegerArray( List,Name,N1,IValues,Proc )
                 ELSE
                   CALL ListAddInteger( List,Name,k,Proc )
                 END IF
               END IF
             ELSE

               IF ( SizeGiven ) THEN
                 IF ( .NOT. ScanOnly ) THEN
                   IF (SIZE(IValues)<n1) THEN
                     DEALLOCATE(IValues)
                     ALLOCATE(IValues(n1))
                   END IF
                 END IF

                 k = 0
                 DO i=1,N1
                   DO WHILE( k <= slen )
                      k = k + 1
                      IF ( str(k:k) == ' ') EXIT
                   END DO

                   DO WHILE( k <= slen )
                     k = k + 1
                     IF ( str(k:k) /= ' ') EXIT
                   END DO

                   IF ( k > slen ) THEN
                     IF ( .NOT. ReadAndTrim(InFileUnit,str,Echo)) &
                       CALL SyntaxError( Section,Name,str )

                     k = 1
                     slen = LEN_TRIM(str)
                   END IF
                   IF ( .NOT. ScanOnly ) THEN
                     READ( str(k:),*,iostat=iostat ) IValues(i)
                     IF( iostat /= 0 ) THEN
                       CALL Fatal('SectionContents','Problem reading integer keyword: '//TRIM(Name)//': '//str(k:)) 
                     END IF
                   END IF

                 END DO
                 IF ( .NOT. ScanOnly ) CALL ListAddIntegerArray( &
                              List,Name,N1,IValues )
               ELSE IF (.NOT.ScanOnly) THEN
                 READ( str(str_beg:),*,iostat=iostat ) k 
                 IF( iostat /= 0 ) THEN
                   CALL Fatal('SectionContents','Problem reading integer keyword: '//TRIM(Name)//': '//str(str_beg:)) 
                 END IF                 
                 CALL ListAddInteger( List,Name,k )
               END IF
             END IF
             EXIT

          CASE('string')

            CALL CheckKeyWord( Name,'string',CheckAbort,FreeNames,Section )

            IF ( .NOT. ScanOnly ) CALL ListAddString( List,Name,str(str_beg:) )
            EXIT

          CASE('file')

            CALL CheckKeyWord( Name,'file',CheckAbort,FreeNames,Section )

            IF ( .NOT. ScanOnly ) CALL ListAddString( List,Name, &
                          str(str_beg:),.FALSE. )
            EXIT

          CASE('variable')

            DO k=LEN(str),1,-1
              IF ( str(k:k) /= ' ' ) EXIT 
            END DO

            Depname = str(str_beg:k)

            TYPE = LIST_TYPE_VARIABLE_SCALAR

          CASE('equals')

            IF ( .NOT. ScanOnly ) THEN
               DO k=LEN(str),1,-1
                 IF ( str(k:k) /= ' ' ) EXIT 
               END DO

               Depname = str(str_beg:k)

               n=1
               IF ( n > SIZE(ATt) ) THEN
                  DEALLOCATE( ATt, ATx )
                  ALLOCATE( ATt(n), ATx(n1,n2,n) )
               END IF

               ATt(1) = 1.0_dp
               ATx(:,:,1) = 1.0_dp

               IF(n1==1.AND.n2==1) THEN
                 CALL ListAddDepReal(List,Name,Depname,n,ATt,ATx(1,1,1))
               ELSE
                 CALL ListAddDepRealArray(List,Name,Depname,n,ATt,n1,n2,ATx)
               END IF
            END IF
            EXIT

          CASE('opposes')

           IF ( .NOT. ScanOnly ) THEN
               DO k=LEN(str),1,-1
                 IF ( str(k:k) /= ' ' ) EXIT 
               END DO

               Depname = str(str_beg:k)

               n = 1
               IF ( n > SIZE( ATt ) ) THEN
                  DEALLOCATE( ATt, ATx )
                  ALLOCATE( ATt(n), ATx(n1,n2,n) )
               END IF

               ATt(1) = 1.0_dp
               ATx(:,:,1) = -1.0_dp

               IF(n1==1.AND.n2==1) THEN
                 CALL ListAddDepReal(List,Name,Depname,n,ATt,ATx(1,1,1))
               ELSE
                 CALL ListAddDepRealArray(List,Name,Depname,n,ATt,n1,n2,ATx)
               END IF
            END IF
            EXIT

          CASE('size')
            N1 = 1
            N2 = 1
            READ( str(str_beg:),*,err=1,END=1) N1,N2

1           CONTINUE

            IF ( .NOT. ScanOnly ) THEN
               IF ( ALLOCATED( ATx ) ) DEALLOCATE( ATx )
               ALLOCATE( ATx(N1,N2,1) )
               IF ( ALLOCATED( ATt ) ) DEALLOCATE( ATt )
               ALLOCATE( ATt(1) )
            END IF

            SizeGiven = .TRUE.

          CASE('-remove')

            IF ( .NOT. ScanOnly ) CALL ListRemove( List, Name )
            EXIT

          CASE DEFAULT

            ReturnType = .TRUE.
            CALL CheckKeyWord( Name, TypeString, CheckAbort, &
                     FreeNames,Section,ReturnType )
            IF ( ReturnType ) THEN 
              str = TRIM(TypeString) // ' ' // str
              GOTO 20
            END IF
            CALL SyntaxError( Section, Name,str )
          END SELECT
!------------------------------------------------------------------------------
        END DO
!------------------------------------------------------------------------------
      END DO
!------------------------------------------------------------------------------
      END SUBROUTINE SectionContents
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
      SUBROUTINE SyntaxError( Section, Name, LastString )
!------------------------------------------------------------------------------
        CHARACTER(LEN=*) :: Section, Name, LastString

         CALL Error( 'Model Input', ' ' )
         WRITE( Message, * ) 'Unknown specifier: [',TRIM(LastString),']'
         CALL Error( 'Model Input', Message )
         WRITE( Message, * ) 'In section: [', TRIM(Section), ']'
         CALL Error( 'Model Input', Message )
         WRITE( Message, * ) 'For property name:[',TRIM(Name),']'
         CALL Fatal( 'Model Input', Message )
!------------------------------------------------------------------------------
      END SUBROUTINE SyntaxError
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  END SUBROUTINE LoadInputFile
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Loads the Gebhardt factors needed in the radiation heat transfer.
!> This may be done externally only when the emissivities are constant.
!> Hence there is a internal computation of Gebhardt factors that is highly optimized.
!> \deprecated
!------------------------------------------------------------------------------
  SUBROUTINE LoadGebhardtFactors( Mesh,FileName )
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
    CHARACTER(LEN=*) FileName
!------------------------------------------------------------------------------

    INTEGER, ALLOCATABLE :: Mapping(:)
    INTEGER :: i,j,k,l,n,m,p
    REAL(KIND=dp) :: s
    CHARACTER(LEN=MAX_STRING_LEN) :: FName
    TYPE(Element_t), POINTER :: elm,celm

!------------------------------------------------------------------------------

    IF ( LEN_TRIM(Mesh % Name) > 0 ) THEN
      FName = TRIM(OutputPath) // '/' // TRIM(Mesh % Name) // '/' // TRIM(FileName)
    ELSE
      FName = TRIM(FileName)
    END IF
    OPEN( 1,file = TRIM(FName),err=10 )

    CALL Info( 'LoadGebhardtFactors', 'Start', Level=5 )

    READ(1,*) n
    ALLOCATE( mapping(n) )
    DO i=1,n
      READ(1,*) j,mapping(i)
    END DO

    DO i=1,n
      READ(1,*) m
      DO j=1,m
        READ(1,*) k,l,s
        k = mapping(k)
        l = mapping(l)
        IF ( .NOT.ASSOCIATED( &
          mesh % elements(k) % boundaryinfo % gebhardtfactors) ) THEN
          ALLOCATE( mesh % elements(k) % boundaryinfo % gebhardtfactors )
          ALLOCATE(  &
          mesh % elements(k) % boundaryinfo % gebhardtfactors % factors(m), &
          mesh % elements(k) % boundaryinfo % gebhardtfactors % elements(m) )
          mesh % elements(k) % boundaryinfo % gebhardtfactors % numberoffactors = m
        ELSE IF ( mesh % elements(k) % boundaryinfo % gebhardtfactors % numberoffactors/=m ) THEN
          DEALLOCATE(  &
          mesh % elements(k) % boundaryinfo % gebhardtfactors % factors, &
          mesh % elements(k) % boundaryinfo % gebhardtfactors % elements )
          ALLOCATE(  &
          mesh % elements(k) % boundaryinfo % gebhardtfactors % factors(m), &
          mesh % elements(k) % boundaryinfo % gebhardtfactors % elements(m) )
          mesh % elements(k) % boundaryinfo % gebhardtfactors % numberoffactors = m
        END IF
        mesh % elements(k) % boundaryinfo % gebhardtfactors % numberofimplicitfactors = &
            mesh % elements(k) % boundaryinfo % gebhardtfactors % numberoffactors
      END DO
    END DO

    REWIND(1)

    READ(1,*) n

    DO i=1,n
      READ(1,*) j,mapping(i)
    END DO

    DO i=1,n
      READ(1,*) m
      DO j=1,m
        READ(1,*) k,l,s
        k = mapping(k)
        l = mapping(l)
        mesh % elements(k) % boundaryinfo % gebhardtfactors % elements(j) = l
        mesh % elements(k) % boundaryinfo % gebhardtfactors % factors(j)  = s
      END DO
    END DO

    DEALLOCATE(mapping)
    CLOSE(1)

    CALL Info( 'LoadGebhardtFactors', '...Done', Level=5 )

    RETURN

10  CONTINUE

    WRITE( Message, * ) 'Can not open file for GebhardtFactors: ',TRIM(FileName)
    CALL Fatal( 'LoadGebhardtFactors', Message )

  END SUBROUTINE LoadGebhardtFactors
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!>    Figure out the requested coordinate system.
!------------------------------------------------------------------------------
  SUBROUTINE SetCoordinateSystem( Model )
!------------------------------------------------------------------------------
     TYPE(Model_t), POINTER :: Model
!------------------------------------------------------------------------------
     LOGICAL :: Found
     TYPE(Mesh_t), POINTER :: Mesh
     REAL(KIND=dp) :: x,y,z
     CHARACTER(LEN=MAX_NAME_LEN) :: csys
     INTEGER :: Mesh_dim, Model_dim
     
     csys = ListGetString( Model % Simulation, 'Coordinate System', Found )
     IF ( .NOT. Found ) Csys = 'cartesian'

     IF ( csys=='cartesian' .OR. csys=='polar' ) THEN
        Mesh => Model % Meshes

        ! Inherit the maximum dimension from the mesh in case
        ! it is not given.
        Model_dim = 0
        DO WHILE( ASSOCIATED( Mesh ) )          
          Mesh_dim = Mesh % MaxDim
          IF( Mesh_dim == 0 ) THEN
            CALL SetMeshDimension( Mesh )
            Mesh_dim = Mesh % MaxDim
          END IF
          Model_dim = MAX( Model_dim, Mesh_dim )
          IF( Model_dim == 3 ) EXIT
          Mesh => Mesh % Next
        END DO

        Model % Dimension = Model_dim
     END IF

     SELECT CASE ( csys )
       CASE( 'cartesian' )
         Coordinates = Cartesian
       CASE( 'cartesian 1d' )
         Model % DIMENSION = 1
         Coordinates = Cartesian
       CASE( 'cartesian 2d' )
         Model % DIMENSION = 2
         Coordinates = Cartesian
       CASE( 'cartesian 3d' )
         Model % DIMENSION = 3
         Coordinates = Cartesian
       CASE( 'axi symmetric' )
         Model % DIMENSION = 2
         Coordinates = AxisSymmetric
       CASE( 'cylindric symmetric' )
         Model % DIMENSION = 2
         Coordinates = CylindricSymmetric
       CASE( 'cylindrical' )
         Model % DIMENSION = 3
         Coordinates = Cylindric
       CASE( 'polar' )
         Coordinates = Polar
       CASE( 'polar 2d' )
         Model % DIMENSION = 2
         Coordinates = Polar
       CASE( 'polar 3d' )
         Model % DIMENSION = 3
         Coordinates = Polar
       CASE DEFAULT
         WRITE( Message, * ) 'Unknown global coordinate system: ', TRIM(csys)
         CALL Fatal( 'SetCoordinateSystem', Message )
     END SELECT
!------------------------------------------------------------------------------
   END SUBROUTINE SetCoordinateSystem
!------------------------------------------------------------------------------

   
!------------------------------------------------------------------------------
!> Function to read the complete Elmer model: sif file and mesh files.
!------------------------------------------------------------------------------
  FUNCTION LoadModel( ModelName,BoundariesOnly,numprocs,mype,MeshIndex) RESULT( Model )
!------------------------------------------------------------------------------
    USE MeshPartition
    USE SParIterGlobals

    IMPLICIT NONE

    CHARACTER(LEN=*) :: ModelName
    LOGICAL :: BoundariesOnly

    INTEGER, OPTIONAL :: numprocs,mype, MeshIndex
 
    TYPE(Model_t), POINTER :: Model

!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh,Mesh1,NewMesh,OldMesh,SerialMesh
    INTEGER :: i,j,k,l,s,nlen,eqn,MeshKeep,MeshLevels,nprocs
    LOGICAL :: GotIt,GotMesh,found,OneMeshName, OpenFile, Transient
    LOGICAL :: stat, single, MeshGrading
    TYPE(Solver_t), POINTER :: Solver
    INTEGER(KIND=AddrInt) :: InitProc
    INTEGER, TARGET :: Def_Dofs(10,6)
    REAL(KIND=dp) :: MeshPower
    REAL(KIND=dp), POINTER :: h(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: Name,ElementDef,ElementDef0
    CHARACTER(LEN=MAX_NAME_LEN) :: MeshDir,MeshName
    TYPE(valuelist_t), POINTER :: lst
    INTEGER, ALLOCATABLE :: EdgeDOFs(:),FaceDOFs(:)
!------------------------------------------------------------------------------

    ALLOCATE( Model )
    CurrentModel => Model
    NULLIFY( Model % Variables )

    MeshDir  = ' '
    MeshName = ' '

    Model % DIMENSION = 0
    Model % NumberOfBoundaries = 0
    Model % NumberOfBodies     = 0
    Model % NumberOfICs        = 0
    Model % NumberOfBCs        = 0
    Model % NumberOfEquations  = 0
    Model % NumberOfSolvers    = 0
    Model % NumberOfMaterials  = 0
    Model % NumberOfBodyForces = 0

#ifdef HAVE_LUA
    BLOCK
      INTEGER :: lstat, ompthread
      CHARACTER(LEN=256) :: txcmd

#if USE_ISO_C_BINDINGS
      character(len=256) :: elmer_home_env
      CALL getenv("ELMER_HOME", elmer_home_env)
#endif

      !$OMP PARALLEL Shared(parenv, ModelName, elmer_home_env) Private(txcmd, ompthread, lstat) Default(none)
      !$OMP CRITICAL
      LuaState = lua_init()
      IF(.NOT. LuaState % Initialized) THEN
        CALL Fatal('ModelDescription', 'Failed to initialize Lua subsystem.')
      END IF

      ! Store mpi task and omp thread ids in a table
      LSTAT = lua_dostring(LuaState, 'ELMER_PARALLEL = {}' // c_null_char)
      write(txcmd,'(A,I0)') 'ELMER_PARALLEL["pe"] = ', parenv % mype
      lstat = lua_dostring(LuaState, txcmd // c_null_char)

      ompthread = 1
      !$ ompthread = omp_get_thread_num()
      WRITE(txcmd,'(A,I0)') 'ELMER_PARALLEL["thread"] = ', ompthread
      lstat = lua_dostring(LuaState, txcmd // c_null_char)
      
      WRITE(txcmd,'(A,I0, A)') 'tx = array.new(', MAX_FNC, ')'

      ! TODO: (2018-09-17) Nowadays ISO_C_BINDINGS are pretty much mandatory to compile elmer
#if USE_ISO_C_BINDINGS
      ! Call defaults.lua using 1) ELMER_HOME environment variable or 2) ELMER_SOLVER_HOME preprocessor macro
      ! TODO: (2018-09-18) ELMER_SOLVER_HOME might be too long

      if (trim(elmer_home_env) == "") then
        lstat = lua_dostring(LuaState, &
            'loadfile("' // &
ELMER_SOLVER_HOME &
                    // '" .. "/lua-scripts/defaults.lua")()'//c_null_char)
      else
#endif
        lstat = lua_dostring(LuaState, &
            'loadfile(os.getenv("ELMER_HOME") .. "/share/elmersolver/lua-scripts/defaults.lua")()'//c_null_char)
#if USE_ISO_C_BINDINGS
      end if
#endif

      ! Execute lua parts 
      lstat = lua_dostring(LuaState, 'loadstring(readsif("'//trim(ModelName)//'"))()' // c_null_char)
      lstat = lua_dostring(LuaState,  trim(txcmd)// c_null_char)
      LuaState % tx => lua_getusertable(LuaState, 'tx'//c_null_char)
      !$OMP END CRITICAL
      !$OMP END PARALLEL
    END BLOCK
#endif

    INQUIRE( Unit=InFileUnit, OPENED=OpenFile )
    IF ( .NOT. OpenFile ) OPEN( Unit=InFileUnit, File=Modelname, STATUS='OLD' )
    CALL LoadInputFile( Model,InFileUnit,ModelName,MeshDir,MeshName, .TRUE., .TRUE. )
    REWIND( InFileUnit )
    CALL LoadInputFile( Model,InFileUnit,ModelName,MeshDir,MeshName, .TRUE., .FALSE. )
    IF ( .NOT. OpenFile ) CLOSE( InFileUnit )


    IF( .NOT. ListCheckPresent( Model % Simulation,'Solver Input File') ) THEN
      CALL ListAddString( Model % Simulation,'Solver Input File',ModelName ) 
    END IF
      
    
    CALL InitializeOutputLevel()

    Transient=ListGetString(Model % Simulation, &
        'Simulation Type',Found)=='transient'

    Def_Dofs = -1; Def_Dofs(:,1)=1

    i = 1
    DO WHILE(i<=Model % NumberOFSolvers)

      Solver => Model % Solvers(i)
      Model % Solver => Solver

      Name = ListGetString( Solver  % Values, 'Procedure', Found )
      IF ( Found ) THEN
        InitProc = GetProcAddr( TRIM(Name)//'_Init0', abort=.FALSE. )
        IF ( InitProc /= 0 ) THEN
          CALL ExecSolver( InitProc, Model, Solver, &
                  Solver % dt, Transient )

          Solver => Model % Solvers(i)
          Model % Solver => Solver
        END IF
      END IF

      GotMesh = ListCheckPresent(Solver % Values, 'Mesh')

      IF(.NOT.ALLOCATED(Solver % Def_Dofs)) THEN
        ALLOCATE(Solver % Def_Dofs(10,Model % NumberOfBodies,6))
        Solver % Def_Dofs = -1; Solver % Def_Dofs(:,:,1)=1
      END IF

      ! Define what kind of element we are working with in this solver
      !-----------------------------------------------------------------
      ElementDef = ListGetString( Solver % Values, 'Element', stat )
   
      IF ( .NOT. stat ) THEN
        IF ( ListGetLogical( Solver % Values, 'Discontinuous Galerkin', stat ) ) THEN
           Solver % Def_Dofs(:,:,4) = 0
           IF ( .NOT. GotMesh ) Def_Dofs(:,4) = MAX(Def_Dofs(:,4),0 )
           i=i+1
           CYCLE
        ELSE
           ElementDef = "n:1"
        END IF
      END IF

      ElementDef0 = ElementDef
      DO WHILE(.TRUE.)
        j = INDEX( ElementDef0, '-' )
        IF (j>0) THEN
          ElementDef = ElementDef0(1:j-1)
        ELSE
          ElementDef = ElementDef0
        END IF
        CALL GetDefs( ElementDef )
        IF(j>0) THEN
          ElementDef0 = ElementDef0(j+1:)
         ELSE
          EXIT
        END IF
      END DO

      i = i + 1
    END DO

    ! Check the mesh 
    !--------------------------------------------------------
    Name = ListGetString( Model % Simulation, 'Mesh', GotIt )
    IF(PRESENT(MeshIndex)) THEN
      IF ( MeshIndex>0 )Name = TRIM(Name)//TRIM(I2S(MeshIndex))
    END IF

    OneMeshName = .FALSE.
    IF ( GotIt ) THEN
      k = 1
      i = 1
      nlen = LEN_TRIM(name)
      DO WHILE( k<=nlen .AND. name(k:k) /= ' ' )
        MeshDir(i:i)  = name(k:k)
        Meshname(i:i) = name(k:k)
        k = k + 1
        i = i + 1
      END DO

      DO WHILE( k<=nlen .AND. Name(k:k) == ' ' )
        k = k + 1
      END DO

      IF ( k<=nlen ) THEN
         MeshName(i:i) = '/'
         i = i + 1
         DO WHILE( name(k:k) /= ' ' )
           MeshName(i:i) = Name(k:k)
           k = k + 1
           i = i + 1
         END DO
      ELSE
         OneMeshName = .TRUE.
         MeshDir = "." // CHAR(0)
      END IF
      MeshName(i:i) = CHAR(0)
    ELSE
      IF(PRESENT(MeshIndex)) THEN
        IF(MeshIndex>0) MeshName = MeshName(1:LEN_TRIM(MeshName)-1) // TRIM(I2S(MeshIndex))//CHAR(0)
      END IF
    END IF

    NULLIFY( Model % Meshes )
    IF ( MeshDir(1:1) /= ' ' ) THEN

      CALL ResetTimer('LoadMesh') 

      Single = ListGetLogical( Model % Simulation,'Partition Mesh', GotIt ) 
      IF ( Single ) THEN
        IF( ParEnv % PEs == 1 ) THEN
          CALL Warn('LoadMesh','Why perform partitioning in serial case?')
        END IF
        IF( ParEnv % MyPe == 0 ) THEN
          SerialMesh => LoadMesh2( Model,MeshDir,MeshName,BoundariesOnly,&
              1,0,def_dofs,LoadOnly = .TRUE. )
          CALL PartitionMeshSerial( Model, SerialMesh, Model % Simulation )
        ELSE
          SerialMesh => AllocateMesh()
        END IF

        IF( ParEnv % PEs > 1) THEN
          Model % Meshes => ReDistributeMesh( Model, SerialMesh, .FALSE., .TRUE. )
        ELSE
          CALL Info('LoadMesh','Only one active partition, using the serial mesh as it is!')
          IF( MAXVAL( SerialMesh % RePartition ) <= 1 ) THEN
            DEALLOCATE( SerialMesh % RePartition ) 
          END IF
          Model % Meshes => SerialMesh
        END IF
          
        CALL PrepareMesh( Model, Model % Meshes, ParEnv % PEs > 1, Def_Dofs )          
      ELSE
        Model % Meshes => LoadMesh2( Model, MeshDir, MeshName, &
            BoundariesOnly, numprocs, mype, Def_Dofs )
      END IF
      

      IF(.NOT.ASSOCIATED(Model % Meshes)) THEN
        CALL FreeModel(Model)
        Model => NULL()
        RETURN
      END IF

      CALL CheckTimer('LoadMesh',Level=5,Delete=.TRUE.)

      CALL SetCoordinateSystem( Model )


      MeshLevels = ListGetInteger( Model % Simulation, 'Mesh Levels', GotIt )
      IF ( .NOT. GotIt ) MeshLevels=1

      IF( MeshLevels > 1 ) THEN
        CALL Info('LoadModel','Creating hierarchy of meshes by mesh multiplication: '&
            //TRIM(I2S(MeshLevels)))
      END IF
      MeshKeep = ListGetInteger( Model % Simulation, 'Mesh keep',  GotIt )
      IF ( .NOT. GotIt ) MeshKeep=MeshLevels

      MeshPower   = ListGetConstReal( Model % Simulation, 'Mesh Grading Power',GotIt)
      MeshGrading = ListGetLogical( Model % Simulation, 'Mesh Keep Grading', GotIt)

      DO i=2,MeshLevels
        OldMesh => Model % Meshes

        IF (MeshGrading) THEN
          ALLOCATE(h(OldMesh % NumberOfNodes))
          Model % Mesh => OldMesh
          CALL GetNodalElementSize(Model,MeshPower,.FALSE.,h)
          NewMesh => SplitMeshEqual(OldMesh,h)
          DEALLOCATE(h)
        ELSE
          NewMesh => SplitMeshEqual(OldMesh)
        END IF

        IF(ASSOCIATED(OldMesh % Faces)) THEN
          CALL FindMeshEdges(NewMesh)

          ALLOCATE( EdgeDOFs(NewMesh % NumberOfBulkElements))
          ALLOCATE( FaceDOFs(NewMesh % NumberOfBulkElements))
          EdgeDOFs = MAX(0,MAXVAL(Def_Dofs(:,2)))
          FaceDOFs = MAX(0,MAXVAL(Def_Dofs(:,3)))
          CALL SetMeshEdgeFaceDofs(NewMesh,EdgeDOFs,FaceDOFs)
          DEALLOCATE(EdgeDOFs,FaceDOFs)

          CALL SetMeshMaxDofs(NewMesh)
          IF(ParEnv % PEs>1) CALL SParEdgeNumbering(NewMesh)
          IF(ParEnv % PEs>1) CALL SParFaceNumbering(NewMesh)
        END IF

        IF ( i>MeshLevels-MeshKeep+1 ) THEN
          NewMesh % Next => OldMesh
          NewMesh % Parent => OldMesh
          OldMesh % Child  => NewMesh
          Newmesh % OutputActive = .TRUE.
          OldMesh % OutputActive = .FALSE.
        ELSE
          CALL ReleaseMesh(OldMesh)
        END IF
        Model % Meshes => NewMesh
      END DO


      IF ( OneMeshName ) THEN
         i = 0
      ELSE
         i = LEN_TRIM(MeshName)
         DO WHILE( i>0 )
           IF (MeshName(i:i) == '/') EXIT 
           i = i-1
         END DO
      END IF

      i = i + 1
      k = 1
      Model % Meshes % Name = ' '
      DO WHILE( MeshName(i:i) /= CHAR(0) )
        Model % Meshes % Name(k:k) = MeshName(i:i)
        k = k + 1
        i = i + 1
      END DO

      ! Ok, give name also to the parent meshes as they might be saved too
      OldMesh => Model % Meshes % Parent
      DO WHILE( ASSOCIATED( OldMesh ) )
        OldMesh % Name = TRIM(OldMesh % Child % Name)//'p'
        OldMesh => OldMesh % Parent
      END DO

      DO i=1,Model % NumberOfSolvers
         Model % Solvers(i) % Mesh => Model % Meshes
      END DO
    END IF

    DO s=1,Model % NumberOfSolvers
      Name = ListGetString( Model % Solvers(s) % Values, 'Mesh', GotIt )
      IF(PRESENT(MeshIndex)) THEN
        IF ( MeshIndex>0 )Name = TRIM(Name)//TRIM(I2S(MeshIndex))
      END IF

      IF( GotIt ) THEN
        WRITE(Message,'(A,I0)') 'Loading solver specific mesh > '//TRIM(Name)// ' < for solver ',s
        CALL Info('LoadModel',Message,Level=7)

        single=.FALSE.
        IF ( SEQL(Name, '-single ') ) THEN
          single=.TRUE.
          Name=Name(9:)
        END IF

        nprocs = numprocs
        IF ( SEQL(Name, '-part ') ) THEN
          READ( Name(7:), * ) nprocs
          i = 7
          DO WHILE(Name(i:i)/=' ')
           i=i+1
          END DO
          Name=Name(i+1:)
        END IF


        OneMeshName = .FALSE.
        k = 1
        i = 1
        nlen = LEN_TRIM(name)
        MeshName = ' '
        DO WHILE( k<=nlen .AND. name(k:k) /= ' ' )
          MeshDir(i:i)  = name(k:k)
          Meshname(i:i) = name(k:k)
          k = k + 1
          i = i + 1
        END DO

        DO WHILE( k<=nlen .AND. Name(k:k) == ' ' )
          k = k + 1
        END DO

        IF ( k<=nlen ) THEN
          MeshName(i:i) = '/'
          i = i + 1
          DO WHILE( name(k:k) /= ' ' )
            MeshName(i:i) = Name(k:k)
            k = k + 1
            i = i + 1
          END DO
        ELSE
          OneMeshName = .TRUE.
          MeshDir = "." // CHAR(0)
        END IF
        MeshName(i:i) = CHAR(0)

        IF ( OneMeshName ) THEN
          i = 0
        ELSE
          DO WHILE( i>0 .AND. MeshName(i:i) /= '/' )
            i = i - 1
          END DO
        END IF

        Mesh => Model % Meshes
        Found = .FALSE.
        DO WHILE( ASSOCIATED( Mesh ) )
           Found = .TRUE.
           k = 1
           j = i+1
           DO WHILE( MeshName(j:j) /= CHAR(0) )
              IF ( Mesh % Name(k:k) /= MeshName(j:j) ) THEN
                Found = .FALSE.
                EXIT
              END IF
              k = k + 1
              j = j + 1
           END DO
           IF ( LEN_TRIM(Mesh % Name) /= k-1 ) Found = .FALSE.
           IF ( Found ) EXIT
           Mesh => Mesh % Next
        END DO

        IF ( Found ) THEN
          Model % Solvers(s) % Mesh => Mesh
          CYCLE
        END IF

        DO i=1,6
          DO j=1,8
            Def_Dofs(j,i) = MAXVAL(Model % Solvers(s) % Def_Dofs(j,:,i))
          END DO
        END DO


        IF ( Single ) THEN
          Model % Solvers(s) % Mesh => &
              LoadMesh2( Model,MeshDir,MeshName,BoundariesOnly,1,0,def_dofs, s )
        ELSE
          IF ( mype < nprocs ) THEN
            Model % Solvers(s) % Mesh => &
                LoadMesh2( Model,MeshDir,MeshName,BoundariesOnly,nprocs,mype,Def_Dofs, s )
          ELSE
            ! There are more partitions than partitions in mesh, just allocate
            Model % Solvers(s) % Mesh => AllocateMesh()
          END IF
        END IF
        Model % Solvers(s) % Mesh % OutputActive = .TRUE.


        MeshLevels = ListGetInteger( Model % Solvers(s) % Values, 'Mesh Levels', GotIt )
        IF ( .NOT. GotIt ) MeshLevels=1

        MeshKeep = ListGetInteger( Model % Solvers(s) % Values, 'Mesh keep',  GotIt )
        IF ( .NOT. GotIt ) MeshKeep=MeshLevels

        MeshPower   = ListGetConstReal( Model % Simulation, 'Mesh Grading Power',GotIt)
        MeshGrading = ListGetLogical( Model % Simulation, 'Mesh Keep Grading', GotIt)


        DO i=2,MeshLevels
          OldMesh => Model % Solvers(s) % Mesh

          IF (MeshGrading) THEN
            ALLOCATE(h(OldMesh % NumberOfNodes))
            Model % Mesh => OldMesh
            CALL GetNodalElementSize(Model,MeshPower,.FALSE.,h)
            NewMesh => SplitMeshEqual(OldMesh,h)
            DEALLOCATE(h)
          ELSE
            NewMesh => SplitMeshEqual(OldMesh)
          END IF

          IF(ASSOCIATED(OldMesh % Faces)) THEN
            CALL FindMeshEdges(NewMesh)

            ALLOCATE( EdgeDOFs(NewMesh % NumberOfBulkElements))
            ALLOCATE( FaceDOFs(NewMesh % NumberOfBulkElements))
            EdgeDOFs = MAX(0,MAXVAL(Def_Dofs(:,2)))
            FaceDOFs = MAX(0,MAXVAL(Def_Dofs(:,3)))
            CALL SetMeshEdgeFaceDofs(NewMesh,EdgeDOFs,FaceDOFs)
            DEALLOCATE(EdgeDOFs,FaceDOFs)

            CALL SetMeshMaxDofs(NewMesh)
            IF(ParEnv % PEs>1) CALL SParEdgeNumbering(NewMesh)
            IF(ParEnv % PEs>1) CALL SParFAceNumbering(NewMesh)
          END IF

          IF ( i>MeshLevels-MeshKeep+1 ) THEN
            NewMesh % Next => OldMesh
            NewMesh % Parent => OldMesh
            OldMesh % Child  => NewMesh
            NewMesh % Name = OldMesh % Name
            Newmesh % OutputActive = .TRUE.
            OldMesh % OutputActive = .FALSE.
          ELSE
            CALL ReleaseMesh(OldMesh)
          END IF
          Model % Solvers(s) % Mesh => NewMesh
        END DO

        IF ( OneMeshName ) i = 0

        k = 1
        i = i + 1
        Model % Solvers(s) % Mesh % Name = ' '
        DO WHILE( MeshName(i:i) /= CHAR(0) )
          Model % Solvers(s) % Mesh % Name(k:k) = MeshName(i:i)
          k = k + 1
          i = i + 1
        END DO

        IF ( ASSOCIATED( Model % Meshes ) ) THEN
          Mesh1 => Model % Meshes
          DO WHILE( ASSOCIATED( Mesh1 % Next ) ) 
            Mesh1 => Mesh1 % Next
          END DO
          Mesh1 % Next => Model % Solvers(s) % Mesh
        ELSE
          Model % Meshes => Model % Solvers(s) % Mesh
        END IF
      END IF
    END DO

    CALL SetCoordinateSystem( Model )
  
    IF ( OutputPath == ' ' ) THEN
      DO i=1,MAX_NAME_LEN
        IF ( MeshDir(i:i) == CHAR(0) ) EXIT
        OutputPath(i:i) = MeshDir(i:i)
      END DO
    END IF

    Mesh => Model % Meshes
    DO WHILE( ASSOCIATED( Mesh ) )
       CALL MeshStabParams( Mesh )
       Mesh => Mesh % Next
    END DO

!------------------------------------------------------------------------------

  CONTAINS

!------------------------------------------------------------------------------
    SUBROUTINE GetDefs(ElementDef)
!------------------------------------------------------------------------------
      CHARACTER(LEN=*) :: ElementDef
!------------------------------------------------------------------------------
      INTEGER  :: ind(8),i,j,n
      INTEGER, POINTER :: gdofs(:,:), sdofs(:,:,:)

      ind = [1,2,3,4,5,6,7,8]

      IF (SEQL(ElementDef, 'point') )     ind=1
      IF (SEQL(ElementDef, 'line') )      ind=2
      IF (SEQL(ElementDef, 'tri') )       ind=3
      IF (SEQL(ElementDef, 'quad') )      ind=4
      IF (SEQL(ElementDef, 'tetra') )     ind=5
      IF (SEQL(ElementDef, 'pyramid') )   ind=6
      IF (SEQL(ElementDef, 'prism') )     ind=7
      IF (SEQL(ElementDef, 'brick') )     ind=8
      IF (SEQL(ElementDef, 'tri_face') )  ind=9
      IF (SEQL(ElementDef, 'quad_face') ) ind=10

      n = INDEX(ElementDef,'-')
      IF (n<=0) n=LEN_TRIM(ElementDef)
          
      j = INDEX( ElementDef(1:n), 'n:' )
      IF ( j>0 ) THEN
        READ( ElementDef(j+2:), * ) l
        Solver % Def_Dofs(ind,:,1) = l
        IF (.NOT. GotMesh ) Def_Dofs(ind,1) = MAX(Def_Dofs(ind,1), l)
      END IF
          
      j = INDEX( ElementDef(1:n), 'e:' )
      IF ( j>0 ) THEN
        READ( ElementDef(j+2:), * ) l
        Solver % Def_Dofs(ind,:,2) = l
        IF ( .NOT. GotMesh ) Def_Dofs(ind,2) = MAX(Def_Dofs(ind,2), l )
      END IF
          
      j = INDEX( ElementDef(1:n), 'f:' )
      IF ( j>0 ) THEN
        READ( ElementDef(j+2:), * ) l
        Solver % Def_Dofs(ind,:,3) = l
        IF ( .NOT. GotMesh ) Def_Dofs(ind,3) = MAX(Def_Dofs(ind,3), l )
      END IF
          
      j = INDEX( ElementDef(1:n), 'd:' )
      IF ( j>0 ) THEN
        READ( ElementDef(j+2:), * ) l
        Solver % Def_Dofs(ind,:,4) = l
        IF ( .NOT. GotMesh ) Def_Dofs(ind,4) = MAX(Def_Dofs(ind,4), l )
      ELSE 
        IF ( ListGetLogical( Solver % Values, &
            'Discontinuous Galerkin', stat ) ) THEN
          Solver % Def_Dofs(ind,:,4) = 0
          IF ( .NOT. GotMesh ) Def_Dofs(ind,4) = MAX(Def_Dofs(ind,4),0 )
        END IF
      END IF
          
      j = INDEX( ElementDef(1:n), 'b:' )
      IF ( j>0 ) THEN
        READ( ElementDef(j+2:), * ) l
        Solver % Def_Dofs(ind,:,5) = l
        IF ( .NOT. GotMesh ) Def_Dofs(ind,5) = MAX(Def_Dofs(ind,5), l )
      END IF
          
      j = INDEX( ElementDef(1:n), 'p:' )
      IF ( j>0 ) THEN
        IF ( ElementDef(j+2:j+2)=='%' ) THEN
          Solver % Def_Dofs(ind,:,6) = 0
        ELSE
          READ( ElementDef(j+2:), * ) l
          Solver % Def_Dofs(ind,:,6) = l
          IF ( .NOT. GotMesh ) Def_Dofs(ind,6) = MAX(Def_Dofs(ind,6), l )
         END IF
      END IF

!------------------------------------------------------------------------------
    END SUBROUTINE GetDefs
!------------------------------------------------------------------------------

    
    !------------------------------------------------------------------------------
    !      Initialize the log file output system for Messages
    !------------------------------------------------------------------------------
    SUBROUTINE InitializeOutputLevel( )
      
      INTEGER, POINTER :: OutputMask(:)
    
       
      MinOutputLevel = ListGetInteger( CurrentModel % Simulation, &
          'Min Output Level', GotIt )
      
      MaxOutputLevel = ListGetInteger( CurrentModel % Simulation, &
          'Max Output Level', GotIt )
      
      IF ( .NOT. GotIt ) MaxOutputLevel = 10
      
      OutputMask => ListGetIntegerArray( CurrentModel % Simulation, &
          'Output Level', GotIt )
      
      IF ( GotIt ) THEN
        DO i=1,SIZE(OutputMask)
          OutputLevelMask(i-1) = OutputMask(i) /= 0
        END DO
      END IF
      
      DO i=0,31
        OutputLevelMask(i) = OutputLevelMask(i) .AND. &
            i >= MinOutputLevel .AND. i <= MaxOutputLevel
      END DO
      
      OutputPrefix = ListGetLogical( CurrentModel % Simulation, &
          'Output Prefix', GotIt )
      IF ( .NOT. GotIt ) OutputPrefix = .FALSE.
      
      OutputCaller = ListGetLogical( CurrentModel % Simulation, &
          'Output Caller', GotIt )
      IF ( .NOT. GotIt ) OutputCaller = .TRUE.
      
      ! By default only on partition is used to show the results
      ! For debugging it may be useful to show several.
      MaxOutputPE = ListGetInteger( CurrentModel % Simulation, &
          'Max Output Partition', GotIt )    
      IF( GotIt ) THEN
        MaxOutputPE = MIN(ParEnv % PEs, MaxOutputPE)        
        MinOutputPE = ListGetInteger( CurrentModel % Simulation, &
            'Min Output Partition', GotIt )    
        MinOutputPE = MAX(0, MinOutputPE)
        
        IF( ParEnv % MyPe >= MinOutputPE .AND. &
            ParEnv % MyPe <= MaxOutputPE ) THEN 
          OutputPE = ParEnv % MyPE
        ELSE
          OutputPE = -1
        END IF
      END IF 
                    
    END SUBROUTINE InitializeOutputLevel
!------------------------------------------------------------------------------
  END FUNCTION LoadModel
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Some keywords automatically require other keywords to be set
!> We could complain on the missing keywords later on, but sometimes 
!> it may be just as simple to add them directly. 
!------------------------------------------------------------------------------
  SUBROUTINE CompleteModelKeywords()

    TYPE(Model_t), POINTER :: Model 
    TYPE(ValueList_t), POINTER :: List, ListB
    INTEGER :: i,j,k,n,nb
    LOGICAL :: Found, Flag
    CHARACTER(LEN=MAX_NAME_LEN) :: Name, NameB
    REAL(KIND=dp) :: Tol = 1.0e-8
    
    CALL Info('CompleteModelKeywords','Completing keywords for mortar BCs',Level=12)

    Model => CurrentModel 

    IF( ListGetLogical( Model % Simulation,'Mortar BCs Rotational',Found ) ) THEN     
      Tol = ListGetConstReal( Model % Simulation,&
          'Mortar BCs Rotational Tolerance',Found )
      IF(.NOT. Found ) Tol = 1.0e-6
      CALL DetectMortarPairs( Model, Model % Meshes, Tol, 4, .TRUE. )
    END IF
    IF( ListGetLogical( Model % Simulation,'Mortar BCs Radial',Found ) ) THEN
      Tol = ListGetConstReal( Model % Simulation,&
          'Mortar BCs Radial Tolerance',Found )
      IF(.NOT. Found ) Tol = 1.0e-3
      CALL DetectMortarPairs( Model, Model % Meshes, Tol, 5, .FALSE. )     
    END IF
    IF( ListGetLogical( Model % Simulation,'Mortar BCs Axial',Found ) ) THEN
      Tol = ListGetConstReal( Model % Simulation,&
          'Mortar BCs Axial Tolerance',Found )
      IF(.NOT. Found ) Tol = 1.0e-6
      CALL DetectMortarPairs( Model, Model % Meshes, Tol, 3, .TRUE. )           
    END IF
      
    
    IF( ListGetLogical( Model % Simulation,'Use Mortar Names',Found ) ) THEN
      DO i=1,Model % NumberOfBCs
        List => Model % BCs(i) % Values
        Name = ListGetString( List,'Name',Found )
        IF(.NOT. Found ) CYCLE
        n = INDEX(Name,'_mortar') - 1
        IF( n <= 0 ) CYCLE
        
        DO j=1,Model % NumberOfBCs
          ListB => Model % BCs(j) % Values
          NameB = ListGetString( ListB,'Name',Found )
          
          IF(.NOT. Found ) CYCLE
          IF( ListCheckPresent( List,'Mortar BC') ) CYCLE

          nb = LEN_TRIM(NameB)
          IF( nb /= n ) CYCLE          
          IF( Name(1:n) == NameB(1:n) ) THEN
            CALL Info('CompleteModelKeywords','Adding > Mortar BC = '&
                //TRIM(I2S(i))//' < to boundary '//TRIM(I2S(j)),Level=5)
            CALL ListAddInteger( ListB,'Mortar BC',i )
            EXIT
          END IF
        END DO
      END DO
    END IF
      
    IF( ListGetLogical( Model % Simulation,'Use Contact Names',Found ) ) THEN
      DO i=1,Model % NumberOfBCs
        List => Model % BCs(i) % Values
        Name = ListGetString( List,'Name',Found )
        IF(.NOT. Found ) CYCLE
        n = INDEX(Name,'_contact') - 1
        IF( n <= 0 ) CYCLE

        DO j=1,Model % NumberOfBCs
          ListB => Model % BCs(j) % Values
          NameB = ListGetString( ListB,'Name',Found )
          
          IF(.NOT. Found ) CYCLE
          IF( ListCheckPresent( List,'Contact BC') ) CYCLE

          nb = LEN_TRIM(NameB)
          IF( nb /= n ) CYCLE
          IF( Name(1:n) == NameB(1:n) ) THEN
            CALL Info('CompleteModelKeywords','Adding > Contact BC = '&
                //TRIM(I2S(i))//' < to boundary '//TRIM(I2S(j)),Level=5)
            CALL ListAddInteger( ListB,'Contact BC',i )
            EXIT
          END IF
        END DO
      END DO
    END IF
      

    DO i=1,Model % NumberOfBCs
      List => Model % BCs(i) % Values
      j = ListGetInteger( List,'Mortar BC',Found )
      IF(j==0) j = ListGetInteger( List,'Contact BC',Found )
      IF(j==0) CYCLE

      IF( j > Model % NumberOfBCs ) CYCLE

      ListB => Model % BCs(j) % Values
      IF(.NOT. ASSOCIATED( ListB ) ) CYCLE

      CALL ListCompareAndCopy( List, ListB,'Mass Consistent Normals',Found )
      IF( Found ) CALL Info('CompleteModelKeywords',&
          'Added > Mass Consistent Normals < to master BC '//TRIM(I2S(j)),Level=10)

      CALL ListCompareAndCopy( List, ListB,'Rotational Normals',Found )
      IF( Found ) CALL Info('CompleteModelKeywords',&
          'Added > Rotational Normals < to master BC '//TRIM(I2S(j)),Level=10)

      CALL ListCompareAndCopy( List, ListB,'Normal-Tangential Displacement',Found )
      IF( Found ) CALL Info('CompleteModelKeywords',&
          'Added > Normal-Tangential Displacement < to master BC '//TRIM(I2S(j)),Level=10)

      CALL ListCompareAndCopy( List, ListB,'Normal-Tangential Velocity',Found )
      IF( Found ) CALL Info('CompleteModelKeywords',&
          'Added > Normal-Tangential Velocity < to master BC '//TRIM(I2S(j)),Level=10)
    END DO


  END SUBROUTINE CompleteModelKeywords
  


!------------------------------------------------------------------------------
!> Save fields in a file that may be used in restarting of Elmer simulation.
!> The data may be saved either in ascii and binary formats.
!------------------------------------------------------------------------------
  FUNCTION SaveResult( Filename,Mesh,Time,SimulationTime,Binary,SaveAll,&
                       FreeSurface ) RESULT(SaveCount)
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER :: Time,SaveCount
    CHARACTER(LEN=*) :: Filename
    REAL(KIND=dp) :: SimulationTime
    LOGICAL :: Binary,SaveAll
    LOGICAL, OPTIONAL :: FreeSurface

!------------------------------------------------------------------------------

    TYPE(Element_t), POINTER :: CurrentElement
    INTEGER :: i,j,k,k2,DOFs, dates(8), n, PermSize,IsVector
    TYPE(Variable_t), POINTER :: Var
    CHARACTER(LEN=MAX_NAME_LEN) :: FName, PosName, DateStr, EqName, VarName
    LOGICAL :: SaveCoordinates, MoveBoundary, GotIt, SaveThis, &
        SaveGlobal, OutputVariableList
    INTEGER, POINTER :: PrevPerm(:) 
    INTEGER(IntOff_k) :: PrevPermPos, Pos
    INTEGER(IntOff_k), SAVE :: VarPos(MAX_OUTPUT_VARS) = 0
    LOGICAL :: Found
    CHARACTER(1) :: E

    SAVE SaveCoordinates

!------------------------------------------------------------------------------
!   If first time here, count number of variables
!------------------------------------------------------------------------------

    FName = FileName
    IF ( .NOT. FileNameQualified(FileName) ) THEN
      IF ( LEN_TRIM(OutputPath) > 0 ) THEN
        FName = TRIM(OutputPath) // '/' // TRIM(FileName)
      END IF
    END IF
    IF( Mesh % SavesDone == 0) THEN
      CALL Info('SaveResult','-----------------------------------------',Level=4)
      CALL Info('SaveResult','Saving results to file: '//TRIM(FName), Level=4 )
    END IF

    PosName = TRIM(FName) // ".pos"

    SaveGlobal = ListGetLogical( CurrentModel % Simulation,'Output Global Variables',&
        Found ) 
    
    OutputVariableList = ListCheckPresent( CurrentModel % Simulation,&
        'Output Variable 1')

    ! The first time we start by writing the header.
    IF ( Mesh % SavesDone == 0 ) THEN

      ! Check whether the coordinates should be saved also 
      IF ( PRESENT( FreeSurface ) ) THEN
        SaveCoordinates = FreeSurface
      ELSE
        SaveCoordinates = .FALSE.
        MoveBoundary    = .FALSE.
        DO i=1,CurrentModel % NumberOfBCs
          SaveCoordinates = ListGetLogical( &
              CurrentModel % BCs(i) % Values,'Free Surface', GotIt )
          IF ( SaveCoordinates ) THEN
            MoveBoundary =  ListGetLogical( &
                CurrentModel % BCs(i) % Values,'Internal Move Boundary', GotIt )         
            IF ( GotIt ) SaveCoordinates = MoveBoundary
          END IF          
          IF ( SaveCoordinates ) EXIT
        END DO
      END IF      
      IF( ListGetLogical( CurrentModel % Simulation,'Output Coordinates',&
          Found ) ) SaveCoordinates = .TRUE.      
      IF( SaveCoordinates ) THEN
        IF ( Mesh % SavesDone==0) CALL Info('SaveResult','Saving also coordinates',Level=12)
      END IF
      
      ! Write the header to file
      IF(Mesh % SavesDone==0) CALL Info('SaveResult','Writing the header part',Level=12)
      CALL InitializeFile( OutputUnit, FName, PosUnit, PosName )

      DateStr = FormatDate()
      WRITE( OutputUnit, '("!File started at: ",A)' ) TRIM(DateStr)
      DOFs = 0
      WRITE( OutputUnit,* ) 'Degrees of freedom: '

      ! Vectors
      DO IsVector = 1,0,-1

        i = 1
        Var => Mesh % Variables
        DO WHILE( ASSOCIATED(Var) )
          IF ( Var % Output ) THEN

            SaveThis = .FALSE.
            IF( SIZE(Var % Values) == Var % Dofs ) THEN
              SaveThis = SaveGlobal
            ELSE
              IF ( .NOT.SEQL(Var % Name,'coordinate') .OR. SaveCoordinates ) THEN
                SaveThis = .TRUE.
              END IF
            END IF
            
            IF( IsVector == 1 .AND. Var % Dofs == 1 ) SaveThis = .FALSE.
            IF( IsVector == 0 .AND. Var % Dofs > 1 ) SaveThis = .FALSE.

            k = LEN_TRIM(Var % Name)
            IF( OutputVariableList .AND. SaveThis ) THEN
              SaveThis = .FALSE.
              DO j=1,1000
                VarName = ListGetString(CurrentModel % Simulation,&
                    'Output Variable '//I2S(j), Found )
                IF( .NOT. Found ) EXIT
                k2 = LEN_TRIM(VarName)
                IF( VarName(1:k2) == Var % Name(1:k2) ) THEN
                  SaveThis = .TRUE.
                  ! This makes it possible to request saving of vectors
                  ! so that also all the corresponding scalar components (1,2,3,...) are saved. 
                  IF( k > k2 ) SaveThis = ( VERIFY( Var % Name(k2+1:k),' 0123456789') == 0 ) 
                  IF( SaveThis ) EXIT
                END IF
              END DO
            END IF

            IF( SaveThis ) THEN
              Found = .FALSE.
              IF( ASSOCIATED( Var % Solver ) ) THEN
                EqName = ListGetString(Var % Solver % Values,'Equation',Found)
              END IF
              IF(.NOT. Found ) EqName = 'no equation'
              
              IF( ASSOCIATED( Var % Perm ) ) THEN
                PermSize = SIZE( Var % Perm ) 
              ELSE
                PermSize = 0
              END IF
              
              IF( k < 25 ) THEN
                WRITE( OutputUnit,'(A,T25,A,I8,I8,I4,A)') Var % Name(1:k),' : ',&
                    SIZE(Var % Values), PermSize, Var % DOFs,' : '//TRIM( EqName )
              ELSE
                WRITE( OutputUnit,'(A,A,I8,I8,I4,A)') Var % Name(1:k),' : ',&
                    SIZE(Var % Values), PermSize, Var % DOFs,' : '//TRIM( EqName )
              END IF

              IF( IsVector == 0 ) THEN
                Dofs = Dofs + 1
                IF ( Binary ) THEN
                  CALL BinWriteString( PosUnit, Var % Name(1:k) )
                  i = i + 1
                END IF
              END IF
            END IF
                        
          END IF
          Var => Var % Next
        END DO
      END DO

      WRITE(OutputUnit,*) 'Total DOFs: ', DOFs
      WRITE(OutputUnit,*) 'Number Of Nodes: ', Mesh % NumberOfNodes


      IF ( Binary ) THEN
         ! Jump to the beginning (remembering the initial L|B!) in the positions
         ! file to fill in the number of variables, then jump back here again.
         Pos = BinFTell( PosUnit )
         CALL BinFSeek( PosUnit, 1_IntOff_k, BIN_SEEK_SET )
         CALL BinWriteInt4( PosUnit, i-1 )
         CALL BinFSeek( PosUnit, Pos, BIN_SEEK_SET )

         CALL SwitchToBinary( OutputUnit,FName,Mesh % NumberOfNodes )
      END IF
    ELSE
      CALL AppendOpen( OutputUnit,FName,PosUnit,PosName )
    END IF

    IF(Mesh % SavesDone==0) CALL Info('SaveResult','Writing data for the current timestep',Level=12)
    CALL WriteTime( OutputUnit,PosUnit,Mesh % SavesDone+1, Time, SimulationTime )

!------------------------------------------------------------------------------
!   Write data to disk
!------------------------------------------------------------------------------
    PrevPerm => NULL()
    Var => Mesh % Variables
    j = 1
    DO WHILE( ASSOCIATED(Var) )
      SaveThis = .FALSE.
      IF ( Var % Output .AND. Var % DOFs==1 ) THEN
        k = LEN_TRIM(Var % Name)
        IF( Var % DOFs == SIZE(Var % Values) ) THEN
          SaveThis = SaveGlobal
        ELSE
          IF ( .NOT.SEQL(Var % Name, 'coordinate') .OR. SaveCoordinates ) THEN
            SaveThis = .TRUE.
          END IF
        END IF
      END IF


      IF( OutputVariableList .AND. SaveThis ) THEN
        SaveThis = .FALSE.
        DO j=1,1000
          VarName = ListGetString(CurrentModel % Simulation,&
              'Output Variable '//I2S(j), Found )
          IF( .NOT. Found ) EXIT
          IF( TRIM( VarName ) == TRIM( Var % Name(1:k) ) ) THEN
            SaveThis = .TRUE.
            EXIT
          END IF
        END DO
      END IF
      
      IF( SaveThis ) THEN
        IF( SaveAll .OR. Var % ValuesChanged ) THEN
          CALL WriteVarName( OutputUnit,PosUnit,Var % Name(1:k),VarPos(j) )
          CALL WritePerm( OutputUnit, Var % Perm, PrevPerm )
          
          IF ( ASSOCIATED(Var % Perm) ) THEN
            n = SIZE(Var % Perm)
          ELSE
            n = SIZE( Var % Values )
          END IF
          
          DO i=1, n
            k = i
            IF ( ASSOCIATED(Var % Perm) ) k = Var % Perm(i)
            IF ( k > 0 ) THEN
              CALL WriteReal( OutputUnit,Var % Values(k) )
            END IF
          END DO
          Var % ValuesChanged = .FALSE.
        ELSE
          IF ( Binary ) CALL BinWriteInt8( PosUnit,INT(VarPos(j),Int8_k) )
        END IF
        j = j + 1
      END IF
      
      Var => Var % Next
    END DO
  
    IF ( Binary ) THEN
      CALL BinClose( OutputUnit )
      CALL BinClose( PosUnit )
    ELSE
      CLOSE( OutputUnit )
    END IF


    IF( Mesh % SavesDone == 0 ) THEN
      CALL Info('SaveResult','Done writing results file',Level=12)
      CALL Info('SaveResult','-----------------------------------------',Level=4)
    END IF

    Mesh % SavesDone = Mesh % SavesDone + 1
    SaveCount = Mesh % SavesDone 


  CONTAINS

      SUBROUTINE WritePerm( OutputUnit, CurrPerm, PrevPerm )
         INTEGER, INTENT(IN) :: OutputUnit
         INTEGER, POINTER :: CurrPerm(:), PrevPerm(:)
         INTEGER :: i, n
         LOGICAL :: SameAsPRev

         IF ( .NOT.ASSOCIATED( CurrPerm ) ) THEN

            IF ( Binary ) THEN
               CALL BinWriteInt4( OutputUnit,0 )
            ELSE
               WRITE( OutputUnit, '(A)' ) 'Perm: NULL'
            END IF

         ELSE
         
            SameAsPrev = .FALSE.
            IF ( ASSOCIATED( CurrPerm, PrevPerm ) ) THEN
               SameAsPrev = .TRUE.
            ELSE IF ( ASSOCIATED( PrevPerm ) ) THEN
               IF ( SIZE(CurrPerm) == SIZE(PrevPerm) ) THEN
                  IF ( ALL( CurrPerm == PrevPerm ) ) SameAsPrev = .TRUE.
               END IF
            END IF

            IF ( SameAsPrev ) THEN

               IF ( Binary ) THEN
                  CALL BinWriteInt4( OutputUnit,-1 )
                  CALL BinWriteInt8( OutputUnit,INT( PrevPermPos,Int8_k ) )
               ELSE
                  WRITE( OutputUnit, '(A)' ) 'Perm: use previous'
               END IF

            ELSE

               PrevPerm => Var % Perm

               n = COUNT( CurrPerm > 0 )
               IF ( Binary ) THEN
                  PrevPermPos = BinFTell( OutputUnit )
                  CALL BinWriteInt4( OutputUnit, SIZE(CurrPerm) )
                  CALL BinWriteInt4( OutputUnit, n )
                  DO i = 1, SIZE( CurrPerm )
                     IF ( CurrPerm(i) > 0 ) THEN
                        CALL BinWriteInt4( OutputUnit,i )
                        CALL BinWriteInt4( OutputUnit,CurrPerm(i) )
                     END IF
                  END DO
               ELSE
                  WRITE( OutputUnit,'(A,i12," ",i12)') 'Perm: ', SIZE(CurrPerm), n
                  DO i = 1, SIZE( CurrPerm )
                     IF ( CurrPerm(i) > 0 ) THEN
                        WRITE( OutputUnit, '(2i11)' ) i,CurrPerm(i)
                     END IF
                  END DO
               END IF

            END IF

         END IF
      END SUBROUTINE WritePerm


      SUBROUTINE WriteVarName( OutputUnit,PosUnit,Name,VarPos )
         INTEGER, INTENT(IN) :: OutputUnit,PosUnit
         CHARACTER(*), INTENT(IN) :: Name
         INTEGER(IntOff_k), INTENT(OUT) :: VarPos

         IF ( Binary ) THEN
            VarPos = BinFTell( OutputUnit )
            CALL BinWriteInt8( PosUnit, INT(VarPos,Int8_k) )

            CALL BinWriteString( OutputUnit,TRIM(Name) )
         ELSE
            WRITE( OutputUnit,'(a)' ) TRIM(Name)
         END IF
      END SUBROUTINE WriteVarName
        

      SUBROUTINE WriteReal( OutputUnit,r )
         INTEGER, INTENT(IN) :: OutputUnit
         REAL(dp), INTENT(IN) :: r

         IF ( Binary ) THEN
            CALL BinWriteDouble( OutputUnit,r )
         ELSE
            WRITE( OutputUnit,* ) r
         END IF
      END SUBROUTINE WriteReal


      SUBROUTINE WriteTime( OutputUnit,PosUnit,SavesDone,nTime,SimulationTime )
         INTEGER, INTENT(IN) :: OutputUnit, PosUnit
         INTEGER, INTENT(IN) :: SavesDone,nTime
         REAL(dp), INTENT(IN) :: SimulationTime
         INTEGER(IntOff_k) :: Pos

         IF ( Binary ) THEN
            Pos = BinFTell( OutputUnit )
            CALL BinWriteInt8( PosUnit, INT(Pos,Int8_k) )

            CALL BinWriteString( OutputUnit, 'Time:' )
            CALL BinWriteInt4( OutputUnit, SavesDone )
            CALL BinWriteInt4( OutputUnit, nTime )
            CALL BinWriteDouble( OutputUnit, SimulationTime )
         ELSE
            WRITE( OutputUnit,'(a,i7,i7,ES17.8E3)' ) 'Time: ',SavesDone,nTime, &
                                                   SimulationTime
         END IF
      END SUBROUTINE WriteTime


      SUBROUTINE AppendOpen( OutputUnit,FName,PosUnit,PosName )
         INTEGER, INTENT(IN) :: OutputUnit
         CHARACTER(*), INTENT(IN) :: FName
         INTEGER, INTENT(IN) :: PosUnit
         CHARACTER(*), INTENT(IN) :: PosName

         IF ( Binary ) THEN
            CALL BinOpen( OutputUnit,FName,"APPEND" )
            CALL BinOpen( PosUnit,PosName,"APPEND" )
         ELSE
            OPEN( OutputUnit,FILE=FName,STATUS="OLD",POSITION="APPEND" )
         END IF
      END SUBROUTINE AppendOpen
      

      SUBROUTINE SwitchToBinary( OutputUnit,FName,nNodes )
         INTEGER, INTENT(IN) :: OutputUnit
         CHARACTER(*), INTENT(IN) :: FName
         INTEGER, INTENT(IN) :: nNodes

         CLOSE( OutputUnit )
         CALL BinOpen( OutputUnit, TRIM(FName), "APPEND" )

         ! The binary part starts with a NULL byte.
         CALL BinWriteString( OutputUnit, "" )
      END SUBROUTINE SwitchToBinary
      

      SUBROUTINE InitializeFile( OutputUnit,FName,PosUnit,PosName )
         INTEGER, INTENT(IN) :: OutputUnit
         CHARACTER(*), INTENT(IN) :: FName
         INTEGER, INTENT(IN) :: PosUnit
         CHARACTER(*), INTENT(IN) :: PosName

         IF ( Binary ) THEN
            CALL BinEndianess( E )
            OPEN( OutputUnit,File=FName,STATUS='UNKNOWN' )
            WRITE( OutputUnit, * ) 'BINARY 3.', E

            CALL BinOpen( PosUnit,PosName,'WRITE' )
            CALL BinWriteChar( PosUnit,E )
            ! Make room for number of variables, we fill it in later
            CALL BinWriteInt4( PosUnit, 0 )
         ELSE
            OPEN( OutputUnit,File=FName,STATUS='UNKNOWN' )
            WRITE( OutputUnit, * ) 'ASCII 3'
         END IF
      END SUBROUTINE InitializeFile

!------------------------------------------------------------------------------
  END FUNCTION SaveResult
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Loads the result file that has been saved by an earlier Elmer simulation.
!> This makes it possible to restart the simulation.
!------------------------------------------------------------------------------
  SUBROUTINE LoadRestartFile( RestartFile,TimeCount,Mesh,Continuous,EOF,RestartList)
    CHARACTER(LEN=*) :: RestartFile
    INTEGER :: TimeCount
    TYPE(Mesh_T), POINTER :: Mesh
    LOGICAL, OPTIONAL :: Continuous,EOF
    TYPE(ValueList_t), POINTER, OPTIONAL :: RestartList
!------------------------------------------------------------------------------
    TYPE(Variable_t),POINTER :: Var, Comp
    CHARACTER(LEN=MAX_NAME_LEN) :: Name,VarName,VarName2,FullName,PosName
    CHARACTER(LEN=:), ALLOCATABLE :: Row
    CHARACTER(LEN=MAX_STRING_LEN) :: FName,Trash
    INTEGER ::i,j,k,k2,n,nt,Node,DOFs,TotalDofs,DofCount,SavedCount,Timestep,NSDOFs,nlen
    INTEGER :: nNodes, Stat, FieldSize, PermSize
    INTEGER, SAVE :: FmtVersion
    INTEGER, ALLOCATABLE :: Perm(:)

    TYPE(Solver_t),   POINTER :: Solver
    TYPE(Variable_t), POINTER :: TimeVar, tStepVar

    LOGICAL :: RestartFileOpen = .FALSE., Cont, Found, LoadThis
    LOGICAL, SAVE :: PosFile = .FALSE.
    LOGICAL, SAVE :: Binary, RestartVariableList, GotPerm, GotIt
    INTEGER, SAVE, ALLOCATABLE :: RestartVariableSizes(:)
    TYPE(ValueList_t), POINTER :: ResList
    
    REAL(KIND=dp) :: Dummy,Val,Time
    REAL(KIND=dp), POINTER :: Component(:), Temp(:)
    REAL(KIND=dp), POINTER :: Velocity1(:),Velocity2(:),Velocity3(:),Pressure(:)
    INTEGER(KIND=IntOff_k) :: Pos
    INTEGER :: iostat, FileCount
    CHARACTER(1) :: E
#ifdef USE_ISO_C_BINDINGS
    REAL(dp) :: tstart, tstop
#else
    REAL(dp) :: tstart, tstop, CPUTime
#endif
    SAVE TotalDOFs

    tstart = CPUTime()
!------------------------------------------------------------------------------
!   Open restart file and search for the right position
!------------------------------------------------------------------------------
    CALL Info( 'LoadRestartFile',' ', Level = 4)
    CALL Info( 'LoadRestartFile','--------------------------------------------', Level= 4 )
    CALL Info( 'LoadRestartFile','Reading data from file: '//TRIM(RestartFile), Level = 4 )

    IF( PRESENT( RestartList ) ) THEN
      ResList => RestartList
    ELSE
      ResList => CurrentModel % Simulation
    END IF
    
    RestartVariableList = ListCheckPresent( ResList,'Restart Variable 1')
    IF( RestartVariableList ) THEN
      CALL Info('LoadRestartFile','Reading only variables given by: Restart Variable i',Level=10)
    ELSE
      CALL Info('LoadRestartFile','Reading all variables (if not wanted use >Restart Variable i< )',Level=10)      
    END IF
    
    Cont = .FALSE.
    IF ( PRESENT( Continuous ) ) Cont = Continuous
    IF ( PRESENT( EOF ) ) EOF = .FALSE.

    IF ( Cont .AND. RestartFileOpen ) GOTO 30

    IF ( .NOT. FileNameQualified(RestartFile) .AND. LEN_TRIM(OutputPath)>0 ) THEN
      FName = TRIM(OutputPath) // '/' // TRIM(RestartFile)
    ELSE
      FName = RestartFile
    END IF
    OPEN( RestartUnit,File=TRIM(FName),STATUS='OLD',IOSTAT=iostat )

    IF( iostat == 0 ) THEN
      FileCount = 1
    ELSE
      FileCount = 0
    END IF
 
    FileCount = NINT( ParallelReduction( 1.0_dp * FileCount ) )
    IF( FileCount == 0 ) THEN
      CALL Error( 'LoadRestartFile','=======================================' )
      CALL Error( 'LoadRestartFile','' )
      CALL Error( 'LoadRestartFile','Could not open file "'//TRIM(FName)//'"' )
      CALL Error( 'LoadRestartFile','No restart possible!' )
      CALL Error( 'LoadRestartFile','' )
      CALL Fatal( 'LoadRestartFile','=======================================' )
    ELSE IF( FileCount < ParEnv % PEs ) THEN
      CALL Info('LoadRestartFile','Succefully opened '//TRIM(I2S(FileCount))//&
          ' restart files out of '//TRIM(I2S(ParEnv % PEs)),Level=6)
      IF( ListGetLogical( ResList,'Restart Error Continue',Found ) ) THEN
        ! This partition does not have a mesh
        IF( iostat /= 0 ) RETURN 
      ELSE
        CALL Error( 'LoadRestartFile','=======================================' )
        CALL Error( 'LoadRestartFile','' )
        CALL Error( 'LoadRestartFile','Expted to find all restart files "'//TRIM(FName)//'"' )
        CALL Error( 'LoadRestartFile','No restart possible!' )
        CALL Error( 'LoadRestartFile','' )
        CALL Fatal( 'LoadRestartFile','=======================================' )
      END IF
    END IF
    
    RestartFileOpen = .TRUE.

    ALLOCATE(CHARACTER(MAX_STRING_LEN)::Row)
    READ( RestartUnit, '(A)', IOSTAT=iostat ) Row
    IF( iostat /= 0 ) THEN
      CALL Fatal('LoadRestartFile','Error reading header line!')
    END IF

    IF ( Row(3:8) == 'BINARY' ) THEN
        Binary = .TRUE.
        FmtVersion = 1
        E = Row(1:1)
        CALL Info( 'LoadRestartFile', TRIM(Row(3:)), Level = 4 )
    ELSE IF ( Row(2:7) == 'BINARY' ) THEN
        Binary = .TRUE.
        FmtVersion = 2
        READ( Row(9:9),'(I1)',IOSTAT=iostat) FmtVersion
        IF( iostat /= 0 ) THEN
          CALL Fatal('LoadRestartFile','Error reading version: '//TRIM(Row))
        END IF
        i = INDEX( Row, "." )
        E = Row(i+1:i+1)
        CALL Info( 'LoadRestartFile', TRIM(Row(2:)), Level = 4 )
    ELSE IF ( Row(2:6) == 'ASCII' ) THEN
        Binary = .FALSE.
        READ( Row(8:8),'(I1)',IOSTAT=iostat) FmtVersion
        IF( iostat /= 0 ) THEN
          CALL Fatal('LoadRestartFile','Error reading version: '//TRIM(Row))
        END IF
        CALL Info( 'LoadRestartFile', TRIM(Row(2:)), Level = 4 )
    ELSE 
        Binary = .FALSE.
        FmtVersion = 0
        CALL Info( 'LoadRestartFile', 'ASCII 0', Level = 4 )
    END IF
    
    CALL Info( 'LoadRestartFile','Reading restart file version '//TRIM(I2S(FmtVersion)), Level = 4)

    ! If we want to skip some of the variables we need to have a list 
    ! of their sizes still. This is particularly true with variables that 
    ! do not have permutation since they could be a field (like coordinate)
    ! or a global variable (like time).
    !----------------------------------------------------------------------
    IF( RestartVariableList ) THEN      
      DO WHILE( ReadAndTrim(RestartUnit,Row) )
        nlen = LEN_TRIM(Row)        
        k = INDEX( Row(1:nlen),'total dofs:',.TRUE.) 
        IF( k /= 0 ) THEN
          READ( Row(k+11:nlen),*,IOSTAT=iostat ) TotalDofs
          IF( iostat /= 0 ) THEN
            CALL Fatal('LoadRestartFile','Unable to load total dofs!')
          END IF
          EXIT
        END IF
      END DO
      ALLOCATE( RestartVariableSizes(TotalDofs) )
      RestartVariableSizes = 0
      REWIND( RestartUnit )
    END IF

    ! Find the start of dof definition part
    ! Here we use the INDEX so that there could be some empty space
    ! also before the keyword. 
    !----------------------------------------------------------------
    DO WHILE( ReadAndTrim(RestartUnit,Row) )
      IF( INDEX( Row(1:20),'degrees of freedom' ) /= 0 ) EXIT
    END DO

!------------------------------------------------------------------------------
    DofCount = 0

    DO WHILE( ReadAndTrim(RestartUnit,Row) )

      nlen = LEN_TRIM(Row)

      k = INDEX( Row(1:nlen),'total dofs:',.TRUE.) 
      IF( k /= 0 ) THEN
        READ( Row(k+11:nlen),*,IOSTAT=iostat ) TotalDofs
        IF( iostat /= 0 ) THEN
          CALL Fatal('LoadRestartFile','Unable to load total dofs!')
        END IF
        CALL Info('LoadRestartFile','Total number of dofs to load: '//I2S(TotalDofs) )
        EXIT
      END IF

      IF( FmtVersion < 3 ) THEN
        ! Figure out what is the solver to which the variable is associated to
        ! this requires that the 'Equation' keyword is unique.
        ! I wonder if this is used at all?
        k = INDEX(Row,']')+1
        
        ! The last colon in the line
        k = k+INDEX(Row(k:),':',.TRUE.)-1
        NULLIFY(Solver)
        DO i = 1,CurrentModel % NumberOfSolvers
          Solver => CurrentModel % Solvers(i)
          IF ( Row(k+1:nlen) == ListGetString(Solver % Values, 'Equation',Found)) EXIT
        END DO
        
        ! Figure out the slot where the number of dofs are given and read them
        ! The rule is to start from ':' and go throug empty space and occupied space
        DO j=k-1,1,-1
          IF ( Row(j:j) /= ' ' ) EXIT
        END DO
        DO k=j,1,-1
          IF ( Row(k:k) == ' ' ) EXIT
        END DO
        READ(Row(k+1:),*,IOSTAT=iostat) DOFs
        IF( iostat /= 0 ) THEN
          CALL Fatal('LoadRestartFile','Error reading DOFs: '//TRIM(Row(k+1:)))
        END IF

        IF( Dofs < 1 ) CALL Fatal('LoadRestartFile','The Dofs should be positive: '//i2s(DOFs))
        
        ! Figure out the name of the variable
        j = INDEX(Row,'[')
        IF( j > 0 ) THEN
          VarName = TRIM(Row(1:j-1))
        ELSE
          VarName = TRIM(Row(1:k-1))
        END IF
        FieldSize = 0
        FullName = VarName
      ELSE IF( FmtVersion == 3 ) THEN

        ! read the field names
        ! the full name includes the possible info behind bracets [...]
        j = INDEX( Row,']')
        IF( j == 0 ) THEN
          ! names are the same
          j = INDEX( Row,':') 
          IF( j > 1 ) THEN
            VarName = TRIM(Row(1:j-1))
            FullName = VarName
          ELSE
            CALL Warn('LoadRestartFile','Cannot read variable name: '//TRIM(Row))
          END IF
        ELSE
          FullName = TRIM(Row(1:j))
          j = INDEX( Row,'[')
          IF( j == 0 ) THEN
            CALL Warn('LoadRestartFile','Missing left parenthesis: '//TRIM(Row))
            VarName = FullName
          ELSE
            VarName = TRIM(Row(1:j-1))
          END IF
        END IF

        CALL Info('LoadRestartFile','Reading variable: '//TRIM(VarName),Level=12)
        
        ! read the size of field, size or perm and number of dofs per node
        !-----------------------------------------------------------------
        j = MAX(INDEX(Row,']'),1)
        k = INDEX( Row(j:),':')
        j = j+k
        READ(Row(j+1:),*,IOSTAT=iostat) FieldSize,PermSize,DOFs
        IF( iostat /= 0 ) THEN
          CALL Fatal('LoadRestartFile','Error reading size information: '//TRIM(Row(j+1:)))
        END IF

        ! read the name of the solver and associate it to existing solver
        !----------------------------------------------------------------
        k = INDEX( Row(j+1:),':')
        j = j+k
        DO k=j+1,nlen
          IF( Row(k:k) /= ' ') EXIT
        END DO
        NULLIFY(Solver)
        Found = .FALSE.
        DO i = 1,CurrentModel % NumberOfSolvers
          Solver => CurrentModel % Solvers(i)
          IF ( Row(k:nlen) == TRIM( ListGetString(Solver % Values, 'Equation',GotIt) ) ) THEN
            Found = .TRUE. 
            EXIT            
          END IF
        END DO

        IF( Found ) THEN
          CALL Info('LoadRestartFile','Associated to existing solver',Level=20)
        ELSE IF( PermSize > 0 ) THEN
          CALL Warn('LoadRestartFile','Could not associate variable to solver: '//TRIM(VarName))
          ! Associated to some solver at least 
          Solver => CurrentModel % Solvers(1)
        END IF
      ELSE
        CALL Fatal('LoadRestartFile','Unknown file version: '//I2S(FmtVersion))
      END IF

      ! The user may optionally give a list of the fields to be loaded from the 
      ! result file. All fields that will be loaded will also be allocated for and
      ! hence this feature could save some resources.
      !---------------------------------------------------------------------------
      LoadThis = .TRUE.
      IF( RestartVariableList ) THEN
        IF( Dofs == 1 ) DofCount = DofCount + 1
        LoadThis = .FALSE.
        k = LEN_TRIM( VarName )
        DO j=1,1000
          VarName2 = ListGetString( ResList,'Restart Variable '//I2S(j), Found )
          IF( .NOT. Found ) EXIT
          k2 = LEN_TRIM(VarName2)

          IF( VarName2(1:k2) == VarName(1:k2) ) THEN
            LoadThis = .TRUE.
            ! This makes it possible to request loading of vectors
            ! so that also all the corresponding scalar components (1,2,3,...) are saved. 
            IF( k > k2 ) LoadThis = ( VERIFY( VarName(k2+1:k),' 0123456789') == 0 ) 
            IF( LoadThis ) EXIT
          END IF
        END DO        
        IF(.NOT. LoadThis ) THEN
          IF( FieldSize == 0 ) THEN
            CALL Fatal('LoadRestartFile','Cannot pick variables with old format!')
          END IF
          IF( Dofs == 1 ) THEN
             IF(PermSize > 0) THEN
                RestartVariableSizes(DofCount) = PermSize
             ELSE
                RestartVariableSizes(DofCount) = FieldSize
             END IF
          END IF
        END IF
      END IF

      IF( LEN_TRIM( VarName ) == 0 ) THEN
        CALL Warn('LoadRestartFile','Could not deduce variable name!')
        LoadThis = .FALSE.
      END IF


      IF(.NOT. LoadThis) CYCLE

      ! Check whether a variable exists or not. If it does not exist then 
      ! create the variable so that it can be filled with the data.
      !------------------------------------------------------------------
      Var => VariableGet( Mesh % Variables, VarName,.TRUE. )      
      IF ( .NOT.ASSOCIATED(Var) ) THEN
        CALL Info('LoadRestartFile','Creating variable: '//TRIM(VarName),Level=6)

        ALLOCATE( Var )

        ! The old format (ver. < 3) does not include information on vector sizes prior to loading
        ! thus make an educated guess.
        !----------------------------------------------------------------------------------------
        IF( FieldSize == 0 .AND. FmtVersion<3 ) THEN
          FieldSize = Mesh % NumberOfNodes * DOFs
          PermSize = Mesh % NumberOfNodes
        END IF

        ALLOCATE( Var % Values(FieldSize) )
        Var % Values = 0.0          

        IF( PermSize > 0 ) THEN
          ALLOCATE( Var % Perm(PermSize) )
          Var % Perm = 0
        END IF

        IF ( SEQL(VarName, 'flow solution ') ) THEN
!------------------------------------------------------------------------------
!         First add components to the variable list separately...
!         (must be done this way for the output routines to work properly...)
!----------------------------------------------------------------------------
!          NSDOFs = CoordinateSystemDimension() + 1
!          IF ( Coordinates == CylindricSymmetric ) NSDOFs = NSDOFs + 1
          NSDOFS = Dofs

          Velocity1 => Var % Values(1:NSDOFs*Mesh % NumberOfNodes:NSDOFs)
          CALL VariableAdd( Mesh % Variables,  Mesh, Solver, 'Velocity 1', &
                       1, Velocity1, Var % Perm )

          Velocity2 => Var % Values(2:NSDOFs*mesh % NumberOfNodes:NSDOFs)
          CALL VariableAdd( Mesh % Variables, Mesh, Solver, 'Velocity 2', &
                       1, Velocity2, Var % Perm )

          IF ( NSDOFs == 3 ) THEN
            Pressure => Var % Values(3:NSDOFs*Mesh % NumberOfNodes:NSDOFs)
            CALL VariableAdd( Mesh % Variables, Mesh, Solver, 'Pressure', &
                       1, Pressure, Var % Perm )
          ELSE
            Velocity3 => Var % Values(3:NSDOFs*Mesh % NumberOfNodes:NSDOFs)
            CALL VariableAdd( Mesh % Variables, Mesh, Solver, 'Velocity 3', &
                       1, Velocity3, Var % Perm )

            Pressure => Var % Values(4:NSDOFs*Mesh % NumberOfNodes:NSDOFs)
             CALL VariableAdd( Mesh % Variables, Mesh, Solver, 'Pressure', &
                         1, Pressure, Var % Perm )
           END IF
!------------------------------------------------------------------------------
!        Then add the thing itself
!------------------------------------------------------------------------------
         CALL VariableAdd( Mesh % Variables, Mesh, Solver, &
                'Flow Solution',NSDOFs,Var % Values,Var % Perm )
        ELSE IF( PermSize == 0 ) THEN
          IF ( DOFs > 1 ) THEN
            DO i=1,DOFs
              Component => Var % Values(i:FieldSize:DOFs)
              name = ComponentName( FullName, i )
              CALL VariableAdd( Mesh % Variables,  Mesh, Solver, name, &
                  1, Component )
            END DO
          END IF
          CALL VariableAdd( Mesh % Variables, Mesh, Solver, &
              FullName,DOFs,Var % Values) 
        ELSE
          IF ( DOFs > 1 ) THEN
            DO i=1,DOFs
              Component => Var % Values(i:FieldSize:DOFs)
              name = ComponentName( FullName, i )
              CALL VariableAdd( Mesh % Variables,  Mesh, Solver, name, &
                           1, Component, Var % Perm )
            END DO
          END IF
          CALL VariableAdd( Mesh % Variables, Mesh, Solver, &
              FullName,DOFs,Var % Values,Var % Perm )
        END IF
      END IF
    END DO

    IF ( Binary ) THEN
      ! Switch to binary reading
      CLOSE( RestartUnit )
      CALL BinOpen( RestartUnit,FName,"read" )
      CALL BinSetInputEndianess( RestartUnit,E )
      CALL BinReadString( RestartUnit,Trash ) !Skip header; we read it already

      PosName = TRIM(FName) // '.pos'
      INQUIRE( FILE=PosName, EXIST=PosFile )
      IF ( PosFile ) THEN
         CALL BinOpen( PosUnit,PosName,'read' )
         CALL BinReadString( PosUnit,E ) 
         CALL BinSetInputEndianess( PosUnit,E )
      END IF
    END IF

30  CONTINUE

!------------------------------------------------------------------------------
!   ...read one timestep to memory...
!------------------------------------------------------------------------------
    IF ( Cont ) THEN
       nt = TimeCount
    ELSE
       nt = 1
    END IF

    IF ( PosFile ) THEN
      ! If TimeCount == 0, we'll get the last time step
      Pos = GetPosition( PosUnit,TimeCount,0,TimeCount )
      nt = TimeCount
      CALL BinFSeek( RestartUnit,Pos,BIN_SEEK_SET )
    END IF

    DO WHILE( nt <= TimeCount .OR. TimeCount == 0 )
      CALL ReadTime( RestartUnit,SavedCount,Timestep,Time,Stat )

      IF ( Stat /= 0 ) THEN
         IF ( PRESENT( EOF ) ) THEN
            EOF = .TRUE.
         ELSE IF ( TimeCount /= 0 ) THEN
            CALL Warn( 'LoadRestartFile','Did not find data at the requested point' )
            WRITE( Message,'(A,I6)') 'Reading from the last existing point: ',nt-1
            CALL Warn( 'LoadRestartFile',Message)
         END IF
         EXIT
      END IF

      TimeVar  => VariableGet( Mesh % Variables, 'Time' )
      tStepVar => VariableGet( Mesh % Variables, 'Timestep' )

      IF ( ASSOCIATED( TimeVar ) )  TimeVar % Values(1)  = Time
      IF ( ASSOCIATED( tStepVar ) ) tStepVar % Values(1) = Timestep

      WRITE( Message,'(A,ES12.3)') 'Reading time sequence: ',Time
      CALL Info( 'LoadRestartFile',Message, Level=4)

      WRITE( Message,'(A,I0)') 'Reading timestep: ',Timestep
      CALL Info( 'LoadRestartFile',Message, Level=4)

      DO i=1,TotalDOFs

        LoadThis = .TRUE.
        IF( RestartVariableList ) THEN
          n = RestartVariableSizes(i)
          LoadThis = ( n == 0 ) 
        END IF

         IF ( PosFile ) THEN
            Pos = GetPosition( PosUnit,TimeCount,i )
            CALL BinFSeek( RestartUnit,Pos,BIN_SEEK_SET )
         END IF

         CALL ReadVariableName( RestartUnit,Row,Stat )

         ! If not all varibales were saved for this time step, and we're not
         ! using a .pos file, we may have reached the end even though i < TotalDOFs.
         IF ( Stat /= 0 ) EXIT
         IF ( SEQL(Row, "Time:") ) THEN
            CALL UnReadLine( RestartUnit, Row )
            EXIT
         END IF

         IF( LoadThis ) THEN
           CALL Info('LoadRestartFile','Reading Variable: '//TRIM(Row),Level=12)
         ELSE
           CALL Info('LoadRestartFile','Skipping Variable: '//TRIM(Row),Level=12)          
         END IF

         ! Note that Var % Perm is the permutation associated with the current field
         ! while Perm will be the permutation associated with the saved field. 
         ! They could be different, even though the usually are not!
         IF ( FmtVersion > 0 ) THEN
           CALL ReadPerm( RestartUnit, Perm, GotPerm )
         END IF

         IF( LoadThis ) THEN
           Var => VariableGet( Mesh % Variables,Row, ThisOnly=.TRUE. )
           IF ( .NOT. ASSOCIATED(Var) ) THEN
             CALL Warn('LoadRestartFile','Variable is not present for reading: '//TRIM(Row))
           END IF
           IF( GotPerm .NEQV. ASSOCIATED( Var % Perm ) ) THEN
             DEALLOCATE( Var % Perm ) ; Var % Perm => NULL()
!            CALL Fatal('LoadRestartFile','Permutation should either exist or not!')
           END IF

           IF ( ASSOCIATED(Var % Perm) .AND. FmtVersion > 0 ) THEN
             n = SIZE(Var % Perm)
           ELSE
             n = SIZE(Var % Values)
           END IF
           ! in case of (.NOT. LoadThis) n has already been set
           CALL Info('LoadRestartFile','Size of variable is '//TRIM(I2S(n)),Level=20)
           
           ! This relies that the "Transient Restart" flag has been used consistently when saving and loading
           IF( ASSOCIATED( Var % Solver ) ) THEN
             IF( ListGetLogical( Var % Solver % Values,'Transient Restart',Found ) ) THEN
               CALL Info('LoadRestartFile','Assuming variable to have transient initialization: '//TRIM(Row),Level=6)
               Var % Solver % DoneTime = Var % Solver % Order
             END IF
           END IF
         END IF
         
         DO j=1, n
           IF ( FmtVersion > 0 ) THEN
             CALL GetValue( RestartUnit, Perm, GotPerm, j, k, Val )
           ELSE
             READ( RestartUnit,* ) Node, k, Val
           END IF
           
           ! One can not really omit reading the lines since otherwise at least the 
           ! ascii format would loose it, but now we can cycle the rest. 
           IF(.NOT. LoadThis ) CYCLE
           
           IF ( .NOT. GotPerm ) THEN
             Var % Values(k) = Val
           ELSE IF ( Var % Perm(j) > 0 ) THEN
             IF ( k > 0 ) Var % Values(Var % Perm(j)) = Val
           ELSE 
             Var % Perm(j) = k
             IF ( k > 0 ) Var % Values(k) = Val
           END IF
         END DO

         IF( InfoActive( 20 ) ) THEN
           PRINT *,'LoadRestartFile range:',ParEnv % MyPe, MINVAL( Var % Values ), MAXVAL( Var % Values )
         END IF
           
         IF( LoadThis ) THEN
           CALL InvalidateVariable( CurrentModel % Meshes, Mesh, Row )
         END IF

      END DO
      nt = nt + 1
    END DO
!------------------------------------------------------------------------------

    IF ( .NOT. Cont ) THEN
      IF ( Binary ) THEN
         CALL BinClose( RestartUnit )
         IF ( PosFile ) CALL BinClose( PosUnit )
       ELSE
         CLOSE( RestartUnit )
       END IF
       RestartFileOpen = .FALSE.
    END IF


    ! This is now obsolite for the new format 
    IF( FmtVersion < 3 ) THEN
      ! Change variable allocations to correct sizes,
      ! first for vectors...
      ! ---------------------------------------------
      Var => Mesh % Variables
      DO WHILE(ASSOCIATED(Var))
        IF ( ASSOCIATED(Var % Perm) .AND. Var % DOFs>1 ) THEN
          n = Var % DOFs*COUNT(Var % Perm>0)
          IF ( SIZE(Var % Values) /= n) THEN
            ALLOCATE(Temp(n))
            Temp = Var % Values(1:n)
            DEALLOCATE(Var % Values)
            Var % Values => Temp
            DO i=1,Var % DOFs
              Comp => VariableGet(Mesh % Variables, &
                  ComponentName(Var,i),ThisOnly=.TRUE.)
              Comp % Values => Var % Values(i::Var % DOFs)
            END DO
          END IF
        END IF
        Var => Var % Next
      END DO

      !... and then for scalars
      ! -----------------------
      Var => Mesh % Variables
      DO WHILE(ASSOCIATED(Var))
        IF ( ASSOCIATED(Var % Perm) .AND. Var % DOFs==1 ) THEN
          n = COUNT(Var % Perm>0)
          IF ( SIZE(Var % Values) /= n) THEN
            ALLOCATE(Temp(n))
            Temp = Var % Values(1:n)
            DEALLOCATE(Var % Values)
            Var % Values => Temp
          END IF
        END IF
        Var => Var % Next
      END DO
    END IF

    
    tstop = CPUTime()
    
    WRITE( Message,'(A,ES15.4)') 'Time spent for restart (s): ', tstop - tstart
    CALL Info( 'LoadRestartFile',Message, Level = 4)
    CALL Info( 'LoadRestartFile', 'All done', Level = 4 )
    CALL Info( 'LoadRestartFile','--------------------------------------------', Level = 4 )


CONTAINS

   SUBROUTINE UnReadLine( Unit,Line )
      INTEGER, INTENT(IN) :: Unit
      CHARACTER(*), INTENT(IN) :: Line
      INTEGER(IntOff_k) :: Offset

      IF ( Binary ) THEN
         Offset = LEN_TRIM(Line) + 1
         CALL BinFSeek( Unit, -Offset, BIN_SEEK_CUR )
      ELSE
         BACKSPACE Unit
      END IF
   END SUBROUTINE UnReadLine


   INTEGER(IntOff_k) FUNCTION GetPosition( PosUnit,TimeStep,VarNr,FoundTStep ) &
                     RESULT(Pos)
      INTEGER, INTENT(IN) :: PosUnit
      INTEGER, INTENT(IN) :: TimeStep
      INTEGER, INTENT(IN) :: VarNr  ! 0 = time, 1 = first var, etc.
      INTEGER, INTENT(OUT),OPTIONAL :: FoundTStep
      CHARACTER(40) :: VarName
      INTEGER :: iTime, nVar, i, Stat
      INTEGER(IntOff_k) :: Offset, Offset2
      INTEGER(Int8_k) :: tmp
      INTEGER(IntOff_k), SAVE :: HeaderEnd, TimeStepSize = 0

      IF ( TimeStepSize == 0 ) THEN
         CALL BinReadInt4( PosUnit, nVar )
         TimeStepSize = (nVar + 1)*8 ! 8 bytes per variable + 8 bytes for time

         DO i = 1, nVar
            CALL BinReadString( PosUnit, VarName )
         END DO

         HeaderEnd = BinFTell( PosUnit )
      END IF

      IF ( TimeStep > 0 ) THEN
         Offset = (TimeStep - 1)*TimeStepSize + VarNr*8
         CALL BinFSeek( PosUnit, HeaderEnd + Offset, BIN_SEEK_SET )
         CALL BinReadInt8( PosUnit, tmp, Stat )

         IF ( Stat == 0 ) THEN
            Pos = tmp
            IF ( PRESENT(FoundTStep) ) FoundTStep = TimeStep
         ELSE
            CALL Warn( 'LoadRestartFile',&
                 'Did not find the the requested timestep in the positions file;' )
            CALL Warn( 'LoadRestartFile','using the last found one instead.')
            Offset2 = -TimeStepSize + VarNr*8
            CALL BinFSeek( PosUnit, Offset2 , BIN_SEEK_END )
            CALL BinReadInt8( PosUnit, tmp )
            Pos = tmp

            IF ( PRESENT(FoundTStep) ) THEN
               CALL BinFSeek( PosUnit, 0_IntOff_k, BIN_SEEK_END )
               Offset = BinFTell( PosUnit )
               FoundTStep = (Offset - HeaderEnd - VarNr*8)/TimeStepSize
            END IF
         END IF
      ELSE
         ! Find last time step
         Offset2 = -TimeStepSize + VarNr*8
         CALL BinFSeek( PosUnit, Offset2 , BIN_SEEK_END )
         CALL BinReadInt8( PosUnit, tmp )
         Pos = tmp
         IF ( PRESENT(FoundTStep) ) THEN
            CALL BinFSeek( PosUnit, 0_IntOff_k, BIN_SEEK_END )
            Offset = BinFTell( PosUnit )
            FoundTStep = (Offset - HeaderEnd - VarNr*8)/TimeStepSize
         END IF
      END IF
   END FUNCTION GetPosition
      

   SUBROUTINE GetValue( RestartUnit, Perm, UsePerm, iNode, iPerm, Val )
   !
   ! Get iPerm = Perm(iNode) and Value from input
   !
      INTEGER, INTENT(IN) :: RestartUnit
      INTEGER, ALLOCATABLE :: Perm(:)
      LOGICAL :: UsePerm
      INTEGER, INTENT(IN) :: iNode
      INTEGER, INTENT(OUT) :: iPerm
      REAL(dp), INTENT(OUT) :: Val

      IF ( UsePerm ) THEN
         iPerm = Perm(iNode)
      ELSE
         iPerm = iNode
      END IF

      IF ( iPerm > 0 ) THEN
         IF ( Binary ) THEN
            CALL BinReadDouble( RestartUnit, Val )
         ELSE
            READ( RestartUnit, * , IOSTAT=iostat ) Val
            IF( iostat /= 0 ) THEN
              CALL Fatal('LoadRestartFile','Error in GetValue')
            END IF
         END IF
      END IF
   END SUBROUTINE GetValue
    

   SUBROUTINE ReadPerm( RestartUnit, Perm, GotPerm )
      INTEGER, INTENT(IN) :: RestartUnit
      INTEGER, ALLOCATABLE :: Perm(:)
      LOGICAL :: GotPerm
      INTEGER :: nPerm, nPositive, i, j, k
      CHARACTER(MAX_NAME_LEN) :: Row
      INTEGER(Int8_k) :: Pos

      GotPerm = .FALSE.
      IF ( Binary ) THEN
         CALL BinReadInt4( RestartUnit, nPerm )
      ELSE
         READ( RestartUnit, '(A)' ) Row
         IF ( Row(7:10) == "NULL" ) THEN
            nPerm = 0
         ELSE IF ( Row(7:18) == "use previous" ) THEN
            nPerm = -1
         ELSE
            READ( Row(7:),*,IOSTAT=iostat) nPerm, nPositive
            IF( iostat /= 0 ) THEN
              CALL Fatal('LoadRestartFile','Error in ReadPerm: '//TRIM(Row))
            END IF
         END IF
      END IF

      IF ( nPerm < 0 ) THEN
         IF ( Binary ) CALL BinReadInt8( RestartUnit, Pos )
         ! At the moment, we always read all variables, and can therefore
         ! safely assume that the "previous" Perm table has been read and is
         ! held in memory at this point. In the future, however, we might be
         ! asked to read only some variables, in which case the previous Perm
         ! table might be yet unread so we need to jump back to 'Pos' and read
         ! it from there.
         GotPerm = .TRUE.
         RETURN
      ELSE IF ( nPerm == 0 ) THEN
!         IF ( ASSOCIATED(Perm) ) DEALLOCATE( Perm )
         RETURN
      ELSE 
         IF ( Binary ) CALL BinReadInt4( RestartUnit, nPositive )
      END IF

      IF( ALLOCATED( Perm ) ) THEN
        IF( SIZE( Perm ) < nPerm ) DEALLOCATE( Perm ) 
      END IF
      IF( .NOT. ALLOCATED( Perm ) ) THEN
        ALLOCATE( Perm(nPerm) )
      END IF
      Perm = 0

      DO i = 1, nPositive
         IF ( Binary ) THEN
            CALL BinReadInt4( RestartUnit, j )
            CALL BinReadInt4( RestartUnit, k )
         ELSE
            READ( RestartUnit, * , IOSTAT=iostat) j,k
            IF( iostat /= 0 ) THEN
              CALL Fatal('LoadRestartFile','Error reading in ReadPerm')
            END IF
         END IF
         Perm(j) = k
      END DO

      GotPerm = .TRUE.

   END SUBROUTINE ReadPerm


   SUBROUTINE ReadVariableName( RestartUnit,VarName,Stat )
      INTEGER, INTENT(IN) :: RestartUnit
      CHARACTER(*), INTENT(OUT) :: VarName
      INTEGER, INTENT(OUT) :: Stat

      IF ( Binary ) THEN
         CALL BinReadString( RestartUnit,VarName,Stat )
      ELSE
         READ( RestartUnit,FMT='(a)',IOSTAT=Stat ) VarName
      END IF
   END SUBROUTINE ReadVariableName


  SUBROUTINE ReadTime( RestartUnit,SavedCount,Timestep,Time,Stat )
  ! Read SavedCount, Timestep and Time.  Stat is set to 0 for success, < 0 for
  ! end-of-file, and > 0 for error.
    INTEGER, INTENT(IN) :: RestartUnit
    INTEGER, INTENT(OUT) :: SavedCount, Timestep
    REAL(dp), INTENT(OUT) :: Time
    INTEGER, INTENT(OUT) :: Stat
    CHARACTER(LEN=:), ALLOCATABLE :: String
    ALLOCATE(CHARACTER(32)::String)

    IF ( Binary ) THEN
      CALL BinReadString( RestartUnit, String, Stat )
      IF (Stat /= 0) RETURN
      CALL BinReadInt4( RestartUnit, SavedCount, Stat )
      CALL BinReadInt4( RestartUnit, Timestep, Stat )
      CALL BinReadDouble( RestartUnit, Time, Stat )
    ELSE
      DO WHILE( ReadAndTrim(RestartUnit,String) )
        IF ( SEQL(String, 'time:') ) THEN
          READ( String(7:),*,IOSTAT=iostat) SavedCount,Timestep,Time
          IF( iostat /= 0 ) THEN
            CALL Fatal('LoadRestartFile','Error in ReadTime!')
          END IF
          Stat = 0
          RETURN
        END IF
      END DO
      Stat = -1
    END IF
  END SUBROUTINE ReadTime

!------------------------------------------------------------------------------
  RECURSIVE SUBROUTINE InvalidateVariable( TopMesh,PrimaryMesh,Name )
!------------------------------------------------------------------------------
    CHARACTER(LEN=*) :: Name
    TYPE(Mesh_t),  POINTER :: TopMesh,PrimaryMesh
!------------------------------------------------------------------------------
    CHARACTER(LEN=MAX_NAME_LEN) :: tmpname
    INTEGER :: i
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Variable_t), POINTER :: Var,Var1
!------------------------------------------------------------------------------
    Mesh => TopMesh

    DO WHILE( ASSOCIATED(Mesh) )
      IF ( .NOT.ASSOCIATED( PrimaryMesh, Mesh) ) THEN
        Var => VariableGet( Mesh % Variables, Name, .TRUE.)
        IF ( ASSOCIATED( Var ) ) THEN
          Var % Valid = .FALSE.
          Var % PrimaryMesh => PrimaryMesh
          IF ( Var % DOFs > 1 ) THEN

! This should not be needed no more
            IF ( .FALSE. .AND. Var % Name == 'flow solution' ) THEN
              Var1 => VariableGet( Mesh % Variables, 'Velocity 1', .TRUE.)
              IF ( ASSOCIATED( Var1 ) ) THEN
                 Var1 % Valid = .FALSE.
                 Var1 % PrimaryMesh => PrimaryMesh
              END IF
              Var1 => VariableGet( Mesh % Variables, 'Velocity 2', .TRUE.)
              IF ( ASSOCIATED( Var1 ) ) THEN
                 Var1 % Valid = .FALSE.
                 Var1 % PrimaryMesh => PrimaryMesh
              END IF
              Var1 => VariableGet( Mesh % Variables, 'Velocity 3', .TRUE.)
              IF ( ASSOCIATED( Var1 ) ) THEN
                 Var1 % Valid = .FALSE.
                 Var1 % PrimaryMesh => PrimaryMesh
              END IF
              Var1 => VariableGet( Mesh % Variables, 'Pressure', .TRUE.)
              IF ( ASSOCIATED( Var1 ) ) THEN
                 Var1 % Valid = .FALSE.
                 Var1 % PrimaryMesh => PrimaryMesh
              END IF
              Var1 => VariableGet( Mesh % Variables, 'Surface', .TRUE.)
              IF ( ASSOCIATED( Var1 ) ) THEN
                 Var1 % Valid = .FALSE.
                 Var1 % PrimaryMesh => PrimaryMesh
              END IF
            ELSE
              DO i=1,Var % DOFs
                tmpname = ComponentName( Name, i )
                Var1 => VariableGet( Mesh % Variables, tmpname, .TRUE. )
                IF ( ASSOCIATED( Var1 ) ) THEN
                   Var1 % Valid = .FALSE.
                   Var1 % PrimaryMesh => PrimaryMesh
                END IF
              END DO
            END IF
          END IF
        END IF
      END IF
!     CALL InvalidateVariable( Mesh % Child, PrimaryMesh, Name )
      Mesh => Mesh % Next
    END DO 
!------------------------------------------------------------------------------
  END SUBROUTINE InvalidateVariable
!------------------------------------------------------------------------------
  END SUBROUTINE LoadRestartFile
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Writes data in ElmerPost format. 
!------------------------------------------------------------------------------
  SUBROUTINE WritePostFile( PostFile,ResultFile,Model,TimeCount,AppendFlag )
!------------------------------------------------------------------------------
    TYPE(Model_t), POINTER :: Model 
    INTEGER :: TimeCount
    LOGICAL, OPTIONAL :: AppendFlag
    CHARACTER(LEN=*) :: PostFile,ResultFile
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: CurrentElement
    TYPE(Variable_t), POINTER :: Var,Var1,Displacement,MeshUpdate,MaskVar

    CHARACTER(LEN=:), ALLOCATABLE :: Row
    CHARACTER(MAX_NAME_LEN) :: Str, DateStr

    LOGICAL :: gotIt, SaveCoordinates, MoveBoundary, MeshMoved, MaskExists

    REAL(KIND=dp) :: Time,Dummy, MeshScale, Coord(3)
    INTEGER :: ii,i,j,k,l,n,m,q,Node,idummy,DOFs,SavedCount,TimeStep, &
         NumberOfNodes, NumberOfElements, ind, nDOFs, MeshDim, Nzeros
    INTEGER, POINTER :: MaskPerm(:), MaskOrder(:)
!------------------------------------------------------------------------------

    IF( Model % Mesh % SavesDone == 0 ) THEN
      CALL Info('WritePostFile','Saving results in ElmerPost format to file '//TRIM(PostFile))
    END IF

    IF ( .NOT. FileNameQualified(PostFile) ) THEN
      IF ( LEN_TRIM(OutputPath) > 0 ) THEN
        IF ( AppendFlag .AND. Model % Mesh % SavesDone /= 0 )  THEN
          OPEN( PostFileUnit,File=TRIM(OutputPath) // '/' // &
             TRIM(PostFile), POSITION='APPEND' )
        ELSE
          OPEN( PostFileUnit,File=TRIM(OutputPath) // '/' // &
             TRIM(PostFile),STATUS='UNKNOWN' )
        END IF
      ELSE
        IF ( AppendFlag .AND. Model % Mesh % SavesDone /= 0 ) THEN
          OPEN( PostFileUnit,File=TRIM(PostFile),POSITION='APPEND' )
        ELSE
          OPEN( PostFileUnit,File=TRIM(PostFile),STATUS='UNKNOWN' )
        ENDIF
      END IF
    ELSE
      IF ( AppendFlag .AND. Model % Mesh % SavesDone /= 0  ) THEN
        OPEN( PostFileUnit,File=TRIM(PostFile),POSITION='APPEND' )
      ELSE
        OPEN( PostFileUnit,File=TRIM(PostFile),STATUS='UNKNOWN' )
      END IF
    END IF

    IF ( .NOT.AppendFlag ) THEN
      IF ( .NOT. FileNameQualified(ResultFile) ) THEN
        IF ( LEN_TRIM(OutputPath) > 0 ) THEN
          OPEN( OutputUnit,File=TRIM(OutputPath) // '/' &
              // TRIM(ResultFile),STATUS='OLD' )
        ELSE
          OPEN( OutputUnit,File=TRIM(ResultFile),STATUS='OLD' )
        END IF
      ELSE
        OPEN( OutputUnit,File=TRIM(ResultFile),STATUS='OLD' )
      END IF
    END IF

    SaveCoordinates = .FALSE.
    MoveBoundary    = .FALSE.
    DO i=1,CurrentModel % NumberOfBCs
      SaveCoordinates = ListGetLogical( &
         CurrentModel % BCs(i) % Values,'Free Surface', GotIt )
      IF ( SaveCoordinates ) THEN
         MoveBoundary =  ListGetLogical( &
             CurrentModel % BCs(i) % Values,'Internal Move Boundary', GotIt )         
         IF ( GotIt ) SaveCoordinates = MoveBoundary
      END IF

      IF ( SaveCoordinates ) EXIT
    END DO

!------------------------------------------------------------------------------
! Initialize stuff for masked saving
!------------------------------------------------------------------------------
    Str = ListGetString( Model % Simulation,'Post File Mask Variable',MaskExists)
    IF(.NOT. MaskExists) THEN
       Str = ListGetString( Model % Simulation,'ElmerPost Mask Variable',MaskExists)      
    END IF

    IF( MaskExists ) THEN
       MaskVar => VariableGet(Model % Variables,TRIM(Str))
       IF( ASSOCIATED(MaskVar)) MaskPerm => MaskVar % Perm
       MaskExists = ASSOCIATED(MaskPerm)
    END IF
    IF(MaskExists) THEN
       CALL Info('WritePostFile','Using '// TRIM(Str) // ' as mask variable')
       NumberOfNodes = MAXVAL(MaskPerm)
       ALLOCATE(MaskOrder(NumberOfNodes))
       DO i=1,SIZE(MaskPerm)
          j = MaskPerm(i)
          IF(j > 0) MaskOrder(j) = i
       END DO
       NumberOfElements = 0
       DO i=1,Model % NumberOfBulkElements + Model % NumberOfBoundaryElements
          CurrentElement => Model % Elements(i)
          IF( ALL(MaskPerm(CurrentElement % NodeIndexes) /= 0)) THEN
             NumberOfElements = NumberOfElements + 1
          END IF
       END DO
    ELSE
       NumberOfNodes = Model % NumberOfNodes
       NumberOfElements =  Model % NumberOfBulkElements + Model % NumberOfBoundaryElements
    END IF
 
!------------------------------------------------------------------------------
!   Count degrees of freedom to be saved
!------------------------------------------------------------------------------
    NULLIFY( Displacement, MeshUpdate )

    DOFs = 0
    Var => Model % Variables
    DO WHILE( ASSOCIATED(Var) )

      IF ( .NOT. Var % Output ) THEN
         Var => Var % Next; CYCLE
      END IF

      IF( SIZE( Var % Values ) == Var % Dofs ) THEN
        Var => Var % Next; CYCLE
      END IF

      IF( Var % TYPE /= Variable_on_nodes ) THEN
        Var => Var % Next; CYCLE
      END IF

      SELECT CASE(Var % Name(1:Var % NameLen))
        CASE( 'mesh update' )
           Var1 => Model % Variables
           DO WHILE( ASSOCIATED( Var1 ) )
             IF ( Var1 % Name == 'displacement' ) EXIT
             Var1 => Var1 % Next
           END DO
           IF ( .NOT. ASSOCIATED( Var1 ) ) THEN
              DOFs = DOFs + 3
           END IF

        CASE('mesh update 1','mesh update 2', 'mesh update 3' )

        CASE( 'displacement' )
          DOFs = DOFs + 3
          IF (ASSOCIATED(Var % Cvalues)) DOFs=DOFs+3

        CASE( 'displacement 1','displacement 2','displacement 3')

        CASE( 'flow solution' )
          DOFs = DOFs + 4

        CASE( 'velocity 1','velocity 2','velocity 3','pressure' )

        CASE( 'magnetic field' )
          DOFs = DOFs + 3

        CASE( 'magnetic field 1','magnetic field 2','magnetic field 3')

        CASE DEFAULT

          nDOFs = 1
          IF (ASSOCIATED(Var % Cvalues)) nDOFs=2

          IF ( Var % DOFs == 1 ) THEN
             DOFs = DOFs + nDOFs

          ELSE
            j=INDEX(Var % Name, '[')
            IF ( j > 0 ) THEN
              str =  ' '
              DO i=1,LEN_TRIM(Var % Name)
                str(i:i) = Var % Name(i:i)
              END DO
              DO WHILE( .TRUE. )
                i = INDEX( str(j+1:), ':' ) + j
                IF ( i<=j ) EXIT
                READ( str(i+1:),'(i1)' ) k
                DOFs = DOFs+nDOFs*k
                IF (k==2) DOFs=DOFs+nDOFs
                j = i + 1
              END DO
              DO k=1,Var % DOFs
                Var => Var % Next
              END DO
            END IF
          END IF

      END SELECT
      Var => Var % Next
    END DO

    IF ( .NOT.SaveCoordinates ) DOFs = DOFs-3

!------------------------------------------------------------------------------
! Write header to output
!------------------------------------------------------------------------------
    IF ( .NOT.AppendFlag .OR. Model % Mesh % SavesDone == 0 ) THEN
      WRITE(PostFileUnit,'(i10,i10,i7,i7)',ADVANCE='NO' ) NumberOfNodes, &
           NumberOfElements, DOFs, TimeCount

      Var => Model % Variables
      DO WHILE( ASSOCIATED( Var ) )

        IF ( .NOT. Var % Output ) THEN
           Var => Var % Next; CYCLE
        END IF

        IF( SIZE( Var % Values ) == Var % Dofs ) THEN
          Var => Var % Next; CYCLE
        END IF

        IF( Var % TYPE /= Variable_on_nodes ) THEN
          Var => Var % Next; CYCLE
        END IF

        SELECT CASE(Var % Name(1:Var % Namelen))
          CASE( 'mesh update' )

             Var1 => Model % Variables
             DO WHILE( ASSOCIATED( Var1 ) )
               IF ( TRIM(Var1 % Name) == 'displacement' ) EXIT
               Var1 => Var1 % Next
             END DO

             IF ( .NOT. ASSOCIATED( Var1 ) ) THEN
                WRITE(PostFileUnit,'(a)',ADVANCE='NO') ' vector: Mesh.Update'
                Displacement => Var
             ELSE
                MeshUpdate   => Var
             END IF

          CASE( 'mesh update 1','mesh update 2', 'mesh update 3' )

          CASE( 'displacement' )
            WRITE(PostFileUnit,'(a)',ADVANCE='NO' ) ' vector: Displacement'
            IF (ASSOCIATED(Var % Cvalues)) &
              WRITE(PostFileUnit,'(a)',ADVANCE='NO' ) ' vector: Displacement.im'
            Displacement => Var

          CASE( 'displacement 1','displacement 2','displacement 3')

          CASE( 'flow solution' )
            WRITE(PostFileUnit,'(a)',ADVANCE='NO') ' vector: Velocity scalar: Pressure'

          CASE( 'velocity 1','velocity 2','velocity 3','pressure' )

          CASE( 'magnetic field' )
            WRITE(PostFileUnit,'(a)',ADVANCE='NO' ) ' vector: MagField'

          CASE( 'magnetic field 1','magnetic field 2','magnetic field 3')

          CASE( 'coordinate 1','coordinate 2','coordinate 3' )

          CASE DEFAULT

            IF ( Var % DOFs == 1 ) THEN
              DO i=1,Var % NameLen
                str(i:i) = Var % Name(i:i)
                IF ( str(i:i) == ' ' ) str(i:i) = '.'
              END DO
              str(1:1) = CHAR(ICHAR(str(1:1))-ICHAR('a')+ICHAR('A'))
              WRITE(PostFileUnit,'(a,a)',ADVANCE='NO' ) ' scalar: ',str(1:Var % NameLen)
              IF (ASSOCIATED(Var % Cvalues)) &
                WRITE(PostFileUnit,'(a,a)',ADVANCE='NO' ) ' scalar: ',str(1:Var % NameLen)//'.im'
            ELSE
               j=INDEX(Var % Name, '[')
               IF ( j > 0 ) THEN
                 str =  ' '
                 DO i=1,LEN_TRIM(Var % Name)
                   str(i:i) = Var % Name(i:i)
                 END DO
                 DOFs = 0
                 DO WHILE( .TRUE. )
                   i = INDEX( str(j+1:), ':' ) + j
                   IF ( i<=j ) EXIT
                   READ( str(i+1:),'(i1)' ) k
                   DOFs = DOFs+k
                   DO WHILE( str(j+1:j+1) == ' ' )
                     j = j + 1
                   END DO
                   str(j+1:j+1) = CHAR(ICHAR(str(j+1:j+1))-ICHAR('a')+ICHAR('A'))
                   DO l=j+1,i-1
                     IF ( str(l:l) == ' ' ) str(l:l) = '.'
                   END DO
                   IF ( k==1 ) THEN
                     WRITE(PostFileUnit,'(a,a)',ADVANCE='NO' ) ' scalar: ',str(j+1:i-1)
                     IF (ASSOCIATED(Var % Cvalues)) &
                      WRITE(PostFileUnit,'(a,a)',ADVANCE='NO' ) ' scalar: ',str(j+1:i-1)//'.im'
                   ELSE IF ( k<= 3 ) THEN
                     WRITE(PostFileUnit,'(a,a)',ADVANCE='NO' ) ' vector: ',str(j+1:i-1)
                     IF (ASSOCIATED(Var % Cvalues)) &
                       WRITE(PostFileUnit,'(a,a)',ADVANCE='NO' ) ' vector: ',str(j+1:i-1)//'.im'
                   ELSE
                     DO l=1,k
                       WRITE(PostFileUnit,'(a,a)',ADVANCE='NO' ) ' scalar: ', &
                           str(j+1:i-1) // '.' // CHAR(l+ICHAR('0'))
                       IF (ASSOCIATED(Var % Cvalues)) &
                         WRITE(PostFileUnit,'(a,a)',ADVANCE='NO' ) ' scalar: ', &
                             str(j+1:i-1) // '.' // CHAR(l+ICHAR('0')) // '.im'
                     END DO
                   END IF
                   j = i + 1
                 END DO
                 DO k=1,Var % DOFs
                    Var => Var % Next
                 END DO
               END IF
            END IF
        END SELECT
        Var => Var % Next
      END DO

      IF ( SaveCoordinates ) THEN
        WRITE(PostFileUnit,'(a)',ADVANCE='NO' ) ' vector: Coordinates'
      END IF

      WRITE(PostFileUnit,'()')
      DateStr = FormatDate()
      WRITE( PostFileUnit, '("#File started at: ",A)' ) TRIM(DateStr)
!------------------------------------------------------------------------------
!   Coordinates
!------------------------------------------------------------------------------

      MeshScale = 1.0_dp
      DO i=1,Model % NumberOfSolvers
        MeshMoved = ListGetLogical( Model % Solvers(i) % Values, &
                    'Displace Mesh', Gotit )
        IF ( Gotit ) THEN
          IF ( .NOT. MeshMoved ) MeshScale = 0.0_dp
        ELSE
          MeshMoved = ListGetLogical( Model % Solvers(i) % Values, &
                 'Output Mesh Deformation', Gotit )
          IF ( GotIt ) THEN
            IF (MeshMoved ) MeshScale = 0.0_dp
          ELSE
            IF (Model % Solvers(i) % NofEigenValues>0) MeshScale=0._dp
          END IF
        END IF
      END DO

      MeshDim = Model % Mesh % MeshDim

      DO ii=1,NumberOfNodes

        i = ii
        IF(MaskExists) i = MaskOrder(i)
        
        Coord(1) = Model % Nodes % x(i)
        Coord(2) = Model % Nodes % y(i)
        Coord(3) = Model % Nodes % z(i)
        
        IF ( ASSOCIATED(Displacement) ) THEN
          k = Displacement % Perm(i)
          
          IF ( k > 0 ) THEN
            DO l=1,Displacement % Dofs
              Coord(l) = Coord(l) - MeshScale * &
                  Displacement % Values( Displacement % Dofs * (k-1) + l )
            END DO
          ELSE IF( ASSOCIATED( MeshUpdate ) ) THEN
            k = MeshUpdate % Perm(i)             
            IF ( k > 0 ) THEN
              DO l=1,MeshUpdate % Dofs
                Coord(l) = Coord(l) - MeshScale * &
                    MeshUpdate % Values( Displacement % Dofs * (k-1) + l )
              END DO
            END IF
          END IF
        END IF
        
        IF( MeshDim == 3 ) THEN
          WRITE(PostFileUnit,'(3ES17.8E3)') Coord(1:MeshDim)
        ELSE IF( MeshDim == 2 ) THEN
          WRITE(PostFileUnit,'(2ES17.8E3,A)') Coord(1:MeshDim),' 0.0'
        ELSE
          WRITE(PostFileUnit,'(ES17.8E3,A)') Coord(1:MeshDim),' 0.0 0.0'
        END IF
      END DO

!------------------------------------------------------------------------------
! Elements
!------------------------------------------------------------------------------
      WRITE(PostFileUnit,'(a)') '#group all'
      DO i=1,Model % NumberOfBulkElements
         CurrentElement => Model % Elements(i)

         IF(MaskExists) THEN
            IF( .NOT. ALL(MaskPerm(CurrentElement % NodeIndexes) /= 0)) CYCLE
         END IF
         
         k = CurrentElement % BodyId
         gotIt = .FALSE.
         IF ( k >= 1 .AND. k <= Model % NumberOfBodies ) THEN
            Str = ListGetString( Model % Bodies(k) % Values,'Name',gotIt )
         END IF
         
         IF ( gotIt ) THEN
            k = LEN_TRIM(Str)
            DO j=1,k
               IF ( Str(j:j) == ' ' ) Str(j:j) = '.'
            END DO
            
            WRITE( PostFileUnit,'(a)',ADVANCE='NO' )  Str(1:k)
         ELSE
            IF ( k > 0 .AND. k < 10 ) THEN
               WRITE(PostFileUnit,'(a,i1,a)',ADVANCE='NO' ) 'body',k,' '
            ELSE IF ( k >= 10 .AND. k < 100 ) THEN
               WRITE(PostFileUnit,'(a,i2,a)',ADVANCE='NO' ) 'body',k,' '
            ELSE
               WRITE(PostFileUnit,'(a,i3,a)',ADVANCE='NO' ) 'body',k,' '
            END IF
         END IF
         
         WRITE(PostFileUnit,'(i5)', ADVANCE='NO') CurrentElement % TYPE % ElementCode
         n = 0
         DO j=1,CurrentElement % TYPE % NumberOfNodes,4
            DO k=1,MIN(4,CurrentElement % TYPE % NumberOfNodes-n)
               n = n + 1
               ind = CurrentElement % NodeIndexes(n)
               IF(MaskExists) ind = MaskPerm(ind)
               WRITE(PostFileUnit, '(i8)', ADVANCE='NO')  ind - 1
            END DO
            WRITE( PostFileUnit,'(a)' ) ''
         END DO
      END DO
      
      DO i=Model % NumberOfBulkElements + 1,Model % NumberOfBulkElements + &
           Model % NumberOfBoundaryElements
         
         CurrentElement => Model % Elements(i)
         
         IF(MaskExists) THEN
            IF( .NOT. ALL(MaskPerm(CurrentElement % NodeIndexes) /= 0)) CYCLE
         END IF
         
         k = CurrentElement % BoundaryInfo % Constraint
         
         gotIt = .FALSE.
         IF ( k >= 1 .AND. k <= Model % NumberOfBCs ) THEN
            Str = ListGetString( Model % BCs(k) % Values,'Name',gotIt )
         END IF
         
         IF ( gotIt ) THEN
            k = LEN_TRIM(Str)
            DO j=1,k
               IF ( Str(j:j) == ' ' ) Str(j:j) = '.'
            END DO
            
            WRITE( PostFileUnit,'(a)',ADVANCE='NO' )  Str(1:k)
         ELSE
            IF ( k < 10 ) THEN
               WRITE( PostFileUnit,'(a,i1,a)',ADVANCE='NO' ) 'Constraint', k, ' '
            ELSE IF ( k < 100 ) THEN
               WRITE( PostFileUnit,'(a,i2,a)',ADVANCE='NO' ) 'Constraint', k, ' '
            ELSE
               WRITE( PostFileUnit,'(a,i3,a)',ADVANCE='NO' ) 'Constraint', k, ' '
            END IF
         END IF
         
         WRITE(PostFileUnit,'(i5)', ADVANCE='NO') CurrentElement % TYPE % ElementCode
         DO k=1,CurrentElement % TYPE % NumberOfNodes
            ind = CurrentElement % NodeIndexes(k)
            IF(MaskExists) ind = MaskPerm(ind)
            WRITE( PostFileUnit, '(i8)', ADVANCE='NO' )  ind-1
         END DO
         WRITE( PostFileUnit,'(a)' ) ''
      END DO
      WRITE(PostFileUnit,'(a)') '#endgroup all'
!------------------------------------------------------------------------------
!   Open result file and go trough it...
!------------------------------------------------------------------------------

      REWIND(OutputUnit)
   END IF ! .NOT.AppendFlag .OR. Model % Mesh % SavesDone == 0
   
   IF ( AppendFlag .AND. Model % Mesh % SavesDone == 0 ) THEN
      CLOSE(PostFileUnit)
      RETURN
   END IF

!------------------------------------------------------------------------------
   ALLOCATE(CHARACTER(MAX_STRING_LEN)::Row)
   DO WHILE( .TRUE. )
      IF ( AppendFlag ) THEN
         SavedCount = Model % Mesh % SavesDone
         TimeStep   = SavedCount
         Var => VariableGet( Model % Variables, 'Time', ThisOnly=.TRUE. )
         Time = 1.0d0
         IF ( ASSOCIATED(Var) ) Time = Var % Values(1)
      ELSE
!------------------------------------------------------------------------------
!   ...read one timestep to memory (if not already there)...
!------------------------------------------------------------------------------
         DO WHILE( ReadAndTrim(OutputUnit,Row) )
            IF ( SEQL(Row, 'total dofs:') ) READ( Row(12:),* ) DOFs
            IF ( SEQL( Row, 'time:') ) EXIT
         END DO
         
         IF ( .NOT.SEQL(Row, 'time:') ) EXIT
         
         READ( Row(7:),* ) SavedCount,Timestep,Time
      END IF
      
      WRITE( PostFileUnit,'(a,i7,i7,ES17.8E3)' ) '#time ',SavedCount,Timestep,Time
      
      IF ( .NOT.AppendFlag ) THEN
         DO i=1,DOFs
            READ(OutputUnit,'(a)' ) Row
            Var => VariableGet( Model % Variables,Row,.TRUE. )
            
            IF ( ASSOCIATED(Var) ) THEN
               DO j=1,NumberOfNodes
                  k = j
                  IF( MaskExists ) k = MaskPerm(k)
                  IF ( ASSOCIATED(Var % Perm) ) k = Var % Perm(k)
                  IF ( k > 0 ) THEN
                     READ(OutputUnit,*) Node,idummy,Var % Values(k)
                  ELSE
                     READ(OutputUnit,*) Node,idummy,Dummy
                  END IF
               END DO
            END IF
        END DO
      END IF
!-----------------------------------------------------------------------------
!     ...then save it to post file.
!------------------------------------------------------------------------------
     DO ii=1,NumberOfNodes
        
        i = ii
        IF(MaskExists) i = MaskOrder(i)

        Var => Model % Variables
        DO WHILE( ASSOCIATED( Var ) )
           IF ( .NOT. Var % Output ) THEN
              Var => Var % Next; CYCLE
           END IF
           
           IF( SIZE( Var % Values ) == Var % DOFs ) THEN
             Var => Var % Next; CYCLE
           END IF

           IF( Var % TYPE /= Variable_on_nodes ) THEN
             Var => Var % Next; CYCLE
           END IF

           SELECT CASE(Var % Name(1:Var % Namelen))
              
           CASE( 'mesh update' )
              Var1 => Model % Variables
              DO WHILE( ASSOCIATED( Var1 ) )
                 IF ( Var1 % Name(1:Var1 % NameLen) == 'displacement' ) EXIT
                 Var1 => Var1 % Next
              END DO
              IF ( .NOT. ASSOCIATED( Var1 ) ) THEN
                 k = i
                 IF ( ASSOCIATED(Var % Perm) ) k = Var % Perm(k)
                 IF ( k > 0 ) THEN
                    DO j=1,Var % DOFs
                       WRITE(PostFileUnit,'(ES17.8E3)',ADVANCE='NO') Var % Values(Var % DOFs*(k-1)+j)
                    END DO
                    IF ( Var % DOFs == 2 ) WRITE(PostFileUnit,'(A)',ADVANCE='NO') ' 0.0'
                 ELSE
                    WRITE(PostFileUnit,'(A)',ADVANCE='NO') ' 0.0 0.0 0.0'
                 END IF
              END IF
              
           CASE(  'mesh update 1','mesh update 2', 'mesh update 3' )
              
           CASE( 'displacement' )
              k = i
              IF ( ASSOCIATED(Var % Perm) ) k = Var % Perm(k)
              
              IF ( k > 0 ) THEN
                 IF (ASSOCIATED(Var % Cvalues)) THEN
                   DO j=1,Var % DOFs
                     WRITE(PostFileUnit,'(ES17.8E3)',ADVANCE='NO') &
                         REAL(Var % Cvalues(Var % DOFs*(k-1)+j))
                   END DO
                   IF (Var % DOFs == 2) WRITE(PostFileUnit,'(A)',ADVANCE='NO') ' 0.0'
                   DO j=1,Var % DOFs
                     WRITE(PostFileUnit,'(ES17.8E3)',ADVANCE='NO') &
                         AIMAG(Var % Cvalues(Var % DOFs*(k-1)+j))
                   END DO
                   IF (Var % DOFs == 2) WRITE(PostFileUnit,'(A)',ADVANCE='NO') ' 0.0'
                 ELSE
                   DO j=1,Var % DOFs
                     WRITE(PostFileUnit,'(ES17.8E3)',ADVANCE='NO') &
                         Var % Values(Var % DOFs*(k-1)+j)
                   END DO
                   IF (Var % DOFs == 2) WRITE(PostFileUnit,'(A)',ADVANCE='NO') ' 0.0'
                 END IF
              ELSE
                 Var1 => Model % Variables
                 DO WHILE( ASSOCIATED( Var1 ) )
                    IF ( Var1 % Name(1:Var1 % Namelen) == 'mesh update' ) EXIT
                    Var1 => Var1 % Next
                 END DO
                 IF ( ASSOCIATED( Var1 ) ) THEN
                    k = i
                    IF ( ASSOCIATED(Var1 % Perm) ) k = Var1 % Perm(k)
                    IF ( k > 0 ) THEN
                      DO j=1,Var1 % DOFs
                        WRITE(PostFileUnit,'(ES17.8E3)',ADVANCE='NO')  &
                            Var1 % Values(Var1 % DOFs*(k-1)+j)
                      END DO
                      IF ( Var1 % DOFs == 2 ) WRITE(PostFileUnit,'(A)',ADVANCE='NO') ' 0.0'
                    ELSE
                      WRITE(PostFileUnit,'(A)',ADVANCE='NO') ' 0.0 0.0 0.0'
                    END IF
                 ELSE
                   IF( ASSOCIATED( Var % Cvalues ) ) THEN
                     WRITE(PostFileUnit,'(A)',ADVANCE='NO') ' 0.0 0.0 0.0 0.0 0.0 0.0'
                   ELSE
                     WRITE(PostFileUnit,'(A)',ADVANCE='NO') ' 0.0 0.0 0.0'
                   END IF
                 END IF
              END IF
              
           CASE( 'displacement 1','displacement 2','displacement 3')
              
           CASE( 'flow solution' )
              k = i
              IF ( ASSOCIATED(Var % Perm) ) k = Var % Perm(k)
              IF ( k > 0 ) THEN
                 DO j=1,Var % DOFs-1
                    WRITE(PostFileUnit,'(ES17.8E3)',ADVANCE='NO') Var % Values(Var % DOFs*(k-1)+j)
                 END DO
                 IF ( Var % DOFs < 4 ) WRITE(PostFileUnit,'(A)',ADVANCE='NO') ' 0.0'
                 
                 WRITE(PostFileUnit,'(ES17.8E3)',ADVANCE='NO')  &
                      Var % Values(Var % DOFs*k)
              ELSE
                 WRITE(PostFileUnit,'(A)',ADVANCE='NO') ' 0.0 0.0 0.0 0.0'
              END IF
              
           CASE( 'velocity 1','velocity 2','velocity 3','pressure' )
              
           CASE( 'magnetic field' )
              k = i
              IF ( ASSOCIATED(Var % Perm) ) k = Var % Perm(k)
              IF ( k > 0 ) THEN
                 DO j=1,Var % DOFs
                    WRITE(PostFileUnit,'(ES17.8E3)',ADVANCE='NO') Var % Values(Var % DOFs*(k-1)+j)
                 END DO
                 IF ( Var % DOFs == 2 ) WRITE(PostFileUnit,'(A)',ADVANCE='NO') ' 0.0'
              ELSE
                 WRITE(PostFileUnit,'(A)',ADVANCE='NO') ' 0.0 0.0 0.0'
              END IF
              
           CASE( 'magnetic field 1','magnetic field 2','magnetic field 3')
              
           CASE( 'coordinate 1','coordinate 2','coordinate 3' )
              
           CASE DEFAULT

              IF ( Var % DOFs == 1 ) THEN
                 k = i
                 IF ( ASSOCIATED(Var % Perm) ) k = Var % Perm(k)
                 IF ( k > 0 ) THEN
                    IF (ASSOCIATED(Var % Cvalues)) THEN
                      WRITE(PostFileUnit,'(ES17.8E3)',ADVANCE='NO') &
                         REAL(Var % Cvalues(k))
                      WRITE(PostFileUnit,'(ES17.8E3)',ADVANCE='NO') &
                         AIMAG(Var % Cvalues(k))
                    ELSE
                      IF(k<=SIZE(Var % Values)) THEN
                        WRITE(PostFileUnit,'(ES17.8E3)',ADVANCE='NO') Var % Values(k)
                      ELSE
                        WRITE(PostFileUnit,'(ES17.8E3)',ADVANCE='NO') 0._dp
                      END IF
                    END IF
                 ELSE
                    WRITE(PostFileUnit,'(A)',ADVANCE='NO') ' 0.0'
                    IF(ASSOCIATED(Var % Cvalues)) &
                      WRITE(PostFileUnit,'(A)',ADVANCE='NO') ' 0.0'
                 END IF
              ELSE
                 l = INDEX( var % name, '[' )
                 IF ( l > 0 ) THEN
                    DOFs = 0
                    DO WHILE( .TRUE. )
                       m = INDEX( Var % Name(l+1:), ':' ) + l
                       IF ( m<=l ) EXIT
                       READ( Var % Name(m+1:),'(i1)' ) q
                       k = i
                       IF ( ASSOCIATED(Var % Perm) ) k = Var % Perm(k)
                       IF ( k > 0 ) THEN
                          IF ( q==2 .OR. q==3 ) THEN
                            IF (ASSOCIATED(Var % Cvalues)) THEN
                               DO j=DOFs+1,DOFs+q
                                  WRITE(PostFileUnit,'(ES17.8E3)',ADVANCE='NO') &
                                     REAL(Var % Cvalues(Var % DOFs*(k-1)+j))
                               END DO
                               IF(q==2) WRITE(PostFileUnit,'(A)',ADVANCE='NO') ' 0.0'
                               DO j=DOFs+1,DOFs+q
                                  WRITE(PostFileUnit,'(ES17.8E3)',ADVANCE='NO') &
                                     AIMAG(Var % Cvalues(Var % DOFs*(k-1)+j))
                               END DO
                               IF(q==2) WRITE(PostFileUnit,'(A)',ADVANCE='NO') ' 0.0'
                            ELSE
                               DO j=DOFs+1,DOFs+q
                                 WRITE(PostFileUnit,'(ES17.8E3)',ADVANCE='NO') &
                                     Var % Values(Var % DOFs*(k-1)+j)
                               END DO
                               IF(q==2) WRITE(PostFileUnit,'(A)',ADVANCE='NO') ' 0.0'
                            END IF
                          ELSE
                            DO j=DOFs+1,DOFs+q
                              IF (ASSOCIATED(Var % Cvalues)) THEN
                                WRITE(PostFileUnit,'(ES17.8E3)',ADVANCE='NO') &
                                   REAL(Var % Cvalues(Var % DOFs*(k-1)+j))
                                WRITE(PostFileUnit,'(ES17.8E3)',ADVANCE='NO') &
                                   AIMAG(Var % Cvalues(Var % DOFs*(k-1)+j))
                              ELSE
                                WRITE(PostFileUnit,'(ES17.8E3)',ADVANCE='NO') &
                                   Var % Values(Var % DOFs*(k-1)+j)
                              END IF
                            END DO
                          END IF
                       ELSE
                         IF( q == 2 ) THEN
                           Nzeros = 3
                         ELSE
                           Nzeros = q
                         END IF
                         IF( ASSOCIATED( Var % CValues ) ) Nzeros = 2 * Nzeros
                         DO j=1,Nzeros
                             WRITE(PostFileUnit,'(A)',ADVANCE='NO') ' 0.0'
                          END DO
                       END IF
                       l=m+1
                       DOFs = DOFs+q
                    END DO
                    DO j=1,Var % DOFs
                       Var => Var % Next
                    END DO
                 END IF
              END IF
           END SELECT
           Var => Var % Next
        END DO

        IF ( SaveCoordinates ) THEN
           Var => VariableGet( Model % Variables,'Coordinate 1' )
           WRITE(PostFileUnit,'(ES17.8E3)',ADVANCE='NO') Var % Values(i)
           
           Var => VariableGet( Model % Variables,'Coordinate 2' )
           WRITE(PostFileUnit,'(ES17.8E3)',ADVANCE='NO') Var % Values(i)
           
           Var => VariableGet( Model % Variables,'Coordinate 3' )
           WRITE(PostFileUnit,'(ES17.8E3)',ADVANCE='NO') Var % Values(i)
        END IF
        
        WRITE(PostFileUnit,'()')
     END DO
     IF (  AppendFlag ) EXIT
!------------------------------------------------------------------------------
  END DO
!------------------------------------------------------------------------------
!   We are done here close the files and deallocate
!------------------------------------------------------------------------------
  CLOSE(PostFileUnit)
  IF ( .NOT. AppendFlag ) CLOSE(OutputUnit)

  IF(MaskExists) DEALLOCATE(MaskOrder)
  
END SUBROUTINE WritePostFile
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>  Solve the ElementSize equation using Galerkin discretization where 
!> at each integration point the element size is the size of the element 
!> in question. The nodal elementsize may be used when defining the 
!> new mesh density after inheritance.
!------------------------------------------------------------------------------
SUBROUTINE GetNodalElementSize(Model,expo,noweight,h)
!------------------------------------------------------------------------------
  USE IterSolve
  USE Integration
  USE ElementDescription

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  LOGICAL :: noweight
  REAL(KIND=dp) :: expo
  REAL(KIND=Dp), POINTER :: h(:)
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Solver_t), POINTER :: Solver
  REAL(KIND=dp) :: dt=1._dp
  LOGICAL :: TransientSimulation=.FALSE.

  TYPE(Element_t),POINTER :: Element

  LOGICAL :: AllocationsDone = .FALSE., Found

  INTEGER :: n, t, istat, active
  REAL(KIND=dp) :: Norm, Power
  TYPE(Mesh_t), POINTER :: Mesh
  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:)
  LOGICAL :: GotIt
  REAL(KIND=dp) :: ElemMin, ElemMax
  TYPE(ValueList_t), POINTER :: Params

  TYPE(Matrix_t), POINTER :: A
  INTEGER, POINTER :: CPerm(:)

  SAVE STIFF, FORCE, AllocationsDone
!------------------------------------------------------------------------------

  CALL Info('GetNodalElementSize','Computing nodal element size indicator')

  ALLOCATE( Solver )
  Solver % Mesh => Model % Mesh
  Mesh => Solver % Mesh

  ALLOCATE(Solver % Def_Dofs(10,Model % NumberOfBodies,6))
  Solver % Def_Dofs = -1; Solver % Def_Dofs(:,:,1)=1

  Solver % Values => ListAllocate()
  Params => Solver % Values
  CALL ListAddString( Params,'Linear System Iterative Method', 'CG' )
  CALL ListAddLogical( Params,'Linear System Symmetric', .TRUE. )
  CALL ListAddInteger( Params, 'Linear System Max Iterations', 5000 )
  CALL ListAddString( Params, 'Linear System Preconditioning', 'ILU0' )
  CALL ListAddInteger( Params, 'Linear System Residual Output', 1 )
  CALL ListAddConstReal( Params, 'Linear System Convergence Tolerance', 1.0d-9 )

  ALLOCATE(CPerm(Mesh % NumberOfNodes+Mesh % NumberOfEdges))

  A => CreateMatrix(Model,Solver, &
           Solver % Mesh,CPerm,1,MATRIX_CRS,.FALSE.,NodalDofsOnly=.TRUE. )

  A % Comm = ELMER_COMM_WORLD
  A % ParMatrix => NULL()
  Solver % Matrix => A
  Model % Solver => Solver

  ALLOCATE(A % RHS(Mesh % NumberOfNodes))

  Solver % TimeOrder = 0

  Mesh % Variables => NULL()
  CALL VariableAdd( Mesh % Variables, Mesh, Solver, &
        'nodal h',1,h,Cperm) !,Output=.FALSE.)
  Solver % Variable=>VariableGet(Mesh % Variables,'nodal h',ThisOnly=.TRUE.)

  IF ( ParEnv % PEs>1 ) THEN
    IF ( ASSOCIATED(Solver % Mesh % ParallelInfo % INTERFACE) ) THEN
      ParEnv % ActiveComm = ELMER_COMM_WORLD

      ALLOCATE(ParEnv % Active(ParEnv % PEs))
      ParEnv % Active=.TRUE.

      CALL ParallelInitMatrix(Solver, Solver % Matrix )

      Solver % Matrix % ParMatrix % ParEnv % ActiveComm = &
                 Solver % Matrix % Comm
      ParEnv = Solver % Matrix % ParMatrix % ParEnv
    END IF
  END IF


  ! Allocate some storage
  !--------------------------------------------------------------
  N = Mesh % MaxElementNodes ! just big enough for elemental arrays
  ALLOCATE( FORCE(N), STIFF(N,N), STAT=istat )
  IF ( istat /= 0 ) THEN
    CALL Fatal( 'GetNodalElementSize', 'Memory allocation error.' )
  END IF
  
  ElemMin =  HUGE(ElemMin)
  ElemMax = -HUGE(ElemMax)

  !Initialize the system and do the assembly:
  !------------------------------------------
  Power = 1.0_dp / expo
  active = Mesh % NumberOfBulkElements

  A % RHS=0._dp
  A % Values=0._dp
  DO t=1,active
    Element => Mesh % Elements(t)
    Model % CurrentElement => Element
    n = Element % TYPE % NumberOfNodes
    
    !Get element local matrix and rhs vector:
    !----------------------------------------
    CALL LocalMatrix(  STIFF, FORCE, Element, n )
    
    !Update global matrix and rhs vector from local matrix & vector:
    !---------------------------------------------------------------
    CALL CRS_GlueLocalMatrix( A,n,1,Element % NodeIndexes,STIFF )
    A % RHS(Element % NodeIndexes) = &
             A % RHS(Element % NodeIndexes)+FORCE(1:n)
  END DO

  h=0
  IF (ASSOCIATED(A % ParMatrix)) THEN
    CALL ParallelIter(A,A % ParallelInfo,1,h,A % RHS,Solver,A % ParMatrix)
  ELSE
    CALL IterSolver(A,h,A % RHS,Solver)
  END IF

  WRITE(Message,'(A,2ES12.4)') 'Minimum Element Size: ',ElemMin, MINVAL(h)
  CALL Info('GetNodalElementSize',Message)
  WRITE(Message,'(A,2ES12.4)') 'Maximum Element Size: ',ElemMax, MAXVAL(h)
  CALL Info('GetNodalElementSize',Message)
  WRITE(Message,'(A,2ES12.4)') 'Element Size Ratio: ',ElemMax / ElemMin, MAXVAL(h)/MINVAL(h)
  CALL Info('GetNodalElementSize',Message)

  Model % Solver=>NULL()
  DEALLOCATE(Mesh % Variables)
  Mesh % Variables => NULL()

  CALL FreeSolver(Solver)
  DEALLOCATE( FORCE, STIFF,Cperm )
 
CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  STIFF, FORCE, Element, n )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:)
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),DetJ,LoadAtIP,Weight
    LOGICAL :: Stat
    INTEGER :: i,j,t
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------

    STIFF = 0.0_dp
    FORCE = 0.0_dp

    ALLOCATE(Nodes % x(n), Nodes % y(n), Nodes % z(n))
    Nodes % x = Mesh  % Nodes % x(Element % NodeIndexes)
    Nodes % y = Mesh  % Nodes % y(Element % NodeIndexes)
    Nodes % z = Mesh  % Nodes % z(Element % NodeIndexes)

    ! Numerical integration:
    !----------------------
    IP = GaussPoints( Element )

    DO t=1,IP % n
       ! Basis function values & derivatives at the integration point:
       !--------------------------------------------------------------
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
               IP % W(t),  detJ, Basis )

       ! The source term at the integration point:
       !------------------------------------------
       LoadAtIP = DetJ ** Power
       IF( NoWeight ) THEN
         Weight = IP % s(t) 
       ELSE
         Weight = IP % s(t) * DetJ
       END IF

       ElemMin = MIN( ElemMin, LoadAtIP )
       ElemMax = MAX( ElemMax, LoadAtIP )
       
       ! Finally, the elemental matrix & vector:
       !----------------------------------------
       DO i = 1, n
         DO j = 1, n
           STIFF(i,j) = STIFF(i,j) + Weight * &
               Basis(i) * Basis(j)
         END DO
         FORCE(i) = FORCE(i) + Weight * Basis(i) * LoadAtIp
       END DO
    END DO
    DEALLOCATE(Nodes % x, Nodes % y, Nodes % z)
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------
END SUBROUTINE GetNodalElementSize
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Release a mesh from the list of meshes.
!------------------------------------------------------------------------------
  RECURSIVE SUBROUTINE FreeMesh(Mesh)
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
!------------------------------------------------------------------------------
    IF (.NOT.ASSOCIATED(Mesh)) RETURN

    CALL FreeMesh(Mesh % Next)

    Mesh % Next   => NULL()
    Mesh % Child  => NULL()
    Mesh % Parent => NULL()

    CALL ReleaseMesh(Mesh)
    DEALLOCATE(Mesh)
!------------------------------------------------------------------------------
  END SUBROUTINE FreeMesh
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Releases structures related to the Solver. 
!------------------------------------------------------------------------------
  SUBROUTINE FreeSolver(Solver)
!------------------------------------------------------------------------------
    TYPE(Solver_t) :: Solver
!------------------------------------------------------------------------------

    CALL Info('FreeSolver','Free solver matrix',Level=20)
    CALL FreeMatrix(Solver % Matrix)

    CALL Info('FreeSolver','Free solver miscallenous',Level=20)
    CALL FreeValueList(Solver % Values)
    IF (ALLOCATED(Solver % Def_Dofs)) DEALLOCATE(Solver % Def_Dofs)
    IF (ASSOCIATED(Solver % ActiveElements)) DEALLOCATE(Solver % ActiveElements)
    IF( ASSOCIATED( Solver % ColourIndexList ) ) THEN
      CALL Graph_Deallocate(Solver % ColourIndexList)
      DEALLOCATE( Solver % ColourIndexList )
    END IF
        
!------------------------------------------------------------------------------
  END SUBROUTINE FreeSolver
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Releases value list which includes all the sif definitions, for example.
!------------------------------------------------------------------------------
  SUBROUTINE FreeValueList(List)
!------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: List
!------------------------------------------------------------------------------
    TYPE(ValueListEntry_t), POINTER :: ptr
   
    IF(.NOT.ASSOCIATED(List)) RETURN
    ptr => List % Head
    DO WHILE(ASSOCIATED(ptr))
      IF (ASSOCIATED(ptr % TValues)) DEALLOCATE(ptr % TValues)
      IF (ASSOCIATED(ptr % FValues)) DEALLOCATE(ptr % FValues)
      IF (ASSOCIATED(ptr % IValues)) DEALLOCATE(ptr % IValues)
      ptr => ptr % Next
    END DO 
    DEALLOCATE(List)
!------------------------------------------------------------------------------
  END SUBROUTINE FreeValueList
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Releases the whole model. 
!------------------------------------------------------------------------------
 SUBROUTINE FreeModel(Model)
!------------------------------------------------------------------------------
   TYPE(Model_t), POINTER :: Model
!------------------------------------------------------------------------------
   TYPE(Matrix_t), POINTER :: A,B
   INTEGER :: i
   IF (.NOT.ASSOCIATED(Model)) RETURN

   CALL Info('FreeModel','Freeing meshes',Level=15)
   CALL FreeMesh(Model % Meshes)

   CALL Info('FreeModel','Freeing constants list',Level=15)
   CALL FreeValueList(Model % Constants)

   CALL Info('FreeModel','Freeing simulation list',Level=15)
   CALL FreeValueList(Model % Simulation)

   IF (ASSOCIATED(Model % BCs)) THEN
     CALL Info('FreeModel','Freeing boundary lists',Level=15)
     DO i=1,Model % NumberOfBCs
#if 0
       A => Model % BCs(i) % PMatrix
       IF (ASSOCIATED(A)) THEN
         DO WHILE( ASSOCIATED(A) )
           B => A % Child
           A % Child => NULL()
           A => B
         END DO
         CALL FreeMatrix(Model % BCs(i) % PMatrix)
       END IF
#endif
       CALL FreeValueList( Model % BCs(i) % Values)
     END DO
     DEALLOCATE(Model % BCs)
   END IF

   CALL Info('FreeModel','Freeing solvers',Level=15)  
   DO i=1,Model % NumberOfSolvers
     CALL Info('FreeModel','Solver: '//TRIM(I2S(i)),Level=20)
     CALL FreeSolver(Model % Solvers(i))
   END DO
   DEALLOCATE(Model % Solvers)

   IF (ASSOCIATED(Model % ICs)) THEN
     CALL Info('FreeModel','Freeing initial conditions lists',Level=15)   
     DO i=1,Model % NumberOfICs
       CALL FreeValueList( Model % ICs(i) % Values)
     END DO
     DEALLOCATE(Model % ICs)
   END IF

   IF (ASSOCIATED(Model % Bodies)) THEN
     CALL Info('FreeModel','Freeing body lists',Level=15)   
     DO i=1,Model % NumberOfBodies
       CALL FreeValueList( Model % Bodies(i) % Values)
     END DO
     DEALLOCATE(Model % Bodies)
   END IF

   IF (ASSOCIATED(Model % Equations)) THEN
     CALL Info('FreeModel','Freeing equations lists',Level=15)    
     DO i=1,Model % NumberOfEquations
       CALL FreeValueList( Model % Equations(i) % Values)
     END DO
     DEALLOCATE(Model % Equations)
   END IF

   IF (ASSOCIATED(Model % BodyForces)) THEN
     CALL Info('FreeModel','Freeing body forces lists',Level=15)   
     DO i=1,Model % NumberOfBodyForces
       CALL FreeValueList( Model % BodyForces(i) % Values)
     END DO
     DEALLOCATE(Model % BodyForces)
   END IF

   Model=>NULL()
!------------------------------------------------------------------------------
 END SUBROUTINE FreeModel
!------------------------------------------------------------------------------


END MODULE ModelDescription

!> \}
