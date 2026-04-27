!/*****************************************************************************/
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

!------------------------------------------------------------------------------
 MODULE VectorHelmholtzUtils 
!------------------------------------------------------------------------------
   USE Types
   USE Lists
   USE ElementUtils, ONLY : SetParentBasis
   USE ElementDescription
   USE ParallelUtils
   IMPLICIT NONE

   COMPLEX(KIND=dp), PARAMETER :: im = (0._dp,1._dp)   
   
 CONTAINS


!------------------------------------------------------------------------------
   FUNCTION ComplexCrossProduct(v1,v2) RESULT(v3)
!------------------------------------------------------------------------------
     COMPLEX(KIND=dp) :: v1(3), v2(3), v3(3)
     v3(1) =  v1(2)*v2(3) - v1(3)*v2(2)
     v3(2) = -v1(1)*v2(3) + v1(3)*v2(1)
     v3(3) =  v1(1)*v2(2) - v1(2)*v2(1)
!------------------------------------------------------------------------------
   END FUNCTION ComplexCrossProduct
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> This routine computes geometric quantities for lumped port making the setting
!> of BCs a little easier.
!------------------------------------------------------------------------------
  SUBROUTINE DefinePortParameters(Model, Mesh)
!------------------------------------------------------------------------------
    IMPLICIT NONE 
    TYPE(Model_t) :: Model 
    TYPE(Mesh_t), POINTER :: Mesh
!------------------------------------------------------------------------------
    TYPE(Nodes_t) :: ElementNodes
    TYPE(Element_t), POINTER :: Element
    TYPE(GaussIntegrationPoints_t) :: IP
    INTEGER ::i,j,k,n,t,t0,bc_id,ierr,PortDir,PortIndex,PortTypeInd
    INTEGER, POINTER :: Indexes(:)
    REAL(KIND=dp) :: s,detJ,Width,Area,Length,RadInner, RadOuter,CenterArray(3,1),Scale
    REAL(KIND=dp), ALLOCATABLE :: Basis(:), LumpVec(:)
    TYPE(ValueList_t), POINTER :: BC, BC0
    CHARACTER(:), ALLOCATABLE :: PortType
    LOGICAL :: Found, stat
    CHARACTER(*), PARAMETER :: Caller = 'DefinePortParameters'

    IF(.NOT. ASSOCIATED( Mesh ) ) THEN
      CALL Fatal(Caller,'Mesh not associated!')
    END IF

    CALL Info(Caller,'Defining geometric port parameters',Level=8)

    n = Mesh % MaxElementNodes    
    t0 = Mesh % NumberOfBulkElements       
    
    ! Check the number of ports and add "port type index" keyword.
    DO bc_id = 1,Model % NumberOfBCs
      BC => Model % BCs(bc_id) % Values

      ! This has already been defined!
      IF( ListCheckPresent( BC,'Port Type Index')) CYCLE

      PortType = ListGetString( BC,'port type',Found)

      IF(.NOT. (Found .OR. ListCheckPresent(BC,'port impedance'))) CYCLE            

      SELECT CASE(PortType)
      CASE('rectangular')
        PortTypeInd = 1

      CASE('coaxial')
        PortTypeInd = 2

      CASE('eigenmode')
        PortTypeInd = 3

      CASE('potential')
        PortTypeInd = 4

      CASE('beta')
        PortTypeInd = 5

      CASE DEFAULT
        CALL Info(Caller,'Port Type "Port Type" defaulted to "rectangular"',Level=4)
        PortTypeInd = 1
      END SELECT
      
      CALL Info(Caller,'Defining parameters for port on BC: '//I2S(bc_id),Level=8)
      CALL ListAddInteger( BC,'Port Type Index',PortTypeInd)

      IF(.NOT. ALLOCATED(Basis) ) THEN
        ALLOCATE(Basis(n), ElementNodes % x(n), ElementNodes % y(n), ElementNodes % z(n), LumpVec(7) )
      END IF
      
      LumpVec(1:3) = HUGE(s)
      LumpVec(4:6) = -HUGE(s)
      LumpVec(7) = 0.0_dp

      
      DO t=1, Mesh % NumberOfBoundaryElements
        Element => Mesh % Elements( t0 + t )               
        IF( Element % BoundaryInfo % Constraint /= Model % BCs(bc_id) % Tag ) CYCLE

        Indexes => Element % NodeIndexes
        n = Element % TYPE % NumberOfNodes

        ElementNodes % x = 0.0_dp
        ElementNodes % y = 0.0_dp
        ElementNodes % z = 0.0_dp

        ElementNodes % x(1:n) = Mesh % Nodes % x(Indexes)
        ElementNodes % y(1:n) = Mesh % Nodes % y(Indexes)
        ElementNodes % z(1:n) = Mesh % Nodes % z(Indexes)

        ! Get min/max range for each coordinate.
        LumpVec(1) = MIN(LumpVec(1),MINVAL(ElementNodes % x(1:n)))
        LumpVec(2) = MIN(LumpVec(2),MINVAL(ElementNodes % y(1:n)))
        LumpVec(3) = MIN(LumpVec(3),MINVAL(ElementNodes % z(1:n)))
        LumpVec(4) = MAX(LumpVec(4),MAXVAL(ElementNodes % x(1:n)))
        LumpVec(5) = MAX(LumpVec(5),MAXVAL(ElementNodes % y(1:n)))
        LumpVec(6) = MAX(LumpVec(6),MAXVAL(ElementNodes % z(1:n)))
        
        ! Integrate over the area.
        IP = GaussPoints( Element, PReferenceElement = .FALSE.)
        DO j=1,IP % n        
          stat = ElementInfo( Element, ElementNodes, IP % U(j), IP % V(j), IP % W(j), detJ, Basis )
          S = DetJ * IP % s(j)
          LumpVec(7) = LumpVec(7) + s
        END DO
      END DO
      
      ! Do parallel communication, if needed.
      IF( ParEnv % PEs > 1 ) THEN
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, LumpVec(1:3), 3, &
            MPI_DOUBLE_PRECISION, MPI_MIN, ELMER_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, LumpVec(4:6), 3, &
            MPI_DOUBLE_PRECISION, MPI_MAX, ELMER_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, LumpVec(7:7), 1, &
            MPI_DOUBLE_PRECISION, MPI_SUM, ELMER_COMM_WORLD, ierr )
      END IF

      ! Area is used in all port models. 
      Area = LumpVec(7)
      
      SELECT CASE(PortTypeInd)
      CASE(1)
        PortDir = ABS( ListGetInteger( BC,'Port Direction',Found) )
        IF(.NOT. Found) PortDir = 3
                
        Length = LumpVec(3+PortDir) - LumpVec(PortDir)        
        Width = Area / Length
        Scale = Width / Length

        !PRINT *,'area:',area, length, width, scale
        
        CALL ListAddConstReal( BC,'Port Length',Length)
        CALL ListAddConstReal( BC,'Port Scale',Scale)
        IF(InfoActive(8)) THEN
          PRINT *,'Setting rectangular port parameters:',Length,Scale
        END IF

      CASE(2)
        RadOuter = 0.0_dp
        DO i=1,3
          RadOuter = MAX(RadOuter,(LumpVec(3+i)-LumpVec(i))/2)
          CenterArray(i,1) = (LumpVec(3+i)+LumpVec(i))/2
        END DO
        RadInner = SQRT(RadOuter**2-Area/PI)        
        Length = (RadInner+RadOuter)/2
        Scale = 2*PI/LOG(RadOuter/RadInner)

        ! PRINT *,'area:',area, radinner, radouter
        
        CALL ListAddConstReal( BC,'Port Length',Length) 
        CALL ListAddConstReal( BC,'Port Scale',Scale)
        CALL ListAddConstRealArray( BC,'Port Center',3,1,CenterArray)
        IF(InfoActive(8)) THEN
          PRINT *,'Setting coaxial port parameters:',Length,Scale,' and center ',CenterArray
        END IF
      END SELECT
          
    END DO

    IF(ALLOCATED(Basis)) THEN
      DEALLOCATE(Basis, ElementNodes % x, ElementNodes % y, ElementNodes % z, LumpVec )
    END IF

  END SUBROUTINE DefinePortParameters
!------------------------------------------------------------------------------

   
!------------------------------------------------------------------------------
  SUBROUTINE ElectricPortModel(Phase,Solver,Element,GotPort,&
      B,L,Basis,dBasisdx,WBasis)
!------------------------------------------------------------------------------
     INTEGER :: Phase
     TYPE(Solver_t) :: Solver
     LOGICAL, OPTIONAL :: GotPort
     TYPE(Element_t), POINTER, OPTIONAL :: Element
     COMPLEX(KIND=dp), OPTIONAL :: B, L(:)
     REAL(KIND=dp), OPTIONAL :: Basis(:), dBasisdx(:,:), WBasis(:,:)

     LOGICAL :: Found
     TYPE(Solver_t), POINTER :: EigenSolver
     TYPE(Variable_t), POINTER :: EigenVar, PotVar
     REAL(KIND=dp), ALLOCATABLE :: Re_Eigenf(:), Im_Eigenf(:), ParentBasis(:)
     INTEGER :: EigenInd, PortDirection, PortTypeIndex, p, n, nd, m, i
     INTEGER, ALLOCATABLE :: DofInds(:)
     COMPLEX(KIND=dp) :: PortZ, PortBeta
     REAL(KIND=dp) :: Omega, PortLength, PortScale, PortCenter(3), mu0inv, muinv, eps0, rob0, &
         epsr, mur
     LOGICAL :: PortPassive !, ThinSheet, GoodConductor, Absorb, EigenSource, EigenWave
     TYPE(ValueHandle_t), SAVE :: EigenInd_h, PortTypeIndex_h, PortZ_h, PortLength_h, PortScale_h, &
         PortDirection_h, PortCenter_h, PortBeta_h, PortPassive_h, MuCoeff_h, EpsCoeff_h
     TYPE(Element_t), POINTER :: Parent
     COMPLEX(KIND=dp), PARAMETER :: im = (0._dp,1._dp)   
     CHARACTER(:), ALLOCATABLE :: str
     CHARACTER(*), PARAMETER :: Caller = 'VectorHelmholtzUtils'
     
     
     SAVE Omega, mu0inv, PortTypeIndex, EigenInd, Re_Eigenf, Im_Eigenf, EigenSolver, EigenVar, &
         PortBeta, PortZ, PortScale, PortDirection, PortLength, PortCenter, PortPassive, &
         DofInds, ParentBasis, m,  n, nd, eps0, rob0, Parent, PotVar
          
     
     SELECT CASE ( Phase ) 
     CASE( 1 )  ! Initialize Handles

       CALL ListInitElementKeyword( PortTypeIndex_h,'Boundary Condition','Port Type Index')
       CALL ListInitElementKeyword( PortZ_h,'Boundary Condition','Port Impedance',InitIm=.TRUE.)
       CALL ListInitElementKeyword( PortLength_h,'Boundary Condition','Port Length')
       CALL ListInitElementKeyword( PortScale_h,'Boundary Condition','Port Scale')
       CALL ListInitElementKeyword( PortDirection_h,'Boundary Condition','Port Direction',DefIValue=3)
       CALL ListInitElementKeyword( PortCenter_h,'Boundary Condition','Port Center',InitVec3D=.TRUE.)
       CALL ListInitElementKeyword( PortBeta_h,'Boundary Condition','Port Beta',InitIm=.TRUE.) 
       CALL ListInitElementKeyword( PortPassive_h,'Boundary Condition','Port Passive')
       CALL ListInitElementKeyword( EigenInd_h,'Boundary Condition','Eigenfunction Index')

       CALL ListInitElementKeyword( MuCoeff_h,'Material','Relative Reluctivity',InitIm=.TRUE.)      
       CALL ListInitElementKeyword( EpsCoeff_h,'Material','Relative Permittivity',InitIm=.TRUE.)      
       
       Omega = ListGetAngularFrequency( Solver % Values )

       Found = .FALSE.
       IF( ASSOCIATED( CurrentModel % Constants ) ) THEN
         mu0inv = ListGetConstReal( CurrentModel % Constants, 'Permeability of Vacuum', Found )
         IF(mu0inv/=0) mu0inv=1/mu0inv
       END IF
       IF(.NOT. Found ) mu0inv = 1.0_dp / ( PI * 4.0d-7 )
       
       Found = .FALSE.
       IF( ASSOCIATED( CurrentModel % Constants ) ) THEN
         eps0 = ListGetConstReal ( CurrentModel % Constants, 'Permittivity of Vacuum', Found )
       END IF
       IF(.NOT. Found ) eps0 = 8.854187817d-12

       rob0 = Omega * SQRT( eps0 / mu0inv )
       
       n = Solver % Mesh % MaxElementDOFs
       IF(.NOT. ALLOCATED(DofInds)) THEN
         ALLOCATE(DofInds(n),ParentBasis(n))
         DofInds = 0
         ParentBasis = 0.0_dp
       END IF
       
       ! If we have eigenfunction BC's then this has been set.
       m = ListGetInteger(Solver % Values, 'Eigensolver Index', Found)     
       IF(m > 0) THEN
         Eigensolver => CurrentModel % Solvers(m)              
         IF(.NOT. ALLOCATED(Re_eigenf) ) THEN
           ALLOCATE(Re_Eigenf(n), Im_Eigenf(n))
         END IF
         EigenVar => EigenSolver % Variable
       END IF       

       str = ListGetString( Solver % Values,'tem potential name',Found)
       IF(.NOT. Found) str = 'potential'
       PotVar => VariableGet( Solver % Mesh % Variables, str, ThisOnly=.TRUE.)

       
     CASE( 2 )  ! Visit new element

       
       PortTypeIndex = ListGetElementInteger(PortTypeIndex_h, Element, GotPort)
       IF(.NOT. GotPort) RETURN
       
       PortZ = ListGetElementComplex( PortZ_h, Element = Element )     
       PortScale = ListGetElementReal( PortScale_h, Element = Element )

       
       IF( PortTypeIndex == 1 ) THEN       ! rectangular
         PortDirection = ListGetElementInteger( PortDirection_h, Element )
         PortLength = ListGetElementReal( PortLength_h, Element = Element )
       ELSE IF( PortTypeIndex == 2 ) THEN  ! coaxial
         PortCenter = ListGetElementReal( PortCenter_h, Element = Element )
         CALL Fatal(Caller,'Unfinished port type: '//I2S(PortTypeIndex))        
       ELSE IF( PortTypeIndex == 3 ) THEN  ! eigenmode
         PortBeta = ListGetElementReal( PortBeta_h, Element = Element )
         EigenInd = MAX(1,ListGetElementInteger(EigenInd_h, Element, Found))

         n = Element % Type % NumberOfNodes
         m = mGetElementDOFs( DofInds, Element, USolver = EigenSolver )
         nd = m - n
         
         Re_eigenf(1:m) = REAL( EigenVar % EigenVectors(EigenInd,EigenVar % Perm(DofInds(1:m))) )
         Im_eigenf(1:m) = AIMAG( EigenVar % EigenVectors(EigenInd,EigenVar % Perm(DofInds(1:m))) )         
       ELSE IF( PortTypeIndex == 4 ) THEN
         Parent => Element % BoundaryInfo % Left
         IF(.NOT. ASSOCIATED(Parent)) THEN
           Parent => Element % BoundaryInfo % Right
         END IF
         IF(.NOT. ASSOCIATED(Parent)) THEN
           CALL Fatal(Caller,'Port model "potential" requires parent element!')
         END IF
         n = Element % Type % NumberOfNodes
       ELSE IF( PortTypeIndex == 5 ) THEN  
         PortBeta = ListGetElementReal( PortBeta_h, Element = Element, Found = Found )
         IF(.NOT. Found) CALL Fatal(Caller,'"Port Beta" not found for port type "beta"')
       ELSE         
         CALL Fatal(Caller,'Uncoded port type: '//I2S(PortTypeIndex))        
       END IF
       PortPassive = ListGetElementLogical( PortPassive_h, Element = Element )

       
     CASE( 3 )  ! Visit new integration point
       
       IF( PortTypeIndex == 1 ) THEN
         B = im * ( omega / mu0inv ) / (PortScale * PortZ ) 
         IF(PRESENT(L)) THEN
           L(ABS(PortDirection)) = SIGN(1,PortDirection) / ( PortLength * SQRT(PortScale) )
         END IF
       ELSE IF( PortTypeIndex == 3 ) THEN
         B = im * PortBeta
         IF( PRESENT(L)) THEN
           DO p=1,nd
             L(:) = L(:) + CMPLX(Re_Eigenf(n+p) * WBasis(p,:), Im_Eigenf(n+p) * WBasis(p,:), kind=dp) 
           END DO
         END IF
       ELSE IF( PortTypeIndex == 4 ) THEN
         epsr = ListGetElementComplex( EpsCoeff_h, Basis, Element, Found )
         IF(.NOT. Found) THEN
           CALL SetParentBasis( Element, n, Basis, Parent, Parent % TYPE % NumberOfNodes, ParentBasis)
           epsr = ListGetElementComplex( EpsCoeff_h, ParentBasis, Parent, Found )
         END IF
         IF(.NOT. Found ) epsr = 1.0_dp
         mur = ListGetElementComplex( MuCoeff_h, Basis, Element, Found )
         IF(.NOT. Found) THEN
           CALL SetParentBasis( Element, n, Basis, Parent, Parent % TYPE % NumberOfNodes, ParentBasis)
           mur = ListGetElementComplex( MuCoeff_h, Basis, Parent, Found )
         END IF
         IF(.NOT. Found ) mur = 1.0_dp
         B = -im * rob0 * SQRT( epsr / mur )          
         IF(PRESENT(L)) THEN
           DO i=1,3
             L(i) = SUM(dBasisdx(1:n,i) * PotVar % Values(PotVar % Perm(Element % NodeIndexes(1:n))))
           END DO
         END IF
       ELSE IF( PortTypeIndex == 5 ) THEN
         B = im * PortBeta
       ELSE
         IF(PRESENT(L)) L(1:3) = 0.0_dp
         B = 0.0_dp
         RETURN
       END IF

       IF( PRESENT(L)) THEN
         L = 2.0_dp * B * L 
         IF( PortPassive) L = 0.0_dp              
       END IF
         
     END SELECT


   END SUBROUTINE ElectricPortModel
  
   
 END MODULE VectorHelmholtzUtils
   
