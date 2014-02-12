!/*****************************************************************************
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
! * Module containing a solver for computing the strain rate Eij and tr(Eij)                               
! * 2D SDOFs = 5 (E11, E22, E33, E12, Eii)
! * 3D SDOFs = 7 (E11, E22, E33, E12, E23, E31, Eii)
! * Keywords : Flow Solver Name (AIFlow, Flow Solution, Porous, ...)
! *
! ******************************************************************************
! *
! *                     Author:       Juha Ruokolainen
! *
! *                    Address: CSC - IT Center for Science Ltd.
! *                            Keilaranta 14, P.O. Box 405
! *                              02101 Espoo, Finland
! *                              Tel. +358 0 457 2723
! *                            Telefax: +358 0 457 2302
! *                          EMail: Juha.Ruokolainen@csc.fi
! *
! *                       Date: 08 Jun 1997
! *
! *                Modified by:  Fabien / OG 
! *
! *       Date of modification: 13/10/05 from version 1.5
! *
! *****************************************************************************/
!------------------------------------------------------------------------------
   RECURSIVE SUBROUTINE ComputeStrainRate( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------

    USE DefUtils

    IMPLICIT NONE

!------------------------------------------------------------------------------
!******************************************************************************
!
!  Solve stress equations for one timestep
!
!  ARGUMENTS:
!
!  TYPE(Model_t) :: Model,  
!     INPUT: All model information (mesh,materials,BCs,etc...)
!
!  TYPE(Solver_t) :: Solver
!     INPUT: Linear equation solver options
!
!  REAL(KIND=dp) :: dt,
!     INPUT: Timestep size for time dependent simulations (NOTE: Not used
!            currently)
!
!******************************************************************************

     TYPE(Model_t)  :: Model
     TYPE(Solver_t), TARGET :: Solver

     LOGICAL ::  TransientSimulation
     REAL(KIND=dp) :: dt
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     TYPE(Solver_t), POINTER :: PSolver

     TYPE(Matrix_t),POINTER :: StiffMatrix

     INTEGER :: i, j, k, l, n, t, iter, NDeg, STDOFs, LocalNodes, istat
     INTEGER :: dim

     TYPE(ValueList_t),POINTER :: Material, BC
     TYPE(Nodes_t) :: ElementNodes
     TYPE(Element_t),POINTER :: CurrentElement

     REAL(KIND=dp) :: RelativeChange, UNorm, PrevUNorm, s
                      

     REAL(KIND=dp), ALLOCATABLE :: Basis(:),ddBasisddx(:,:,:)
     REAL(KIND=dp), ALLOCATABLE :: dBasisdx(:,:)
     REAL(KIND=dp) :: u, v, w, detJ
     
     LOGICAL :: stat, CSymmetry 

     INTEGER :: NewtonIter, NonlinearIter

     TYPE(Variable_t), POINTER :: EiiSol, FlowVariable 

     REAL(KIND=dp), POINTER :: Eii(:), FlowValues(:), &
           ForceVector(:),  NodalEii(:), EiiComp(:)

     INTEGER, POINTER :: EiiPerm(:), NodeIndexes(:), &
                         FlowPerm(:)

     LOGICAL :: GotIt, AllocationsDone = .FALSE.

     REAL(KIND=dp), ALLOCATABLE:: LocalMassMatrix(:,:), &
       LocalStiffMatrix(:,:), LocalForce(:), &
       LocalVelo(:,:)
            
     CHARACTER(LEN=MAX_NAME_LEN) :: FlowSolverName

     REAL(KIND=dp) :: at, at0, CPUTime, RealTime


!------------------------------------------------------------------------------
     SAVE Basis, dBasisdx, ddBasisddx
     SAVE LocalMassMatrix, LocalStiffMatrix, LocalForce, &
          ElementNodes, AllocationsDone  
     SAVE LocalVelo,  dim

!------------------------------------------------------------------------------
!  Read the name of the Flow Solver (NS, AIFlow, Porous, ...)
!------------------------------------------------------------------------------
      
     FlowSolverName = GetString( Solver % Values, 'Flow Solver Name', GotIt )    
     IF (.NOT.Gotit) FlowSolverName = 'aiflow'
     FlowVariable => VariableGet( Solver % Mesh % Variables, FlowSolverName )
     IF ( ASSOCIATED( FlowVariable ) ) THEN
       FlowPerm    => FlowVariable % Perm
       FlowValues  => FlowVariable % Values
     ELSE
       CALL Info('ComputeStrainRate', &
                      & 'No variable for velocity associated.', Level=4)
     END IF
!              
!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
      IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN

      EiiSol => Solver % Variable
      EiiPerm => EiiSol % Perm
      STDOFs =  EiiSol % DOFs
      Eii => EiiSol % Values

      LocalNodes = COUNT( EiiPerm > 0 )
      IF ( LocalNodes <= 0 ) RETURN


      StiffMatrix => Solver % Matrix
      ForceVector => StiffMatrix % RHS
      UNorm = Solver % Variable % Norm

!------------------------------------------------------------------------------
!     Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
      IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed) THEN
        N = Model % MaxElementNodes
        dim = CoordinateSystemDimension()

       IF ( AllocationsDone ) THEN
         DEALLOCATE( ElementNodes % x,     &
                     ElementNodes % y,     &
                     ElementNodes % z,     &
                     LocalVelo,            &                      
                     Basis, ddBasisddx,    &
                     dBasisdx,             &
                     LocalMassMatrix,      &
                     LocalStiffMatrix,     &
                     LocalForce )
       END IF

       ALLOCATE( ElementNodes % x( N ), &
                 ElementNodes % y( N ), &
                 ElementNodes % z( N ), &
                 LocalVelo( 3,N ), &                                    
                 Basis( 2*N ),ddBasisddx(1,1,1), dBasisdx( 2*N,3 ), &
                 LocalMassMatrix( 2*STDOFs*N,2*STDOFs*N ),  &
                 LocalStiffMatrix( 2*STDOFs*N,2*STDOFs*N ),  &
                 LocalForce( 2*STDOFs*N ),  STAT=istat )

       IF ( istat /= 0 ) THEN
          CALL Fatal( 'ComputeStrainRate', 'Memory allocation error.' )
       END IF
!------------------------------------------------------------------------------
       AllocationsDone = .TRUE.
      END IF


!------------------------------------------------------------------------------
      NonlinearIter = 1
      DO iter=1,NonlinearIter

       at  = CPUTime()
       at0 = RealTime()

       CALL Info( 'ComputeStrainRate', ' ', Level=4 )
       CALL Info( 'ComputeStrainRate', ' ', Level=4 )
       CALL Info( 'ComputeStrainRate', ' ', Level=4 )
       CALL Info( 'ComputeStrainRate', ' ', Level=4 )
       CALL Info( 'ComputeStrainRate', 'Starting assembly...',Level=4 )

!------------------------------------------------------------------------------
       CALL DefaultInitialize()
!------------------------------------------------------------------------------
       DO t=1,Solver % NumberOFActiveElements

         IF ( RealTime() - at0 > 1.0 ) THEN
           WRITE(Message,'(a,i3,a)' ) '   Assembly: ',  &
             INT(100.0 - 100.0 * (Solver % NumberOfActiveElements-t) / &
             (1.0*Solver % NumberOfActiveElements)), ' % done'
           CALL Info( 'ComputeStrainRate', Message, Level=5 )
           at0 = RealTime()
         END IF

         CurrentElement => GetActiveElement(t)
         n = GetElementNOFNodes()
         NodeIndexes => CurrentElement % NodeIndexes

         ElementNodes % x(1:n) = Model % Nodes % x(NodeIndexes(1:n))
         ElementNodes % y(1:n) = Model % Nodes % y(NodeIndexes(1:n))
         ElementNodes % z(1:n) = Model % Nodes % z(NodeIndexes(1:n))

         Material => GetMaterial()


!------------------------------------------------------------------------------
!        Get element local stiffness & mass matrices
!------------------------------------------------------------------------------

          LocalVelo = 0.0d0
          DO i=1, dim
             LocalVelo(i,1:n) = FlowValues((dim+1)*(FlowPerm(NodeIndexes(1:n))-1) + i)
          END DO

          CALL LocalMatrix( LocalMassMatrix, LocalStiffMatrix, &
              LocalForce, LocalVelo, CurrentElement, n, ElementNodes )
              

!------------------------------------------------------------------------------
!        Update global matrices from local matrices 
!------------------------------------------------------------------------------
         CALL DefaultUpdateEquations( LocalStiffMatrix, LocalForce )
 
      END DO

      CALL Info( 'ComputeStrainRate', 'Assembly done', Level=4 )


      CALL DefaultFinishAssembly()

!------------------------------------------------------------------------------
!     Dirichlet boundary conditions
!------------------------------------------------------------------------------
      CALL DefaultDirichletBCs()

!------------------------------------------------------------------------------

      CALL Info( 'ComputeStrainRate', 'Set boundaries done', Level=4 )

!------------------------------------------------------------------------------
!     Solve the system and check for convergence
!------------------------------------------------------------------------------
      PrevUNorm = UNorm

      UNorm = DefaultSolve()

      unorm = 0
      n = Solver % Variable % DOFs
      DO i=1,n-1
        unorm = unorm + SUM( solver % variable % values(i::n)**2 )
      END DO
      unorm = SQRT( unorm )
      solver % variable % norm = unorm

      IF ( PrevUNorm + UNorm /= 0.0d0 ) THEN
         RelativeChange = 2.0d0 * ABS( PrevUNorm - UNorm) / ( PrevUnorm + UNorm)
      ELSE
         RelativeChange = 0.0d0
      END IF

      WRITE( Message, * ) 'Result Norm   : ',UNorm, PrevUNorm
      CALL Info( 'ComputeStrainRate', Message, Level=4 )
      WRITE( Message, * ) 'Relative Change : ',RelativeChange
      CALL Info( 'ComputeStrainRate', Message, Level=4 )


!------------------------------------------------------------------------------
    END DO ! of nonlinear iter
!------------------------------------------------------------------------------

      
CONTAINS


!------------------------------------------------------------------------------
      SUBROUTINE LocalMatrix( MassMatrix, StiffMatrix, ForceVector, &
              NodalVelo, Element, n, Nodes )
              
!------------------------------------------------------------------------------

     REAL(KIND=dp) :: StiffMatrix(:,:), MassMatrix(:,:)
     REAL(KIND=dp) ::  NodalVelo(:,:), ForceVector(:)
     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t) :: Element
     INTEGER :: n
!------------------------------------------------------------------------------
!
     REAL(KIND=dp) :: Basis(2*n), ddBasisddx(1,1,1)
     REAL(KIND=dp) :: dBasisdx(2*n,3), detJ

     REAL(KIND=dp) :: LGrad(3,3), SR(3,3),  Eij(7)

     INTEGER :: i, j, k, p, q, t, dim, cc

     REAL(KIND=dp) :: s, u, v, w, Radius
  
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     INTEGER :: N_Integ

     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ, V_Integ, W_Integ, S_Integ

     LOGICAL :: stat, CSymmetry

!------------------------------------------------------------------------------
      dim = CoordinateSystemDimension()
      cc = 2*dim + 1

      ForceVector = 0.0D0
      StiffMatrix = 0.0D0
      MassMatrix  = 0.0D0

      IntegStuff = GaussPoints( Element )

      U_Integ => IntegStuff % u
      V_Integ => IntegStuff % v
      W_Integ => IntegStuff % w
      S_Integ => IntegStuff % s
      N_Integ =  IntegStuff % n
!
!   Now we start integrating
!
      DO t=1,N_Integ

      u = U_Integ(t)
      v = V_Integ(t)
      w = W_Integ(t)

!------------------------------------------------------------------------------
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo(Element,Nodes,u,v,w,detJ, &
                 Basis,dBasisdx,ddBasisddx,.FALSE.,.FALSE.)


       s = detJ * S_Integ(t)


       Radius = SUM( Nodes % x(1:n) * Basis(1:n) )
       CSymmetry = CurrentCoordinateSystem() == AxisSymmetric
       IF ( CSymmetry ) s = s * Radius
!
! Strain-Rate and Eii = tr(Eij)
!
        SR = 0.0
        Eij = 0.0
        LGrad = MATMUL( NodalVelo(:,1:n), dBasisdx(1:n,:) )
        SR = 0.5 * ( LGrad + TRANSPOSE(LGrad) )
        IF ( CSymmetry ) THEN
          SR(1,3) = 0.0
          SR(2,3) = 0.0
          SR(3,1) = 0.0
          SR(3,2) = 0.0
          SR(3,3) = 0.0
          IF ( Radius > 10*AEPS ) THEN
            SR(3,3) = SUM( Nodalvelo(1,1:n) * Basis(1:n) ) /Radius
                 
          END IF
        END IF
        Eij(1) = SR(1,1)        
        Eij(2) = SR(2,2)        
        Eij(3) = SR(3,3)        
        Eij(4) = SR(1,2)        
        IF (dim > 2) THEN
          Eij(5) = SR(2,3)
          Eij(6) = SR(3,1)
        END IF
        Eij(cc) = SR(1,1) + SR(2,2) + SR(3,3)

        DO i=1,cc 
          DO p=1,n         
            DO q=1,n        
            StiffMatrix(cc*(p-1)+i,cc*(q-1)+i ) =  &
               StiffMatrix(cc*(p-1)+i,cc*(q-1)+i ) + s*Basis(q)*Basis(p)
            END DO

             ForceVector(cc*(p-1)+i) =  &
                     ForceVector(cc*(p-1)+i) + s*Eij(i)*Basis(p) 
          END DO
        END DO
      END DO 

!------------------------------------------------------------------------------
      END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
      END SUBROUTINE ComputeStrainRate
!------------------------------------------------------------------------------
