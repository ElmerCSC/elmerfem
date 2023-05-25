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
! *  Module for ODY solver with global dofs.
! *
! *  Authors: Peter Råback
! *  Email:   Peter.Raback@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 4.5.2015
! *
! *****************************************************************************/

!> \ingroup Solvers
!> \{
 

!------------------------------------------------------------------------------
!> Initialization for the primary solver
!------------------------------------------------------------------------------
SUBROUTINE OdeSolver_init( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver       !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model         !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt            !< Timestep size for time dependent simulations
  LOGICAL :: TransientSimulation !< Steady state or transient simulation
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: Params

  Params => GetSolverParams()

  IF(.NOT. ListCheckPresent( Params,'Time Derivative Order') ) THEN
    CALL ListAddInteger( Params, 'Time derivative order', 2 )
  END IF
  CALL ListAddLogical( Params, 'Use Global Mass Matrix', .TRUE. )
  CALL ListAddLogical( Params, 'Ode Matrix', .TRUE. )
  CALL ListAddLogical( Params, 'Variable Global', .TRUE. )

!------------------------------------------------------------------------------
END SUBROUTINE OdeSolver_init
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> This is a simple solver for ordinary differential equations (ODE) 
!> utilizing many of the same subroutines as the default PDE solution. 
!> The idea is to easily allow the use of same time integration schemes
!> etc. There is plenty of room for improvement and polishing.
!------------------------------------------------------------------------------
SUBROUTINE OdeSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver       !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model         !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt            !< Timestep size for time dependent simulations
  LOGICAL :: TransientSimulation !< Steady state or transient simulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  LOGICAL :: Found
  REAL(KIND=dp) :: Norm  
  INTEGER :: iter, MaxIter,TimeOrder
  CHARACTER(LEN=MAX_NAME_LEN) :: str
  TYPE(ValueList_t), POINTER :: Params
!------------------------------------------------------------------------------

  CALL Info('OdeSolver','Solving ordinary differential equation',Level=4)


  Params => GetSolverParams()
 
  MaxIter = GetInteger( Params,'Nonlinear System Max Iterations',Found )
  IF(.NOT. Found) MaxIter = 1

  TimeOrder = GetInteger( Params,'Time Derivative Order' )

  DO iter = 1, MaxIter 
    IF( MaxIter > 1 ) THEN
      CALL Info('OdeSolver','Nonlinear iteration: '//I2S(iter),Level=5)
    END IF

    CALL DefaultInitialize()
    
    CALL OdeMatrix()

    CALL DefaultFinishBulkAssembly()
    CALL DefaultFinishBoundaryAssembly()
    
    ! This sets the time-integration and is therefore imperative
    CALL DefaultFinishAssembly()
    
    ! This is probably not needed...
    !CALL DefaultDirichletBCs()
    
    Norm = DefaultSolve()

    IF( Solver % Variable % NonlinConverged == 1 ) EXIT
  END DO

  CALL Info('OdeSolver','All done',Level=5)
 

CONTAINS


  ! Creates the local matrix equation for the ODY before time integration. 
  ! This optionally uses > Active Components < where the component entry
  ! could be a natural place for parameters of the ODE system.
  !-----------------------------------------------------------------------

  SUBROUTINE OdeMatrix()
    
    TYPE(Matrix_t), POINTER :: A
    INTEGER :: i,j,Dofs
    REAL(KIND=dp), POINTER CONTIG :: SaveValues(:)
    REAL(KIND=dp) :: val
    TYPE(ValueList_t), POINTER :: OdeList
    INTEGER, POINTER :: ActiveComponents(:)


    ActiveComponents => ListGetIntegerArray( Params, &        
        'Active Components', Found )
    IF( Found ) THEN
      IF( SIZE( ActiveComponents ) > 1 ) THEN
        CALL Fatal('OdeSolver','Currently implemented only for one component!')
      END IF
      i = ActiveComponents(1) 
      IF( i > CurrentModel % NumberOfComponents ) THEN
        CALL Fatal('OdeSolver','Active Component index out of range')
      END IF
      OdeList => CurrentModel % Components(i) % Values 
      CALL Info('OdeSolver','Using active component: '//I2S(i),Level=10)
    ELSE
      OdeList => Params
      CALL Info('OdeSolver','Using solver section',Level=10)
    END IF

    A => Solver % Matrix
    Dofs = Solver % Variable % Dofs
    
    IF ( .NOT. ASSOCIATED( A % MassValues ) ) THEN
      CALL Info('OdeSolver','Allocating mass matrix',Level=10)
      ALLOCATE( A % MassValues(SIZE(A % Values)) )
      A % MassValues = 0.0_dp
    END IF
    
    IF( TimeOrder >= 2 ) THEN
      IF ( .NOT. ASSOCIATED( A % DampValues ) ) THEN
        CALL Info('OdeSolver','Allocating damping matrix',Level=10)
        ALLOCATE( A % DampValues(SIZE(A % Values)) )
        A % DampValues = 0.0_dp
      END IF
    END IF

    ! Set stiffness matrix values
    DO i=1,Dofs
      DO j=1,Dofs
        str = 'Stiffness Matrix '//I2S(i)//I2S(j)
        val = ListGetCReal( OdeList, str, Found ) 
        CALL CRS_SetMatrixElement( A,i,j,val )
      END DO
    END DO

    ! Set mass matrix values
    SaveValues => A % Values
    A % Values => A % MassValues
    DO i=1,Dofs
      DO j=1,Dofs

        ! In Elmer the convention is to call the highest time-derivatie always mass
        ! Even for 1st order PDEs. However, in this solver we like to call the coefficients
        ! mass, damping and spring.
        IF( TimeOrder >= 2 ) THEN
          str = 'Mass Matrix '//I2S(i)//I2S(j)
        ELSE
          str = 'Damping Matrix '//I2S(i)//I2S(j)
        END IF
        val = ListGetCReal( OdeList, str, Found ) 
        CALL CRS_SetMatrixElement( A,i,j,val )
      END DO
    END DO
    A % Values => SaveValues

    IF( TimeOrder >= 2 ) THEN
      ! Set damp matrix values
      SaveValues => A % Values
      A % Values => A % DampValues
      DO i=1,Dofs
        DO j=1,Dofs
          str = 'Damping Matrix '//I2S(i)//I2S(j)
          val = ListGetCReal( OdeList, str, Found ) 
          CALL CRS_SetMatrixElement( A,i,j,val )
        END DO
      END DO
      A % Values => SaveValues
    END IF

    DO i=1,Dofs
      str = 'Force '//I2S(i)
      val = ListGetCReal( OdeList, str, Found ) 
      A % Rhs(i) = val
    END DO

    IF(.FALSE.) THEN
      PRINT *,'Dofs:',Dofs
      PRINT *,'A % Values',A % values
      PRINT *,'A % DampValues',A % Dampvalues
      IF( TimeOrder >= 2 ) THEN
        PRINT *,'A % MassValues',A % Massvalues
      END IF
      PRINT *,'A % Rhs',A % Rhs
    END IF

  END SUBROUTINE OdeMatrix

!------------------------------------------------------------------------------
END SUBROUTINE OdeSolver
!------------------------------------------------------------------------------
