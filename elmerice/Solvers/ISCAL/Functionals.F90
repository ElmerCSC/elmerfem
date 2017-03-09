!
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
!
!------------------------------------------------------------------------------
!>  Module computing functional values used to estimate errors
!------------------------------------------------------------------------------


MODULE Functionals

CONTAINS

  SUBROUTINE FluxAcrossPoint ( Model,Solver,dt,TransientSimulation,a  )
    !******************************************************************************
    !
    !  Write a function for computing the horizontal flux over a vertical line,
    !  in the form of a vector, a. a = integral of vx over the line. The 
    !  trapezoidal method is used to compute the integral. 
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
    !     INPUT: Timestep size for time dependent simulations
    !
    !  REAL(KIND=dp) :: a,
    !     OUTPUT: value of functional 
    !
    !
    !******************************************************************************
    !------------------------------------------------------------------------------


    USE DefUtils
    IMPLICIT NONE

    TYPE(Model_t) :: Model
    TYPE(Solver_t), TARGET :: Solver

    REAL(KIND=dp) :: dt
    LOGICAL :: TransientSimulation

    TYPE(Variable_t), POINTER :: FlowSol 
    INTEGER, POINTER :: FlowPerm(:) 
    REAL(KIND=dp), POINTER :: FlowSolution(:) 

    REAL(KIND=dp) :: xcoord
    REAL(KIND=dp), POINTER :: a(:)
    INTEGER :: i, mini
    REAL(KIND=dp) :: dist, olddist

    !-----------Variables needed for integration --------------------------------
    TYPE(Solver_t), POINTER :: PSolver
    INTEGER :: j,k,l,Dofs,dof,nsize,TopNodes,BotNodes
    INTEGER, POINTER :: TopPointer(:),BotPointer(:),UpPointer(:),DownPointer(:)
    LOGICAL :: Initialized = .FALSE.,GotVar
    REAL(KIND=dp) :: dx,Level,q
    REAL(KIND=dp), POINTER :: Coord(:)
    TYPE(Variable_t), POINTER :: Var, OldVar

    SAVE :: BotPointer, TopPointer, UpPointer, DownPointer, Coord, TopNodes, &
         BotNodes

    !------------------------------------------------------------------------------
    !    Get the solution  vx,vy,p
    !------------------------------------------------------------------------------
    FlowSol => Solver % Variable
    FlowPerm       => FlowSol % Perm
    FlowSolution   => FlowSol % Values


    !------------------------------------------------------------------------------
    !   Initialize the pointers to top and bottom nodes, needed for the integration 
    !------------------------------------------------------------------------------
    IF( .NOT. Initialized ) THEN

       ! Choose active direction coordinate and set corresponding unit vector
       !---------------------------------------------------------------------
       PSolver => Solver
       CALL DetectExtrudedStructure( Solver % Mesh, PSolver, Var, &
            TopNodePointer = TopPointer, BotNodePointer = BotPointer, &
            UpNodePointer = UpPointer, DownNodePointer = DownPointer )

       Coord => Var % Values
       nsize = SIZE( Coord )
       Initialized = .TRUE.
    END IF

    !------------------------------------------------------------------------------
    !    Figure out what line you wanna integrate
    !------------------------------------------------------------------------------
    !get the coordinate  
    xcoord = GetConstReal( Solver % Values, 'Point x-coord', GotVar)    
    IF (.NOT. GotVar) THEN
       CALL FATAL( 'Error Estimation: ','Point x-coord not set')
    END IF


    olddist = 10e10 !big number
    !figure out which line is closest
    DO i=1,Model % Mesh % NumberOfNodes
       dist = ABS(Model % Nodes % x(i) - xcoord)
       IF (dist < olddist) THEN
          mini=i
          olddist=dist
       END IF
    END DO


    CALL Info( 'Error Estimation: ',Message, Level=4 )
    !------------------------------------------------------------------------------
    !    Compute a
    !------------------------------------------------------------------------------

    a=0.0      

    !Get the bottom pointer for the node you found above
    i = BotPointer(mini)

    dx = (Coord(UpPointer(i)) - Coord(i))

    a(3*(FlowPerm(i)-1)+1)= 0.5*dx

    DO WHILE (i /= TopPointer(i))
       i = UpPointer(i)
       dx = (Coord(UpPointer(i)) - Coord(DownPointer(i)))
       a(3*(FlowPerm(i)-1)+1)=0.5*dx
    END DO

    dx = (Coord(i) - Coord(DownPointer(i)))
    a(3*(FlowPerm(i)-1)+1)= 0.5*dx

  END SUBROUTINE FluxAcrossPoint

  SUBROUTINE FluxAcrossLine ( Model,Solver,dt,TransientSimulation, &
       a  )
    !******************************************************************************
    !
    !  Write a function for computing the horizontal flux over a vertical line,
    !  in the form of a vector, a. a = integral of vx over the line. The 
    !  trapezoidal method is used to compute the integral. 
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
    !     INPUT: Timestep size for time dependent simulations
    !
    !  REAL(KIND=dp) :: a,
    !     OUTPUT: value of functional 
    !
    !
    !******************************************************************************
    !------------------------------------------------------------------------------


    USE DefUtils
    IMPLICIT NONE

    TYPE(Model_t) :: Model
    TYPE(Solver_t), TARGET :: Solver

    REAL(KIND=dp) :: dt,ds,R, nx, ny
    LOGICAL :: TransientSimulation

    TYPE(Variable_t), POINTER :: FlowSol 
    INTEGER, POINTER :: FlowPerm(:) 
    REAL(KIND=dp), POINTER :: FlowSolution(:) 

    CHARACTER(LEN=MAX_NAME_LEN) :: filename
    INTEGER :: i,linelength,istat
    INTEGER, ALLOCATABLE :: linenode(:)
    REAL(KIND=dp), POINTER :: a(:)

    !-----------Variables needed for integration --------------------------------
    TYPE(Solver_t), POINTER :: PSolver
    INTEGER :: j,k,l,Dofs,dof,nsize,TopNodes,BotNodes
    INTEGER, POINTER :: TopPointer(:),BotPointer(:),UpPointer(:),DownPointer(:)
    LOGICAL :: Initialized = .FALSE.,GotVar
    REAL(KIND=dp) :: dx,Level,q
    REAL(KIND=dp), POINTER :: Coord(:)
    TYPE(Variable_t), POINTER :: Var, OldVar

    SAVE :: BotPointer, TopPointer, UpPointer, DownPointer, Coord, TopNodes, &
         BotNodes

    !------------------------------------------------------------------------------
    !    Get the solution  vx,vy,p
    !------------------------------------------------------------------------------
    FlowSol => Solver % Variable
    FlowPerm       => FlowSol % Perm
    FlowSolution   => FlowSol % Values


    !------------------------------------------------------------------------------
    !   Initialize the pointers to top and bottom nodes, needed for the integration 
    !------------------------------------------------------------------------------
    IF( .NOT. Initialized ) THEN

       ! Choose active direction coordinate and set corresponding unit vector
       !---------------------------------------------------------------------
       PSolver => Solver
       CALL DetectExtrudedStructure( Solver % Mesh, PSolver, Var, &
            TopNodePointer = TopPointer, BotNodePointer = BotPointer, &
            UpNodePointer = UpPointer, DownNodePointer = DownPointer )

       Coord => Var % Values
       nsize = SIZE( Coord )
       Initialized = .TRUE.
    END IF

    !------------------------------------------------------------------------------
    !    Figure out what line you wanna integrate over
    !------------------------------------------------------------------------------

    !get the description of the line
    filename = GetString( Solver % Values, 'Line Description', GotVar)    
    IF (.NOT. GotVar) THEN
       CALL FATAL( 'Error Estimation: ','File describing line not found')
    END IF

    OPEN(unit = 102, file = FileName, status = 'old', action = 'read')
    read(102,*) linelength

    ALLOCATE( linenode(linelength), STAT=istat )

    IF ( istat /= 0 ) THEN
       CALL Fatal( 'Error Estimation','Memory allocation error, Aborting.' )
    END IF

    DO i=1,linelength
       READ(102,*) linenode(i)
    END DO
    
    READ(102,*) ds !cirkelbågelängd
    READ(102,*) R !radie


    CLOSE(102)
    

    !------------------------------------------------------------------------------
    !    Compute a
    !------------------------------------------------------------------------------

    a=0.0      

    DO k=1,linelength !looping over the nodes in the line, computing an integral from bottom to top
    
    nx =  Solver % Mesh % Nodes % x(linenode(k))/R         
    ny =  Solver % Mesh % Nodes % y(linenode(k))/R 

    i = BotPointer(linenode(k)) !bottom node
    
!-----------
    dx = (Coord(UpPointer(i)) - Coord(i)) 

    a(4*(FlowPerm(i)-1)+1)= 0.5*nx*dx*ds
    a(4*(FlowPerm(i)-1)+2)= 0.5*ny*dx*ds

    DO WHILE (i /= TopPointer(i))
       i = UpPointer(i)
       dx = 0.5*(Coord(UpPointer(i)) - Coord(DownPointer(i)))
       a(4*(FlowPerm(i)-1)+1)=nx*dx*ds
       a(4*(FlowPerm(i)-1)+2)=ny*dx*ds
    END DO

    dx = (Coord(i) - Coord(DownPointer(i)))
    a(4*(FlowPerm(i)-1)+1)= 0.5*nx*dx*ds    
    a(4*(FlowPerm(i)-1)+2)= 0.5*ny*dx*ds

!----------
END DO

    DEALLOCATE(linenode)

  END SUBROUTINE FluxAcrossLine


END MODULE Functionals
