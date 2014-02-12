
  FUNCTION LeftWall( Model, n, t ) RESULT(f)

    USE Types
    USE SolverUtils
    USE Lists
    USE Integration
    USE ElementDescription
    USE ElementUtils
    USE DefUtils

    IMPLICIT NONE
    
    TYPE(Model_t) :: Model
    INTEGER :: n
    REAL(KIND=dp) :: t,f
!-----------------------------------------------------
    LOGICAL :: Visited = .FALSE.,Found, FixedArea, IsLeftWall
    INTEGER :: i,j,Nopt,Power1,Nleft,Nboth,nstart,nfin
    REAL(KIND=dp) :: x,y,fi
    TYPE(Variable_t), POINTER :: OptSol
    REAL(KIND=dp), POINTER :: Opt(:)

    SAVE Visited, Opt, Nopt, Power1, Nleft, Nboth, nstart, nfin, FixedArea
    
    IsLeftWall = .TRUE.

    IF(.NOT. Visited) THEN
      FixedArea = ListGetLogical( Model % Simulation,'Fixed Area',Found)
      OptSol => VariableGet( Model % Variables, 'opt' )
      IF( ASSOCIATED( OptSol ) ) THEN
        Opt => OptSol % Values
        Nopt = SIZE( Opt )
      ELSE
        ALLOCATE( Opt(10) )
        
        ! Shorter names may be needed for some uses
        IF(.FALSE.) THEN
           Opt(1) = ListGetConstReal( Model % Simulation,'Left Wall B1',Found)
           Opt(3) = ListGetConstReal( Model % Simulation,'Left Wall B2',Found)
           Opt(5) = ListGetConstReal( Model % Simulation,'Left Wall B3',Found)
           Opt(7) = ListGetConstReal( Model % Simulation,'Left Wall B4',Found)
           Opt(9) = ListGetConstReal( Model % Simulation,'Left Wall B5',Found)
           
           Opt(2) = ListGetConstReal( Model % Simulation,'Right Wall B1',Found)
           Opt(4) = ListGetConstReal( Model % Simulation,'Right Wall B2',Found)
           Opt(6) = ListGetConstReal( Model % Simulation,'Right Wall B3',Found)
           Opt(8) = ListGetConstReal( Model % Simulation,'Right Wall B4',Found)
           Opt(10) = ListGetConstReal( Model % Simulation,'Right Wall B5',Found)
        ELSE
           Opt(1) = ListGetConstReal( Model % Simulation,'LB1',Found)
           Opt(3) = ListGetConstReal( Model % Simulation,'LB2',Found)
           Opt(5) = ListGetConstReal( Model % Simulation,'LB3',Found)
           Opt(7) = ListGetConstReal( Model % Simulation,'LB4',Found)
           Opt(9) = ListGetConstReal( Model % Simulation,'LB5',Found)
           
           Opt(2) = ListGetConstReal( Model % Simulation,'RB1',Found)
           Opt(4) = ListGetConstReal( Model % Simulation,'RB2',Found)
           Opt(6) = ListGetConstReal( Model % Simulation,'RB3',Found)
           Opt(8) = ListGetConstReal( Model % Simulation,'RB4',Found)
           Opt(10) = ListGetConstReal( Model % Simulation,'RB5',Found)  
        END IF



        IF( FixedArea ) THEN
          Nopt = 9
          Opt(2:9) = Opt(3:10)
        ELSE
          Nopt = 10
        END IF
      END IF

      Visited = .TRUE.

      Power1 = 2
      IF( FixedArea ) THEN
        nstart = 2
      ELSE
        nstart = 3
      END IF

      Nfin = Nopt
    END IF

    x = Model % Nodes % x(n)
    y = Model % Nodes % y(n)    
    y = 1 - y
        
    f = 0.0_dp    
    IF( IsLeftWall .OR. FixedArea ) THEN
      f = f + Opt(1) * y 
    ELSE
      f = f + Opt(2) * y
    END IF
    
    i = Power1
    IF( IsLeftWall ) THEN
      DO j=nstart,nfin,2
        fi = (i+1)*y**i - i*y**(i-1)
        f = f + Opt(j) * fi
        i = i + 1
      END DO
    ELSE
      DO j=nstart,nfin,2
        fi = (i+1)*y**i - i*y**(i-1)
        f = f + Opt(j+1) * fi
        i = i + 1
      END DO
    END IF


!    PRINT *,'y',y,'dx',f,'Opt',Power1,Opt

  END FUNCTION LeftWall



  FUNCTION RightWall( Model, n, t ) RESULT(f)

    USE Types
    USE SolverUtils
    USE Lists
    USE Integration
    USE ElementDescription
    USE ElementUtils
    USE DefUtils

    IMPLICIT NONE
    
    TYPE(Model_t) :: Model
    INTEGER :: n
    REAL(KIND=dp) :: t,f
!-----------------------------------------------------
    LOGICAL :: Visited = .FALSE.,Found, FixedArea, IsLeftWall
    INTEGER :: i,j,Nopt,Power1,Nleft,Nboth,nstart,nfin
    REAL(KIND=dp) :: x,y,fi
    TYPE(Variable_t), POINTER :: OptSol
    REAL(KIND=dp), POINTER :: Opt(:)

    SAVE Visited, Opt, Nopt, Power1, Nleft, Nboth, nstart, nfin, FixedArea
    
    IsLeftWall = .FALSE.

    IF(.NOT. Visited) THEN
      FixedArea = ListGetLogical( Model % Simulation,'Fixed Area',Found)
      OptSol => VariableGet( Model % Variables, 'opt' )
  
      IF( ASSOCIATED( OptSol ) ) THEN
        Opt => OptSol % Values
        Nopt = SIZE( Opt )
      ELSE
        ALLOCATE( Opt(10) )
        
        IF(.FALSE.) THEN
           Opt(1) = ListGetConstReal( Model % Simulation,'Left Wall B1',Found)
           Opt(3) = ListGetConstReal( Model % Simulation,'Left Wall B2',Found)
           Opt(5) = ListGetConstReal( Model % Simulation,'Left Wall B3',Found)
           Opt(7) = ListGetConstReal( Model % Simulation,'Left Wall B4',Found)
           Opt(9) = ListGetConstReal( Model % Simulation,'Left Wall B5',Found)
           
           Opt(2) = ListGetConstReal( Model % Simulation,'Right Wall B1',Found)
           Opt(4) = ListGetConstReal( Model % Simulation,'Right Wall B2',Found)
           Opt(6) = ListGetConstReal( Model % Simulation,'Right Wall B3',Found)
           Opt(8) = ListGetConstReal( Model % Simulation,'Right Wall B4',Found)
           Opt(10) = ListGetConstReal( Model % Simulation,'Right Wall B5',Found)
        ELSE
           Opt(1) = ListGetConstReal( Model % Simulation,'LB1',Found)
           Opt(3) = ListGetConstReal( Model % Simulation,'LB2',Found)
           Opt(5) = ListGetConstReal( Model % Simulation,'LB3',Found)
           Opt(7) = ListGetConstReal( Model % Simulation,'LB4',Found)
           Opt(9) = ListGetConstReal( Model % Simulation,'LB5',Found)
           
           Opt(2) = ListGetConstReal( Model % Simulation,'RB1',Found)
           Opt(4) = ListGetConstReal( Model % Simulation,'RB2',Found)
           Opt(6) = ListGetConstReal( Model % Simulation,'RB3',Found)
           Opt(8) = ListGetConstReal( Model % Simulation,'RB4',Found)
           Opt(10) = ListGetConstReal( Model % Simulation,'RB5',Found)  
        END IF



        IF( FixedArea ) THEN
          Nopt = 9
          Opt(2:9) = Opt(3:10)
        ELSE
          Nopt = 10
        END IF
      END IF

      Visited = .TRUE.

      Power1 = 2
      IF( FixedArea ) THEN
        nstart = 2
      ELSE
        nstart = 3
      END IF

      Nfin = Nopt
    END IF

    x = Model % Nodes % x(n)
    y = Model % Nodes % y(n)    
    y = 1 - y
        
    f = 0.0_dp    
    IF( IsLeftWall .OR. FixedArea ) THEN
      f = f + Opt(1) * y 
    ELSE
      f = f + Opt(2) * y
    END IF
    
    i = Power1
    IF( IsLeftWall ) THEN
      DO j=nstart,nfin,2
        fi = (i+1)*y**i - i*y**(i-1)
        f = f + Opt(j) * fi
        i = i + 1
      END DO
    ELSE
      DO j=nstart,nfin,2
        fi = (i+1)*y**i - i*y**(i-1)
        f = f + Opt(j+1) * fi
        i = i + 1
      END DO
    END IF


!    PRINT *,'y',y,'dx',f,'Opt',Power1,Opt

  END FUNCTION RightWall



