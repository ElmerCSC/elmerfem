
!------------------------------------------------------------------------------
SUBROUTINE BulkAssembly( Model,Solver,dt,Transient, &
    Mass, Damp, Stiff, Force, InElement, nrow, ncol )
!------------------------------------------------------------------------------

   USE DefUtils

 
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t), TARGET :: Solver
  LOGICAL ::  Transient
  REAL(KIND=dp) :: dt
  REAL(KIND=dp) :: Stiff(:,:), Damp(:,:), Mass(:,:), Force(:)
  TYPE(Element_t), TARGET :: InElement
  INTEGER :: Nrow, Ncol 
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------

  TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
  TYPE(Nodes_t) :: Nodes
  REAL(KIND=dp), ALLOCATABLE :: Source(:), Cond(:), React(:)
  REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:)
  REAL(KIND=dp) :: SourceAtIP, CondAtIP, ReactAtIp, Weight, DetJ
  INTEGER :: i,j,t,p,q,n
  LOGICAL :: Visited = .FALSE., Found
  TYPE(Element_t), POINTER :: Element
  TYPE(ValueList_t), POINTER :: Material, Params
  CHARACTER(LEN=MAX_NAME_LEN) :: str


  
  SAVE Visited, Source, Cond, React, Nodes, Basis, dBasisdx, Params


  IF( .NOT. Visited ) THEN
     N = Solver % Mesh % MaxElementNodes    
     ALLOCATE( Source(n), Cond(n), React(n), Basis(n), dBasisdx(n,3), &
         Nodes % x(n), Nodes % y(n), Nodes % z(n) )     
     CALL info('PoissonBulkAssembly','1st time')
     Params => CurrentModel % Solver % Values
     Visited = .TRUE.
  END IF

  n = Nrow 
  CALL GetElementNodes( Nodes ) 
  Material => GetMaterial()  

  i = ListGetInteger( Params,'Block Matrix Row')
  j = ListGetInteger( Params,'Block Matrix Column')
  
  WRITE (str,'(A,I1,I1)') 'Conductivity ',i,j
  Cond(1:n) = GetReal( Material,str,Found)
  
  WRITE (str,'(A,I1,I1)') 'Reaction ',i,j
  React(1:n) = GetReal( Material,str,Found )

  IF( i == j ) THEN
    WRITE (str,'(A,I1)') 'Source ',i
    Source(1:n) = GetReal( Material,str,Found)
  ELSE
    Source(1:n) = 0.0_dp
  END IF

  IF( .FALSE. ) THEN
    PRINT *,'ij:',i,j
    PRINT *,'coeffs:',Cond(1),React(1),Source(1)
  END IF



  Element => GetCurrentElement()
  IntegStuff = GaussPoints( Element )

  DO t=1,IntegStuff % n
    
    Found = ElementInfo( Element, Nodes, IntegStuff % u(t), &
        IntegStuff % v(t), IntegStuff % w(t), detJ, Basis, dBasisdx )

    Weight = IntegStuff % s(t) * detJ

    CondAtIP = SUM( Basis(1:Nrow) * Cond(1:Nrow) )
    SourceAtIP = SUM( Basis(1:Nrow) * Source(1:Nrow) )
    ReactAtIP = SUM( Basis(1:Nrow) * React(1:Nrow) )
    
    DO p=1,n
      DO q=1,n
        STIFF(p,q) = STIFF(p,q) + Weight * CondAtIP * &
            SUM(dBasisdx(p,1:3)*dBasisdx(q,1:3))
        STIFF(p,q) = STIFF(p,q) + Weight * ReactAtIP * &
            Basis(p)*Basis(q)
      END DO
      FORCE(p) = FORCE(p) + Weight * SourceAtIP * Basis(p)
    END DO
  END DO
  
!------------------------------------------------------------------------------
END SUBROUTINE BulkAssembly
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
SUBROUTINE BoundaryAssembly( Model,Solver,dt,Transient, &
            Mass, Damp, Stiff, Force, InElement, nrow, ncol )
!------------------------------------------------------------------------------

  USE DefUtils
 
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t), TARGET :: Solver
  LOGICAL ::  Transient
  REAL(KIND=dp) :: dt
  REAL(KIND=dp) :: Stiff(:,:), Damp(:,:), Mass(:,:), Force(:)
  TYPE(Element_t), TARGET :: InElement
  INTEGER :: Nrow, Ncol 
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------

  TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
  TYPE(Nodes_t) :: Nodes
  REAL(KIND=dp), ALLOCATABLE :: Flux(:), Text(:), Coeff(:)
  REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:)
  REAL(KIND=dp) :: FluxAtIP, TextAtIP, CoeffAtIP, Weight, DetJ
  INTEGER :: i,j,t,p,q,n
  LOGICAL :: Visited = .FALSE., Found, Found2
  TYPE(ValueList_t), POINTER :: BC
  TYPE(Element_t), POINTER :: Element

  
  SAVE Visited, Flux, Text, Coeff, Nodes, Basis, dBasisdx


  IF( .NOT. Visited ) THEN
    N = Solver % Mesh % MaxElementNodes    
    ALLOCATE( Flux(n), Text(n), Coeff(n), Basis(n), dBasisdx(n,3), &
        Nodes % x(n), Nodes % y(n), Nodes % z(n) )     
    CALL info('PoissonBoundaryAssembly','1st time')
    Visited = .TRUE.
  END IF
  
  Element => GetCurrentElement()
  n = Nrow 
  CALL GetElementNodes( Nodes ) 
  
  BC => GetBC( )
  
  Flux(1:n) = GetReal( BC,'Flux',Found )
  Coeff(1:n) = GetReal( BC,'Coeffient',Found2 )
  IF(.NOT. (Found .OR. Found2)) RETURN
  Text(1:n) = GetReal( BC,'External',Found )


  IntegStuff = GaussPoints( Element )
  
  DO t=1,IntegStuff % n
    
    Found = ElementInfo( Element, Nodes, IntegStuff % u(t), &
        IntegStuff % v(t), IntegStuff % w(t), detJ, Basis, dBasisdx )
    
    Weight = IntegStuff % s(t) * detJ
    
    FluxAtIP = SUM( Basis(1:Nrow) * Flux(1:Nrow) )
    CoeffAtIP = SUM( Basis(1:Nrow) * Coeff(1:Nrow) )
    TextAtIP = SUM( Basis(1:Nrow) * Text(1:Nrow) )
    
    DO p=1,n
      DO q=1,n
        STIFF(p,q) = STIFF(p,q) + Weight * CoeffAtIP * &
            Basis(p)*Basis(q)
      END DO
      FORCE(p) = FORCE(p) + Weight * Basis(p) * &
          ( FluxAtIP + CoeffAtIP * TextAtIP )
    END DO
  END DO

!------------------------------------------------------------------------------
END SUBROUTINE BoundaryAssembly
!------------------------------------------------------------------------------

