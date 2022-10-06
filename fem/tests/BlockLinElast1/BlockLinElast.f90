
!------------------------------------------------------------------------------
SUBROUTINE BulkAssembly( Model,Solver,dt,Transient, &
    Mass, Damp, Stiff, Force, InElement, nrow, ncol )
!------------------------------------------------------------------------------

   USE DefUtils

 
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t) :: Solver
  LOGICAL ::  Transient
  REAL(KIND=dp) :: dt
  TYPE(Element_t):: InElement
  INTEGER :: Nrow, Ncol 
  REAL(KIND=dp) :: Stiff(:,:), Damp(:,:), Mass(:,:), Force(:)
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------

  TYPE(GaussIntegrationPoints_t) :: IntegStuff
  TYPE(Nodes_t) :: Nodes
  REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:)
  REAL(KIND=dp), ALLOCATABLE :: Source(:)
  REAL(KIND=dp) :: Young, Poisson, Lame1, Lame2
  REAL(KIND=dp) :: SourceAtIP, Weight, DetJ
  INTEGER :: i,j,k,t,p,q,n, nd,dim
  LOGICAL :: Visited = .FALSE., Found
  TYPE(Element_t), POINTER :: Element
  TYPE(ValueList_t), POINTER :: Material, Params

  SAVE Visited, Source, Nodes, Basis, dBasisdx, Params

  dim = CoordinateSystemDimension()

  IF( .NOT. Visited ) THEN
     n = Solver % Mesh % MaxElementDOFs
     ALLOCATE( Source(n), Basis(n), dBasisdx(n,3) )
     CALL info('PoissonBulkAssembly','1st time')
     Params => GetSolverParams()
     Visited = .TRUE.
  END IF

  n  = GetElementNOFNodes()
  nd = GetElementNOFDOFs()
  CALL GetElementNodes( Nodes ) 
  Material => GetMaterial()  

  i = GetInteger( Params,'Block Matrix Row')
  j = GetInteger( Params,'Block Matrix Column')
  Source(1:n) = GetReal( Material,'Source',Found)

  Young   = GetCReal( Material, 'Youngs modulus', Found )
  Poisson = GetCReal( Material, 'Poisson Ratio',  Found )

  Lame1 = Young*Poisson/((1+Poisson)*(1-2*Poisson))
  Lame2 = Young/(2*(1+Poisson))

  Element => GetCurrentElement()
  IntegStuff = GaussPoints( Element )


  DO t=1,IntegStuff % n
    Found = ElementInfo( Element, Nodes, IntegStuff % u(t), &
        IntegStuff % v(t), IntegStuff % w(t), detJ, Basis, dBasisdx )

    Weight = IntegStuff % s(t) * detJ

    SourceAtIP = SUM( Basis(1:n) * Source(1:n) )

    IF ( i==j ) THEN
      STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + &
            Weight*Lame2*MATMUL(dBasisdx,TRANSPOSE(dBasisdx))
      FORCE(1:nd) = FORCE(1:nd) + Weight * SourceAtIP * Basis(1:nd)
    END IF
    
    DO p=1,nd
      DO q=1,nd
        STIFF(p,q) = STIFF(p,q) + &
              Weight*Lame1*dBasisdx(q,j)*dBasisdx(p,i)
        STIFF(p,q) = STIFF(p,q) + &
              Weight*Lame2*dBasisdx(q,i)*dBasisdx(p,j)
      END DO
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
  Coeff(1:n) = GetReal( BC,'Coefficient',Found2 )
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

