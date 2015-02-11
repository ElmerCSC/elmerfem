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
! *  Authors: Peter Råback, Juha Ruokolainen
! *  Email:   Peter.Raback@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 5.9.2013
! *
! *****************************************************************************/



!------------------------------------------------------------------------------
!> Subroutine that helps to give current sources for closed coils. Normally there
!> is the difficulty that the potential in such a coil should be discontinuous. 
!> In this solver two different potential fields are solved (say a left and a right 
!> field). The two fields have different boundary conditions that are automatically 
!> set on the bulk nodes so that a current is induced. The solution settles typically
!> on the other side such that the union of the solutions is always reasonable.
!> The main assumption is that the coil axis is aligned with the z-axis. 
!
!> The module includes a special user function that returns the potential 
!> at the good side. 
!
!> \ingroup Solvers
!------------------------------------------------------------------------------


SUBROUTINE CoilSolver_init( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: Params 
  INTEGER :: dim
  LOGICAL :: Found  

  dim = CoordinateSystemDimension()
  Params => Solver % Values

  IF( .NOT. ListCheckPresent( Params,'Variable') ) THEN
    CALL ListAddString( Params,'Variable','-nooutput CoilTmp')
  END IF
  CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params),&
      'CoilPot')
  CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params),&
      'CoilFix')

  IF( GetLogical( Params,'Coil Closed', Found ) ) THEN
     CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params),&
          'CoilPotB')
  END IF

  CALL ListAddString( Params,&
       NextFreeKeyword('Exported Variable',Params),&
       'CoilCurrent[CoilCurrent:'//TRIM(I2S(dim))//']')

! Loads are needed to compute the induced currents in a numerically optimal way
  CALL ListAddLogical( Params,'Calculate Loads',.TRUE.)

END SUBROUTINE CoilSolver_init



!------------------------------------------------------------------------------
SUBROUTINE CoilSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Element_t),POINTER :: Element
  REAL(KIND=dp) :: Norm, s
  INTEGER :: i,j,k,n, nb, nd, t, active, iter, Part, sgn, nsize, CoilParts, &
       MaxNonlinIter, dim,dimi
  INTEGER, POINTER :: Perm(:), Set(:)
  INTEGER, ALLOCATABLE, TARGET :: SetA(:), SetB(:)
  TYPE(Matrix_t), POINTER :: StiffMatrix
  REAL(KIND=dp), POINTER :: ForceVector(:)
  TYPE(Variable_t), POINTER :: PotVar, FixVar, SolVar, FluxVar, LoadVar, DistVar
  TYPE(Variable_t), POINTER :: PotVarA,PotVarB
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: Params 
  REAL(KIND=dp) :: possum, negsum, DesiredCoilCurrent, DesiredCurrentDensity, &
      CoilCrossSection,InitialCurrent, Coeff, val, x0
  LOGICAL :: Found, CoilClosed, CoilAnisotropic, UseDistance, FixConductivity, &
      NormalizeCurrent, GotCurr, GotDens

 !------------------------------------------------------------------------------

  CALL Info('CoilSolver','--------------------------------------')
  CALL Info('CoilSolver','Solving current distribution in a coil')
  CALL Info('CoilSolver','--------------------------------------')

  Params => Solver % Values

  nsize = SIZE( Solver % Variable % Values ) 
  ALLOCATE( SetA(nsize), SetB(nsize) )

  StiffMatrix => Solver % Matrix
  ForceVector => Solver % Matrix % Rhs
  Mesh => Solver % Mesh
  SolVar => Solver % Variable
  Perm => SolVar % Perm

  dim = CoordinateSystemDimension()

  FixConductivity = GetLogical( Params,'Coil Conductivity Fix', Found ) 

  CoilAnisotropic = GetLogical( Params,'Coil Anisotropic', Found )

  UseDistance = GetLogical( Params,'Use Wall Distance',Found)
  IF( UseDistance ) THEN
     DistVar => VariableGet( Mesh % Variables,'Wall Distance' )
     IF( .NOT. ASSOCIATED( DistVar ) ) THEN
        CALL Fatal('CoilSolver','> Wall Distance < not associated!')
     END IF
  END IF

  CoilClosed = GetLogical( Params,'Coil Closed', Found )
  IF( CoilClosed ) THEN
     CoilParts = 2
  ELSE
     IF( .NOT. ListGetLogicalAnyBC( Model,'Coil Start') ) THEN
       CALL Info('CoilSolver','Assuming coil that is not closed',Level=3)
       CALL Fatal('CoilSolver','> Coil Start < must be defined on some BC')
     END IF
     IF( .NOT. ListGetLogicalAnyBC( Model,'Coil End') ) THEN
       CALL Info('CoilSolver','Assuming coil that is not closed',Level=3)
       CALL Fatal('CoilSolver','> Coil End < must be defined on some BC')
     END IF
     CoilParts = 1
  END IF

  DesiredCoilCurrent = ListGetCReal( Params,'Desired Coil Current',GotCurr )
  IF(.NOT. GotCurr ) DesiredCoilCurrent = 1.0_dp

  DesiredCurrentDensity = ListGetCReal( Params,'Desired Current Density',GotDens)
  IF(.NOT. GotDens ) DesiredCurrentDensity = 1.0_dp

  ! If we know the coil cross section we can relate the wishes in total 
  ! current through cross section and current density. 
  CoilCrossSection = ListGetCReal( Params,'Coil Cross Section',Found ) 
  IF( Found ) THEN
    IF( GotCurr ) THEN
      DesiredCurrentDensity = DesiredCoilCurrent / CoilCrossSection
      GotDens = .TRUE.
    ELSE IF( GotDens ) THEN
      DesiredCoilCurrent = DesiredCurrentDensity * CoilCrossSection 
      GotCurr = .TRUE.
    END IF
  END IF

  NormalizeCurrent = ListGetLogical( Params,'Normalize Coil Current',Found ) 

  MaxNonlinIter = GetInteger( Params,'Nonlinear System Max Iterations',Found )
  IF(.NOT. Found ) THEN
    IF( FixConductivity ) THEN
      IF( CoilAnisotropic ) THEN
        MaxNonlinIter = 3
      ELSE
        MaxNonlinIter = 2
      END IF
    ELSE
      MaxNonlinIter = 1
    END IF
  END IF

  IF( MaxNonlinIter > 1 .AND. .NOT. FixConductivity ) THEN
    CALL Warn('CoilSolver','No conductivity fixing, setting nonlinear iterations to one!')
  END IF

  PotVarA => VariableGet( Mesh % Variables,'CoilPot' )
  IF( .NOT. ASSOCIATED( PotVarA ) ) THEN
    CALL Fatal('CoilSolver','CoilPot not associated!')
  END IF
  IF( CoilParts == 2 ) THEN
    PotVarB => VariableGet( Mesh % Variables,'CoilPotB' )
    IF( .NOT. ASSOCIATED( PotVarB ) ) THEN
      CALL Fatal('CoilSolver','CoilPotB not associated!')
    END IF
  END IF

  FixVar => VariableGet( Mesh % Variables,'CoilFix' )
  IF( .NOT. ASSOCIATED( FixVar ) ) THEN
    CALL Fatal('CoilSolver','CoilFix not associated!')
  END IF

  ! Get the loads
  LoadVar => VariableGet( Mesh % Variables,&
      TRIM(SolVar % Name)//' Loads' )
  IF( .NOT. ASSOCIATED( LoadVar ) ) THEN
    CALL Fatal('CoilSolver','> '//TRIM(SolVar % Name)//' < Loads not associated!')
  END IF
  

  ! Choose nodes where the Dirichlet values are set. 
  IF( CoilClosed ) THEN
    Set => SetA
    CALL ChooseFixedBulkNodes(Set,1)
    Set => SetB
    CALL ChooseFixedBulkNodes(Set,2)
    x0 = ListGetCReal( Params,'Coil x0',Found )
  ELSE
    Set => SetA
    CALL ChooseFixedEndNodes(Set)
  END IF

     
  DO iter=1,MaxNonlinIter
      
    IF( iter == 1 ) THEN
      FixVar % Values = 1.0_dp
    ELSE
      CALL Info('CoilSolver','Fixing the conductivity field')
      
      CALL DefaultInitialize()
      
      Active = GetNOFActive()
      DO t=1,Active
        Element => GetActiveElement(t)
        n  = GetElementNOFNodes()
        nd = GetElementNOFDOFs()           
        CALL LocalFixMatrix(  Element, n, nd )
      END DO
      
      ! Solve the fixing field
      !--------------------------
      CALL ListAddLogical( Params,'Skip Compute Nonlinear Change',.TRUE.) 
      Norm = DefaultSolve()
      FixVar % Values = SolVar % Values
      
      PRINT *,'Fix Range:',MINVAL( FixVar % Values), MAXVAL( FixVar % Values )
    END IF

    
    CALL Info('CoilSolver','Computing the dummy potential field')

    ! For closed coils the solution is computed in two parts
    DO Part = 1,CoilParts
      
      ! For closed coils the solution is computed in two parts
      IF( Part == 1 ) THEN        
        Set => SetA
        PotVar => PotVarA
      ELSE
        Set => SetB
        PotVar => PotVarB
      END IF
      
      CALL DefaultInitialize()

      Active = GetNOFActive()          
      DO t=1,Active
        Element => GetActiveElement(t)
        n  = GetElementNOFNodes()
        nd = GetElementNOFDOFs()
        CALL LocalPotMatrix(  Element, n, nd )
      END DO
      
      ! This routine is needed to save the BulkValues
      ! Bulk values are needed to compute the currents
      CALL DefaultFinishBulkAssembly()
      CALL DefaultFinishAssembly()
      
      ! Set the potential values on the nodes
      ! Values 0/1 are used as the potential is later scaled. 
      DO i = 1,nsize
        sgn = Set(i)
        IF( sgn == 0 ) CYCLE
        
        IF( sgn > 0 ) THEN
          val = 1.0_dp
        ELSE
          val = 0.0_dp
        END IF
        
        s = StiffMatrix % Values(StiffMatrix % Diag(i))
        ForceVector(i) = val * s
        CALL ZeroRow( StiffMatrix,i )
        CALL SetMatrixElement( StiffMatrix,i,i,1.0d0*s )     
      END DO
      
      ! Solve the potential field
      !--------------------------
      CALL ListAddLogical( Params,'Skip Compute Nonlinear Change',.FALSE.) 

      ! Here solve actually the potential
      !--------------------------------------
      Norm = DefaultSolve()
      
      PotVar % Values = SolVar % Values
      PRINT *,'Pot Range:',MINVAL( PotVar % Values), MAXVAL( PotVar % Values )


      ! Get the nodal loads i.e. nodal currents in this case
      LoadVar => VariableGet( Mesh % Variables,&
          TRIM(SolVar % Name)//' Loads' )
      IF( .NOT. ASSOCIATED( LoadVar ) ) THEN
        CALL Fatal('CoilSolver','> '//TRIM(SolVar % Name)//' < Loads not associated!')
      END IF
      
      ! Evaluate the positive and negative currents
      ! As the total current vanishes these should be roughly the same 
      ! but with different signs. 
      possum = 0.0_dp
      negsum = 0.0_dp
      DO i=1,nsize
        sgn = Set(i)
        
        ! Note that the narrow part of the gap is omitted and here only 
        ! currents related to the outher parts of the gap (with |sgn|=2) 
        ! are accounted for. 
        IF( sgn == 2 ) THEN
          possum = possum + LoadVar % Values(i)
        ELSE IF( sgn == -2 ) THEN
          negsum = negsum + LoadVar % Values(i)
        END IF
      END DO
    
      PRINT *,'current estimate:',Part,negsum,possum
      InitialCurrent = ( possum - negsum ) / 2.0
    
      ! Scale the potential such that the current is as desired
      ! The current scales linearly with the potential.
      Coeff = DesiredCoilCurrent / InitialCurrent
      PRINT *,'multiplying with:',Coeff,InitialCurrent
      
!      CorrFactor( Part ) = Coeff 
      PotVar % Values = Coeff * PotVar % Values     

    END DO
  END DO


!  PotVarA % Values = CorrFactor(1) * PotVarA % Values
!  PotVarB % Values = CorrFactor(2) * PotVarB % Values


  ! Compute the current
  !---------------------------------
  DO Part = 1,CoilParts
    
    ! The solution is computed separately for each component
    !--------------------------------------------------------
    DO dimi = 1,dim
      CALL Info('CoilSolver','Computing the current components')
      
      CALL DefaultInitialize()
      Active = GetNOFActive()          
      DO t=1,Active
        Element => GetActiveElement(t)
        n  = GetElementNOFNodes()
        nd = GetElementNOFDOFs()           
        CALL LocalFluxMatrix(  Element, n, nd, dimi )
      END DO
      
      ! Solve the flux in direction dimi
      !--------------------------------------
      CALL ListAddLogical( Params,'Skip Compute Nonlinear Change',.TRUE.) 
      Norm = DefaultSolve()
      
      ! If the coil is computed in two parts then 
      ! for the other part pick only the values which are on the better half
      FluxVar => VariableGet( Mesh % Variables,'CoilCurrent '//TRIM(I2S(dimi)) )
      IF( .NOT. ASSOCIATED( FluxVar ) ) THEN
        CALL Fatal('CoilSolver','CoilCurrent not associated!')
      END IF

      FluxVar % Values = SolVar % Values

    END DO    
  END DO
  
  IF( NormalizeCurrent ) THEN
    CALL NormalizeCurrentDensity() 
  END IF


  DEALLOCATE( SetA, SetB ) 

  CALL Info('CoilSolver','All done')
   
 

CONTAINS 


  ! Chooses bulk nodes which are used to set the artificial boundary conditions
  ! in the middle of the coil.
  !----------------------------------------------------------------------------  
  SUBROUTINE ChooseFixedBulkNodes( Set, SetNo )
    
    INTEGER :: SetNo
    INTEGER, POINTER :: Set(:)

    LOGICAL :: Mirror 
    TYPE(Mesh_t), POINTER :: Mesh
    REAL(KIND=dp) :: x,y,z,x0,y0,rad2deg,fii,dfii,dy
    REAL(KIND=dp) :: MinCoord(3),MaxCoord(3)
    INTEGER :: i,j,k,nminus,nplus
    LOGICAL :: Found

    Mirror = ( SetNo == 2 )

    rad2deg = 180.0_dp / PI

    Mesh => Solver % Mesh

    ! The angle of acceptable nodes 
    ! 90 degs effectively chooses the right half
    dfii = ListGetCReal( Params,'Coil dfii',Found)
    IF(.NOT. Found ) dfii = 90.0

    ! The maximum coordinate difference for an acceptable node
    ! The larger the value the more there will be nodes in the set.
    ! There should be enough nodes so that the BC is good one,
    ! but not too many either. 

    MinCoord = HUGE( MinCoord )
    MaxCoord = -HUGE( MaxCoord )

    DO i=1,Mesh % NumberOfNodes
      IF( Perm(i) == 0 ) CYCLE

      MinCoord(1) = MIN( MinCoord(1), Mesh % Nodes % x(i) )      
      MaxCoord(1) = MAX( MaxCoord(1), Mesh % Nodes % x(i) )

      MinCoord(2) = MIN( MinCoord(2), Mesh % Nodes % y(i) )      
      MaxCoord(2) = MAX( MaxCoord(2), Mesh % Nodes % y(i) )

      MinCoord(3) = MIN( MinCoord(3), Mesh % Nodes % z(i) )      
      MaxCoord(3) = MAX( MaxCoord(3), Mesh % Nodes % z(i) )
    END DO

    PRINT *,'Min Coords:',MinCoord
    PRINT *,'Max Coords:',MaxCoord


    dy = ListGetCReal( Params,'Coil dy',Found)
    IF(.NOT. Found ) THEN
      dy = 0.2 * ( MaxCoord(2) - MinCoord(2) )
    END IF

    x0 = ListGetCReal( Params,'Coil x0',Found)
    IF(.NOT. Found ) x0 = ( MaxCoord(1) + MinCoord(1) ) / 2

    y0 = ListGetCReal( Params,'Coil y0',Found)
    IF(.NOT. Found ) y0 = ( MaxCoord(2) + MinCoord(2) ) / 2

    PRINT *,'Center:',X0,Y0

    ! Add this also to the Simulation section as the UDF needs it
    IF( Found ) CALL ListAddConstReal( Model % Simulation,'Coil x0',y0)
    Set = 0
    
    DO i=1,Mesh % NumberOfNodes
      j = Perm(i)

      IF( j == 0 ) CYCLE
      
      x = Mesh % Nodes % x(i) - x0
      y = Mesh % Nodes % y(i) - y0
      z = Mesh % Nodes % z(i)

      IF( mirror ) THEN      
        x = -x
        y = -y
      END IF

      fii = rad2deg * ATAN2( y, x )
      IF( fii > 180.0 ) fii = fii - 360.0
      
      IF( ABS( fii ) > dfii ) CYCLE
      
      IF( ABS( y ) > dy ) CYCLE

      ! Values with abs 1 indicate the narrow band that is omitted when computing the currents
      ! Values with abs 2 indicate the wide band
      IF( y > 0 ) THEN
        Set(j) = 1
        IF( y > dy / 2 ) Set(j) = 2
      ELSE
        Set(j) = -1
        IF( y < -dy / 2 ) Set(j) = -2 
      END IF      
    END DO
    
    ! Just count the nodes set 
    nplus = COUNT( Set > 0 ) 
    nminus = COUNT( Set < 0 ) 
    CALL Info('CoilSolver','Set'//TRIM(I2S(SetNo))//' : '&
        //TRIM(I2S(nplus))//' +nodes and ' &
        //TRIM(I2S(nminus))//' -nodes')

    IF( nplus == 0 .OR. nminus == 0 ) THEN
      CALL Fatal('CoilSolver','Cannot set Dirichlet conditions with this set')
    END IF


  END SUBROUTINE ChooseFixedBulkNodes


  ! Choose end nodes as assingled by "Coil Start" and "Coil End" flags.
  ! The result of this imitate the previous routine in order to be able
  ! to use the same way to computed the resulting currents.
  !--------------------------------------------------------------------
  SUBROUTINE ChooseFixedEndNodes( Set )
    
    INTEGER, POINTER :: Set(:)

    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER, POINTER :: Indexes(:)
    INTEGER :: t,i,j,k,nminus,nplus
    TYPE(ValueList_t), POINTER :: BC
    LOGICAL :: Found

    Set = 0

    Mesh => Solver % Mesh
    
    DO t=Mesh % NumberOfBulkElements + 1, &
         Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
       
       Element => Mesh % Elements(t)
       Model % CurrentElement => Element
       Indexes => Element % NodeIndexes

       IF( ANY( Perm( Indexes ) == 0 ) ) CYCLE

       BC => GetBC(Element)
       IF( ListGetLogical( BC,'Coil Start',Found ) ) THEN
          Set( Perm( Indexes ) ) = 2
!          PRINT *,'Coord1: ',Mesh % Nodes % x( Indexes(1) ), &
!              Mesh % Nodes % y( Indexes(1) )
       ELSE IF( ListGetLogical( BC,'Coil End',Found ) ) THEN
          Set( Perm( Indexes ) ) = -2 
!          PRINT *,'Coord2: ',Mesh % Nodes % x( Indexes(1) ), &
!              Mesh % Nodes % y( Indexes(1) )
       END IF
    END DO
          
    ! Just count the nodes set 
    nplus = COUNT( Set > 0 ) 
    nminus = COUNT( Set < 0 ) 
    CALL Info('CoilSolver','Found '&
        //TRIM(I2S(nplus))//' start nodes and ' &
        //TRIM(I2S(nminus))//' end nodes')

    IF( nplus == 0 .OR. nminus == 0 ) THEN
      CALL Fatal('CoilSolver','Cannot set Dirichlet conditions with this set')
    END IF

  END SUBROUTINE ChooseFixedEndNodes



! Assembly the potential equation related to the dummy potential 
!------------------------------------------------------------------------------
  SUBROUTINE LocalPotMatrix( Element, n, nd )
!------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,Weight,LoadAtIP,GradAtIp(3), &
         FixAtIp, AbsGradAtIp, CondAtIp(3), DistGradAtIp(3),AbsDistGradAtIp,AbsCondAtIp
    REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd), NodalPot(n), &
         NodalFix(n), NodalDist(n),DotProd
    LOGICAL :: Stat,Found
    INTEGER :: i,t,p,q
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------

    CALL GetElementNodes( Nodes )
    STIFF = 0._dp
    FORCE = 0._dp

    IF( CoilParts == 1 ) THEN
      NodalPot(1:n) = PotVar % Values( Perm( Element % NodeIndexes ) )
    ELSE      
      IF( MINVAL( Nodes % x(1:n)) - x0 > 0.0 ) THEN
        NodalPot(1:n) = PotVarB % Values( Perm( Element % NodeIndexes ) )
      ELSE
        NodalPot(1:n) = PotVarA % Values( Perm( Element % NodeIndexes ) )
      END IF
    END IF

    IF( iter > 1 ) THEN
      NodalFix(1:n) = FixVar % Values( Perm( Element % NodeIndexes(1:n) ) )
    END IF
    
    IF( UseDistance ) THEN
      NodalDist(1:n) = DistVar % Values( DistVar % Perm( Element % NodeIndexes(1:n) ) )
    END IF
    
    ! Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis, dBasisdx )

      IF( iter > 1 ) THEN
         FixAtIp = SUM( Basis(1:n) * NodalFix(1:n) )
         IF( CoilAnisotropic ) THEN
            DO i=1,3
               GradAtIP(i) = SUM( dBasisdx(1:n,i) * NodalPot(1:n) )
            END DO
            AbsGradAtIp = SQRT( SUM( GradAtIp ** 2 ) )
            CondAtIp = FixAtIp * ABS( GradAtIp ) / AbsGradAtIp
         ELSE
            CondAtIp = FixAtIp
         END IF
      ELSE
         CondAtIp = 1.0_dp
      END IF
 
      ! If distance from the tangential walls is used then
      ! remove from conductivity the normal component.
      ! The normal vector is defined as the gradient of distance.
      ! As it can in center of coil be less than zero because it is 
      ! poorly defined, normalize the vector only if its larger. 
      ! At center the model could therefore be more isotropic.
      !-------------------------------------------------------------
      IF( UseDistance ) THEN
         DO i=1,3
            DistGradAtIP(i) = SUM( dBasisdx(1:n,i) * NodalDist(1:n) )
         END DO
         AbsDistGradAtIp = SQRT( SUM( DistGradAtIp ** 2 ) )

         AbsCondAtIp = SQRT(SUM( CondAtIP**2 ))
         DotProd = SUM( CondAtIp * DistGradAtIp ) / ( AbsDistGradAtIp * SQRT(3.0) )

         IF( AbsDistGradAtIp > 1.0_dp ) THEN
            DistGradAtIp = DistGradAtIp / AbsDistGradAtIp
         END IF

         CondAtIp = CondAtIp - DotProd * DistGradAtIp

         IF( ANY( CondAtIp < 0 ) ) THEN
            PRINT *,'DistGrad',DistGradAtIp,AbsDistGradAtIp
            PRINT *,'Cond',CondAtIp,DotProd,AbsCondAtIp
         END IF

      END IF

      Weight = IP % s(t) * DetJ

      ! diffusion term (Cond*grad(u),grad(v)):
      ! -----------------------------------
      DO p=1,nd
         DO q=1,nd
            DO i=1,dim
               STIFF(p,q) = STIFF(p,q) + Weight * &
                    CondAtIp(i) * dBasisdx(p,i) * dBasisdx(q,i) 
            END DO
         END DO
      END DO
    END DO

    CALL DefaultUpdateEquations(STIFF,FORCE)
!------------------------------------------------------------------------------
  END SUBROUTINE LocalPotMatrix
!------------------------------------------------------------------------------


! Assembly the equation for the fixing coefficient
!------------------------------------------------------------------------------
  SUBROUTINE LocalFixMatrix( Element, n, nd )
!------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: diff_coeff(n), AbsGradAtIp, Weight, Dreg, AbsCondAtIp
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP,GradAtIp(3),DistGradAtIp(3)
    REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd),NodalPot(nd),NodalDist(nd)
    LOGICAL :: Stat,Found
    INTEGER :: i,t,p,q
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------

    CALL GetElementNodes( Nodes )
    STIFF = 0._dp
    FORCE = 0._dp

    IF( CoilParts == 1 ) THEN
      NodalPot(1:n) = PotVar % Values( Perm( Element % NodeIndexes ) )
    ELSE
      IF( MINVAL( Nodes % x(1:n)) - x0 > 0.0 ) THEN
        NodalPot(1:n) = PotVarB % Values( Perm( Element % NodeIndexes ) )
      ELSE
        NodalPot(1:n) = PotVarA % Values( Perm( Element % NodeIndexes ) )
      END IF
    END IF

    IF( ALL( ABS( NodalPot(1:n) ) < 1.0e-8 ) ) THEN
      PRINT *,'NodalPot',NodalPot
      CALL Fatal('','')
    END IF

    IF( UseDistance ) THEN
       NodalDist(1:n) = DistVar % Values( DistVar % Perm( Element % NodeIndexes ) )
    END IF
    
    diff_coeff(1:n) = GetReal(Params,'CFix Diffusion',Found)

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis, dBasisdx )

      ! By construction the r.h.s. should be one
      !------------------------------------------
      LoadAtIP = 1.0_dp

      DO i=1,3
         GradAtIP(i) = SUM( dBasisdx( 1:nd, i) * NodalPot(1:nd) )
      END DO

      IF( UseDistance ) THEN
         DO i=1,3
            DistGradAtIP(i) = SUM( dBasisdx(1:n,i) * NodalDist(1:n) )
         END DO
         GradAtIp = GradAtIp - SUM( GradAtIp * DistGradAtIp ) * DistGradAtIp 
      END IF

      AbsGradAtIp = SQRT( SUM( GradAtIp**2) )      

      Dreg = SUM( Basis(1:nd) * Diff_Coeff(1:nd) )
      
      Weight = IP % s(t) * DetJ

      ! diffusion term (D*grad(u),grad(v)):
      ! -----------------------------------
      STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + Weight * &
             Dreg * MATMUL( dBasisdx, TRANSPOSE( dBasisdx ) )

      DO p=1,nd
        DO q=1,nd
          STIFF(p,q) = STIFF(p,q) + Weight * AbsGradAtIp * Basis(q) * Basis(p)
        END DO
      END DO

      FORCE(1:nd) = FORCE(1:nd) + Weight * LoadAtIP * Basis(1:nd)
    END DO

!    CALL LCondensate( nd-nb, nb, STIFF, FORCE )
    CALL DefaultUpdateEquations(STIFF,FORCE)
!------------------------------------------------------------------------------
  END SUBROUTINE LocalFixMatrix
!------------------------------------------------------------------------------


! Assembly the equation for the flux component
!------------------------------------------------------------------------------
  SUBROUTINE LocalFluxMatrix( Element, n, nd, dimi  )
!------------------------------------------------------------------------------
    INTEGER :: n, nd, dimi
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),Weight,DetJ,LoadAtIP,FixAtIp, &
         GradAtIp(3),AbsGradAtIp,CondAtIp(3),DistGradAtIp(3),AbsDistGradAtIp, &
         AbsCondAtIp,DotProd
    REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd),NodalPot(nd),&
         NodalFix(nd),NodalDist(nd)
    LOGICAL :: Stat
    INTEGER :: i,t,p,q
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t) :: Nodes

    SAVE Nodes
!------------------------------------------------------------------------------

    CALL GetElementNodes( Nodes )
    STIFF = 0._dp
    FORCE = 0._dp

    IF( CoilParts == 1 ) THEN
      NodalPot(1:n) = PotVar % Values( Perm( Element % NodeIndexes ) )
    ELSE
      IF( MINVAL( Nodes % x(1:n)) - x0 > 0.0 ) THEN
        NodalPot(1:n) = PotVarB % Values( Perm( Element % NodeIndexes ) )
      ELSE
        NodalPot(1:n) = PotVarA % Values( Perm( Element % NodeIndexes ) )
      END IF
    END IF

    IF( iter > 1 ) THEN
      NodalFix(1:n) = FixVar % Values( Perm( Element % NodeIndexes ) )
    END IF

    IF( UseDistance ) THEN
       NodalDist(1:n) = DistVar % Values( DistVar % Perm( Element % NodeIndexes ) )
    END IF

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis, dBasisdx )

      DO i=1,3
         GradAtIP(i) = SUM( dBasisdx( 1:nd, i) * NodalPot(1:nd) )
      END DO
      AbsGradAtIp = SQRT( SUM( GradAtIp**2) )            

      IF( iter == 1 ) THEN
        CondAtIp = 1.0_dp
      ELSE
        FixAtIP = SUM( Basis( 1:nd) * NodalFix(1:nd) )
        IF( CoilAnisotropic ) THEN
          CondAtIp = FixAtIp * ABS( GradAtIp ) / AbsGradAtIp
        ELSE
          CondAtIp = FixAtIp 
        END IF
      END IF


      IF( UseDistance ) THEN
         DO i=1,3
            DistGradAtIP(i) = SUM( dBasisdx(1:n,i) * NodalDist(1:n) )
         END DO
         AbsDistGradAtIp = SQRT( SUM( DistGradAtIp ** 2 ) )

         AbsCondAtIp = SQRT(SUM( CondAtIP**2 ))
         DotProd = SUM( CondAtIp * DistGradAtIp ) / ( AbsDistGradAtIp * SQRT(3.0) )

         IF( AbsDistGradAtIp > 1.0_dp ) THEN
            DistGradAtIp = DistGradAtIp / AbsDistGradAtIp
         END IF

         CondAtIp = CondAtIp - DotProd * DistGradAtIp
      END IF

      LoadAtIp = CondAtIp(dimi) * GradAtIp(dimi) 
      
      Weight = IP % s(t) * DetJ

      DO p=1,nd
        DO q=1,nd
          STIFF(p,q) = STIFF(p,q) + Weight * Basis(q) * Basis(p)
        END DO
      END DO

      FORCE(1:nd) = FORCE(1:nd) + Weight * LoadAtIP * Basis(1:nd)
    END DO

!    CALL LCondensate( nd-nb, nb, STIFF, FORCE )
    CALL DefaultUpdateEquations(STIFF,FORCE)
!------------------------------------------------------------------------------
  END SUBROUTINE LocalFluxMatrix
!------------------------------------------------------------------------------


! Normalize the current density to a given length
!------------------------------------------------------------------------------
  SUBROUTINE NormalizeCurrentDensity( )
!------------------------------------------------------------------------------
    INTEGER :: n, nd, dimi
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Weight,DetJ,LocalCurr(3),AbsCurr,ScaleCurr,TotCurr,TotVol
    REAL(KIND=dp), ALLOCATABLE :: Basis(:),NodalCurr(:,:)
    LOGICAL :: Stat
    INTEGER :: i,elem,t
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t) :: Nodes

    SAVE Nodes
!------------------------------------------------------------------------------

    CALL Info('CoilSolver','Normalizing current density to a constant value')

    n = Mesh % MaxElementNodes
    ALLOCATE( Basis(n), NodalCurr(3,n) )

    
    FluxVar => VariableGet( Mesh % Variables,'CoilCurrent')
    IF( .NOT. ASSOCIATED( FluxVar ) ) THEN
      CALL Fatal('CoilSolver','CoilCurrent not associated!')
    END IF
    
    ! Current is already scaled, but we don't know the desired current density
    ! Let's assume that the average current density is a good candidate
    !--------------------------------------------------------------------------
    IF( .NOT. GotDens .AND. GotCurr ) THEN

      TotVol = 0.0_dp
      TotCurr = 0.0_dp
    
      Active = GetNOFActive()          
      DO elem=1,Active
        Element => GetActiveElement(elem)
        n  = GetElementNOFNodes()
        nd = GetElementNOFDOFs()           
        
        CALL GetElementNodes( Nodes )
        
        DO dimi=1,dim
          NodalCurr(dimi,1:n) = FluxVar % Values( Perm( Element % NodeIndexes ) )
        END DO
        
        !Numerical integration:
        !----------------------
        IP = GaussPoints( Element )
        DO t=1,IP % n
          stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis )
          
          DO dimi=1,dim
            LocalCurr(dimi) = SUM( Basis(1:n) * NodalCurr(dimi,1:n) )
          END DO
          AbsCurr = SQRT( SUM( LocalCurr(1:dim)**2 ) )
          
          Weight = IP % s(t) * DetJ
          
          TotVol = TotVol + Weight
          TotCurr = TotCurr + Weight * AbsCurr
        END DO
      END DO

      DesiredCurrentDensity = TotCurr / TotVol
      
      PRINT *,'Average current density:',DesiredCurrentDensity      
    END IF
    

    DO i=1, SIZE(FluxVar % Values ) / dim
      LocalCurr(1:dim) = FluxVar % Values( dim*(i-1)+1: dim*(i-1)+dim )
      AbsCurr = SQRT( SUM( LocalCurr(1:dim) ** 2 ) )
      IF( AbsCurr > TINY( AbsCurr ) ) THEN
        ScaleCurr = DesiredCurrentDensity / AbsCurr 
        FluxVar % Values( dim*(i-1)+1: dim*(i-1)+dim ) = &
            ScaleCurr * LocalCurr(1:dim)
      END IF
    END DO
    
!------------------------------------------------------------------------------
  END SUBROUTINE NormalizeCurrentDensity
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE CoilSolver
!------------------------------------------------------------------------------



! Returns the coil potential in the case of closed coil if it should be rather 
! used as a potential value than a current. The potential consists of two parts
! and is not continuous but it's gradient should be continuous. 
!------------------------------------------------------------------------------

  FUNCTION UnitedPotential( Model, n, t ) RESULT(f)

    USE DefUtils

    IMPLICIT NONE
    
    TYPE(Model_t) :: Model
    INTEGER :: n
    REAL(KIND=dp) :: t,f

    LOGICAL :: Visited = .FALSE., Found
    TYPE( Variable_t), POINTER :: PotA, PotB
    TYPE( Nodes_t), POINTER :: Nodes
    INTEGER :: m
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: xmin, x0

    SAVE Visited, PotA, PotB, Nodes, x0

    IF( .NOT. Visited ) THEN
      PotB => VariableGet( Model % Mesh % Variables,'CoilPotB' )
      PotA => VariableGet( Model % Mesh % Variables,'CoilPot' )
      Nodes => Model % Mesh % Nodes
      x0 = ListGetCReal( Model % Simulation,'Coil x0',Found )
      Visited = .TRUE.
    END IF

    Element => Model % CurrentElement
    m = GetElementNOFNodes()

    ! One could use as well max or mean, for example
    ! Consistancy is most important
    xmin = MINVAL( Nodes % x(Element % NodeIndexes) )
    IF( xmin - x0 > 0.0 ) THEN
      f = PotB % Values( PotB % Perm(n) )
    ELSE
      f = PotA % Values( PotA % Perm(n) )
    END IF

  END FUNCTION UnitedPotential
