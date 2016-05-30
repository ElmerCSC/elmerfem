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
! *  Authors: Juha Ruokolainen, Peter Råback
! *  Email:   Peter.Raback@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 16.11.2005
! *
! *****************************************************************************/
!------------------------------------------------------------------------------
!>  Compute the volume and area in 3D or area and line intgral in 2D over the 
!>  levelset function. This is better done within a dedicated solver since 
!>  it is crucial for the accuracy that the Heaviside and Delta function are
!>  computed at Gaussian integration points.
!> \ingroup Solvers
!------------------------------------------------------------------------------
   SUBROUTINE LevelSetIntegrate( Model,Solver,Timestep,TransientSimulation )
!------------------------------------------------------------------------------
     USE DefUtils
     USE SolverUtils
     USE MaterialModels
     USE Integration

     IMPLICIT NONE
!------------------------------------------------------------------------------
 
     TYPE(Model_t), TARGET :: Model
     TYPE(Solver_t) :: Solver
     REAL(KIND=dp) :: Timestep
     LOGICAL :: TransientSimulation
 
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     TYPE(Nodes_t)   :: ElementNodes
     TYPE(Element_t),POINTER :: CurrentElement
     TYPE(ValueList_t), POINTER :: Material
     TYPE(Variable_t), POINTER :: SurfSol

     INTEGER :: i,j,k,n,t,istat,bf_id
     INTEGER, POINTER :: NodeIndexes(:), SurfPerm(:)
     LOGICAL :: Visited = .FALSE., GotIt
     REAL(KIND=dp), POINTER :: Surface(:), NodalSurf(:)
     INTEGER :: body_id, dim
     REAL(KIND=dp) :: TotVolume, TotArea, Alpha, Relax, InitVolume, dSurface, dt, &
         Moment(3)
     CHARACTER(LEN=MAX_NAME_LEN) :: LevelSetVariableName
     TYPE(ValueList_t), POINTER :: Params

     SAVE ElementNodes, Visited, NodalSurf, InitVolume

!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------

     Params => GetSolverParams()

     ! The variable that should be renormalized
     LevelSetVariableName = ListGetString(Params,'Level Set Variable',GotIt) 
     IF(GotIt) THEN
       SurfSol => VariableGet( Solver % Mesh % Variables, TRIM(LevelSetVariableName) )
     ELSE  
       SurfSol => VariableGet( Solver % Mesh % Variables, 'Surface' )
     END IF
     IF(ASSOCIATED(SurfSol)) THEN
       Surface  => SurfSol % Values
       SurfPerm => SurfSol % Perm
     ELSE
       CALL Warn('LevelSetIntegrate','Surface variable does not exist')
     END IF

     IF ( ALL( SurfPerm == 0) ) THEN
       CALL Warn('LevelSetIntegrate','Nothing to compute')
       RETURN
     END IF
 
     dim = CoordinateSystemDimension()
 
!------------------------------------------------------------------------------
!    Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
     IF ( .NOT. Visited ) THEN
       N = Solver % Mesh % MaxElementNodes

       ALLOCATE( NodalSurf( N ), &
           ElementNodes % x( N ),   &
           ElementNodes % y( N ),   &
           ElementNodes % z( N ),   &
           STAT=istat )
 
       IF ( istat /= 0 ) THEN
         CALL Fatal( 'LevelSetIntegrate', 'Memory allocation error.' )
       END IF
     END IF
!------------------------------------------------------------------------------
!    Do some additional initialization, and go for it
!------------------------------------------------------------------------------

     TotVolume = 0.0d0
     TotArea = 0.0d0
     Moment = 0.0d0
     
     Alpha = ListGetConstReal(Model % Simulation,'Levelset Bandwidth',GotIt) 
     IF(.NOT. GotIt) Alpha = ListGetConstReal(Params,'Levelset Bandwidth')      
       
     CALL Info( 'LevelSetIntegrate','-------------------------------------', Level=4 )
     CALL Info( 'LevelSetIntegrate', 'Integrating over levelset function', Level=4 )
     CALL Info( 'LevelSetIntegrate','-------------------------------------', Level=4 )

     DO t=1,Solver % Mesh % NumberOfBulkElements
       
       CurrentElement => Solver % Mesh % Elements(t)
       n = CurrentElement % TYPE % NumberOfNodes
       NodeIndexes => CurrentElement % NodeIndexes
       IF( ANY(SurfPerm(NodeIndexes) == 0)) CYCLE
       
       Model % CurrentElement => CurrentElement
       body_id = CurrentElement % Bodyid    
       k = ListGetInteger( Model % Bodies( body_id ) % Values, 'Material' )
       Material => Model % Materials(k) % Values
       
!-----------------------------------------------------------------------------
!        Get element nodal coordinates
!------------------------------------------------------------------------------
       ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
       ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
       ElementNodes % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes)
         
       NodalSurf = Surface( SurfPerm(NodeIndexes) )

       CALL HeavisideIntegrate( NodalSurf, CurrentElement, n, ElementNodes, &
           Alpha, TotVolume, TotArea, Moment)      
     END DO
     
     Moment = Moment / TotVolume
          
     IF(dim == 3) THEN
       WRITE(Message,'(a,ES12.3)') 'Center 3',Moment(3)
       CALL Info( 'LevelSetIntegrate',Message, Level=4 )

       WRITE(Message,'(a,ES12.3)') 'Inside Volume',TotVolume
       CALL Info( 'LevelSetIntegrate',Message, Level=4 )
     
       WRITE(Message,'(a,ES12.3)') 'Interface Area',TotArea
       CALL Info( 'LevelSetIntegrate',Message, Level=4 )
       
       CALL ListAddConstReal(Model % Simulation,'res: LevelSet Center 3',Moment(3))
       CALL ListAddConstReal(Model % Simulation,'res: LevelSet Volume',TotVolume)
       CALL ListAddConstReal(Model % Simulation,'res: LevelSet Area',TotArea)
     ELSE       
       WRITE(Message,'(a,ES12.3)') 'Inside Area',TotVolume
       CALL Info( 'LevelSetIntegrate',Message, Level=4 )
     
       WRITE(Message,'(a,ES12.3)') 'Interface length',TotArea
       CALL Info( 'LevelSetIntegrate',Message, Level=4 )
       
       CALL ListAddConstReal(Model % Simulation,'res: LevelSet Area',TotVolume)
       CALL ListAddConstReal(Model % Simulation,'res: LevelSet Length',TotArea)       
     END IF

     WRITE(Message,'(a,ES12.3)') 'Center 2',Moment(2)
     CALL Info( 'LevelSetIntegrate',Message, Level=4 )

     WRITE(Message,'(a,ES12.3)') 'Center 1',Moment(1)
     CALL Info( 'LevelSetIntegrate',Message, Level=4 )
     
     CALL ListAddConstReal(Model % Simulation,'res: LevelSet Center 2',Moment(2))
     CALL ListAddConstReal(Model % Simulation,'res: LevelSet Center 1',Moment(1))

     IF ( ListGetLogical(Params,'Conserve Volume',GotIt) ) THEN
       IF(.NOT. Visited) THEN
         InitVolume = ListGetConstReal(Params,'Initial Volume',GotIt)
         IF(.NOT. GotIt) InitVolume = TotVolume
       END IF

       Relax = ListGetConstReal(Params,'Conserve Volume Relaxation',GotIt)
       IF(.NOT. GotIt) Relax = 1.0d0
      
       dSurface = Relax * (InitVolume - TotVolume) / TotArea
       IF(ABS(dSurface) > AEPS) THEN
          Surface = Surface + dSurface
       END IF

       WRITE(Message,'(a,ES12.3)') 'Levelset correction',dSurface
       CALL Info( 'LevelSetIntegrate',Message, Level=4 )
       CALL ListAddConstReal(Model % Simulation,'res: Levelset Correction',dSurface)
     END IF

     Visited = .TRUE.


!------------------------------------------------------------------------------

 CONTAINS

!------------------------------------------------------------------------------
   SUBROUTINE HeavisideIntegrate( Surf,Element,n,Nodes,Alpha,Volume,Area,Moment)
!------------------------------------------------------------------------------

     REAL(KIND=dp), DIMENSION(:)   :: Surf
     INTEGER :: n
     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t), POINTER :: Element
     REAL(KIND=dp) :: Alpha, Volume, Area, Moment(3)

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------

     REAL(KIND=dp) :: Basis(n)
     REAL(KIND=dp) :: dBasisdx(n,3),detJ
     REAL(KIND=dp) :: Val,Grad(3),Velo(3),NormalVelo,GradAbs,Heavi,Delta
     INTEGER :: t,N_Integ, CoordinateSystem
     REAL(KIND=dp) :: s,u,v,w,x,y,z
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ
     LOGICAL :: stat

!------------------------------------------------------------------------------

     dim = CoordinateSystemDimension()
     CoordinateSystem = CurrentCoordinateSystem()

!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
     IntegStuff = GaussPoints( element )

     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------
!    Now we start integrating
!------------------------------------------------------------------------------
     DO t=1,N_Integ

       u = U_Integ(t)
       v = V_Integ(t)
       w = W_Integ(t)

!------------------------------------------------------------------------------
!      Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo( Element,Nodes,u,v,w,detJ, Basis,dBasisdx)

       s = detJ * S_Integ(t)

       x = SUM( Nodes % x(1:n) * Basis(1:n) )
       y = SUM( Nodes % y(1:n) * Basis(1:n) )
       IF(dim == 3) z = SUM( Nodes % z(1:n) * Basis(1:n) )

!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
       IF ( CoordinateSystem == AxisSymmetric .OR. CoordinateSystem == CylindricSymmetric ) THEN
         s = s * x * 2.0d0 * PI
       END IF

!------------------------------------------------------------------------------
       Val = SUM( Basis(1:n) * Surf(1:n) )

       IF( Val < -Alpha) THEN
         Heavi = 0.0d0
       ELSE IF(Val > Alpha) THEN
         Heavi = 1.0d0
       ELSE
         Heavi = (1.0d0 + SIN( (Val/Alpha) * (PI/2) ) ) / 2.0d0
         Delta = (1.0d0 + COS( (Val/Alpha) * PI ) ) / (2.0d0 * Alpha)
         
         DO i=1,dim
           Grad(i) = SUM( dBasisdx(1:n,i) * Surf(1:n) )
         END DO
         GradAbs = SQRT( SUM( Grad(1:dim) * Grad(1:dim) ) )          
         
         Area = Area + s * Delta * GradAbs
       END IF

       Volume = Volume + s * Heavi
       Moment(1) = Moment(1) + s * Heavi * x
       Moment(2) = Moment(2) + s * Heavi * y
       IF(dim == 3) Moment(3) = Moment(3) + s * Heavi * z       

     END DO
     
!------------------------------------------------------------------------------
   END SUBROUTINE HeavisideIntegrate
!------------------------------------------------------------------------------

 END SUBROUTINE LevelSetIntegrate
!------------------------------------------------------------------------------
