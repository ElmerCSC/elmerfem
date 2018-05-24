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
!>  Renormalizes the levelset function using straight-forward geometric search
!>  Also includes an option to do the covection at the same time as an alternavtive
!>  for using a separate solver for the advection.
!> \ingroup Solvers
!------------------------------------------------------------------------------
   SUBROUTINE LevelSetDistance( Model,Solver,Timestep,TransientSimulation )
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
     TYPE(Element_t),POINTER :: CurrentElement, Element
     TYPE(ValueList_t), POINTER :: Material
     TYPE(Variable_t), POINTER :: SurfSol, DistanceSol
 
     INTEGER :: i,j,k,l,n,t,iter,istat,body_id,mat_id,bf_id,&
          CoordinateSystem, TimesVisited = 0
     REAL(KIND=dp) :: Norm,x,y,s,a,b,c,d,x0,x1,y0,y1
     INTEGER, POINTER :: NodeIndexes(:)
     LOGICAL :: GotIt, Convect, ExtractAllocated = .FALSE., DistanceAllocated=.FALSE.
     INTEGER, POINTER :: SurfPerm(:)
     REAL(KIND=dp), POINTER :: Surface(:),Distance(:), Surf(:)
     REAL(KIND=dp), ALLOCATABLE :: ZeroNodes(:,:,:), Direction(:)
     INTEGER, POINTER :: DistancePerm(:)
     INTEGER :: ZeroLevels, ReinitializeInterval, ExtractInterval
     LOGICAL :: Reinitialize, Extrct
     REAL(KIND=dp) :: Relax, dt, r, NarrowBand, DsMax
     REAL(KIND=dp), ALLOCATABLE :: ElemVelo(:,:), SurfaceFlux(:)
#ifdef USE_ISO_C_BINDINGS
     REAL(KIND=dp) :: at,totat,st,totst
#else
     REAL(KIND=dp) :: at,totat,st,totst,CPUTime
#endif
     CHARACTER(LEN=MAX_NAME_LEN) :: LevelSetVariableName

     SAVE ElementNodes, ElemVelo, Direction, ZeroNodes, TimesVisited, &
         Distance, DistancePerm, ExtractAllocated, DistanceAllocated

!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------

     TimesVisited = TimesVisited + 1
     ReinitializeInterval = ListGetInteger(Solver % Values,&
          'Reinitialize Interval',GotIt) 
     IF(.NOT. GotIt) ReinitializeInterval = 1

     ExtractInterval = ListGetInteger(Solver % Values,&
         'Extract Interval',GotIt) 
     IF(.NOT. GotIt) ExtractInterval = ReinitializeInterval
     
     IF( ReinitializeInterval == 0) THEN
       Reinitialize = .FALSE.
     ELSE       
       Reinitialize = ( MOD(TimesVisited, ReinitializeInterval) == 0 )
     END IF

     IF( ExtractInterval == 0) THEN
       Extrct = Reinitialize
     ELSE
       Extrct = Reinitialize .OR. ( MOD(TimesVisited, ExtractInterval) == 0 )
     END IF
     
     IF(.NOT. Extrct) THEN
       CALL Info( 'LevelSetDistance','--------------------------------------', Level=4 )
       CALL Info( 'LevelSetDistance','Doing nothing this time', Level=4 )
       CALL Info( 'LevelSetDistance','--------------------------------------', Level=4 )          
       RETURN
     END IF

     ! The variable that should be reinitialized
     LevelSetVariableName = ListGetString(Solver % Values,'LevelSet Variable',GotIt) 
     IF(.NOT. GotIT) LevelSetVariableName = 'Surface'
     SurfSol => VariableGet( Solver % Mesh % Variables, TRIM(LevelSetVariableName) )
     IF(ASSOCIATED(SurfSol)) THEN
       Surface  => SurfSol % Values
       SurfPerm => SurfSol % Perm
     ELSE
       CALL Warn('LevelSetDistance','SurfSol does not exist: '//TRIM(LevelSetVariableName))
       RETURN
     END IF

     CoordinateSystem = CurrentCoordinateSystem() 
     Convect = ListGetLogical(Solver % Values,'Levelset Convect',GotIt)
     NarrowBand = ListGetConstReal(Solver % Values,'Narrow Band',GotIt)
     IF(.NOT. GotIt) NarrowBand = HUGE(NarrowBand)
     dsMax = 0.0d0
     dt = Timestep
     
 
!------------------------------------------------------------------------------
!    Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
     IF ( .NOT. ExtractAllocated ) THEN
       N = Solver % Mesh % MaxElementNodes
       ALLOCATE( ElementNodes % x( N ), ElementNodes % y( N ), ElementNodes % z( N ),   &
           ElemVelo( 2, N), ZeroNodes(Solver % Mesh % NumberOfBulkElements,2,2), &
           STAT=istat )
       IF ( istat /= 0 ) THEN
         CALL Fatal( 'LevelSetDistance', 'Memory allocation error 1.' )
       END IF
       IF( Convect ) THEN
         ALLOCATE( Direction( Solver % Mesh % NumberOfBulkElements), STAT=istat )
       END IF
       IF ( istat /= 0 ) THEN
         CALL Fatal( 'RerormalizeSolver', 'Memory allocation error 2.' )
       END IF
       ExtractAllocated = .TRUE.
     END IF

!------------------------------------------------------------------------------
!    Extract the zero levelset 
!------------------------------------------------------------------------------

     CALL Info( 'LevelSetDistance','--------------------------------------', Level=4 )
     CALL Info( 'LevelSetDistance','Extracting the zero levelset', Level=4 )
     CALL Info( 'LevelSetDistance','--------------------------------------', Level=4 )

     st = CPUTime()

     CALL ExtractZeroLevel()

     st = CPUTIme()-st
     WRITE(Message,'(a,F8.2)') 'Zero level extracted in time (s):',st
     CALL Info( 'LevelSetDistance',Message, Level=4 )

     IF( ZeroLevels == 0) THEN
       CALL Warn('LevelSetDistance','The does not seem to be a zero level-set present, exiting...')
       RETURN
     END IF

     IF(.NOT. Reinitialize) THEN
       CALL Info('LevelSetDistance','Exiting without reinitialization')
       RETURN
     END IF

     CALL Info( 'LevelSetDistance','--------------------------------------', Level=4 )
     CALL Info( 'LevelSetDistance','Computing the signed distance function', Level=4 )
     CALL Info( 'LevelSetDistance','--------------------------------------', Level=4 )

!------------------------------------------------------------------------------
!    Allocate some permanent storage for computing the signed distance
!------------------------------------------------------------------------------
     IF ( .NOT. DistanceAllocated ) THEN
       
       ! The variable for computing the distance
       DistanceSol => Solver % Variable
       IF(ASSOCIATED(DistanceSol)) THEN
         DistancePerm => DistanceSol % Perm
         Distance => DistanceSol % Values
       ELSE
         ALLOCATE(Distance (SIZE(Surface)),STAT=istat)
         IF ( istat /= 0 ) THEN
           CALL Fatal( 'LevelSetDistance', 'Memory allocation error 1.' )
         END IF
         DistancePerm => SurfPerm
       END IF
       DistanceAllocated = .TRUE.
     END IF

!------------------------------------------------------------------------------
!    Compute the signed distance
!------------------------------------------------------------------------------
     st = CPUTIme()
     IF(Convect) THEN
       DO i=1,Solver % Mesh % NumberOfNodes
         Distance(DistancePerm(i)) =  ComputeDistanceWithDirection( Solver % Mesh % Nodes % x(i), &
             Solver % Mesh % Nodes % y(i), 0.0d0, Surface(i) )
       END DO
     ELSE
       DO i=1,Solver % Mesh % NumberOfNodes
         Distance(DistancePerm(i)) =  ComputeDistance( Solver % Mesh % Nodes % x(i), &
             Solver % Mesh % Nodes % y(i), 0.0d0 )
       END DO
       WHERE( Surface < 0 ) Distance = -Distance
     END IF

     n = Solver % Mesh % NumberOfNodes
     Solver % Variable % Norm = SQRT( SUM(Distance**2)/n )

!------------------------------------------------------------------------------
!    Apply the reinitialization to the primary levelset field
!------------------------------------------------------------------------------

     IF( ListGetLogical(Solver % Values,'Reinitialize Passive',GotIt) ) THEN
       CALL Info('LevelSetDistance','Reinitialization not applied to Levelset function')
     ELSE
       ! Update also the previous timesteps so that the differentials remain
       ! unchanged. Otherwise spurious effects are introduced. 
       IF(ASSOCIATED(SurfSol % PrevValues)) THEN        
         j = MIN(2, SIZE(SurfSol % PrevValues,2) )
         IF( ReinitializeInterval > j) THEN                 
           DO i=1,j
             SurfSol % PrevValues(:,i) = SurfSol % PrevValues(:,i) + Distance - Surface
           END DO
         END IF
       END IF
       Surface = Distance
     END IF

     st = CPUTIme()-st
     WRITE(Message,'(a,F8.2)') 'Reinitialization done in time (s):',st
     CALL Info( 'LevelSetNormalize',Message, Level=4 )
 
     IF(Convect) THEN
       WRITE(Message,'(a,ES12.3)') 'Maximum Levelset Change',dsmax
       CALL Info( 'LevelSetDistance',Message, Level=4 )     
       CALL ListAddConstReal(Model % Simulation,'res: LevelSet Max Change',dsmax)
     END IF

!------------------------------------------------------------------------------

CONTAINS


!------------------------------------------------------------------------------
!> Extract the zero levelset as defined by the levelset function.
!------------------------------------------------------------------------------
  SUBROUTINE ExtractZeroLevel()
!------------------------------------------------------------------------------

    INTEGER :: i,j,k,l,m,n,div,onetwo,corners, &
        TriangleIndexes(3), mat_id, body_id, LocalInd(3)
    TYPE(Variable_t), POINTER :: Var
    INTEGER, POINTER :: NodeIndexes(:)
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: x0,y0,x1,y1,t,nx(3),ny(3),nz(3),srf(3),w0(3),w1(3),fval, &
        r1x, r1y, r2x, r2y, aid, Maxsrf, ds1, ds2
    LOGICAL :: FileCreated = .FALSE., FileAppend, GotIt, FileSave, Found, FileNumber
    TYPE(ValueList_t),POINTER :: Material
    CHARACTER(LEN=MAX_NAME_LEN) :: Filename, Filename2
    INTEGER :: VisitedTimes = 0, NumberOfFields=0
    TYPE(ValueList_t), POINTER :: Params

    SAVE VisitedTimes, FileAppend, NumberOfFields
!------------------------------------------------------------------------------

    VisitedTimes = VisitedTimes + 1
    Params => GetSolverParams()
    Filename = ListGetString(Params,'Filename',FileSave )

    IF(FileSave) THEN         
      FileNumber = ListGetLogical(Params,'Filename Numbering',GotIt)     
      FileAppend = ListGetLogical(Params,'File Append',GotIt)

      IF( FileNumber ) THEN
        WRITE( Filename,'(A,I0)') TRIM(Filename),VisitedTimes
        OPEN (10,FILE=Filename)
      ELSE IF(FileAppend .AND. VisitedTimes > 1) THEN 
        OPEN (10, FILE=Filename, POSITION='APPEND')
      ELSE 
        OPEN (10,FILE=Filename)
      END IF      
    END IF
    
    Surface => SurfSol % Values
    SurfPerm => SurfSol % Perm
    
    ZeroLevels = 0
    DO i=1,Solver % Mesh % NumberOfBulkElements
      
      Element => Solver % Mesh % Elements(i)
      n = Element % TYPE % NumberOfNodes
      NodeIndexes => Element % NodeIndexes
      
      IF ( ALL( Surface(SurfPerm(NodeIndexes)) < 0) .OR. &
          ALL( Surface(SurfPerm(NodeIndexes)) > 0) ) CYCLE
      
      corners = Element % TYPE % ElementCode / 100 
      IF(corners < 3 .OR. corners > 4) THEN
        CALL Warn('ExtractZeroLevel','Implemented only for triangles and quads')
      END IF

      IF( Convect ) THEN
        Model % CurrentElement => Element

        body_id = Element % BodyId
        mat_id = ListGetInteger( Model % Bodies(body_id) % Values, 'Material', &
            minv=1, maxv=Model % NumberOfMaterials )
        Material => Model % Materials(mat_id) % Values

        ElemVelo(1,1:n) = ListGetReal( Material,'Levelset Velocity 1',n,NodeIndexes,GotIt)
        ElemVelo(2,1:n) = ListGetReal( Material,'Levelset Velocity 2',n,NodeIndexes,GotIt)
      END IF


      DO div = 1,corners-2

        SELECT CASE (corners) 
        CASE (3)
          LocalInd(1) = 1
          LocalInd(2) = 2
          LocalInd(3) = 3
          
        CASE(4)               
          IF(div == 1) THEN
            LocalInd(1) = 1 
            LocalInd(2) = 2
            LocalInd(3) = 4
          ELSE
            LocalInd(1) = 2
            LocalInd(2) = 3
            LocalInd(3) = 4
          END IF

        END SELECT
        
        TriangleIndexes = NodeIndexes(LocalInd)
        srf = Surface(SurfPerm(TriangleIndexes))                       
        IF ( ALL(srf < 0) .OR. ALL( srf > 0) ) CYCLE

        nx = Solver % Mesh % Nodes % x(TriangleIndexes)
        ny = Solver % Mesh % Nodes % y(TriangleIndexes)
        nz = Solver % Mesh % Nodes % z(TriangleIndexes)

        CALL TriangleIsoLineWeights( nx,ny,nz,srf,w0,w1,found)
        IF( .NOT. Found) CYCLE

        x0 = SUM(w0 * nx)
        y0 = SUM(w0 * ny)
        x1 = SUM(w1 * nx)
        y1 = SUM(w1 * ny)

        r1x = x1 - x0
        r1y = y1 - y0 

        ds1 = SQRT( r1x*r1x + r1y*r1y)
        IF(ds1 < AEPS) CYCLE

        ZeroLevels = ZeroLevels + 1

        IF(Convect) THEN

          ! Find the value that differs most from the zero levelset
          j = 1
          MaxSrf = srf(1)
          DO k=2,3
            IF(ABS(srf(k)) > ABS(Maxsrf)) THEN
              Maxsrf = srf(k)
              j = k
            END IF
          END DO
          
          r2x = nx(j) - x0
          r2y = ny(j) - y0
          
          aid = r1x * r2y - r2x * r1y
          ds2 = SQRT( r2x*r2x + r2y*r2y)
          
          Direction( ZeroLevels ) = aid / (ds1 * ds2)
          IF( Maxsrf < 0.0) THEN
            Direction( ZeroLevels ) = -Direction( ZeroLevels )
          END IF
          
          r1x = SUM(w0 * ElemVelo(1,LocalInd)) * dt
          r1y = SUM(w0 * ElemVelo(2,LocalInd)) * dt
          r2x = SUM(w1 * ElemVelo(1,LocalInd)) * dt
          r2y = SUM(w1 * ElemVelo(2,LocalInd)) * dt

          ds1 = SQRT( r1x*r1x + r1y*r1y)
          ds2 = SQRT( r2x*r2x + r2y*r2y)
          dsmax = MAX(dsmax, MAX(ds1, ds2) )
            
          ZeroNodes(ZeroLevels,1,1) = x0 + r1x
          ZeroNodes(ZeroLevels,1,2) = y0 + r1y
          ZeroNodes(ZeroLevels,2,1) = x1 + r2x
          ZeroNodes(ZeroLevels,2,2) = y1 + r2y
        ELSE
          ZeroNodes(ZeroLevels,1,1) = x0
          ZeroNodes(ZeroLevels,1,2) = y0
          ZeroNodes(ZeroLevels,2,1) = x1
          ZeroNodes(ZeroLevels,2,2) = y1
        END IF


        IF(FileSave) THEN
          
          DO onetwo = 1,2
            
            IF( FileAppend ) THEN
              WRITE(10,'(I4)',ADVANCE='NO') Solver % DoneTime
            END IF            

	    m = 0
            Var => Model % Variables
            DO WHILE( ASSOCIATED( Var ) )
              
              IF ( .NOT. Var % Output .OR. SIZE(Var % Values) == 1 .OR. (Var % DOFs /= 1) ) THEN
                Var => Var % Next        
                CYCLE
              END IF
	      m = m + 1              

              fval = 0.0d0
              DO k=1,3
                l = TriangleIndexes(k)
                IF ( ASSOCIATED(Var % Perm) ) l = Var % Perm(l)
                IF(l > 0) THEN
                  IF(onetwo == 1) THEN
                    fval = fval + w0(k) * (Var % Values(l))
                  ELSE
                    fval = fval + w1(k) * (Var % Values(l))
                  END IF
                END IF
              END DO
              
              WRITE(10,'(ES20.11E3)',ADVANCE='NO') fval
              Var => Var % Next          
            END DO
            WRITE(10,'(A)') ' '
            
          END DO
	END IF        

      END DO

    END DO ! of elements

    


    IF(FileSave) THEN
      CLOSE(10)

      IF( m /= NumberOfFields ) THEN
	IF( NumberOfFields > 0 ) THEN
  	  CALL Warn('ExtractZeroLevel','Mismacth in number of fields')
          PRINT *,j,' vs. ',NumberOfFields
	END IF
	NumberOfFields = m

        OPEN (10, FILE=TRIM(Filename)//TRIM(".names") )
        WRITE(10,'(A,A)') 'Variables in file: ',TRIM(Filename)
        j = 1
        WRITE(10,'(I3,": ",A)') j,'timestep'
        
        Var => Model % Variables
        DO WHILE( ASSOCIATED( Var ) )          
          IF ( .NOT. Var % Output .OR. SIZE(Var % Values) == 1 .OR. (Var % DOFs /= 1) ) THEN
            Var => Var % Next        
            CYCLE 
          END IF          
          j = j + 1
          WRITE(10,'(I3,": ",A)') j,TRIM(Var % Name)
          Var => Var % Next          
        END DO
        CLOSE(10)
      END IF  
    END  IF

!------------------------------------------------------------------------------
   END SUBROUTINE ExtractZeroLevel
!------------------------------------------------------------------------------
 

!------------------------------------------------------------------------------
!> This subroutine extracts the zero line of one triangular element. 
!------------------------------------------------------------------------------
   SUBROUTINE TriangleIsoLineWeights( NX,NY,NZ,S,w0,w1,Found )
!------------------------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND=dp) :: NX(:),NY(:),NZ(:),S(:),w0(3),w1(3)
      LOGICAL :: Found
      REAL(KIND=dp) :: t
!------------------------------------------------------------------------------

      Found = .TRUE.

      w0 = 0.0d0
      w1 = 0.0d0

      IF ( ABS(S(1)) < AEPS .AND. ABS(S(2)) < AEPS ) THEN
        w0(1) = 1.0d0
        w1(2) = 1.0d0
      ELSE IF ( ABS(S(1)) < AEPS .AND. ABS(S(3)) < AEPS ) THEN
        w0(1) = 1.0d0
        w1(3) = 1.0d0
      ELSE IF ( ABS(S(2)) < AEPS .AND. ABS(S(3)) < AEPS ) THEN
        w0(2) = 1.0d0
        w1(3) = 1.0d0
      ELSE IF ( ALL(S <= 0) .OR. ALL( S >= 0) ) THEN
        Found = .FALSE.
      ELSE
        IF ( S(1) >= 0 .AND. S(2) >= 0 .OR. &
            S(1) <= 0 .AND. S(2) <= 0 ) THEN          
          t = -S(1) / ( S(3) - S(1) )
          w0(3) = t
          w0(1) = 1-t
          t = -S(2) / ( S(3) - S(2) )
          w1(3) = t
          w1(2) = 1-t
        ELSE IF ( S(1) >= 0 .AND. S(3) >= 0 .OR. &
            S(1) <= 0 .AND. S(3) <= 0 ) THEN
          t = -S(1) / ( S(2) - S(1) )
          w0(2) = t
          w0(1) = 1-t
          t = -S(3) / ( S(2) - S(3) )
          w1(2) = t
          w1(3) = 1-t
          
        ELSE IF ( S(2) >= 0 .AND. S(3) >= 0 .OR. &
            S(2) <= 0 .AND. S(3) <= 0 ) THEN          
          t = -S(2) / ( S(1) - S(2) )
          w0(1) = t
          w0(2) = 1-t
          t = -S(3) / ( S(1) - S(3) )
          w1(1) = t
          w1(3) = 1-t
        ELSE 
          PRINT *,'TriangleIsoLineWeights: this should not occur'
          PRINT *,s(1),s(2),s(3)
          STOP
        END IF
      END IF
!------------------------------------------------------------------------------
    END SUBROUTINE TriangleIsoLineWeights
!------------------------------------------------------------------------------

 
!------------------------------------------------------------------------------
!> Computes the distance from the given zero levelset given by ZeroNodes. 
!------------------------------------------------------------------------------
   FUNCTION ComputeDistance(xp,yp,zp) RESULT(dist)
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: xp,yp,zp,dist
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: x0,y0,x1,y1,a,b,c,d,s
     INTEGER :: i,j,k,n
!------------------------------------------------------------------------------
     dist = HUGE(dist)
     DO i=1,ZeroLevels
       x0 = ZeroNodes(i,1,1)
       y0 = ZeroNodes(i,1,2)
       
       x1 = ZeroNodes(i,2,1)
       y1 = ZeroNodes(i,2,2)
       
       a = xp - x0
       b = x0 - x1
       d = y0 - y1
       c = yp - y0
       s = b**2 + d**2
       
       x = x0
       y = y0
       IF ( s > 10*AEPS ) THEN
         s = MIN( MAX( -(a*b + c*d) / s, 0.0d0), 1.0d0 )
         x = (1-s) * x0 + s * x1
         y = (1-s) * y0 + s * y1
       END IF
       
       dist = MIN( dist, SQRT( (xp - x)**2 + (yp - y)**2 ) )
     END DO
!------------------------------------------------------------------------------
   END FUNCTION ComputeDistance
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Computes the signed distance from the given zero levelset given by ZeroNodes. 
!------------------------------------------------------------------------------
   FUNCTION ComputeDistanceWithDirection(xp,yp,zp,prevdist) RESULT(mindist)
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: xp,yp,zp,prevdist,mindist
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: x0,y0,x1,y1,a,b,c,d,s,dist,r1x,r1y,r2x,r2y,angle,angle0
     INTEGER :: i,j,k,n
!------------------------------------------------------------------------------

     IF(prevdist - dsmax > narrowband) THEN
       mindist = prevdist - dsmax
       RETURN
     ELSE IF(prevdist + dsmax < -narrowband) THEN
       mindist = prevdist + dsmax
       RETURN
     END IF

     mindist = HUGE(mindist)

     DO i=1,ZeroLevels
       x0 = ZeroNodes(i,1,1)
       y0 = ZeroNodes(i,1,2)
       
       x1 = ZeroNodes(i,2,1)
       y1 = ZeroNodes(i,2,2)
       
       a = xp - x0
       b = x0 - x1
       d = y0 - y1
       c = yp - y0
       s = b**2 + d**2
       
       x = x0
       y = y0
       IF ( s > 10*AEPS ) THEN
         s = MIN( MAX( -(a*b + c*d) / s, 0.0d0), 1.0d0 )
         x = (1-s) * x0 + s * x1
         y = (1-s) * y0 + s * y1
       END IF
       
       dist = SQRT( (xp - x)**2 + (yp - y)**2 ) 
       

       IF(dist <= (ABS(mindist) + AEPS) ) THEN
         
         r1x = x1 - x0
         r1y = y1 - y0
         r2x = xp - x0
         r2y = yp - y0
         
         angle = r1x * r2y - r2x * r1y
         
         ! Favor parents with clear angles
         IF( dist < (ABS(mindist) - AEPS) .OR. (ABS(angle) > ABS(angle0)) ) THEN
           IF(Direction(i) * angle < 0.0) THEN
             mindist = -dist
           ELSE
             mindist = dist
           END IF
           angle0 = angle
         END IF
           
       END IF

      END DO
!------------------------------------------------------------------------------
    END FUNCTION ComputeDistanceWithDirection
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
 END SUBROUTINE LevelSetDistance
!------------------------------------------------------------------------------

