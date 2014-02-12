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
! * CavitySolve                                           
! * To be compiled with two functions fbed = bed(x) and dbed = d bed(x) / dx
! *    move the nodes of a Free surface by -dt(u.n)n 
! *    n = outside normal to the free surface      
! *    u = velocity from either NS solver or AIFlow solver      
! *    dt =  time increment 
! *
! ******************************************************************************
! *
! *                    Author:  Olivier Gagliardini
! *
! *                       Date: 27 Sep 2005
! *
! *                Modified by: 
! *
! *
!/******************************************************************************
! *
! *       Modified by: 
! *
! *       Date of modification: 
! *
! *****************************************************************************/
!---------------------------------------------------------------------
!---------------------------------------------------------------------
       SUBROUTINE CavitySolve( Model, Solver, dt, TransientSimulation )

       USE types
       USE DefUtils
       USE CoordinateSystems
       USE SolverUtils
       USE ElementDescription
!-----------------------------------------------------------
       IMPLICIT NONE
!-----------------------------------------------------------
       TYPE(Solver_t), TARGET :: Solver
       TYPE(Model_t) :: Model
       REAL(KIND=dp) :: dt
       LOGICAL :: TransientSimulation
!-----------------------------------------------------------
       TYPE(Variable_t), POINTER :: moveFSSolution
       REAL(KIND=dp), POINTER :: moveFSvector(:)
       INTEGER, POINTER :: moveFSPerm(:), NodeIndexes(:)
       TYPE(Element_t), POINTER ::  BCElement
       TYPE(Nodes_t) :: Nodes
       TYPE(Variable_t), POINTER :: TimeVar
       TYPE(Variable_t), POINTER :: DevStressVar, & 
              FlowVariable, NormalVariable, MHVariable
       TYPE(ValueList_t), POINTER :: BC, SolverParams
       REAL(KIND=dp), POINTER ::  MHPrev(:,:), DSValues(:), & 
              FlowValues(:), NormalValues(:)
       INTEGER, POINTER :: FlowPerm(:), NormalPerm(:), MHPerm(:), DSPerm(:)
       CHARACTER(LEN=MAX_NAME_LEN)  ::  FlowSolverName, BCTypeName, & 
                                        SaveFileName
       INTEGER :: DIM, n, i, j, k, p, istat, Nn, nbed2, &
                                 nbed1,nbed10, NMAX, Nn2, Node
       REAL(KIND=dp) :: x, y, u, v, Velocity(3),lbed10, lbed20, &
                     LastTime= 0.0, xEnd, xFirst, xf, xd, xbed2, xbed1, &
                     dMULN, dMUFN, mu1, mu2, Large, Normal(3)
       REAL(KIND=dp) ::  dt1, dt2, dt3, tanbed(2), dx, TrueTime
       REAL(KIND=dp) :: aa, bb, cc, dd, s
       REAL(KIND=dp) :: EPS, pwater, ds, Nvelo, sc           
       REAL(KIND=dp) :: taub, ub, fx, fy, pice, deltaFS, MaxdeltaFS, Bu, Bv   
       REAL(KIND=dp) :: Sxx(3),Syy(3),Sxy(3),Velo(3,3), LocalNormal(3,3), &
                         lc, ltot, xa, xb, slope
       REAL(KIND=dp), ALLOCATABLE ::  xn(:), yn(:), xn0(:), &
               yn0(:), xnew(:), ynew(:), y2new(:), intersept(:),  &
               dMU(:), x0(:), y0(:), y1new(:), y3new(:), Nvector(:)
       INTEGER, ALLOCATABLE :: ordre(:), ind(:), BCType(:)
        
       LOGICAL :: AllocationsDone = .FALSE. , FirstTime = .TRUE., &
                 FirstTime2 = .TRUE., SaveFile = .FALSE.
       LOGICAL :: GotIt, succes, FirstNodeMove
       
       COMMON /interbed/aa,bb,cc,dd

       INTERFACE 
         SUBROUTINE indexx(arr,ind) 
         USE types
         USE DefUtils
         IMPLICIT NONE 
         REAL(KIND=dp) :: arr(:) 
         INTEGER :: ind(:) 
         END SUBROUTINE indexx
         SUBROUTINE spline(x,y,yp1,ypn,y2) 
         USE types
         USE DefUtils
         IMPLICIT NONE 
         REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: x,y 
         REAL(KIND=dp), INTENT(IN) :: yp1,ypn 
         REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: y2 
         END SUBROUTINE spline
         FUNCTION splint(xa,ya,y2a,x) 
         USE types
         USE DefUtils
         IMPLICIT NONE 
         REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: xa,ya,y2a 
         REAL(KIND=dp), INTENT(IN) :: x 
         REAL(KIND=dp) :: splint
         END FUNCTION splint
         FUNCTION zbrent(func,x1,x2,tol)
         USE types
         USE DefUtils
         IMPLICIT NONE 
         REAL(KIND=dp), INTENT(IN) :: x1,x2,tol 
         REAL(KIND=dp) :: zbrent 
         INTERFACE 
         FUNCTION func(x) 
         USE types
         USE DefUtils
         IMPLICIT NONE 
         REAL(KIND=dp), INTENT(IN) :: x 
         REAL(KIND=dp) :: func 
         END FUNCTION func 
         END INTERFACE 
         END FUNCTION zbrent
         FUNCTION intbed(x)
         USE types
         USE DefUtils
         REAL(KIND=dp) :: intbed,x
         END FUNCTION intbed
         FUNCTION intbed1(x)
         USE types
         USE DefUtils
         REAL(KIND=dp) :: intbed1,x
         END FUNCTION intbed1
         FUNCTION intbed2(x)
         USE types
         USE DefUtils
         REAL(KIND=dp) :: intbed2,x
         END FUNCTION intbed2
         FUNCTION curvi(x1,x2,func,N)
         USE types
         USE DefUtils
         IMPLICIT NONE 
         REAL(KIND=dp), INTENT(IN) :: x1,x2
         REAL(KIND=dp) :: curvi  
         INTEGER :: N
         INTERFACE 
         FUNCTION func(x) 
         USE types
         USE DefUtils
         IMPLICIT NONE 
         REAL(KIND=dp), INTENT(IN) :: x 
         REAL(KIND=dp) :: func 
         END FUNCTION func 
         END INTERFACE 
         END FUNCTION curvi 
         FUNCTION fbed(x)
         USE types
         USE DefUtils
         REAL(KIND=dp) :: fbed,x
         END FUNCTION fbed
         FUNCTION dbed(x)
         USE types
         USE DefUtils
         REAL(KIND=dp) :: dbed,x
         END FUNCTION dbed
       END INTERFACE

       Save DIM, Nodes,  FlowSolverName,  &
                 LastTime, xEnd, xFirst, xd, xf, Nn, &
                 FirstNodeMove, EPS, TrueTime, intersept, &
                 Large, Nn2, dMU, x0, y0, SaveFile, SaveFileName
       Save BCType, ordre, xn, xn0, yn, yn0, dMULN, dMUFN, & 
             lbed10, lbed20, Nvector



       moveFSSolution => VariableGet( Solver % Mesh % Variables, &
&                                          'moveFS' ) 
       IF ( ASSOCIATED( moveFSSolution ) ) THEN 
          moveFSvector => moveFSSolution % Values
          moveFSPerm => moveFSSolution % Perm
       ELSE
          PRINT *,'FATAL: Unable to set pointer to the current solution'
          STOP
       END IF
!
! Read data from Solver
! 
! IF FIRSTIME ------------------------------------------------------
      IF (FirstTime) THEN
         FirstTime=.FALSE.
         TrueTime = 0.0_dp 
         DIM = CoordinateSystemDimension()
         SolverParams => GetSolverParams()

         FlowSolverName = GetString( SolverParams, &
&                                    'Flow Solver Name', GotIt )    
         IF (.NOT.Gotit) THEN
           CALL INFO('Cavity Solver','No Flow Solver Name given,&
&                   default Value = Flow Solution    ', Level=4)
            FlowSolverName = 'Flow Solution'
         END IF

         FirstNodeMove = GetLogical( SolverParams, &
&                                    'First Node Move', GotIt )    
         IF (.NOT.Gotit) THEN
           CALL INFO('Cavity Solver','No First Node Move given,&
&                   default Value = False    ', Level=4)
            FirstNodeMove = .False. 
         END IF

         SaveFileName = GetString( SolverParams, &
&                                    'Save File Name', SaveFile )    

        IF (SaveFile) THEN
        OPEN(22,File=SaveFileName)
        WRITE(22,'("# 1  : Time")')
        WRITE(22,'("# 2  : TrueTime")')
        WRITE(22,'("# 3  : dt")')
        WRITE(22,'("# 4  : Tau_b")')
        WRITE(22,'("# 5  : p_ice")')
        WRITE(22,'("# 6  : p_water")')
        WRITE(22,'("# 7  : U_b")')
        WRITE(22,'("# 8  : Tau_b / N")')
        WRITE(22,'("# 9  : U_b / N")')
        WRITE(22,'("# 10 : f_x")')
        WRITE(22,'("# 11 : f_y")')
        WRITE(22,'("# 12 : U.n")')
        WRITE(22,'("# 13 : Max(U.n)")')
        WRITE(22,'("# 14 : Max(slope)")')
        WRITE(22,'("# 15 : xbed1")')
        WRITE(22,'("# 16 : xbed2")')
        WRITE(22,'("# 17 : s_cavity")')
        CLOSE(22)
        END IF

        xf = -1.0e20_dp
        xd = 1.0e20_dp
        EPS = 100.0_dp*EPSILON(xf)
        Large = 0.1_dp*HUGE(xf) 
        Nn = 1
        NMAX = 0
        DO p = 1, Model % NumberOfBoundaryElements
          BCElement => GetBoundaryElement(p) 

!         IF( .NOT. ActiveBoundaryElement(BCElement) )  CYCLE
          IF( GetElementFamily(BCElement) == 1 ) CYCLE
          
          BC => GetBC( BCElement ) 

          BCTypeName  = GetString(BC, 'Cavity Type', GotIt)

          IF (GotIt.AND.(BCTypeName=='cavity')) THEN
            n = BCElement % Type % NumberOfNodes
            IF (n>NMAX) NMAX = n
            Nn = Nn + n - 1
          END IF
          
          IF (MAXVAL(Model % Nodes % x(BCElement % NodeIndexes))>xf)  &
                 xf = MAXVAL(Model % Nodes % x(BCElement % NodeIndexes))

          IF (MINVAL(Model % Nodes % x(BCElement % NodeIndexes))<xd)  &
                 xd = MINVAL(Model % Nodes % x(BCElement % NodeIndexes))

        END DO
        Nn2 = FLOOR(Nn/2.0)
        ALLOCATE( Nodes % x( NMAX ), Nodes % y( NMAX ), Nodes % z( NMAX ))

        NMAX = Model % NumberOfNodes
        ALLOCATE( dMU(2*NMAX), x0(NMAX), y0(NMAX), NVector(2*NMAX) ) 
                 
        NMAX = Model % NumberOfBoundaryElements
        ALLOCATE( BCType(NMAX) )

        ALLOCATE( ordre(Nn), ind(Nn), xn(Nn), yn(Nn), &
                    xn0(Nn), yn0(Nn), intersept(Nn) )
        
        x0 = Model % Nodes % x
        y0 = Model % Nodes % y

! Get the nodes number of the free surface
! and the bed and the ramp

        j = 0 
        ordre = 0
        BCType = 0
        DO p = 1, Model % NumberOfBoundaryElements
          BCElement => GetBoundaryElement(p) 

!         IF( .NOT. ActiveBoundaryElement(BCElement) )  CYCLE
          IF( GetElementFamily(BCElement) == 1 ) CYCLE
          BC => GetBC( BCElement ) 
          BCTypeName  = GetString(BC, 'Cavity Type', GotIt)
            

          IF (GotIt) THEN
            n = BCElement % Type % NumberOfNodes
              IF (BCTypeName=='cavity') THEN 
                BCType(p) = 1
                DO i = 1, n
                  IF (.NOT.ANY(ordre== BCElement % NodeIndexes( i ))) THEN
                    j = j + 1
                    ordre(j) = BCElement % NodeIndexes( i )
                    xn(j) = Model % Nodes % x( BCElement % NodeIndexes( i ) )
                  END IF 
                END DO 
              ELSE IF (BCTypeName=='bed') THEN 
                BCType(p) = 2
              ELSE IF (BCTypeName=='ramp') THEN 
                BCType(p) = 3
              END IF
          END IF
        END DO

! Get ordre(i) ordered with increasing x 

        ind = 0
        CALL indexx(xn,ind)
        ordre=ordre(ind)
        xFirst = Model % Nodes % x( ordre(1) ) 
        xEnd = Model % Nodes % x( ordre(Nn) ) 
        lbed10 = curvi(xd,xFirst,fbed,100)
        lbed20 = curvi(xEnd,xf,fbed,100)
         

        DEALLOCATE(ind) 
       END IF
! END IF FIRSTIME --------------------------------------------------
!-------------------------------------------------------------------
!    Get Velocity from (Navier)-Stokes-Solver or AIFlow Solver
!------------------------------------------------------------------------
      TimeVar => VariableGet( Model % Variables, 'Time' )
      FlowVariable => VariableGet( Solver % Mesh % Variables, FlowSolverName )
      IF ( ASSOCIATED( FlowVariable ) ) THEN
       FlowPerm    => FlowVariable % Perm
       FlowValues  => FlowVariable % Values
      ELSE
      CALL Info('moveFS Solver', 'No variable for velocity &
&                                  associated.', Level=4)
      END IF

!------------------------------------------------
!    Get MeshUpdate computed from MeshUpdateSolver              
!------------------------------------------------

     MHVariable => VariableGet( Solver % Mesh % Variables, &
&                                                  "Mesh Update" )
     IF ( ASSOCIATED( MHVariable ) ) THEN
       MHPerm    => MHVariable % Perm
       MHPrev  => MHVariable % PrevValues
     ELSE
       CALL Info('moveFS Solver', 'No variable for Mesh &
&                   Update associated.', Level=4)
     END IF
      If ( FirstTime2 ) MHPrev=0.0_dp
      FirstTime2=.FALSE.
!--------------------------
!    Compute Normal Vector                                      
!--------------------------
      Nvector = 0.0_dp
      DO p = 1, Model % NumberOfBoundaryElements
        BCElement => GetBoundaryElement(p) 
        IF (GetElementFamily(BCElement) == 1 ) CYCLE
        IF (BCType(p)==0) CYCLE 
        n = BCElement % Type % NumberOfNodes
        CALL GetElementNodes( Nodes, BCElement )
        NodeIndexes => BCElement % NodeIndexes

        DO i=1,n
          Bu = BCElement % Type % NodeU(i)
          Bv = 0.0D0
          k = NodeIndexes( i )
          Normal =  NormalVector(BCElement, Nodes, Bu, Bv, .TRUE.) 
          Nvector(DIM*(k-1)+1:DIM*k) = Nvector(DIM*(k-1)+1:DIM*k) + & 
                   Normal(1:DIM)
        END DO
      END DO
! Norm 
      DO p = 1, Model % NumberOfBoundaryElements
        BCElement => GetBoundaryElement(p) 
        IF (GetElementFamily(BCElement) == 1 ) CYCLE
        IF (BCType(p)==0) CYCLE 
        n = BCElement % Type % NumberOfNodes
        NodeIndexes => BCElement % NodeIndexes
        DO i=1,n
          k = NodeIndexes( i )
          s = SQRT( SUM( Nvector(DIM*(k-1)+1:DIM*k)**2 ) )
          IF ( s /= 0.0D0 ) THEN
            Nvector(DIM*(k-1)+1:DIM*k) = Nvector(DIM*(k-1)+1:DIM*k)/s 
          END IF
        END DO
      END DO


! ---------------------------------------
! Get the MeshUpdate of the cavity  Nodes 
! ---------------------------------------

           LastTime = TimeVar % Values( 1 )
           CALL Info( 'CavitySolve', ' ', Level=4 )
           CALL Info( 'CavitySolve', '---------------',Level=4 )
           CALL Info( 'CavitySolve', ' CAVITY SOLVER ',Level=4 )
           CALL Info( 'CavitySolve', '---------------',Level=4 )
           CALL Info( 'CavitySolve', ' ', Level=4 )

          intersept = Large
          k = FlowVariable % DOFs
          deltaFS=0.0_dp
          MaxdeltaFS=0.0_dp
          DO p = 1, Nn         
            xn0(p) = x0( ordre(p) ) + &
                    MHPrev( DIM*(MHPerm( ordre(p) ) - 1 ) + 1 , 1 )
            yn0(p) = y0( ordre(p) ) + &
                    MHPrev( DIM*(MHPerm( ordre(p) ) - 1 ) + 2 , 1 )

! Get velocity  node ordre(p)   
            Velocity = 0.0_dp
            Nvelo = 0.0_dp
            DO i=1,k-1
            Velocity(i) = FlowValues( k*( FlowPerm( ordre(p) )-1 ) +i ) 
            Nvelo = Nvelo + Velocity(i)*Nvector(DIM*(ordre(p)-1)+i)
            END DO
            xn(p) = xn0(p) + dt * Velocity(1)                          
            yn(p) = yn0(p) + dt * Velocity(2)                          

            deltaFS = deltaFS + Nvelo**2 
            IF (Nvelo**2 > MaxdeltaFS) MaxdeltaFS = Nvelo**2

! test if node intersept bed
            IF ((yn(p)<=fbed(xn(p))).AND.(p/=1).AND.(p/=Nn))  THEN
              bb = xn0(p) - xn(p)
              aa = (yn0(p) - yn(p))/bb
              bb = (yn(p)*xn0(p)-xn(p)*yn0(p))/bb
              x = zbrent(intbed,xn0(p),xn(p),1.0e-08_dp)
              intersept(p) = dt * (x-xn0(p))/(xn(p)-xn0(p))
            END IF
          END DO
! should be 
          intersept(1) = 0.0_dp
          intersept(Nn) = 0.0_dp

! Find first node intersepting bed
          dt1 = dt           
          dt2 = dt            
          dt3 = dt
          nbed1 = 1
          nbed2 = Nn
          DO p = 2, Nn2
            IF (FirstNodeMove) THEN
              IF (intersept(p)<=dt1) THEN
                 dt1 = intersept(p)
                 nbed1 = p
              ENDIF
            END IF
            IF (intersept(Nn-p+1)<=dt2) THEN
               dt2 = intersept(Nn-p+1)
               nbed2 = Nn-p+1
            ENDIF
          END DO

! Find in which side it intersept first 
          IF ((nbed1/=1).AND.(nbed2/=Nn)) THEN
            IF (dt1 < dt2) THEN
! First Node intersepting is on left
              nbed2 = Nn
            ELSEIF (dt2 < dt1) THEN
! First Node intersepting is on right
              nbed1 = 1 
            END IF           
          END IF           
!
! Find dt3, xbed1 and xbed2
!
         xbed1 = xn0(1)
         nbed10=nbed1
!
! CASE 1 : No intersecting node 
!
         IF ((nbed10==1).AND.(nbed2==Nn)) THEN
           dt3 = dt
           IF (FirstNodeMove) THEN
             aa = xn(2)
             bb = yn(2)
             cc = FlowValues( k*( FlowPerm( ordre(2) )-1 )+1) 
             dd = FlowValues( k*( FlowPerm( ordre(2) )-1 )+2) 
             xa = xn0(1)
             xb = aa      
             IF (intbed2(xa)*intbed2(xb)<=0.0_dp) THEN
                xbed1 = zbrent(intbed2,xa,xb,1.0e-8_dp)
             ELSE
                xbed1 = xn(1)
             END IF
           END IF
           aa = xn(Nn-1)
           bb = yn(Nn-1) 
           cc = FlowValues( k*( FlowPerm( ordre(Nn-1) )-1 )+1) 
           dd = FlowValues( k*( FlowPerm( ordre(Nn-1) )-1 )+2) 
           xa = aa     
           xb = xn(Nn)
           IF (intbed2(xa)*intbed2(xb)<=0.0_dp) THEN
              xbed2 = zbrent(intbed2,xa,xb,1.0e-8_dp)
           ELSE
              xbed2 = xn(Nn)
           END IF
           write(*,*)'Case 1',dt3,xbed1,xbed2
           
!
! CASE 2 : intersecting node on right 
!
         ELSE IF ((nbed10==1).AND.(nbed2<Nn)) THEN
           dt3 = dt2
           IF (dt3 > EPS) THEN
             IF (FirstNodeMove) THEN
               aa = xn(2)
               bb = yn(2)
               cc = FlowValues( k*( FlowPerm( ordre(2) )-1 )+1) 
               dd = FlowValues( k*( FlowPerm( ordre(2) )-1 )+2) 
               xa = xn0(1)
               xb = aa    
               IF (intbed2(xa)*intbed2(xb)<=0.0) THEN
                 xbed1 = zbrent(intbed2,xa,xb,1.0e-8_dp)
               ELSE
                  xbed1 = xn(1)
               END IF
             END IF

             bb = xn0(nbed2) - xn(nbed2)
             aa = (yn0(nbed2) - yn(nbed2))/bb
             bb = (yn(nbed2)*xn0(nbed2)-xn(nbed2)*yn0(nbed2))/bb
             xbed2 = zbrent(intbed,xn0(nbed2),xn(nbed2),1.0e-08_dp)
           
           ELSE
            xbed1 = xn0(1)
            xbed2 = xn0(nbed2)
            nbed1 = nbed1 + 1
           END IF

           write(*,*)'Case 2',dt3,xbed1,xbed2

!
! CASE 3 : intersecting node on left 
!
         ELSE IF ((nbed10>1).AND.(nbed2==Nn)) THEN
           dt3 = dt1
            
           IF (dt3 > EPS) THEN
             IF (FirstNodeMove) THEN
               bb = xn0(nbed1) - xn(nbed1)
               aa = (yn0(nbed1) - yn(nbed1))/bb
               bb = (yn(nbed1)*xn0(nbed1)-xn(nbed1)*yn0(nbed1))/bb
               xbed1 = zbrent(intbed,xn0(nbed1),xn(nbed1),1.0e-08_dp)
             END IF
             aa = xn(Nn-1)
             bb = yn(Nn-1) 
             cc = FlowValues( k*( FlowPerm( ordre(Nn-1) )-1 )+1) 
             dd = FlowValues( k*( FlowPerm( ordre(Nn-1) )-1 )+2) 
             xa = aa        
             xb = xn(Nn)
             IF (intbed2(xa)*intbed2(xb)<=0.0) THEN
                xbed2 = zbrent(intbed2,xa,xb,1.0e-8_dp)
             ELSE
                xbed2 = xn(Nn)
             END IF
           ELSE
             xbed1 = xn0(nbed1)
             xbed2 = xn0(Nn)
             nbed1 = nbed1 + 1
           END IF

           write(*,*)'Case 3',dt3,xbed1,xbed2
!
! CASE 4 : intersecting node on left and right 
!
         ELSE  ! nbed1 > 1 and nbed2 < Nn
           dt3 = MIN(dt1,dt2) 
           
           IF (dt3 < EPS) THEN 
             IF (FirstNodeMove) THEN 
               bb = xn0(nbed1) - xn(nbed1)
               aa = (yn0(nbed1) - yn(nbed1))/bb
               bb = (yn(nbed1)*xn0(nbed1)-xn(nbed1)*yn0(nbed1))/bb
               xbed1 = zbrent(intbed,xn0(nbed1),xn(nbed1),1.0e-08_dp)
             END IF

             bb = xn0(nbed2) - xn(nbed2)
             aa = (yn0(nbed2) - yn(nbed2))/bb
             bb = (yn(nbed2)*xn0(nbed2)-xn(nbed2)*yn0(nbed2))/bb
             xbed2 = zbrent(intbed,xn0(nbed2),xn(nbed2),1.0e-08_dp)

           ELSE
             xbed1 = xn0(nbed1)
             xbed2 = xn0(nbed2)
             nbed1 = nbed1 + 1
           END IF

           write(*,*)'Case 4',dt3,xbed1,xbed2

         END IF

         IF (nbed1==1) nbed1=0 


          ALLOCATE( xnew(nbed2+1-nbed1),ynew(nbed2+1-nbed1), &
                   y1new(nbed2+1-nbed1), y2new(nbed2+1-nbed1), &
                   y3new(nbed2+1-nbed1) )
! new position of nbed nodes
          xnew(1) = xbed1 
          ynew(1) = fbed(xbed1) 
          k = FlowVariable % DOFs
          n = 2
          lc =0.0_dp
          DO p = 1+nbed1, nbed2
            Velocity = 0.0_dp
            DO i=1,k-1
            Velocity(i) = FlowValues( k*( FlowPerm( ordre(p) )-1 ) +i ) 
            END DO
            
            xnew(n) = xn0(p) + dt3 * Velocity(1)
            ynew(n) = yn0(p) + dt3 * Velocity(2)



! IF nbed=Nn -> move last node exactly on the bed             
            IF (p==nbed2)  THEN
               ynew(n)=fbed(xbed2)
            END IF

! Cavity length 
            lc = lc + SQRT((xnew(n)-xnew(n-1))**2 + &
                                 (ynew(n)-ynew(n-1))**2 ) 

            n = n +1
          END DO

! First derivative of the cavity
          y1new(1) = dbed(xbed1)         
          y1new(nbed2+1-nbed1) = dbed(xbed2)         
          DO n = 2, nbed2-nbed1
            y1new(n) = (ynew(n+1) - ynew(n-1))/(xnew(n+1) - xnew(n-1))
          END DO

! New cavity lenght and space between nodes
          
          ds = lc / (Nn - 1)

! New position of the nodes on the cavity surface
          dMU = 0.0
          CALL spline(xnew,ynew,dbed(xbed1),dbed(xbed2),y2new)
          CALL spline(xnew,y1new,1.0e40_dp,1.0e40_dp,y3new)
          x = xbed1 
          y = fbed(x)           
          dMU(2*ordre(1)-1) = x - xn0(1)
          dMU(2*ordre(1)) = y - yn0(1)
          

          DO p = 2, Nn 

            tanbed(1) =  1.0_dp            
            tanbed(2) =  splint(xnew,y1new,y3new,x)
            tanbed = tanbed / SQRT(tanbed(1)**2 + tanbed(2)**2)
            tanbed(1) = ABS(tanbed(1)) 
            IF (tanbed(1) < EPS) tanbed(1) = EPS/ds
            x = x + ds*tanbed(1)

! Again take care that node Nn is exactly on bed
             
            IF (p==Nn) THEN
              x = xbed2
              y = fbed(x)           
            ELSE
              y = splint(xnew,ynew,y2new,x) 
            ENDIF

! Take care that no node are bellow the bed

            IF (y < fbed(x)) THEN
               write(*,*)'y < fbed',p,x,y, fbed(x)
               y = fbed(x) + EPS
            END IF

            dMU(2*ordre(p)-1) = x - xn0(p)
            dMU(2*ordre(p)) = y - yn0(p)
         
          END DO

          dMUFN = curvi(xn0(1),xbed1,fbed,20)
          dMULN = curvi(xn0(Nn),xbed2,fbed,20)


          TrueTime = TrueTime + dt3
!
! Change time to get next time right
!
          TimeVar % Values(1) =  LastTime -dt + dt3 

          DEALLOCATE(xnew,ynew,y2new)


! ----------------------------------
! Compute the Mesh Update value
! ----------------------------------
!
        DO p = 1, Model % NumberOfBoundaryElements
          BCElement => GetBoundaryElement(p) 
          IF (GetElementFamily(BCElement) == 1 ) CYCLE
          IF (BCType(p)==0) CYCLE 

          n = BCElement % Type % NumberOfNodes
          
          DO i=1, n
            mu1 = 0.0
            mu2 = 0.0 
            Node = BCElement % NodeIndexes(i)
! Nodes on the Free Surface
            IF (BCType(p)==1) THEN
               ! Read the water pressure (constant in space) 
               IF (Node==ordre(1)) THEN
                 BC => GetBC( BCElement ) 
                 IF (FlowSolverName=='aiflow') THEN
                  Sxx(1:n) =  GetReal(BC, 'Normal Force', GotIt)
                 ELSE IF (FlowSolverName=='flow solution') THEN
                  Sxx(1:n) = GetReal(BC, 'External Pressure', GotIt)
                 END IF
                  pwater = -Sxx(1)
               END IF

               mu1 = dMU(2*Node-1) + &
                        MHPrev( DIM*MHPerm( Node ) - 1 , 1 )
               mu2 = dMU(2*Node) + &
                        MHPrev( DIM*MHPerm( Node ) , 1 )
              IF ((.Not.FirstNodeMove).AND.(Node==ordre(1))) THEN
                mu1 = 0.0_dp  
                mu2 = 0.0_dp            
              END IF 

! Nodes before and after the free surface
            ELSE IF ((BCType(p)==2).OR.(BCType(p)==3)) THEN
              IF ((x0(Node)>=xEnd).And.(xf-x0(Node)>EPS)) THEN
                mu1 = MHPrev( DIM*( MHPerm( Node ) -1 ) + 1 , 1 )
                x = Model % Nodes % x( Node )
                tanbed(1) = 1.0_dp
                tanbed(2) = dbed(x)
                tanbed = tanbed / sqrt(tanbed(1)**2 + tanbed(2)**2)
                tanbed(1) = ABS(tanbed(1)) 
                IF (tanbed(1)>EPS) THEN          
                  mu1 = mu1 + dMULN * tanbed(1) * &
                                    curvi(x0(Node),xf,fbed,100)/lbed20 
                ELSE
                  mu1 = mu1 + EPS
                END IF
                x = x0(Node) + mu1
                y = fbed(x)
                mu2 = y - y0( Node ) 
                IF (Node==ordre(Nn)) THEN
                  mu1 = MHPrev( DIM*( MHPerm( Node ) -1 ) + 1, 1 ) + &
                              dMU( 2*(Node -1) + 1) 
                  mu2 = MHPrev( DIM*( MHPerm( Node ) -1 ) + 2, 1 ) + &
                              dMU( 2*(Node -1) + 2) 
                END IF

! Nodes before the free surface
              ELSE IF ((x0(Node)<=xFirst).AND.(x0(Node)-xd>EPS) &
                                        .AND.(FirstNodeMove)) THEN
                mu1 = MHPrev( DIM*( MHPerm( Node ) -1 ) + 1 , 1 )
                x = Model % Nodes % x( Node )
                tanbed(1) = 1.0_dp 
                tanbed(2) = dbed(x)
                tanbed = tanbed / sqrt(tanbed(1)**2 + tanbed(2)**2)
                tanbed(1) = ABS(tanbed(1)) 
                IF (tanbed(1)<EPS) tanbed(1)=EPS
                mu1 = mu1 + dMUFN * tanbed(1) * & 
                                 curvi(xd,x0(Node),fbed,100)/lbed10 
                x = x0(Node) + mu1
                y = fbed(x)
                mu2 = y - y0( Node ) 

                IF (Node==ordre(1)) THEN
                  mu1 = MHPrev( DIM*( MHPerm( Node ) -1 ) + 1, 1 ) + &
                           dMU( 2*(Node -1) + 1) 
                  mu2 = MHPrev( DIM*( MHPerm( Node ) -1 ) + 2, 1 ) + &
                           dMU( 2*(Node -1) + 2) 
                END IF
              END IF  
            END IF
              moveFSVector(DIM*moveFSPerm(Node) - 1) = mu1 
              moveFSVector(DIM*moveFSPerm(Node) ) = mu2 

          END DO

        END DO

! ----------------------------------
! Compute output (taub, us, pi, ...)
! ----------------------------------
!
      DevStressVar => &
             VariableGet(Solver % Mesh % Variables,'DeviatoricStress')
      IF ( ASSOCIATED( DevStressVar ) ) THEN
      DSPerm => DevStressVar % Perm    
      DSValues => DevStressVar % Values  
      END IF

      sc = 0.0_dp
      taub = 0.0_dp
      ub = 0.0_dp
      pice = 0.0_dp
      fx = 0.0_dp
      fy = 0.0_dp
      ltot=0.0_dp 
      slope = 0.0_dp
      DO p = 1, Model % NumberOfBoundaryElements
        BCElement => GetBoundaryElement(p) 
        IF (GetElementFamily(BCElement) == 1 ) CYCLE
        IF (BCType(p)==0) CYCLE 
        n = BCElement % Type % NumberOfNodes

        CALL GetElementNodes( Nodes, BCElement )

        NodeIndexes => BCElement % NodeIndexes

        Nodes % x(1:n) = Model % Nodes % x(NodeIndexes(1:n))
        Nodes % y(1:n) = Model % Nodes % y(NodeIndexes(1:n))
        Nodes % z(1:n) = Model % Nodes % z(NodeIndexes(1:n))


        Sxx(1:n) = DSValues(2*dim*(DSPerm(NodeIndexes(1:n))-1) + 1) -&
                  FlowValues( (dim+1)*FlowPerm(NodeIndexes(1:n)) ) 
        Syy(1:n) = DSValues(2*dim*(DSPerm(NodeIndexes(1:n))-1) + 2) -&
                  FlowValues( (dim+1)*FlowPerm(NodeIndexes(1:n)) ) 
        Sxy(1:n) = DSValues(2*dim*(DSPerm(NodeIndexes(1:n))-1) + 4) 

        Velo = 0.0d0
        LocalNormal=0.0_dp
        DO i=1,DIM         
         LocalNormal(i,1:n)=Nvector(DIM*(NodeIndexes(1:n)-1)+i)
         Velo(i,1:n) = FlowValues((DIM+1)*(FlowPerm(NodeIndexes(1:n))-1)+i) 
        END DO
        Call BCIntegrals(ub, taub, pice, fx, fy, ltot, sc, &
                     Sxx, Syy, Sxy, Velo, LocalNormal,  BCType(p),  &
                     BCElement, n, Nodes)
                            

        DO i = 1, n
          k = NodeIndexes( i )
          IF ( SQRT(SUM(Nvector(DIM*(k-1)+1:DIM*k)**2)) /= 0.0_dp ) THEN
            IF (Nvector(2*k)/=0.0_dp) THEN
              IF (-Nvector(2*k-1)/Nvector(2*k)>slope) &
                           slope = -Nvector(2*k-1)/Nvector(2*k)
            ELSE
              IF (Nvector(2*k-1)>0.0) slope = Large
            END IF
          END IF
        END DO

      END DO
      taub = taub / ltot          
      pice = pice / ltot
      ub = ub / ltot       
      deltaFS = SQRT(deltaFS) / ( Nn * ub ) 
      MaxdeltaFS = SQRT(MaxdeltaFS) / ub

      IF (SaveFile) THEN
        OPEN(22,File=SaveFileName,POSITION='APPEND') 
        WRITE(22,'(17(e14.8,2x))')LastTime, TrueTime, dt3, taub, & 
             pice, pwater, ub, taub/(pice-pwater),  ub/(pice - pwater),&
             fx, fy, deltaFS, MaxdeltaFS, slope, xbed1, xbed2, sc
                 
        CLOSE(22)
      CALL Info( 'CavitySolve', '---------------------',Level=4 )
      WRITE( Message, * ) ' dFS = V.n , Max(V.n)  : ',deltaFS, MaxdeltaFS
      CALL Info( 'CavitySolve', Message, Level=4 )
      CALL Info( 'CavitySolve', '---------------------',Level=4 )
      END IF
     
! ----------------------------------

CONTAINS

!
! --------------------------------------------------------------------------------------
! Integration over the Boundary Element of some data
! --------------------------------------------------------------------------------------
!
         SUBROUTINE BCIntegrals(ub, tau, pice, fx, fy, l, sc, &
             Sxx, Syy, Sxy, Velo, LocalNormal, BCn, Element, n, ElementNodes)  
                                
 
!          USE types
!          USE DefUtils
!          USE CoordinateSystems
!          USE SolverUtils
!          USE ElementDescription
!          IMPLICIT NONE
           REAL(KIND=dp) :: ub, tau, pice,  fx, fy, l, sc
           REAL(KIND=dp) :: Sxx(3), Syy(3), Sxy(3), Velo(3,3), &
                            LocalNormal(3,3)
           TYPE(Element_t),POINTER  :: Element
           TYPE(Nodes_t)    :: ElementNodes
           INTEGER :: BCn, n
!------------------------------------------------------------------------------
           REAL(KIND=dp) :: Basis(n),ddBasisddx(1,1,1)
           REAL(KIND=dp) :: dBasisdx(n,3),SqrtElementMetric
           REAL(KIND=dp) :: u,v,w,s
           REAL(KIND=dp) :: Sig(3,3)
           REAL(KIND=dp), POINTER :: U_Integ(:),V_Integ(:),W_Integ(:),S_Integ(:)
           REAL(KIND=dp) :: Velocity(3), Snn, norma(3)
           INTEGER :: i,j,t,q,p,DIM,N_Integ
           LOGICAL :: stat
           TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
           
!------------------------------------------------------------------------------
           DIM = CoordinateSystemDimension()
!
!  Integration stuff
!
           IntegStuff = GaussPoints( Element )
           U_Integ => IntegStuff % u
           V_Integ => IntegStuff % v
           W_Integ => IntegStuff % w
           S_Integ => IntegStuff % s
           N_Integ =  IntegStuff % n
!
!  Now we start integrating
!
           DO t=1,N_Integ

             u = U_Integ(t)
             v = V_Integ(t)
             w = W_Integ(t)

!------------------------------------------------------------------------------
!    Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
              stat = ElementInfo( Element, ElementNodes, u, v, w, &
                SqrtElementMetric, Basis, dBasisdx, ddBasisddx, .FALSE.)
                

              s = SqrtElementMetric * S_Integ(t)
           
              Sig = 0.0
              Sig(1,1) = SUM(Sxx(1:n)*Basis(1:n)) 
              Sig(2,2) = SUM(Syy(1:n)*Basis(1:n)) 
              Sig(1,2) = SUM(Sxy(1:n)*Basis(1:n)) 
              Sig(2,1) = Sig(1,2)
! Normal Upward-pointing normal
!             Normal = - NormalVector( Element,ElementNodes,u,v,.TRUE. )
              Snn = 0.0
              Velocity = 0.0
              normal = 0.0
              DO i=1, DIM
                Velocity(i) = SUM(Velo(i,1:n)*Basis(1:n) )
                normal(i) = -SUM(LocalNormal(i,1:n)*Basis(1:n) )
              END DO

              DO i=1, DIM
                DO j=1, DIM
                  Snn = Snn + Sig(i,j)*Normal(i)*Normal(j) 
                END DO
              END DO
              
              ub = ub + Velocity(1) * Normal(2) * s
              l = l + Normal(2) * s 
              tau  = tau + Snn * Normal(1) * s 
              pice = pice - Snn * Normal(2) * s 
             
              IF (BCn==1)  THEN
                sc = sc + s
              END IF
              
              IF (BCn==3) THEN
                fx = fx + Snn * Normal(1) * s 
                fy = fy - Snn * Normal(2) * s 
              END IF

         END DO
         END SUBROUTINE BCINtegrals
!
! --------------------------------------------------------------------------------------


END SUBROUTINE CavitySolve



! --------------------------------------------------------
     FUNCTION intbed(x)
       USE types
       USE DefUtils
     REAL(KIND=dp) :: intbed,x, fbed
     REAL(KIND=dp) :: aa, bb, cc, dd
     COMMON /interbed/aa,bb,cc,dd

     intbed = fbed(x) - aa*x -bb 
      

     END FUNCTION intbed
!
      FUNCTION intbed2(x)
       USE types
       USE DefUtils
       REAL(KIND=dp) :: intbed2, x, fbed, dbed
       REAL(KIND=dp) :: aa, bb, cc, dd
       REAL(KIND=dp) :: EPS, dx
       COMMON /interbed/aa,bb,cc,dd
! last node function
      EPS = 100.0_dp*EPSILON(aa)
      IF (cc < 0.0) THEN
         cc = -cc
         dd = -dd
      ENDIF
      IF (cc < EPS) cc = EPS
      dx = x - aa
      IF (ABS(dx) < EPS) dx = EPS
      intbed2 = 2.0_dp/dx*(fbed(x)-bb)-dd/cc-dbed(x) 
      END FUNCTION intbed2

! --------------------------------------------------------
      FUNCTION curvi(x0,x1,func,N)
       USE types
       USE DefUtils
       REAL(KIND=dp) :: x0, x1, x,  dx, l
       REAL(KIND=dp) ::  dl1, dl2, curvi
       INTEGER :: i, N
       INTERFACE 
         FUNCTION func(x) 
         USE types
         USE DefUtils
         IMPLICIT NONE 
         REAL(KIND=dp), INTENT(IN) :: x 
         REAL(KIND=dp) :: func 
         END FUNCTION func 
       END INTERFACE 
!
! Calculate the length of a curve
! l > 0 if x0 < x1, else l < 0 
!
!
       dx = (x1-x0)/N
       l = 0.0_dp
       x = x0
       y = func(x) 
         DO i=2,N+1
         dl1 = x 
         dl2 = y
         x = x + dx
         y = func(x)
         dl1 = x - dl1
         dl2 = y - dl2
         l = l + SQRT(dl1**2 + dl2**2)
         END DO
      IF (x0 < x1) THEN
        curvi = l
      ELSE
        curvi = -l
      END IF
      END FUNCTION curvi

! --------------------------------------------------------
      FUNCTION intbed1(x)
       USE types
       USE DefUtils
     REAL(KIND=dp) :: intbed1,x, dbed, fbed, Norm
     REAL(KIND=dp) :: normal(2),XM(2)
     REAL(KIND=dp) :: aa, bb, cc, dd
     COMMON /interbed/aa,bb,cc,dd
     COMMON /interbed/a,b

! function that return f(x) = XM.n X(x,fbed), M(a,b), n (dbed,-1)
    
     ! Norm of XM 
      XM(1) = aa - x
      XM(2) = bb - fbed(x)
      Norm = sqrt(XM(1)**2+XM(2)**2)
      XM = XM/Norm

      normal(1) = dbed(x) 
      normal(2) = -1.0       
      Norm = sqrt(normal(1)**2+normal(2)**2)
      normal = normal/Norm

      intbed1 = SUM(normal*XM) 

      END FUNCTION intbed1
! -------------------------------------------------------------
! Numerical recipes subroutines  and functions          -------
! -------------------------------------------------------------
       SUBROUTINE indexx(arr,ind) 
       USE types
       USE DefUtils
       IMPLICIT NONE 
       REAL(KIND=dp) :: arr(:) 
       INTEGER :: ind(:) 
       INTEGER, PARAMETER :: NN=15, NSTACK=50
       REAL(KIND=dp) :: a 
       INTEGER :: n,k,i,j,indext,jstack,l,r 
       INTEGER, DIMENSION(NSTACK) :: istack 

        n = size(arr)
        Do i = 1, n
        ind( i ) = i 
        END DO
jstack=0 
l=1 
r=n 
do 
if (r-l < NN) then 
do j=l+1,r 
indext=ind(j) 
a=arr(indext) 
do i=j-1,l,-1 
if (arr(ind(i)) <= a) exit 
ind(i+1)=ind(i) 
end do 
ind(i+1)=indext 
end do 
if (jstack == 0) RETURN 
r=istack(jstack) 
l=istack(jstack-1) 
jstack=jstack-2 
else 
k=(l+r)/2 
call swap(ind(k),ind(l+1)) 
call icomp_xchg(ind(l),ind(r)) 
call icomp_xchg(ind(l+1),ind(r)) 
call icomp_xchg(ind(l),ind(l+1)) 
i=l+1 
j=r 
indext=ind(l+1) 
a=arr(indext) 
do 
do 
i=i+1 
if (arr(ind(i)) >= a) exit 
end do 
do 
j=j-1 
if (arr(ind(j)) <= a) exit 
end do 
if (j < i) exit 
call swap(ind(i),ind(j)) 
end do 
ind(l+1)=ind(j) 
ind(j)=indext 
jstack=jstack+2 
if (jstack > NSTACK) Write(*,*)' indexx: NSTACK too small'
if (r-i+1 >= j-l) then 
istack(jstack)=r
istack(jstack-1)=i 
r=j-1 
else 
istack(jstack)=j-1 
istack(jstack-1)=l 
l=i 
end if 
end if 
end do 
CONTAINS 
SUBROUTINE icomp_xchg(i,j) 
INTEGER :: i,j 
INTEGER :: swp 
if (arr(j) < arr(i)) then 
swp=i 
i=j 
j=swp 
end if 
END SUBROUTINE icomp_xchg 
SUBROUTINE swap(a,b) 
INTEGER :: a, b
INTEGER :: dum
dum = a
a = b
b = dum
END SUBROUTINE swap
END SUBROUTINE indexx
! --------------------------------------------------------
FUNCTION splint(xa,ya,y2a,x) 
       USE types
       USE DefUtils
IMPLICIT NONE 
REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: xa,ya,y2a 
REAL(KIND=dp), INTENT(IN) :: x 
REAL(KIND=dp) :: splint
INTEGER :: khi,klo,n
REAL(KIND=dp) :: a,b,h 
n = size(xa)
klo=max(min(locate(xa,x),n-1),1)
khi=klo+1 
h=xa(khi)-xa(klo) 
if (h == 0.0) THEN
    Write(*,*)' bad xa input in splint'  
    Write(*,*)khi,klo,xa(khi),xa(klo),x
End IF
a=(xa(khi)-x)/h 
b=(x-xa(klo))/h 
splint=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0_dp 
CONTAINS
! ----------------------------------------------------------------------------

FUNCTION locate(xx,x) 
       USE types
       USE DefUtils
IMPLICIT NONE 
REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: xx 
REAL(KIND=dp), INTENT(IN) :: x 
INTEGER :: locate 
INTEGER :: n,jl,jm,ju 
LOGICAL :: ascnd
n=size(xx) 
ascnd = (xx(n) >= xx(1)) 
jl=0 
ju=n+1  
do 
if (ju-jl <= 1) exit 
jm=(ju+jl)/2 
if (ascnd .eqv. (x >= xx(jm))) then 
jl=jm 
else 
ju=jm 
end if 
end do 
if (x == xx(1)) then 
locate=1 
else if (x == xx(n)) then 
locate=n-1 
else 
locate=jl 
end if 
END FUNCTION locate

END FUNCTION splint
! ----------------------------------------------------------------------------

SUBROUTINE spline(x,y,yp1,ypn,y2) 
       USE types
       USE DefUtils
IMPLICIT NONE 
REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: x,y 
REAL(KIND=dp), INTENT(IN) :: yp1,ypn 
REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: y2 
INTEGER :: n 
REAL(KIND=dp), DIMENSION(size(x)) :: a,b,c,r 
INTERFACE
SUBROUTINE tridag(a,b,c,r,u) 
       USE types
       USE DefUtils
IMPLICIT NONE 
REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: a,b,c,r 
REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: u 
END SUBROUTINE tridag
END INTERFACE
n = size(x)
c(1:n-1)=x(2:n)-x(1:n-1) 
r(1:n-1)=6.0_dp*((y(2:n)-y(1:n-1))/c(1:n-1)) 
r(2:n-1)=r(2:n-1)-r(1:n-2) 
a(2:n-1)=c(1:n-2) 
b(2:n-1)=2.0_dp*(c(2:n-1)+a(2:n-1)) 
b(1)=1.0 
b(n)=1.0 
if (yp1 > 0.99e30_dp) then 
r(1)=0.0 
c(1)=0.0 
else 
r(1)=(3.0_dp/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
c(1)=0.5 
end if 
if (ypn > 0.99e30_dp) then 
r(n)=0.0 
a(n)=0.0 
else 
r(n)=(-3.0_dp/(x(n)-x(n-1)))*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn) 
a(n)=0.5 
end if 
call tridag(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n)) 
END SUBROUTINE spline

! ----------------------------------------------------------------------------

SUBROUTINE tridag(a,b,c,r,u) 
       USE types
       USE DefUtils
IMPLICIT NONE 
REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: a,b,c,r 
REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: u 
REAL(KIND=dp), DIMENSION(size(b)) :: gam
INTEGER :: n,j 
REAL(KIND=dp) :: bet 
n = size(b)
bet=b(1) 
if (bet == 0.0) Write(*,*)' tridag_ser: Error at code stage 1 ' 
u(1)=r(1)/bet 
do j=2,n 
gam(j)=c(j-1)/bet 
bet=b(j)-a(j-1)*gam(j) 
if (bet == 0.0) Write(*,*)' tridag_ser: Error at code stage 2 ' 
u(j)=(r(j)-a(j-1)*u(j-1))/bet 
end do 
do j=n-1,1,-1 
u(j)=u(j)-gam(j+1)*u(j+1) 
end do 
END SUBROUTINE tridag

! ----------------------------------------------------------------------------
FUNCTION zbrent(func,x1,x2,tol)
       USE types
       USE DefUtils
IMPLICIT NONE 
REAL(KIND=dp), INTENT(IN) :: x1,x2,tol 
REAL(KIND=dp) :: zbrent 
INTERFACE 
FUNCTION func(x) 
       USE types
       USE DefUtils
IMPLICIT NONE 
REAL(KIND=dp), INTENT(IN) :: x 
REAL(KIND=dp) :: func 
END FUNCTION func 
END INTERFACE 
INTEGER, PARAMETER :: ITMAX=100 
REAL(KIND=dp), PARAMETER :: EPS=epsilon(x1)
INTEGER :: iter 
REAL(KIND=dp) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm 
a=x1 
b=x2 
fa=func(a) 
fb=func(b) 
if ((fa > 0.0 .and. fb > 0.0) .or. (fa < 0.0 .and. fb < 0.0)) THEN 
            write(*,*)' root must be bracketed for zbrent r' 
            write(*,*)a,b,fa,fb
END IF
c=b 
fc=fb 
do iter=1,ITMAX 
if ((fb > 0.0 .and. fc > 0.0) .or. (fb < 0.0 .and. fc < 0.0)) then 
c=a 
fc=fa 
d=b-a 
e=d 
end if 
if (abs(fc) < abs(fb)) then 
a=b 
b=c 
c=a 
fa=fb 
fb=fc 
fc=fa 
end if 
tol1=2.0_dp*EPS*abs(b)+0.5_dp*tol 
xm=0.5_dp*(c-b) 
if (abs(xm) <= tol1 .or. fb == 0.0) then 
zbrent=b 
RETURN 
end if 
if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then 
s=fb/fa  
if (a == c) then 
p=2.0_dp*xm*s 
q=1.0_dp-s
else 
q=fa/fc 
r=fb/fc 
p=s*(2.0_dp*xm*q*(q-r)-(b-a)*(r-1.0_dp)) 
q=(q-1.0_dp)*(r-1.0_dp)*(s-1.0_dp) 
end if 
if (p > 0.0) q=-q 
p=abs(p) 
if (2.0_dp*p < min(3.0_dp*xm*q-abs(tol1*q),abs(e*q))) then 
e=d 
d=p/q
else
d=xm
e=d 
end if 
else 
d=xm 
e=d 
end if 
a=b 
fa=fb 
b=b+merge(d,sign(tol1,xm), abs(d) > tol1 ) 
fb=func(b) 
end do 
Write(*,*)' zbrent: exceeded maximum iterations ' 
zbrent=b 

END FUNCTION zbrent
