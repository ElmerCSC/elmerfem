!/*****************************************************************************/
! *
! *  Elmer/Ice, a glaciological add-on to Elmer
! *  http://elmerice.elmerfem.org
! *
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
! ******************************************************************************
! *
! *  Authors: Olivier Gagliardini
! *  Email:   
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 
! * 
! *****************************************************************************
!>   Shape Factor user functions
!>   
!>  return the gravity force -g + (1-f)(g.t).t
!>  where t is the tangent vector of the surface
!>   end f is the shape factor
!>   work only in 2D (no sense in 3D)
!>   work for non-structured mesh
!>   f can be given or calculated as a function of H(x,t) (height)
!>   for different shape (rectangle, parabola, shelf)
FUNCTION ShapeFactorGravity_x ( Model, nodenumber, x) RESULT(gx)
   USE Types
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   INTEGER :: nodenumber
   REAL(KIND=dp) :: x, gx, ShapeFactorGravity

   gx = ShapeFactorGravity ( Model, nodenumber, x, 1 )
    

END FUNCTION ShapeFactorGravity_x

FUNCTION ShapeFactorGravity_y ( Model, nodenumber, x) RESULT(gy)
   USE Types
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   INTEGER :: nodenumber
   REAL(KIND=dp) :: x, gy, ShapeFactorGravity

   gy = ShapeFactorGravity ( Model, nodenumber, x, 2 )

END FUNCTION ShapeFactorGravity_y




FUNCTION ShapeFactorGravity ( Model, nodenumber, x, axis ) RESULT(gi)
   USE types
   USE GeneralUtils
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   INTEGER :: nodenumber, axis
   REAL(KIND=dp) :: x, gi

   TYPE(Nodes_t), SAVE :: Nodes
   TYPE(variable_t), POINTER :: Timevar
   TYPE(Element_t), POINTER ::  BoundaryElement, BCElement, CurElement, ParentElement
   TYPE(ValueList_t), POINTER :: BodyForce, BC


   REAL(KIND=dp), ALLOCATABLE :: zs(:), xs(:), ts(:,:), zb(:), xb(:)
   REAL(KIND=dp), TARGET, ALLOCATABLE :: z2s(:), z2b(:), ts2x(:), ts2y(:), para2(:)
   REAL(KIND=dp), DIMENSION(:), POINTER :: z2s_p, z2b_p, ts2x_p, ts2y_p, para2_p
   REAL(KIND=dp), ALLOCATABLE :: auxReal(:), para(:)
   REAL(KIND=dp) :: t, t0, told, Bu, Bv, f, afactor, width, H, w, NVelo
   REAL(KIND=dp) :: g(2), tangent(2),  normal(3), Velo(2)

   INTEGER, ALLOCATABLE :: ind(:), NodesOnSurface(:), NodesOnBed(:)
   INTEGER, ALLOCATABLE :: TempS(:), TempB(:)
   INTEGER :: Ns, Nb 
   INTEGER ::i, j, n, k, p, DIM 

   LOGICAL :: FirstTime = .TRUE., NewTime, GotIt  
   LOGICAL :: Parabola, Rectangle, Computed=.TRUE.
   LOGICAL :: Bedrock, Surface
       
!--------------------------------------------
! Parameter for the shape function definition
! f(w) = 2/ (n pi) sum_i=1^n ATAN(a(i) x w^i)
!--------------------------------------------
   INTEGER, PARAMETER :: nr=6, np=1                            
   REAL(KIND=dp), PARAMETER :: ap(np) = (/0.814589_dp/), &
                               ar(nr) = (/5.30904_dp,0.0934244_dp,3.94177_dp,74.6014_dp,0.0475342_dp,1.35021_dp/)

   SAVE told, t0, FirstTime, NewTime, g
   SAVE zs, xs, ts, zb, xb, Ns, Nb, para, auxReal
   SAVE z2s, z2b, ts2x, ts2y, para2 
   SAVE Parabola, Rectangle, Computed
   SAVE NodesOnSurface, NodesOnBed, TempS, TempB

   Timevar => VariableGet( Model % Variables,'Time')
   t = TimeVar % Values(1)

   IF (FirstTime) THEN
! STOP if we are not solving a 2D flow line problem
      DIM = CoordinateSystemDimension()
      IF ( DIM /= 2 ) THEN
         CALL FATAL('ShapeFactorGravity','Work only for 2D flow line problem')
      END IF

      FirstTime = .FALSE.
      NewTime = .TRUE.
      told = t
      t0 = t

      n = Model % MaxElementNodes 
      ALLOCATE( auxReal(n) )
      
      BodyForce => GetBodyForce()  

!-------------------------
! Get the gravity vector g                    
!-------------------------
      g(1) = GetConstReal( BodyForce, 'Shape Gravity 1', GotIt )
      g(2) = GetConstReal( BodyForce, 'Shape Gravity 2', GotIt )
         IF ( .Not.GotIt ) CALL FATAL('ShapeFactorGravity', &
                'Gravity vector must be specified')

!---------------------------------------------
! Determine how is calculated the shape factor                    
!---------------------------------------------
      Parabola = GetLogical( BodyForce, 'Parabola Shape', GotIt)
      IF (.Not.GotIt) Parabola = .FALSE.
      Rectangle = GetLogical( BodyForce, 'Rectangle Shape', GotIt)
      IF (.Not.GotIt) Rectangle = .FALSE.
      IF (Rectangle.AND.Parabola) CALL FATAL('ShapeFactorGravity', & 
                'Make a choice between Parabola or Rectangle shapes')
      IF ((.Not.Rectangle).AND.(.Not.Parabola)) THEN
          Computed = .FALSE.
          CurElement => Model % CurrentElement
          n = CurElement % Type % NumberOfNodes   
          auxReal(1:n) = GetReal( BodyForce, &
                             'Shape Factor', GotIt )
          IF ( .Not.GotIt ) CALL FATAL('ShapeFactorGravity', &
                'Shape Factor must be given or Computed')
      END IF


!-------------------------------------------
! Number of nodes on the bed and the surface  
!-------------------------------------------
      CurElement => Model % CurrentElement
      Ns = 1
      Nb = 1
      DO p = 1, Model % NumberOfBoundaryElements
         BCElement => GetBoundaryElement(p) 
         BC => GetBC( BCElement ) 
         IF( GetElementFamily(BCElement) == 1 ) CYCLE
         Surface = .FALSE.
         Bedrock = .FALSE.
         Surface = GetLogical( BC, 'Shape Surface', GotIt)
         Bedrock = GetLogical( BC, 'Shape Bedrock', GotIt)
         IF (Surface) THEN
            n = BCElement % Type % NumberOfNodes
            Ns = Ns + n - 1
         ELSE IF (Bedrock) THEN
            n = BCElement % Type % NumberOfNodes
            Nb = Nb + n - 1
         END IF
      END DO
      Model % CurrentElement => CurElement 
      
      IF ((Ns==0).OR.(Nb==0)) CALL FATAL('ShapeFactorGravity', &
                'Shape Surface and/or Shape Bedrock not defined')

      ALLOCATE ( zs(Ns), ts(2,Ns), xs(Ns), xb(Nb), zb(Nb) ) 
      ALLOCATE (z2s(Ns), z2b(Ns), ts2x(Ns), ts2y(Ns) ) 
      ALLOCATE ( NodesOnSurface(Ns), NodesOnBed(Nb), TempS(Ns), TempB(Nb) ) 
      IF (Computed) ALLOCATE( para(Ns), para2(Ns) ) 


!-----------------------------------------------------
! Store the nodes number on the bed and on the surface  
!-----------------------------------------------------
      CurElement => Model % CurrentElement
      Ns = 0
      Nb = 0
      NodesOnSurface = 0
      NodesOnBed = 0
      DO p = 1, Model % NumberOfBoundaryElements
         BCElement => GetBoundaryElement(p) 
         BC => GetBC( BCElement ) 
         IF( GetElementFamily(BCElement) == 1 ) CYCLE
         Surface = .FALSE.
         Bedrock = .FALSE.
         Surface = GetLogical( BC, 'Shape Surface', GotIt)
         Bedrock = GetLogical( BC, 'Shape Bedrock', GotIt)
         IF (Surface) THEN
            n = BCElement % Type % NumberOfNodes
            DO i = 1, n
               j = BCElement % NodeIndexes(i)
               IF (.NOT.ANY(NodesOnSurface == j )) THEN
                  Ns = Ns + 1
                  NodesOnSurface(Ns) = j                           
                  xs(Ns) = Model % Nodes % x( j )
               END IF
            END DO
         ELSE IF (Bedrock) THEN
            n = BCElement % Type % NumberOfNodes
            DO i = 1, n
               j = BCElement % NodeIndexes(i)
               IF (.NOT.ANY(NodesOnBed == j )) THEN
                  Nb = Nb + 1
                  NodesOnBed(Nb) = j                           
                  xb(Nb) = Model % Nodes % x( j )
               END IF
            END DO
         END IF
      END DO
      Model % CurrentElement => CurElement 
       
!------------------------------------------
! Get ordered NodesOnSurface and NodesOnBed  
! SortD : sort in increase order
!------------------------------------------
      ALLOCATE (ind(Ns))
      ind = (/(i,i=1,Ns)/) 
      CALL SortD(Ns,xs,ind)
      NodesOnSurface = NodesOnSurface(ind)
      DEALLOCATE(ind)
      ALLOCATE (ind(Nb))
      ind = (/(i,i=1,Nb)/) 
      CALL SortD(Nb,xb,ind)
      NodesOnBed = NodesOnBed(ind)
      DEALLOCATE(ind)
      
!----------------------------------------
! If Computed
! Get the width of the rectangle channel       
! Get the afactor of the parabola channel
!----------------------------------------
      IF (Computed) THEN
         CurElement => Model % CurrentElement
         DO p = 1, Model % NumberOfBoundaryElements
            BCElement => GetBoundaryElement(p) 
            BC => GetBC( BCElement ) 
            IF( GetElementFamily(BCElement) == 1 ) CYCLE
            Surface = .FALSE.
            Surface = GetLogical( BC, 'Shape Surface', GotIt)
            IF (Surface) THEN
                n = BCElement % Type % NumberOfNodes
                IF (Parabola) THEN
                   auxReal(1:n) = GetReal( BodyForce, &
                                 'Parabola aFactor', GotIt ) 
                   IF (.Not.GotIt) CALL FATAL('ShapeFactorGravity', &
                        'Parabola aFactor must be given for Parabola shape') 
                ELSE IF (Rectangle) THEN
                   auxReal(1:n) = GetReal( BodyForce, &
                                   'Rectangle Width', GotIt ) 
                   IF ((.Not.GotIt).AND.Rectangle) CALL FATAL('ShapeFactorGravity', &
                        'Rectangle Width must be given for Rectangle shape') 
                END IF
                DO i = 1, n
                   j = BCElement % NodeIndexes(i)
                   TempS = 0
                   WHERE (NodesOnSurface==j) TempS = 1
                   k = MAXLOC(TempS, dim=1)
                   para(k) = auxReal(i)
                END DO
            END IF
         END DO
         Model % CurrentElement => CurElement 
      END IF

   ELSE
      IF (t > told) THEN
         NewTime = .TRUE.
         told = t
      END IF
   ENDIF
   
   IF (NewTime) THEN
      NewTime = .FALSE.
!---------------------------------------------
! Nodal value of zs, zb, tangent for that time
!---------------------------------------------
      ts = 0.0_dp
      CurElement => Model % CurrentElement
      DO p = 1, Model % NumberOfBoundaryElements
         BCElement => GetBoundaryElement(p) 
         BC => GetBC( BCElement ) 
         IF( GetElementFamily(BCElement) == 1 ) CYCLE
         Surface = .FALSE.
         Bedrock = .FALSE.
         Surface = GetLogical( BC, 'Shape Surface', GotIt)
         Bedrock = GetLogical( BC, 'Shape Bedrock', GotIt)
         IF (Surface) THEN
           n = BCElement % Type % NumberOfNodes
           CALL GetElementNodes( Nodes , BCElement )
           DO i = 1,n
             j = BCElement % NodeIndexes( i )
             Bu = BCElement % Type % NodeU(i)
             Bv = 0.0D0
             Normal  = NormalVector(BCElement, Nodes, Bu, Bv, .TRUE.)
             TempS = 0
             WHERE (NodesOnSurface==j) TempS = 1
             k = MAXLOC(TempS, dim=1)
             ts(1,k) = ts(1,k) + Normal(2)
             ts(2,k) = ts(2,k) - Normal(1)
             xs(k) = Model % Nodes % x(j)
             zs(k) = Model % Nodes % y(j)
           END DO
         ELSE IF(Bedrock) THEN
           n = BCElement % Type % NumberOfNodes
           DO i = 1,n
             j = BCElement % NodeIndexes( i )
             TempB = 0
             WHERE (NodesOnBed==j) TempB = 1
             k = MAXLOC(TempB, dim=1)
             xb(k) = Model % Nodes % x(j)
             zb(k) = Model % Nodes % y(j)
           END DO
         END IF
      END DO

      DO i=1, Ns                        
         ts(:,i) = ts(:,i) / SQRT(SUM(ts(:,i)**2.0)) 
      END DO
      Model % CurrentElement => CurElement 

!-----------------------------------
! Construct the spline interpolation 
!-----------------------------------
       CALL CubicSpline(Ns,xs,zs,z2s)
       CALL CubicSpline(Nb,xb,zb,z2b)
       CALL CubicSpline(Ns,xs,ts(1,:),ts2x)
       CALL CubicSpline(Ns,xs,ts(2,:),ts2y)
       IF (Computed) CALL CubicSpline(Ns,xs,para,para2)

   ENDIF  ! new dt
    

     BodyForce => GetBodyForce()
     CurElement => Model % CurrentElement
     n = CurElement % Type % NumberOfNodes   
!-----------------------
! Interpolate for that x
!-----------------------
    x = Model % Nodes % x ( nodenumber)
    IF (Computed) THEN
       z2s_p => z2s
       z2b_p => z2b
       H = InterpolateCurve(xs,zs,x,z2s_p) - InterpolateCurve(xb,zb,x,z2b_p)
       IF (Parabola) THEN
               para2_p => para2
               afactor = InterpolateCurve(xs,para,x,para2_p)
               IF ((H>0.0).AND.(afactor>0.0)) THEN 
                   w = 1.0/SQRT(H*afactor)
                   f = 0.0_dp
                   DO i=1, np
                      f = f + ATAN(ap(i)*w**i)
                   END DO
                   f = 2.0*f/(np*pi) 
               ELSE
                   f = 1.0
               END IF    
       ELSE IF (Rectangle) THEN
               para2_p => para2
               width = InterpolateCurve(xs,para,x,para2_p)
               IF (H>0.0) THEN 
                  w = width / (2.0*H)
                  f = 0.0_dp
                  DO i=1, nr
                     f = f + ATAN(ar(i)*w**i)
                  END DO
                  f = 2.0*f/(nr*pi) 
               ELSE
                  f = 1.0
               END IF
       END IF
    ELSE
!------------------------------------
! Get the Shape factor for that nodes
!------------------------------------
         auxReal(1:n) = GetReal( BodyForce, &
                             'Shape Factor', GotIt )
         DO i=1, n
            j = CurElement % NodeIndexes (i)  
            IF (nodenumber == j) EXIT  
         END DO 
         f = auxReal(i) 
    END IF

    ts2x_p => ts2x
    ts2y_p => ts2y
    tangent(1) = InterpolateCurve(xs,ts(1,:),x,ts2x_p)
    tangent(2) = InterpolateCurve(xs,ts(2,:),x,ts2y_p)
    tangent = tangent / SQRT(SUM(tangent**2.0))
   
    gi = g(axis) - (1.0 - f) * SUM(g*tangent) * tangent(axis) 



END FUNCTION ShapeFactorGravity

