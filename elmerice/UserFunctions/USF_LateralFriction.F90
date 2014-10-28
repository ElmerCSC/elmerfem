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
! *  Date modifications:
! *
! * 
! *****************************************************************************
!>   Lateral Friction user functions
!>   
!>  return the gravity force -g + K * ||u||^(m-1) u  
!>  where K is the a lateral friction coefficient, 
!>        m the lateral friction exponent, 
!>   end  u is the velocity vector
!>   work only in 2D (no sense in 3D)
!>   work for non-structured mesh
FUNCTION LateralFriction_x ( Model, nodenumber, x) RESULT(gx)
   USE Types
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   INTEGER :: nodenumber
   REAL(KIND=dp) :: x, gx, LateralFriction

   gx = LateralFriction ( Model, nodenumber, x, 1 )
    

END FUNCTION LateralFriction_x

FUNCTION LateralFriction_y ( Model, nodenumber, x) RESULT(gy)
   USE Types
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   INTEGER :: nodenumber
   REAL(KIND=dp) :: x, gy, LateralFriction

   gy = LateralFriction ( Model, nodenumber, x, 2 )

END FUNCTION LateralFriction_y




FUNCTION LateralFriction ( Model, nodenumber, x, axis ) RESULT(gi)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   INTEGER :: nodenumber, axis
   REAL(KIND=dp) :: x, gi, Kspring, mm

   TYPE(Nodes_t), SAVE :: Nodes
   TYPE(Element_t), POINTER ::  CurElement
   TYPE(variable_t), POINTER :: FlowVariable 
   TYPE(ValueList_t), POINTER :: BodyForce, Material

   REAL(KIND=dp), POINTER :: FlowValues(:)
   INTEGER, POINTER :: FlowPerm(:)
    
   REAL(KIND=dp), ALLOCATABLE :: auxReal(:)
  
   REAL(KIND=dp) :: g(2), Velo(2), NVelo
   INTEGER :: i, j, n


   LOGICAL :: FirstTime = .TRUE., GotIt
       
   CHARACTER(LEN=MAX_NAME_LEN) :: FlowSolverName

   SAVE FirstTime, g, mm, FlowSolverName, auxReal

   BodyForce => GetBodyForce()  
   material => GetMaterial()

   IF (FirstTime) THEN

      FirstTime = .FALSE.

      n = Model % MaxElementNodes 
      ALLOCATE( auxReal(n) )

      !---------------------
      ! Get the gravity vector g                    
      !-------------------------
      g(1) = GetConstReal( BodyForce, 'Lateral Friction Gravity 1', GotIt )
      g(2) = GetConstReal( BodyForce, 'Lateral Friction Gravity 2', GotIt )
      IF (.Not.GotIt ) CALL FATAL('LateralFriction', &
                'Gravity vector must be specified')


      mm = GetConstReal( BodyForce, 'Lateral Friction Exponent', GotIt )
      IF (.Not.GotIt ) CALL FATAL('LateralFriction', &
                'Lateral Friction Exponent must be defined')

      FlowSolverName = GetString( BodyForce, 'Flow Solver Name', GotIt )    
      IF (.NOT.Gotit) FlowSolverName = 'Flow Solution'

   ENDIF  !FirstTime


      FlowVariable => VariableGet( Model % Variables, FlowSolverName )
      IF ( ASSOCIATED( FlowVariable ) ) THEN
         FlowPerm    => FlowVariable % Perm
         FlowValues  => FlowVariable % Values
      ELSE
         CALL Info('ShapeFactorGravity', &
                    & 'No variable for velocity associated.', Level=4)
      END IF
      NVelo = 0.0
      DO i=1, 2
        Velo(i) = FlowValues(3*(FlowPerm(nodenumber)-1) + i)
        NVelo = NVelo + Velo(i)**2.0
      END DO
      NVelo = SQRT(NVelo)
      
         !------------------------------------
         ! Get K coefficient for that nodes
         !------------------------------------
         CurElement => Model % CurrentElement
         n = CurElement % Type % NumberOfNodes   
         auxReal(1:n) = GetReal( BodyForce, &
              'Lateral Friction Coefficient', GotIt )
         DO i=1, n
            j = CurElement % NodeIndexes (i)  
            IF (nodenumber == j) EXIT  
         END DO
         Kspring = auxReal(i) 
         IF ((NVelo > 1.0e-6).AND.(ABS(mm-1.0)>1.0e-6))  Kspring = Kspring * Nvelo**(mm-1.0) 
        
         gi = g(axis) - Kspring*velo(axis) 

END FUNCTION LateralFriction



