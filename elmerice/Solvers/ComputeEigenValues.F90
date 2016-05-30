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
! *  Authors: Fabien / OG 
! *  Email:   
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 13/10/05 from version 1.5
! * 
! *****************************************************************************/
!------------------------------------------------------------------------------
   RECURSIVE SUBROUTINE ComputeEigenValues( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------

    USE DefUtils

    IMPLICIT NONE

!------------------------------------------------------------------------------
!******************************************************************************
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
!     INPUT: Timestep size for time dependent simulations (NOTE: Not used
!            currently)
!
!******************************************************************************

     TYPE(Model_t)  :: Model
     TYPE(Solver_t), TARGET :: Solver

     LOGICAL ::  TransientSimulation
     REAL(KIND=dp) :: dt
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------


    
     INTEGER :: dim, StressDOFs
     TYPE(Variable_t), POINTER :: StressSol,EigenSol,EigenV1,EigenV2,EigenV3
     REAL(KIND=dp), POINTER ::  Stress(:),Eigen(:)
     INTEGER, POINTER :: StressPerm(:),EigenPerm(:)

     REAL(KIND=dp),dimension(3,3) :: localM,EigenVec
     REAL(KIND=dp),dimension(3) :: EigValues, EI, ki
     REAL(KIND=dp) :: WORK(24),Dumy(1)
     Real(KIND=dp) :: a 
     INTEGER :: i, j, t, ordre(3),infor
     LOGICAL :: GotIt,UnFoundFatal=.TRUE.
     CHARACTER(LEN=MAX_NAME_LEN) :: TensorVarName, EigenVarName


!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
      dim = CoordinateSystemDimension()

      TensorVarName = GetString( Solver % Values, 'Tensor Variable Name', GotIt )    
      IF (.NOT.Gotit) TensorVarName = 'Stress'
      StressSol => VariableGet( Solver % Mesh % Variables, TensorVarName )
      StressPerm => StressSol % Perm
      StressDOFs = StressSol % DOFs
      Stress => StressSol % Values

     ! Eigen Values (dimension 3)
      EigenVarName = GetString( Solver % Values, 'EigenValue Variable Name', GotIt )    
      IF (.NOT.Gotit) EigenVarName = 'EigenStress'
      EigenSol => VariableGet( Solver % Mesh % Variables, EigenVarName,UnFoundFatal=UnFoundFatal)
      EigenPerm => EigenSol % Perm
      Eigen => EigenSol % Values

     ! Eigen Vector (dimension 3)
      EigenV1 => VariableGet( Solver % Mesh % Variables, 'EigenVector1' )
      EigenV2 => VariableGet( Solver % Mesh % Variables, 'EigenVector2' )
      EigenV3 => VariableGet( Solver % Mesh % Variables, 'EigenVector3' )


     Do t=1,Solver % Mesh % NumberOfNodes

          ! the Stress components [Sxx, Syy, Szz, Sxy, Syz, Szx] 
          localM=0.0_dp
          localM(1,1)=Stress( StressDOFs * (StressPerm(t) - 1 ) + 1)
          localM(2,2)=Stress( StressDOFs * (StressPerm(t) - 1 ) + 2)
          localM(3,3)=Stress( StressDOFs * (StressPerm(t) - 1 ) + 3)
          localM(1,2)=Stress( StressDOFs * (StressPerm(t) - 1 ) + 4)
          localM(2,1)=localM(1,2)
          if (dim.eq.3) then
                  localM(2,3)=Stress( StressDOFs * (StressPerm(t) - 1 ) + 5)
                  localM(1,3)=Stress( StressDOFs * (StressPerm(t) - 1 ) + 6)
                  localM(3,2)=localM(2,3)
                  localM(3,1)=localM(1,3)
           end if

! Compute EigenValues using lapack DGEEV subroutine
           CALL DGEEV('N','V',3,localM,3,EigValues,EI,Dumy,1,EigenVec,3,Work,24,infor )
           IF (infor.ne.0) &
               CALL FATAL('Compute EigenValues', 'Failed to compute EigenValues')
           localM(1:3,1:3)=EigenVec(1:3,1:3)
! Ordered value Ev1 < Ev2 < Ev3
         
           Do i=1,3
             ki(i)=EigValues(i)
             ordre(i)=i
           End Do
           Do j=2,3
             a=ki(j)
             Do i=j-1,1,-1
               If (ki(i).LE.a) Goto 20
               ki(i+1)=ki(i)
               ordre(i+1)=ordre(i)
             End Do
  20         Continue
             ki(i+1)=a
             ordre(i+1)=j
           End Do
                       
           Eigen( 3 * ( EigenPerm(t) - 1) + 1) = EigValues(ordre(1))
           Eigen( 3 * ( EigenPerm(t) - 1) + 2) = EigValues(ordre(2))
           Eigen( 3 * ( EigenPerm(t) - 1) + 3) = EigValues(ordre(3))

           If (associated(EigenV1)) then
              EigenV1 % Values( 3 * (EigenV1 % Perm(t) - 1) + 1 ) = localM(ordre(1),1)
              EigenV1 % Values( 3 * (EigenV1 % Perm(t) - 1) + 2 ) = localM(ordre(1),2)
              EigenV1 % Values( 3 * (EigenV1 % Perm(t) - 1) + 3 ) = localM(ordre(1),3)
           endif
           If (associated(EigenV2)) then
              EigenV2 % Values( 3 * (EigenV2 % Perm(t) - 1) + 1 ) = localM(ordre(2),1)
              EigenV2 % Values( 3 * (EigenV2 % Perm(t) - 1) + 2 ) = localM(ordre(2),2)
              EigenV2 % Values( 3 * (EigenV2 % Perm(t) - 1) + 3 ) = localM(ordre(2),3)
           endif
           If (associated(EigenV3)) then
              EigenV3 % Values( 3 * (EigenV3 % Perm(t) - 1) + 1 ) = localM(ordre(3),1)
              EigenV3 % Values( 3 * (EigenV3 % Perm(t) - 1) + 2 ) = localM(ordre(3),2)
              EigenV3 % Values( 3 * (EigenV3 % Perm(t) - 1) + 3 ) = localM(ordre(3),3)
           Endif

     End do

!------------------------------------------------------------------------------
      END SUBROUTINE ComputeEigenValues
!------------------------------------------------------------------------------

