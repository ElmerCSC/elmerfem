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
! ******************************************************************************
! ******************************************************************************
! *
! *  Authors: Basile de Fleurian
! *  Email:   basile.defleurian@uci.edu; gagliar@lgge.obs.ujf-grenoble.fr
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 02 Jully 2013
!------------------------------------------------------------------------------
! Transfer of the water between two drainage systems,transfer is computed to go 
! from the EPL layer to the sediment one giving positive value of the transfer 
! when water flows in this way and negative in the other;
! This function requires two hydrological solvers (inefficient and efficient)
!-----------------------------------------------------------------------------

FUNCTION EPLToIDS(Model,nodenumber,x) RESULT(Transfer)
  USE DefUtils
  USE types
  USE CoordinateSystems
  USE SolverUtils
  USE ElementDescription

  IMPLICIT NONE
  !-----------------------------------------------
  !    External variables
  !-----------------------------------------------
  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL(KIND=dp) :: Transfer, x
  !----------------------------------------------- 
  !    Local variables
  !----------------------------------------------- 

  TYPE(Element_t),POINTER :: Element
  TYPE(Variable_t), POINTER :: IDSHead, EPLHead, PipingMask
  TYPE(ValueList_t), POINTER :: Constants, Material, BodyForce

  REAL(KIND=dp), POINTER :: IDSValues(:), EPLValues(:), MaskValues(:)
  REAL(KIND=dp), ALLOCATABLE :: g(:,:)
  REAL(KIND=dp), ALLOCATABLE ::  Density(:), Gravity(:), &
       EPLComp(:), EPLPorosity(:), EPLThick(:), EPLStoring(:),&
       IDSComp(:), IDSPorosity(:), IDSThick(:), IDSStoring(:),&
       IDSTrans(:),UpperLimit(:), LeakFact(:)
  REAL(KIND=dp)::WatComp

  INTEGER, POINTER :: IDSPerm(:), EPLPerm(:), MaskPerm(:)
  INTEGER :: material_id, bf_id
  INTEGER :: N, istat, i, k

  LOGICAL :: FirstTime=.TRUE., &
       Found= .FALSE., &
       AllocationsDone= .FALSE., &
       UnFoundFatal=.TRUE.

  SAVE g, &
       Density, &     
       Gravity, & 
       EPLComp, &
       EPLPorosity, &
       EPLThick, &
       EPLStoring, &
       IDSComp, &
       IDSPorosity, &
       IDSThick, &
       IDSStoring, &
       IDSTrans, &
       UpperLimit, &
       LeakFact, &
       FirstTime, &
       Found, &
       AllocationsDone

  IF (FirstTime)THEN
     FirstTime = .FALSE.
     IF ( AllocationsDone ) THEN
        DEALLOCATE(g, & 
             Density, &     
             Gravity, & 
             EPLComp, &
             EPLPorosity, &
             EPLThick, &
             EPLStoring, &
             IDSComp, &
             IDSPorosity, &
             IDSThick, &
             IDSStoring, &
             IDSTrans, &
             UpperLimit, &
             LeakFact)
        
     END IF

     N = Model % MaxElementNodes
     ALLOCATE(g( 3,N ), &
          Density( N ), &     
          Gravity( N ), & 
          EPLComp( N ), &
          EPLPorosity( N ), &
          EPLThick( N ), &
          EPLStoring( N ), &
          IDSComp( N ), &
          IDSPorosity( N ), &
          IDSThick( N ), &
          IDSStoring( N ), &
          IDSTrans( N ), &
          UpperLimit( N ), &
          LeakFact( N ), &
          STAT=istat)

     IF ( istat /= 0 ) THEN
        CALL FATAL(  'USF_WaterTransfer', 'Memory allocation error' )
     ELSE
        WRITE(Message,'(a)') 'Memory allocation done'
        CALL INFO("WaterTransfer",Message,Level=4)
     END IF
     AllocationsDone = .TRUE. 
  END IF


  !Get the Mask
  !------------
  PipingMask => VariableGet( Model % Variables, 'Open EPL',UnFoundFatal=UnFoundFatal)
  MaskPerm    => PipingMask % Perm
  MaskValues  => PipingMask % Values

  IF(MaskValues(MaskPerm(nodenumber)).GE.0.0)THEN
     !EPL is not active, no transfer
     Transfer = 0.0
  ELSE

     !Get the EPL and IDS water Heads
     !-------------------------------------
     EPLHead => VariableGet( Model % Variables, 'EPLHead',UnFoundFatal=UnFoundFatal)
     EPLPerm    => EPLHead % Perm
     EPLValues  => EPLHead % Values

     IDSHead => VariableGet( Model % Variables, 'IDSHead',UnFoundFatal=UnFoundFatal)
     IDSPerm    => IDSHead % Perm
     IDSValues  => IDSHead % Values

     !Get Parameters needed to compute the storing coeficient
     !-------------------------------------------------------
     Element => Model % CurrentElement
     Material => GetMaterial(Element)
     N = GetElementNOFNodes(Element)
     Constants => GetConstants()
     BodyForce => GetBodyForce()

     IF (.NOT.ASSOCIATED(Material)) THEN
        WRITE (Message,'(A)') 'No Material found for boundary element no. '
        CALL FATAL("WaterTransfer",Message)
     ELSE
        material_id = GetMaterialId( Element, Found)
        IF(.NOT.Found) THEN
           WRITE (Message,'(A)') 'No Material ID found for boundary element no. '
           CALL FATAL("WaterTransfer",Message)
        END IF
     END IF

     !Generic parameters
     !------------------
     WatComp = GetConstReal( Constants, &
          'Water Compressibility', Found)
     IF ( .NOT.Found ) THEN
        WatComp = 5.04e-4
     END IF

     Density(1:N) = ListGetReal( Material, 'Water Density',  N, Element % NodeIndexes, Found,&
          UnfoundFatal=UnFoundFatal)
        !Previous default value: Density = 1.0055e-18

     IF ( ASSOCIATED( BodyForce ) ) THEN
        bf_id = GetBodyForceId()
        g = 0.0_dp  
        g(1,1:N) = GetConstReal( Model % BodyForces(bf_id) % Values, 'Flow BodyForce 1',Found)
        g(2,1:N) = GetConstReal( Model % BodyForces(bf_id) % Values, 'Flow BodyForce 2',Found)
        g(3,1:N) = GetConstReal( Model % BodyForces(bf_id) % Values, 'Flow BodyForce 3',Found)
        Gravity(1:N) = SQRT(SUM(g**2.0/N))
     END IF

     !EPL Parameters
     !--------------
     EPLComp(1:N) = listGetReal( Material,'EPL Compressibility', N, Element % NodeIndexes, Found,&
          UnfoundFatal=UnFoundFatal)
        !Previous default value: EPLComp(1:N) = 1.0D-2

     EPLPorosity(1:N) = listGetReal( Material,'EPL Porosity', N, Element % NodeIndexes, Found,&
          UnfoundFatal=UnFoundFatal)
        !Previous default value: EPLPorosity(1:N) = 0.4D00

     EPLThick(1:N) = listGetReal( Material,'EPL Thickness', N, Element % NodeIndexes, Found,&
          UnfoundFatal=UnFoundFatal)
        !Previous default value: EPLThick(1:N) = 10.0D00

     !Inefficient Drainage System Parameters
     !--------------------------------------
     IDSComp(1:N) = listGetReal( Material,'IDS Compressibility', N, Element % NodeIndexes, Found,&
          UnfoundFatal=UnFoundFatal)
        !Previous default value: IDSComp(1:N) = 1.0D-2

     IDSPorosity(1:N) = listGetReal( Material,'IDS Porosity', N, Element % NodeIndexes, Found,&
          UnfoundFatal=UnFoundFatal)
        !Previous default value: IDSPorosity(1:N) = 0.4D00

     IDSThick(1:N) = listGetReal( Material,'IDS Thickness', N, Element % NodeIndexes, Found,&
          UnfoundFatal=UnFoundFatal)
        !Previous default value: IDSThick(1:N) = 10.0D00

     !Computing the Storing coeficient of the EPL
     !-------------------------------------------
     EPLStoring(1:N) = EPLThick(1:N) * &
          Gravity(1:N) * EPLPorosity(1:N) * Density(1:N) * &
          (WatComp + EPLComp(1:N)/EPLPorosity(1:N))

     !Computing the Storing coeficient of the Inefficient Drainage System 
     !-------------------------------------------------------------------
     IDSStoring(1:N) = IDSThick(1:N) * &
          Gravity(1:N) * IDSPorosity(1:N) * Density(1:N) * &
          (WatComp + IDSComp(1:N)/IDSPorosity(1:N))

     !Get Inefficient Drainage System Transmitivity
     !---------------------------------------------
      IDSTrans(1:N) = listGetReal( Material,'IDS Transmitivity', N, Element % NodeIndexes, Found,&
          UnfoundFatal=UnFoundFatal)
        !Previous default value: IDSTrans(1:N) = 1.5D04

     !Get Inefficient Drainage System Upper Limit
     !-------------------------------------------
     UpperLimit(1:n) = ListGetReal(Material,'IDSHead Upper Limit',&
          N, Element % NodeIndexes, Found)
     IF (.NOT. Found) THEN
        WRITE(Message,'(a)') 'No upper limit of solution for element no.' 
        CALL INFO("WaterTransfer", Message, level=10)
     END IF

     !Getting parameters needed by the Transfer Equation
     !--------------------------------------------------
     LeakFact(1:n) = ListGetReal( Material, 'Leakage Factor', N, Element % NodeIndexes, Found,&
          UnfoundFatal=UnFoundFatal)
        !Previous default value: LeakFact(1:N) = 50.0D00

     !Computation of the Transfer
     !---------------------------
     DO i=1,N
        k = Element % NodeIndexes(i)
        IF(nodenumber.EQ. k)THEN
           IF((IDSValues(IDSPerm(k)).GE.UpperLimit(i)).AND. &
                (EPLValues(EPLPerm(k)).GT.IDSValues(IDSPerm(k))))THEN
              !If Inefficient Drainage System Head is at upper limit and EPL Head greater than Inefficient Drainage System Head
              ! Then do nothing
              Transfer = 0.0
           ELSE
              Transfer =  IDSTrans(i) &
                   *(EPLValues(EPLPerm(k)) - IDSValues(IDSPerm(k))) & 
                   /(IDSThick(i) * LeakFact(i))
              IF(Transfer.GT.0.0)THEN               !Positive Transfer is from EPL to IDS
                 Transfer = Transfer * EPLStoring(i)
              ELSEIF(Transfer.LT.0.0)THEN           !Negative Transfer from IDS to EPL
                 Transfer = Transfer * IDSStoring(i)
              END IF
           END IF
           EXIT
        END IF
     END DO
  END IF

END FUNCTION EPLToIDS
