!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This library is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU Lesser General Public
! *  License as published by the Free Software Foundation; either
! *  version 2.1 of the License, or (at your option) any later version.
! *
! *  This library is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! *  Lesser General Public License for more details.
! * 
! *  You should have received a copy of the GNU Lesser General Public
! *  License along with this library (in file ../LGPL-2.1); if not, write 
! *  to the Free Software Foundation, Inc., 51 Franklin Street, 
! *  Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Authors: fabien Gillet-Chaulet, Rupert Gladstone
! *  Email:   fabien.gillet-chaulet@univ-grenoble-alpes.fr
! *           rupertgladstone1972@gmail.com
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: May 2022
! *
! ******************************************************************************/
!--------------------------------------------------------------------------------
!>  Module containing utility routines for the SSA
!--------------------------------------------------------------------------------
MODULE SSAMaterialModels
  
  USE DefUtils
  
  IMPLICIT NONE

  CONTAINS

!--------------------------------------------------------------------------------
!>  Return the effective friction coefficient
!--------------------------------------------------------------------------------
   FUNCTION SSAEffectiveFriction(Element,n,Basis,ub,SEP,PartlyGrounded,h,rho,rhow,sealevel,SlipDer) RESULT(Slip)
   IMPLICIT NONE
   REAL(KIND=dp) :: Slip ! the effective friction coefficient
   TYPE(Element_t), POINTER :: Element ! the current element
   INTEGER :: n ! number of nodes
   REAL(KIND=dp) :: Basis(:) ! basis functions
   REAL(KIND=dp) :: ub  ! the velocity for non-linear friction laws
   LOGICAL :: SEP ! Sub-Element Parametrisation of the friction
   LOGICAL :: PartlyGrounded ! is the GL within the current element?
   REAL(KIND=dp) :: h ! for SEP: the ice thickness at current location
   REAL(KIND=dp) :: rho,rhow,sealevel ! density, sea-water density, sea-level
   REAL(KIND=dp),OPTIONAL :: SlipDer ! dSlip/du=dSlip/dv if ub=(u^2+v^2)^1/2 ! required to compute the Jacobian



   INTEGER            :: iFriction
   INTEGER, PARAMETER :: LINEAR = 1
   INTEGER, PARAMETER :: WEERTMAN = 2
   INTEGER, PARAMETER :: BUDD = 5
   INTEGER, PARAMETER :: REG_COULOMB_GAG = 3 ! Schoof 2005 & Gagliardini 2007
   INTEGER, PARAMETER :: REG_COULOMB_JOU = 4 ! Joughin 2019
   
   TYPE(ValueList_t), POINTER :: Material, Constants, pMaterial
   TYPE(Variable_t), POINTER :: GMSol,BedrockSol,NSol
   INTEGER, POINTER :: NodeIndexes(:)
   CHARACTER(LEN=MAX_NAME_LEN) :: Friction
   REAL(KIND=dp) :: Slip2, gravity, qq, hafq
   REAL(KIND=dp) :: fm,fq,MinN,U0
   REAL(KIND=dp) :: alpha,beta,fB
   INTEGER :: GLnIP

   REAL(KIND=dp),DIMENSION(n) :: NodalBeta, NodalGM, NodalBed, NodalLinVelo,NodalC,NodalN
   REAL(KIND=dp) :: bedrock,Hf,fC,fN,LinVelo
   LOGICAL :: Found
   TYPE(Solver_t), POINTER :: pSolver => NULL()

   SAVE :: GLnIP, GMSol, BedRockSol, pSolver
   
!  Sub - element GL parameterisation
   IF (SEP) THEN
     ! If we have adaptive rule then GLnIP keyword is not given. 
     IF(.NOT. ASSOCIATED(pSolver, CurrentModel % Solver ) ) THEN
       pSolver => CurrentModel % Solver 
       GLnIP=ListGetInteger( pSolver % Values, &
           'GL integration points number',Found)
       IF(GLnIP > 0) THEN
         GMSol => VariableGet( CurrentModel % Variables, 'GroundedMask',UnFoundFatal=.TRUE. )
       END IF
       BedrockSol => VariableGet( CurrentModel % Variables, 'bedrock',UnFoundFatal=.TRUE. )
     END IF
     IF(GLnIP > 0) THEN
       CALL GetLocalSolution( NodalGM,UElement=Element,UVariable=GMSol)
     END IF
     CALL GetLocalSolution( NodalBed,UElement=Element,UVariable= BedrockSol)
   END IF

! Friction law
   Material => GetMaterial(Element)
   NodeIndexes => Element % NodeIndexes

   Friction = ListGetString(Material, 'SSA Friction Law',Found, UnFoundFatal=.TRUE.)
   SELECT CASE(Friction)
     CASE('linear')
       iFriction = LINEAR
       fm = 1.0_dp
     CASE('weertman')
       iFriction = WEERTMAN
     CASE('budd')
       iFriction = BUDD
     CASE('coulomb')
       iFriction = REG_COULOMB_GAG
     CASE('regularized coulomb')
       iFriction = REG_COULOMB_JOU
     CASE DEFAULT
       CALL FATAL("SSAEffectiveFriction",'Friction choice not recognised')
    END SELECT

    ! coefficient for all friction parameterisations
    NodalBeta(1:n) = ListGetReal( &
           Material, 'SSA Friction Parameter', n, NodeIndexes(1:n), Found,&
           UnFoundFatal=.TRUE.)

    ! for nonlinear powers of sliding velocity
    SELECT CASE (iFriction)
    CASE(REG_COULOMB_JOU,REG_COULOMB_GAG,WEERTMAN,BUDD)
      fm = ListGetConstReal( Material, 'SSA Friction Exponent', Found , UnFoundFatal=.TRUE.)
      NodalLinVelo(1:n) = ListGetReal( &
           Material, 'SSA Friction Linear Velocity', n, NodeIndexes(1:n), Found,&
           UnFoundFatal=.TRUE.)
    CASE DEFAULT
    END SELECT

    ! where explicit dependence on effective pressure is present...
    SELECT CASE (iFriction)
    CASE(REG_COULOMB_GAG,BUDD)
      NSol => VariableGet( CurrentModel % Variables, 'Effective Pressure', UnFoundFatal=.TRUE. )
      CALL GetLocalSolution( NodalN,UElement=Element, UVariable=NSol)
      MinN = ListGetConstReal( Material, 'SSA Min Effective Pressure', Found, UnFoundFatal=.TRUE.)
      fN = SUM( NodalN(1:n) * Basis(1:n) )
      ! Effective pressure should be >0 (for the friction law)
      fN = MAX(fN, MinN)
    END SELECT
    
    ! parameters unique to one sliding parameterisation
    SELECT CASE (iFriction)

    CASE(BUDD)
      Constants => GetConstants()
      gravity = ListGetConstReal( Constants, 'Gravity Norm', UnFoundFatal=.TRUE. )
      ! calculate haf from N = rho_i g z*
      qq = ListGetConstReal( Material, 'SSA Haf Exponent', Found, UnFoundFatal=.TRUE.)
      hafq = ( fN / (gravity * rho) ) ** qq
      
    CASE(REG_COULOMB_GAG)
      fq = ListGetConstReal( Material, 'SSA Friction Post-Peak', Found, UnFoundFatal=.TRUE. )
      NodalC = 0.0_dp
      NodalC(1:n) = ListGetReal( &
          Material, 'SSA Friction Maximum Value', n, NodeIndexes(1:n), Found,&
          UnFoundFatal=.TRUE.)
      fC = SUM( NodalC(1:n) * Basis(1:n) )

    CASE(REG_COULOMB_JOU)
      U0 = ListGetConstReal( Material, 'SSA Friction Threshold Velocity', Found, UnFoundFatal=.TRUE.)

    END SELECT

    Beta=SUM(Basis(1:n)*NodalBeta(1:n))

    IF (SEP) THEN
      ! Floating
      IF (PartlyGrounded .OR. GLnIP==0) THEN
        bedrock = SUM( NodalBed(1:n) * Basis(1:n) )
        Hf = rhow * (sealevel-bedrock) / rho
        if (h < Hf) beta = 0._dp
      ELSE IF (ALL(NodalGM(1:n) < 0._dp)) THEN
        beta = 0._dp
      END IF
   END IF

   Slip2=0.0_dp
   IF (iFriction .NE. LINEAR) THEN
     LinVelo = SUM( NodalLinVelo(1:n) * Basis(1:n) )
     IF ((iFriction == WEERTMAN).AND.(fm==1.0_dp)) iFriction=LINEAR
     Slip2=1.0_dp
     IF (ub < LinVelo) then
       ub = LinVelo
       Slip2=0.0_dp
     ENDIF
   END IF

   SELECT CASE (iFriction)

   CASE(LINEAR)
     Slip = beta
     IF (PRESENT(SlipDer)) SlipDer = 0._dp

   CASE(WEERTMAN)
     Slip = beta * ub**(fm-1.0_dp)
     IF (PRESENT(SlipDer)) SlipDer = Slip2*Slip*(fm-1.0_dp)/(ub*ub)

   CASE(BUDD)
     Slip = beta * hafq * ub**(fm-1.0_dp)
     ! TODO:
     !     IF (PRESENT(SlipDer)) SlipDer = 

   CASE(REG_COULOMB_GAG)
     IF (fq.NE.1.0_dp) THEN
       alpha = (fq-1.0_dp)**(fq-1.0_dp) / fq**fq
     ELSE
       alpha = 1.0_dp
     END IF
     fB = alpha * (beta / (fC*fN))**(fq/fm)
     Slip = beta * ub**(fm-1.0_dp) / (1.0_dp + fB * ub**fq)**fm
     IF (PRESENT(SlipDer)) SlipDer  = Slip2 * Slip * ((fm-1.0_dp) / (ub*ub) - &
         fm*fq*fB*ub**(fq-2.0_dp)/(1.0_dp+fB*ub**fq))

   CASE(REG_COULOMB_JOU)
     Slip = beta * ub**(fm-1.0_dp) / (ub + U0)**fm
     IF (PRESENT(SlipDer)) SlipDer = Slip2 * Slip * ((fm-1.0_dp) / (ub*ub) - &
         fm*ub**(-1.0_dp)/(ub+U0))

   END SELECT

  END FUNCTION SSAEffectiveFriction

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Compute the element averaged friction
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION ComputeMeanFriction(Element,n,ElementNodes,STDOFs,NodalU,NodalV,NodalZs,NodalZb,MinH, &
                               NodalDensity,SEP,GLnIP,sealevel,rhow) RESULT(strbasemag)
  REAL(KIND=dp) :: strbasemag
  TYPE(Element_t), POINTER :: Element
  INTEGER :: n
  TYPE(Nodes_t) :: ElementNodes
  INTEGER :: STDOFs
  REAL(KIND=dp) :: NodalU(n),NodalV(n),NodalZs(n),NodalZb(n),NodalDensity(n)
  REAL(KIND=dp) :: MinH
  LOGICAL :: SEP
  INTEGER :: GLnIP
  REAL(KIND=dp) :: sealevel,rhow

  LOGICAL :: PartlyGroundedElement
  TYPE(Variable_t),POINTER :: GMSol
  REAL(KIND=dp) :: NodalGM(n)
  TYPE(GaussIntegrationPoints_t) :: IP
  REAL(KIND=dp) :: Basis(n), detJ
  REAL(KIND=dp) :: h,ub,rho,Velo(2)
  REAL(KIND=dp) :: area,tb
  REAL(KIND=dp) :: Ceff
  LOGICAL :: stat
  INTEGER :: t

  strbasemag=0._dp
  IF (SEP) THEN
     GMSol => VariableGet( CurrentModel % Variables, 'GroundedMask',UnFoundFatal=.TRUE. )
     CALL GetLocalSolution( NodalGM,UElement=Element,UVariable=GMSol)
     PartlyGroundedElement=(ANY(NodalGM(1:n).GE.0._dp).AND.ANY(NodalGM(1:n) < 0._dp))
     IF (PartlyGroundedElement .AND. GLnIP > 0 ) THEN
        IP = GaussPoints( Element , np=GLnIP )
     ELSE
        IP = GaussPoints( Element )
     ENDIF
   ELSE
     IP = GaussPoints( Element )
   ENDIF

   area=0._dp
   tb=0._dp
   DO t=1,IP % n
      stat = ElementInfo( Element, ElementNodes, IP % U(t), IP % V(t), &
           IP % W(t),  detJ, Basis )

     h = SUM( (NodalZs(1:n)-NodalZb(1:n)) * Basis(1:n) )
     h=max(h,MinH)

     rho = SUM( NodalDensity(1:n) * Basis(1:n) )

     Velo = 0.0_dp
     Velo(1) = SUM(NodalU(1:n) * Basis(1:n))
     IF (STDOFs == 2) Velo(2) = SUM(NodalV(1:n) * Basis(1:n))
     ub = SQRT(Velo(1)*Velo(1)+Velo(2)*Velo(2))

     Ceff=SSAEffectiveFriction(Element,n,Basis,ub,SEP,PartlyGroundedElement,h,rho,rhow,sealevel)

     area=area+detJ*IP % s(t)
     tb=tb+Ceff*ub*detJ*IP % s(t)
   END DO

   strbasemag=tb/area

   END FUNCTION ComputeMeanFriction

!--------------------------------------------------------------------------------
!>  Return the effective basal mass balance (to be called separately for each IP in
!>  a partly grounded element)
!--------------------------------------------------------------------------------
   FUNCTION SSAEffectiveBMB(Element,nn,Basis,SEM,BMB,hh,FIPcount,rho,rhow,sealevel,FAF) RESULT(BMBatIP)
     
     IMPLICIT NONE
     
     REAL(KIND=dp)              :: BMBatIP ! the effective basal melt rate at integration point
     
     INTEGER,INTENT(IN)         :: nn ! element number of nodes
     REAL(KIND=dp), INTENT(IN)  :: BMB(:) ! basal mass balance
     LOGICAL, INTENT(IN)        :: SEM ! Sub-Element Parametrisation (requires interpolation of floatation on IPs) 
     REAL(KIND=dp),INTENT(IN)   :: hh ! the ice thickness at current location
     REAL(KIND=dp),INTENT(IN)   :: Basis(:)
     TYPE(Element_t),POINTER,INTENT(IN) :: Element
     
     ! optional arguments, depending on melt param
     INTEGER,INTENT(INOUT),OPTIONAL    :: FIPcount
     REAL(KIND=dp),INTENT(IN),OPTIONAL :: rho,rhow,sealevel ! to calculate floatation for SEM3
     REAL(KIND=dp),INTENT(IN),OPTIONAL :: FAF ! Floating area fraction for SEM1
     
     TYPE(ValueList_t), POINTER  :: Material
     TYPE(Variable_t), POINTER   :: GMSol,BedrockSol
     CHARACTER(LEN=MAX_NAME_LEN) :: MeltParam
     
     REAL(KIND=dp),DIMENSION(nn) :: NodalBeta, NodalGM, NodalBed, NodalLinVelo,NodalC
     REAL(KIND=dp) :: bedrock,Hf
     
     LOGICAL :: Found
     
     !  Sub - element GL parameterisation
     IF (SEM) THEN
        GMSol => VariableGet( CurrentModel % Variables, 'GroundedMask',UnFoundFatal=.TRUE. )
        CALL GetLocalSolution( NodalGM,UElement=Element,UVariable=GMSol )
        BedrockSol => VariableGet( CurrentModel % Variables, 'bedrock',UnFoundFatal=.TRUE. )
        CALL GetLocalSolution( NodalBed,UElement=Element,UVariable= BedrockSol )
     END IF
     
     Material => GetMaterial(Element)     
     MeltParam = ListGetString(Material, 'SSA Melt Param',Found, UnFoundFatal=.TRUE.)

     BMBatIP=SUM(Basis(1:nn)*BMB(1:nn))
     
     SELECT CASE(MeltParam)
        
     CASE('FMP','fmp')
        
     CASE('NMP','nmp')
        BMBatIP = 0.0_dp
        
     CASE('SEM1','sem1')
        ! check element type is triangular (would need to modify
        ! CalcFloatingAreaFraction to allow other element types)
        IF (element % type % ElementCode .NE. 303) THEN
           CALL Fatal('SSAEffectiveBMB','Expecting element type 303!')
        END IF

        IF (PRESENT(FAF)) THEN
           BMBatIP = BMBatIP * FAF
        ELSE
           CALL Fatal('SSAEffectiveBMB','FAF (floating area fraction) not present!')
        END IF
        
     CASE('SEM3','sem3')
        bedrock = SUM( NodalBed(1:nn) * Basis(1:nn) )
        Hf= rhow * (sealevel-bedrock) / rho
        IF (hh > Hf) THEN
           BMBatIP = 0.0_dp
        ELSE
           IF (PRESENT(FIPcount)) FIPcount = FIPcount + 1
        END IF
           
     CASE DEFAULT
        WRITE( Message, * ) 'SSA Melt Param not recognised:', MeltParam
        CALL FATAL("SSAEffectiveBMB",Message)
        
     END SELECT
     
   END FUNCTION SSAEffectiveBMB
   
!--------------------------------------------------------------------------------
!> Calculate the fractional floating area of a partly grounded element for SEM1.
!> For implementing SEP1 use (1-FAF) for grounded area fraction.
!> Written for element type 303.
!> To be called per element from ThicknessSolver for SEM1.
!> See Helene Seroussi TC papers from 2014 and 2018.
!> More notes here: https://www.overleaf.com/read/chpfpgzhwvjr
!--------------------------------------------------------------------------------
   FUNCTION CalcFloatingAreaFraction(element,NodalGM,hhVar,sealevel,rho,rhow) RESULT(FAF)

     IMPLICIT NONE
     
     REAL(KIND=dp)              :: FAF ! the area fraction of floating ice

     TYPE(Element_t),POINTER,INTENT(IN)     :: Element
     REAL(KIND=dp),INTENT(IN)               :: NodalGM(:)
     TYPE(Variable_t),POINTER,INTENT(IN)    :: hhVar
     REAL(KIND=dp), INTENT(IN)              :: rho,rhow,sealevel ! to calculate floatation 

     TYPE(Variable_t),POINTER               :: bedVar
     INTEGER,POINTER                        :: hhPerm(:),bedPerm(:)
     REAL(KIND=dp),POINTER                  :: hh(:),bed(:)
     
     ! The following real vars use Yu's terminology (see overleaf linked above)
     REAL(KIND=dp) :: ss,tt,A1,A2,x1,x2,x3,x4,x5,y1,y2,y3,y4,y5
     REAL(KIND=dp) :: B1,B2,B3,Zb1,Zb2,Zb3

     ! NI refers to Node Index.
     ! NI2 is the sole node of its category (floating or GL, see Yu's notes);
     ! NI1 and NI3 are both in the other category.
     ! nn is the number of nodes for the current element (3 for triangles, which is assumed here)
     INTEGER :: NumFLoatingNodes, nn, NI1, NI2, NI3

     bedVar => VariableGet( CurrentModel % Variables, 'bedrock', UnFoundFatal=.TRUE. )
     bedPerm => bedVar % Perm
     bed     => bedVar % Values

     hhPerm  => hhVar % Perm
     hh      => hhVar % Values
     
     nn = GetElementNOFNodes()
     NumFLoatingNodes = -SUM(NodalGM(1:nn))

     IF (ANY(NodalGM(1:nn) > 0._dp)) THEN
        CALL Fatal('CalcFloatingAreaFraction','Fully grounded nodes found!')
     END IF

     IF (NumFLoatingNodes < 1) THEN
        CALL Fatal('CalcFloatingAreaFraction','Not enough floating nodes!')

     ELSEIF (NumFLoatingNodes == 1) THEN
        IF (NodalGM(1) == -1) THEN
           NI2 = element % NodeIndexes(1)
           NI1 = element % NodeIndexes(2)
           NI3 = element % NodeIndexes(3)
        ELSEIF (NodalGM(2) == -1) THEN
           NI2 = element % NodeIndexes(2)
           NI1 = element % NodeIndexes(1)
           NI3 = element % NodeIndexes(3)
        ELSEIF (NodalGM(3) == -1) THEN
           NI2 = element % NodeIndexes(3)
           NI1 = element % NodeIndexes(1)
           NI3 = element % NodeIndexes(2)
        END IF

     ELSEIF (NumFLoatingNodes == 2) THEN
        IF (NodalGM(1) == 0) THEN
           NI2 = element % NodeIndexes(1)
           NI1 = element % NodeIndexes(2)
           NI3 = element % NodeIndexes(3)
        ELSEIF (NodalGM(2) == 0) THEN
           NI2 = element % NodeIndexes(2)
           NI1 = element % NodeIndexes(1)
           NI3 = element % NodeIndexes(3)
        ELSEIF (NodalGM(3) == 0) THEN
           NI2 = element % NodeIndexes(3)
           NI1 = element % NodeIndexes(1)
           NI3 = element % NodeIndexes(2)
        END IF

     ELSEIF (NumFLoatingNodes > 2) THEN
        CALL Fatal('CalcFloatingAreaFraction','Too many floating nodes!')

     END IF

     x1 = CurrentModel % Mesh % Nodes % x(NI1)
     x2 = CurrentModel % Mesh % Nodes % x(NI2)
     x3 = CurrentModel % Mesh % Nodes % x(NI3)

     y1 = CurrentModel % Mesh % Nodes % y(NI1)
     y2 = CurrentModel % Mesh % Nodes % y(NI2)
     y3 = CurrentModel % Mesh % Nodes % y(NI3)

     B1 = bed(bedPerm(NI1))
     B2 = bed(bedPerm(NI2))
     B3 = bed(bedPerm(NI3))

     Zb1= sealevel - hh(hhPerm(NI1)) * rho/rhow
     Zb2= sealevel - hh(hhPerm(NI2)) * rho/rhow
     Zb3= sealevel - hh(hhPerm(NI3)) * rho/rhow
     
     ss = (Zb2-B2)/(B3-B2-Zb3+Zb2)
     tt = (Zb1-B1)/(B2-B1-Zb2+Zb1)

     x4 = x1 + tt*(x2-x1)
     x5 = x2 + ss*(x3-x2)
     y4 = y1 + tt*(y2-y1)
     y5 = y2 + ss*(y3-y2)
     
     A1 = 0.5_dp * ABS( x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2) )
     A2 = 0.5_dp * ABS( x4*(y5-y2) + x5*(y2-y4) + x2*(y4-y5) )

     FAF = A2/A1
     
     IF (NumFLoatingNodes == 2) FAF = 1.0_dp - FAF  ! Needed because FAF was grounded fraction

   END FUNCTION CalcFloatingAreaFraction
   
END MODULE SSAMaterialModels
 
