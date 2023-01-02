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
! * This module solves electrochemical reactions for two phase (solid+electrolyte)
! * battery. The model uses the Galerkin method for all equations. The robustness
! * of the code still could be improved. The continuous model itself follows very
! * closely the references while the way how the discrete system is solved is almost
! * completely different.
! *
! * If you're interested in continuing the work we would appreciate contacting
! * the authors (Peter & Timo). 
! * 
! *  References:
! *  [1] Ashwin S. Borakhadikar, "One Dimensional Computer Modeling of a Lithium-Ion
! *      Battery". Master's Thesis 2017.
! *  [2] Lu Xia et al. "A Computationally Efficient Implementation of an
! *      Electrochemistry-Based Model for Lithium-Ion Batteries", IFAC Proceedings,
! *      2017.
! *
! ******************************************************************************
! *
! *  Authors: Peter Råback, CSC
! *           Timo Uimonen, University of Vaasa, timo.uimonen@gmail.com
! *  Email:   elmeradm@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 12 March 2019
! *
! *****************************************************************************/


!-------------------------------------------------------------------------------
!> This module includes some common functinalities needed by the four different
!> PDE solvers that constitute the battery model.
!-------------------------------------------------------------------------------
MODULE BatteryModule

  USE DefUtils
  IMPLICIT NONE

  LOGICAL :: BatteryInitialized = .FALSE.
  LOGICAL :: UseTimeAveFlux, UseTimeAveDiff, UseMeanFlux 
  REAL(KIND=dp) :: FaradayConstant, GasConstant
  TYPE(Mesh_t), POINTER :: MainMesh, SubMesh
  TYPE(Variable_t), POINTER :: PhisVar, PhieVar, CsVar, CeVar, EtaVar,JliVar,&
      SOCVar, CsAveVar, CsInitVar, CsErrVar, Cs1DVar, CsDiffVar, JliIntegVar, &
      CsFullVar => NULL()
  REAL(KIND=dp), ALLOCATABLE :: AnodeWeight(:),MassMult1D(:)
  REAL(KIND=dp) :: AnodeWeightSum, CathodeWeightSum
  REAL(KIND=dp), POINTER :: Jli0(:) => NULL()
  

CONTAINS 

  ! Initializes some data for the battery simulation.
  ! This may be called by any solver who gets there first - except the 1d solver for Cs.
  !-------------------------------------------------------------------------------------
  SUBROUTINE InitializeBattery()
    CHARACTER(*), PARAMETER :: Caller = 'InitializeBattery'
    TYPE( Mesh_t), POINTER :: Mesh
    LOGICAL :: Found
    INTEGER :: i
    
    IF( BatteryInitialized ) RETURN
    
    IF( ListCheckPresent( CurrentModel % Constants,'Universal Gas Constant') ) THEN
      CALL Fatal(Caller,'Use keyword "Gas Constant", no "Universal"')
    END IF

    GasConstant = ListGetConstReal( CurrentModel % Constants,'Gas Constant',Found)
    IF(.NOT. Found ) GasConstant = 8.31446261815324
    
    FaradayConstant = ListGetConstReal( CurrentModel % Constants,'Faraday Constant',Found)
    IF(.NOT. Found ) FaradayConstant = 96485.33212 
    
    UseTimeAveFlux = ListGetLogicalAnySolver( CurrentModel,'Use Time Average Flux')
    IF( UseTimeAveFlux ) THEN
      CALL Info(Caller,'Using time-averaged Jli flux',Level=7)
    END IF

    UseTimeAveDiff = ListGetLogicalAnySolver( CurrentModel,'Use Time Average Diffusion')
    IF( UseTimeAveDiff ) THEN
      CALL Info(Caller,'Using time-averaged Ce diffusion term',Level=7)
    END IF
    
    ! Possible field for previous flux
    UseMeanFlux = ListGetLogicalAnySolver( CurrentModel,'Use Mean Flux')
    IF( UseMeanFlux ) THEN
      CALL Info(Caller,'Using iteration-averaged relaxed Jli flux',Level=7)
    END IF
          
    MainMesh => GetMesh()
    
    DO i=1,CurrentModel % NumberOfSolvers
      IF( ListGetString( CurrentModel % Solvers(i) % Values,'Variable',Found ) == 'cs onedim' ) THEN
        SubMesh => CurrentModel % Solvers(i) % Mesh
        EXIT
      END IF
    END DO
    IF(.NOT. ASSOCIATED( SubMesh ) ) THEN
      CALL Fatal(Caller,'Could not define the 1D submesh!')
    END IF
    
    Mesh => MainMesh
    
    CALL SetVariablePointers()    
    CALL CalculateAnodeWeight()

    ! Using some "inertia" in the flux helps in convergence.
    ! We may either use natural inertia coming from time,
    ! or just iteration coming from inertia. 
    IF( UseTimeAveFlux ) THEN      
      JLi0 => JliVar % PrevValues(:,1)    
    ELSE IF( UseMeanFlux ) THEN
      ALLOCATE( Jli0( SIZE( JliVar % Values ) ) )
      Jli0 = 0.0_dp
    END IF
    
    BatteryInitialized = .TRUE.
    

  CONTAINS

    SUBROUTINE SetVariablePointers()
      PhisVar => VariableGet( Mesh % Variables,'Phis')
      IF(.NOT. ASSOCIATED( PhisVar ) ) THEN
        CALL Fatal(Caller,'Variable "Phis" not present!')
      END IF
      PhieVar => VariableGet( Mesh % Variables,'Phie')
      IF(.NOT. ASSOCIATED( PhieVar ) ) THEN
        CALL Fatal(Caller,'Variable "Phie" not present!')
      END IF
      CsVar => VariableGet( Mesh % Variables,'Cs')
      IF(.NOT. ASSOCIATED( CsVar ) ) THEN
        CALL Fatal(Caller,'Variable "Cs" not present!')
      END IF
      CeVar => VariableGet( Mesh % Variables,'Ce')
      IF(.NOT. ASSOCIATED( CeVar ) ) THEN
        CALL Fatal(Caller,'Variable "Ce" not present!')
      END IF
      IF( UseTimeAveDiff ) THEN
        CeVar % PrevValues(:,1) = CeVar % Values 
      END IF

      EtaVar => VariableGet( Mesh % Variables,'Eta')
      IF(.NOT. ASSOCIATED( EtaVar ) ) THEN
        CALL Fatal(Caller,'Variable "Eta" not present!')
      END IF
      JliVar => VariableGet( Mesh % Variables,'Jli')

      ! Enforce history for JliVar
      IF( .NOT. ASSOCIATED( JliVar % PrevValues ) ) THEN 
        ALLOCATE( JliVar % PrevValues(SIZE(JliVar % Values), 1 ) )
      END IF
      JliVar % PrevValues(:,1) = JliVar % Values 
      
      IF(.NOT. ASSOCIATED( JliVar ) ) THEN
        CALL Fatal(Caller,'Variable "Jli" not present!')
      END IF
      CsAveVar => VariableGet( Mesh % Variables,'Cs Ave')
      IF(.NOT. ASSOCIATED( CsAveVar ) ) THEN
        CALL Fatal(Caller,'Variable "Cs Ave" not present!')
      END IF
      SOCVar => VariableGet( Mesh % Variables,'SOC')
      IF(.NOT. ASSOCIATED( SOCVar ) ) THEN
        CALL Fatal(Caller,'Variable "SOC" not present!')
      END IF
      Cs1dVar => VariableGet( SubMesh % Variables,'Cs OneDim')
      IF(.NOT. ASSOCIATED( SOCVar ) ) THEN
        CALL Fatal(Caller,'Variable "Cs Onedim" not present!')
      END IF
      ALLOCATE( MassMult1D(SIZE(Cs1dVar % Values) ) )
      MassMult1D = 0.0_dp
     
      CsDiffVar => VariableGet( Mesh % Variables,'Cs Diff')
      JliIntegVar => VariableGet( Mesh % Variables,'Jli Integral')            
      CSInitVar => VariableGet( Mesh % Variables,'Cs Init')
      CSErrVar => VariableGet( Mesh % Variables,'Cs Err')
      
    END SUBROUTINE SetVariablePointers
      

    ! Calculates the vector that includes nodal weights w_i that can be used
    ! in integrating or averaging over field variables. 
    ! For anode w_i<0 and cathode w_i>0.
    !---------------------------------------------------------------------
    SUBROUTINE CalculateAnodeWeight()

      REAL(KIND=dp) :: Basis(8), DetJ, Weight
      TYPE(Nodes_t), SAVE :: Nodes      
      TYPE(GaussIntegrationPoints_t) :: IP
      TYPE(ValueList_t), POINTER :: Material, PrevMaterial
      TYPE(Element_t), POINTER :: Element
      INTEGER :: t, n, elem
      LOGICAL :: Anode, Stat, AnodeOrCathode           
     
      ALLOCATE( AnodeWeight( SIZE( JliVar % Values ) ) )
      AnodeWeight = 0.0_dp

      PrevMaterial => NULL()

      DO elem=1,Mesh % NumberOfBulkElements
        Element => Mesh % Elements(elem)

        AnodeOrCathode = ( ALL( JliVar % Perm(Element % NodeIndexes) /= 0) )
        IF(.NOT. AnodeOrCathode ) CYCLE

        Material => GetMaterial(Element)
        IF( .NOT. ASSOCIATED( Material, PrevMaterial ) ) THEN
          Anode = ListGetLogical(Material,'Anode',Found)
          PrevMaterial => Material
        END IF
          
        n = GetElementNOFNodes( Element )        
        CALL GetElementNodes( Nodes, UElement=Element, UMesh=Mesh )
        
        IP = GaussPoints( Element )

        DO t=1,IP % n
          stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis )
          Weight = IP % s(t) * DetJ
          
          IF( .NOT. Anode ) Weight = -Weight 
          AnodeWeight( JliVar % Perm( Element % NodeIndexes ) ) = &
              AnodeWeight( JliVar % Perm( Element % NodeIndexes ) ) + Weight * Basis(1:n)
        END DO
      END DO

      AnodeWeightSum = SUM( AnodeWeight, AnodeWeight > 0 )
      CathodeWeightSum = SUM( -AnodeWeight, AnodeWeight < 0 )

    END SUBROUTINE CalculateAnodeWeight

    
  END SUBROUTINE InitializeBattery


  ! Print on output information on variables.
  ! This may be given priority to save time.
  ! There are two modes
  ! 1) Min/max + plain number average
  ! 2) Min/Max + weighted integral and average for Anode/Cathode
  !--------------------------------------------------------------------
  SUBROUTINE VariableRange( Var, Level, AWeight ) 
    TYPE(Variable_t), POINTER :: Var
    INTEGER :: Level
    REAL(KIND=dp), OPTIONAL :: AWeight(:)
    REAL(KIND=dp) :: RangeVals(8)

    IF( .NOT. InfoActive( Level ) ) RETURN

    IF( PRESENT( AWeight ) ) THEN
      RangeVals(1) = MINVAL( Var % Values * AWeight, AWeight > 0 )
      RangeVals(2) = MAxVAL( Var % Values * AWeight, AWeight > 0 )
      RangeVals(3) = SUM( Var % Values * AWeight, AWeight > 0 )
      RangeVals(4) = RangeVals(3) / SUM( AWeight, AWeight > 0 )
      PRINT *,'Anode range '//TRIM(Var % Name)//' (min,max,int,ave):',RangeVals(1:4)
      
      RangeVals(5) = MINVAL( Var % Values * AWeight, AWeight < 0 )
      RangeVals(6) = MAXVAL( Var % Values * AWeight, AWeight < 0 )
      RangeVals(7) = -SUM( Var % Values * AWeight, AWeight < 0 )
      RangeVals(8) = -RangeVals(7) / SUM( AWeight, AWeight < 0 )
      PRINT *,'Cathode range for '//TRIM(Var % Name)//' (int, ave):',RangeVals(5:8)
    ELSE
      RangeVals(1) = MINVAL(Var % Values)
      RangeVals(2) = MAXVAL(Var % Values)
      RangeVals(3) = SUM( Var % Values) / SIZE( Var % Values )
      PRINT *,'Range for '//TRIM(Var % Name)//' (min, max, ave):',RangeVals(1:3)
    END IF
    
  END SUBROUTINE VariableRange


  ! Get extreme left node from 1D mesh using coordinate
  ! This is used to set BC condition for 1D solid phase mesh.
  !--------------------------------------------------------
  FUNCTION ExtremeLeftNode( Mesh ) RESULT ( node ) 
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER :: node
    !--------------------------------------------------------
    REAL(KIND=dp), POINTER :: x(:)    
    INTEGER i    
    x => Mesh % Nodes % x
    node = 1
    DO i = 2, Mesh % NumberOfNodes
      IF( x(i) < x(node) ) node = i
    END DO
  END FUNCTION ExtremeLeftNode

  ! Get extreme right node from 1D mesh using coordinate
  ! This is used to set BC condition for 1D solid phase mesh.
  !--------------------------------------------------------
  FUNCTION ExtremeRightNode( Mesh ) RESULT ( node ) 
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER :: node
    !--------------------------------------------------------
    REAL(KIND=dp), POINTER :: x(:)    
    INTEGER i    
    x => Mesh % Nodes % x
    node = 1
    DO i = 2, Mesh % NumberOfNodes
      IF( x(i) > x(node) ) node = i
    END DO
  END FUNCTION ExtremeRightNode

  ! The relation for the equilibrium potential for the negative electrode     
  ! Eq. (3.20) in [1]
  !-------------------------------------------------------------------------
  FUNCTION Uneg(x) RESULT ( Un ) 
    REAL(KIND=dp) :: x, Un
    REAL(KIND=dp) :: sqrtx

    sqrtx = SQRT(x)  

    Un = 8.0029 + 5.0647*x - 12.578*sqrtx - 8.6322e-4/x &
        + 2.1765e-5*x*sqrtx - 0.46016*EXP(15*(0.06-x)) - 0.55364*EXP(-2.4326*(x-0.92))

  END FUNCTION Uneg

  ! The relation for the equilibrium potential for the positive electrode     
  ! Eq. (3.21) in [1]
  !-------------------------------------------------------------------------
  FUNCTION Upos(y) RESULT ( Up ) 
    REAL(KIND=dp) :: y, y2, y4, y6, Up

    ! We use these to make the computation slightly faster, I hope
    y2 = y*y
    y4 = y2*y2
    y6 = y2*y4
    
    Up = 85.6781*y6 - 357.7*y*y4 + 613.89*y4 - 555.65*y*y2 + 281.06*y2 &
        - 76.648*y - 0.30987*EXP(5.657*y**115) + 13.1983
  END FUNCTION Upos
  
  
  ! Compute ButlerVolmer equation for one node.
  ! Optionally also compute sensitivity to given variations using numerical
  ! or symbolid differentiation. Numerical differentiation uses recursion. 
  !----------------------------------------------------------------------------
  RECURSIVE FUNCTION ButlerVolmer(Material, node, Phis, Phie, Cs, Ce, &
      Eta, EtaFixed, djdPhis, djdPhie, djdCs, djdCe ) RESULT ( j_Li )
    IMPLICIT NONE
    TYPE( ValueList_t), POINTER :: Material
    INTEGER :: node
    REAL(KIND=dp) :: Phis, Phie, Cs, Ce, Eta, j_Li, dx, TimeStep
    LOGICAL :: EtaFixed
    REAL(KIND=dp), OPTIONAL ::  djdPhis, djdPhie, djdCs, djdCe

    LOGICAL :: Found
    LOGICAL, SAVE :: cathode, anode, Visited = .FALSE.
    REAL(KIND=dp), SAVE :: Cs_max, alpha_c, alpha_a, beta_c, beta_a, &
        a_s, k_0, T, eps, PotScale
    TYPE( ValueList_t), POINTER, SAVE :: PrevMaterial => NULL()
    REAL(KIND=dp) :: U,c_0,f1,f2,Eta0
    LOGICAL, SAVE :: NewMaterial, FirstTime = .TRUE.

    IF(.NOT. Visited ) THEN
      eps = 1.0d-3
      Visited = .TRUE.
    END IF

    NewMaterial = .NOT. ASSOCIATED( Material, PrevMaterial ) 

    ! Fetching values from keywords takes some time so only do it
    ! if the material has changed. 
    IF( NewMaterial ) THEN
      T = ListGetCReal( CurrentModel % Constants, 'Ambient Temperature')

      PotScale = ListGetCReal( CurrentModel % Solver % Values,'Overpotential Scaling',Found )
      IF(.NOT. Found ) PotScale = 1.0_dp
            
      Anode = ListGetLogical(Material,'Anode',Found)
      Cathode = ListGetLogical(Material,'Cathode',Found)

      Cs_max = ListGetCReal( Material,'Maximum solid phase concentration')
      alpha_a = ListGetCReal( Material,'Anodic charge transfer coefficient',Found)
      IF(.NOT. Found ) alpha_a = 0.5_dp
      alpha_c = ListGetCReal( Material,'Cathodic charge transfer coefficient',Found)
      IF(.NOT. Found ) alpha_c = 0.5_dp
      k_0 = ListGetCReal( Material,'Kinetic constant')

      ! Eq. (3.4) in [1]   
      a_s = InterfaceSurfaceArea( Material )
      
      ! Lets combine some constants to save time
      beta_a = alpha_a * FaradayConstant / (GasConstant*T) 
      beta_c = alpha_c * FaradayConstant / (GasConstant*T)
      PrevMaterial => Material

      Timestep = GetTimeStep()
      FirstTime = ( TimeStep == 1 ) 
    END IF

    ! IF( FirstTime )
    ! we could do here some exceptions 
    ! END IF
    
    ! Eq. (3.22) in [1] 
    IF( .NOT. EtaFixed ) THEN    
      IF( Anode ) THEN
        U = Uneg( Cs / Cs_max )
      ELSE IF( Cathode ) THEN
        U = Upos( Cs / Cs_max )
      END IF
      
      ! Eq. (3.19) in [1]
      eta = Phis - Phie - U
    END IF
      
    IF( NewMaterial ) THEN
      !PRINT *,'eta:',Anode,eta,Phis,Phie,U
    END IF
      
    ! Eq. (3.18) in [1]
    c_0 = a_s * k_0 * Ce**alpha_a * (Cs_max - Cs)**alpha_a * Cs**alpha_c
    
    ! Eq. (3.17) in [1] with R_SEI set to zero. 
    j_Li = c_0 * ( EXP( eta * beta_a ) - EXP( -eta * beta_c ) )
    
    ! Optional numerical derivatives in case we want to make
    ! Newton linearization with respect to some variable.

    ! This basically tests for NaN's
    IF( J_li /= J_li ) THEN
      PRINT *,'jli',a_s,c_0,k_0,Cs_max - Cs,Ce,Cs,j_li
      STOP
    END IF
      
    ! For potentials we use analytical derivatives and chain rule
    IF( PRESENT(djdPhis) ) THEN
      djdPhis = c_0 * ( beta_a * EXP( eta * beta_a ) + beta_c * EXP( -eta * beta_c ) )
    END IF

    IF( PRESENT(djdPhie) ) THEN
      djdPhie = -c_0 * ( beta_a * EXP( eta * beta_a ) + beta_c * EXP( -eta * beta_c ) )
    END IF

    ! For concentrations we use numerical central differences
    IF( PRESENT(djdCs) ) THEN
      ! Does not include derivative of U
      IF( .FALSE.) THEN
        djdCs = j_Li * ( alpha_c / Cs - alpha_a / (Cs_max - Cs ) )
      ELSE
        dx = eps*Cs
        Eta0 = Eta
        f1 = ButlerVolmer(Material, node, Phis, Phie, Cs+dx, Ce, Eta, EtaFixed)
        f2 = ButlerVolmer(Material, node, Phis, Phie, Cs-dx, Ce, Eta, EtaFixed)
        Eta = Eta0
        djdCs = (f1-f2)/(2*dx)
      END IF
    END IF

    IF( PRESENT(djdCe) ) THEN
      ! Does not include derivative of U
      IF( .FALSE. ) THEN
        djdCe = j_Li * ( alpha_a / Ce )
      ELSE
        dx = eps*Ce
        Eta0 = Eta
        f1 = ButlerVolmer(Material, node, Phis, Phie, Cs, Ce+dx, Eta, EtaFixed)
        f2 = ButlerVolmer(Material, node, Phis, Phie, Cs, Ce-dx, Eta, EtaFixed)
        Eta = Eta0
        djdCe = (f1-f2)/(2*dx)
      END IF
    END IF

  END FUNCTION ButlerVolmer


  ! Given homogeneous flux J_Li convert it to flux for the solid phase spheres
  ! The idea is that we can hide scaling of Jli such that what remains is
  ! rather standard PDE and BCs for (3.2) for Cs.
  !----------------------------------------------------------------------------
  FUNCTION SolidFluxScaling( Material ) RESULT ( FluxScale ) 
    TYPE( ValueList_t), POINTER :: Material
    REAL(KIND=dp) :: FluxScale

    REAL(KIND=dp), SAVE :: R_s, Eps_s, FluxScale0
    TYPE( ValueList_t), POINTER, SAVE :: PrevMaterial => NULL()
    LOGICAL :: Found
        
    IF(.NOT. ASSOCIATED( Material, PrevMaterial ) ) THEN
      R_s = ListGetCReal( Material,'Particle Radius')
      Eps_s = ListGetCReal( Material,'Active particle volume fraction')       
      PrevMaterial => Material

      ! Scaling coefficient 1/(a_s*F) as defined by Eqs. (3.4) and (3.5) in [1]
      ! Note that we further scale by R_s**2 since in the weak form the BC is multiplied by area
      ! and divide by R_s**3 since the mesh is a unit mesh. 
      FluxScale0 = 1.0_dp / ( 3 * Eps_s * FaradayConstant) 
    END IF

    FluxScale = FluxScale0

  END FUNCTION SolidFluxScaling


  FUNCTION InterfaceSurfaceArea( Material ) RESULT ( A_s )
    TYPE( ValueList_t), POINTER :: Material
    REAL(KIND=dp) :: A_s

    REAL(KIND=dp), SAVE :: R_s, Eps_s
    TYPE( ValueList_t), POINTER, SAVE :: PrevMaterial => NULL()

    IF(.NOT. ASSOCIATED( Material, PrevMaterial ) ) THEN
      R_s = ListGetCReal( Material,'Particle Radius')
      Eps_s = ListGetCReal( Material,'Active particle volume fraction')
      PrevMaterial => Material
    END IF

    ! Eq. (3.5) in [1]
    A_s = 3 * Eps_s / R_s

  END FUNCTION InterfaceSurfaceArea
 

  ! Given homogeneous flux J_Li convert it to flux for electrolyte phase
  ! The idea is that we can hide scaling of Jli such that what remains is
  ! rather standard PDE for (3.6) for Cs.
  !-----------------------------------------------------------------------
  FUNCTION ElectrolyteFluxScaling( Material ) RESULT ( FluxScale ) 
    TYPE( ValueList_t), POINTER :: Material
    REAL(KIND=dp) :: FluxScale
    LOGICAL :: Found

    REAL(KIND=dp), SAVE :: t_plus
    LOGICAL, SAVE :: Visited = .FALSE.

    IF(.NOT. Visited ) THEN
      t_plus = ListGetCReal( CurrentModel % Constants,'Transference Number',Found)
      IF(.NOT. Found ) t_plus = 0.363_dp 
      Visited = .TRUE.
    END IF

    ! Eq. (3.6) in [1]
    FluxScale = (1-t_plus) / FaradayConstant

  END FUNCTION ElectrolyteFluxScaling


  FUNCTION EffIonConductivity( Material, Ce ) RESULT ( KeffIon ) 
    TYPE( ValueList_t), POINTER :: Material
    REAL(KIND=dp) :: Ce, KeffIon

    REAL(KIND=dp), SAVE :: Eps_e
    TYPE( ValueList_t), POINTER, SAVE :: PrevMaterial => NULL()

    IF( .NOT. ASSOCIATED( Material, PrevMaterial ) ) THEN
      Eps_e = ListGetCReal( Material,'Electrolyte volume fraction')
      PrevMaterial => Material
    END IF
    
    ! Eq. (3.14) and (3.15) in [1] for electrolyte phase ionic conductivity
    KeffIon = Eps_e * 15.8e-4 * Ce * EXP( 0.85*(1.0e-3 * Ce)**1.4)
    
  END FUNCTION EffIonConductivity


  FUNCTION EffDiffConductivity( Material, Ce ) RESULT ( KeffDiff )   
    TYPE( ValueList_t), POINTER :: Material
    REAL(KIND=dp) :: Ce, KeffDiff
    LOGICAL :: Found

    REAL(KIND=dp), SAVE :: t_plus, T
    REAL(KIND=dp) :: KeffIon
    LOGICAL, SAVE :: Visited = .FALSE.

    IF(.NOT. Visited ) THEN
      T = ListGetConstReal( CurrentModel % Constants,'Ambient Temperature')
      t_plus = ListGetCReal( CurrentModel % Constants,'Transference Number',Found)
      IF(.NOT. Found ) t_plus = 0.363_dp
      Visited = .TRUE.
    END IF

    KeffIon = EffIonConductivity( Material, Ce ) 

    ! Eqs. (3.13) in [1] for effective diffusive conductivity assuming constant activity f
    KeffDiff = KeffIon * ( 2*T*GasConstant / FaradayConstant ) * (t_plus - 1 ) 

  END FUNCTION EffDiffConductivity

  
  ! This updates the flux from Butler-Volmer equation to present the current state
  ! Optionally one may also obtain the sensitivity which may be desirable for
  ! Newton's linearization. Note that this may be called by any of the four solvers
  ! and hence the keywords in the solver section may be used to control behavior
  ! of each solver separately. 
  !--------------------------------------------------------------------------------
  SUBROUTINE  ButlerVolmerUpdate( Solver )
    TYPE(Solver_t), TARGET :: Solver

    CHARACTER(LEN=MAX_NAME_LEN) :: str    
    REAL(KIND=dp) :: Phie, Phis, Cs, Ce, Flux, Eta, FluxDer, Relax, EtaRelax    
    TYPE(Variable_t), POINTER :: SensVar
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Element_t), POINTER :: Element
    TYPE( ValueList_t), POINTER :: Material
    LOGICAL, ALLOCATABLE :: NodeDone(:)
    INTEGER :: t, i, j, k, n, vari
    LOGICAL :: Found, EtaFixed
    TYPE(Solver_t), POINTER :: PSolver
    TYPE(ValueList_t), POINTER :: Params
    LOGICAL :: CalcSens(4)
    TYPE(VariableTable_t) :: SensVars(4)
    REAL(KIND=dp), POINTER :: TmpVals(:)
    INTEGER :: timestep, prevtimestep = 0

    SAVE prevtimestep

    Params => Solver % Values
    IF( ListGetLogical( Params,'Skip Butler Volmer',Found ) ) THEN
      CALL Info('ButlerVolmerUpdate','Skipping update for this solver',Level=8)
      RETURN
    END IF
    
    CALL Info('ButlerVolmerUpdate','Precomputing flux "Jli"',Level=8)

    ! Use native mesh for 'phis' to obtain the fields,
    ! not the 1D mesh for 'cs'
    DO i=1,CurrentModel % NumberOfSolvers
      IF( ListGetString( CurrentModel % Solvers(i) % Values,'Variable',Found ) == 'phis' ) THEN
        Mesh => CurrentModel % Solvers(i) % Mesh
        PSolver => CurrentModel % Solvers(i)
        EXIT
      END IF
    END DO

    ! When visiting this routine for first time, memorize the previes values.
    timestep = GetTimestep()
    IF( timestep /= prevtimestep ) THEN
      IF( ASSOCIATED( JliVar % PrevValues ) ) THEN
        JliVar % PrevValues(:,1) = JliVar % Values 
      END IF
      prevtimestep = timestep
    END IF

    CALL Info('ButlerVolmerUpdate','Solver index related to variable "phis" is: '&
        //I2S(i),Level=10)       
    
    EtaFixed = ListGetLogical( Params,'Fixed Overpotential',Found )
    IF( EtaFixed ) THEN
      CALL Info('ButlerVolmerUpdate','Using previously computed overpotential',Level=8)
    END IF
    
    EtaRelax = ListGetCReal( Params,'Overpotential Relaxation Factor',Found)
    IF(.NOT. Found) EtaRelax = 1.0_dp

    Relax = ListGetCReal( Params,'Butler Volmer Relaxation Factor',Found)
    IF(.NOT. Found) Relax = 1.0_dp
    
    PSolver => Solver
    n = SIZE( JliVar % Values ) 

    DO vari=1,4
      SELECT CASE( vari )
      CASE( 1 )
        str = 'Cs'
      CASE( 2 )
        str = 'Ce'
      CASE( 3 )
        str = 'Phis'
      CASE( 4 )
        str = 'Phie'        
      END SELECT

      CalcSens(vari) = ListGetLogical( Params,'Calculate '//TRIM(str)//' sensitivity',Found )

      ! Check for which variables to compute the sensitivity.
      ! If the sensitivity variable does not exist, create it on-the-fly.
      IF( CalcSens(vari) ) THEN
        SensVar => VariableGet( Mesh % Variables,'dJli d'//TRIM(str))        
        IF(.NOT. ASSOCIATED( SensVar ) ) THEN
          TmpVals => NULL()
          ALLOCATE( TmpVals( n ) )
          TmpVals = 0.0_dp
          CALL VariableAdd( Mesh % Variables, Mesh, PSolver,'dJli d'//TRIM(str),&
              1, TmpVals, JliVar % Perm, Secondary = .TRUE. )        
          SensVar => VariableGet( Mesh % Variables,'dJli d'//TRIM(str))
          CALL Info('ButlerVolmer','Created sensitivity variable: '//TRIM(SensVar % Name))
        END IF
        SensVars(vari) % Variable => SensVar
      END IF
    END DO

    ! We loop over elements since element have the materials.
    ! However, we really have the variables on the nodes.
    ! This is used to tag the nodes that have already been visited. 
    ALLOCATE( NodeDone( Mesh % NumberOfNodes ) )
    NodeDone = .FALSE.

    DO t=1,Mesh % NumberOfBulkElements
      Element => Mesh % Elements(t)

      IF( ANY( CsVar % Perm(Element % NodeIndexes) == 0) ) CYCLE        
      Material => GetMaterial(Element)
      
      n = GetElementNOFNodes(Element)
      DO i=1,n
        j = Element % NodeIndexes(i)          
        IF( NodeDone(j) ) CYCLE           
        NodeDone(j) = .TRUE.

        ! Calculate flux from Butler-Volmer equations
        k = JliVar % Perm(j) 

        Phis = PhisVar % Values( k )
        Phie = PhieVar % Values( PhieVar % Perm(j) )
        Cs = CsVar % Values( k )
        Ce = CeVar % Values( CeVar % Perm(j) )
        Eta = EtaVar % Values( k )
        
        IF( .NOT. ANY( CalcSens ) ) THEN
          Flux = ButlerVolmer( Material, j, Phis, Phie, Cs, Ce, Eta, EtaFixed ) 
          JliVar % Values( k ) = (1-Relax)*JLiVar % Values(k) + Relax * Flux
          IF(.NOT. EtaFixed) EtaVar % Values( k ) = &
              (1-EtaRelax) * EtaVar % Values( k ) + EtaRelax * Eta
          CYCLE
        END IF
        
        ! Calculate the required sensitivities.
        DO vari = 1, 4
          IF( .NOT. CalcSens(vari) ) CYCLE

          IF( vari == 1 ) THEN ! 'Cs'
            Flux = ButlerVolmer( Material, j, Phis, Phie, Cs, Ce, &
                Eta, EtaFixed, djdCs = FluxDer )               
          END IF
          IF( vari == 2 ) THEN ! 'Ce'
            Flux = ButlerVolmer( Material, j, Phis, Phie, Cs, Ce, &
                Eta, EtaFixed, djdCe = FluxDer )               
          END IF
          IF( vari == 3 ) THEN ! 'Phis'
            Flux = ButlerVolmer( Material, j, Phis, Phie, Cs, Ce, &
                Eta, EtaFixed, djdPhis = FluxDer )               
          END IF
          IF( vari == 4 ) THEN ! 'Phie'
            Flux = ButlerVolmer( Material, j, Phis, Phie, Cs, Ce, &
                Eta, EtaFixed, djdPhie = FluxDer )               
          END IF
          SensVars(vari) % Variable % Values(k) = &
              (1-Relax)*SensVars(vari) % Variable % Values(k) + Relax * FluxDer
        END DO

        JliVar % Values( k ) = (1-Relax)*JLiVar % Values(k) + Relax * Flux
        IF(.NOT. EtaFixed) EtaVar % Values( k ) = &
            (1-EtaRelax) * EtaVar % Values( k ) + EtaRelax * Eta

      END DO
    END DO

    CALL ButlerVolmerImbalance() 

    
  CONTAINS

    ! Correct Butler-Volmer potential such that the integral over potential goes
    ! to zero. Assumes that we have "anode" and "cathode" with fluxes of different
    ! signs. 
    !------------------------------------------------------------------------------
    SUBROUTINE ButlerVolmerImbalance() 

      REAL(KIND=dp) :: TotFluxA, TotFluxC, Corr
      LOGICAL :: CorrectFluxes
      
      !n = COUNT( JliVar % Values /= JliVar % Values ) 
      !IF( n > 0 ) THEN
      !  CALL Fatal('ButlerVolmerUpdate',&
      !      'There are NaNs in the flux solution: '//I2S(n))
      !END IF

      ! This is only performed if "Max Output Level" is 7 or more. 
      CALL VariableRange( EtaVar, 7 )

      CorrectFluxes = ListGetLogical( Params,'Correct Butler Volmer Fluxes',Found ) 
      
      IF(.NOT. ( InfoActive(7) .OR. CorrectFluxes ) ) RETURN
      
      TotFluxA = SUM( JliVar % Values * AnodeWeight, AnodeWeight > 0 )
      TotFluxC = SUM( -JliVar % Values * AnodeWeight, AnodeWeight < 0 )

      str = ListGetString(Params,'Equation')
      PRINT *,'Imbalance:'//TRIM(str)//':',&
          (TotFluxA+TotFluxC)/(ABS(TotFLuxA)+ABS(TotFluxC)), TotFluxA, TotFluxC, CalcSens
      
      ! Optionally we may enforce cathode and anode fluxes to be exactly same.
      ! This just fixes disbalance during iteration. It is not a proper fix
      ! but might still be useful. 
      IF( CorrectFluxes ) THEN
        IF( TotFluxA * TotFluxC < -1.0e-8 )  THEN
          Corr = SQRT( -TotFluxC / TotFluxA ) 
        ELSE
          PRINT *,'Cannot correct flux disbalance:',TotFluxA, TotFluxC
          RETURN
        END IF

        PRINT *,'Correcting fluxes with coefficients:',corr
        WHERE( AnodeWeight > 0 ) 
          JliVar % Values = JliVar % Values * Corr
        ELSE WHERE
          JliVar % Values = JliVar % Values / Corr
        END WHERE
      END IF      
      
    END SUBROUTINE ButlerVolmerImbalance
             
    
  END SUBROUTINE ButlerVolmerUpdate



  ! Function to compute cell voltage of given potential solver
  !------------------------------------------------------------------
  FUNCTION CellVoltage( Solver ) RESULT ( VCell ) 
    USE DefUtils
    IMPLICIT NONE
    TYPE(Solver_t) :: Solver
     
    INTEGER, SAVE :: iLeft, iRight
    INTEGER :: jLeft, jRight, bc
    REAL(KIND=dp) :: Vcell, Rf, CurrentDensity, PhiS0, PhiSL
    LOGICAL, SAVE :: Visited = .FALSE.
    TYPE(Mesh_t), POINTER, SAVE :: Mesh
    LOGICAL :: Found, AverageCellVoltage
    
    IF(.NOT. Visited ) THEN
      iLeft = ExtremeLeftNode( MainMesh ) 
      iRight = ExtremeRightNode( MainMesh ) 
      Visited = .TRUE.
    END IF
      
    jLeft = PhisVar % Perm(iLeft)
    jRight = PhisVar % Perm(iRight)

    Rf = ListGetCReal( CurrentModel % Constants,&
        'Current Collector Resistance',UnfoundFatal = .TRUE.)

    ! Here we assume that the current density is constant
    Found = .FALSE.
    DO bc = 1,CurrentModel % NumberOfBCs
      CurrentDensity = ListGetCReal( CurrentModel % BCs(bc) % Values,&
          'Current Density', Found )
      IF( Found ) EXIT
    END DO
    IF(.NOT. Found ) THEN
      CALL Fatal('CellVoltage','Could not find "Current Density" in any BC')
    END IF
    CurrentDensity = ABS( CurrentDensity )

    AverageCellVoltage = ListGetLogical( Solver % Values,'Use Average Cell Voltage',Found ) 
    
    IF( AverageCellVoltage ) THEN
      Phis0 = SUM( PhisVar % Values * AnodeWeight, AnodeWeight > 0 ) / AnodeWeightSum
      PhisL = -SUM( PhisVar % Values * AnodeWeight, AnodeWeight < 0 ) / CathodeWeightSum
    ELSE          
      Phis0 = PhisVar % Values( jLeft )
      PhisL = PhisVar % Values( jRight )
    END IF
      
    ! Cell voltage = PhiS(at x=L) - PhiS(at x=0) - Rf/A*I
    ! eq. 3.25 [1]
    !---------------------------------------------------------
    VCell = PhisL - Phis0 - Rf * CurrentDensity
    
  END FUNCTION CellVoltage
  
END MODULE BatteryModule
