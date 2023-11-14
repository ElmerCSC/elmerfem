!------------------------------------------------------------------------------
! For copyrights see the BatteryUtils.F90 file in this directory!
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> This routine computes information on the battery status.
!> This is placed in a separate solver such that it can be called only when really
!> needed, e.g. just before saving using the "exec solver = before saving" slot.
!------------------------------------------------------------------------------
SUBROUTINE BatteryPost( Model,Solver,dt,Transient )
  !------------------------------------------------------------------------------
  USE DefUtils
  USE BatteryModule
  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t):: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
  !------------------------------------------------------------------------------
  CHARACTER(*), PARAMETER :: Caller = 'BatteryPost'
  CHARACTER(len = 20) :: GetModel
  TYPE(ValueList_t), POINTER :: Params
  TYPE(ValueList_t), POINTER :: Material
  INTEGER :: t,i,j,k,l,n,dofs,CsNodes,iLeft,iRight,jLeft,jRight,matid,NUnderStoich, &
      SocModel
  INTEGER, POINTER :: CsPerm(:)  
  REAL(KIND=dp) :: Cs_max,pFull,pNill,cs_avg,p1,p2
  LOGICAL :: Anode, Cathode, SurfaceSOC
  REAL(KIND=dp) :: CathodeSOC, AnodeSOC, MeanSOC, FluxCoeff, mincserr, maxcserr, &
      CellV, Vlimit, AnodeQ, CathodeQ, ElectrolyteQ, AnodeU, CathodeU, &
      ElectrolyteU, EndArea, pActive, StoichLimit, DischTime, CurrApp, &
      Capacity, SOCAlt, CellCapacity, InitSOC, PrevSOC
  LOGICAL :: Found, Visited = .FALSE.
  TYPE(Variable_t), POINTER :: AveVarPtr 


  SAVE Visited, iLeft, iRight, jLeft, jRight, PrevSOC, InitSOC

  !------------------------------------------------------------------------------
  !    Get variables needed for solution
  !------------------------------------------------------------------------------

  CALL Info(Caller,'-------------------------------------------------------')
  CALL Info(Caller,'Postprocessing information on battery status')
  CALL Info(Caller,'-------------------------------------------------------')  

  Params => GetSolverParams()

  EndArea = ListGetCReal( Model % Constants,'Electrode Plate Area',Found )
  IF(.NOT. Found ) EndArea = 1.0_dp

  IF(.NOT. Visited) THEN  
    iLeft = ExtremeLeftNode( SubMesh ) 
    iRight = ExtremeRightNode( SubMesh )     
    jLeft = Cs1DVar % Perm(iLeft)
    jRight = Cs1DVar % Perm(iRight)
    InitSOC =  ListGetCReal( Params, 'Initial SOC', Found)
    IF ( .NOT. Found ) InitSOC = 1
    PrevSOC = InitSOC
    Visited = .TRUE.
  END IF

  CsPerm => CsVar % Perm  
  CsNodes = SIZE( CsVar % Values ) 

  CsFullVar => VariableGet( MainMesh % Variables,'Cs profile')
  dofs = CsFullVar % dofs

  ! Integrate over the surface flux
  IF( ASSOCIATED( JliIntegVar ) ) THEN
    IF( UseTimeAveFlux ) THEN
      JliIntegVar % Values = JliIntegVar % Values + &
          0.5_dp * ( JliVar % Values + Jli0 ) * dt
    ELSE
      JliIntegVar % Values = JliIntegVar % Values + JliVar % Values * dt
    END IF
  END IF

  ! Compute the average of solid phase potential
  ! It is not enough to take a number average since we are integrating over
  ! a sphere. Instead use averaging over the whole volume integral.
  !----------------------------------------------------------------------  
  DO j=1, MainMesh % NumberOfNodes
    k = CsPerm(j)
    IF( k == 0 ) CYCLE

    ! Copy the 1D concentration values related to a node
    Cs1DVar % Values(1:dofs) = CsFullVar % Values(dofs*(k-1)+1:dofs*k)    
    CsVar % Values(k) = Cs1DVar % Values(jRight)

    ! Calculate average of 1D mesh using the weight vector precomputed by 1D solver.
    ! Note: volume of unit mesh is 1.0/3.0
    cs_avg = 3.0_dp * SUM( MassMult1D * Cs1DVar % Values )

    CsAveVar % Values(k) = cs_avg

    ! The difference in the extreme outer and inner locations of the 1d sphere    
    IF( ASSOCIATED( CsDiffVar ) ) THEN
      p1 = Cs1dVar % Values(jright)
      p2 = Cs1dVar % Values(jleft)
      CsDiffVar % Values(CsDiffVar % Perm(j)) = p2 - p1
    END IF

  END DO

  SurfaceSOC = ListGetLogical( Params,&
      'SOC on Surface',Found ) 
  IF( SurfaceSOC ) THEN
    CALL Info(Caller,'Using surface concentration to evaluate SOC',Level=7)
    AveVarPtr => CsVar
  ELSE
    AveVarPtr => CsAveVar
  END IF

  ! Checks which SOC model user have chosen
  GetModel = GetString( Solver % Values, 'Soc Model',Found)
  IF( .NOT. Found ) GetModel = 'default'

  SELECT CASE( GetModel )
  CASE('default')
    SocModel = 1
  CASE('simplified')
    SocModel = 2
  CASE('coulomb')
    SocModel = 3
  END SELECT

  DO matid = 1, CurrentModel % NumberOfMaterials
    Material => CurrentModel % Materials(matid) % Values
    Anode = ListGetLogical( Material,'Anode',Found )
    Cathode = ListGetLogical( Material,'Cathode',Found ) 
    IF(.NOT. ( Anode .OR. Cathode ) ) CYCLE

    pFull = ListGetCReal( Material,'Stoichiometry at Full Charge' ) 
    pNill = ListGetCReal( Material,'Stoichiometry at Nill Charge' )
    Cs_Max = ListGetCReal( Material,'Maximum solid phase concentration')
    Anode = ListGetLogical( Material,'Anode',Found )
    Cathode = ListGetLogical( Material,'Cathode',Found ) 
    pActive = ListGetCReal( Material,'Active Particle Volume Fraction')

    ! Computing State of Charge for anode and cathode
    ! SoC = ((Cs_avg/Cs_max)-p0%)/(p100%-p0%)
    IF( SocModel /= 3 ) THEN
      IF( Anode ) THEN
        IF ( SocModel == 2 ) THEN
          WHERE  ( AnodeWeight > 0 )
            SOCVar % Values = (AveVarPtr % Values/Cs_max)
          END WHERE
        ELSE
          WHERE  ( AnodeWeight > 0 )
            SOCVar % Values = ((AveVarPtr % Values/Cs_max)-pNill)/(pFull-pNill)
          END WHERE
        END IF
      ELSE
        IF ( SocModel == 2 ) THEN
          WHERE  ( AnodeWeight < 0 )
            SOCVar % Values = (AveVarPtr % Values/Cs_max)
          END WHERE
        ELSE
          WHERE  ( AnodeWeight < 0 )
            SOCVar % Values = ((AveVarPtr % Values/Cs_max)-pNill)/(pFull-pNill)
          END WHERE
        END IF
      END IF
    END IF

    IF( Anode ) THEN
      AnodeQ = FaradayConstant * EndArea * pActive * &
          SUM( AveVarPtr % Values * AnodeWeight, AnodeWeight > 0 )
    ELSE
      CathodeQ = FaradayConstant * EndArea * pActive * &
          SUM( -AveVarPtr % Values * AnodeWeight, AnodeWeight < 0 )
    END IF

    IF( ASSOCIATED( CsErrVar ) ) THEN
      FluxCoeff = SolidFluxScaling( Material )        
      IF( Anode ) THEN
        WHERE ( AnodeWeight > 0 )                    
          CsErrVar % Values = (CsInitVar % Values - CsAveVar % Values ) / &
              ( 3 * FLuxCoeff * JliIntegVar % Values ) - 1
        END WHERE
      ELSE
        WHERE ( AnodeWeight < 0 )          
          CsErrVar % Values = (CsInitVar % Values - CsAveVar % Values ) / &
              ( 3 * FLuxCoeff * JliIntegVar % Values ) - 1
        END WHERE
      END IF
      mincserr = MINVAL( CsErrVar % Values )
      maxcserr = MAXVAL( CsErrVar % Values )
      CALL ListAddConstReal( Model % Simulation,'res: min cs error',mincserr)
      CALL ListAddConstReal( Model % Simulation,'res: max cs error',maxcserr)      
    END IF

  END DO

  IF( SocModel /= 3 ) THEN
    AnodeSOC = SUM( SOCVar % Values * AnodeWeight, AnodeWeight > 0 ) 
    CathodeSOC = -SUM( SOCVar % Values * AnodeWeight, AnodeWeight < 0 ) 
    MeanSOC = ( AnodeSOC + CathodeSOC ) / ( AnodeWeightSum + CathodeWeightSum )

    AnodeSOC = AnodeSOC / AnodeWeightSum
    CathodeSOC = CathodeSOC / CathodeWeightSum

    CALL VariableRange( SOCVar, 6, AnodeWeight )

    CALL ListAddConstReal( Model % Simulation,'res: Anode SoC',AnodeSOC )
    CALL ListAddConstReal( Model % Simulation,'res: Cathode SoC',CathodeSOC )
    CALL ListAddConstReal( Model % Simulation,'res: Mean SoC',MeanSOC )
  END IF

  IF( ListGetLogical( Params,'Calculate Charges',Found ) ) THEN
    CALL IntegrateElectrolyte()  

    CALL ListAddConstReal( Model % Simulation,'res: Anode charge',AnodeQ )
    CALL ListAddConstReal( Model % Simulation,'res: Cathode charge',CathodeQ )
    CALL ListAddConstReal( Model % Simulation,'res: Electrolyte charge',ElectrolyteQ )

    AnodeU = SUM( AnodeWeight * PhiSVar % Values, AnodeWeight > 0 ) / AnodeWeightSum
    CathodeU = SUM( -AnodeWeight * PhiSVar % Values, AnodeWeight < 0 ) / CathodeWeightSum

    CALL ListAddConstReal( Model % Simulation,'res: Anode Mean Potential',AnodeU )
    CALL ListAddConstReal( Model % Simulation,'res: Cathode Mean Potential',CathodeU )   
    CALL ListAddConstReal( Model % Simulation,'res: Electrolyte Mean Potential',ElectrolyteU )
  END IF

  ! Checks that stoichiometric limit is not passed
  !----------------------------------------------------
  StoichLimit = ListGetCReal( Params, 'Stoichiometric limit', Found)
  NUnderStoich = COUNT( CsVar % Values <= StoichLimit )
  IF ( Found .AND. NUnderStoich > 0 ) THEN
    CALL Warn(Caller,'Discharge in '//I2S(NUnderStoich)//&
        ' nodes has reached the stoichiometric limit!')  
  END IF

  ! Check that there are no negative concentrations.
  ! If there are then we need to do something about it!
  !---------------------------------------------------
  i = COUNT( CeVar % Values < 0.0_dp )
  IF( i > 0 ) THEN
    n = SIZE( CeVar % Values )
    CALL Warn(Caller,'Number of negative concentrations (out of '&
        //I2S(n)//') for Ce: '//I2S(i))
  END IF

  i = COUNT( CsVar % Values < 0.0_dp )
  IF( i > 0 ) THEN
    n = SIZE( CsVar % Values )
    CALL Warn(Caller,'Number of negative concentrations (out of '&
        //I2S(n)//') for Cs: '//I2S(i))
  END IF

  ! Solves the Cell Voltage
  !---------------------------------------------------
  CellV = CellVoltage( Solver )

  WRITE(Message,'(A,ES12.5)') 'Cell voltage is: ', CellV
  CALL Info(Caller,Message,Level=7)
  CALL ListAddConstReal( Model % Simulation,'res: cell voltage',CellV)

  VLimit = ListGetCReal( Params, 'Cutoff voltage', Found)
  IF (Found .AND. CellV <= Vlimit) THEN
    CALL Warn(Caller,'Cell Voltage has reached the cut-off Voltage!')
    CALL Info(Caller,'Setting exit condition for early termination!')
    CALL ListAddConstReal( Model % Simulation,'Exit Condition',1.0_dp)
  END IF

  ! Calculate cell capacity in Ah/m^2
  !----------------------------------------------------------------------
  CurrApp = ListGetCReal( Params, 'Applied Current', Found)
  IF ( found ) THEN
    DischTime = GetTime()
    Capacity = ABS( CurrApp ) * (DischTime/3600)
    CALL ListAddConstReal( Model % Simulation,'res: capacity',Capacity)
    CALL ListAddConstReal( Model % Simulation,'res: applied current',CurrApp)

    ! Alternative way to calculate SoC with Coulomb counting
    ! Faradaic efficiency is assumed as 1
    !----------------------------------------------------------------------
    IF( SocModel == 3 ) THEN
      CellCapacity = ListGetCReal( Params, 'Cell Capacity', Found)
      IF (found) THEN
        SOCAlt = InitSOC - ( (1/( CellCapacity * 3600 ) ) * CurrApp * dt + (1 - PrevSOC) )
        PrevSOC = SOCAlt
        CALL ListAddConstReal( Model % Simulation,'res: Coulomb SoC',SOCAlt)
      END IF
    END IF
  END IF


CONTAINS

  ! Study electrolyte info also.
  ! There are conflicting material parameters at the interfaces so we
  ! cannot precompute the weights as in the case of solid phase.
  !---------------------------------------------------------------------
  SUBROUTINE IntegrateElectrolyte()

    REAL(KIND=dp) :: Basis(8), ElemPot(8), DetJ, Weight, SumWeight
    TYPE(Nodes_t), SAVE :: Nodes      
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: Material, PrevMaterial
    TYPE(Element_t), POINTER :: Element
    INTEGER :: t, n, elem
    LOGICAL :: Anode, Stat, AnodeOrCathode           
    TYPE(Mesh_t), POINTER :: Mesh

    Mesh => GetMesh()

    PrevMaterial => NULL()
    ElectrolyteQ = 0.0_dp
    ElectrolyteU = 0.0_dp
    SumWeight = 0.0_dp

    DO elem=1,Mesh % NumberOfBulkElements
      Element => Mesh % Elements(elem)

      Material => GetMaterial(Element)
      IF( .NOT. ASSOCIATED( Material, PrevMaterial ) ) THEN
        PrevMaterial => Material
        pActive = ListGetCReal( Material,'Electrolyte Volume Fraction')
      END IF

      n = GetElementNOFNodes( Element )        
      CALL GetElementNodes( Nodes, UElement=Element )

      IP = GaussPoints( Element )

      ElemPot(1:n) = PhiEVar % Values( PhiEVar % Perm( Element % NodeIndexes ) ) 

      DO t=1,IP % n
        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t), detJ, Basis )
        Weight = IP % s(t) * DetJ

        SumWeight = SumWeight + Weight
        ElectrolyteQ = ElectrolyteQ + Weight * pActive
        ElectrolyteU = ElectrolyteU + Weight * SUM( Basis(1:n) * ElemPot(1:n) )
      END DO
    END DO

    ElectrolyteQ = FaradayConstant * EndArea * ElectrolyteQ
    ElectrolyteU = ElectrolyteU / SumWeight

  END SUBROUTINE IntegrateElectrolyte

END SUBROUTINE BatteryPost






