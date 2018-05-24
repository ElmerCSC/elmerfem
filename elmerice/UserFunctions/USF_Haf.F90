!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! USF_Haf.f90
! Calculates height above floatation at the lower boundary based on ice thickness, density, 
! lower surface height and grounded mask (latter is optional).
!
! Call this function from the lower boundary condition (BC) in the sif like this:
!  Haf =  real coord 3
!    Real Procedure "ElmerIceUSF" "Calculate_Haf"
! (where coord 3 is NOT a dummy argument - we actually use this as height of current node 
! relative to sea level!)
!
! Relevant parameters can be given in the lower BC like this:
! Haf GroundedMask = String GroundedMask
! Haf Thickness = String Depth
! Haf rhoRatio = Variable Coordinate 3
!    Real MATC "1.0 * rhoi/rhow"
!
! The second and third are required, the first is optional.  The third, rhoRatio, assumes that 
! rhow and rhoi are defined at the start of the sif, e.g.
! $yearinsec = 365.25*24*60*60
! $rhoi = 910.0/(1.0e6*yearinsec^2)
! $rhow = 1000.0/(1.0e6*yearinsec^2)                                      
!
! Note that if depth is a variable giving the vertical distance of each node from the upper 
! surface then it is equal to the ice thickness when taken at the lower surface.
!
!
! In order to calculate volume above floatation, this user function can be used in conjunction 
! with the SaveScalars Solver like this:
!
! Solver x
!  Equation = SaveScalars
!  Procedure = "SaveData" "SaveScalars"
!  Filename = f.dat
!  Variable 1 = Haf
!  Operator 1 = Int
!  Variable 2 = dummy
!  Operator 2 = Area
!  Exported Variable 1 = Haf
! End
!
! The SaveScalars solver needs, in the lower surface boundary condition:
!  
!  Save Scalars = True
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION Calculate_Haf ( Model, nodenumber, nodeHeight) RESULT(Haf)
  
  USE types
  USE DefUtils
  
  IMPLICIT NONE
  
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: nodeHeight
  INTEGER       :: nodenumber

  REAL(KIND=dp) :: Haf
  
  REAL(KIND=dp) :: rhoRatio, Thickness
  LOGICAL       :: FirstTime=.TRUE., GotIt, UseGroundedMask, UnFoundFatal=.TRUE.
  
  CHARACTER(LEN=MAX_NAME_LEN) :: FunctionName = "USF_Haf", thickName, GroundedMaskName
  TYPE(Variable_t), POINTER   :: GroundedMask, thick
  TYPE(ValueList_t), POINTER  :: BC
  
  SAVE FirstTime, UseGroundedMask, UnFoundFatal, thickName, GroundedMaskName
  
  IF (FirstTime) THEN
     FirstTime = .False.
     
     BC => GetBC(Model % CurrentElement)
     IF (.NOT.ASSOCIATED(BC))THEN
        CALL Fatal(FunctionName, 'No BC Found')
     END IF
     
     GroundedMaskName = GetString( BC, 'Haf GroundedMask', GotIt )
     IF (Gotit) THEN
        UseGroundedMask = .TRUE.
     ELSE
        UseGroundedMask = .FALSE.       
     END IF
     
     thickName = GetString( BC, 'Haf Thickness', GotIt )
     IF (.NOT.Gotit) THEN
        CALL Fatal(FunctionName, 'Cant find >Haf Thickness< in BC')           
     END IF
     
     rhoRatio = GetConstReal( BC, 'Haf rhoRatio', GotIt )
     IF (.NOT.Gotit) THEN
        CALL Fatal(FunctionName, 'Cant find >Haf rhoRatio< in BC')           
     END IF
     
  END IF
  
  GroundedMask => VariableGet(Model % Variables,GroundedMaskName,UnFoundFatal=UnFoundFatal)
  thick => VariableGet(Model % Variables,thickName,UnFoundFatal=UnFoundFatal)
  
  Haf = 0.0_dp
  
  thickness = thick%values(thick%perm( nodenumber ))
  !   nodeHeight = Model%Nodes%z
  
  ! we assume sea level is at zero height.  Hence:
  ! (-nodeHeight/rhoRatio) = 
  ! (vertical distance from lower surface to sea level)*rho_ocean / rho_ice = 
  ! is the thicknes of ice up to floatation. 
  IF (nodeHeight.GT.0.0_dp)  THEN
     Haf = thickness
  ELSE
     Haf = thickness + nodeHeight/rhoRatio
  END IF
  
  ! We assume mask value for ungrounded ice is -1
  IF (UseGroundedMask) THEN
     IF (GroundedMask%values(GroundedMask%perm( nodenumber )) .LT.-0.5 ) THEN 
        Haf = 0.0_dp
     END IF
  END IF
  
END FUNCTION Calculate_Haf
