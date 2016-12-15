!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! USF_Damage.f90
! Last modif : 12 December 2014
! Computes the evolution of damage and computes the Enhancement factor of the Glen's law
! as a function of damage evolution 
!
!
! (1) Enhancement Factor 
! Need some inputs in the sif file.
! Parameters: 
! Glen Exponent
!             
!
! (2) SourceDamage
! Need some inputs in the sif file.
! Parameters: 
! Damage Enhancement Factor
! Damage Parameter sigmath
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! EnhancementFactor
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION EnhancementFactor ( Model, nodenumber, D) RESULT(E)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   TYPE(ValueList_t), POINTER :: Material
   TYPE(Solver_t), TARGET :: Solver
   REAL(KIND=dp) :: D, E, n  
   INTEGER :: nodenumber
   LOGICAL :: FirstTime=.TRUE., GotIt

   SAVE FirstTime, n 

   IF (FirstTime) THEN
   FirstTime = .False.
    
      Material => GetMaterial()
      n = GetConstReal( Material, 'Glen Exponent', GotIt )
      IF (.NOT.GotIt) THEN
         WRITE(Message,'(A)') 'Variable Glen Exponent not found. &
              &Setting to 3.0'
         CALL INFO('Damage EnhancementFactor', Message, level=2)
         n = 3.0_dp
      ELSE
         WRITE(Message,'(A,F10.4)') 'n = ', n 
         CALL INFO('Damage EnhancementFactor', Message, level=2)
      END IF
   END IF

       E = (1.0 - D)**(-n) 
    
  ! write(*,*) D
  ! write(*,*)'E', E
END FUNCTION EnhancementFactor


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SourceDamage 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION SourceDamage (Model, nodenumber, D) RESULT(Source)


   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   USE GeneralUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   REAL (KIND=dp) :: D, Source          
   INTEGER :: nodenumber
 
   TYPE(Solver_t):: Solver 
   TYPE(ValueList_t), POINTER :: Material, Constants
   TYPE(Variable_t), POINTER :: StressVariable, FlowVariable, ChiVariable, PSeaDVariable
   REAL(KIND=dp), POINTER :: StressValues(:), FlowValues(:), ChiValues(:), PSeaDValues(:)
   INTEGER, POINTER :: StressPerm(:), FlowPerm(:), ChiPerm(:), PSeaDPerm(:)

   INTEGER :: Ind(3,3), DIM, i, j, indice(3), infor
   REAL (KIND=dp) :: Sig(3,3), SigDev(3,3), EigVect(3,3), EigValues(3),tmp 
   REAL (KIND=dp) :: SigmaI, SigmaII, Chi, B, sigmath, pwater
   REAL (KIND=DP) :: EI(3),Dumy(1),Work(24), sigmath_var, stress_threshold
   REAL (KIND=DP) :: u, v, nbrPi, s
   LOGICAL :: GotIt, FirstTime = .TRUE., Cauchy, UnFoundFatal=.TRUE.
   CHARACTER*20 :: USF_Name='SourceDamage'


   SAVE :: Ind, DIM, sigmath, B
   SAVE :: FirstTime, Cauchy
   

   IF (FirstTime) THEN
      FirstTime = .FALSE.  
      DIM = CoordinateSystemDimension()

      DO i=1, 3
         Ind(i,i) = i
      END DO
      Ind(1,2) = 4
      Ind(2,1) = 4
      Ind(2,3) = 5
      Ind(3,2) = 5
      Ind(3,1) = 6
      Ind(1,3) = 6

    Material => GetMaterial()
    !Read the coefficients B and sigmath  

    B = GetConstReal( Material, 'Damage Enhancement Factor', GotIt )
      IF (.NOT.GotIt) THEN
      CALL FATAL('Damage Source', 'Damage Enhancement Factor B not Found')
      ELSE
         WRITE(Message,'(A,F10.4)') 'Damage Enhancement Factor = ', B
         CALL INFO('Damage Source', Message, level=2)
      ENDIF

      sigmath = GetConstReal( Material, 'Damage Parameter sigmath', GotIt )
      IF (.NOT.GotIt) THEN
      CALL FATAL('Damage Source', 'Damage Parameter Sigmath not Found')
      ELSE
         WRITE(Message,'(A,F10.4)') 'Damage Parameter sigmath = ', sigmath 
         CALL INFO('Damage Source', Message, level=2)
      ENDIF

   ! Cauchy or deviatoric stresses ?
      Cauchy = ListGetLogical( Material , 'Cauchy', Gotit )
      WRITE(Message,'(A,L1)') 'Cauchy stress tensor computed ? ', Cauchy 
         CALL INFO('Damage Source', Message, level=2)
   END IF ! FirstTime


   ! Get the Stress                     
   StressVariable => VariableGet( Model % Variables, 'Stress',UnFoundFatal=UnFoundFatal)
   StressPerm    => StressVariable % Perm
   StressValues  => StressVariable % Values

   ! Get Chi variable (positive where damage increases)
   ChiVariable => VariableGet( Model % Variables, 'Chi',UnFoundFatal=UnFoundFatal)
   ChiPerm    => ChiVariable % Perm
   ChiValues  => ChiVariable % Values
   
   ! Get the Sea Damage Pressure
   PSeaDVariable => VariableGet( Model % Variables, 'PSeaD' )
   IF ( ASSOCIATED( PSeaDVariable ) ) THEN
      PSeaDPerm   => PSeaDVariable % Perm
      PSeaDValues => PSeaDVariable % Values
   ELSE
      CALL WARN('Damage Source', 'PSeaD not associated, basal pressure not taken into account in damage formation' )
      CALL WARN('Damage Source', 'Taking default value PSeaD=0.0')
   END IF

   ! Get the variables to compute the hydrostatic pressure  
   FlowVariable => VariableGet( Model % Variables, 'Flow Solution',UnFoundFatal=UnFoundFatal)
   FlowPerm    => FlowVariable % Perm
   FlowValues  => FlowVariable % Values

   Sig = 0.0
   DO i=1, DIM
      DO j= 1, DIM
         Sig(i,j) =  &
              StressValues( 2*DIM *(StressPerm(Nodenumber)-1) + Ind(i,j) )
      END DO
   END DO
   IF (DIM==2) Sig(3,3) = StressValues( 2*DIM *(StressPerm(Nodenumber)-1) + Ind(3,3))

  ! S = sigma + p
  ! Need Cauchy Stress and Deviatoric Stress 
   IF (.NOT.Cauchy) THEN ! If Deviatoric Stress is computed, then, get the
                         ! Cauchy Stress
       SigDev = Sig
       DO i=1,3  
           Sig(i,i) = SigDev(i,i) - FlowValues((DIM+1)*FlowPerm(Nodenumber))
       END DO
   ELSE ! If the Cauchy Stress is computed, then get the Deviatoric Stress 
       DO i=1,3  
           SigDev(i,i) = Sig(i,i) + FlowValues((DIM+1)*FlowPerm(Nodenumber))
            !write(*,*)Sig(i,1), Sig(i,2), Sig(i,3)
       END DO
   END IF


   ! Compute the principal stresses:

   ! Get the principal stresses
   CALL DGEEV('N','N',3,Sig,3,EigValues,EI,Dumy,1,Dumy,1,Work,24,infor )
   IF (infor.ne.0) &
   CALL FATAL('Compute EigenValues', 'Failed to compute EigenValues') 

   ! Get the eigenvectors (if necessary)
   CALL DGEEV('N','V',3,Sig,3,EigValues,EI,Dumy,1,EigVect,3,Work,24,infor ) 
   IF (infor.ne.0) &
   CALL FATAL('Compute EigenVectors', 'Failed to compute EigenVectors') 
   
   indice = (/(i,i=1,3)/)
   CALL sortd(3,EigValues,indice)
 
   SigmaI = EigValues(3)
   SigmaII = EigValues(2) 

   !SigmaI = MAXVAL(EigValues)
   !write(*,*)'SigI',SigmaI,EigValues
   nbrPi = 3.141592_dp

   u = EvenRandom()
   v = EvenRandom()

   ! Determination of the stress threshold   
      Constants => GetConstants()
   s = GetConstReal( Constants, 'Dev Tensile Strength Modifier', GotIt) 
    IF (.NOT.GotIt) THEN
      CALL FATAL('USF_Damage','No "Dev tensile strength modifier" given, set &
      to 0.05')
    END IF
  
   ! Get a normal distribution of mean 0 and std dev "s" using two random numbers
   ! "u" and "v" 
   sigmath_var = ABS(0+s*SQRT(-2.0_dp*LOG((1.0_dp-u)))*COS(2.0_dp*nbrPi*v))

   stress_threshold = sigmath*(1+sigmath_var)
  
   ! Get the Sea pressure at the node
   
   IF ( ASSOCIATED( PSeaDVariable ) ) THEN
   pwater = PseaDValues ( PSeaDPerm (nodenumber) )
   ELSE
   pwater = 0.0_dp
   END IF

   ! Damage Criterion
   Chi = 1.0_dp / (1.0_dp - D) * (SigmaI + pwater) - stress_threshold

   ! Save Chi in ChiVariable
   ChiValues(ChiPerm(nodenumber)) = Chi
 
   ! Advection-Reaction Source term :
   Source = B * MAX(Chi,0.0_dp)

END FUNCTION SourceDamage   
