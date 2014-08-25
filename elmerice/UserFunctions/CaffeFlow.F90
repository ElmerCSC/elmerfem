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
! *  Authors:  Juha Ruokolainen, Hakime Seddik
! *  Email:   hakime@pop.lowtem.hokudai.ac.jp
! *  Address: Institute of Low Temperature Science  
! *     Hokkaido University                   
! *     Sapporo-shi Kita-ku, Kita 19, nishi 6 
! *     Hokkaido ; Japan                      
! *     EMail: hakime@pop.lowtem.hokudai.ac.jp
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 9 July 1997
! *  Modification Date: 
! * 
! *****************************************************************************
!> CaffeFlow.f90  anisotropic model for ice flow
!> Caffe model: viscosity factor as a function of ice anisotropy and temperature
FUNCTION caffeGetViscosity ( Model, nodenumber, temperature ) RESULT(anisoVisFact)
   USE types
   USE DefUtils
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
!-----------------------------------------------------------
   IMPLICIT NONE
!------------ external variables ---------------------------
   TYPE(Model_t) :: Model
   INTEGER :: nodenumber
    REAL(KIND=dp) :: temperature, anisoVisFact
!------------ internal variables----------------------------
   TYPE(Element_t), POINTER :: CurrentElement
   TYPE(Nodes_t) :: ElementNodes
   TYPE(Variable_t), POINTER :: FlowVariable, FabricVariable, TempVariable
   TYPE(ValueList_t), POINTER :: Material
   REAL(KIND=dp), POINTER :: Hwrk(:,:,:)
   REAL(KIND=dp) :: LocalVelo(3, Model % MaxElementNodes)
   REAL(KIND=dp), DIMENSION(Model % MaxElementNodes) :: Expo
   REAL(KIND=dp) ::&
       rateFactor, aToMinusOneThird, gasconst, temphom
   REAL (KIND=dp), ALLOCATABLE :: activationEnergy(:,:), arrheniusFactor(:,:),&
       aniEnhancementFact(:), viscosityExponent(:), PressureMeltingPoint(:),&
       LimitTemp(:)
   REAL(KIND=dp), POINTER ::  FlowValues(:), FabricValues(:), TempValues(:)           
   REAL(KIND=dp) :: Basis( Model % MaxElementNodes ), ddBasisddx(1,1,1), dBasisdx( Model % MaxElementNodes , 3 )
   REAL(KIND=dp) :: u, v, w, detJ, LGrad(3,3), Dij(3,3)
   REAL (KIND=dp) :: a1133, a2233, a3333, a3331, a3332, a3312
   REAL (KIND=dp) :: AfirstArg,AsecondArg, Aconst, Eps
   REAL (KIND=dp) :: Emin, MinGamma
   REAL(KIND=dp) :: a2(6), a4(9), A, E, m, Temp
   INTEGER, POINTER :: FlowPerm(:), FabricPerm(:), TempPerm(:), NodeIndexes(:)
   INTEGER :: DIM, nMax, t, i, j, k, STDOFs, n
   INTEGER :: material_id, body_id,  elementNbNodes, nodeInElement, istat
   LOGICAL :: stat, FirstTime=.TRUE., GotIt
   CHARACTER(LEN=MAX_NAME_LEN) :: TempName
  
!------------ remember this -------------------------------
   SAVE ElementNodes, DIM, FirstTime, gasconst, activationEnergy, arrheniusFactor,&
         aniEnhancementFact, viscosityExponent, Hwrk, PressureMeltingPoint, &
        LimitTemp

!-----------------------------------------------------------
! Read in constants from SIF file and do some allocations
!-----------------------------------------------------------
IF (FirstTime) THEN   
 ! inquire coordinate system dimensions  and degrees of freedom from NS-Solver
 ! ---------------------------------------------------------------------------

        ALLOCATE( ElementNodes % x ( Model % MaxElementNodes ), &
                  ElementNodes % y ( Model % MaxElementNodes ), &
                  ElementNodes % z ( Model % MaxElementNodes ) ) 

        DIM = CoordinateSystemDimension()

        gasconst = ListGetConstReal( Model % Constants,'Gas Constant',GotIt)
        IF (.NOT. GotIt) THEN
        gasconst = 8.314D00 ! m-k-s
        WRITE(Message,'(a,e10.4,a)') 'No entry for Gas Constant (Constants) in input file found. Setting to ',&
             gasconst,' (J/mol)'
        CALL INFO('CAFFE (CAFFEViscosity)', Message, level=4)
        END IF
        nMax = Model % MaxElementNodes
        ALLOCATE(activationEnergy(2,nMax),&
        arrheniusFactor(2,nMax),&
        LimitTemp( nMax),&
        aniEnhancementFact(nMax),&
        PressureMeltingPoint( nMax ),&
        viscosityExponent(nMax),&
        STAT=istat)
        IF ( istat /= 0 ) THEN
            CALL Fatal('CAFFE (CAFFEViscosity)','Memory allocation error, Aborting.')
        END IF
        NULLIFY( Hwrk )
        FirstTime = .FALSE.
     CALL Info('CAFFE (CAFFEViscosity)','Memory allocations done', Level=3)           
END IF

!---------------------------------------------
! get element properties and solver variables
!--------------------------------------------- 
       FlowVariable => VariableGet( Model % Solver % Mesh % Variables, "flow solution" )
       IF ( ASSOCIATED( FlowVariable ) ) THEN
         FlowPerm    => FlowVariable % Perm
         FlowValues  => FlowVariable % Values
       ELSE
         CALL Info('CAFFE', &
                      & 'No variable for velocity associated.', Level=4)
       END IF

       STDOFs =  FlowVariable % DOFs
      
       FabricVariable => VariableGet( Model % Solver % Mesh % Variables, "fabric" )
       IF ( ASSOCIATED( FabricVariable ) ) THEN
         FabricPerm    => FabricVariable % Perm
         FabricValues  => FabricVariable % Values
       ELSE
         CALL Info('CAFFE', &
                      & 'No variable for fabric associated.', Level=4)
       END IF

        
!      TempVariable => VariableGet( Model % Solver % Mesh % Variables, 'Temperature' )
!      IF ( ASSOCIATED( TempVariable) ) THEN
!        TempPerm    => TempVariable % Perm
!        TempValues => TempVariable % Values
!      END IF

 body_id = Model % CurrentElement % BodyId
           
 material_id = ListGetInteger(Model % Bodies(body_id) % Values, 'Material', GotIt)
 IF (.NOT.GotIt) CALL FATAL('CAFFE (CAFFEViscosity)','No Material ID found')

 CurrentElement => Model % CurrentElement 
 elementNbNodes = GetElementNOFNodes()
 IF (.NOT. GotIt) THEN
    WRITE(Message,'(a,I2,a,I2,a)') 'No material id for current element of node ',nodenumber,', body ',body_id,' found'
     CALL FATAL('CAFFE (CAFFEViscosity)', Message)
 END IF
 NodeIndexes => CurrentElement % NodeIndexes


Material => Model % Materials(material_id) % Values
IF (.NOT.ASSOCIATED(Material)) THEN  
   WRITE(Message,'(a,I2,a,I2,a)') 'No Material for current element of node ',nodenumber,', body ',body_id,' found'
   CALL FATAL('CAFFE (CAFFEViscosity)',Message)
END IF
       
  ElementNodes % x(1:elementNbNodes) = Model % Nodes % x(NodeIndexes(1:elementNbNodes))
  ElementNodes % y(1:elementNbNodes) = Model % Nodes % y(NodeIndexes(1:elementNbNodes))
  ElementNodes % z(1:elementNbNodes) = Model % Nodes % z(NodeIndexes(1:elementNbNodes))
 
! 2D U,V,p    STDOFs=3
! 3D U,V,W,p  STDOFs=4
  LocalVelo = 0.0_dp
  DO i = 1, STDOFs - 1
     LocalVelo(i,1:elementNbNodes) = FlowValues( STDOFs*(FlowPerm(NodeIndexes(1:elementNbNodes))-1) + i)
  END DO

!          Locala11(1:n) = FabricValues( 5 * (FabricPerm(NodeIndexes(1:n))-1) + 1 )
!          Locala22(1:n) = FabricValues( 5 * (FabricPerm(NodeIndexes(1:n))-1) + 2 )
!          Locala12(1:n) = FabricValues( 5 * (FabricPerm(NodeIndexes(1:n))-1) + 3 )
!          Locala23(1:n) = FabricValues( 5 * (FabricPerm(NodeIndexes(1:n))-1) + 4 )
!          Locala31(1:n) = FabricValues( 5 * (FabricPerm(NodeIndexes(1:n))-1) + 5 )
           

DO  i=1,elementNbNodes
     IF (nodenumber == NodeIndexes( i )) THEN

!-------------------------------------------
! get material properties
!-------------------------------------------
  ! activation energies
  !--------------------
  CALL ListGetRealArray( Material,'Activation Energies',Hwrk, elementNbNodes, &
  Model % CurrentElement % NodeIndexes, GotIt )
  IF (.NOT. GotIt) THEN
  WRITE(Message,'(a,I2,a,I2)') 'No Value for Activation Energy  found in Material ', material_id,' for node ', nodenumber
  CALL FATAL(' CAFFE (CAFFEViscosity)',Message)
  END IF
  IF ( SIZE(Hwrk,2) == 1 ) THEN
     DO j=1,MIN(3,SIZE(Hwrk,1))
        activationEnergy(j,1:elementNbNodes) = Hwrk(j,1,1:elementNbNodes)
     END DO
  ELSE
     WRITE(Message,'(a,I2,a,I2)') 'Incorrect array size for Activation Energy in Material ', material_id,' for node ', nodenumber
     CALL FATAL('CAFFE (CAFFEViscosity)',Message)
  END IF

! Arrhenius Factors
! ------------------
  CALL ListGetRealArray( Material,'Arrhenius Factors',Hwrk,elementNbNodes, &
       Model % CurrentElement % NodeIndexes, GotIt )
  IF (.NOT. GotIt) THEN
     WRITE(Message,'(a,I2,a,I2)') 'No Value for Arrhenius Factors  found in Material ', material_id,' for node ', nodenumber
     CALL FATAL('CAFFE (CAFFEViscosity)',Message)
  END IF
  IF ( SIZE(Hwrk,2) == 1 ) THEN
     DO j=1,MIN(3,SIZE(Hwrk,1))
        arrheniusFactor(j,1:elementNbNodes) = Hwrk(j,1,1:elementNbNodes)
     END DO
  ELSE
     WRITE(Message,'(a,I2,a,I2)') 'Incorrect array size for Arrhenius Factors in Material ', material_id,' for node ', nodenumber
     CALL FATAL('CAFFE (CAFFEViscosity)',Message)
  END IF

! Threshold temperature for switching activation energies and Arrhenius factors
!------------------------------------------------------------------------------
  LimitTemp(1:elementNbNodes) = ListGetReal( Material,'Limit Temperature', elementNbNodes,&
  Model % CurrentElement % NodeIndexes, GotIt )
  IF (.NOT. GotIt) THEN
     LimitTemp(1:elementNbNodes) = -1.0D01
     WRITE(Message,'(a,I2,a,I2,a)') 'No keyword >Limit Temperature< found in Material ',&
          material_id,' for node ', nodenumber, '.setting to -10'
     CALL INFO('CAFFE (CAFFEViscosity)', Message, level=4)
  END IF

! Viscosity Exponent
!-------------------
  viscosityExponent(1:elementNbNodes) = ListGetReal( Material,'Viscosity Exponent', elementNbNodes, &
Model % CurrentElement % NodeIndexes, GotIt )
  IF (.NOT. GotIt) THEN
     viscosityExponent(1:elementNbNodes) = 1.0D00/3.0D00
     WRITE(Message,'(a,I2,a,I2,a)') 'No Viscosity Exponent found in Material ', material_id,&
' for node ', nodenumber, '.setting k=1/3'
     CALL INFO('CAFFE (CAFFEViscosity)', Message, level=4)
  END IF
           
! Enhancement factor for anisotropic flow law 
! -------------------------------------------
  aniEnhancementFact(1:elementNbNodes) = ListGetReal( Material,'Anisotropic Enhancement Factor', &
elementNbNodes, Model % CurrentElement % NodeIndexes, GotIt )
  IF (.NOT. GotIt) THEN
     WRITE(Message,'(a,I2,a,I2,a)') 'Value for Maximum Enhancement factor not found in Matrial', &
material_id,' for node ', nodenumber
     CALL Fatal('CAFFE (CAFFEViscosity)', Message)
  END IF

! Critical Enhancement factor to avoid singularities in anisotropic flow law
! --------------------------------------------------------------------------
  Emin = GetConstReal(Material, 'Critical Enhancement factor', GotIt)
      IF(.NOT. GotIt) THEN 
          Emin=0.0001
          WRITE(Message,'(a,I2,a,I2,a)') 'Value for Critical Enhancement factor not found in Matrial', &
material_id,' for node ', nodenumber, '.setting Emin=1.Oe-4'
          CALL INFO('CAFFE (CAFFEViscosity)', Message, level=3)
      END IF

! Critical parameter to avoid singulairties in computing the polycrystall deformability
! -------------------------------------------------------------------------------------
  MinGamma = GetConstReal(Material, 'Critical Strain Rate', GotIt)
      IF(.NOT. GotIt) THEN 
         MinGamma=1.0e-15
         WRITE(Message,'(a,I2,a,I2,a)') 'Critical Strain rate not found in Matrial', material_id, &
' for node ', nodenumber, '.setting MinGamma=1.0e-15'
         Call INFO('CAFFE (CAFFEViscosity)', Message, level=3)
      END IF

! Pressure Melting Point and homologous temperature
! --------------------------------------------------
  TempName =  GetString(Material  ,'Temperature Name', GotIt)
  IF (.NOT.GotIt) CALL FATAL('CAFFE (CAFFEViscosity)','No Temperature Name found')
  PressureMeltingPoint(1:elementNbNodes) =&
       ListGetReal( Material, TRIM(TempName) // ' Upper Limit',&
       elementNbNodes, Model % CurrentElement % NodeIndexes, GotIt)
  IF (.NOT.GotIt) THEN
     temphom = 0.0d00
     WRITE(Message,'(A,A,A,i3,A)') 'No entry for ',TRIM(TempName) // ' Upper Limit',&
          ' found in material no. ', material_id,'. Using 273.16 K.'
     CALL WARN('iceproperties (getViscosityFactor)',Message)
  ELSE
     temphom = MIN(temperature - PressureMeltingPoint(i), 0.0d00)
  END IF
  !----------------------------------------------------
  ! homologous Temperature is below -10 degrees celsius
  !----------------------------------------------------
  IF (temphom < LimitTemp(i)) THEN
     k=1
     !----------------------------------------------------
     ! homologous Temperature is above -10 degrees celsius
     !----------------------------------------------------
  ELSE
     k=2
  END IF
  !----------------------------------------------
  ! Comput the rate factor
  !----------------------------------------------
  rateFactor =&
       arrheniusFactor(k,i)*exp(-1.0D00*activationEnergy(k,i)/(gasconst*(2.7315D02 + temphom)))
                                            
! u, v, w local coord of node i
                u = CurrentElement % TYPE % NodeU(i)
                v = CurrentElement % TYPE % NodeV(i)
                w = CurrentElement % TYPE % NodeW(i)
                                                                                    
               stat = ElementInfo(CurrentElement, ElementNodes, u, v, w, detJ, &
               Basis, dBasisdx, ddBasisddx, .FALSE., .FALSE.)
             
               LGrad = MATMUL( LocalVelo(:,1:elementNbNodes), dBasisdx(1:elementNbNodes,:))
               
               Dij = 0.5_dp * ( LGrad + TRANSPOSE(LGrad) )

!               Temp = TempValues( TempPerm( NodeIndexes(i) ) )
!               m = Expo( i ) 
                
                IF ( ASSOCIATED( FabricVariable ) ) THEN
                a2(1) = FabricValues( 5 * (FabricPerm(NodeIndexes(i))-1) + 1 )
                a2(2) = FabricValues( 5 * (FabricPerm(NodeIndexes(i))-1) + 2 )
                a2(3) = 1.0_dp - a2(1) - a2(2)
                a2(4) = FabricValues( 5 * (FabricPerm(NodeIndexes(i))-1) + 3 )
                a2(5) = FabricValues( 5 * (FabricPerm(NodeIndexes(i))-1) + 4 )
                a2(6) = FabricValues( 5 * (FabricPerm(NodeIndexes(i))-1) + 5 )
!               a2(1) = SUM( Locala11(1:n) * Basis(1:n) )
                ELSE
                  a2(4:6) = 0.0
                  a2(1:3) = 1.0/3.0_dp
                END IF
 
                CALL IBOF( a2, a4 )

!Compute the rest of the components of the a4 orientation tensor

!        a1133 = -a1111+a11-a1122
         a1133 = -a4(1)+a2(1)-a4(3)

!        a2233 = -a2222+a22-a1122
         a2233 = -a4(2)+a2(2)-a4(3)

!        a3333 = a33-a1133-a2233
         a3333 = a2(3)-a1133-a2233
    
!        a3331=a1333=a13-a1311-a1322=a13-a1131-a2231
         a3331 = a2(6)-a4(6)-a4(5)

!        a3332=a2333=a23-a2311-a2322=a23-a1123-a2223
         a3332 = a2(5)-a4(4)-a4(8)

!        a3312=a1233=a12-a1211-a1222=a12-a1112-a2212
         a3312 = a2(4)-a4(7)-a4(9)

!If necessary correct D11, D22, D33 to preserve incompressibility
Eps = Dij(1,1)+Dij(2,2)+Dij(3,3) 
IF (Eps > 0 .OR. Eps < 0 ) THEN
Eps = Eps/3.0D00

Dij(1,1)= Dij(1,1)-Eps
Dij(2,2)= Dij(2,2)-Eps
Dij(3,3)= Dij(3,3)-Eps
END IF 

!        Computes the arguments of the equation of A                   
         AfirstArg = a2(1)*((Dij(1,1))**2.0+(Dij(1,2))**2.0+(Dij(1,3))**2.0) &
                     +a2(2)*((Dij(1,2))**2.0+(Dij(2,2))**2.0+(Dij(2,3))**2.0) &
                     +a2(3)*((Dij(1,3))**2.0+(Dij(2,3))**2.0+(Dij(3,3))**2.0) &
                     +a2(4)*(2.0*Dij(1,1)*Dij(1,2)+2.0*Dij(2,2)*Dij(1,2)+2.0*Dij(2,3)*Dij(1,3)) &
                     +a2(5)*(2.0*Dij(1,3)*Dij(1,2)+2.0*Dij(2,2)*Dij(2,3)+2.0*Dij(3,3)*Dij(2,3)) &
                     +a2(6)*(2.0*Dij(1,1)*Dij(1,3)+2.0*Dij(2,3)*Dij(1,2)+2.0*Dij(3,3)*Dij(1,3))

         AsecondArg = a4(1)*(Dij(1,1))**2.0+a4(2)*(Dij(2,2))**2.0+a3333*(Dij(3,3))**2.0 &
                     +a4(3)*(2*Dij(1,1)*Dij(2,2)+4.0*(Dij(1,2))**2.0) &
                     +a1133*(2*Dij(1,1)*Dij(3,3)+4.0*(Dij(1,3))**2.0) &
                     +a2233*(2*Dij(2,2)*Dij(3,3)+4.0*(Dij(2,3))**2.0) &
                     +4.0*a4(4)*Dij(1,1)*Dij(2,3)+8.0*a4(4)*Dij(1,2)*Dij(1,3) &
                     +4.0*a4(5)*Dij(2,2)*Dij(1,3)+8.0*a4(5)*Dij(1,2)*Dij(2,3) &
                     +4.0*a3312*Dij(3,3)*Dij(1,2)+8.0*a3312*Dij(1,3)*Dij(2,3) &
                     +4.0*a4(7)*Dij(1,1)*Dij(1,2)+4.0*a4(6)*Dij(1,1)*Dij(1,3) &
                     +4.0*a4(9)*Dij(2,2)*Dij(2,1)+4.0*a4(8)*Dij(2,2)*Dij(2,3) &
                     +4.0*a3331*Dij(3,3)*Dij(1,3)+4.0*a3332*Dij(3,3)*Dij(2,3) 
         
                 
!        Compute the 5/tr(D)^2
         Aconst = (Dij(1,1))**2.0+(Dij(2,2))**2.0+(Dij(3,3))**2.0+2.0*(Dij(1,2))**2.0+ &
                   2.0*(Dij(3,1))**2.0+2.0*(Dij(2,3))**2.0
              
         IF (Aconst < MinGamma) Aconst = MinGamma

!        IF (sqrt(2.0*Aconst) < MinGamma) THEN
!        A = 1.0_dp              
!        ELSE
                  
!        Compute the value of A=(5/tr(D)^2)*(D.a2.D-(a4:D):D)
         A = 5.0_dp * (AfirstArg-AsecondArg) / Aconst
!        END IF

!        Compute the Enhancement factor function
         IF (A >= 1.0 ) THEN
             E = (4.0*A**2.0*(aniEnhancementFact(i)-1.0)+25.0-4.0*aniEnhancementFact(i))/21.0
         ELSE
             E = (1.0-Emin)*A**((8.0/21.0)*((aniEnhancementFact(i)-1.0)/(1.0-Emin)))+Emin
         END IF
            
!        Compute the viscosity
         anisoVisFact = 1.0D00 / (( 2.0D00 * E * rateFactor )**viscosityExponent(i))

         EXIT           
     END IF 
  END DO
! write(*,*)A,E,anisoVisFact   

CONTAINS       

!!! Compute fourth order tensor a4 from a2 with closure function IBOF (Chung, 2002)
!!! a2 enters in the order : 11, 22, 33, 12, 23 ,13
!!! Output for a4 is in the order : 1111, 2222, 1122, 1123, 2231, 1131, 1112, 2223, 2212
!!! Code modified from Gillet-Chaullet source
!------------------------------------------------------------------------------
      SUBROUTINE IBOF(a2,a4)

      USE Types
       
       implicit none
       Real(dp),dimension(6),intent(in):: a2  
       Real(dp),dimension(9),intent(out):: a4  
       Real(dp):: a_11,a_22,a_33,a_12,a_13,a_23
       Real(dp):: b_11,b_22,b_12,b_13,b_23
       Real(dp):: aPlusa

       Real(dp),dimension(21) :: vec
       Real(dp),dimension(3,21) :: Mat
       Real(dp),dimension(6) :: beta
       Real(dp) :: Inv2,Inv3
       integer :: i,j      

       
       a_11=a2(1)
       a_22=a2(2)
       a_33=a2(3)
       a_12=a2(4)
       a_23=a2(5)
       a_13=a2(6)


      !Coefficiants 

      Mat(1,1)=0.217774509809788e+02_dp
      Mat(1,2)=-.297570854171128e+03_dp
      Mat(1,3)=0.188686077307885e+04_dp
      Mat(1,4)=-.272941724578513e+03_dp
      Mat(1,5)=0.417148493642195e+03_dp
      Mat(1,6)=0.152038182241196e+04_dp
      Mat(1,7)=-.137643852992708e+04_dp
      Mat(1,8)=-.628895857556395e+03_dp
      Mat(1,9)=-.526081007711996e+04_dp
      Mat(1,10)=-.266096234984017e+03_dp
      Mat(1,11)=-.196278098216953e+04_dp
      Mat(1,12)=-.505266963449819e+03_dp
      Mat(1,13)=-.110483041928547e+03_dp
      Mat(1,14)=0.430488193758786e+04_dp
      Mat(1,15)=-.139197970442470e+02_dp
      Mat(1,16)=-.144351781922013e+04_dp
      Mat(1,17)=-.265701301773249e+03_dp
      Mat(1,18)=-.428821699139210e+02_dp
      Mat(1,19)=-.443236656693991e+01_dp
      Mat(1,20)=0.309742340203200e+04_dp
      Mat(1,21)=0.386473912295113e+00_dp
      Mat(2,1)=-.514850598717222e+00_dp
      Mat(2,2)=0.213316362570669e+02_dp
      Mat(2,3)=-.302865564916568e+03_dp
      Mat(2,4)=-.198569416607029e+02_dp
      Mat(2,5)=-.460306750911640e+02_dp
      Mat(2,6)=0.270825710321281e+01_dp
      Mat(2,7)=0.184510695601404e+03_dp
      Mat(2,8)=0.156537424620061e+03_dp
      Mat(2,9)=0.190613131168980e+04_dp
      Mat(2,10)=0.277006550460850e+03_dp
      Mat(2,11)=-.568117055198608e+02_dp
      Mat(2,12)=0.428921546783467e+03_dp
      Mat(2,13)=0.142494945404341e+03_dp
      Mat(2,14)=-.541945228489881e+04_dp
      Mat(2,15)=0.233351898912768e+02_dp
      Mat(2,16)=0.104183218654671e+04_dp
      Mat(2,17)=0.331489412844667e+03_dp
      Mat(2,18)=0.660002154209991e+02_dp
      Mat(2,19)=0.997500770521877e+01_dp
      Mat(2,20)=0.560508628472486e+04_dp
      Mat(2,21)=0.209909225990756e+01_dp
      Mat(3,1)=0.203814051719994e+02_dp
      Mat(3,2)=-.283958093739548e+03_dp
      Mat(3,3)=0.173908241235198e+04_dp
      Mat(3,4)=-.195566197110461e+03_dp
      Mat(3,5)=-.138012943339611e+03_dp
      Mat(3,6)=0.523629892715050e+03_dp
      Mat(3,7)=0.859266451736379e+03_dp
      Mat(3,8)=-.805606471979730e+02_dp
      Mat(3,9)=-.468711180560599e+04_dp
      Mat(3,10)=0.889580760829066e+01_dp
      Mat(3,11)=-.782994158054881e+02_dp
      Mat(3,12)=-.437214580089117e+02_dp
      Mat(3,13)=0.112996386047623e+01_dp
      Mat(3,14)=0.401746416262936e+04_dp
      Mat(3,15)=0.104927789918320e+01_dp
      Mat(3,16)=-.139340154288711e+03_dp
      Mat(3,17)=-.170995948015951e+02_dp
      Mat(3,18)=0.545784716783902e+00_dp
      Mat(3,19)=0.971126767581517e+00_dp
      Mat(3,20)=0.141909512967882e+04_dp
      Mat(3,21)=0.994142892628410e+00_dp

       
      ! Compute the invariants
      Inv2=0.5_dp*(1._dp-(a_11*a_11+a_22*a_22+a_33*a_33+ &
            2._dp*(a_12*a_12+a_13*a_13+a_23*a_23)))
            
       Inv3=a_11*(a_22*a_33-a_23*a_23)+a_12*(a_23*a_13-a_12*a_33)+ &
             a_13*(a_12*a_23-a_22*a_13)
       
     ! complete polynome of degree 5 for the 2 invariants.
         vec(1)=1._dp
         vec(2)=Inv2
         vec(3)=vec(2)*vec(2)
         vec(4)=Inv3
         vec(5)=vec(4)*vec(4)
         vec(6)=vec(2)*vec(4)
         vec(7)=vec(3)*vec(4)
         vec(8)=vec(2)*vec(5)
         vec(9)=vec(2)*vec(3)
         vec(10)=vec(5)*vec(4)
         vec(11)=vec(9)*vec(4)
         vec(12)=vec(3)*vec(5)
         vec(13)=vec(2)*vec(10)
         vec(14)=vec(3)*vec(3)
         vec(15)=vec(5)*vec(5)
         vec(16)=vec(14)*vec(4)
         vec(17)=vec(12)*vec(2)
         vec(18)=vec(12)*vec(4)
         vec(19)=vec(2)*vec(15)
         vec(20)=vec(14)*vec(2)
         vec(21)=vec(15)*vec(4)

       ! Compites beta_bar (cf annexe C Chung)
       ! Warning: beta(1)=beta_bar_3 (Chung); beta(2)=beta_bar_4; beta(3)=beta_bar_6
       !           beta(4)=beta_bar_1        ; beta(5)=beta_bar_2; beta(6)=beta_bar_5

       ! calcul the three betas in terms of the polynomes
         beta(:)=0._dp
         Do i=1,3
          Do j=1,21
            beta(i)=beta(i)+Mat(i,j)*vec(j)
          End do
         End do
          
       ! calcul the other 3 to get the normalisation
         beta(4)=3._dp*(-1._dp/7._dp+beta(1)*(1._dp/7._dp+4._dp*Inv2/7._dp+8._dp*Inv3/3._dp)/5._dp- &
                  beta(2)*(0.2_dp-8._dp*Inv2/15._dp-14._dp*Inv3/15._dp)- &
                  beta(3)*(1._dp/35._dp-24._dp*Inv3/105._dp-4._dp*Inv2/35._dp+ &
                  16._dp*Inv2*Inv3/15._dp+8._dp*Inv2*Inv2/35._dp))/5._dp

         beta(5)=6._dp*(1._dp-0.2_dp*beta(1)*(1._dp+4._dp*Inv2)+ &
                  7._dp*beta(2)*(1._dp/6._dp-Inv2)/5._dp- &
                  beta(3)*(-0.2_dp+2._dp*Inv3/3._dp+4._dp*Inv2/5._dp- &
                  8._dp*Inv2*Inv2/5._dp))/7._dp

         beta(6)=-4._dp*beta(1)/5._dp-7._dp*beta(2)/5._dp- &
                   6._dp*beta(3)*(1._dp-4._dp*Inv2/3._dp)/5._dp

        !beta_bar
        Do i=1,6
         beta(i)=beta(i)/3._dp
        End do
         beta(2)=beta(2)/2._dp
         beta(5)=beta(5)/2._dp
         beta(6)=beta(6)/2._dp

        !! Compute 5 b=a.a
        b_11=a_11*a_11+a_12*a_12+a_13*a_13
        b_22=a_22*a_22+a_12*a_12+a_23*a_23
        b_12=a_11*a_12+a_12*a_22+a_13*a_23
        b_13=a_11*a_13+a_12*a_23+a_13*a_33
        b_23=a_12*a_13+a_22*a_23+a_23*a_33

        !Compute the  9 terms of a4

        a4(1)=3._dp*beta(4)+6._dp*beta(5)*a_11+3._dp*beta(1)*a_11*a_11+&
         6._dp*beta(2)*b_11+6._dp*beta(6)*a_11*b_11+3._dp*beta(3)*b_11*b_11
        a4(2)=3._dp*beta(4)+6._dp*beta(5)*a_22+3._dp*beta(1)*a_22*a_22+&
         6._dp*beta(2)*b_22+6._dp*beta(6)*a_22*b_22+3._dp*beta(3)*b_22*b_22

        a4(3)=beta(4)+beta(5)*(a_22+a_11)+beta(1)*(a_11*a_22+2._dp*a_12*a_12)+&
         beta(2)*(b_22+b_11)+beta(6)*(a_11*b_22+a_22*b_11+4._dp*a_12*b_12)+&
         beta(3)*(b_11*b_22+2._dp*b_12*b_12)


         a4(4)=beta(5)*a_23+beta(1)*(a_11*a_23+2._dp*a_12*a_13)+beta(2)*b_23+&
          beta(6)*(a_11*b_23+a_23*b_11+2._dp*(a_12*b_13+a_13*b_12))+beta(3)*&
          (b_11*b_23+2._dp*b_12*b_13)
         a4(5)=beta(5)*a_13+beta(1)*(a_22*a_13+2._dp*a_12*a_23)+beta(2)*b_13+&
          beta(6)*(a_22*b_13+a_13*b_22+2._dp*(a_12*b_23+a_23*b_12))+beta(3)*&
          (b_22*b_13+2._dp*b_12*b_23)


         a4(6)=3._dp*beta(5)*a_13+3._dp*beta(1)*a_11*a_13+3._dp*beta(2)*b_13+&
          3._dp*beta(6)*(a_11*b_13+a_13*b_11)+3._dp*beta(3)*b_11*b_13
         a4(7)=3._dp*beta(5)*a_12+3._dp*beta(1)*a_11*a_12+3._dp*beta(2)*b_12+&
          3._dp*beta(6)*(a_11*b_12+a_12*b_11)+3._dp*beta(3)*b_11*b_12
         a4(8)=3._dp*beta(5)*a_23+3._dp*beta(1)*a_22*a_23+3._dp*beta(2)*b_23+&
          3._dp*beta(6)*(a_22*b_23+a_23*b_22)+3._dp*beta(3)*b_22*b_23
         a4(9)=3._dp*beta(5)*a_12+3._dp*beta(1)*a_22*a_12+3._dp*beta(2)*b_12+&
          3._dp*beta(6)*(a_22*b_12+a_12*b_22)+3._dp*beta(3)*b_22*b_12

         END SUBROUTINE IBOF
!**************************************************************************************

END FUNCTION caffeGetViscosity
