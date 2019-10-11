!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation; either
! * version 2.1 of the License, or (at your option) any later version.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! * 
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library (in file ../LGPL-2.1); if not, write 
! * to the Free Software Foundation, Inc., 51 Franklin Street, 
! * Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Authors: Juha Ruokolainen, Thomas Zwinger
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 24 Apr 1997
! *
! *****************************************************************************/

!> \ingroup ElmerLib
!> \}

!---------------------------------------------------------------------
!> Module for material models of fluids mainly.
!---------------------------------------------------------------------


MODULE MaterialModels
 
   USE DefUtils

   IMPLICIT NONE
   INTEGER, PARAMETER :: Incompressible = 0, UserDefined1 = 1,UserDefined2 = 2
   INTEGER, PARAMETER :: PerfectGas1 = 3, PerfectGas2 = 4, PerfectGas3 = 5, Thermal = 6


CONTAINS


 
!------------------------------------------------------------------------------
!> Return second invariant.
!> Note: Actually SQUARE of the second invariant of velocity is returned
!------------------------------------------------------------------------------
   FUNCTION SecondInvariant( Velo,dVelodx,CtrMetric,Symb ) RESULT(SecInv)
!------------------------------------------------------------------------------

     REAL(KIND=dp), OPTIONAL :: CtrMetric(3,3),Symb(3,3,3)
     REAL(KIND=dp) :: Velo(3),dVelodx(3,3),SecInv

!------------------------------------------------------------------------------

     INTEGER :: i,j,k,l
     REAL(KIND=dp) :: CovMetric(3,3),s,t

     SecInv = 0.0D0

     IF ( CurrentCoordinateSystem() == Cartesian ) THEN
!------------------------------------------------------------------------------

       DO i=1,3
         DO j=1,3
           s = dVelodx(i,j) + dVelodx(j,i)
           SecInv = SecInv + s * s
         END DO
       END DO
!------------------------------------------------------------------------------
     ELSE IF ( CurrentCoordinateSystem() == AxisSymmetric ) THEN

        SecInv = (2*dVelodx(1,1))**2 + (2*dVelodx(2,2))**2 + &
          2*(dVelodx(1,2) + dVelodx(2,1))**2 + (2*Velo(1)*symb(1,3,3))**2

     ELSE

!------------------------------------------------------------------------------
       CovMetric = CtrMetric
       CALL InvertMatrix( CovMetric,3 )

       DO i=1,3
         DO j=1,3
            s = 0.0d0
            t = 0.0d0

            DO k=1,3
               s = s + CovMetric(i,k) * dVelodx(k,j) + &
                       CovMetric(j,k) * dVelodx(k,i)

              t = t + CtrMetric(j,k) * dVelodx(i,k) + &
                      CtrMetric(i,k) * dVelodx(j,k)

              DO l=1,3
                s = s - CovMetric(i,k) * Symb(l,j,k) * Velo(l)
                s = s - CovMetric(j,k) * Symb(l,i,k) * Velo(l)

                t = t - CtrMetric(j,k) * Symb(l,k,i) * Velo(l)
                t = t - CtrMetric(i,k) * Symb(l,k,j) * Velo(l)
              END DO
           END DO
           SecInv = SecInv + s * t
         END DO
       END DO
!------------------------------------------------------------------------------

     END IF
!------------------------------------------------------------------------------
   END FUNCTION SecondInvariant
!------------------------------------------------------------------------------



#if 0
this ise not in USE
!------------------------------------------------------------------------------
   SUBROUTINE FrictionHeat( Heat,Viscosity,Ux,Uy,Uz,Element,Nodes )
!------------------------------------------------------------------------------
     REAL(KIND=dp)  :: Heat(:),Viscosity(:),Ux(:),Uy(:),Uz(:)
     TYPE(Nodes_t)     :: Nodes
     TYPE(Element_t)   :: Element
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Basis(MAX_ELEMENT_NODES),dBasisdx(MAX_ELEMENT_NODES,3)
     REAL(KIND=dp) :: s,u,v,w,SqrtMetric,SqrtElementMetric,Velo(3)
     REAL(KIND=dp) :: Metric(3,3),dVelodx(3,3), &
                           CtrMetric(3,3),Symb(3,3,3),dSymb(3,3,3,3)

     LOGICAL :: stat
     INTEGER :: i,j,n
!------------------------------------------------------------------------------

     n = Element % TYPE % NumberOfNodes
     DO i=1,n
        u = Element % TYPE % NodeU(i)
        v = Element % TYPE % NodeV(i)
        w = Element % TYPE % NodeW(i)
!------------------------------------------------------------------------------
!       Basis function values & derivatives at the calculation point
!------------------------------------------------------------------------------
        stat = ElementInfo( Element,Nodes, u, v, w, &
          SqrtElementMetric, Basis,dBasisdx )
!------------------------------------------------------------------------------
!       Coordinate system dependent information
!------------------------------------------------------------------------------
        CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb, &
               Nodes % x(i),Nodes % y(i),Nodes % z(i) )
!------------------------------------------------------------------------------

        DO j=1,3
          dVelodx(1,j) = SUM( Ux(1:n) * dBasisdx(1:n,j) )
          dVelodx(2,j) = SUM( Uy(1:n) * dBasisdx(1:n,j) )
          dVelodx(3,j) = SUM( Uz(1:n) * dBasisdx(1:n,j) )
        END DO

        Velo(1) = Ux(i)
        Velo(2) = Uy(i)
        Velo(3) = Uz(i)
        Heat(i) = 0.5d0*Viscosity(i)*SecondInvariant(Velo,dVelodx,Metric,Symb)
     END DO

!------------------------------------------------------------------------------
   END SUBROUTINE FrictionHeat
!------------------------------------------------------------------------------
#endif



!------------------------------------------------------------------------------
!> Returns effective viscosity for Navier-Stokes equation. 
!> The viscosity model may be either some nonnewtonian material law, 
!> or from turbulence models, but not from both at the same time.
!------------------------------------------------------------------------------
   FUNCTION EffectiveViscosity( Viscosity,Density,Ux,Uy,Uz,Element, &
        Nodes,n,nd,u,v,w, muder, LocalIP ) RESULT(mu)
     !------------------------------------------------------------------------------

     USE ModelDescription

     REAL(KIND=dp)  :: Viscosity,Density,u,v,w,mu,Ux(:),Uy(:),Uz(:)
     REAL(KIND=dp), OPTIONAL :: muder
     TYPE(Nodes_t)  :: Nodes
     INTEGER :: n,nd
     INTEGER, OPTIONAL :: LocalIP
     TYPE(Element_t),POINTER :: Element

     !------------------------------------------------------------------------------
     REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3)
     REAL(KIND=dp) :: ss,s,SqrtMetric,SqrtElementMetric,Velo(3)
     REAL(KIND=dp) :: Metric(3,3), dVelodx(3,3), CtrMetric(3,3), &
          Symb(3,3,3), dSymb(3,3,3,3)

     INTEGER :: i,j,k
     LOGICAL :: stat,GotIt,UseEUsrf=.FALSE.

     CHARACTER(LEN=MAX_NAME_LEN) :: ViscosityFlag, TemperatureName, EnhcmntFactFlag
     TYPE(ValueList_t), POINTER :: Material
     REAL(KIND=dp) :: x, y, z, c1n(n), c2n(n), c3n(n), c4n(n), &
          c1, c2, c3, c4, c5, c6, c7, Temp, NodalTemperature(n), Tlimit, TempCoeff, &
          h, A1, A2, Q1, Q2, R, NodalEhF(n), EhF, ArrheniusFactor

     ! Temperature is needed for thermal models
     TYPE(Variable_t), POINTER :: TempSol 
     REAL(KIND=dp), POINTER :: Temperature(:)
     INTEGER, POINTER :: TempPerm(:)

     INTEGER(KIND=AddrInt) :: Fnc

     TYPE(Variable_t), POINTER :: Var

     REAL(KIND=dp) :: dist,F2,F3
     REAL(KIND=dp) :: KE_K, KE_E, KE_Z, CT, TimeScale,Clip, Cmu, Vals(n)

     CHARACTER(LEN=MAX_NAME_LEN) :: str

     LOGICAL :: SetArrheniusFactor=.FALSE.

#ifndef USE_ISO_C_BINDINGS
     INTERFACE
        FUNCTION MaterialUserFunction( Proc,Model,Element,Nodes,n,nd, &
             Basis,dBasisdx,Viscosity,Velo, dVelodx ) RESULT(s)
          USE Types
          INTEGER(KIND=AddrInt) :: Proc
          TYPE(Model_t) :: Model
          TYPE(Nodes_t) :: Nodes
          TYPE(Element_t), POINTER :: Element
          INTEGER :: n,nd
          REAL(KIND=dp) :: Basis(:),dBasisdx(:,:),Viscosity, &
               Velo(:), dVelodx(:,:), s
        END FUNCTION MaterialUserFunction
        FUNCTION EnhancementFactorUserFunction( Proc,Model,Element,Nodes,n,nd, &
             Basis,dBasisdx,Viscosity,Velo,dVelodx,SecondInvariantSqr,LocalIP ) RESULT(Ehf)
          USE Types
          INTEGER(KIND=AddrInt) :: Proc
          TYPE(Model_t) :: Model
          TYPE(Nodes_t) :: Nodes
          TYPE(Element_t), POINTER :: Element
          INTEGER :: n,nd,LocalIP 
          REAL(KIND=dp) :: Basis(:),dBasisdx(:,:),Viscosity, &
               Velo(:), dVelodx(:,:), SecondInvariantSqr, Ehf
        END FUNCTION EnhancementFactorUserFunction
     END INTERFACE     
#endif
     !------------------------------------------------------------------------------
     mu = Viscosity
     IF ( PRESENT(muder) ) muder=0

     k = ListGetInteger( CurrentModel % Bodies(Element % BodyId) % Values, 'Material', &
          minv=1, maxv=CurrentModel % NumberOFMaterials )

     Material => CurrentModel % Materials(k) % Values

     ViscosityFlag = ListGetString( Material,'Viscosity Model', GotIt)
     
     IF(.NOT. gotIt) RETURN
     !------------------------------------------------------------------------------
     !    Basis function values & derivatives at the calculation point
     !------------------------------------------------------------------------------
     stat = ElementInfo( Element,Nodes,u,v,w, &
          SqrtElementMetric, Basis,dBasisdx )
     !------------------------------------------------------------------------------
     !   Coordinate system dependent information
     !------------------------------------------------------------------------------
     x = SUM( Nodes % x(1:n) * Basis(1:n) )
     y = SUM( Nodes % y(1:n) * Basis(1:n) )
     z = SUM( Nodes % z(1:n) * Basis(1:n) )
     CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,x,y,z )
     !------------------------------------------------------------------------------
     DO j=1,3
        dVelodx(1,j) = SUM( Ux(1:nd)*dBasisdx(1:nd,j) )
        dVelodx(2,j) = SUM( Uy(1:nd)*dBasisdx(1:nd,j) )
        dVelodx(3,j) = SUM( Uz(1:nd)*dBasisdx(1:nd,j) )
     END DO

     Velo(1) = SUM( Basis(1:nd) * Ux(1:nd) )
     Velo(2) = SUM( Basis(1:nd) * Uy(1:nd) )
     Velo(3) = SUM( Basis(1:nd) * Uz(1:nd) )

     ! This is the square of shearrate which results to 1/2 in exponent 
     ! Also the derivative is taken with respect to the square
     !-------------------------------------------------------------------
     ss = 0.5_dp * SecondInvariant(Velo,dVelodx,Metric,Symb)


     SELECT CASE( ViscosityFlag )


     CASE('glen')
        c2n = ListGetReal( Material, 'Glen Exponent', n, Element % NodeIndexes, GotIt ) ! this is the real exponent, n, not 1/n
        IF (.NOT.GotIt) c2n(1:n) = 3.0_dp
        c2 = SUM( Basis(1:n) * c2n(1:n) )
        s = ss/4.0_dp ! the second invariant is not taken from the strain rate tensor, but rather 2*strain rate tensor (that's why we divide by 4 = 2**2)        
        
        SetArrheniusFactor = GetLogical(Material, 'Set Arrhenius Factor', GotIt)
        IF ( (.NOT.GotIt) .OR. .NOT.(SetArrheniusFactor)) THEN
           NodalTemperature(1:n) = ListGetReal(Material, 'Constant Temperature', n, Element % NodeIndexes, GotIt) !we are happy as is
           IF(.NOT.GotIt) THEN !we have to find a temperature field

              TemperatureName = GetString(Material, 'Temperature Field Variable', GotIt)
              IF (.NOT.GotIt) WRITE(TemperatureName,'(A)') 'Temperature'
              TempSol => VariableGet( CurrentModel % Variables,TRIM(TemperatureName))
              IF ( ASSOCIATED( TempSol) ) THEN
                 TempPerm    => TempSol % Perm
                 Temperature => TempSol % Values   
                 Temp =  SUM(Basis(1:n) * Temperature(TempPerm(Element % NodeIndexes(1:n))))
              ELSE
                 WRITE(Message, '(A,A,A)') 'Could not find variable ',&
                      TRIM(TemperatureName),' to inquire temperature field for Glen'
                 CALL FATAL('EffectiveViscosity',Message)
              END IF
           
           ELSE
              Temp = SUM(Basis(1:n) * NodalTemperature(1:n))
           END IF
        
           R = GetConstReal( CurrentModel % Constants,'Gas Constant',GotIt)
           IF (.NOT.GotIt) R = 8.314_dp
           ! lets for the time being have this hardcoded
           Tlimit = GetConstReal(Material, 'Limit Temperature', GotIt)
           IF (.NOT.GotIt) THEN
              Tlimit = -10.0_dp
              CALL INFO('EffectiveViscosity','Limit Temperature not found. Setting to -10', Level=5)
           END IF
           A1 = GetConstReal(Material, 'Rate Factor 1', GotIt)
           IF (.NOT.GotIt) THEN
              A1 = 3.985d-13
              CALL INFO('EffectiveViscosity','Rate Factor 1 not found. Setting to 3.985e-13', Level=5)
           END IF
           A2 = GetConstReal(Material, 'Rate Factor 2', GotIt)
           IF (.NOT.GotIt) THEN
              A2 = 1.916d03
              CALL INFO('EffectiveViscosity','Rate Factor 2 not found. Setting to 1.916E03', Level=5)
           END IF
           Q1 = GetConstReal(Material, 'Activation Energy 1', GotIt)
           IF (.NOT.GotIt) THEN
              Q1 = 60.0d03
              CALL INFO('EffectiveViscosity','Activation Energy 1 not found. Setting to 60.0E03', Level=5)
           END IF
           Q2 = GetConstReal(Material, 'Activation Energy 2', GotIt)
           IF (.NOT.GotIt) THEN
              Q2 = 139.0d03
              CALL INFO('EffectiveViscosity','Activation Energy 2 not found. Setting to 139.0d03', Level=5)
           END IF
        
           IF (Temp.LE. Tlimit) THEN
              ArrheniusFactor = A1 * EXP( -Q1/(R * (273.15_dp + Temp)))
           ELSE IF((Tlimit<Temp) .AND. (Temp .LE. 0.0_dp)) THEN
              ArrheniusFactor = A2 * EXP( -Q2/(R * (273.15_dp + Temp)))
           ELSE
              ArrheniusFactor = A2 * EXP( -Q2/(R * (273.15_dp)))
              CALL INFO('EffectiveViscosity',&
                   'Positive Temperature detected in Glen - limiting to zero!', Level = 5)
           END IF
        ELSE
          ArrheniusFactor = GetConstReal(Material,'Arrhenius Factor', GotIt)
          IF (.NOT.(GotIt)) THEN 
            CALL FATAL('EffectiveViscosity',&
                 '<Set Arrhenius Factor> is TRUE, but no value <Arrhenius Factor> found')
          END IF
        END IF
        Ehf = 1.0_dp
        EnhcmntFactFlag = ListGetString( Material,'Glen Enhancement Factor Function', UseEUsrf )
        IF (UseEUsrf) THEN
          IF (.NOT.PRESENT(LocalIP)) CALL FATAL('EffectiveViscosity',&
               'Chose "Glen Enhancement Factor Function" but no LocalIP provided by calling function')
          Fnc = GetProcAddr( EnhcmntFactFlag, Quiet=.TRUE. )
          EhF = EnhancementFactorUserFunction( Fnc, CurrentModel, Element, Nodes, n, nd, &
               Basis, dBasisdx, Viscosity, Velo, dVelodx, s , LocalIP)
        ELSE
          NodalEhF(1:n) =  ListGetReal( Material, 'Glen Enhancement Factor',&
               n, Element % NodeIndexes, GotIt )
          IF (GotIt) &
               EhF = SUM(Basis(1:n) * NodalEhF(1:n))
        END IF
        
        IF (PRESENT(muder)) muder = 0.5_dp * (  EhF * ArrheniusFactor)**(-1.0_dp/c2) &
             * ((1.0_dp/c2)-1.0_dp)/2.0_dp * s**(((1.0_dp/c2)-1.0_dp)/2.0_dp - 1.0_dp)/4.0_dp

        c3n = ListGetReal( Material, 'Critical Shear Rate',n, Element % NodeIndexes,GotIt )
        IF (GotIt) THEN
           c3 = SUM( Basis(1:n) * c3n(1:n) )
           IF(s < c3**2) THEN
              s = c3**2
              IF (PRESENT(muder)) muder = 0._dp
           END IF
        END IF

        ! compute the effective viscosity
        mu = 0.5_dp * (EhF * ArrheniusFactor)**(-1.0_dp/c2) * s**(((1.0_dp/c2)-1.0_dp)/2.0_dp);

     CASE('power law')
        c2n = ListGetReal( Material, 'Viscosity Exponent', n, Element % NodeIndexes )
        c2 = SUM( Basis(1:n) * c2n(1:n) )

        s = ss
        IF (PRESENT(muder)) THEN
           IF (s /= 0) THEN
              muder = Viscosity * (c2-1)/2 * s**((c2-1)/2-1)
           ELSE
              muder = 0.0_dp
           END IF
        END IF

        c3n = ListGetReal( Material, 'Critical Shear Rate',n, Element % NodeIndexes,gotIt )
        IF (GotIt) THEN
           c3 = SUM( Basis(1:n) * c3n(1:n) )
           IF(s < c3**2) THEN
              s = c3**2
              IF (PRESENT(muder)) muder = 0._dp
           END IF
        END IF
        mu = Viscosity * s**((c2-1)/2)

        c4n = ListGetReal( Material, 'Nominal Shear Rate',n, Element % NodeIndexes,gotIt )
        IF (GotIt) THEN
           c4 = SUM( Basis(1:n) * c4n(1:n) )
           mu = mu / c4**(c2-1)
           IF (PRESENT(muder)) muder = muder / c4**(c2-1)
        END IF

     CASE('power law too')
        c2n = ListGetReal( Material, 'Viscosity Exponent', n, Element % NodeIndexes )
        c2 = SUM( Basis(1:n) * c2n(1:n) )
        mu = Viscosity **(-1/c2)* ss**(-(c2-1)/(2*c2)) / 2
        IF ( PRESENT(muder) )  muder = &
             Viscosity**(-1/c2)*(-(c2-1)/(2*c2))*ss*(-(c2-1)/(2*c2)-1) / 2

     CASE ('carreau')
        c1n = ListGetReal( Material, 'Viscosity Difference',n,Element % NodeIndexes )
        c1 = SUM( Basis(1:n) * c1n(1:n) )
        c2n = ListGetReal( Material, 'Viscosity Exponent', n, Element % NodeIndexes )
        c2 = SUM( Basis(1:n) * c2n(1:n) )
        c3n = ListGetReal( Material, 'Viscosity Transition',n,Element % NodeIndexes )
        c3 = SUM( Basis(1:n) * c3n(1:n) )
        c4 = ListGetConstReal( Material, 'Yasuda Exponent',gotIt)
        IF(gotIt) THEN
           s = SQRT(ss)
           mu = Viscosity + c1 * (1 + c3**c4*ss**(c4/2))**((c2-1)/c4) 
           IF ( PRESENT(muder ) ) muder =  &
                c1*(1+c3**c4*ss**(c4/2))**((c2-1)/c4-1)*(c2-1)/2*c3**c4*ss**(c4/2-1)
        ELSE
           mu = Viscosity + c1 * (1 + c3*c3*ss)**((c2-1)/2) 
           IF ( PRESENT(muder) ) muder = &
                c1*(c2-1)/2*c3**2*(1+c3**2*ss)**((c2-1)/2-1)
        END IF

     CASE ('cross')
        c1n = ListGetReal( Material, 'Viscosity Difference',n,Element % NodeIndexes )
        c1 = SUM( Basis(1:n) * c1n(1:n) )
        c2n = ListGetReal( Material, 'Viscosity Exponent', n, Element % NodeIndexes )
        c2 = SUM( Basis(1:n) * c2n(1:n) )
        c3n = ListGetReal( Material, 'Viscosity Transition',n,Element % NodeIndexes )
        c3 = SUM( Basis(1:n) * c3n(1:n) )
        mu = Viscosity + c1 / (1 + c3*ss**(c2/2))
        IF ( PRESENT(muder) ) muder = &
             -c1*c3*ss**(c2/2)*c2 / (2*(1+c3*ss**(c2/2))**2*ss)

     CASE ('powell eyring')
        c1n = ListGetReal( Material, 'Viscosity Difference',n,Element % NodeIndexes)
        c1 = SUM( Basis(1:n) * c1n(1:n) )
        c2 = ListGetConstReal( Material, 'Viscosity Transition')
        s = SQRT(ss)
        IF(c2*s < 1.0d-5) THEN
           mu = Viscosity + c1
        ELSE
           mu = Viscosity + c1 * LOG(c2*s+SQRT(c2*c2*ss+1))/(c2*s)
           IF ( PRESENT(muder) ) muder = &
                c1*(c2/(2*s)+c2**2/(2*SQRT(c2**2*ss+1)))/((c2*s+SQRT(c2*ss+1))*c2*s) - &
                c1*LOG(c2*s+SQRT(c2**2*ss+1))/(c2*s**3)/2
        END IF

     CASE( 'smagorinsky' )
        c2n = ListGetReal( Material, 'Smagorinsky Constant', &
             n, Element % NodeIndexes,gotit )
        c2 = SUM( Basis(1:n) * c2n(1:n) )
        h  = ElementDiameter( Element, Nodes )
        Viscosity = Viscosity + Density * c2 * h**2 * SQRT(2*ss) / 2
        IF ( PRESENT(muder) ) muder = &
             Density*c2*h**2*SQRT(2._dp)/(4*SQRT(ss))

     CASE( 'ke','k-epsilon' )
        IF (ListGetString(Material,'KE Model',gotIt)/='v2-f' ) THEN
           Var => VariableGet( CurrentModel % Variables, 'Kinetic Energy' )
           IF ( .NOT. ASSOCIATED( Var ) ) &
                CALL Fatal( 'Viscosity Model', 'The kinetic energy variable not defined?' )
           KE_K = SUM(Basis(1:n) * Var % Values(Var % Perm(Element % NodeIndexes)))

           Var => VariableGet( CurrentModel % Variables, 'Kinetic Dissipation' )
           IF ( .NOT. ASSOCIATED( Var ) ) &
                CALL Fatal( 'Viscosity Model', 'The kinetic dissipation rate variable not defined?' )
           KE_E = SUM(Basis(1:n) * Var % Values(Var % Perm(Element % NodeIndexes)))

           Vals(1:n) = ListGetReal( Material, 'KE Cmu',n,Element % NodeIndexes,gotIt )
           IF ( .NOT. GotIt ) THEN
              Cmu = SUM( Basis(1:n) * Vals(1:n) )
           ELSE
              Cmu = 0.09_dp 
           END IF
           mu = Viscosity + Cmu*Density*KE_K**2 / KE_E
        ELSE
           Var => VariableGet( CurrentModel % Variables, 'Kinetic Energy' )
           IF ( .NOT. ASSOCIATED( Var ) ) &
                CALL Fatal( 'Viscosity Model', 'The kinetic energy variable not defined?' )
           KE_K = SUM(Basis(1:n) * Var % Values(Var % Perm(Element % NodeIndexes)))

           Var => VariableGet( CurrentModel % Variables, 'Kinetic Dissipation' )
           IF ( .NOT. ASSOCIATED( Var ) ) &
                CALL Fatal( 'Viscosity Model', 'The kinetic dissipation rate variable not defined?' )
           KE_E = SUM(Basis(1:n) * Var % Values(Var % Perm(Element % NodeIndexes)))

           Var => VariableGet( CurrentModel % Variables, 'V2' )
           IF ( .NOT. ASSOCIATED( Var ) ) &
                CALL Fatal( 'Viscosity Model', 'The V2 variable not defined?' )
           KE_Z = SUM(Basis(1:n) * Var % Values(Var % Perm(Element % NodeIndexes)))

           Vals(1:n) = ListGetReal( Material, 'V2-F CT',n,Element % NodeIndexes )
           CT = SUM( Basis(1:n) * Vals(1:n) )
           TimeScale = MAX( KE_K/KE_E, CT*SQRT(Viscosity/Density/KE_E) )

           Vals(1:n) = ListGetReal( Material, 'KE Cmu',n,Element % NodeIndexes )
           Cmu = SUM( Basis(1:n) * Vals(1:n) )

           mu = Viscosity + Cmu*Density*KE_Z*TimeScale
        END IF

     CASE( 'rng k-epsilon' )
        Var => VariableGet( CurrentModel % Variables, 'Effective Viscosity')
        mu = SUM( Basis(1:n) * Var % Values( Var % Perm( Element % NodeIndexes )))

     CASE( 'spalart-allmaras' )
        Var => VariableGet( CurrentModel % Variables, 'Turbulent Viscosity')
        IF ( .NOT. ASSOCIATED( Var ) ) &
             CALL Fatal( 'Viscosity Model', 'The turbulent viscosity variable not defined?' )
        mu = SUM( Basis(1:n) * Var % Values( Var % Perm( Element % NodeIndexes )))
        c1 = mu/(Viscosity/Density)
        c1 = c1**3 / (c1**3 + 7.1_dp**3) 
        mu = Viscosity + mu*Density*c1

     CASE( 'k-omega' )
        Var => VariableGet( CurrentModel % Variables, 'Kinetic Energy' )
        IF ( .NOT. ASSOCIATED( Var ) ) &
             CALL Fatal( 'Viscosity Model', 'The kinetic energy variable not defined?' )
        KE_K = SUM(Basis(1:n) * Var % Values(Var % Perm(Element % NodeIndexes)))

        Var => VariableGet( CurrentModel % Variables, 'Kinetic Dissipation' )
        IF ( .NOT. ASSOCIATED( Var ) ) &
             CALL Fatal( 'Viscosity Model', 'The kinetic dissipation rate variable not defined?' )
        KE_E = SUM(Basis(1:n) * Var % Values(Var % Perm(Element % NodeIndexes)))

        mu = Viscosity + Density * KE_K / KE_E

     CASE( 'sst k-omega' )
        Var => VariableGet( CurrentModel % Variables, 'Kinetic Energy' )
        IF ( .NOT. ASSOCIATED( Var ) ) &
             CALL Fatal( 'Viscosity Model', 'The kinetic energy variable not defined?' )
        KE_K = SUM(Basis(1:n) * Var % Values(Var % Perm(Element % NodeIndexes)))

        Var => VariableGet( CurrentModel % Variables, 'Kinetic Dissipation' )
        IF ( .NOT. ASSOCIATED( Var ) ) &
             CALL Fatal( 'Viscosity Model', 'The kinetic dissipation rate variable not defined?' )
        KE_E = SUM(Basis(1:n) * Var % Values(Var % Perm(Element % NodeIndexes)))

        Var => VariableGet( CurrentModel % Variables, 'Wall distance' )
        IF ( .NOT. ASSOCIATED( Var ) ) &
             CALL Fatal( 'Viscosity Model', 'The wall distance variable not defined?' )
        Dist = SUM(Basis(1:n) * Var % Values(Var % Perm(Element % NodeIndexes)))

        F2 = TANH( MAX(2*SQRT(KE_K)/(0.09_dp*KE_E*Dist), &
             500._dp*Viscosity/(Density*KE_E*Dist**2))**2)

        !        F3 = 1-TANH((150*Viscosity/Density/KE_E/Dist**2)**4)
        F3 = 1

        mu = Viscosity+0.31_dp*Density*KE_K/MAX(0.31_dp*KE_E,SQRT(ss)*F2*F3)

     CASE( 'levelset' )
        TempSol => VariableGet( CurrentModel % Variables, 'Surface' )
        IF ( ASSOCIATED( TempSol) ) THEN
           TempPerm    => TempSol % Perm
           Temperature => TempSol % Values
        ELSE
           CALL Warn('EffectiveViscosity','variable Surface needed in levelset viscosity model')
        END IF

        c1n = ListGetReal( Material, 'Viscosity Difference',n,Element % NodeIndexes)
        c1 = SUM( Basis(1:n) * c1n(1:n) )
        c2 = ListGetConstReal( Material, 'Levelset bandwidth')

        Temp = SUM(Basis(1:n) * Temperature(TempPerm(Element % NodeIndexes(1:n))))
        Temp = Temp / c2
        IF(Temp < -1.0) THEN
           c3 = 0.0d0
        ELSE IF(Temp > 1.0) THEN
           c3 = 1.0d0
        ELSE
           c3 = 0.75d0 * (Temp - Temp**3/3) + 0.5d0
        END IF

        mu = Viscosity + c3 * c1 

     CASE( 'user function' )
        str = ListGetString( Material, 'Viscosity Function' )
        Fnc = GetProcAddr( str, Quiet=.TRUE. )
        mu = MaterialUserFunction( Fnc, CurrentModel, Element, Nodes, n, nd, &
             Basis, dBasisdx, Viscosity, Velo, dVelodx )

     CASE DEFAULT 
        CALL WARN('EffectiveViscosity','Unknown material model')

     END SELECT


     ! Add a generic temperature coefficient at the integration point
     ! for backward compatibility this is activated by an existing keyword
     !--------------------------------------------------------------------
     c1 = ListGetConstReal( Material, 'Viscosity Temp Exp',GotIt)
     IF( GotIt ) THEN 	
        TempSol => VariableGet( CurrentModel % Variables, 'Temperature' )
        IF ( ASSOCIATED( TempSol) ) THEN
           TempPerm    => TempSol % Perm
           Temperature => TempSol % Values
        ELSE
           CALL Warn('EffectiveViscosity','variable Temperature needed for thermal viscosity model')
        END IF
        c2 = ListGetConstReal( Material, 'Viscosity Temp Ref')
        c3 = ListGetConstReal( Material, 'Viscosity Temp Offset',GotIt)
        Temp = SUM(Basis(1:n) * Temperature(TempPerm(Element % NodeIndexes(1:n))))
        TempCoeff = EXP(c1*(1/(Temp+c3)-1/c2))

        mu = TempCoeff * mu
     END IF


     !------------------------------------------------------------------------------
   END FUNCTION EffectiveViscosity
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Returns effective heat conductivity mainly related to turbulence models.
!------------------------------------------------------------------------------
   FUNCTION EffectiveConductivity( Conductivity,Density,Element, &
        Temperature,Ux,Uy,Uz,Nodes,n,nd,u,v,w ) RESULT(PCond)
!------------------------------------------------------------------------------
     USE ModelDescription

     REAL(KIND=dp)  :: Conductivity,Density,u,v,w,PCond, &
              Ux(:),Uy(:),Uz(:), Temperature(:), NodalViscosity(n)
     TYPE(Nodes_t)  :: Nodes
     INTEGER :: n,nd
     TYPE(Element_t),POINTER :: Element

#ifndef USE_ISO_C_BINDINGS
     INTERFACE
       FUNCTION MaterialUserFunction( Proc,Model,Element,Nodes,n,nd, &
          Basis,dBasisdx,Conductivity,Temp, dTempdx ) RESULT(s)
       USE Types
       INTEGER(KIND=AddrInt) :: Proc
       TYPE(Model_t) :: Model
       TYPE(Nodes_t) :: Nodes
       TYPE(Element_t), POINTER :: Element
       INTEGER :: n,nd
       REAL(KIND=dp) :: Basis(:),dBasisdx(:,:),Conductivity, &
                    Temp(:), dTempdx(:,:), s
       END FUNCTION MaterialUserFunction
     END INTERFACE
#endif
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3)
     LOGICAL :: stat,GotIt

     INTEGER :: i

     CHARACTER(LEN=MAX_NAME_LEN) :: ConductivityFlag
     TYPE(ValueList_t), POINTER :: Material
     REAL(KIND=dp) :: x, y, z, c1n(n), Temp(1), dTempdx(3,1), Pr_t, c_p,&
                      mu, Tmu, DetJ

     ! Temperature is needed for thermal models
     TYPE(Variable_t), POINTER :: TempSol 

     INTEGER(KIND=AddrInt) :: Fnc
     CHARACTER(LEN=MAX_NAME_LEN) :: str
!------------------------------------------------------------------------------
     PCond = Conductivity

     Material => GetMaterial( Element )
     ConductivityFlag = GetString( Material,'Heat Conductivity Model', GotIt)
     IF(.NOT. gotIt) RETURN

!------------------------------------------------------------------------------
     SELECT CASE( ConductivityFlag )
       CASE( 'ke','k-epsilon', 'turbulent' )
          stat = ElementInfo( Element,Nodes,u,v,w,detJ,Basis )

          c1n(1:n) = GetReal( Material, 'Heat Capacity' )
          c_p = SUM( Basis(1:n) * c1n(1:n) )

          c1n(1:n) = GetReal( Material, 'Viscosity' )
          mu = SUM( Basis(1:n) * c1n(1:n) )

          c1n(1:n) = GetReal( Material, 'Turbulent Prandtl Number',GotIt )
          IF ( GotIt ) THEN
            Pr_t = SUM( Basis(1:n) * c1n(1:n) )
          ELSE
            Pr_t = 0.85_dp
          END IF

          Tmu = EffectiveViscosity( mu,Density,Ux,Uy,Uz,Element, &
                        Nodes,n,nd,u,v,w ) - mu
          PCond = Conductivity + c_p * Tmu / Pr_t

       CASE( 'user function' )
         str = ListGetString( Material, 'Heat Conductivity Function' )
         Fnc = GetProcAddr( str, Quiet=.TRUE. )

         stat = ElementInfo( Element,Nodes,u,v,w, detJ, Basis,dBasisdx )
         Temp(1) = SUM( Basis(1:nd) * Temperature(1:nd) )
         DO i=1,3
            dTempdx(i,1) = SUM( dBasisdx(1:nd,i) * Temperature(1:nd) )
         END DO
       
         PCond = MaterialUserFunction( Fnc, CurrentModel, Element, Nodes, n, nd, &
              Basis, dBasisdx, Conductivity, Temp, dTempdx )

     CASE DEFAULT 
       CALL WARN('EffectiveConductivity','Unknown material model')

     END SELECT
!------------------------------------------------------------------------------

   END FUNCTION EffectiveConductivity
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Returns density stemming from various equations of state.
!------------------------------------------------------------------------------
   SUBROUTINE ElementDensity( Density, n )
!------------------------------------------------------------------------------
     INTEGER :: n
     REAL(KIND=dp) :: Density(:)
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: HeatCapacity(n), SpecificHeatRatio, GasConstant(n), &
      ReferencePressure, Pressure(n), Temperature(n), ReferenceTemperature(n), &
      HeatExpansionCoeff(n)
     LOGICAL :: Found
     TYPE(ValueList_t), POINTER :: Material

     Material => GetMaterial()

     SELECT CASE( GetString( Material, 'Compressibility Model', Found) )
     CASE( 'perfect gas', 'ideal gas' )
        HeatCapacity(1:n) = GetReal( Material, 'Heat Capacity' )
        SpecificHeatRatio = ListGetConstReal( Material, &
                 'Specific Heat Ratio', Found )
        IF ( .NOT.Found ) SpecificHeatRatio = 5.d0/3.d0

        GasConstant(1:n) = ( SpecificHeatRatio - 1.d0 ) *  &
            HeatCapacity(1:n) / SpecificHeatRatio

        ReferencePressure = GetCReal( Material, 'Reference Pressure', Found )
        IF ( .NOT.Found ) ReferencePressure = 0.0d0

        CALL GetScalarLocalSolution( Pressure, 'Pressure' )
        CALL GetScalarLocalSolution( Temperature, 'Temperature' )

        Density(1:n) = ( Pressure(1:n) + ReferencePressure ) / &
            ( GasConstant(1:n) * Temperature(1:n) )

     CASE( 'thermal' )
        HeatExpansionCoeff(1:n) = GetReal( Material, &
              'Heat Expansion Coefficient' )

        ReferenceTemperature(1:n) = GetReal( Material, &
              'Reference Temperature' )
        CALL GetScalarLocalSolution( Temperature, 'Temperature' )

        Density(1:n) = GetReal( Material,'Density' )
        Density(1:n) = Density(1:n) * ( 1 - HeatExpansionCoeff(1:n)  * &
           (  Temperature(1:n) - ReferenceTemperature(1:n) ) )

     CASE( 'user defined' )
        CALL GetScalarLocalSolution( Density, 'Density' )

     CASE DEFAULT
        Density(1:n) = GetReal( Material, 'Density' )
     END SELECT
   END SUBROUTINE ElementDensity

END MODULE MaterialModels

!> \} ElmerLib
