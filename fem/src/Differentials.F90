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
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 02 Jun 1997
! *
! *****************************************************************************/

!> \ingroup ElmerLib 
!> \{

!-----------------------------------------------------------------------------
!>  This module contains some built-in material laws, and also some 
!> vector utilities, curl, dot, cross, etc. Some of these may be 
!> of no use currently.
!-----------------------------------------------------------------------------
MODULE Differentials

  USE Types
  USE Lists
  USE LinearAlgebra
  USE ElementDescription

  IMPLICIT NONE

CONTAINS

!------------------------------------------------------------------------------
!> Computes the Lorentz force resulting from a magnetic field at 
!> integration point (u,v,w).
!------------------------------------------------------------------------------
  FUNCTION LorentzForce( Element,Nodes,u,v,w,n ) RESULT(L)
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t) :: Nodes
    INTEGER :: n
    REAL(KIND=dp) :: L(3),u,v,w,x,y,z
!------------------------------------------------------------------------------
    TYPE(Variable_t), POINTER :: Mx,My,Mz,MFx,MFy,MFz
    INTEGER :: i,j,k,bfId
    LOGICAL :: stat,GotIt
    INTEGER, POINTER :: NodeIndexes(:)

    TYPE(ValueList_t), POINTER :: Material

    REAL(KIND=dp) :: B(3),dHdx(3,3)
    REAL(KIND=dp) :: dBasisdx(n,3),SqrtElementMetric
    REAL(KIND=dp) :: Basis(n),Permeability(n),mu

    REAL(KIND=dp) :: ExtMx(n),ExtMy(n),ExtMz(n)

    REAL(KIND=dp) :: SqrtMetric,Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3)
!------------------------------------------------------------------------------
    L = 0.0D0

!     bfId = ListGetInteger( CurrentModel % Bodies( Element % BodyId ) % Values, &
!                 'Body Force', GotIt, 1, CurrentModel % NumberOFBodyForces )

!     IF ( .NOT.GotIt ) RETURN

!     IF ( .NOT.ListGetLogical( CurrentModel % BodyForces( &
!           bfId ) % Values, 'Lorentz Force' , GotIt ) ) RETURN
!------------------------------------------------------------------------------
    NodeIndexes => Element % NodeIndexes

    Mx => VariableGet( CurrentModel % Variables, 'Magnetic Field 1' )
    My => VariableGet( CurrentModel % Variables, 'Magnetic Field 2' )
    Mz => VariableGet( CurrentModel % Variables, 'Magnetic Field 3' )
    IF ( .NOT.ASSOCIATED( Mx ) ) RETURN

    IF ( ANY(Mx % Perm(NodeIndexes)<=0) ) RETURN

    k = ListGetInteger( CurrentModel % Bodies &
                    (Element % BodyId) % Values, 'Material', &
                     minv=1, maxv=CurrentModel % NumberOFMaterials )
    Material => CurrentModel % Materials(k) % Values

    Permeability(1:n) = ListGetReal( Material, 'Magnetic Permeability', &
                              n, NodeIndexes ) 
!------------------------------------------------------------------------------
    ExtMx(1:n) = ListGetReal( Material, 'Applied Magnetic Field 1', &
                    n,NodeIndexes, Gotit )

    ExtMy(1:n) = ListGetReal( Material, 'Applied Magnetic Field 2', &
                  n,NodeIndexes, Gotit )

    ExtMz(1:n) = ListGetReal( Material, 'Applied Magnetic Field 3', &
                  n,NodeIndexes, Gotit )

! If you want to use time-domain solution for high-frequendy part,
! leave external field out. Better to use frequency-domain solver!
#if 1
    MFx => VariableGet( CurrentModel % Variables, 'Magnetic Flux Density 1' )
    MFy => VariableGet( CurrentModel % Variables, 'Magnetic Flux Density 2' )
    MFz => VariableGet( CurrentModel % Variables, 'Magnetic Flux Density 3' )
    IF ( ASSOCIATED( MFx ) ) THEN
      ExtMx(1:n) = ExtMx(1:n) + MFx % Values(MFx % Perm(NodeIndexes))
      ExtMy(1:n) = ExtMy(1:n) + MFy % Values(MFy % Perm(NodeIndexes))
      ExtMz(1:n) = ExtMz(1:n) + MFz % Values(MFz % Perm(NodeIndexes))
    END IF
#endif

!------------------------------------------------------------------------------
!   Get element info 
!------------------------------------------------------------------------------
    stat = ElementInfo( Element,Nodes,u,v,w,SqrtElementMetric, &
               Basis,dBasisdx )
!------------------------------------------------------------------------------
    B(1) = SUM( Basis(1:n)*Mx % Values(Mx % Perm(NodeIndexes)) )
    B(2) = SUM( Basis(1:n)*My % Values(My % Perm(NodeIndexes)) )
    B(3) = SUM( Basis(1:n)*Mz % Values(Mz % Perm(NodeIndexes)) )

    B(1) = B(1) + SUM( Basis(1:n)*ExtMx(1:n) )
    B(2) = B(2) + SUM( Basis(1:n)*ExtMy(1:n) )
    B(3) = B(3) + SUM( Basis(1:n)*ExtMz(1:n) )

    DO i=1,3
      dHdx(1,i) = SUM( dBasisdx(1:n,i)* &
           Mx % Values(Mx % Perm(NodeIndexes)) / Permeability(1:n) )
      dHdx(2,i) = SUM( dBasisdx(1:n,i)* &
           My % Values(My % Perm(NodeIndexes)) / Permeability(1:n) )
      dHdx(3,i) = SUM( dBasisdx(1:n,i)* &
           Mz % Values(Mz % Perm(NodeIndexes)) / Permeability(1:n) )
    END DO
!------------------------------------------------------------------------------
!       Get coordinate system info
!------------------------------------------------------------------------------
    x = SUM( Nodes % x(1:n) * Basis(1:n) )
    y = SUM( Nodes % y(1:n) * Basis(1:n) )
    z = SUM( Nodes % z(1:n) * Basis(1:n) )
    CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,x,y,z )
    IF ( CurrentCoordinateSystem() /= Cartesian ) CALL InvertMatrix( Metric,3 )

    mu = SUM( Permeability(1:n)*Basis(1:n) )
    L = ComputeLorentz( B,dHdx,mu,SqrtMetric,Metric,Symb )

CONTAINS

!------------------------------------------------------------------------------
  FUNCTION ComputeLorentz( B,dHdx,mu,SqrtMetric,Metric,Symb ) RESULT(LF)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: B(:),dHdx(:,:),mu,LF(3),SqrtMetric,Metric(:,:),Symb(:,:,:)
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,l,m
    REAL(KIND=dp) :: Bc(3),Ji(3),Jc(3),s,Perm(3,3,3),r
!------------------------------------------------------------------------------

    IF ( CurrentCoordinateSystem() == Cartesian ) THEN
      Ji(1) = dHdx(3,2) - dHdx(2,3)
      Ji(2) = dHdx(1,3) - dHdx(3,1)
      Ji(3) = dHdx(2,1) - dHdx(1,2)
      LF(1) = Ji(2)*B(3) - Ji(3)*B(2)
      LF(2) = Ji(3)*B(1) - Ji(1)*B(3)
      LF(3) = Ji(1)*B(2) - Ji(2)*B(1)
      RETURN
    END IF

    r = SqrtMetric

    IF ( CurrentCoordinateSystem()  == CylindricSymmetric ) THEN
      Ji(1) = -dHdx(3,2)
      Ji(2) =  dHdx(3,1)
      IF (r > 1.0d-10) THEN
         Ji(2) = Ji(2) + B(3)/(r*mu)
      ELSE
         Ji(2) = Ji(2) + Ji(2)
      END IF
      Ji(3) = dHdx(1,2) - dHdx(2,1)

      LF(1) = Ji(3)*B(2) - Ji(2)*B(3)
      LF(2) = Ji(1)*B(3) - Ji(3)*B(1)
! You might want to use SI units for the azimuthal component,
! if you compute Lorentz force at nodal points and symmetry axis,
! otherwise you divide by zero.
#ifdef SI_UNITS
      LF(3) = Ji(2)*B(1) - Ji(1)*B(2)
#else
      IF (r > 1.0d-10) THEN
         LF(3) = ( Ji(2)*B(1) - Ji(1)*B(2) ) / r
      ELSE
         LF(3) = 0.d0
      END IF
#endif
      RETURN
    END IF

    Perm = 0
    Perm(1,2,3) = -1.0d0 / SqrtMetric
    Perm(1,3,2) =  1.0d0 / SqrtMetric
    Perm(2,1,3) =  1.0d0 / SqrtMetric
    Perm(2,3,1) = -1.0d0 / SqrtMetric
    Perm(3,1,2) = -1.0d0 / SqrtMetric
    Perm(3,2,1) =  1.0d0 / SqrtMetric
!------------------------------------------------------------------------------

    Bc = 0.0d0
    DO i=1,3
      DO j=1,3
        Bc(i) = Bc(i) + Metric(i,j)*B(j)
      END DO
    END DO

!------------------------------------------------------------------------------

    Ji = 0.0d0
    DO i=1,3
      s = 0.0D0
      DO j=1,3
        DO k=1,3
          IF ( Perm(i,j,k) /= 0 ) THEN
            DO l=1,3
              s = s + Perm(i,j,k)*Metric(j,l)*dHdx(l,k)
              DO m=1,3
                s = s + Perm(i,j,k)*Metric(j,l)*Symb(k,m,l)*B(m)/mu
              END DO
            END DO
          END IF
        END DO
      END DO
      Ji(i) = s
    END DO
 
    Jc = 0.0d0
    DO i=1,3
      DO j=1,3
        Jc(i) = Jc(i) + Metric(i,j)*Ji(j)
      END DO
    END DO
!------------------------------------------------------------------------------

    LF = 0.0d0
    DO i=1,3
      s = 0.0D0
      DO j=1,3
        DO k=1,3
          IF ( Perm(i,j,k) /= 0 ) THEN
            s = s + Perm(i,j,k)*Jc(k)*Bc(j)
          END IF
        END DO
      END DO
      LF(i) = s
    END DO
!------------------------------------------------------------------------------
  END FUNCTION ComputeLorentz
!------------------------------------------------------------------------------
  END FUNCTION LorentzForce
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Compute the Joule heating at integration point (u,v,w) given the 
!> appropriate electrostatic or magnetic field that initiates the 
!> current through a conductor. 
!------------------------------------------------------------------------------
  FUNCTION JouleHeat( Element,Nodes,u,v,w,n ) RESULT(JouleH)
!------------------------------------------------------------------------------
    TYPE(Element_t) :: Element
    TYPE(Nodes_t) :: Nodes
    INTEGER :: n
    REAL(KIND=dp) :: JouleH,u,v,w,x,y,z
!------------------------------------------------------------------------------
    TYPE(Variable_t), POINTER :: Mx,My,Mz,MFx,MFy,MFz
    INTEGER :: i,j,k,bfId
    LOGICAL :: stat,GotIt
    INTEGER, POINTER :: NodeIndexes(:)

    TYPE(ValueList_t), POINTER :: Material

    REAL(KIND=dp) :: B(3),dHdx(3,3)
    REAL(KIND=dp) :: dBasisdx(n,3),SqrtElementMetric
    REAL(KIND=dp) :: Basis(n),Permeability(n), &
             ElectricConductivity(n)

    REAL(KIND=dp) :: ExtMx(n),ExtMy(n),ExtMz(n)
    REAL(KIND=dp) :: mu,elcond,SqrtMetric,Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3)

    INTEGER, SAVE :: JouleMode = 0 
    TYPE(Variable_t), POINTER, SAVE :: Jvar
!------------------------------------------------------------------------------
    JouleH = 0.0_dp

    bfId = ListGetInteger( CurrentModel % Bodies( Element % BodyId ) % &
         Values, 'Body Force', GotIt, 1, CurrentModel % NumberOfBodyForces )

    IF ( .NOT.GotIt ) RETURN
    
    IF ( .NOT.ListGetLogical( CurrentModel % BodyForces( &
         bfId ) % Values, 'Joule Heat' , GotIt ) ) RETURN
!------------------------------------------------------------------------------

    IF( JouleMode == 0 ) THEN
      Jvar => VariableGet( CurrentModel % Variables, 'Joule Heating e' )
      IF ( ASSOCIATED( Jvar ) ) JouleMode = 1 

      IF( JouleMode == 0 ) THEN
        Jvar => VariableGet( CurrentModel % Variables, 'Joule Field' )
        IF ( ASSOCIATED( Jvar ) ) JouleMode = 2
      END IF

      IF( JouleMode == 0 ) THEN
        Jvar => VariableGet( CurrentModel % Variables, 'Potential' )
        IF ( ASSOCIATED( Jvar ) ) JouleMode = 3
      END IF

      IF( JouleMode == 0 ) THEN
        Jvar => VariableGet( CurrentModel % Variables, 'Magnetic Field 1' )
        IF ( ASSOCIATED( Jvar ) ) JouleMode = 4 
      END IF

      IF( JouleMode == 0 ) THEN
        CALL Warn('JouleHeat','Joule heating requested but no field to compute it!')
        RETURN
      END IF
    END IF

    IF( JouleMode == 1 ) THEN
      NodeIndexes => Element % DgIndexes 
    ELSE
      NodeIndexes => Element % NodeIndexes
    END IF
    IF( ANY( Jvar % Perm( NodeIndexes ) == 0 ) ) RETURN

    
    !------------------------------------------------------------------------------
    ! The simplest model just evaluates precomputed elemental heating at integration point
    !------------------------------------------------------------------------------
    IF( JouleMode == 1 ) THEN
      stat = ElementInfo( Element,Nodes,u,v,w,SqrtElementMetric,Basis )
      JouleH = SUM( Basis(1:n) * Jvar % Values( Jvar % Perm( NodeIndexes ) ) )

      ! Make an early exit since we don't need conductivity 
      RETURN
    END IF


    !------------------------------------------------------------------------------
    !  All other models require electric conductivity
    !------------------------------------------------------------------------------
    k = ListGetInteger( CurrentModel % Bodies &
         (Element % BodyId) % Values, 'Material')
    Material => CurrentModel % Materials(k) % Values
    
    ElectricConductivity(1:n) = ListGetReal( Material, &
        'Electric Conductivity',n,NodeIndexes,GotIt )

    IF( JouleMode == 2 ) THEN
      ! This model uses a "joule field" and multiplies it with conductivity at the integration point
      stat = ElementInfo( Element,Nodes,u,v,w,SqrtElementMetric,Basis )     
      elcond = SUM( ElectricConductivity(1:n) * Basis(1:n) )
      JouleH = elcond * SUM( Basis(1:n) * Jvar % Values(Jvar % Perm(NodeIndexes)) )

    ELSE IF( JouleMode == 3 ) THEN
      ! This model uses potential and evaluates the current from its gradient
      stat = ElementInfo( Element,Nodes,u,v,w,SqrtElementMetric,Basis,dBasisdx )
      elcond = SUM( ElectricConductivity(1:n) * Basis(1:n) )

      B(1) = SUM( dBasisdx(1:n,1) * Jvar % Values(Jvar % Perm(NodeIndexes)) )
      B(2) = SUM( dBasisdx(1:n,2) * Jvar % Values(Jvar % Perm(NodeIndexes)) )
      B(3) = SUM( dBasisdx(1:n,3) * Jvar % Values(Jvar % Perm(NodeIndexes)) )     
      JouleH = elcond * SUM( B * B )

    ELSE IF( JouleMode == 4 ) THEN
      !------------------------------------------------------------------------------
      !  Magnetic induction equation, might be obsolete
      !------------------------------------------------------------------------------
      stat = ElementInfo( Element,Nodes,u,v,w,SqrtElementMetric,Basis,dBasisdx )
      elcond = SUM( ElectricConductivity(1:n) * Basis(1:n) )
      IF( elcond < TINY( elcond ) ) RETURN

      Permeability(1:n) = ListGetReal( Material, 'Magnetic Permeability', &
          n, NodeIndexes ) 
      
      Mx => VariableGet( CurrentModel % Variables, 'Magnetic Field 1' )
      My => VariableGet( CurrentModel % Variables, 'Magnetic Field 2' )
      Mz => VariableGet( CurrentModel % Variables, 'Magnetic Field 3' )
      
      !------------------------------------------------------------------------------
      ExtMx(1:n) = ListGetReal( Material, 'Applied Magnetic Field 1', &
          n,NodeIndexes, Gotit )       
      ExtMy(1:n) = ListGetReal( Material, 'Applied Magnetic Field 2', &
          n,NodeIndexes, Gotit )       
      ExtMz(1:n) = ListGetReal( Material, 'Applied Magnetic Field 3', &
          n,NodeIndexes, Gotit )
      !
      MFx => VariableGet( CurrentModel % Variables, 'Magnetic Flux Density 1' )
      MFy => VariableGet( CurrentModel % Variables, 'Magnetic Flux Density 2' )
      MFz => VariableGet( CurrentModel % Variables, 'Magnetic Flux Density 3' )
      IF ( ASSOCIATED( MFx ) ) THEN
        ExtMx(1:n) = ExtMx(1:n) + MFx % Values(MFx % Perm(NodeIndexes))
        ExtMy(1:n) = ExtMy(1:n) + MFy % Values(MFy % Perm(NodeIndexes))
        ExtMz(1:n) = ExtMz(1:n) + MFz % Values(MFz % Perm(NodeIndexes))
      END IF
      
      !------------------------------------------------------------------------------
      B(1) = SUM( Basis(1:n) * Mx % Values(Mx % Perm(NodeIndexes)) )
      B(2) = SUM( Basis(1:n) * My % Values(My % Perm(NodeIndexes)) )
      B(3) = SUM( Basis(1:n) * Mz % Values(Mz % Perm(NodeIndexes)) )
      
      B(1) = B(1) + SUM( Basis(1:n) * ExtMx(1:n) )
      B(2) = B(2) + SUM( Basis(1:n) * ExtMy(1:n) )
      B(3) = B(3) + SUM( Basis(1:n) * ExtMz(1:n) )
      
      mu = SUM( Basis(1:n) * Permeability(1:n) )
      DO i=1,3
        dHdx(1,i) = SUM( dBasisdx(1:n,i)* &
            Mx % Values(Mx % Perm(NodeIndexes))/Permeability(1:n) )
        dHdx(2,i) = SUM( dBasisdx(1:n,i)* &
            My % Values(My % Perm(NodeIndexes))/Permeability(1:n) )
        dHdx(3,i) = SUM( dBasisdx(1:n,i)* &
            Mz % Values(Mz % Perm(NodeIndexes))/Permeability(1:n) )
      END DO
      
      IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
        x = SUM( Nodes % x(1:n) * Basis(1:n))
        y = SUM( Nodes % y(1:n) * Basis(1:n))
        z = SUM( Nodes % z(1:n) * Basis(1:n))
        CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,x,y,z )
        CALL InvertMatrix( Metric,3 )
      END IF
      JouleH = ComputeMagneticHeat( B,dHdx,mu,SqrtMetric,Metric,Symb ) / &
          elcond            
    END IF
    
CONTAINS

!------------------------------------------------------------------------------
  FUNCTION ComputeMagneticHeat( B,dHdx,mu,SqrtMetric,Metric,Symb ) RESULT(JH)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: B(:),dHdx(:,:),mu,JH,SqrtMetric,Metric(:,:),Symb(:,:,:)
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,l,m
    REAL(KIND=dp) :: Bc(3),Ji(3),Jc(3),s,Perm(3,3,3),r
!------------------------------------------------------------------------------

    IF ( CurrentCoordinateSystem() == Cartesian ) THEN
      Ji(1) = dHdx(3,2) - dHdx(2,3)
      Ji(2) = dHdx(1,3) - dHdx(3,1)
      Ji(3) = dHdx(2,1) - dHdx(1,2)
      JH = Ji(1)*Ji(1) + Ji(2)*Ji(2) + Ji(3)*Ji(3)
      RETURN
    END IF

    IF ( CurrentCoordinateSystem() == CylindricSymmetric ) THEN
      r = SqrtMetric
      Ji(1) = -dHdx(3,2)
      Ji(2) = B(3)/(r*mu) + dHdx(3,1)
      Ji(3) = dHdx(1,2) - dHdx(2,1)
      JH = Ji(1)*Ji(1) + Ji(2)*Ji(2) + Ji(3)*Ji(3)
      RETURN
    END IF

    Perm = 0
    Perm(1,2,3) = -1.0d0 / SqrtMetric
    Perm(1,3,2) =  1.0d0 / SqrtMetric
    Perm(2,1,3) =  1.0d0 / SqrtMetric
    Perm(2,3,1) = -1.0d0 / SqrtMetric
    Perm(3,1,2) = -1.0d0 / SqrtMetric
    Perm(3,2,1) =  1.0d0 / SqrtMetric
!------------------------------------------------------------------------------

    Bc = 0.0d0
    DO i=1,3
      DO j=1,3
        Bc(i) = Bc(i) + Metric(i,j)*B(j)
      END DO
    END DO

!------------------------------------------------------------------------------

    Ji = 0.0d0
    DO i=1,3
      s = 0.0D0
      DO j=1,3
        DO k=1,3
          IF ( Perm(i,j,k) /= 0 ) THEN
            DO l=1,3
              s = s + Perm(i,j,k)*Metric(j,l)*dHdx(l,k)
              DO m=1,3
                s = s + Perm(i,j,k)*Metric(j,l)*Symb(k,m,l)*B(m)/mu
              END DO
            END DO
          END IF
        END DO
      END DO
      Ji(i) = s
    END DO
 
    Jc = 0.0d0
    DO i=1,3
      DO j=1,3
        Jc(i) = Jc(i) + Metric(i,j)*Ji(j)
      END DO
    END DO

!------------------------------------------------------------------------------

    JH = 0.0d0
    DO i=1,3
      JH = JH + Ji(i) * Jc(i)
    END DO
!------------------------------------------------------------------------------
  END FUNCTION ComputeMagneticHeat
!------------------------------------------------------------------------------
  END FUNCTION JouleHeat
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Compute the curl vector B = curl(A) at all nodes.
!> \deprecated Is this used anywhere?
!------------------------------------------------------------------------------
  SUBROUTINE Curl( Ax,Ay,Az,Bx,By,Bz,Reorder )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Ax(:),Ay(:),Az(:),Bx(:),By(:),Bz(:)
    INTEGER :: Reorder(:)
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t) :: Nodes 

    LOGICAL :: Stat
    INTEGER :: i,j,k,l,m,n,p,q,t
    REAL(KIND=dp) :: x,y,z,u,v,w,s,A(3),B(3),dx(3,3)

    INTEGER :: Perm(3,3,3)
    INTEGER, POINTER :: NodeIndexes(:),Visited(:)

    REAL(KIND=dp), ALLOCATABLE :: dBasisdx(:,:)
    REAL(KIND=dp), ALLOCATABLE :: Basis(:),aaz(:)

    REAL(KIND=dp) :: SqrtElementMetric
    REAL(KIND=dp) :: SqrtMetric,Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3)
!------------------------------------------------------------------------------

    ALLOCATE( Visited(CurrentModel % NumberOfNodes) )
    Visited = 0

    Perm = 0
    Perm(1,2,3) = -1
    Perm(1,3,2) =  1
    Perm(2,1,3) =  1
    Perm(2,3,1) = -1
    Perm(3,1,2) = -1
    Perm(3,2,1) =  1

    n = CurrentModel % Mesh % MaxElementNodes
    ALLOCATE( dBasisdx(n,3), Basis(n), aaz(n) )
    ALLOCATE( Nodes % x(n), Nodes % y(n), Nodes % z(n) )

    Bx = 0.0D0
    By = 0.0D0
    Bz = 0.0D0
!------------------------------------------------------------------------------
!   Go through model elements, we will compute on average of elementwise
!   curls on nodes of the model
!------------------------------------------------------------------------------
    DO t=1,CurrentModel % NumberOfBulkElements
!------------------------------------------------------------------------------
      Element => CurrentModel % Elements(t)
      n = Element % TYPE % NumberOfNodes
      NodeIndexes => Element % NodeIndexes

      Nodes % x(1:n) = CurrentModel % Nodes % x( NodeIndexes )
      Nodes % y(1:n) = CurrentModel % Nodes % y( NodeIndexes )
      Nodes % z(1:n) = CurrentModel % Nodes % z( NodeIndexes )
!------------------------------------------------------------------------------
!     Through element nodes
!------------------------------------------------------------------------------

      IF (MINVAL(Reorder(NodeIndexes)) > 0) THEN

      DO p=1,n
        q = Reorder(NodeIndexes(p))
        u = Element % TYPE % NodeU(p)
        v = Element % TYPE % NodeV(p)

        IF ( Element % TYPE % DIMENSION == 3 ) THEN
          w = Element % TYPE % NodeW(p)
        ELSE
          w = 0.0D0
        END IF
!------------------------------------------------------------------------------
!       Get element basis functions, basis function derivatives, etc,
!       and compute partials derivatives of the vector A with respect
!       to global coordinates.
!------------------------------------------------------------------------------
        stat = ElementInfo( Element,Nodes,u,v,w,SqrtElementMetric, &
                  Basis,dBasisdx )
!------------------------------------------------------------------------------
!       Get coordinate system info
!------------------------------------------------------------------------------
        IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
          x = SUM( Nodes % x(1:n) * Basis(1:n))
          y = SUM( Nodes % y(1:n) * Basis(1:n))
          z = SUM( Nodes % z(1:n) * Basis(1:n))
          CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,x,y,z )
          CALL InvertMatrix( Metric,3 )
        END IF

!        print *, '***',p
!        print *, '***',x,y,z
!        print *, '***',NodeIndexes
!        print *, '***',Reorder(NodeIndexes)

        DO k=1,3
          dx(1,k) = SUM( dBasisdx(1:n,k) * Ax(Reorder(NodeIndexes)) )
          dx(2,k) = SUM( dBasisdx(1:n,k) * Ay(Reorder(NodeIndexes)) )
          dx(3,k) = SUM( dBasisdx(1:n,k) * Az(Reorder(NodeIndexes)) )
        END DO
!------------------------------------------------------------------------------
!       And compute the curl for the node of the current element
!------------------------------------------------------------------------------
        IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
          A(1) = Ax(q)
          A(2) = Ay(q)
          A(3) = Az(q)
          B = 0.0D0
          IF ( ABS(SqrtMetric) > 1.0d-15 ) THEN
            DO i=1,3
              s = 0.0D0
              DO j=1,3
                DO k=1,3
                  IF ( Perm(i,j,k) /= 0 ) THEN
                    DO l=1,3
                      s = s + Perm(i,j,k)*Metric(j,l)*dx(l,k)
                      DO m=1,3
                        s = s + Perm(i,j,k)*Metric(j,l)*Symb(k,m,l)*A(m)
                      END DO
                    END DO
                  END IF
                END DO
              END DO
              B(i) = s
            END DO
       
            Bx(q) = Bx(q) + B(1) / SqrtMetric
            By(q) = By(q) + B(2) / SqrtMetric
            Bz(q) = Bz(q) + B(3) / SqrtMetric
          END IF
        ELSE

          Bx(q) = Bx(q) + dx(3,2) - dx(2,3)
          By(q) = By(q) + dx(1,3) - dx(3,1)
          Bz(q) = Bz(q) + dx(2,1) - dx(1,2)

        END IF
        Visited(q) = Visited(q) + 1
      END DO

    END IF

!------------------------------------------------------------------------------
    END DO
!------------------------------------------------------------------------------
!   Finally, compute average of the the curls at nodes
!------------------------------------------------------------------------------
    DO i=1,CurrentModel % NumberOfNodes
      IF ( Visited(i) > 0 ) THEN
        Bx(i) = Bx(i) / Visited(i)
        By(i) = By(i) / Visited(i)
        Bz(i) = Bz(i) / Visited(i)
      END IF
    END DO

    DEALLOCATE( Visited, Nodes % x, Nodes % y, Nodes % z, Basis, dBasisdx, aaz )
!------------------------------------------------------------------------------
  END SUBROUTINE Curl
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Compute the curl in axisymmetric coordinates.
!> \deprecated Is this used anywhere?
!------------------------------------------------------------------------------
SUBROUTINE AxiSCurl( Ar,Az,Ap,Br,Bz,Bp,Reorder )
!------------------------------------------------------------------------------
  IMPLICIT NONE
  REAL(KIND=dp) :: Ar(:),Az(:),Ap(:),Br(:),Bz(:),Bp(:)
  INTEGER :: Reorder(:)

  TYPE(Element_t), POINTER :: Element
  TYPE(Nodes_t) :: Nodes 

  LOGICAL :: Stat

  INTEGER, POINTER :: NodeIndexes(:),Visited(:)
  INTEGER :: p,q,i,t,n

  REAL(KIND=dp) :: u,v,w,r

  REAL(KIND=dp) :: SqrtElementMetric
  REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:)
!------------------------------------------
  ALLOCATE( Visited(CurrentModel % NumberOfNodes) )

  n = CurrentModel % Mesh % MaxElementDOFs
  ALLOCATE( Basis(n), dBasisdx(n,3) )
  ALLOCATE(Nodes % x(n),Nodes % y(n),Nodes % z(n))

  Visited = 0

  Br = 0.0d0
  Bz = 0.0d0
  Bp = 0.0d0

  DO t=1,CurrentModel % NumberOfBulkElements

     Element => CurrentModel % Elements(t)
     n = Element % TYPE % NumberOfNodes
     NodeIndexes => Element % NodeIndexes

     Nodes % x(1:n) = CurrentModel % Nodes % x( NodeIndexes )
     Nodes % y(1:n) = CurrentModel % Nodes % y( NodeIndexes )
     Nodes % z(1:n) = CurrentModel % Nodes % z( NodeIndexes )

     IF ( MINVAL(Reorder(NodeIndexes)) > 0 ) THEN
        DO p=1,n
           q = Reorder(NodeIndexes(p))
           u = Element % TYPE % NodeU(p)
           v = Element % TYPE % NodeV(p)

           IF ( Element % TYPE % DIMENSION == 3 ) THEN
              w = Element % TYPE % NodeW(p)
           ELSE
              w = 0.0D0
           END IF

           stat = ElementInfo( Element, Nodes, u, v, w, SqrtElementMetric, &
                          Basis, dBasisdx )

           r = SUM( Basis(1:n) * Nodes % x(1:n) )

           Br(q) = Br(q) - SUM( dBasisdx(1:n,2)*Ap(Reorder(NodeIndexes)) )

           Bp(q) = Bp(q) + SUM( dBasisdx(1:n,2) * Ar(Reorder(NodeIndexes)) ) &
                - SUM( dBasisdx(1:n,1) * Az(Reorder(NodeIndexes)) )

           Bz(q) = Bz(q) + SUM( dBasisdx(1:n,1) * Ap(Reorder(NodeIndexes)) )

           IF (r > 1.0d-10) THEN
              Bz(q) = Bz(q) + SUM( Basis(1:n)*Ap(Reorder(NodeIndexes)) ) / r
           ELSE
              Bz(q) = Bz(q) + SUM( dBasisdx(1:n,1)*Ap(Reorder(NodeIndexes)) )
           END IF

           Visited(q) = Visited(q) + 1           
        END DO
     END IF
  END DO

  DO i=1,CurrentModel % NumberOfNodes
     IF ( Visited(i) > 0 ) THEN
        Br(i) = Br(i) / Visited(i)
        Bp(i) = Bp(i) / Visited(i)
        Bz(i) = Bz(i) / Visited(i)
     END IF
  END DO

  DEALLOCATE( Visited, Basis, dBasisdx )
  DEALLOCATE( Nodes % x, Nodes % y, Nodes % z )

!------------------------------------------------------------------------------
END SUBROUTINE AxiSCurl
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>  Compute cross product of given vectors: C = A x B in generalized coordinates.
!> \deprecated Is this used anywhere?
!------------------------------------------------------------------------------
  SUBROUTINE Cross( Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,n )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,n
    REAL(KIND=dp) :: SqrtMetric,Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3),x,y,z
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!   Compute the cross product
!------------------------------------------------------------------------------
    Cx = Ay * Bz - Az * By
    Cy = Az * Bx - Ax * Bz
    Cz = Ax * By - Ay * Bx
!------------------------------------------------------------------------------
!   Make contravariant
!------------------------------------------------------------------------------
    IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
      x = CurrentModel % Nodes % x(n)
      y = CurrentModel % Nodes % y(n)
      z = CurrentModel % Nodes % z(n)
      CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,x,y,z )

      x = SqrtMetric * Cx
      y = SqrtMetric * Cy
      z = SqrtMetric * Cz

      Cx = Metric(1,1)*x + Metric(1,2)*y + Metric(1,3)*z
      Cy = Metric(2,1)*x + Metric(2,2)*y + Metric(2,3)*z
      Cz = Metric(3,1)*x + Metric(3,2)*y + Metric(3,3)*z
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE Cross
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>  Compute dot product of given vectors: L=A \cdot B in orthogonal 
!> coordinate system.
!> \deprecated Is this used anywhere?
!------------------------------------------------------------------------------
  FUNCTION Dot( Ax,Ay,Az,Bx,By,Bz,n ) RESULT(L)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Ax,Ay,Az,Bx,By,Bz
    REAL(KIND=dp) :: L
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,n
    REAL(KIND=dp) :: SqrtMetric,Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3),x,y,z
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!     Compute the dot product
!------------------------------------------------------------------------------
    IF ( CurrentCoordinateSystem() == Cartesian ) THEN
      L =  Ax*Bx + Ay*By + Az*Bz
    ELSE
      x = CurrentModel % Nodes % x(n)
      y = CurrentModel % Nodes % y(n)
      z = CurrentModel % Nodes % z(n)
      CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,x,y,z )
!
!  NOTE: this is for orthogonal coordinates only
!
      L = Ax*Bx / Metric(1,1) + Ay*By / Metric(2,2) + Az*Bz / Metric(3,3)
    END IF
!------------------------------------------------------------------------------
  END FUNCTION Dot
!------------------------------------------------------------------------------

END MODULE Differentials

!> \} ElmerLib

