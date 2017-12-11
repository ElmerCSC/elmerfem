!#########################################################
!# User functions for the bodies in rotation test case
!#  Ux: return the Ux velocity (Ux=2Pi(y-0.5))
!#  Uy: return the Uy velocity (Uy=2Pi(0.5-x))
!#  HIni: Return the initial shape of the bodies
!#########################################################
      FUNCTION Ux(Model,nodenumber,VarIn) RESULT(VarOut)
      USE DefUtils
       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: VarIn
       REAL(kind=dp) :: VarOut
       !-----------------
       REAL(kind=dp) :: y

       y=Model%Mesh%Nodes%y(nodenumber)

       VarOut=2.*Pi*(y-0.5)
      End FUNCTION Ux

      FUNCTION Uy(Model,nodenumber,VarIn) RESULT(VarOut)
      USE DefUtils
       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: VarIn
       REAL(kind=dp) :: VarOut
       !-----------------
       REAL(kind=dp) :: x

       x=Model%Mesh%Nodes%x(nodenumber)
       VarOut=2.*Pi*(0.5-x)
      End FUNCTION Uy

      FUNCTION HIni(Model,nodenumber,VarIn) RESULT(VarOut)
      USE DefUtils
       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: VarIn(2) !x,y
       REAL(kind=dp) :: VarOut
       !-----------------
       TYPE(ValueList_t), POINTER :: Constants
       REAL(kind=dp) :: x,y
       Real(kind=dp) :: Cone,Bump,Cylinder
       LOGICAL,SAVE :: HasCone,HasBump,HasCylinder
       LOGICAL :: Found
       LOGICAL,SAVE :: FirstTime=.True.

        IF (FirstTime) THEN
         Constants => Model % Constants
         IF (.NOT.ASSOCIATED(Constants)) THEN
           HasCone=.True.
           HasBump=.True.
           HasCylinder=.True.
         ELSE
           HasCone=ListGetLogical(Constants, 'WithCone',Found)
           HasBump =ListGetLogical(Constants, 'WithBump',Found)
           HasCylinder=ListGetLogical(Constants, 'WithCylinder',Found)
         END IF
          PRINT *,"HasCone",HasCone
          PRINT *,"HasBump",HasBump
          PRINT *,"HasCylinder",HasCylinder
          FirstTime=.FALSE.
        END IF

        x=VarIn(1)
        y=VarIn(2)

        VarOut=0._dp
        IF (HasCone) VarOut=VarOut+Cone(x,y)
        IF (HasBump) VarOut=VarOut+Bump(x,y)
        IF (HasCylinder) VarOut=VarOut+Cylinder(x,y)
      End FUNCTION HIni

      FUNCTION Cone(x,y) RESULT(VarOut)
      USE TYPES
       implicit none
       !-----------------
       REAL(kind=dp) :: VarOut
       !-----------------
       REAL(kind=dp) :: x,y,r
       REAL(kind=dp),parameter :: x0=0.5,y0=0.25,r0=0.15

        r=sqrt((x-x0)**2.0+(y-y0)**2.0)
        VarOut=max(1.0-r/r0,0.0)

      End FUNCTION Cone

      FUNCTION Bump(x,y) RESULT(VarOut)
      USE TYPES
       implicit none
       !-----------------
       REAL(kind=dp) :: VarOut
       !-----------------
       REAL(kind=dp) :: x,y,r
       REAL(kind=dp),parameter :: x0=0.25,y0=0.5,r0=0.15

        r=sqrt((x-x0)**2.0+(y-y0)**2.0)/r0

        VarOut=0.5*(1.0+cos(Pi*min(r,1.0_dp)))


      END FUNCTION Bump

      FUNCTION Cylinder(x,y) RESULT(VarOut)
      USE TYPES
       implicit none
       !-----------------
       REAL(kind=dp) :: VarOut
       !-----------------
       REAL(kind=dp) :: x,y,r
       REAL(kind=dp),parameter :: x0=0.5,y0=0.75,r0=0.15

       r=sqrt((x-x0)**2.0+(y-y0)**2.0)/r0

       if ((r.LT.1.0_dp).AND.((abs(x-x0).GT.0.0225).OR.(y.GT.0.85))) THEN
          VarOut=1.0_dp
       else
          VarOut=0.0_dp
       endif
      END FUNCTION Cylinder


