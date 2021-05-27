# Utility {#utility}


## Change of Variable

The Elmer/Ice User Function [USF_CoV.F90](../../UserFunctions/USF_CoV.F90)
provides a collection of user function for standard change of varaibles used
with the optimisation, *e.g.* to guarantee that a varaible stay positive.

Below are the available functions:

```
!#  => B=10^alpha
*Parameter* = Variable Alpha
  REAL procedure "ElmerIceUSF" "TenPowerA"

!# the derivative of the function above
!# => DB/Dalpha = ln(10)*10.0^(alpha)
*Parameter Der* = Variable Alpha
  Real procedure "ElmerIceUSF" "TenPowerA_d"

!# DJDAlpha=DJ/DB*DB/Dalpha
!# DJDAlpha=DJ/DB * ln(10)*10.0^(alpha)
DJDAlpha = Variable DJDB, alpha
  Real procedure "ElmerIceUSF" "Derivative_TenPowerA"

!# alpha=log10(B=10^alpha)
Alpha = Variable B
  Real procedure "ElmerIceUSF" "Log10A"

```

```
!#  => alpha^2
*Parameter* = Variable Alpha
  REAL procedure "ElmerIceUSF" "Asquare"

!# the derivative of the function above
!# => DB/Dalpha = 2 * alpha
*Parameter Der* = Variable alpha
  Real procedure "ElmerIceUSF" "Asquare_d"

!# DJDAlpha=DJ/DB*DB/Dalpha
!# DJDAlpha=DJ/DB * 2 * alpha
DJDAlpha = Variable DJDB, alpha
  Real procedure "ElmerIceUSF" "Derivative_Asque"

!# alpha=sqrt(B=alpha^2)
Alpha = Variable B
  Real procedure "ElmerIceUSF" "SQRTA"

```
