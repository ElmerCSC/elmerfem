# User Function USF_LateralFriction
## General Information
- **USF Fortran File:** USF_LateralFriction.f90
- **USF Name:** LateralFriction_x and LateralFriction_y
- **Required Input Variable(s):** A Flow solution in FlowSolverName

## General Description
This user function return the modified gravity force *-g + K * |u|^(m-1) u*
where *K* is the a lateral friction coefficient and *m* the lateral friction exponent, as in Gagliardini et al. (2010). When *m*=0, this is equivalent to the [shape factor](./ShapeFactor.md) (no dependency to the velocity). How to evaluate the lateral friction coefficient *K* as a function of the glacier width *W* is given in the supplementary material of Gagliardini et al. (2010). Noting that *K* as defined here is equal to *K/rho_i* in the paper, it becomes *K = (n+1)^(1/n) / [rho_i W^(1+1/n) (2A)^(1/n)]* and *m=1/n* (*A* and *n* are Glen's law parameter and exponent, respectively).

This solver works only in 2D (no sense in 3D). It works for non-structured mesh.

## SIF contents
The required keywords in the SIF file for these user functions are:

```
$yearinsec = 365.25*24*60*60
$rhoi = 900.0/(1.0e6*yearinsec^2)
$gravity = 9.81*yearinsec^2
$Afactor = 80.0
$n = 3.0
$etai = 1.0/(2*Afactor)^(1.0/n)
$W = 10.0e3
$Kspring = etai* (n+1)^(1/n) / (rhoi * W^(1+(1/n)))

Body Force 1
  Flow BodyForce 1 = Variable Coordinate 1
     Real Procedure "ElmerIceUSF" "LateralFriction_x"
  Flow BodyForce 2 = Variable Coordinate 1
     Real Procedure "ElmerIceUSF" "LateralFriction_y"

  Lateral Friction Gravity 1 = Real 0.0    
  Lateral Friction Gravity 2 = Real $-gravity
  Lateral Friction Coefficient = Real $Kspring
  Lateral Friction Exponent = Real $1.0/n
  Flow Solver Name = String Flow Solution
End
```

## Examples
An example using the user function USF_LateralFriction can be found in [ELMER_TRUNK]/elmerice/Tests/LateralFriction.

## Reference
When used this solver can be cited using the following references:
Gagliardini O., G. Durand, T. Zwinger, R. Hindmarsh and E. Le Meur, 2010. Coupling of ice-shelf melting and buttressing is a key process in ice-sheets dynamics, Geophys. Res. Lett., 37, L14501, doi:10.1029/2010GL043334.
