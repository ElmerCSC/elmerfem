# User Function Weertman Friction Law
## General Information
- **USF Fortran File:** USF_Sliding.f90
- **USF Name:** Sliding_Weertman
- **Required Input Variable(s):** A Flow Solution in Flow Solution Name, Normal Vector

## General Description
The file USF_Sliding.f90 contains user functions to apply non-linear friction at the base of glacier.

The first user function (Sliding_Weertman) is a non-linear Weertman-type friction law and is described in [this page](./Weertman.md). The second user function (Friction_Coulomb) is a non-linear water-pressure-dependent friction law, as proposed by Schoof (2005) and Gagliardini et al. (2007), and is presented [here](./Coulomb.md). The third user function (Sliding_Budd) is described here and is from Budd et al 1984 (Annals of Glaciology 5, page 29-36).

The friction law in Weertman_Sliding is of the form:
*tau_b = C.{u_b}^{m - 1} . u_b*
The Slip Coefficient in Elmer is then given as
*C.{u_b}^{m - 1}*
When *u_b < u_{t0}*, *u_b* in the previous equation is replaced by *u_{t0}* (linearisation for small velocity).

The parameters to be given are:
- Weertman Friction Coefficient → *C*
- Weertman Exponent → *m*
- Weertman Linear Velocity → *u_{t0}*

## SIF contents
The required keywords in the SIF file for this user function are:

```
!!! Bedrock Boundary Condition 
Boundary Condition 1
  Target Boundaries = 1

  Normal-Tangential Velocity = Logical True
  Flow Force BC = Logical True

  Velocity 1 = Real 0.0
  
  Slip Coefficient 2 =  Variable Coordinate 1
    Real Procedure "ElmerIceUSF" "Sliding_Weertman"
  Slip Coefficient 3 =  Variable Coordinate 1
    Real Procedure "ElmerIceUSF" "Sliding_Weertman"
    
  Weertman Friction Coefficient = Real 2.412579e-2        
  Weertman Exponent = Real $1.0/3.0
  Weertman Linear Velocity = Real 0.00001

End
```

## Examples
The Weertman friction law is used in the tests [ELMER_TRUNK]/elmerice/Tests/GL_MISMIP and [ELMER_TRUNK]/elmerice/Tests/Contact and especially tested in [ELMER_TRUNK]/elmerice/Tests/Friction_Weertman.
