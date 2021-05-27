# User Function Coulomb Friction Law
## General Information
- **USF Fortran File:** USF_Sliding.f90
- **USF Name:** Friction_Coulomb
- **Required Input Variable(s):** A Flow Solution in Flow Solution Name, Normal Vector, Stress or the Effective Pressure variable.

## General Description
The file USF_Sliding.f90 contains three user functions to apply non-linear friction at the base of glacier.

The first user function (Sliding_Weertman) is a non-linear Weertman-type friction law and is described [here](./Weertman.md). The second user function (Friction_Coulomb) is a non-linear water-pressure-dependent friction law, as proposed by Schoof (2005) and Gagliardini et al. (2007), and is presented in this page. The third user function (Sliding_Budd) is described [here](./Budd.md) and is from Budd et al 1984 (Annals of Glaciology 5, page 29-36).

The friction law in Friction_Coulomb is of the form:
*tau_b = C.N {[{ {chi . {u_b}^{-n} }/ {(1 + a . chi^q)} }]}^{1/n} . u_b*
where
*a = {(q - 1)^{q-1}}/{q^q}*
and
*chi = {u_b}/{C^n N^n A_s}*
The Slip Coefficient in Elmer is then given as
*C.N {[{ {chi . {u_b}^{-n} }/ {(1 + a . chi^q)} }]}^{1/n}*
When *u_b < u_{t0}*, *u_b* in the previous equation is replaced by *u_{t0}*.

The parameters to be given are:
- Friction Law Sliding Coefficient → *A_s*
- Friction Law Post-Peak Exponent → *q* >= 1
- Friction Law Maximum Value → *C* ~ max bed slope
- Friction Law Exponent → *m* = (n Glen's law)
- Friction Law Linear Velocity → *u_{t0}*
The effective pressure is defined as *N = -sigma_{nn} -p_w*, where *sigma_{nn}* is the normal Cauchy stress and *p_w* the water pressure. If a variable Effective Pressure exists, it is used to evaluate directly *N*. Else, the normal Cauchy stress is estimated from the stress computed at previous timestep. The water pressure is prescribed as an External Pressure (Negative - Compressive convention, and therefore 'External Pressure' should be equal to the opposite of the water pressure in the sif).

## SIF contents
The required keywords in the SIF file for this user function are:

```
!!! Bedrock Boundary Condition 
Boundary Condition 1
  Target Boundaries = 1

  Normal-Tangential Velocity = Logical True
  Flow Force BC = Logical True
  
  !! Water pressure given through the Stokes 'External Pressure' parameter 
  !! (Negative = Compressive)
  External Pressure = Equals Water Pressure
   
  Velocity 1 = Real 0.0
  
  Slip Coefficient 2 =  Variable Coordinate 1
    Real Procedure "ElmerIceUSF" "Friction_Coulomb"
  Slip Coefficient 3 =  Variable Coordinate 1
    Real Procedure "ElmerIceUSF" "Friction_Coulomb"
    
  !! Parameters needed for the Coulomb Friction Law
  Friction Law Sliding Coefficient = Real 4.1613e5  
  Friction Law Post-Peak Exponent  = Real 1.0      !(q=1)
  Friction Law Maximum Value = Real 1.0            !(C=1)
  Friction Law PowerLaw Exponent = Real 3.0        !(m = n = 3 Glen's law) 
  Friction Law Linear Velocity = Real 0.01         
End
```

## Examples
The Coulomb friction law is tested in [ELMER_TRUNK]/elmerice/Tests/Friction_Coulomb with a direct input of the effective pressure and [ELMER_TRUNK]/elmerice/Tests/Friction_Coulomb_Pw with the effective pressure computed from the stress and a prescribed water pressure.

## Reference
When this friction law is used, it can be cited using the following reference:
Gagliardini O., D. Cohen, P. Råback and T. Zwinger, 2007. Finite-Element Modeling of Subglacial Cavities and Related Friction Law. J. of Geophys. Res., Earth Surface, 112, F02027.
