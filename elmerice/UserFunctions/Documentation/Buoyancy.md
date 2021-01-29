# User Function Buoyancy
## General Information
- **USF Fortran File:** Buoyancy.f90
- **USF Name:** SeaPressure and SeaSpring
- **Required Input Variable(s):** None

## General Description
The aim of these user functions is to apply the water hydrostatic pressure induced by the ocean at the base and at the front of floating ice mass.

The hydrostatic water pressure exerted by the ocean on the floating ice mass is:
*sigma_{nn} = -rho_w. g . (h_{sl} - z_b(t))*
To get a stable solution, one needs to account for the fact that the bottom free surface elevation *z_b* is evolving with time, as:
*z_b(t) = z_b(t-dt) + (u_n + b).N_s.dt*
where *u_n* is the normal velocity, *b* the normal basal melt/accretion rate (positive for melting) and *(u_n + b).N_s* the vertical projection of these quantities.

Accounting for the bottom free surface displacement, the normal stress reads:
*sigma_{nn} = -rho_w. g . (h_{sl} - b.N_s.dt - z_b(t-dt) ) + (rho_w. g.N_s.dt).u_n*
where the term *rho_w. g.N_s.dt* can be assimilated as a normal viscous spring accounting for any shift of the free surface from the hydrostatic equilibrium.

The first user function (SeaPressure) is used to apply the hydrostatic water pressure for a given sea level. Basal melt is only accounted for where the viscous spring (see below) is also applied.

The second user function (SeaSpring) evaluates the viscous spring induced by any shift of the free surface from the hydrostatic equilibrium.

In case of basal melting, the value of the basal melting is read from the Accumulation keyword of the bottom free surface.

## SIF contents
The required keywords in the SIF file for these user functions are:

```
$yearinsec = 365.25*24*60*60
$rhoi = 900.0/(1.0e6*yearinsec^2)
$rhow = 1000.0/(1.0e6*yearinsec^2)

Constants
! For the Buoyancy User function
  Buoyancy Use Basal Melt = Logical True
  Bottom Surface Name = String "Zs Bottom"
  Water Density = Real $rhow 
End

! Body force for the bottom free surface
Body Force 3
  ...
!! melting/accretion under ice/shelf
!! positive for melting
!! negative for accretion
  Zs Bottom Accumulation = Real 0.5e0
End

Material 1
  ...
  Sea level = Real 0.0
End

!! vertical front (air and sea contact)
Boundary Condition 2
  Name = "front"
  Target Boundaries = 2

  Flow Force BC = Logical True
  External Pressure = Variable Coordinate 3
     Real Procedure "ElmerIceUSF" "SeaPressure"

  Compute Sea Pressure = Logical True
End

!! Bottom BC (Sea contact)
Boundary Condition 1
  Name = "bottom"
  Target Boundaries = 1
  Body Id = 3

  Normal-Tangential Velocity = Logical True
  Flow Force BC = Logical True
  External Pressure = Variable Coordinate 3
   Real Procedure "ElmerIceUSF" "SeaPressure"

  Slip Coefficient 1 = Variable Coordinate 3
    Real Procedure "ElmerIceUSF" "SeaSpring"

  Compute Sea Pressure = Logical True
  Compute Sea Spring = Logical True
End
```

## Examples
An example of the usage of the user function SeaPressure can be found in the the Tête Rousse application of the Elmer/Ice course material.

Another example can be found in [ELMER_TRUNK]/elmerice/Tests/Buoyancy. This example is simply a floating iceberg with basal melting equals to the surface accumulation. Obviously, if the solution is correct, this iceberg should stay at the same elevation, with a vertical ice velocity equals to the melt/accumulation.
