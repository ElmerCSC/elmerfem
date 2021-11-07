# User Function IceProperties
## General Information
- **USF Fortran File:** USF_IceProperties.f90
- **USF Name:** IceConductivity, IceCapacity and IcePressureMeltingPoint
- **Required Input Variable(s):** Temperature (IceConductivity, IceCapacity), Pressure (IcePressureMeltingPoint)

## General Description
The aim of these user functions is to provide a Fortran version of the else-as-MATC-functions-prescribed material parameters for ice (except for viscosity, which is handled by Glen's flow law). Fortran functions are way faster in execution time, which, in a run that repeatedly calls those parameters, can lead to tremendous speed-ups. Hence, if computing thermo-mechanically coupled problems, rather stick to this USF.

All input is expected to be in SI units (Temperature in Kelvin). All outputs by default are in SI units, except the user provides scaling factors (see below).

### IceConductivity
The heat conductivity of ice as a function of temperature (T) is defined (in SI units) as:
 <img src="https://latex.codecogs.com/svg.latex?\Large&space;\kappa_{ice} = 9.828 \cdot exp(-5.7^{-3} \cdot T) [W m^{−1} K^{−1}]}" title="heat conductivity of ice" />
### IceCapacity
The capacity of ice as a function of temperature (T) is defined (in SI units) as:
 <img src="https://latex.codecogs.com/svg.latex?\Large&space;c = 146.3 + (7.253 \cdot T) [J kg^{−1} K^{−1}]" title="Heat capacity of ice" />
### IcePressureMeltingPoint
The pressure melting point of ice as a function of pressure (p) is defined as (in Kelvin):
<img src="https://latex.codecogs.com/svg.latex?\Large&space;T_{pmp} = 273.15 - C_{cc} \cdot max(p, 0) [K]" title="Clausius Clapeyron relation" />
where *C_{cc}* is the Clausius Clapeyron constant. In case of negative ice pressures (actually, any below atmospheric pressures), the function uses the reference value at atmospheric pressure.

## SIF contents
The required keywords in the SIF file for these user functions are:

```
$secondsperyear = 365.25 * 24.0 * 3600.0

Constants
   Clausius Clapeyron = Real 9.8e-08
End

Material 1
  Name = "ice"

  ! Heat transfer stuff (converted to MPa-m-a system)
  Temp Heat Capacity = Variable Temp
    Real Procedure "ElmerIceUSF" "IceCapacity"
  Heat Capacity Scaling Factor = Real $(secondsperyear)^(2.0)

  Temp Heat Conductivity = Variable Temp
    Real Procedure "ElmerIceUSF" "IceConductivity"  
  Heat Conductivity Scaling Factor = Real $(secondsperyear)*1.0E-06 
  
  Temp Upper Limit = Variable HydroPressure
    Real Procedure "ElmerIceUSF" "IcePressureMeltingPoint"
  Pressure Scaling Factor = Real 1.0E06 ! from MPa to Pa
End
```

## Examples
An example demonstrating the use of the thermal properties of ice can be found in [ELMER_TRUNK]/elmerice/Tests/TemperateIceTestFct.

## References
Ritz, C. 1987. Time dependent boundary conditions for calculation of temperature fields in ice sheets. In: E. D. Waddington and J. S. Walder (Eds.), The Physical Basis of Ice Sheet Modelling, IAHS Publication No. 170, pp. 207–216. IAHS Press, Wallingford, UK.
