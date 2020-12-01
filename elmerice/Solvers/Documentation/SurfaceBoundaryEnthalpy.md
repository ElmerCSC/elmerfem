#Surface Boundary Condition for steady state thermal regime
##General Information
- **Solver Fortran File:** SurfaceBoundaryEnthalpy.f90
- **Solver Name:** SurfEnthBoundarySolver
- **Required Output Variable(s):** Surf Enth
- **Optional Output Variable(s):** Mass Balance, Densi, Firn, Melting, Refreeze, Accu, Rad_Fact, Rain, PotRad
- **Required Input Variable(s):** Depth and SurfGrad from the FlowDepth Solver
- **Input Data:** Daily air temperature time series
- **Optional Input Data:** Daily precipitation time series
For vertically structured 3D mesh only. Works in serial and parallel.

##General Description
SurfEnthBoundarySolver is a pseudo-solver which compute surface mass balance and Dirichlet surface boundary condition for the Enthalpy solver. It takes into account firn heating processes by solving vertical melt-water percolation and refreezing.

The solver uses the provided air temperature (and precipitation) daily record to compute the associated mean surface characteristic of the glacier over the period covered by the provided data time series. It can output the following variables:

- Mass Balance (m w.eq./yr): Mean surface mass balance
- Surf Enth (J/kg) : Enthalpy value bellow active layer. Can be use as a Dirichlet condition in the Enthalpy Solver
- Densi (kg/m3): Density field in 3D
- Firn (m w. eq.) : Firn thickness
- Melting (m w.eq./yr) : Surface melting
- Refreeze (m w.eq./yr) : Amount of refreezing (superimposed ice)
- Accu (m w.eq./yr) : Snow accumulation
- Rad_fact (m w.eq./(W/m2)) : Melting factor for radiation
- Rain (m w.eq./yr) : Amount of rain
- PotRad (W/m2) : Potential solar radiation
The mass balance model is based on a degree day model that takes into account potential solar radiation. Mean Enthalpy at 10m-depth (below active layer), is computed by solving the heat equation on a 1D vertical profile forced by a mean annual cycle of air temperature and precipitation determined from the data. This is done for each surface node using a Crank-Nicholson scheme on a 6 cm-resolution grid at a daily time-step. It takes into account the seasonal change of the density profile and allows percolation of water only where density is lower than 800 kg/m3. More details about the model can be found in:

Gilbert, A., Sinisalo, A., Gurung, T. R., Fujita, K., Maharjan, S. B., Sherpa, T. C., & Fukuda, T. (2020). The influence of water percolation through crevasses on the thermal regime of a Himalayan mountain glacier. The Cryosphere, 14(4), 1273–1288. [https://doi.org/10.5194/tc-14-1273-2020](https://doi.org/10.5194/tc-14-1273-2020)

##SIF contents
The parameters are set in the constant section of the sif file :

```
Constants

  rho_surf = real 350.0		! Snow surface density
  rho_ice = real 917.0 		! Ice density
  rho_w = real 1000.0 		! Water density
  Sr = real 0.005		! Residual water saturation in Snow/Firn
  T_ref_enthalpy = real 200.0   ! Use to compute Surf Enth (see Enthalpy solver)
  L_heat = real 334000.0        ! Latent heat of fusion
	
  AirTemperatureFile = File "YourTempFile.dat"  ! Contain daily temperature data
  PrecipFile = File "YourPrecipFile.dat"  ! Contain daily precipitation data (optional)
    
  Precip = real 0.300  		!Mean annual precipitation if PrecipFile not provided
  TempCorrec= real -0.12	!Possibility of shifting temperature to get steady state mass balance for example (optional)
  PrecipCorrec = real 1.0	!Possibility of adding a correcting factor on precipitation if PrecipFile provided (optional)
  
  GradTemp = real 0.0065	!Air temperature Lapse Rate (K m^-1)
  GradPrecip= real 0.001	!Precipitation Lapse Rate (% m^-1)
  z_temp = real 5310.0		!Elevation of temperature measurement from AirTemperatureFile (m)
  z_precip = real 5310.0	!Elevation of Precipitation measurement (m)
  
  RadFact_ice = real 0.0000925		! Melting factor for ice from radiation
  RadFact_snow = real $0.0000925/2.0	! Melting factor for snow from radiation
  Deg_jour = real 0.0114		! Melting factor from air temperature
  
  seuil_precip = real 2.0		!Rain/Snow air temperature threshold (degree C)
  seuil_fonte = real 0.0		!Melting air temperature threshold (degree C)
  
  firn_param = real 30.0		! Firn densification factor (yr)
  super_ice = real 0.15			! Superimposed ice factor
  
  Latitude = real 28.82			!Latitude (degree) to compute Potential Solar Radiation
  
  !Possibility to export 1D profile simulation at one node of coordinate (X_output1D,Y_output1D) (optional)
  X_output1D = real x_coordinate
  Y_output1D = real y_coordinate
  

End
```
The solver needs output from the FlowDepth Solver :

```
Solver 1

  Equation = "Flowdepth"
   Procedure = File "ElmerIceSolvers" "FlowDepthSolver"
   Variable = String "Depth"
   Variable DOFs = 1
   Linear System Solver = "Direct"
   ! this sets the direction
   ! -1 is negative z-direction (upside down)
   ! +1 is positive (downside up)
   Gradient = Real -1.0E00
   Calc Free Surface = Logical True
   Freesurf Name = String "Surf"
End
```
The solver is called with the following:

```
Solver 2

  Equation = SurfBoundary
  Variable = Surf Enth
  Variable DOFs = 1
  procedure =  "ElmerIceSolvers" "SurfEnthBoundarySolver"

  ! The following variables declaration are optional:
  
  Exported Variable 1 = String "Mass Balance"
  Exported Variable 1 DOFs = 1

  Exported Variable 2 = String "Densi"
  Exported Variable 2 DOFs = 1

  Exported Variable 3 = String "Firn"
  Exported Variable 3 DOFs = 1
  
  Exported Variable 4 = String "Melting"
  Exported Variable 4 DOFs = 1

  Exported Variable 5 = String "Refreeze"
  Exported Variable 5 DOFs = 1

  Exported Variable 6 = String "Accu"
  Exported Variable 6 DOFs = 1

  Exported Variable 7 = String "Rad_Fact"
  Exported Variable 7 DOFs = 1
  
  Exported Variable 8 = String "Rain"
  Exported Variable 8 DOFs = 1
  
  Exported Variable 9 = String "PotRad"
  Exported Variable 9 DOFs = 1

End
```
Boundary Condition to setup Dirichlet condition for the Enthalpy Solver:

```
! Upper Surface
Boundary Condition 2
  Target Boundaries = 2
   
  Depth = real 0.0
  Enthalpy_h = Equals Surf Enth 

End
```

##Examples
An example solving for the enthalpy within the Rika Samba Glacier using the SurfEnthBoundarySolver can be found in [ELMER_TRUNK]/elmerice/examples/Test_SurfaceBoundaryEnth.

##Reference
Gilbert, A., Sinisalo, A., Gurung, T. R., Fujita, K., Maharjan, S. B., Sherpa, T. C., & Fukuda, T. (2020). The influence of water percolation through crevasses on the thermal regime of a Himalayan mountain glacier. The Cryosphere, 14(4), 1273–1288. [https://doi.org/10.5194/tc-14-1273-2020](https://doi.org/10.5194/tc-14-1273-2020)
