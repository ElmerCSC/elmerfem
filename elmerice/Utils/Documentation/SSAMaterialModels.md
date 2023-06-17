# Module  SSAMaterialModels

## General Information  
- **Module Fortran File:** SSAMaterialModels.F90
 
## History
- Creation date: May 2022
	- Move SSA friction laws in a separate module to ease use by other piece of code
- Rev. afcafe865: 
	- add ComputeMeanFriction
- 4th June 2023
      	- add SSAEffectiveBMB

## General Description  
Module containing utility functions to define the Material/Frictions laws for the SSA.

## Available routines

- SSAEffectiveFriction:
	- Return the effective friction coefficient at current location

- ComputeMeanFriction:
	- compute the element-averaged basal friction

- SSAEffectiveBMB:
	- Return the basal mass balance at current location.  Requires the following.

In material section:

SSA Melt Param = String "sem"

Valid values include sem for sub element melt parameterisation (i.e. melt at each IP is allowed or not based on floatation), fmp for full melt parameterisation and nmp for no melt parameterisation.

In Thickness Solver:

GL integration points number = Integer 20
