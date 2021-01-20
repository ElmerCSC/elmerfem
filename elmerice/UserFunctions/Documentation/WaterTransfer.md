# User Function USF_WaterTransfer
## General Information
- **USF Fortran File:** USF_WaterTransfer.f90
- **USF Name:** EPLToIDS
- **Required Input Variable(s):** EPLHead, IDSHead and IDSHead Upper Limit
- **Optional Input Variable(s):** None

## General Description
The aim of this function is to compute a water transfer between the two layers of the hydrological Double Continuum Equivalent model.

## SIF contents
The required keywords in the SIF file for these user functions are given in the presentation of the efficient drainage system solver (EPLSolver).

## Examples
This transfer function is used in test in [ELMER_TRUNK]/elmerice/Tests/Hydro_Coupled.
