# Solver Grounded Solver
## General Information
- **Solver Fortran File:** GroundedSolver.f90
- **Solver Name:** GroundedSolver
- **Required Output Variable(s):** GroundedMask (user defined)
- **Required Input Variable(s):** Need to read the bedrock, which can be either a variable or a material parameter (default is the material parameter 'Min Zs Bottom')
- **Optional Output Variable(s):** None
- **Optional Input Variable(s):** None

## Versions
The solver GroundedSolverInit has been suppressed (version 7078). Instead, the solver GroundedSolver should be executed twice (one before the simulation and one at each time step).

## General Description
The aim of these solvers is to compute a mask to qualify which part of the bed is grounded, which one is floating and which one belongs to the grounding line. The GroundedMask value is +1 where grounded, -1 where floating and 0 on the grounding line.

## SIF contents
The required keywords in the SIF file for this solver are:

To initialise the GroundedMask variable

```
Solver 1
  Exec Solver = Before All
  Equation = "GroundedMask"
  Variable = "GroundedMask"
  Variable DOFs = 1
  Procedure = "ElmerIceSolvers" "GroundedSolver"
  ! Give a tolerance for the bedrock
  Toler = Real 1.0e-03
  
  ! DEFAULT: Bedrock is read in the material parameter "Min Zs Bottom"
  !! OR Use this keyword if the bedrock is given as a variable
  ! Bedrock Variable = String "bedrock"
  !! OR use this keyword if the bedrock is given as a material parameter
  ! Bedrock Material = String "Min Zb"  
End
```
To get the new value of GroundedMask. To be executed just after the mesh update.

```
Solver 15
  Equation = "GroundedMask"
  Variable = "GroundedMask"
  Variable DOFs = 1
  Procedure = "ElmerIceSolvers" "GroundedSolver"

  Toler = Real 1.0e-3

  ! DEFAULT: Bedrock is read in the material parameter "Min Zs Bottom"
  !! OR Use this keyword if the bedrock is given as a variable
  ! Bedrock Variable = String "bedrock"
  !! OR use this keyword if the bedrock is given as a material parameter
  ! Bedrock Material = String "Min Zb" 
End
```
The bedrock elevation is read through the keyword Min Zs Bottom which is in the Material section. Here bedrock is a variable which contains the bedrock elevation.

```
Material 1
  ...
  Min Zs Bottom = Equals bedrock
  Max Zs Bottom = Real 1.0e6
End
```

## Examples
An example testing the Grounded solver can be found in [ELMER_TRUNK]/elmerice/Tests/Grounded. See also GL_MISMIP and Contact (also in the Tests directory).
