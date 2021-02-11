# Solver/User Function {Name of Solver/USF}
## General Information
- **Solver/User Function Fortran File:** {Name of .F90 source file}
- **Solver/User Function Name:** {Name of actual solver/USF subroutine}
- **Required Output Variable(s):** {list all variables that the solver mandatorily outputs}
- **Required Input Variable(s):** {list all required inputs and the associated keywords}
  - {e.g. Flow solution (Flow Solution Name = String...)}
- **Optional Output Variable(s):** {list any optional output variables and associated keywords}
  - {e.g. SomeDerivedVar (Output Optional Variables = Logical True)}
- **Optional Input Variable(s):** {list any optional inputs and associated keywords}
  - {e.g. SomeInputVar (Extra Input Name = String...)}
- **Solver Keywords:** {list all other solver-specific keywords with brief explanations}
  - {e.g. LoseGame = Logical True (You have lost The Game)}
  
## General Description
{Include a brief summary of what the solver/USF does, with equations if they help. For some solvers/USFs this will be briefer than for others. If there are specific things that need to be addressed or explained, add additional subsections}

### Additional subsection 1 {e.g. Geographical Restriction}
{In Soviet Russia, this solver solves YOU}

## Known Bugs and Limitations
{If there are any known issues or limitations that could be improved, detail them here}
- {e.g. 2019-02-22 Lunar phase issue fixed}
- {e.g. 2017-05-14 This solver only works at full moon}
- {e.g. Does not work in parallel}

## SIF Contents
{Provide examples of all the required SIF contents for a typical use of the solver/USF. Please ensure you include examples of lines that need to be added to other parts of the SIF (BCs, Material, etc.), not just the solver/USF block}
The required keywords in the SIF file for this solver/USF are:

```
Solver 1
  Equation = ...
  Variable = ...
  Variable DOFs = ...

  Procedure = ...
End

Material 1
...
End

Boundary Condition 1
...
End
```
## Examples
{Show where a relevant example in the Elmer distribution can be found}
An example in which the ... can be found here [ELMER_TRUNK]/elmerice/Tests/...

## References
{If the solver is derived from a published source or is the main subject of one, provide a reference here}
