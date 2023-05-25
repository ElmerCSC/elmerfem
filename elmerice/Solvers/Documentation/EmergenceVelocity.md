# Retrieve Emergence Velocity
## General Information
- **Solver Fortran File:** Emergence.F90
- **Solver Name:** GetEmergenceVelocity
- **Required Output Variable(s):**
  - (1) EmergenceVelocity
- **Required Input Variable(s):**
  - (1) Normal Vector
  - (2) Flow Solution
- **Optional Output Variable(s):** None
- **Optional Input Variable(s):** None

## General Description
This is a pseudo solver (i.e., it only composes a new variable from given, not solving a matrix) to retrieve the emergence velocity from a given surface normal vector and a velocity field. It uses the fact that the scalar-product between surface velocity and surface normal gives (apart from a factor that is very close to unity) the emergence velocity

*v_{em} =  - u  {{\partial h}/{\partial x}} - v  {{\partial h}/{\partial y}} + w,*

as given in the kinematic boundary condition of the free surface.

*{{\partial h}/{\partial t}} - v_{em} = a ||grad F_h||,*

where *h* is the z-coordinate of the free surface, *(u,v,w)* the components of the ice velocity vector *\vec{u},* a the net normal accumulation/ablation. The gradient of the implicit free surface function *F_h = z - h*

*||grad F_h|| = \sqrt{({{\partial h}/{\partial x}})^2 + ({{\partial h}/{\partial y}})^2 + 1},*

usually is approximated by unity, as the derivatives of the free surface equation are of the order of the aspect ratio (usually small). Consequently, the surface normal is given by

*\vec{n} = {{grad F_h}/{||grad F_h||}} \approx{ {grad F_h}}*

and hence *\vec{n} =( -{{\partial h}/{\partial x}}, -{{\partial h}/{\partial y}}, 1 ),*

and the emergence velocity can be approximated by

*v_{em} = \vec{u} . \vec{n}*

## SIF Contents
The following SIF excerpt additionally contains solvers needed for the surface Normal Vector (Solver 2) and the *Flow Solution* (not shown here).

```
!==============================================================================
! /// Compute Normals ///
!==============================================================================
Solver 2
   Exec Solver = "Before Simulation"
   !Exec Solver = Never
   Equation = "NormalVector"
   Procedure = "ElmerIceSolvers" "ComputeNormalSolver"
   Variable = String "Normal Vector"
   Variable DOFs = 3
   Optimize Bandwidth = Logical False
   ComputeAll = Logical False
End
!==============================================================================
! /// Computing emergence velocity ///
!==============================================================================
Solver 3
  Equation = "SMB"
  Procedure = "ElmerIceSolvers" "GetEmergenceVelocity"
  Variable = -dofs 1 EmergenceVelocity
End
```
GetEmergenceVelocity (Solver 4) is - in contrary to the other Solvers - executed only on the body declared at the free surface boundary, where ComputeNormal = Logical True has to be set. The corresponding Equation 2 has to contain a keyword Convection = “Computed” as well as the name of the variable contining the velocities (usually Flow Solution)

```
Equation 2
  Name = "Surface Equations"
  Active Solvers(1) = 4 
  Convection = String "Computed"
  Flow Solution Name = String "Flow Solution"
End
```
## Example
An example solving the emergence velocity for a surface velocity and shape distribution computed on a Bueler-profile can be found in [ELMER_TRUNK]/elmerice/Tests/Emergence.

## References
Välisuo, I., T. Zwinger, and J. Kohler, 2017. Inverse solution of surface mass balance of Midtre Lovénbreen, Svalbard, Journal of Glaciology, 1-10, doi:10.1017/jog.2017.26
