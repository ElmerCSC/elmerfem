# Adaptive timestepping in ElmerSolver

Adaptive timestepping controls the local time discretization error by
dividing each user timestep (`Timestep Sizes`) into substeps of varying size.
The user timestep is the output interval: results are written at its end. The
substeps happen inside it, and `After Timestep` solvers are run after each
accepted substep.

A substep (a *try*) is solved, its error estimated and the try either accepted
or rejected. A rejected try is undone and repeated with a smaller step. The
size of the next try is computed from the error,

    ddt_new = ddt * clamp( Safety * err^(-1/(p+1)), Decrease, Increase )

where `p` is the order of the time integration method and `err` is the error
scaled so that `err <= 1` means "accept".

Undoing a try restores the complete state of the simulation at the start of
the step:
- the solver variables and their time histories;
- all other mesh variables, including the exported ones;
- the node coordinates of moving meshes;
- the UMAT state variables of ElasticSolve.


## Minimal setup

    Simulation
      Simulation Type = Transient
      Timestepping Method = BDF
      BDF Order = 2
      Timestep Sizes = 0.1
      Timestep Intervals = 10

      Adaptive Timestepping = True
      Adaptive Time Error = 1.0e-3         ! relative tolerance, required
      Adaptive Time Error Absolute = 1e-6  ! recommended when a field starts from zero
    End

All other keywords have defaults suitable for most cases.


## Error estimates

`Adaptive Error Type` selects how the error of a try is estimated.

| Value | Cost per try | Description |
|---|---|---|
| `predictor` (default) | 1 solve | The solution is compared to its extrapolation from the previous accepted steps (Milne's device). Works for BDF (and implicit Euler). Step doubling is used for the first steps, until the history is long enough, and for methods other than BDF. |
| `weighted` | 3 solves | Step doubling: one full step and two half steps, compared with the weighted norm below. |
| `norm` | 3 solves | The original estimate: step doubling, with the error the largest relative change of the solver *norms*. Kept for compatibility, it tracks the true error poorly. |

With `predictor` and `weighted`, the error is a weighted RMS norm over the
degrees of freedom of each solver:

    err = RMS( (u_a - u_b) / (atol + rtol * |u|) ) / (2^p - 1)

Here `u_a - u_b` is the difference of the two solutions compared, and
`rtol = Adaptive Time Error`. The largest error over the solvers decides.

Which solvers are included:
- Solvers that have a time derivative term (mass or damping matrix)
  are included. Quasi-static solvers such as MeshSolve or postprocessing
  solvers are not.
- Pressure components (e.g. of Navier-Stokes) are excluded by default. They
  are algebraic, so their error does not decrease with the step.
- With explicit Runge-Kutta 2, its embedded error estimate is combined in.
  This catches the onset of instability, which step doubling alone does not.

In parallel runs each shared degree of freedom is counted once, and all
partitions take the same steps as a serial run.


## Keywords

### Simulation section

| Keyword | Type | Default | Meaning |
|---|---|---|---|
| `Adaptive Timestepping` | Logical | False | Turns adaptive timestepping on. |
| `Adaptive Time Error` | Real | (required) | Relative error tolerance `rtol`. |
| `Adaptive Time Error Absolute` | Real | see below | Absolute error tolerance `atol`, in the units of the fields. |
| `Adaptive Error Type` | String | `predictor` | `predictor`, `weighted` or `norm`, see above. |
| `Adaptive Safety Factor` | Real | 0.9 | Safety factor of the step size formula. |
| `Adaptive Increase Coefficient` | Real | 2.0 | Largest growth factor of the step. |
| `Adaptive Decrease Coefficient` | Real | 0.5 | Smallest reduction factor of the step. |
| `Adaptive Max Timestep` | Real | `dt` | Largest substep. It is never larger than the user timestep `dt`. |
| `Adaptive Min Timestep` | Real | 1e-8 × max | Smallest substep. A try at the minimum size is accepted whatever its error, with a warning. |
| `Adaptive Keep Smallest` | Integer | 0 | After a rejection, keep the reduced size for this many steps before letting it grow. |
| `Adaptive Initial Timestep` | Real | – | Size of the first substep. If not given, it is found by the probe step. |
| `Adaptive Initial Probe` | Logical | True | Estimate the first substep by a discarded trial step of 1e-3 × max size. False starts from the max size instead. |
| `Adaptive Timestep Stretch` | Real | 0.1 | How much a step may be stretched to reach the end of a user timestep (see below). |
| `Adaptive Smart Timestep` | Logical | – | Old ways at the end of a user timestep: True splits the rest into two equal steps, False takes it as is. Giving it at all turns off the default look-ahead. |
| `Adaptive Rough Timestep` | Logical | False | Do not shorten the last substep to end exactly at the user timestep. The solution is labelled with its true time, and the next timestep is shortened by the overshoot. |
| `Adaptive Reject Not Converged` | Logical | False | Reject a try also when a solver's nonlinear iteration ends at its iteration limit without converging (see below). A diverged iteration always rejects the try. |

### Solver section

| Keyword | Type | Default | Meaning |
|---|---|---|---|
| `Adaptive Time Error` | Real | Simulation value | Relative tolerance for this solver. |
| `Adaptive Time Error Absolute` | Real | Simulation value | Absolute tolerance for this solver. |
| `Adaptive Error Exclude` | Logical | – | True leaves the solver out of the error estimate. False forces it in, for solvers that integrate in time on their own. |
| `Adaptive Error Components` | Integer array | all but pressure | The components of a vector variable that are included. |

### Set by the code (not by the user)

| Keyword | Section | Meaning |
|---|---|---|
| `Adaptive Error Measure` | Simulation | Error estimate supplied by a solver (Runge-Kutta 2 sets it). |
| `Adaptive Timestep Request` | Simulation | Largest step a solver accepts for the current try (absolute time). If the try was longer it is rejected, and the next steps are limited to it. ElasticSolve sets it from the `PNEWDT` of a UMAT. |


## Absolute tolerance

If `Adaptive Time Error Absolute` is not given, `atol` is `rtol` times the
largest magnitude of each component, with a floor of 1e-6 times the largest
component. This works for fields of a definite size, such as temperature.

For a field that starts from zero, such as velocity starting from rest, give
`atol` explicitly. A field growing from zero changes by a large *relative*
amount on every step at the start. Without an `atol` the steps then shrink
towards the minimum size.


## End of a user timestep

The last substep is shortened to end exactly at the user timestep (unless
`Adaptive Rough Timestep` is set). To avoid a sliver of a step there, the
`predictor` and `weighted` types look ahead:
- If the rest of the timestep is at most `1 + Adaptive Timestep Stretch` times
  the step, it is taken as one slightly longer step. This is not done right
  after a rejection, nor beyond `Adaptive Max Timestep` or a step size request.
- If the rest is less than two steps, it is taken as two equal steps.
- A step that was shortened only to reach the end of the timestep does not
  limit the next one. The next timestep continues from the step size before
  the shortening (growth still limited by `Adaptive Increase Coefficient`).


## Nonlinear divergence

A try is also rejected when the nonlinear iteration of a solver diverges. The
rest of the try's solves are skipped, and the step is reduced by `Adaptive
Decrease Coefficient`. The log line of the try is marked `(diverged)`. A
solver reports divergence:
- through `Nonlinear System Divergence Limit` (any solver using the standard
  convergence test);
- in ElasticSolve with a UMAT, when the residual is not finite or an element
  is inverted (`det F <= 0`). The iteration stops there, before the next
  assembly. Without adaptive timestepping or `Convergence Control` this is
  fatal;
- when an iterative linear solver diverges and `Linear System Abort Not
  Converged = False` (otherwise this is fatal already).

An iteration that only reaches `Nonlinear System Max Iterations` without
converging is accepted with a warning. With `Adaptive Reject Not Converged =
True` it is rejected too. This is useful when the step size is limited by the
nonlinear iteration rather than by the time discretization error. Examples are
quasi-static plasticity driven through time, or any case without mass, where
the error estimate is close to zero. Use it only with solvers that have a
`Nonlinear System Convergence Tolerance`: without one, a solver never reports
convergence.

If a try diverges even at `Adaptive Min Timestep`, it is accepted with a
warning.

The log line of a try shows its number of nonlinear iterations as `iter n`.

A limit on the number of nonlinear iterations, as in Abaqus, does not work
for such cases: at the onset of yield or of necking even tiny steps take as
many iterations as long ones, and the step would be held small for nothing.
For the same reason, do not set `Nonlinear System Max Iterations` low to make
failures cheaper. In a necking tensile test the steps through the neck needed
10-19 iterations whatever their size; with a limit of 10 the steps shrank
tenfold and the run took twice as long as with 25.

### Line search with a UMAT

Near a limit load a full Newton step of a plastic material often overshoots,
typically into inverted elements, and the try diverges although a solution
exists. Then the steps shrink far below what the error tolerance would allow.
In ElasticSolve with a UMAT, a line search avoids this:

    Nonlinear System Linesearch = True

After each Newton solve the full step `du` is tried first. It is judged by the
slope `s(alpha) = du . R(u + alpha*du)`, the derivative of the potential
energy along the step, as in Crisfield's method (and Abaqus). The residual
norm would be a poor measure: a Newton step of a bending structure raises
large residuals across `du`, such as axial forces, which the iteration then
removes by itself. The step is shortened, and the iterate assembled again
without a new linear solve, when
- an element is inverted or the residual is not finite: `alpha` is reduced by
  `Linesearch Factor`;
- the step went beyond the energy minimum along `du`, `s(alpha) < -Limit *
  s(0)`: `alpha` is set to the root of the secant of `s`.

A step short of the minimum is accepted, and after the given number of
reductions the last step is taken as it is. The UMAT always starts from the
state of the beginning of the step, so a rejected trial leaves nothing behind.

| Keyword (Solver section) | Type | Default | Meaning |
|---|---|---|---|
| `Nonlinear System Linesearch` | Logical | False | Turns the line search on. |
| `Nonlinear System Linesearch Iterations` | Integer | 10 | Largest number of step reductions per Newton iteration. |
| `Nonlinear System Linesearch Factor` | Real | 0.5 | Reduction of the step length after an inverted element or a non-finite residual. |
| `Nonlinear System Linesearch Limit` | Real | 0.5 | How far beyond the energy minimum a step may go, as a fraction of the initial slope. |

In a quasi-static plasticity case driven through time this allowed steps
about 60 times longer, with the same results. The test `UMAT_linesearch`
applies a load in one step that plain Newton cannot.

The residual of the UMAT Newton iteration is relative to the external load
vector. When there is none, as in a case driven by prescribed displacements
only, it is relative to the reaction forces of the current iterate instead
(the residual at the constrained degrees of freedom). `Nonlinear System
Convergence Tolerance` is then a fraction of the forces carried, whatever
the units. `Calculate Loads` does not change this. Near a limit load the residual stops decreasing at some
level, from the line search and from points switching between loading and
unloading. A tolerance below that level is never met, and every try is
rejected as not converged. In a tensile test through necking, 1e-2 passed
smoothly where 1e-4 and stricter stalled; CalculiX accepts 0.5 % of the
average force by default.

`Nonlinear System Relaxation Factor` has no effect with a UMAT: ElasticSolve
adds the Newton increment to the total solution itself. If the factor is
given with a value other than 1, a warning is printed once. Use the line
search instead.

### UMAT tangent with large deflection

With `Large Deflection = True` (the default of ElasticSolve) the UMAT gets the
deformation gradients `DFGRD0` and `DFGRD1`. It returns the Cauchy stress and
its tangent `DDSDDE`. The tangent can be defined in two ways, and ElasticSolve
must be told which one the UMAT uses. It is given in the Material section:

| Keyword (Material section) | Type | Default | Meaning |
|---|---|---|---|
| `UMAT Tangent` | String | `Material` | `Material` or `Jaumann` (alias `Abaqus`), see below. |

- `Material`: `DDSDDE` is the derivative of the Cauchy stress with respect to
  the Green-Lagrange strain. This is the original convention of ElasticSolve.
- `Jaumann`: `DDSDDE` follows the Abaqus convention. `J * DDSDDE` relates
  the Jaumann rate of the Kirchhoff stress `tau = J * sigma` to the rate of
  deformation `d = sym(L)`:

      tau_dot = J * DDSDDE : d + w*tau - tau*w,   L = d + w

  ElasticSolve adds the spin terms `w*tau - tau*w` and the geometric terms of
  the stress itself. Use this for UMATs written for Abaqus.

The stress is the same either way. Only the Newton iteration changes. With the
wrong convention the tangent does not match the stress, and the iteration
converges only linearly. With large rotations or near the limit load it may
not converge at all. In a test with a plastic beam stretched to 8 % strain
(an Abaqus-convention UMAT), `Material` failed after a few steps, while
`Jaumann` converged in 9-15 iterations per step. For a purely elastic step of
0.5 % strain, `Material` took 20 iterations and `Jaumann` 6.

Without `Large Deflection` the two conventions are the same, and the keyword
has no effect.

### Locking of linear elements: UMAT Bbar

Plastic flow, and a nearly incompressible material, keeps the volume constant.
A fully integrated linear element has to keep it at each integration point,
which it can only do by hardly deforming. The mesh is then too stiff
(volumetric locking). Quadratic elements, such as 820 with
`Element Integration Points = "-brick 8"` (the C3D20R of Abaqus and CalculiX),
do not lock. For linear elements, set in the Material section:

| Keyword (Material section) | Type | Default | Meaning |
|---|---|---|---|
| `UMAT Bbar` | Logical | False | Take the volume change of the linear elements (404, 706, 808) from the whole element. |

- Small strain (`Large Deflection = False`): B-bar, the mean dilatation of
  Hughes (1980). The volumetric part of the strain is replaced by its average
  over the element. The tangent stays symmetric.
- Large deflection: F-bar (de Souza Neto et al., 1996). At each integration
  point the UMAT gets `F_bar = (J_c/J)^(1/3) F`, where `J_c` is the volume
  change at the element centre. The tangent gets a term coupling to the
  centre and is not symmetric. Both `UMAT Tangent` conventions are supported.

This is what the Abaqus C3D8 element does as well. It is not reduced
integration: `Element Integration Points = "-brick 1"` on 808 leaves the
hourglass modes free, and the iteration fails.

Results with 808 bricks:

| Case | Without | `UMAT Bbar` | 820 reference |
|---|---|---|---|
| Cantilever tip, nu = 0.4999, small strain | -0.0079 | -0.0366 | -0.0383 |
| Same, distorted mesh | -0.0028 | -0.0358 | |
| Plastic bar at 10 % elongation, force | 461.7 | 449.7 | 450.1 |

In an axisymmetric thick cylinder with nu = 0.4999 (404 elements), the inner
displacement was 28 % of the exact value without and 99.9 % with it. The
Newton iteration stays quadratic. The patch test (a uniform strain on a
distorted mesh) is passed exactly.

Notes:

- Wedges (706) gain little. Their triangular cross section locks like a linear
  triangle. In the bar above the force went from 463.7 to 460.1.
- Linear triangles and tetrahedra have one integration point, and quadratic
  elements do not lock. For them the keyword has no effect.
- Shear locking in bending is not affected. Use several linear elements
  through the thickness, or quadratic elements.
- With F-bar the first Newton iterates of a very large step can be rougher.
  A 15 % uniform axisymmetric stretch in one step from rest did not converge
  with `UMAT Tangent = Material`. In four steps it converged.

### Quasi-static cases

For a case without inertia, give the solver

    Time Derivative Order = 1
    Quasi Stationary = True

Then the time derivative is of first order and the `predictor` estimate can
be used, at one solve per try. The estimate then measures how well the
solution follows its extrapolation along the load path. Without these the
solver is of second order and step doubling is used, at three solves per try.


## Notes and limitations

- The `predictor` estimate needs BDF (or implicit Euler). With Newmark, Bossak
  or Crank-Nicolson the step doubling estimate is used.
- Every try with step doubling costs three solves, with `predictor` one solve.
  Coupled systems repeat their coupling iterations on every solve.
- ElasticSolve with a UMAT supports transient and adaptive runs. The UMAT may
  ask for a smaller step via `PNEWDT < 1`. Without adaptive timestepping this
  only gives a warning.
- Related but separate mechanisms, not covered here: `Convergence Control`
  (step reduction on divergence) and `Predictor-Corrector Control`.


## Test cases

| Test | Shows |
|---|---|
| `TimeAdapt` | the legacy `norm` estimate |
| `TimeAdaptBDF2`, `TimeAdaptBDF2_mpi` | `weighted` step doubling with BDF2, serial and parallel |
| `TimeAdaptPredictor`, `TimeAdaptPredictor_mpi` | the default `predictor` estimate |
| `TimeAdaptRK2` | explicit Runge-Kutta 2 with its embedded estimate |
| `TimeAdaptRough` | `Adaptive Rough Timestep` |
| `TimeAdaptMeshUpdate` | a moving mesh, state restored on rejection |
| `UMAT_viscoelastic_adaptive` | a UMAT with state variables and `PNEWDT` |
| `UMAT_linesearch` | the line search of the UMAT Newton iteration |
