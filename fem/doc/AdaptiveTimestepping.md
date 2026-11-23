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
| `predictor` (default) | 1 solve | The solution is compared to its extrapolation from the previous
   accepted steps (Milne's device). Works for BDF (and implicit Euler). Step doubling is used
    for the first steps, until the history is long enough, and for methods other than BDF. |
| `weighted` | 3 solves | Step doubling: one full step and two half steps, compared with the weighted norm below. |
| `norm` | 3 solves | The original estimate: step doubling, with the error the largest relative
   change of the solver *norms*. Kept for compatibility, it tracks the true error poorly. |

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
| `Adaptive Min Timestep` | Real | 1e-8 × max | Smallest substep. A try at the minimum size is accepted whatever
   its error, with a warning. |
| `Adaptive Keep Smallest` | Integer | 0 | After a rejection, keep the reduced size for this many steps before letting it grow. |
| `Adaptive Initial Timestep` | Real | – | Size of the first substep. If not given, it is found by the probe step. |
| `Adaptive Initial Probe` | Logical | True | Estimate the first substep by a discarded trial step of 1e-3 × max size.
    False starts from the max size instead. |
| `Adaptive Timestep Stretch` | Real | 0.1 | How much a step may be stretched to reach the end of a user timestep (see below). |
| `Adaptive Smart Timestep` | Logical | – | Old ways at the end of a user timestep: True splits the rest into two equal steps,
   False takes it as is. Giving it at all turns off the default look-ahead. |
| `Adaptive Rough Timestep` | Logical | False | Do not shorten the last substep to end exactly at the user timestep. The solution
  is labelled with its true time, and the next timestep is shortened by the overshoot. |

### Solver section

| Keyword | Type | Default | Meaning |
|---|---|---|---|
| `Adaptive Time Error` | Real | Simulation value | Relative tolerance for this solver. |
| `Adaptive Time Error Absolute` | Real | Simulation value | Absolute tolerance for this solver. |
| `Adaptive Error Exclude` | Logical | – | True leaves the solver out of the error estimate. False forces it in, for
   solvers that integrate in time on their own. |
| `Adaptive Error Components` | Integer array | all but pressure | The components of a vector variable that are included. |

### Set by the code (not by the user)

| Keyword | Section | Meaning |
|---|---|---|
| `Adaptive Error Measure` | Simulation | Error estimate supplied by a solver (Runge-Kutta 2 sets it). |
| `Adaptive Timestep Request` | Simulation | Largest step a solver accepts for the current try (absolute time).
   If the try was longer it is rejected, and the next steps are limited to it. ElasticSolve sets it from the `PNEWDT` of a UMAT. |


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
