# Elmer NVHPC Compatibility Fixes

This document describes the changes made to the Elmer source tree to allow compilation with the NVIDIA HPC SDK Fortran compiler (`nvhpc 24.7`, as available on Levante (DKRZ) as well as `nvhpc 25.3`, as available on Dolpung (DKRZ / MPI-M)).
All fixes have been written to preserve full backward compatibility with the `gcc 11.2.0` compiler (as available on Levante (DKRZ)).


## Change 1: Module-level pointer initialization of derived types with pointer components

**Files changed:** `fem/src/SParIterGlobals.F90`, `fem/src/SParIterComm.F90`

### Compiler Error

```
NVFORTRAN-S-0034-Syntax error at or near end of line
    (fem/src/SParIterGlobals.F90: 69)
NVFORTRAN-S-0134-Illegal attribute - duplicate pointer
    (fem/src/SParIterGlobals.F90: 70)
NVFORTRAN-S-0038-Symbol, parenv, has not been explicitly declared
    (fem/src/SParIterGlobals.F90: 45)
0 inform, 0 warnings, 3 severes, 0 fatal for spariterglobals
```

### Root Cause

`nvfortran` does not support module-level pointer initialization (the `=>` syntax in a declaration) when the derived type itself contains `DIMENSION(:), POINTER` components (here `Active` and `IsNeighbour` inside `ParEnv_t`). 
`gfortran` and Intel compilers accept this as an extension, but it is not strictly standard Fortran and `nvfortran` rejects it.

### Fix

**`fem/src/SParIterGlobals.F90`**
Remove the `=> ParEnv_Common` initialiser from the `ParEnv` pointer declaration:

```fortran
! BEFORE
TYPE (ParEnv_t), SAVE, TARGET :: ParEnv_Common
TYPE (ParEnv_t), POINTER, SAVE :: ParEnv => ParEnv_Common

! AFTER
TYPE (ParEnv_t), SAVE, TARGET :: ParEnv_Common
! nvhpc cannot initialise a pointer to a derived type with pointer components
! at module level (NVFORTRAN-S-0034). ParEnv is initialised explicitly in
! ParCommInit (SParIterComm.F90) before first use.
TYPE (ParEnv_t), POINTER, SAVE :: ParEnv
```

**`fem/src/SParIterComm.F90`**
Add an explicit association guard at the top of `ParCommInit`, before `ParEnv` is first accessed:

```fortran
! BEFORE
ParallelEnv => ParEnv

! AFTER
! Ensure ParEnv points to the module-level ParEnv_Common storage.
! (With gfortran/Intel this is done via module-level pointer initialization;
!  nvhpc does not support that for derived types with pointer components.)
IF (.NOT.ASSOCIATED(ParEnv)) ParEnv => ParEnv_Common

ParallelEnv => ParEnv
```

## Change 2: Internal compiler error (ICE) on derived-type structure assignment and array operations

**File changed:** `fem/src/SParIterSolver.F90`

### Compiler Error

```
Lowering Error: array upper bound is not a symbol for datatype 36322
Lowering Error: array extnt is not a symbol for datatype 36322
NVFORTRAN-F-0000-Internal compiler error. Errors in Lowering 2
    (fem/src/SParIterSolver.F90: 104)
NVFORTRAN/x86-64 Linux 24.7-0: compilation aborted
```

As fixes were applied, the same class of ICE resurfaced at further locations in the same file (lines 113, 3148, 3150, and in the `CombineCRSMatIndices` subroutine), always with the same root cause.

### Root Cause

`nvhpc 24.7` has an internal compiler error in its IR-lowering pass triggered by:

1. **Whole-structure assignment** (`A = B`) where the derived type contains `DIMENSION(:), POINTER` components.
2. **Whole-array assignment** to `ALLOCATABLE` or `POINTER` array components where the array bounds are derived-type component access expressions (e.g., `SMat2 % NumberOfRows`) rather than plain scalar symbols.

The compiler backend requires a concrete, simple integer symbol to resolve array extents at lowering time. `gfortran` handles these cases correctly.

### Fix: Structure assignment (line ~104)

Replace the single structure assignment with explicit field-by-field copies of scalar components. The pointer array components (`Active`, `IsNeighbour`) are assigned on the immediately following lines anyway, so nothing is lost.

```fortran
! BEFORE
pes = ParEnv % PEs
ALLOCATE( SParMatrixDesc )
SParMatrixDesc % ParEnv = ParEnv

ALLOCATE(SParMatrixDesc % ParEnv % Active(ParEnv % PEs))
SParMatrixDesc % ParEnv % Active = ParEnv % Active

! AFTER
pes = ParEnv % PEs
ALLOCATE( SParMatrixDesc )
! nvhpc ICE workaround: structure assignment of a derived type with
! pointer array components triggers "Lowering Error" in the compiler backend.
! Copy scalar fields explicitly; Active and IsNeighbour are handled below.
SParMatrixDesc % ParEnv % PEs             = ParEnv % PEs
SParMatrixDesc % ParEnv % MyPE            = ParEnv % MyPE
SParMatrixDesc % ParEnv % Initialized     = ParEnv % Initialized
SParMatrixDesc % ParEnv % ActiveComm      = ParEnv % ActiveComm
SParMatrixDesc % ParEnv % NumOfNeighbours = ParEnv % NumOfNeighbours
SParMatrixDesc % ParEnv % NumberOfThreads = ParEnv % NumberOfThreads
SParMatrixDesc % ParEnv % ExternalInit    = ParEnv % ExternalInit

ALLOCATE(SParMatrixDesc % ParEnv % Active(ParEnv % PEs))
! nvhpc ICE workaround: whole-array assignment of a DIMENSION(:),POINTER
! component triggers "array upper bound is not a symbol" Lowering Error.
! SIZE() produces a concrete integer symbol the compiler backend can handle.
BLOCK
  INTEGER :: npes, i
  npes = SIZE(ParEnv % Active)
  DO i = 1, npes
    SParMatrixDesc % ParEnv % Active(i) = ParEnv % Active(i)
  END DO
END BLOCK
```

### Fix: Array initialization and `CombineCRSMatIndices` (lines ~3148 ff.)

Change local `POINTER` array variables to `ALLOCATABLE` (same semantics here, avoids one class of ICE), introduce a local scalar copy of derived-type component bounds and replace all whole-array assignments with explicit element-wise loops. 
The `BLOCK` construct is used to scope the temporary scalar so it does not pollute the outer scope.

```fortran
! BEFORE
INTEGER :: i, j, k, i1, i2, j1, j2, ind, ind1, ind2, DRows, DCols, row, col
INTEGER, POINTER :: cols(:)
LOGICAL, POINTER :: done(:)
...
ALLOCATE( Done( Smat2 % NumberOfRows ) )
done = .FALSE.
...
ALLOCATE( DMat % Rows(  SMat2 % NumberOfRows + 1) )
...
DMat % Rows = SMat2 % Rows(1:SMat2 % NumberOfRows+1)
DMat % GRows = SMat2 % GRows(1:SMat2 % NumberOfRows)
...

! AFTER
! nvhpcICE workaround: nrows2 is a plain local scalar used instead
! of Smat2 % NumberOfRows as an array bound — derived-type component accesses
! as bounds trigger "array numelm is not a symbol" Lowering Error.
INTEGER :: i, j, k, i1, i2, j1, j2, ind, ind1, ind2, DRows, DCols, row, col, nrows2
! nvhpc ICE workaround: ALLOCATE of a POINTER(:) with a non-trivial
! extent expression triggers "array numelm is not a symbol" Lowering Error.
! ALLOCATABLE has the same semantics here and avoids the compiler bug.
INTEGER, ALLOCATABLE :: cols(:)
LOGICAL, ALLOCATABLE :: done(:)
...
! nvhpc ICE workaround: any whole-array assignment to an ALLOCATABLE
! triggers "array numelm is not a symbol" Lowering Error. Use a scalar loop.
nrows2 = Smat2 % NumberOfRows
ALLOCATE( Done( nrows2 ) )
DO i = 1, nrows2
  Done(i) = .FALSE.
END DO
...
! nvhpc ICE workaround: whole-array assignments to derived-type
! POINTER array components trigger "array numelm" Lowering Error.
! Use local scalars and element-wise loops instead.
nrows2 = SMat2 % NumberOfRows
ALLOCATE( DMat % Rows(  nrows2 + 1) )
ALLOCATE( DMat % GRows( nrows2 ) )
ALLOCATE( DMat % RowOwner( nrows2 ) )
ALLOCATE( DMat % Cols( SMat2 % Rows(nrows2 + 1)-1 ) )
DMat % NumberOfRows = nrows2
DO i = 1, nrows2+1
  DMat % Rows(i) = SMat2 % Rows(i)
END DO
DO i = 1, nrows2
  DMat % GRows(i)    = SMat2 % GRows(i)
  DMat % RowOwner(i) = SMat2 % RowOwner(i)
END DO
DO i = 1, SIZE(DMat % Cols)
  DMat % Cols(i) = SMat2 % Cols(i)
END DO
```

The same pattern is applied symmetrically to the `SMat2 % NumberOfRows == 0` branch.


## Change 3: Procedure pointer component calls require explicit pass-object argument

**File changed:** `fem/src/ZirkaHysteresis.F90`

### Compiler Error

`nvfortran` rejected calls of the form:

```fortran
H = rc_p % simple_eval(B)
```

where `simple_eval` is a `PROCEDURE POINTER` component of a derived type. `nvfortran` correctly requires the object to be passed explicitly as the first argument because `NOPASS` is implied for procedure pointer components. `gfortran` silently passes the object implicitly, which is non-standard behavior.

An additional complication is that `gfortran 11` has an internal compiler error (ICE) when `PROCEDURE` pointer variables are declared inside a `BLOCK` construct within a subroutine that uses `CONTAINS`. To work around both issues simultaneously, procedure pointer temporaries are hoisted to the enclosing scoping unit.

### Fix

Introduce a local `PROCEDURE` pointer variable and extract the component pointer before calling, passing the object explicitly. For chained pointer accesses (`x % parent % parent % simple_eval`), introduce a `CLASS(*), POINTER` temporary first to avoid a second component-access ICE.

```fortran
! BEFORE
H = rc_p % simple_eval(B)

! AFTER
! Portable fix: extract procedure pointer and call directly (no % syntax)
! so neither gfortran nor nvhpc adds an implicit pass-object.
PROCEDURE(SimpleEvalRevCurve), POINTER :: eval_fn
...
eval_fn => rc_p % simple_eval
H = eval_fn(rc_p, B)
```

For chained accesses in `AddStack`:

```fortran
! BEFORE
Hpp = x % parent % parent % simple_eval(B)
Hp  = x % parent % simple_eval(B)

! AFTER
! Portable fix: extract procedure pointer and call directly (no % syntax)
! so neither gfortran nor nvhpc adds an implicit pass-object.
CLASS(RevCurve_t), POINTER :: tmp
PROCEDURE(SimpleEvalRevCurve), POINTER :: eval_fn
! (hoisted to subroutine scope to avoid GCC 11 ICE in lower_nested_functions
!  when PROCEDURE pointers appear inside BLOCK in a CONTAINS subroutine)
...
tmp => x % parent % parent
eval_fn => tmp % simple_eval
Hpp = eval_fn(tmp, B)
tmp => x % parent
eval_fn => tmp % simple_eval
Hp = eval_fn(tmp, B)
```

The same pattern is applied in `rc_printeval` wherever `% simple_eval(B)` appeared.

## Change 4: `VariableAdd` called with positional `NULL()` solver argument

**File changed:** `elmerice/Solvers/ValvingGeometry.F90`

### Compiler Error

`nvfortran` rejected the call because passing `NULL()` as a positional argument for a `TYPE(Solver_t), POINTER` dummy argument is ambiguous or non-conforming when the interface requires a typed pointer. Using keyword arguments makes the intended `OPTIONAL` / null-pointer usage unambiguous to all compilers.

### Fix

```fortran
! BEFORE
CALL VariableAdd(Mesh % Variables, Mesh, NULL(), "isoline id", 1, WorkReal, WorkPerm)

! AFTER
CALL VariableAdd(Mesh % Variables, Mesh, Name="isoline id", DOFs=1, Values=WorkReal, Perm=WorkPerm)
```

The solver argument is simply omitted (relying on its `OPTIONAL` or default-pointer status in the interface), which is cleaner and portable across all compilers.

## Change 5: `ISNAN` intrinsic not recognised by nvfortran

**Files changed:** `elmerice/Solvers/GlaDSCoupledSolver.F90`, `fem/src/modules/MagnetoDynamics2D.F90`, `fem/src/modules/CoordinateTransform.F90`

### Compiler Error

`nvfortran` does not provide `ISNAN` as a built-in intrinsic. `ISNAN` is a GNU extension available in `gfortran` but not guaranteed by the Fortran standard and not implemented in `nvfortran`. 
The standard-conforming replacement is `IEEE_IS_NAN` from the `IEEE_ARITHMETIC` intrinsic module (Fortran 2003+), which is supported by all modern compilers including `gfortran`, Intel, and `nvfortran`.

### Fix: `GlaDSCoupledSolver.F90`

Add the module `USE` at the top of the subroutine and replace the call:

```fortran
! BEFORE
IF(ISNAN(AreaSolution(k))) AreaSolution(k) = 0.0

! AFTER
USE IEEE_ARITHMETIC, ONLY: IEEE_IS_NAN
...
IF(IEEE_IS_NAN(AreaSolution(k))) AreaSolution(k) = 0.0
```

### Fix: `MagnetoDynamics2D.F90`

Add the module `USE` inside `BulkAssembly` and replace all occurrences:

```fortran
! BEFORE
IF (ISNAN(BodyLorentzForcesRe(i, j))) THEN
  BodyLorentzForcesRe(i, j) = 0._dp
END IF
IF (ISNAN(BodyLorentzForcesIm(i, j))) THEN
  BodyLorentzForcesIm(i, j) = 0._dp
END IF

! AFTER
USE IEEE_ARITHMETIC, ONLY: IEEE_IS_NAN
...
IF (IEEE_IS_NAN(BodyLorentzForcesRe(i, j))) THEN
  BodyLorentzForcesRe(i, j) = 0._dp
END IF
IF (IEEE_IS_NAN(BodyLorentzForcesIm(i, j))) THEN
  BodyLorentzForcesIm(i, j) = 0._dp
END IF
```

## Change 6: Missing explicit `USE` for `GetString` and `GaussPointsAdapt` in nested subroutines

**File changed:** `fem/src/modules/HeatSolveVec.F90`

### Compiler Error

`nvfortran` could not resolve `GetString` and `GaussPointsAdapt` inside the internal subroutines `LocalMatrixVec` and `LocalMatrix` because `nvfortran` is stricter than `gfortran` about host-association of `USE`-associated symbols into nested scoping units.
An explicit `USE` statement in each subroutine makes the dependency unambiguous.

### Fix

Add an explicit `USE` statement inside each affected internal subroutine:

```fortran
! In LocalMatrixVec:
! BEFORE
SUBROUTINE LocalMatrixVec( Element, n, nd, nb, VecAsm, InitHandles )
  USE LinearForms
  IMPLICIT NONE

! AFTER
SUBROUTINE LocalMatrixVec( Element, n, nd, nb, VecAsm, InitHandles )
  USE LinearForms
  USE DefUtils, ONLY: GetString, GaussPointsAdapt
  IMPLICIT NONE

! In LocalMatrix:
! BEFORE
SUBROUTINE LocalMatrix( Element, n, nd, nb, InitHandles )
  IMPLICIT NONE

! AFTER
SUBROUTINE LocalMatrix( Element, n, nd, nb, InitHandles )
  USE DefUtils, ONLY: GetString, GaussPointsAdapt
  IMPLICIT NONE
```

## Change 7: `SetMagneticFluxDensityBC` uses host-associated `Solver` and `Mesh` implicitly

**File changed:** `fem/src/modules/MagnetoDynamics2D.F90`

### Compiler Error

`nvfortran` rejected the use of `Solver` and `Mesh` inside `SetMagneticFluxDensityBC` because the subroutine is a `CONTAINS`-internal routine and `nvfortran` does not silently host-associate the outer `Solver` and `Mesh` objects in this context. 
Additionally, the fixed-size local arrays `Bx(Solver % Mesh % MaxElementDofs)` are not valid in `nvfortran` because the size is a derived-type component expression rather than a compile-time constant.

The fix makes the dependency explicit by adding `LSolver` and `LMesh` as dummy arguments, converting the fixed-size local arrays to `ALLOCATABLE` and updating all two call sites. The same fix is applied to both the real and harmonic variants of `SetMagneticFluxDensityBC` in the file.

### Fix

```fortran
! BEFORE
SUBROUTINE SetMagneticFluxDensityBC()
  IMPLICIT NONE
  ...
  REAL(KIND=dp) :: Bx(Solver % Mesh % MaxElementDofs), &
                   By(Solver % Mesh % MaxElementDofs)
  Perm => Solver % Variable % Perm
  A    => Solver % Matrix
  ...
  x = Mesh % Nodes % x(k)
  y = Mesh % Nodes % y(k)

! AFTER
SUBROUTINE SetMagneticFluxDensityBC(LSolver, LMesh)
  IMPLICIT NONE
  TYPE(Solver_t) :: LSolver
  TYPE(Mesh_t), POINTER :: LMesh
  ...
  REAL(KIND=dp), ALLOCATABLE :: Bx(:), By(:)
  ...
  ALLOCATE(Bx(LSolver % Mesh % MaxElementDofs), By(LSolver % Mesh % MaxElementDofs))
  Perm => LSolver % Variable % Perm
  A    => LSolver % Matrix
  ...
  x = LMesh % Nodes % x(k)
  y = LMesh % Nodes % y(k)
```

Call sites updated accordingly:

```fortran
! BEFORE
CALL SetMagneticFluxDensityBC()

! AFTER
CALL SetMagneticFluxDensityBC(Solver, Mesh)
```

The same change (with the additional `Bxim`/`Byim` arrays) is applied to the harmonic variant of `SetMagneticFluxDensityBC` further down in the same file.

## Note on Quad Precision (`HAVE_QP`): Build Flag vs. Source Change

**Relevant files (exemplarily):** `fem/src/ElementDescription.F90`, `fem/src/LinearAlgebra.F90` and others using `SELECTED_REAL_KIND(24)`

### Compiler Error

```
NVFORTRAN-S-0081-Illegal selector - KIND value must be non-negative
    (fem/src/LinearAlgebra.F90: 565)
NVFORTRAN-S-0081-Illegal selector - KIND value must be non-negative
    (fem/src/LinearAlgebra.F90: 566)
0 inform, 0 warnings, 2 severes, 0 fatal for invertmatrix3x3qp
```

### Root Cause

`nvfortran` does not support 128-bit quad precision floating point. `SELECTED_REAL_KIND(24)` returns `-5` (no such kind available) under `nvfortran`, which is then passed as a `KIND` selector, an illegal negative value that causes a hard compiler error.

### Two Approaches to Fix This

**Approach A: Source-level workaround (as done in branch `nvhpc-mods`):**

In each affected file, the quad precision kind is aliased to double precision (`dp`) directly:

```fortran
! BEFORE
INTEGER, PARAMETER :: qp = SELECTED_REAL_KIND(24)

! AFTER
!INTEGER, PARAMETER :: qp = SELECTED_REAL_KIND(24)   ! not supported by nvfortran
INTEGER, PARAMETER :: qp = dp                         ! fall back to double precision
```

**Approach B: CMake build flag (used in this fork):**

Quad precision support is disabled at the CMake configuration level:

```cmake
-DHAVE_QP=OFF
```

## Summary Table

| # | File | Nature of fix | Compiler(s) affected |
|---|------|---------------|----------------------|
| 1 | `fem/src/SParIterGlobals.F90` | Remove module-level `=>` init for derived type with pointer components | nvhpc |
| 1 | `fem/src/SParIterComm.F90` | Add explicit `IF (.NOT.ASSOCIATED) ParEnv => ParEnv_Common` guard | nvhpc |
| 2 | `fem/src/SParIterSolver.F90` | Replace structure assignment and whole-array ops with explicit scalar/element loops to avoid ICE in lowering pass | nvhpc 24.7 |
| 3 | `fem/src/ZirkaHysteresis.F90` | Extract procedure pointer components before calling; pass object explicitly | nvhpc, gfortran 11 |
| 4 | `elmerice/Solvers/ValvingGeometry.F90` | Use keyword arguments to `VariableAdd` instead of positional `NULL()` | nvhpc |
| 5 | `elmerice/Solvers/GlaDSCoupledSolver.F90` | Replace `ISNAN` (GNU extension) with `IEEE_IS_NAN` | nvhpc |
| 5 | `fem/src/modules/MagnetoDynamics2D.F90` | Replace `ISNAN` with `IEEE_IS_NAN`; add `USE IEEE_ARITHMETIC` | nvhpc |
| 6 | `fem/src/modules/HeatSolveVec.F90` | Add explicit `USE DefUtils` in nested subroutines | nvhpc |
| 7 | `fem/src/modules/MagnetoDynamics2D.F90` | Pass `Solver`/`Mesh` explicitly to `SetMagneticFluxDensityBC`; convert fixed-size arrays to `ALLOCATABLE` | nvhpc |
