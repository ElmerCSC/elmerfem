# User Function Budd Sliding Relation
## General Information
- **USF Fortran File:** USF_Sliding.f90
- **USF Name:** Sliding_Budd
- **Required Input Variable(s):** Vertical coord, z

## General Description
This is a summary description. See comments in the code module for a more complete description.

Basal shear stress is given by:

*tau_b = C.{u_b}^{m - 1} . u_b . Zab^{q}*
where *Zab* is height above buoyancy, *C* is a friction coefficient, and *m* and *q* are exponents.

The parameters to be given (in the bed boundary condition) are:
- Budd Friction Coefficient → *C*
- Budd Velocity Exponent → *m*
- Budd Zab Exponent → *q*
- Budd Linear Velocity → *u_{t0}*
- Budd Floatation
If Budd Floatation is set to true then the height above buoyancy will be based on the floatation condition instead of inferred from the effective pressure (i.e. depth is used instead of normal stress). Default is false.

## References
Bill Budd's original work:

Budd, W., Keage, P. L., and Blundy, N. A.: Empirical studies of ice sliding, J. Glaciol., 23, 157–170, 1979.

Budd, W., Jenssen, D., and Smith, I.: A 3-dimensional time- dependent model of the West Antarctic Ice-Sheet, Ann. Glaciol., 5, 29–36, 1984.

Recent implementation in Elmer/Ice:

Gladstone, R.M., R.C. Warner, B.K. Galton-Fenzi, O. Gagliardini, T. Zwinger and R. Greve, 2017. Marine ice sheet model performance depends on basal sliding physics and sub-shelf melting, The Cryosphere, 11, 319-329, doi:10.5194/tc-11-319-2017 pdf
