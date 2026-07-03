Test for SubShelfMeltLinearTF solver
====================================

This test now checks branch-sensitive behaviour, not only the base melt formula.

Base formula:

  m0 = gammaT * (rhow * cp / (rhoi * Lf)) * (Tw - Tf)  [m/yr]

where:

  Tf = lambda1*S + lambda2 + cc*z

Common setup:
- 2D mesh (1000 m x 500 m)
- Uniform ocean temperature: temp_oce_post = -1.0 degC
- Uniform ocean salinity: sal_oce_post = 34.0 PSU
- Ice base elevation: Zb = -500 m

Scenario 1 (floating + water column scaling):
- GroundedMask = -1
- water column scaling = True
- bedrock = -510 m, scaling factor = 75 m
- Expected:
  Tf      = -2.2455 degC
  Tw - Tf = 1.2455 degC
  m0      = 52.19234573 m/yr
  wct     = Zb - bedrock = 10 m
  scaling = tanh(10 / (75/e)) = 0.3473593091
  m       = 18.12949715 m/yr

Scenario 2 (grounding-line masking):
- GroundedMask = 0 everywhere
- grounding line melt = False
- Expected melt = 0 everywhere

Scenario 3 (grounded masking):
- GroundedMask = 1 everywhere
- grounding line melt = True
- Expected melt = 0 everywhere

Limitations (important)
-----------------------

This test is not exhaustive.

- It validates solver norms only. It does not compare full nodal fields, so some spatially local mistakes can remain undetected if they do not change the norm enough.
- It does not exercise the branch `IF (Solver % Variable % Perm(ii) .LE. 0) CYCLE` because this test runs in a simple serial configuration where node permutations are active.
- It does not directly validate the element flux output `<variable>_flux`; the element integration block is not checked by this test.
- The `water column scaling = False` path is currently only used in fully masked scenarios (expected zero melt). This means a bug in the non-masked `meltScaling = 1.0` assignment can still escape detection in this specific test.
- The setup uses uniform ocean forcing and mask values per scenario; mixed masks and spatially varying forcing are not covered.
