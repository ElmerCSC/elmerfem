This is a simple 2-D test case where the UMAT subroutine is called to obtain 
the definition of the constitutive law of St Venant-Kirchhoff type (the case of 
plane strain). The case has been checked to yield the same result as obtained 
by using the earlier implementation of the same model within this solver. Note 
that the nonlinear iterations of the two versions may however perform 
differently as the umat version employs an inexact Newton iteration.
