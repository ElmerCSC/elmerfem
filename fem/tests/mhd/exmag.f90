! The additional terms in the induction equation
! due to a time harmonic vertical field B_ac = B0 sin(wt) e_y

FUNCTION appliedB( model,n,t ) RESULT(c1)
   USE types

   TYPE(model_t) :: model
   INTEGER :: n
   REAL(KIND=dp) :: t,c1
   REAL(KIND=dp), PARAMETER :: B0=0.0001_dp, f=3_dp

! B_ac for nabla x ( v x B )
   c1 = B0 * SIN(2.0_dp*PI*f*t)

 END FUNCTION appliedB

FUNCTION applieddBdt( model,n,t ) RESULT(c2)
   USE types

   TYPE(model_t) :: model
   INTEGER :: n
   REAL(KIND=dp) :: t,c2,w
   REAL(KIND=dp), PARAMETER :: B0=0.0001_dp, f=3_dp

   w = 2.0_dp * PI * f
! -@B_ac/@t
   c2 = - B0 * w * COS(w*t)

 END FUNCTION applieddBdt
