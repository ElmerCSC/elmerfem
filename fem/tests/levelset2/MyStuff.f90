
  FUNCTION Viscosity( Model, n, fii ) RESULT(r)
   USE Types
   USE SolverUtils

   IMPLICIT NONE

   TYPE(Model_t) :: Model
   INTEGER :: i,j,k,n, body_id

   REAL(KIND=dp) :: x,y,r,fii,dfii,t,visc0,visc1,fstep
   LOGICAL :: Visited = .FALSE., GotIt


   SAVE Visited, visc1, visc0, dfii

   IF(.NOT. Visited) THEN
     body_id = Model % CurrentElement % Bodyid    
     k = ListGetInteger( Model % Bodies( body_id ) % Values, 'Material' )
     visc0 = ListGetConstReal(Model % Materials(k) % Values,'Inside Viscosity')
     visc1 = ListGetConstReal(Model % Materials(k) % Values,'Outside Viscosity')
     dfii = ListGetConstReal(Model % Materials(k) % Values,'Levelset Bandwidth')
     Visited = .TRUE.
   END IF

   t = fii / dfii
   IF(t < -1.0) THEN
     fstep = 0.0d0
   ELSE IF(t > 1.0) THEN
     fstep = 1.0d0
   ELSE
     fstep = 0.75d0 * (t - t**3.0d0/3) + 0.5d0
   END IF

   r = (visc0-visc1) * fstep + visc1

 END FUNCTION Viscosity



  FUNCTION Density( Model, n, fii ) RESULT(r)
   USE Types
   USE SolverUtils

   IMPLICIT NONE

   TYPE(Model_t) :: Model
   INTEGER :: i,j,k,n, body_id

   REAL(KIND=dp) :: x,y,r,fii,dfii,t,dens0,dens1,fstep
   LOGICAL :: Visited = .FALSE., GotIt


   SAVE Visited, dens1, dens0, dfii

   IF(.NOT. Visited) THEN
     body_id = Model % CurrentElement % Bodyid    
     k = ListGetInteger( Model % Bodies( body_id ) % Values, 'Material' )
     dens0 = ListGetConstReal(Model % Materials(k) % Values,'Inside Density')
     dens1 = ListGetConstReal(Model % Materials(k) % Values,'Outside Density')
     dfii = ListGetConstReal(Model % Materials(k) % Values,'Levelset Bandwidth')
     Visited = .TRUE.
   END IF

   t = fii / dfii
   IF(t < -1.0) THEN
     fstep = 0.0d0
   ELSE IF(t > 1.0) THEN
     fstep = 1.0d0
   ELSE
     fstep = 0.75d0 * (t - t**3.0d0/3) + 0.5d0
   END IF

   r = (dens0-dens1) * fstep + dens1

 END FUNCTION Density

