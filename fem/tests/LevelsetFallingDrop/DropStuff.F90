
  FUNCTION InitCircle( Model, n, t ) RESULT(f)
    USE Types
    USE SolverUtils
    IMPLICIT NONE
    
    TYPE(Model_t) :: Model
    LOGICAL :: Visited = .FALSE., GotIt
    INTEGER :: n,i,j
    REAL(KIND=dp) :: t,x,y,r,Lx,Ly,f,a1,b1,r1
    
    SAVE Visited, a1, b1, r1  

    IF(.NOT. Visited) THEN
      Lx = 1.0
      Ly = 1.0
      
      a1 = 0.5 * Lx 
      b1 = 0.5 * Ly
      r1 = 0.25 * Lx
      Visited = .TRUE.
    END IF
    
    x = Model % Nodes % x(n)
    y = Model % Nodes % y(n)    
    

    ! positive inside the circle
    f = r1 - SQRT((x-a1)**2.0 + (y-b1)**2.0) 

  END FUNCTION InitCircle




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



 FUNCTION FrameVelo( Model, n, dt ) RESULT(v)
   USE Types
   USE SolverUtils

   IMPLICIT NONE

   TYPE(Model_t) :: Model
   INTEGER :: i,j,k,n,DoneTime

   REAL(KIND=dp) :: x,y,dy,dt,v,v1,relax,dist=0.0d0
   LOGICAL :: GotIt, Visited 

   SAVE v1, dist, DoneTime

   IF( Model % Solver % DoneTime /= DoneTime) THEN
     relax = ListGetConstReal( Model % Simulation,'Frame Velocity Relaxation Factor')
     DoneTime = Model % Solver % DoneTime

     dy = ListGetConstReal(Model % Simulation,"res: levelset center 2",GotIt) 
     v1 = relax * (0.5-dy)/dt

     dist = dist + v1 * dy
     CALL ListAddConstReal(Model % Simulation,"res: cumulative distance",dist) 

     Visited = .TRUE.
   END IF

   v = v1

 END FUNCTION FrameVelo
