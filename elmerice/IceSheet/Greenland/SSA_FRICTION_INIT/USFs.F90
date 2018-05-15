!########################################################################      
!### A user function to define elements as passive if VarIn(=H)<Hmin
!#######################################################################
      FUNCTION PassiveCond_H(Model,nodenumber,VarIn)  RESULT(VarOut)
       USE Types
       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: VarIn
       REAL(kind=dp) :: VarOut
       !-----------------
       REAL(kind=dp) :: H
       REAL(kind=dp),parameter :: Hmin=1.0

       H=VarIn
       IF (H.LT.Hmin) THEN
          VarOut=1.0_dp
       ELSE
          VarOut=-1.0_dp
       END IF

      END FUNCTION PassiveCond_H
