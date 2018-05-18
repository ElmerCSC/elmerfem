      FUNCTION DistanceCond(Model,nodenumber,VarIn) RESULT(VarOut)
      USE DefUtils
      implicit none
      !-----------------
      TYPE(Model_t) :: Model
      INTEGER :: nodenumber
      REAL(kind=dp) :: VarIn,VarOut

      IF (VarIn.LT.0.1) THEN
        VarOut = +1.0
      ELSE
        VarOut = -1.0
      END IF
      End FUNCTION DistanceCond

      FUNCTION HMax(Model,nodenumber,VarIn) RESULT(VarOut)
      USE DefUtils
      implicit none
      !-----------------
      TYPE(Model_t) :: Model
      INTEGER :: nodenumber
      REAL(kind=dp) :: VarIn !Distance
      REAL(kind=dp) :: VarOut

      REAL(kind=dp),SAVE :: Extent,HMaxIN,HMaxOut
      LOGICAL,SAVE :: FirstVisit=.TRUE.

      IF (FirstVisit) THEN
        Extent=ListGetConstReal( Model % Constants,'Hmax margin extent',UnFoundFatal=.TRUE.)
        HMaxIN=ListGetConstReal( Model % Constants,'Hmax within margin',UnFoundFatal=.TRUE.)
       HMaxOut=ListGetConstReal( Model % Constants,'Hmax outside margin',UnFoundFatal=.TRUE.)
      END IF

      IF ((VarIn.GT.1.0).AND.(VarIn.LT.Extent)) THEN
         VarOut = HMaxIN
      ELSE
         VarOut = HMaxOut
      END IF

      End FUNCTION HMax
