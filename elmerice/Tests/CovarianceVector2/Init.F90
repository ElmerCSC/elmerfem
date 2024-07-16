! *****************************************************************************
! Initialize an impulse; i.e.:
! r=1 if distance to prescribed "center" i < 100 APES
! r=0 otherwise
! There should be a node at the center!!
!#############################################################################
      FUNCTION Init(Model,nodenumber,xy) RESULT(r)
       USE DefUtils
       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp),INTENT(IN) :: xy(4)
       REAL(kind=dp) :: r
       REAL(kind=dp) :: d


       d=sqrt((xy(1)-xy(3))**2+(xy(2)-xy(4))**2)
       If (d.LT.100*AEPS) THEN
               r=1._dp
       Else
               r=0._dp
       End if

       End FUNCTION Init


