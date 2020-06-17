!#
!# A User function for Ronne-Filchner test case;
!# to impose Conditionnal Dirichlet conditions only for Grounded ice and GL
!#   VarIn=GroundedMask (-1 : for floating ice; 0 : for GL; 1: for grounded ice)
!#   Varout return a positive value for Grounded ice and GL
       FUNCTION GM_CONDITION(Model,nodenumber,VarIn) RESULT(VarOut)
       USE DefUtils
       implicit none
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: VarIn,VarOut

        VarOut=VarIn+0.05
       END FUNCTION GM_CONDITION

!# A User function for Ronne-Filchner test case;
!# to impose passive condition for regulrisation of beta in floating
!elements
!#   VarIn=GroundedMask (-1 : for floating ice; 0 : for GL; 1: for grounded ice)
!#   Varout return a positive value for floating ice
     FUNCTION BetaRegPassive(Model,nodenumber,VarIn) RESULT(VarOut)
       USE DefUtils
       implicit none
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: VarIn,VarOut

        VarOut=-VarIn-0.5_dp
       END FUNCTION BetaRegPassive
