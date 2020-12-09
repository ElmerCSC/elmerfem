      FUNCTION Random(Model,nodenumber,x) RESULT(r)
       USE DefUtils
       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp),INTENT(IN) :: x(2)
       REAL(kind=dp) :: r

       INTEGER :: ssize
       INTEGER,ALLOCATABLE :: seed(:)
       LOGICAL, SAVE :: Firsttime=.TRUE.
       LOGICAL :: Found

       CHARACTER(LEN=MAX_NAME_LEN),SAVE :: Rtype

       IF (Firsttime) THEN
         CALL random_seed(size=ssize)
         allocate(seed(ssize))
         
         seed = ListGetInteger( Model % Constants, 'Random Seed',Found )
         IF (Found)  call random_seed( put=seed )
         CALL random_seed(get=seed)

         Rtype = ListGetString( Model % Constants, 'Random Function',UnFoundFatal=.TRUE.)
         deallocate(seed)
         Firsttime=.FALSE.
       ENDIF


       SELECT CASE(Rtype)
       CASE('even')
         r = EvenRandom()
       CASE('normal')
         r = NormalRandom()
       CASE DEFAULT
         CALL FATAL('Random Function','Random Function should be <even> or <normal>')
       END SELECT

       r = r*x(1)+x(2)

       End FUNCTION Random


