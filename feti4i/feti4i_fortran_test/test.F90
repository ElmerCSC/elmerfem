PROGRAM test

  USE FETI4I

  TYPE(C_PTR) :: instance
  TYPE(C_PTR) :: matrix

  INTEGER(KIND=FETI4IInt) :: n = 4
  INTEGER(KIND=FETI4IInt) :: nelems = 5
  INTEGER(KIND=FETI4IInt) :: size = 0 !TODO when this is known
  INTEGER(KIND=FETI4IInt), ALLOCATABLE :: inds(:)
  REAL(KIND=FETI4IReal), ALLOCATABLE :: vals(:)
  REAL(KIND=FETI4IReal), ALLOCATABLE :: b(:)
  INTEGER(KIND=FETI4IInt), ALLOCATABLE :: l2g(:)
  INTEGER(KIND=FETI4IMPIInt) :: neighbours_size = 3
  INTEGER(KIND=FETI4IMPIInt), ALLOCATABLE :: neighbours(:)
  INTEGER(KIND=FETI4IInt) :: dirichlet_size = 2
  INTEGER(KIND=FETI4IInt), ALLOCATABLE :: dirichlet_indices(:)
  REAL(KIND=FETI4IReal), ALLOCATABLE :: dirichlet_values(:)
  REAL(KIND=FETI4IReal), ALLOCATABLE :: x(:)

  INTEGER :: ALLOC_ERR

  CALL FETI4ICreateStiffnessMatrix(matrix, 1)

  ALLOCATE(vals(n*n), inds(n))
  DO i=1,n
    DO j=1,n
     vals((i-1)*n+j) = (i-1)*n+j
    END DO
    inds(i) = i
  END DO

  size = 1
  DO i=1,nelems
    CALL FETI4IAddElement(matrix, n, inds, vals)
    size = size + (n-1)
  END DO
  
  ALLOCATE(b(size), l2g(size))
  DO i=1,size
    b(i) = -i
    l2g(i) = 100+i
  END DO

  write (*,*) "b: ", b

  ALLOCATE(neighbours(neighbours_size))
  DO i=1,neighbours_size
    neighbours(i) = i-1
  END DO

  ALLOCATE(dirichlet_indices(dirichlet_size), dirichlet_values(dirichlet_size))
  DO i=1,dirichlet_size
    dirichlet_indices(i) = 2*(i+1)
    dirichlet_values(i)  = -2*(i+1)
  END DO

  CALL FETI4ICreateInstance(instance, matrix, size, b, l2g, &
    neighbours_size, neighbours, &
    dirichlet_size, dirichlet_indices, dirichlet_values)

  ALLOCATE(x(size), STAT = ALLOC_ERR)
  write (*,*) "ALLOC_ERR=",ALLOC_ERR

  x(:) = 1.0

  write (*,*) "x: ", x

  CALL FETI4ISolve(instance, size, x)

  CALL FETI4IDestroy(instance)

END PROGRAM test 
