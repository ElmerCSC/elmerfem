!/******************************************************************************
! *
! *  Authors: Peter R�back & Juha Ruokolainen
! *  Email:   Peter.Raback@csc.fi & Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 26.5. 2010
! *
! *****************************************************************************/


!-------------------------------------------------------------------------
!> This subroutine provides a shared library handle to the internal subroutine
!> in order to enable more flexible use of the routine.
!> \ingroup Solvers
!-------------------------------------------------------------------------
SUBROUTINE ParticleOutputSolver( Model,Solver,dt,TransientSimulation )
  USE DefUtils
  USE ParticleUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation

  CALL Info('ParticleOutputSolver','-----------------------------------------')
  CALL Info('ParticleOutputSolver','Saving particle data')

  CALL SaveParticleData( Model,Solver,dt,TransientSimulation )

  CALL Info('ParticleOutputSolver','All done')
  CALL Info('ParticleOutputSolver','-----------------------------------------')

END SUBROUTINE ParticleOutputSolver


!------------------------------------------------------------------------------
!>  This subroutine is used to save results in a uniform grid.
!>  This is achieved by first creating particles at the positions,
!>  and then saving them with more generic routines.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE SaveGridData( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  USE ParticleUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver  !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model    !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt       !< Timestep size for time dependent simulations
  LOGICAL :: TransientSimulation      !< Steady state or transient simulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------

  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: Params
  LOGICAL :: GotIt, Structured, Visited = .FALSE.
  CHARACTER(LEN=MAX_NAME_LEN) :: FileFormat
  LOGICAL :: Found,VtiFormat, VtuFormat, TableFormat, AnyFormat, &
      RecreateGrid
  INTEGER :: Extent(6)
  REAL(KIND=dp) :: Origin(3), Dx(3)
  INTEGER, POINTER :: GridIndex(:,:,:)
  TYPE( Particle_t), POINTER :: Particles  

  SAVE Visited, GridIndex, Origin, Dx, Extent

!------------------------------------------------------------------------------

  CALL Info('SaveGridData','-----------------------------------------', Level=4 )
  CALL Info('SaveGridData','Saving data an uniform grid point        ', Level=4 )

  Particles => GlobalParticles
  Params => GetSolverParams()
  Mesh => Solver % Mesh

  TableFormat = ListGetLogical( Params,'Table Format',Found)
  VtuFormat = ListGetLogical( Params,'Vtu Format',Found)
  VtiFormat = GetLogical( Params,'Vti Format',GotIt )

  FileFormat = ListGetString( Params,'Output Format',Found) 
  IF( Found ) THEN
    IF( FileFormat == 'vtu') VtuFormat = .TRUE.
    IF( FileFormat == 'table') TableFormat = .TRUE.
    IF( FileFormat == 'vti') VtiFormat = .TRUE.
  END IF

  AnyFormat = VtuFormat .OR. TableFormat .OR. VtiFormat 
  IF( .NOT. AnyFormat ) THEN
    CALL Warn('SaveGridData','No active file format given!')
    RETURN
  END IF

  RecreateGrid = ListGetLogical( Params,'Recreate Grid', Found ) 

  ! Initialize the particles on the first calling
  !------------------------------------------------------------------------
  Structured = VtiFormat
  IF( .NOT. Visited .OR. RecreateGrid ) THEN
    Particles % TimeOrder = 0
    Particles % dim = CoordinateSystemDimension()
    CALL CreateGridParticles( Particles ) 
    CALL CreateListForSaving( Model, Solver % Values,.TRUE. )    
    Visited = .TRUE.
  END IF
  
  ! The calling is split to two since the 1st one requires totally different 
  ! data structure. 
  !--------------------------------------------------------------------------
  IF( VtiFormat ) CALL ParticleOutputVti( Particles, Extent, Origin, Dx, GridIndex )
  IF( TableFormat ) CALL ParticleOutputTable( Particles )
  IF( VtuFormat ) CALL ParticleOutputVtu( Particles ) 

  IF( RecreateGrid ) CALL DestroyParticles( Particles ) 


  CALL Info('SaveGridData','All done',Level=4)
  CALL Info('SaveGridData', '-----------------------------------------', Level=4 )
  

CONTAINS


  !---------------------------------------------------------
  ! Allocate particles for visualization
  !---------------------------------------------------------
  SUBROUTINE AllocateGridParticles(Particles,NoParticles)
    
    TYPE(Particle_t), POINTER :: Particles
    INTEGER :: NoParticles   
    REAL(KIND=dp), POINTER :: Coordinate(:,:), UVW(:,:)
    INTEGER, POINTER :: Status(:), ElementIndex(:)
    INTEGER :: PrevNoParticles, dim, No, n
    

    IF( NoParticles <= Particles % MaxNumberOfParticles ) THEN
      CALL Info('AllocateParticles','There are already enough particles',Level=12)
      RETURN
    ELSE
      CALL Info('AllocateParticles','Allocating number of particles: '// &
          TRIM(I2S(NoParticles)),Level=12)    
    END IF
    
    dim = Particles % dim 
    
    IF( Particles % MaxNumberOfParticles == 0 ) THEN
      ALLOCATE( Particles % Coordinate(NoParticles,dim))
      ALLOCATE( Particles % uvw(NoParticles,dim))
      ALLOCATE( Particles % ElementIndex(NoParticles))
      ALLOCATE( Particles % Status(NoParticles))
      
      Particles % Coordinate = 0.0_dp
      Particles % uvw = 0.0_dp
      Particles % ElementIndex = 0
      Particles % Status = PARTICLE_ALLOCATED

      Particles % MaxNumberOfParticles = NoParticles
    ELSE
      Coordinate => Particles % Coordinate
      UVW => Particles % UVW
      Status => Particles % Status
      ElementIndex => Particles % ElementIndex
      
      ALLOCATE( Particles % Coordinate(NoParticles,dim) )
      ALLOCATE( Particles % UVW(NoParticles,dim) )
      ALLOCATE( Particles % Status(NoParticles) )
      ALLOCATE( Particles % ElementIndex(NoParticles) )
      
      ! ------------------------

      PrevNoParticles = Particles % MaxNumberOfParticles
      Particles % NumberOfParticles = NoParticles
      

      Particles % Coordinate(1:PrevNoParticles,:) = Coordinate
      Particles % UVW(1:PrevNoParticles,:) = UVW
      Particles % ElementIndex(1:PrevNoParticles) = ElementIndex
      Particles % Status(1:PrevNoParticles) = Status


      DEALLOCATE(Coordinate, uvw, Status, ElementIndex)
      
      Particles % Coordinate(PrevNoParticles+1:NoParticles,:) = 0.0_dp
      Particles % UVW(PrevNoParticles+1:NoParticles,:) = 0.0_dp
      Particles % ElementIndex(PrevNoParticles+1:NoParticles) = 0
      Particles % Status(PrevNoParticles+1:NoParticles) = PARTICLE_ALLOCATED
      
      Particles % MaxNumberOfParticles = NoParticles
    END IF
  END SUBROUTINE AllocateGridParticles
  

  !-----------------------------------------------------------------
  ! Add one single particle to list, allocating more space if needed
  !------------------------------------------------------------------

  SUBROUTINE AddGridParticle(Particles,ElementIndex,GlobalCoords,LocalCoords)

    TYPE(Particle_t), POINTER :: Particles
    INTEGER :: ElementIndex
    REAL(KIND=dp) :: GlobalCoords(3),LocalCoords(3)
    
    INTEGER :: i,j,n,m,dim

    n = Particles % NumberOfParticles 
    IF( n == Particles % MaxNumberOfParticles ) THEN
      m = MAX( 1000, n / 2 ) 
      CALL AllocateGridParticles( Particles, n + m )
    END IF
    
    dim = Particles % dim
    n = n + 1
    
    Particles % NumberOfParticles = n
    Particles % ElementIndex(n) = ElementIndex
    Particles % Status(n) = PARTICLE_READY
    
    Particles % Coordinate(n,1:dim) = GlobalCoords(1:dim)
    Particles % uvw(n,1:dim) = LocalCoords(1:dim)

    
  END SUBROUTINE AddGridParticle


  !-----------------------------------------------------------------
  ! Find grid particles in a uniform grid
  !------------------------------------------------------------------     
  SUBROUTINE CreateGridParticles(Particles)
     
    TYPE(Particle_t), POINTER :: Particles

    TYPE(Mesh_t), POINTER :: Mesh
    REAL(KIND=dp) :: MinCoord(3), MaxCoord(3)
    REAL(KIND=dp) :: LocalCoords(3), GlobalCoords(3)
    REAL(KIND=dp) :: x,y,z,u,v,w
    INTEGER, POINTER :: MaskPerm(:)
    INTEGER :: t,i,j,k,n,nx,ny,nz,imin,imax,jmin,jmax,kmin,kmax, &
        imintot,imaxtot,jmintot,jmaxtot,kmintot,kmaxtot,&
        meshdim, griddim, ActiveCoordinate,ntot
    INTEGER :: ioff,joff,koff,cands1, cands2, ierr, totcount(3),tmpcount(3)
    INTEGER :: ParallelNodes, NumberOfNodes
    TYPE(Nodes_t) :: Nodes
    TYPE(ValueList_t), POINTER :: Params
    TYPE(ELement_t), POINTER :: Element
    LOGICAL, POINTER :: GridPointActive(:,:,:)
    LOGICAL :: CheckForDuplicates, LowerDimensional, MaskExist, MaskOnBulk, GotIt
    INTEGER :: ElemStart, ElemFin
    CHARACTER(MAX_NAME_LEN) :: Str


    CALL Info('SaveGridData','Saving data on uniform grid',Level=4)

    Mesh => GetMesh()
    IF( .NOT. ASSOCIATED( Mesh ) ) THEN
      CALL Fatal('FindGridParticles','No Mesh associated')
    END IF

    IF( Visited .AND. Structured) THEN
      DEALLOCATE( GridIndex )
    END IF


    meshdim = Mesh % MeshDim
    Params => GetSolverParams()

    ! Create a mask for saving only part of data
    !---------------------------------------------------------------
    MaskExist = .FALSE.
    LowerDimensional = .FALSE.
    Str = ListGetString( Params,'Mask Name',GotIt) 
    IF( GotIt ) THEN
      ALLOCATE( MaskPerm( Model % NumberOfNodes ) ) 
      CALL MakePermUsingMask( Model,Solver,Mesh,Str, &
          .FALSE., MaskPerm, NumberOfNodes, MaskOnBulk )
      ParallelNodes = NINT( ParallelReduction( 1.0_dp * NumberOfNodes ) )
      IF( ParallelNodes == 0 ) THEN
        CALL Fatal('SaveGridData','Given mask not active: '//TRIM(Str) )
      ELSE
        MaskExist = .TRUE.
        CALL Info('SaveGridData','Using > '// TRIM(Str) // ' < as mask variable',Level=5)
        LowerDimensional = .NOT. MaskOnBulk
      END IF
    END IF
    
    IF( LowerDimensional ) THEN
      griddim = meshdim - 1
      ElemStart = Mesh % NumberOfBulkElements + 1
      ElemFin = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
      CALL Info('SaveGridData','Assuming reduced coordinate to be: '//TRIM(I2S(meshdim)),Level=5)
    ELSE
      griddim = meshdim
      ElemStart = 1
      ElemFin = Mesh % NumberOfBulkElements 
      ActiveCoordinate = 0
    END IF
    CALL Info('SaveGridData','Saving data on '//TRIM(I2S(griddim))//'D grid',Level=5)


    ! The bounding box may be given, otherwise it is taken to include the whole mesh
    !-------------------------------------------------------------------------------
    MinCoord(1) = GetCReal( Params,'Min Coordinate 1',GotIt) 
    IF(.NOT. GotIt) MinCoord(1) = MINVAL(Mesh % Nodes % x )
    MinCoord(2) = GetCReal( Params,'Min Coordinate 2',GotIt) 
    IF(.NOT. GotIt) MinCoord(2) = MINVAL(Mesh % Nodes % y )
    MinCoord(3) = GetCReal( Params,'Min Coordinate 3',GotIt) 
    IF(.NOT. GotIt) MinCoord(3) = MINVAL(Mesh % Nodes % z )

    MaxCoord(1) = GetCReal( Params,'Max Coordinate 1',GotIt) 
    IF(.NOT. GotIt) MaxCoord(1) = MAXVAL(Mesh % Nodes % x )
    MaxCoord(2) = GetCReal( Params,'Max Coordinate 2',GotIt) 
    IF(.NOT. GotIt) MaxCoord(2) = MAXVAL(Mesh % Nodes % y )    
    MaxCoord(3) = GetCReal( Params,'Max Coordinate 3',GotIt) 
    IF(.NOT. GotIt) MaxCoord(3) = MAXVAL(Mesh % Nodes % z )

    ! print *,'Bounding box min:',MinCoord
    ! print *,'Bounding box max:',MaxCoord

    ! Optionally the mesh origin may be moved to guarantee that there is 
    ! a node at (x0,y0,z0) always.
    !--------------------------------------------------------------------
    IF( GetLogical( Params,'Grid Origin At Corner',GotIt ) ) THEN
      Origin(1:3) = MinCoord(1:3)
    ELSE
      Origin(1) = GetCReal( Params,'Grid Origin 1',GotIt) 
      Origin(2) = GetCReal( Params,'Grid Origin 2',GotIt) 
      Origin(3) = GetCReal( Params,'Grid Origin 3',GotIt) 
    END IF
    ! print *,'Origin:',Origin


    ! Get the grid resolution assuming that the grid is cartesian and uniform.
    !----------------------------------------------------------------------------
    dx(1) = GetCReal( Params,'Grid dx',GotIt) 
    IF(.NOT. GotIt ) THEN
      nx = GetInteger( Params,'Grid nx',GotIt) 
      IF( GotIt) THEN
        dx(1) = ( MaxCoord(1) - MinCoord(1) ) / nx 
      ELSE
        CALL Fatal('FindGridParticles','Give either > Grid dx < or > Grid nx <')
      END IF
    END IF

    IF( griddim >= 2 ) THEN
      dx(2) = GetCReal( Params,'Grid dy',GotIt) 
      IF(.NOT. GotIt ) THEN
        nx = GetInteger( Params,'Grid ny',GotIt) 
        IF( GotIt) THEN
          dx(2) = ( MaxCoord(2) - MinCoord(2) ) / nx
        ELSE
          dx(2) = dx(1)
        END IF
      END IF
    END IF
      
    IF( griddim == 3 ) THEN
      dx(3) = GetCReal( Params,'Grid dz',GotIt) 
      IF(.NOT. GotIt ) THEN
        nx = GetInteger( Params,'Grid nz',GotIt) 
        IF( GotIt) THEN
          dx(3) = ( MaxCoord(3) - MinCoord(3) ) / nx
        ELSE
          dx(3) = dx(1)
        END IF
      END IF
    END IF


    ! Set limits for the global indexes. These are used particularly if the 
    ! bounding box has been manually reduced. 
    !----------------------------------------------------------------------------
    imintot = CEILING( ( MinCoord(1) - Origin(1) ) / dx(1) ) 
    imaxtot = FLOOR( ( MaxCoord(1) - Origin(1) ) / dx(1) ) 
    
    IF( griddim >= 2 ) THEN
      jmintot = CEILING( ( MinCoord(2) - Origin(2) ) / dx(2) ) 
      jmaxtot = FLOOR( ( MaxCoord(2) - Origin(2) ) / dx(2) ) 
    ELSE
      jmintot = 0
      jmaxtot = 0
    END IF

    IF( griddim == 3 ) THEN 
      kmintot = CEILING( ( MinCoord(3) - Origin(3) ) / dx(3) ) 
      kmaxtot = FLOOR( ( MaxCoord(3) - Origin(3) ) / dx(3) ) 
    ELSE
      kmintot = 0
      kmaxtot = 0
    END IF

    CALL Info('SaveGridData','Index i range: '&
        //TRIM(I2S(imintot))//' - '//TRIM(I2S(imaxtot)),Level=12)
    CALL Info('SaveGridData','Index j range: '&
        //TRIM(I2S(jmintot))//' - '//TRIM(I2S(jmaxtot)),Level=12)
    CALL Info('SaveGridData','Index k range: '&
        //TRIM(I2S(kmintot))//' - '//TRIM(I2S(kmaxtot)),Level=12)

    ioff = imintot-1
    joff = jmintot-1
    koff = kmintot-1
    
    ! Create a table for checking active gridpoints
    !----------------------------------------------------------------------------
    CheckForDuplicates = Structured .OR. GetLogical( Params,'Check for Duplicates')   

    IF( CheckForDuplicates ) THEN
      ALLOCATE( GridPointActive(imaxtot-ioff,jmaxtot-joff,kmaxtot-koff) )
      GridPointActive = .FALSE.

      IF( Structured ) THEN
        ALLOCATE( GridIndex(imaxtot-ioff,jmaxtot-joff,kmaxtot-koff) )
        GridIndex = 0
      END IF
    END IF

    ! It is most convenient to allocate enough at the start but this could 
    ! mean excessive memory usage 
    IF( .NOT. ListGetLogical( Params,'Adaptive Allocation',Found ) ) THEN
      ntot = (imaxtot-ioff)*(jmaxtot-joff)*(kmaxtot-koff)
      CALL AllocateGridParticles( Particles, ntot )
    END IF


    Extent(1) = imintot
    Extent(2) = imaxtot 
    Extent(3) = jmintot
    Extent(4) = jmaxtot
    Extent(5) = kmintot
    Extent(6) = kmaxtot

    

    ! Create particles in the uniform grid
    !----------------------------------------------------------------------------
    
    cands1 = 0
    cands2 = 0


    DO t = ElemStart, ElemFin 

      Element => Mesh % Elements(t)
      n = GetElementNOFNodes(Element)
      CALL GetElementNodes(Nodes,Element)

      IF( MaskExist ) THEN
        IF( ANY( MaskPerm( Element % NodeIndexes ) == 0 ) ) CYCLE
      END IF

      ! Only use the correct dimensional elements for interpolation!
      IF( GetElementDim(Element) /= griddim ) CYCLE


      imin = CEILING( ( MINVAL( Nodes % x(1:n) ) - Origin(1) ) / dx(1) )
      imax = FLOOR( ( MAXVAL( Nodes % x(1:n) ) - Origin(1) ) / dx(1) )
      imin = MAX( imin, imintot )
      imax = MIN( imax, imaxtot )

      IF( griddim >= 2 ) THEN
        jmin = CEILING( ( MINVAL( Nodes % y(1:n) ) - Origin(2) ) / dx(2) )
        jmax = FLOOR( ( MAXVAL( Nodes % y(1:n) ) - Origin(2) ) / dx(2) )
        jmin = MAX( jmin, jmintot )
        jmax = MIN( jmax, jmaxtot )
      ELSE
        jmin = jmintot
        jmax = jmintot

        ! If the element is of reduced order the flatten it in order to make 
        ! PointInElement function better.
        !-------------------------------------------------------------------
        IF( meshdim >= 2 ) THEN
          Nodes % y(1:n) = 0.0_dp
          Nodes % z(1:n) = 0.0_dp
        END IF        
      END IF

      IF( griddim == 3 ) THEN 
        kmin = CEILING( ( MINVAL( Nodes % z(1:n) ) - Origin(3) ) / dx(3) )
        kmax = FLOOR( ( MAXVAL( Nodes % z(1:n) ) - Origin(3) ) / dx(3) )
        kmin = MAX( kmin, kmintot )
        kmax = MIN( kmax, kmaxtot )     
      ELSE
        kmin = kmintot
        kmax = kmaxtot
        IF( meshdim == 3 ) THEN
          Nodes % z(1:n) = 0.0_dp
        END IF
      END IF

      ! The loop is ordered in this way since more often 
      ! nz < ny < nx than any other way minimizing the cost
      !-----------------------------------------------------
      GlobalCoords = 0.0_dp

      DO k=kmin,kmax
        IF( griddim == 3 ) GlobalCoords(3) = k * dx(3) + Origin(3)
        
        DO j=jmin,jmax
          IF( griddim >= 2 ) GlobalCoords(2) = j * dx(2) + Origin(2)
          
          DO i=imin,imax
            GlobalCoords(1) = i * dx(1) + Origin(1)
 
            cands1 = cands1 + 1

            IF( CheckForDuplicates ) THEN
              IF( GridPointActive(i-ioff,j-joff,k-koff) ) CYCLE
              cands2 = cands2 + 1
            END IF
            
            IF ( PointInElement( Element, Nodes, &
                GlobalCoords, LocalCoords ) ) THEN
             
              CALL AddGridParticle(Particles,t,GlobalCoords,LocalCoords)

              IF( CheckForDuplicates ) THEN
                GridPointActive(i-ioff,j-joff,k-koff) = .TRUE.
                IF( Structured ) THEN
                  GridIndex(i-ioff,j-joff,k-koff) = Particles % NumberOfParticles
                END IF
              END IF
            END IF
            
          END DO
        END DO
      END DO
    END DO

    IF( MaskExist ) THEN
      DEALLOCATE( MaskPerm )
    END IF

    IF( CheckForDuplicates ) THEN
      DEALLOCATE( GridPointActive )
    END IF

    totcount(1) = cands1
    totcount(2) = cands2
    totcount(3) = Particles % NumberOfParticles

    IF( ParEnv % PEs > 1 ) THEN
      tmpcount = totcount
      CALL MPI_ALLREDUCE( tmpcount, totcount, 3, MPI_INTEGER, &
          MPI_SUM, ELMER_COMM_WORLD, ierr ) 
    END IF

    WRITE( Message,'(A,I8)') 'Number of candidate nodes:',totcount(1)
    CALL Info('CreateGridParticles',Message,Level=6)
       
    IF( CheckForDuplicates ) THEN
      WRITE( Message,'(A,I8)') 'Number of duplicate nodes:',totcount(1)-totcount(2)
      CALL Info('CreateGridParticles',Message,Level=6)
    END IF

    WRITE( Message,'(A,I8)') 'Number of created nodes:',totcount(3)
    CALL Info('CreateGridParticles',Message,Level=5)
    
    IF( totcount(3) > 0 ) THEN
      WRITE( Message,'(A,F8.2)') 'Search hit fraction:',1.0_dp * totcount(3) / totcount(1)
      CALL Info('CreateGridParticles',Message,Level=6)
    END IF


  END SUBROUTINE CreateGridParticles
     
!------------------------------------------------------------------------------
END SUBROUTINE SaveGridData
!------------------------------------------------------------------------------




