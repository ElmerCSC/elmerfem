!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! * 
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the 
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Authors: Peter Råback, Juha Ruokolainen, 
! *            Samuel Cook, Fabien Gillet-Chaulet , Mondher Chekki 
! *  Email:   elmeradm@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 26.05.2010
! *  Updated to include NetCDF: 05.2021
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
  USE Types

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
  LOGICAL :: Found,VtiFormat, VtuFormat, TableFormat, NetCDFFormat, &
      AnyFormat, RecreateGrid
  INTEGER :: Extent(6)
  REAL(KIND=dp) :: Origin(3), Dx(3)
  INTEGER, POINTER :: GridIndex(:,:,:)
  TYPE( Particle_t), POINTER :: Particles  

  SAVE Visited, GridIndex, Origin, Dx, Extent

#ifdef HAVE_NETCDF
  INTERFACE
    SUBROUTINE ParticleOutputNetCDF( Particles, GridExtent, GridOrigin, GridDx, GridIndex )
      USE Types
      USE ParticleUtils
      USE NetCDF
      TYPE(Particle_t), POINTER :: Particles  
      INTEGER :: GridExtent(6)
      REAL(KIND=dp) :: GridOrigin(3), GridDx(3)
      INTEGER, POINTER :: GridIndex(:,:,:)
    END SUBROUTINE ParticleOutputNetCDF
  END INTERFACE
#endif

!------------------------------------------------------------------------------

  CALL Info('SaveGridData','-----------------------------------------', Level=4 )
  CALL Info('SaveGridData','Saving data on uniform grid point        ', Level=4 )

  Particles => GlobalParticles
  Params => GetSolverParams()
  Mesh => Solver % Mesh

  TableFormat = ListGetLogical( Params,'Table Format',Found)
  VtuFormat = ListGetLogical( Params,'Vtu Format',Found)
  VtiFormat = GetLogical( Params,'Vti Format',GotIt )
  NetCDFFormat =  GetLogical( Params,'NetCDF Format',GotIt )

  FileFormat = ListGetString( Params,'Output Format',Found) 
  IF( Found ) THEN
    IF( FileFormat == 'vtu') VtuFormat = .TRUE.
    IF( FileFormat == 'table') TableFormat = .TRUE.
    IF( FileFormat == 'vti') VtiFormat = .TRUE.
    IF( FileFormat == 'netcdf')  NetCDFFormat = .TRUE. 
  END IF

#ifndef HAVE_NETCDF
  IF( NetCDFFormat ) THEN
    CALL Warn('SaveGridData','Please recompile Elmer with Netcdf library or choose another file format !')
    NetCDFFormat = .FALSE.
  ENDIF 
#endif  

  AnyFormat = VtuFormat .OR. TableFormat .OR. VtiFormat .OR. NetCDFFormat
  IF( .NOT. AnyFormat ) THEN
    CALL Warn('SaveGridData','No active file format given, nothing to do!')
    RETURN
  END IF

  RecreateGrid = ListGetLogical( Params,'Recreate Grid', Found ) 

  ! Initialize the particles on the first calling
  !------------------------------------------------------------------------
  Structured = VtiFormat .OR. NetCDFFormat

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
#ifdef HAVE_NETCDF
  IF( NetCDFFormat ) CALL ParticleOutputNetCDF( Particles, Extent, Origin, Dx, GridIndex )
#endif

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
          I2S(NoParticles),Level=12)    
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
    REAL(KIND=dp) :: MinCoord(3), MaxCoord(3), gMinCoord(3), gMaxCoord(3)
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
    LOGICAL :: CheckForDuplicates, LowerDimensional, MaskExist, MaskOnBulk, GotIt, Parallel=.FALSE.
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

    IF(ParEnv % PEs > 1) Parallel = .TRUE.

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
      ParallelNodes = ParallelReduction( NumberOfNodes ) 
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
      CALL Info('SaveGridData','Assuming reduced coordinate to be: '//I2S(meshdim),Level=5)
    ELSE
      griddim = meshdim
      ElemStart = 1
      ElemFin = Mesh % NumberOfBulkElements 
      ActiveCoordinate = 0
    END IF
    CALL Info('SaveGridData','Saving data on '//I2S(griddim)//'D grid',Level=5)
    
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

    ! We need separately global range (with "g") for determining nx, ny, nz etc.
    ! and the local range to not allocate too much memory. 
    IF( Parallel ) THEN
      DO i=1,3
        gMinCoord(i) = ParallelReduction(MinCoord(i),1)
        gMaxCoord(i) = ParallelReduction(MaxCoord(i),2)
      END DO      
#ifdef HAVE_NETCDF
      IF(NetCDFFormat) THEN
        MinCoord = gMinCoord
        MaxCoord = gMaxCoord
      END IF
#endif 
    ELSE
      gMinCoord = MinCoord
      gMaxCoord = MaxCoord
    END IF

     !print *,'Bounding box min:',MinCoord,ParEnv % myPE
     !print *,'Bounding box max:',MaxCoord,ParEnv % myPE

    ! Optionally the mesh origin may be moved to guarantee that there is 
    ! a node at (x0,y0,z0) always.
    !--------------------------------------------------------------------
    IF( GetLogical( Params,'Grid Origin At Corner',GotIt ) ) THEN
      Origin(1:3) = gMinCoord(1:3)
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
        dx(1) = ( gMaxCoord(1) - gMinCoord(1) ) / nx 
      ELSE
        CALL Fatal('FindGridParticles','Give either > Grid dx < or > Grid nx <')
      END IF
    END IF

    IF( griddim >= 2 ) THEN
      dx(2) = GetCReal( Params,'Grid dy',GotIt) 
      IF(.NOT. GotIt ) THEN
        nx = GetInteger( Params,'Grid ny',GotIt) 
        IF( GotIt) THEN
          dx(2) = ( gMaxCoord(2) - gMinCoord(2) ) / nx
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
          dx(3) = ( gMaxCoord(3) - gMinCoord(3) ) / nx
        ELSE
          dx(3) = dx(1)
        END IF
      END IF
    END IF


    ! Set limits for the global indexes. These are used particularly if the 
    ! bounding box has been manually reduced. Note use of local bounding box
    ! in parallel too
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
        //I2S(imintot)//' - '//I2S(imaxtot),Level=12)
    CALL Info('SaveGridData','Index j range: '&
        //I2S(jmintot)//' - '//I2S(jmaxtot),Level=12)
    CALL Info('SaveGridData','Index k range: '&
        //I2S(kmintot)//' - '//I2S(kmaxtot),Level=12)

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


#ifdef HAVE_NETCDF
!------------------------------------------------------------------------------
!> Writes data out in NetCDF format which assumes a uniform grid where
!> the position of each point is defined by the origin and the grid density. 
!> Also single precision is supported. 
!------------------------------------------------------------------------------
  SUBROUTINE ParticleOutputNetCDF( Particles, GridExtent, GridOrigin, GridDx, GridIndex )
!------------------------------------------------------------------------------

!    USE DefUtils 
!    USE MeshUtils
!    USE ElementDescription
!    USE AscBinOutputUtils
    USE NetCDF
    USE Types
    USE ParticleUtils    

    IMPLICIT NONE
    TYPE(Particle_t), POINTER :: Particles  
    INTEGER :: GridExtent(6)
    REAL(KIND=dp) :: GridOrigin(3), GridDx(3)
    INTEGER, POINTER :: GridIndex(:,:,:)

    !NetCDF variables
    INTEGER :: FileId = -1
    INTEGER, SAVE :: DimId(4)
    INTEGER, SAVE :: VarId(100)
    INTEGER :: NetCDFStatus

    TYPE(ValueList_t),POINTER :: Params
    INTEGER, SAVE :: nTime = 0
    LOGICAL :: GotIt, Parallel, FixedMeshend
    
    CHARACTER(MAX_NAME_LEN), SAVE :: FilePrefix
    CHARACTER(MAX_NAME_LEN) :: NetCDFFile 
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Variable_t), POINTER :: Var
    INTEGER :: i, j, k, Partitions, Part, ExtCount, FileindexOffSet, iTime
    CHARACTER(MAX_NAME_LEN) :: Dir
    REAL(KIND=dp) :: SaveNodeFraction, Time
    LOGICAL :: Found,SinglePrec,NoFileIndex,SuppressDim
    LOGICAL, SAVE :: AllocationDone=.False.
    INTEGER :: NFTYPE
    REAL(Kind=dp) :: FillValue
    REAL :: FillValue_sp
  
    CHARACTER(MAX_NAME_LEN) :: Str
    INTEGER :: NumberOfNodes, ParallelNodes, Dim, ierr
    
    Params => ListGetSolverParams()
    Mesh => GetMesh()
    Time = GetTime()
    
    ExtCount = ListGetInteger( Params,'Output Count',GotIt)
    IF( GotIt ) THEN
      nTime = ExtCount
    ELSE
      nTime = nTime + 1
    END IF
    FileIndexOffset = ListGetInteger( Params,'Fileindex offset',GotIt)
    iTime = nTime + FileIndexOffset

    !Option for user to specify a no-data fill value for particles that fall
    !outside the mesh
    FillValue = ListGetCReal( Params,'No Data Fill Value',GotIt)
    IF(.NOT. GotIt) FillValue = -9999.9_dp

    SinglePrec = GetLogical( Params,'Single Precision',GotIt) 
    IF (SinglePrec) THEN 
       NFTYPE=NF90_FLOAT
    ELSE
       NFTYPE=NF90_DOUBLE
    ENDIF
    IF (SinglePrec) FillValue_sp=FillValue
    NoFileindex = GetLogical( Params,'No Fileindex',GotIt)

    IF ( nTime == 1 ) THEN
      FilePrefix = ListGetString( Params,'Filename Prefix')
      CALL Info('ParticleOutputNetCDF','Saving in NetCDF format to file: ' &
	//TRIM(FilePrefix)//'.nc')
    END IF
    
    Partitions = ParEnv % PEs
    Part = ParEnv % MyPE
    Parallel = (Partitions > 1) .OR. ListGetLogical(Params,'Enforce Parallel format',GotIt)
    
    Dim = Particles % dim
    !This switch will prevent the coordinates for the unused dimension (if there
    !is one) from being written. Not linked up to masks or anything, because you
    !might (not) want the full set of coordinates independently of running on a
    !boundary of a higher-dimensional mesh
    SuppressDim = GetLogical( Params,'Suppress Extra Dimension ',GotIt)
    IF(.NOT. GotIt) SuppressDim = .FALSE.
    IF((GridExtent(2*Dim)-GridExtent(2*Dim-1)) == 0 .AND. SuppressDim) Dim = Dim-1
    
    NumberOfNodes = Particles % NumberOfParticles
     
    IF (LEN_TRIM(Mesh % Name) > 0 ) THEN
      Dir = TRIM(Mesh % Name) // "/"
    ELSE
      Dir = "./"
    END IF
   
    !Default NetCDF behaviour is to write just one file with an unlimited time
    !dimension, so each call to the solver results in all the arrays being saved
    !at the next temporal increment, rather than writing a separate file per
    !timestep 
    !IF( NoFileIndex ) THEN
      WRITE( NetCDFFile,'(A,A,".nc")' ) TRIM(Dir),TRIM(FilePrefix)
    !ELSE
      !WRITE( NetCDFFile,'(A,A,I4.4,".nc")' ) TRIM(Dir),TRIM(FilePrefix),iTime
    !END IF

    CALL WriteNetCDFFile( NetCDFFile )

  CONTAINS


    SUBROUTINE WriteNetCDFFile( NetCDFFile )
      CHARACTER(LEN=*), INTENT(IN) :: NetCDFFile
      TYPE(Variable_t), POINTER :: Var, Solution
      CHARACTER(LEN=512) :: str, DimName
      CHARACTER(LEN=1) :: WorkChar,WorkChar2
      INTEGER :: i,j,k,l,dofs,Rank,cumn,n,vari,sdofs,ind,IsVector,IsAppend,&
                 GridPoints,Offset,DimLen,NumVars2,AllInd,FieldLength,&
                 ValidPart,MPIVP,status
      INTEGER, SAVE :: NumVars
      CHARACTER(LEN=1024) :: Txt, ScalarFieldName, VectorFieldName, FieldName, &
          FieldName2, BaseString, OutStr, WorkString, WorkString2
      LOGICAL :: ScalarsExist, VectorsExist, Found, ParticleMode, ComponentVector, &
          ComplementExists, ThisOnly, Stat, WriteData, WriteXML, NotInBoss
      INTEGER, POINTER :: Perm(:), Perm2(:), Indexes(:)
      INTEGER, ALLOCATABLE,SAVE :: ElemInd(:),ElemInd2(:)
      REAL(KIND=dp), POINTER :: ScalarValues(:), VectorValues(:,:),Values(:),Values2(:),&
          Values3(:), Values4(:), Values5(:), Values6(:), Values7(:),&
          Values8(:), Values9(:)
      REAL(KIND=dp) :: x,y,z,u,v,w,DetJ,val
      REAL :: fvalue
      TYPE(Nodes_t),SAVE :: Nodes      
      TYPE(Element_t), POINTER :: Element
      REAL(KIND=dp),SAVE,ALLOCATABLE :: Array(:,:,:),PArray(:,:,:),Basis(:)
      REAL(KIND=dp) :: rt,rt0,rtc
      INTEGER :: nx,ny,nz

      ! Initialize the NetCDF file for writing (only in boss partition)
      ! Or, if this not first call, open the existing file
      !--------------------------------------------------------------
      IF(Part == 0 .OR. .NOT. Parallel) THEN
        IF(nTime == 1) THEN
          NetCDFStatus = NF90_CREATE(NetCDFFile, 0, FileId)
          IF ( NetCDFStatus /= 0 ) THEN
            CALL Fatal( 'WriteNetCDFFile', 'NetCDF file could not be created: '//TRIM(NetCDFFile))
          END IF
        ELSE
          NetCDFStatus = NF90_OPEN(NetCDFFile, NF90_WRITE, FileId)
          IF ( NetCDFStatus /= 0 ) THEN
            CALL Fatal( 'WriteNetCDFFile', 'NetCDF file could not be opened: '//TRIM(NetCDFFile))
          END IF

        END IF
      END IF
      
      nx=(GridExtent(2)-GridExtent(1))+1
      ny=(GridExtent(4)-GridExtent(3))+1
      nz=(GridExtent(6)-GridExtent(5))+1
      
      IF (.NOT.AllocationDone) THEN
        Allocate(Array(nx,ny,nz))
        IF (Parallel.AND.(Part == 0)) Allocate(PArray(nx,ny,nz))
        n = Mesh % MaxElementNodes
        ALLOCATE( Basis(n), Nodes % x(n), Nodes % y(n), Nodes % z(n) )

        n = Mesh % MaxElementDOFS
        ALLOCATE( ElemInd(n), ElemInd2(n) )
        AllocationDone=.TRUE.
      ENDIF

      ThisOnly = .TRUE.

      !Set up dims here beforehand if first time
      IF(nTime==1 .AND. (Part == 0 .OR. .NOT. Parallel)) THEN

        IF(Dim==2) DimId(4) = 0
        DO i=1,Dim+1
          IF(i==1) THEN
            DimName = 'Time'
          ELSE IF(i==2) THEN
            DimName = 'x'
            DimLen = nx
          ELSE IF(i==3) THEN
            DimName = 'y'
            DimLen = ny
          ELSE IF(i==4) THEN
            DimName = 'z'
            DimLen = nz
          ELSE
            CALL Fatal( 'WriteNetCDFFile', 'Are you sure your glacier has more than 3 dimensions?')
          END IF
          IF(i==1) THEN
            NetCDFStatus = NF90_DEF_DIM(FileId, DimName, nf90_unlimited, DimId(i))
            IF ( NetCDFStatus /= 0 ) THEN
              CALL Fatal( 'WriteNetCDFFile', 'NetCDF dimension could not be created: '//TRIM(NetCDFFile))
            END IF
          ELSE
            NetCDFStatus = NF90_DEF_DIM(FileId, DimName, DimLen, DimId(i))
            IF ( NetCDFStatus /= 0 ) THEN
              CALL Fatal( 'WriteNetCDFFile', 'NetCDF dimension could not be created: '//TRIM(NetCDFFile))
            END IF
          END IF
        END DO
      END IF

      !Check there are actually some variables
      ScalarFieldName = ListGetString( Params,'Scalar Field 1',ScalarsExist)
      VectorFieldName = ListGetString( Params,'Vector Field 1',VectorsExist)
      IF( .NOT. ScalarsExist .AND. .NOT. VectorsExist) THEN
        CALL Warn('WriteNetCDFFile','Are there really no scalars or vectors?')
        WriteData = .FALSE.
        RETURN
      END IF
      WriteData = .TRUE.

      !Create coordinate variables if first time so you know where the grid
      !came from      
      IF(nTime==1) NumVars = 1

      IF(nTime==1 .AND. (Part == 0 .OR. .NOT. Parallel)) THEN
        DO i=1,Dim+1
          IF(i==1) DimName = 'Time'
          IF(i==2) DimName = 'x'
          IF(i==3) DimName = 'y'
          IF(i==4) DimName = 'z'

          NetCDFStatus = NF90_DEF_VAR(FileId, DimName, NF90_DOUBLE, (/ DimId(i) /),VarId(NumVars))
          IF ( NetCDFStatus /= 0 ) THEN
             CALL Fatal( 'WriteNetCDFFile', 'NetCDF variable could not be created: '//TRIM(DimName))
          END IF
          NumVars = NumVars + 1
        END DO

        !Create all the variables in the NetCDF if first time. Vectors will have
        !an unknown number of components, so just loop over them until found
        !them all  
        DO i=1,2
          IF(i==1) BaseString = 'Scalar Field'
          IF(i==2) BaseString = 'Vector Field'
          DO Vari= 1, 99
            WRITE(Txt,'(A)') TRIM(BaseString)//' '//I2S(Vari)
            IF(i==1) THEN          
              FieldName = ListGetString( Params, TRIM(Txt), Found )
              IF(.NOT. Found) EXIT
              IF(Dim==2) THEN
                NetCDFStatus = NF90_DEF_VAR(FileId, TRIM(FieldName), NFTYPE,&
                             (/ DimId(2), DimId(3), DimId(1) /),VarId(NumVars))
                IF ( NetCDFStatus /= 0 ) THEN
                  CALL Fatal( 'WriteNetCDFFile', 'NetCDF variable could not be created: '//TRIM(FieldName))
                END IF
                IF (SinglePrec) THEN
                  NetCDFStatus = NF90_DEF_VAR_Fill(FileId, VarId(NumVars), 0, FillValue_sp)
                ELSE
                  NetCDFStatus = NF90_DEF_VAR_Fill(FileId, VarId(NumVars), 0, FillValue)
                ENDIF
                IF ( NetCDFStatus /= NF90_NOERR ) THEN
                  CALL Fatal( 'WriteNetCDFFile', 'NetCDF no-data fill value could not be defined: '//TRIM(FieldName))
                END IF
              ELSE IF(Dim==3) THEN
                NetCDFStatus = NF90_DEF_VAR(FileId, TRIM(FieldName), NFTYPE,&
                             (/ DimId(2), DimId(3), DimId(4), DimId(1) /),VarId(NumVars))
                IF ( NetCDFStatus /= 0 ) THEN
                CALL Fatal( 'WriteNetCDFFile', 'NetCDF variable could not be created: '//TRIM(FieldName))
                END IF
                IF (SinglePrec) THEN
                  NetCDFStatus = NF90_DEF_VAR_Fill(FileId, VarId(NumVars), 0, FillValue_sp)
                ELSE
                  NetCDFStatus = NF90_DEF_VAR_Fill(FileId, VarId(NumVars), 0, FillValue)
                ENDIF
                IF ( NetCDFStatus /= 0 ) THEN
                  CALL Fatal( 'WriteNetCDFFile', 'NetCDF no-data fill value could not be defined: '//TRIM(FieldName))
                END IF
              END IF
              NumVars = NumVars + 1
            ELSE IF(i==2) THEN
              FieldName = ListGetString( Params, TRIM(Txt), Found )
              IF(.NOT. Found) EXIT
              DO j=1,10
                Solution => VariableGet( Mesh % Variables, TRIM(FieldName)//' '//I2S(j),ThisOnly )
                IF( .NOT. ASSOCIATED(Solution)) THEN
                  EXIT
                END IF
                IF(Dim==2) THEN
                  NetCDFStatus =  NF90_DEF_VAR(FileId,&
                                TRIM(FieldName)//' '//TRIM(I2S(j)), NFTYPE,&
                                (/ DimId(2), DimId(3), DimId(1) /),VarId(NumVars))
                  IF ( NetCDFStatus /= 0 ) THEN
                    CALL Fatal( 'WriteNetCDFFile', 'NetCDF variable could not be created: '//TRIM(FieldName))
                  END IF
                  IF (SinglePrec) THEN
                    NetCDFStatus = NF90_DEF_VAR_Fill(FileId, VarId(NumVars), 0, FillValue_sp)
                  ELSE
                    NetCDFStatus = NF90_DEF_VAR_Fill(FileId, VarId(NumVars), 0, FillValue)
                  ENDIF
                  IF ( NetCDFStatus /= 0 ) THEN
                    CALL Fatal( 'WriteNetCDFFile', 'NetCDF no-data fill value could not be defined: '//TRIM(FieldName))
                  END IF
                ELSE IF(Dim==3) THEN
                  NetCDFStatus =  NF90_DEF_VAR(FileId,&
                                TRIM(FieldName)//' '//I2S(j), 6,&
                                (/ DimId(2), DimId(3), DimId(4), DimId(1) /),VarId(NumVars))
                  IF ( NetCDFStatus /= 0 ) THEN
                    CALL Fatal( 'WriteNetCDFFile', 'NetCDF variable could not be created: '//TRIM(FieldName))
                  END IF
                  NetCDFStatus = NF90_DEF_VAR_Fill(FileId, VarId(NumVars), 0, FillValue)
                  IF ( NetCDFStatus /= 0 ) THEN
                    CALL Fatal( 'WriteNetCDFFile', 'NetCDF no-data fill value could not be defined: '//TRIM(FieldName))
                  END IF
                END IF
                NumVars = NumVars + 1
              END DO
            END IF

          END DO
        END DO
      END IF

      !All definitions done, on to writing values
      IF(nTime==1 .AND. (Part == 0 .OR. .NOT. Parallel)) THEN
        NetCDFStatus = NF90_ENDDEF(FileId)
        IF ( NetCDFStatus /= 0 ) THEN
          CALL Fatal( 'WriteNetCDFFile', 'NetCDF file could not be put in data-entry mode: '//TRIM(NetCDFFile))
        END IF
      END IF

      NumVars2 = 1

      !Add coordinate values first. Only time gets updated if not first call
      IF(Part == 0 .OR. .NOT. Parallel) THEN
        DO j=1,Dim+1
          IF(j==1) THEN
              NetCDFStatus = NF90_PUT_VAR(FileId, VarId(NumVars2), Time, start=(/ nTime /))
              IF ( NetCDFStatus /= 0 ) THEN
                CALL Fatal( 'WriteNetCDFFile', 'NetCDF variable could not be written: '//TRIM(FieldName))
              END IF
            NumVars2 = NumVars2 + 1
          ELSE
            IF(nTime==1) THEN
              DO i=1,INT((GridExtent(2*(j-1))-GridExtent((2*(j-1))-1))+1)
                val = GridExtent((2*(j-1))-1) + (GridDx(j-1) * (i-1)) + GridOrigin(j-1)
                NetCDFStatus = NF90_PUT_VAR(FileId, VarId(NumVars2), val, start=(/ i /))
                IF ( NetCDFStatus /= 0 ) THEN
                  CALL Fatal( 'WriteNetCDFFile', 'NetCDF variable could not be written: '//TRIM(FieldName))
                END IF
              END DO
            END IF
            NumVars2 = NumVars2 + 1
          END IF
        END DO
      END IF
 
      !Send NumVars and NumVars2 from boss to all partitions so that everyone
      !runs the same sized loop
      IF(Parallel) THEN
        CALL MPI_BCAST(NumVars, 1, MPI_INTEGER, 0, ELMER_COMM_WORLD, ierr)
        CALL MPI_BCAST(NumVars2, 1, MPI_INTEGER, 0, ELMER_COMM_WORLD, ierr)
      END IF      

      l = NumVars2

      DO Vari = l, NumVars-1

        !---------------------------------------------------------------------
        ! We've already read the variable list when defining the variables in
        ! the NetCDF file, so here we can just read back the variable names
        ! from the file and go get the variable in the usual manner, without
        ! having to distinguish between scalars and vectors
        !--------------------------------- -----------------------------------
        !---------------------------------------------------------------------
        ! Find the variable with the given name in the normal manner 
        !---------------------------------------------------------------------
        
        FieldLength = 0
        WorkChar2 = 'x'
        IF(NumVars2 > NumVars-1) EXIT
        IF(Part == 0 .OR. .NOT. Parallel) THEN
          NetCDFStatus = NF90_INQUIRE_VARIABLE(FileId, VarId(NumVars2), FieldName)
          IF ( NetCDFStatus /= 0 ) THEN
            CALL Fatal( 'WriteNetCDFFile', 'NetCDF variable could not be read: '//TRIM(FieldName))
          END IF
          FieldLength = LEN(TRIM(FieldName))
        END IF

        !To avoid having to open the same file in multiple partitions and
        !possibly cause all kinds of errors, this will broadcast each separate
        !character of fieldname and then re-concatenate them so all processes
        !know what to look for. This feels very stupid.
        IF(Parallel) THEN
          CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)
          CALL MPI_BCAST(FieldLength, 1, MPI_INTEGER, 0, ELMER_COMM_WORLD, ierr)
          DO i=1,FieldLength
            IF(Part == 0) WorkChar = FieldName(i:i)
            CALL MPI_BCAST(WorkChar, 1, MPI_CHARACTER, 0, ELMER_COMM_WORLD, ierr)
            IF(Part .NE. 0 .AND. i==1) WorkString = WorkChar

            !We have to consider, certainly in vectors, that there will be a
            !space in the middle of the filename that TRIM will remove. These
            !IFs ensure it's re-inserted
            IF(WorkChar == ' ') WorkChar2 = WorkChar
            IF(Part .NE. 0 .AND. i .NE. 1) THEN
              IF(WorkChar2 == ' ' .AND. WorkChar .NE. ' ') THEN
                WorkString = TRIM(WorkString)//WorkChar2//WorkChar
                WorkChar2 = 'x'
              ELSE
                WorkString = TRIM(WorkString)//WorkChar
              END IF
            END IF
          END DO
          IF(Part .NE. 0) FieldName = TRIM(WorkString)
        END IF

        !Actually get the variable!
        Solution => VariableGet( Mesh % Variables,TRIM(FieldName),ThisOnly )
        IF( .NOT. ASSOCIATED( Solution ) ) THEN
          WRITE(Txt, '(A,A)') 'Nonexistent variable: ',TRIM(FieldName)
          CALL Warn('WriteNetCDFFile', Txt)
          CYCLE
        END IF

        !Because we've already separated out the vector components, every
        !variable will have dofs = 1
        Perm => Solution % Perm
        dofs = 1 !Solution % DOFs
        Values => Solution % Values
          
        !---------------------------------------------------------------------
        ! Eigenmodes have not yet been implemented
        !---------------------------------------------------------------------
        IF( ASSOCIATED(Solution % EigenVectors)) THEN
          CALL Warn('WriteNetCDFFile','Do the eigen values')
        END IF

        !---------------------------------------------------------------------
        ! There may be special complementary variables such as 
        ! displacement & mesh update.
        ! I think this should still function in NetCDF, but I haven't tested it,
        ! so be alert. 
        !---------------------------------------------------------------------
        ComplementExists = .FALSE.

        FieldName2 = ListGetString( Params, TRIM(FieldName)//' Complement', Found )
        IF( Found ) THEN
          Solution => VariableGet( Mesh % Variables, &
              TRIM(FieldName2), ThisOnly )
          IF( ASSOCIATED(Solution)) THEN 
            Values2 => Solution % Values
            Perm2 => Solution % Perm 
            ComplementExists = .TRUE.
          ELSE
            CALL Warn('WriteNetCDFFile','Complement does not exist:'//TRIM(FieldName2))
          END IF
        END IF
              
        !---------------------------------------------------------------------
        ! Finally save the field values for scalars and vectors
        !---------------------------------------------------------------------
        IF( WriteData ) THEN
          Array=-HUGE(1.0_dp)
          !PArray=Array
          DO k = 1,nz
            DO j = 1,ny
              DO i = 1,nx

                ind = GridIndex( i, j, k ) 
                IF(ind.GT.0) THEN

                    Element => Mesh % Elements( Particles % ElementIndex(ind) )            
                    IF ( Solution % TYPE == Variable_on_elements ) THEN
                      val = Values(Perm(Element % ElementIndex))
                    ELSE
                      Indexes => Element % NodeIndexes
                      n = Element % TYPE % NumberOfNodes
                    
                      Nodes % x(1:n) = Mesh % Nodes % x( Indexes ) 
                      Nodes % y(1:n) = Mesh % Nodes % y( Indexes ) 
                      Nodes % z(1:n) = Mesh % Nodes % z( Indexes ) 
                    
                      u = Particles % uvw(ind,1)
                      v = Particles % uvw(ind,2)
                      IF( dim == 3 ) THEN
                        w = Particles % uvw(ind,3)
                      ELSE
                        w = 0.0_dp
                      END IF
                      stat = ElementInfo( Element,Nodes,u,v,w,detJ,Basis)
                                        
                      IF( Solution % TYPE == Variable_on_nodes_on_elements ) THEN
                        ElemInd(1:n) = Perm( Element % DGIndexes(1:n) )
                        IF( ComplementExists ) THEN
                          ElemInd2(1:n) = Perm2( Element % DGIndexes(1:n) )
                        END IF
                      ELSE 
                        ElemInd(1:n) = Perm( Indexes(1:n) )
                        IF( ComplementExists ) THEN
                          ElemInd2(1:n) = Perm2( Indexes(1:n) )
                        END IF
                      END IF

                      val = 0.0_dp
                      IF( ALL( ElemInd(1:n) > 0 ) ) THEN
                        val = SUM( Basis(1:n) * Values(dofs*(ElemInd(1:n)-1)+1) )
                      ELSE IF( ComplementExists ) THEN
                        IF( ALL( ElemInd2(1:n) > 0 ) ) THEN
                          val = SUM( Basis(1:n) * Values2(dofs*(ElemInd2(1:n)-1)+1) )
                        END IF
                      END IF
                    END IF

                    Array(i,j,k)=val
                  END IF

              END DO ! i
            END DO ! j
          END DO ! k

          IF(Parallel) THEN
            CALL MPI_REDUCE(Array,PArray,nx*ny*nz,MPI_DOUBLE,MPI_MAX,0,ELMER_COMM_WORLD, ierr)
            IF(Part == 0) Array=PArray
          END IF
        
          IF(Part == 0 .OR. (.NOT.Parallel)) THEN
            !Array=PArray
            WHERE(Array.EQ.-HUGE(1.0_dp)) Array=FillValue
            IF(Dim == 2) THEN
               NetCDFStatus = NF90_PUT_VAR(FileId, VarId(NumVars2), Array(:,:,1), start=(/ 1,1,nTime /))
            ELSE IF(Dim == 3) THEN
               NetCDFStatus = NF90_PUT_VAR(FileId, VarId(NumVars2), Array(:,:,:), start=(/ 1,1,1,nTime /))
            END IF
            IF ( NetCDFStatus /= 0 ) THEN
               CALL Fatal( 'WriteNetCDFFile', 'NetCDF variable could not be written: '//TRIM(FieldName))
            END IF
          END IF

        END IF
        NumVars2 = NumVars2 + 1
      END DO
        
      !Boss only
      IF(Part == 0 .OR. .NOT. Parallel) THEN
        NetCDFStatus = NF90_CLOSE(FileId)
        IF ( NetCDFStatus /= 0 ) THEN
          CALL Fatal( 'WriteNetCDFFile', 'NetCDF file could not be closed: '//TRIM(NetCDFFile))
        END IF
      END IF

    END SUBROUTINE WriteNetCDFFile
      
!----------------------------------------------------------------------------
  END SUBROUTINE ParticleOutputNetCDF
!----------------------------------------------------------------------------
#endif


