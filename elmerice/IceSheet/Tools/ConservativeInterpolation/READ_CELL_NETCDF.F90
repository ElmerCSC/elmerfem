!/*****************************************************************************/
! *
! *  Elmer/Ice, a glaciological add-on to Elmer
! *  http://elmerice.elmerfem.org
! *
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
! ******************************************************************************
! *
! *  Authors: F. Gillet-Chaulet (IGE-France)
! *  Web:     http://elmerice.elmerfem.org
! *  Original Date: 04/2019
! * 
! *****************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read a cell variable stored in a netcdf file
!  The netcdf should contains all the cell of the serial mesh
!  IF used in parallel, the parallel mesh should have the same global
!  cell ordering as the serial mesh
!  
!  IF netcdf contains a time dimension, the current simulation time is
!  used as time index : if t = ]0._dp,1._dp] => index=1 etc...
!
!  Required input parameters:
!   File Name = File <netcdf file>
!   VarName = File <name of the netcdf variable>
!   Target Variable Name = String OPTIONAL <name of the Elmer variable>
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      SUBROUTINE READ_CELL_NETCDF( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
      USE DefUtils
      USE NETCDF
      IMPLICIT NONE
!------------------------------------------------------------------------------
      TYPE(Solver_t), TARGET :: Solver
      TYPE(Model_t) :: Model
      REAL(KIND=dp) :: dt
      LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
      TYPE(ValueList_t), POINTER :: SolverParams
      TYPE(Variable_t),POINTER :: Var
      TYPE(Element_t), POINTER :: Element
      INTEGER :: i
      INTEGER :: NElements
      CHARACTER (len=100) :: FName
      CHARACTER (len=100) :: VarName,TVarName
      INTEGER :: varid,ncells,ncid,ndims,ntime,TVarID
      INTEGER :: NetCDFstatus
      REAL(KIND=dp), ALLOCATABLE :: Values(:)
      REAL(KIND=dp) :: time
      INTEGER :: TimeIndex
      INTEGER,SAVE :: SavedTime=-1
      INTEGER :: EIndex
      CHARACTER(LEN=MAX_NAME_LEN) :: SolverName="READ_NETCDF_CELL"
      LOGICAL :: Parallel,Found

! get parameters
      SolverParams => GetSolverParams()
      FName = ListGetString(SolverParams,'File Name',UnFoundFatal=.TRUE.)

      write(Message,'(a,a)') 'File name: ',Trim(FName)
      CALL INFO(SolverName,Message,Level=4)

      VarName = ListGetString(SolverParams,'Variable Name',UnFoundFatal=.TRUE.)
      TVarName = ListGetString(SolverParams,'Target Variable Name',Found)
      IF (.NOT.Found) TVarName=VarName

! check if this is a paralell run
      Parallel=(ParEnv % PEs > 1)

! get variable
      Var => VariableGet( Model % Mesh % Variables,TRIM(TVarName),UnFoundFatal=.TRUE.)
      IF(Var % TYPE /= Variable_on_elements) &
        CALL FATAL(SolverName,'Wrong variable type; use -elem ')

      NetCDFstatus = NF90_OPEN(trim(FName),NF90_NOWRITE,ncid)
      IF ( NetCDFstatus /= NF90_NOERR ) &
         CALL FATAL(SolverName,"file open failed")

      NetCDFstatus = nf90_inq_dimid(ncid, 'ncells' , varid)
      IF ( NetCDFstatus /= NF90_NOERR ) &
          CALL FATAL(SolverName,"unable to get ncells dim")

      NetCDFstatus = nf90_inquire_dimension(ncid,varid,len=ncells)
      IF ( NetCDFstatus /= NF90_NOERR ) &
         CALL FATAL(SolverName,"unable to get ncells")

      allocate(Values(ncells))

      ! get variable ID
      NetCDFstatus = nf90_inq_varid(ncid,trim(VarName),TVarId)
      IF ( NetCDFstatus /= NF90_NOERR ) &
          CALL FATAL(SolverName,"unable to get varid")

      ! variable dimensions
      NetCDFstatus = nf90_inquire_variable(ncid, TVarId, ndims=ndims)
      IF ( NetCDFstatus /= NF90_NOERR ) &
          CALL FATAL(SolverName,"unable to get variable dimensions")
      
      ! if ndim > 1 we should have a time dimension
      IF (ndims.GT.1) THEN
        NetCDFstatus = nf90_inq_dimid(ncid, 'time' , varid)
        IF ( NetCDFstatus /= NF90_NOERR ) &
           CALL FATAL(SolverName,"unable to get time dimension")

        NetCDFstatus = nf90_inquire_dimension(ncid,varid,len=ntime)
        IF ( NetCDFstatus /= NF90_NOERR ) &
           CALL FATAL(SolverName,"unable to get ntime")
        
        ! get time index
        time = GetTime()
        dt = GetTimeStepSize()
        TimeIndex = floor(time-dt/2) + 1
        TimeIndex = max(1,min(TimeIndex,ntime))

        ! we have already done this; return
        IF (SavedTime.EQ.TimeIndex) THEN 
          NetCDFstatus=nf90_close(ncid)
          CALL INFO(SolverName,"Nothing to do; return",level=4)
          RETURN
        ENDIF
        SavedTime=TimeIndex
        
        WRITE(Message,'(a,i0)') 'reading time step: ',TimeIndex
        CALL INFO(SolverName,Message,level=4)

        NetCDFstatus=nf90_get_var(ncid,TVarId,Values,start=(/1,TimeIndex/),count=(/ncells,1/))
        IF ( NetCDFstatus /= NF90_NOERR ) &
          CALL FATAL(SolverName,"unable to get variable")

      ELSE
        ! we have already done this; return
        IF (SavedTime.GT.0) THEN
          NetCDFstatus=nf90_close(ncid)
          CALL INFO(SolverName,"Nothing to do; return",level=4)
          RETURN
        ENDIF
        SavedTime=+1

        NetCDFstatus=nf90_get_var(ncid,TVarId,Values,start=(/1/),count=(/ncells/))
        IF ( NetCDFstatus /= NF90_NOERR ) & 
          CALL FATAL(SolverName,"unable to get variable")

      END IF

      ! close file
      NetCDFstatus=nf90_close(ncid)

      NElements = GetNOFActive()
      DO i=1,NElements
        Element => GetActiveElement(i)
        IF (Parallel) THEN 
          EIndex=Element % GElementIndex
        ELSE
          EIndex=Element % ElementIndex
        ENDIF
        Var % Values(Var % Perm(Element % ElementIndex))=Values(EIndex)
      END DO

      DEALLOCATE(Values)


      END SUBROUTINE READ_CELL_NETCDF
