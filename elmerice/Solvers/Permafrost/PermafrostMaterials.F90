!/*****************************************************************************/
! *
! *  Elmer/Ice, a glaciological add-on to Elmer
! *  http://elmerice.elmerfem.org
! *
! * 
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) asny later version.
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
! *  Authors: Thomas Zwinger, Denis Cohen, Juha Hartikainen
! *  Email:  thomas Zwinger [at] csc.fi 
! *  Web:     http://elmerice.elmerfem.org
! *  Address: CSC - Scientific Computing Ltd.  
! *               Keilaranta 14                    
! *               02101 Espoo, Finland             
! *                                                 
! *       Original Date:  January 2017  -               
! * 
! *****************************************************************************
!>  Module containing material declarations and material functions
!>  for enhanced permafrost problem
! *****************************************************************************
MODULE PermafrostMaterials
  
  USE Types
  USE DefUtils
  USE SolverUtils
  IMPLICIT NONE
  !---------------------------------
  ! type for solvent (water and ice)
  !---------------------------------
  TYPE SolventMaterial_t
     REAL(KIND=dp) :: &
          Mw,rhow0,rhoi0,hw0,hi0,vi0,bccw0,&
          Ei0, nui0, betai, &
          kw0th,ki0th,bw,bi, &
          cw0,acw(0:5),bcw(0:5), &
          ci0,aci(0:5),&
          aw0,kw0,zw0,aaw(0:5),bzw(0:5),ckw(0:5), &
          ai0,ki0,aai(0:5),cki(0:5),&
          muw0,nu10,anw(0:5),bnw(0:5)
     INTEGER :: &
          acwl,bcwl,aawl,bzwl,ckwl,&
          acil,aail,ckil,anwl,bnwl 
  END type SolventMaterial_t
  
  !---------------------------------
  ! type for solute (ions)
  !---------------------------------
  TYPE SoluteMaterial_t
     REAL(KIND=dp) ::  Mc,vc0,kc0th,&
          d1,d2,bc,&
          cc0,acc(0:5),bcc(0:5),&
          rhoc0,ac0,kc0,zc0,aac(0:5),ckc(0:5),bzc(0:5),&
          nu20,anc(0:5),bnc(0:5)
     INTEGER :: accl, bccl,aacl,ckcl,bzcl,ancl,bncl
     CHARACTER(LEN=MAX_NAME_LEN) :: SoluteName
  END type SoluteMaterial_t
  
  !---------------------------------
  ! type for rock material
  !---------------------------------
  TYPE RockMaterial_t
     INTEGER :: NumberOfRockRecords
     REAL(KIND=dp), ALLOCATABLE :: ks0th(:),e1(:),bs(:),rhos0(:),&
          Xi0(:),eta0(:),etak(:),hs0(:),Kgwh0(:,:,:),qexp(:),alphaL(:),alphaT(:),RadGen(:),&
          cs0(:),acs(:,:),as0(:),aas(:,:),ks0(:),cks(:,:),Es0(:),nus0(:),betas(:)
     INTEGER, ALLOCATABLE :: acsl(:),aasl(:),cksl(:)
     CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: VariableBaseName(:)
  END TYPE RockMaterial_t

  TYPE(SolventMaterial_t), TARGET :: GlobalSolventMaterial
  TYPE(SoluteMaterial_t), TARGET :: GlobalSoluteMaterial
  TYPE(RockMaterial_t) :: GlobalRockmaterial
  
  
CONTAINS
  !---------------------------------------------------------------------------------------------
  !---------------------------------------------------------------------------------------------
  !---------------------------------------------------------------------------------------------

  !-------------------------------------------------
  ! I/O related functions
  !-------------------------------------------------

  SUBROUTINE SetPermafrostSolventMaterial( CurrentSolventMaterial)
    IMPLICIT NONE
    TYPE(ValueList_t), POINTER :: Params, Constants
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    ! ----------- local
    TYPE(SolventMaterial_t), TARGET :: LocalSolventMaterial
    LOGICAL :: FirstTime=.TRUE.
    CHARACTER(LEN=MAX_NAME_LEN) :: SubroutineName='SetPermafrostSolventMaterial'
    SAVE LocalSolventMaterial, FirstTime

    IF (FirstTime) THEN
      !------------------------------------------------------------------------------
      ! set constants for water and ice 
      ! Mw,rhow0,rhoi0,hw0,hi0,vi0,cw0,ci0,acw(3),bcw(0:2),aci(0:1),kw0th,ki0th, bi, bw
      !------------------------------------------------------------------------------
      LocalSolventMaterial % Mw =    1.8015d-2      
      LocalSolventMaterial % hw0 =   0.0_dp      
      LocalSolventMaterial % hi0 =  -333360.0_dp !!
      
      ! --------------------- polynomials

      !!! water !!!
      
      ! heat capacity water      
      LocalSolventMaterial % cw0  = 4207.7_dp
      LocalSolventMaterial % acw(0:5) = &
           RESHAPE([1.0_dp,-0.0887_dp,0.2859_dp,0.0_dp,0.0_dp,0.0_dp], &
           SHAPE(LocalSolventMaterial % acw))
      LocalSolventMaterial % acwl=2      
      LocalSolventMaterial % bcw(0:5) = &
           RESHAPE([1.0_dp,1.5852_dp,8.0686_dp,0.0_dp,0.0_dp,0.0_dp],&
           SHAPE(LocalSolventMaterial % bcw))
      LocalSolventMaterial % bcwl=2

      !heat conductivity of water
      LocalSolventMaterial % kw0th = 0.56_dp 
      LocalSolventMaterial % bw = 0.0_dp
      
      ! density water
      LocalSolventMaterial % rhow0 = 999.9_dp ! density at reference temperature
      LocalSolventMaterial % kw0  = 4.4534d-10 ! Isothermal compressibility

      LocalSolventMaterial % ckw(0:5) = &
           RESHAPE([1.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp],&
           SHAPE(LocalSolventMaterial % ckw))
      LocalSolventMaterial % ckwl=0      
      LocalSolventMaterial % aw0  = -5.3358d-05 ! Isobaric thermal expansion
      LocalSolventMaterial % aaw(0:5) = &
           RESHAPE([1.0_dp,-79.1305_dp,207.4836_dp,-403.8270_dp,395.5347_dp,-166.1466_dp],&
           SHAPE(LocalSolventMaterial % aaw))
      LocalSolventMaterial % aawl=5
      LocalSolventMaterial % zw0  = -2.0217d-01 ! Isothermal chemical compaction
      LocalSolventMaterial % bzw(0:5) = &
           RESHAPE([1.0_dp,12.8298_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp],&
           SHAPE(LocalSolventMaterial % bzw))
      LocalSolventMaterial % bzwl=1

      !viscosity water
      LocalSolventMaterial % muw0 = 1.7914d-03   ! viscosity at reference temperature
      LocalSolventMaterial % nu10 = -0.034376_dp ! temperature dependence of viscosity
      LocalSolventMaterial % anw(0:5) = &
           RESHAPE([1.0_dp,-2.3302_dp,4.0084_dp,-2.9697_dp,0.0_dp,0.0_dp],&
           SHAPE(LocalSolventMaterial % anw))
      LocalSolventMaterial % anwl=3
      LocalSolventMaterial % bnw(0:5) = &
           RESHAPE([1.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp],&
           SHAPE(LocalSolventMaterial % bnw))
      LocalSolventMaterial % bnwl=0

      !!!!! ice !!!!!!
      
      ! density ice
      LocalSolventMaterial % rhoi0 = 916.8_dp !reference density
      LocalSolventMaterial % ki0  = 1.1417d-10 ! Isothermal compressibility
      LocalSolventMaterial % cki(0:5) = &
           RESHAPE([1.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp],&
	    SHAPE(LocalSolventMaterial % cki))
      LocalSolventMaterial % ckil=0
      LocalSolventMaterial % ai0  = 1.6781d-04 ! Isobaric thermal expansion
      LocalSolventMaterial % aai(0:5) = &
           RESHAPE([1.0_dp,1.1923_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp], &
	   SHAPE(LocalSolventMaterial % aai))
      LocalSolventMaterial % aail=1
      
      ! specivic volume ice
      LocalSolventMaterial % vi0 = 1.0_dp/(LocalSolventMaterial % rhoi0) ! reference specific volume
      
      ! heat conductivity ice
      LocalSolventMaterial % ki0th = 2.24_dp!!      
      LocalSolventMaterial % bi = 0.0_dp

      ! heat capacity     
      LocalSolventMaterial % ci0  = 2088.8_dp ! reference value
      LocalSolventMaterial % aci(0:5) = &
           RESHAPE([1.0_dp,0.9557_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp], &
	   SHAPE(LocalSolventMaterial % aci))
      LocalSolventMaterial % acil=1

      ! elastic properties
      LocalSolventMaterial % Ei0 = 9.33d09
      LocalSolventMaterial % nui0= 0.325_dp
      LocalSolventMaterial % betai = 0.0_dp  ! CHANGE

      CALL INFO(SubroutineName,"-----------------------------------------------------------------",Level=9)
      CALL INFO(SubroutineName,"Solvent related constants",Level=9)
      WRITE(Message,*) "Mw",LocalSolventMaterial % Mw,"rhow0",LocalSolventMaterial % rhow0,"rhoi0",LocalSolventMaterial % rhoi0,&
           "hw0",LocalSolventMaterial % hw0,"hi0",LocalSolventMaterial % hi0,"vi0",LocalSolventMaterial % vi0,&
           "cw0",LocalSolventMaterial % cw0,"ci0",LocalSolventMaterial % ci0,"acw(3)",LocalSolventMaterial % acw(0:2)
      CALL INFO(SubroutineName,Message,Level=9)
      WRITE(Message,*) "bcw(0:2)",LocalSolventMaterial % bcw(0:2),"aci(0:1)",LocalSolventMaterial % aci(0:1),&
           "kw0th",LocalSolventMaterial % kw0th,"ki0th",LocalSolventMaterial % ki0th," bi",LocalSolventMaterial % bi,&
           "bw",LocalSolventMaterial % bw
      CALL INFO(SubroutineName,Message,Level=9)
      CALL INFO(SubroutineName,"-----------------------------------------------------------------",Level=9)
      FirstTime = .FALSE.
    END IF
    CurrentSolventMaterial => LocalSolventMaterial
  END SUBROUTINE SetPermafrostSolventMaterial
  
  !---------------------------------------------------------------------------------------------
  SUBROUTINE ReadPermafrostSoluteMaterial( Params,Constants,CurrentSoluteMaterial )
    IMPLICIT NONE
    TYPE(ValueList_t), POINTER :: Params, Constants
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    ! ----------- local
    TYPE(SoluteMaterial_t), TARGET :: LocalSoluteMaterial
    INTEGER :: i,j,k,l, n,t, active, DIM, ok,InitialNumerOfSoluteRecords, EntryNumber
    INTEGER,parameter :: io=20
    LOGICAL :: Found, DataRead=.FALSE.
    CHARACTER(LEN=MAX_NAME_LEN) ::  SoluteFileName, Comment
    CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SubroutineName='ReadPermafrostSoluteMaterial'

    SAVE DataRead,SoluteFileName, LocalSoluteMaterial

    IF (DataRead) THEN
      CurrentSoluteMaterial => LocalSoluteMaterial   
      RETURN
    ELSE
      DIM = CoordinateSystemDimension()
      !------------------------------------------------------------------------------
      ! Inquire and open file
      !------------------------------------------------------------------------------
      ! give preference to a defined material database
      SoluteFileName = GetString( Params, 'Solute Material File', Found )
      IF (.NOT.Found) THEN
        CALL INFO(SubroutineName," 'Solute Material File' keyword not found - assigning default values!",Level=1)
        DataRead=.TRUE.
        LocalSoluteMaterial % SoluteName = TRIM('Sea Salt')
        LocalSoluteMaterial % Mc    = 0.031404_dp
        LocalSoluteMaterial % rhoc0 = 2206.6_dp 
        LocalSoluteMaterial % vc0   = 1.0/(LocalSoluteMaterial % rhoc0)
        LocalSoluteMaterial % kc0th = 0.56_dp      
        LocalSoluteMaterial % d1    = 0.87_dp		       
        LocalSoluteMaterial % d2    = 2.00_dp
        LocalSoluteMaterial % bc    = 0.0_dp
        ! heat capacity polynomials
        LocalSoluteMaterial % cc0   = 1906.6_dp
        LocalSoluteMaterial % acc(0:5) = &
             RESHAPE([1.0_dp,-0.0887_dp,0.2859_dp,0.0_dp,0.0_dp,0.0_dp],&
	      SHAPE(LocalSoluteMaterial % acc))
        LocalSoluteMaterial % accl = 2
        LocalSoluteMaterial % bcc(0:5) = &
             RESHAPE([1.0_dp,-1.5852_dp,8.0686_dp,0.0_dp,0.0_dp,0.0_dp],&
	      SHAPE(LocalSoluteMaterial % bcc))
        LocalSoluteMaterial % bccl = 2
        ! density        
        LocalSoluteMaterial % ac0 =   -5.3358d-05 ! thermal expansion
        LocalSoluteMaterial % kc0 =   4.4534d-10  ! compressibility
        LocalSoluteMaterial % zc0 =   2.0217d-01  ! chemical compaction
        LocalSoluteMaterial % aac = &
             RESHAPE([1.0_dp, -79.1305_dp, 207.4835_dp, -403.827_dp, 395.5347_dp, -166.1466_dp],&
	      SHAPE(LocalSoluteMaterial % aac))
        LocalSoluteMaterial % aacl = 5
        LocalSoluteMaterial % ckc = &
             RESHAPE([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp],&
	      SHAPE(LocalSoluteMaterial % ckc))
        LocalSoluteMaterial % ckcl = 0 
        LocalSoluteMaterial % bzc = &
             RESHAPE([1.0_dp, -12.8298_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp],&
	      SHAPE(LocalSoluteMaterial % bzc))
        LocalSoluteMaterial % bzcl = 1
        ! viscosity
        LocalSoluteMaterial % nu20 = 2.6870_dp !influence of salinity on viscosity
        LocalSoluteMaterial % anc = &
             RESHAPE([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp],&
	      SHAPE(LocalSoluteMaterial % anc))
        LocalSoluteMaterial % ancl = 0
        LocalSoluteMaterial % bnc= &
             RESHAPE([1.0_dp, 6.7992_dp, -31.3293_dp, 44.7717_dp, 0.0_dp, 0.0_dp],&
	      SHAPE(LocalSoluteMaterial % bnc))
        LocalSoluteMaterial % bncl = 3
        DataRead=.TRUE.
      ELSE      
        OPEN(unit = io, file = TRIM(SoluteFileName), status = 'old',iostat = ok)
        IF (ok /= 0) THEN
          WRITE(Message,'(A,A)') 'Unable to open file ',TRIM(SoluteFileName)
          CALL FATAL(Trim(SubroutineName),Trim(message))
        ELSE
          !------------------------------------------------------------------------------
          ! Read in the number of records in file (first line integer)
          !------------------------------------------------------------------------------
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % SoluteName
          WRITE(Message,'(A,A)') 'Reading entry', TRIM(LocalSoluteMaterial % SoluteName)
          CALL INFO(Trim(SubroutineName),Trim(Message),Level=3)
          !------------------------------------------------------------------------------
          ! Read in information for (currently fixed) solute (= salts)
          ! 
          !     Mc, rhoc0
          !     vc0,cc0,acc(0:2),bcc(0:2)
          !     kc0th,muw0,d1,d2,bc
          !------------------------------------------------------------------------------
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % Mc, Comment 
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % vc0, Comment
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % kc0th, Comment
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % d1, Comment
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % d2, Comment
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % bc, Comment
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % cc0, Comment
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % acc(0:5), Comment
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % accl, Comment	
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % bcc(0:5), Comment
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % bccl, Comment	
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % rhoc0, Comment
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % ac0, Comment	
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % kc0, Comment	
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % zc0 , Comment 
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % aac(0:5), Comment
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % aacl, Comment
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % ckc(0:5), Comment
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % ckcl, Comment
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % bzc(0:5), Comment
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % bzcl , Comment
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % nu20, Comment
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % anc(0:5), Comment
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % ancl , Comment
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % bnc(0:5), Comment
          READ (io, *, END=10, IOSTAT=OK, ERR=20) LocalSoluteMaterial % bncl, Comment
          DataRead=.TRUE.
10        CLOSE(io)
          IF (.NOT.DataRead) THEN
            WRITE(Message,'(A,A,A)')  'Not all entries in "Solute material File" ',TRIM(SoluteFileName),' found.'
            CALL FATAL(Trim(SubroutineName),Trim(message))
          END IF
        END IF
      END IF
      CurrentSoluteMaterial => LocalSoluteMaterial
    END IF
    CALL INFO(SubroutineName,"-----------------------------------------------------------------",Level=9)
    CALL INFO(SubroutineName,"Solute related constants",Level=9)
    WRITE(Message,*) "Mc",CurrentSoluteMaterial % Mc,"vc0",CurrentSoluteMaterial % vc0,&
         "kc0th", CurrentSoluteMaterial %kc0th
    CALL INFO(SubroutineName,Message,Level=9)    
    WRITE(Message,*) "d1",CurrentSoluteMaterial % d1,"d2",&         
         CurrentSoluteMaterial % d2,"bc",CurrentSoluteMaterial % bc
    CALL INFO(SubroutineName,Message,Level=9)
    WRITE(Message,*) "cc0",CurrentSoluteMaterial % cc0,"acc(0:5)",CurrentSoluteMaterial % acc(0:5),&
         "bcc(0:5)",CurrentSoluteMaterial % bcc(0:5)
    CALL INFO(SubroutineName,Message,Level=9)
    WRITE(Message,*) "rhoc0",CurrentSoluteMaterial % rhoc0,"ac0",CurrentSoluteMaterial % ac0,"kc0",&
         CurrentSoluteMaterial % kc0,"zc0",CurrentSoluteMaterial % zc0
    WRITE(Message,*)  "aac(0:5)",CurrentSoluteMaterial % aac(0:5),"ckc(0:5)",CurrentSoluteMaterial % ckc(0:5),&
         "bzc(0:5)",CurrentSoluteMaterial % bzc(0:5)
    CALL INFO(SubroutineName,Message,Level=9)
    WRITE(Message,*)  "bzc(0:5)",CurrentSoluteMaterial % bzc(0:5),"bnc(0:5)",CurrentSoluteMaterial % bnc(0:5),&
         "nu20",CurrentSoluteMaterial % nu20    
    CALL INFO(SubroutineName,Message,Level=9)
    WRITE(Message,*)  "aacl",LocalSoluteMaterial % aacl,"ckcl", LocalSoluteMaterial % ckcl    
    CALL INFO(SubroutineName,Message,Level=9)
    WRITE(Message,*)  "bzcl",LocalSoluteMaterial % bzcl,"ancl", LocalSoluteMaterial % ancl,&
         "bncl",LocalSoluteMaterial % bncl
    CALL INFO(SubroutineName,Message,Level=9)
    CALL INFO(SubroutineName,"-----------------------------------------------------------------",Level=9)
    RETURN
20  WRITE(Message,'(A,A,A)')  'Not all entries in "Solute material File" ',TRIM(SoluteFileName),' found.'
    CLOSE(io)
    CALL FATAL(Trim(SubroutineName),Trim(message))
  END SUBROUTINE ReadPermafrostSoluteMaterial

  !---------------------------------------------------------------------------------------------
  FUNCTION ReadPermafrostRockMaterial( Params ) RESULT(NumberOfRockRecords)
    IMPLICIT NONE
    TYPE(ValueList_t), POINTER :: Params
    INTEGER :: NumberOfRockRecords
    !--------------
    INTEGER :: I,J,K, n,t, RockMaterialID, active, DIM, ok,InitialNumberOfRockRecords, EntryNumber
    INTEGER,parameter :: io=21
    LOGICAL :: Found, fexist, FirstTime=.TRUE., AllocationsDone=.FALSE., DataRead=.FALSE.,&
         dimensionalTensorReduction = .FALSE.
    CHARACTER(LEN=MAX_NAME_LEN) ::  MaterialFileName, NewMaterialFileName, str, Comment
    CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName='ReadPermafrostRockMaterial'

    SAVE AllocationsDone,DataRead,InitialNumberOfRockRecords,MaterialFileName

    IF (DataRead) THEN
      NumberOfRockRecords = InitialNumberOfRockRecords
      RETURN
    ELSE ! we read Data from file database
      DIM = CoordinateSystemDimension()
      dimensionalTensorReduction = GetLogical( Params, 'Swap Tensor', Found)
      !------------------------------------------------------------------------------
      ! Inquire and open file
      !------------------------------------------------------------------------------
      ! give preference to a defined material database
      MaterialFileName = GetString( Params, 'Rock Material File', Found )
      IF (.NOT.Found) THEN
        CALL INFO(FunctionName," 'Rock Material File' keyword not found - looking for default DB!")
        fexist = .FALSE.
#ifdef USE_ISO_C_BINDINGS
        str = 'ELMER_LIB'
#else
        str = 'ELMER_LIB'//CHAR(0)
#endif
        CALL envir( str,MaterialFileName,k ) 
        IF ( k > 0  ) THEN
          MaterialFileName = MaterialFileName(1:k) // '/permafrostmaterialdb.dat'
          INQUIRE(FILE=TRIM(MaterialFileName), EXIST=fexist)
        END IF
        IF (.NOT. fexist) THEN
#ifdef USE_ISO_C_BINDINGS
          str = 'ELMER_HOME'
#else
          str = 'ELMER_HOME'//CHAR(0)
#endif
          CALL envir( str,MaterialFileName,k ) 
          IF ( k > 0 ) THEN
            MaterialFileName = MaterialFileName(1:k) // '/share/elmersolver/lib/' // 'permafrostmaterialdb.dat'
            INQUIRE(FILE=TRIM(MaterialFileName), EXIST=fexist)
          END IF
          IF ((.NOT. fexist) .AND. k>0) THEN
            MaterialFileName = MaterialFileName(1:k) // '/permafrostmaterialdb.dat'
            INQUIRE(FILE=TRIM(MaterialFileName), EXIST=fexist)
          END IF
        END IF
        IF (.NOT. fexist) THEN
          CALL Fatal('CheckKeyWord', 'permafrostmaterialdb.dat not found')
        END IF
      END IF

      ! if we are still here, we open the file (what ever it may be)
      OPEN(unit = io, file = TRIM(MaterialFileName), status = 'old',iostat = ok)
      IF (ok /= 0) THEN
        WRITE(Message,'(A,A)') 'Unable to open file ',TRIM(MaterialFileName)
        CALL FATAL(Trim(FunctionName),Trim(message))
      ELSE
        !------------------------------------------------------------------------------
        ! Read in the number of records in file (first line integer)
        !------------------------------------------------------------------------------
        READ (io, *, END=30, IOSTAT=OK, ERR=40) NumberOfRockRecords, Comment
        WRITE (Message,*) "Attempting to read ",NumberOfRockRecords," ",&
             TRIM(Comment)," records from data file ",TRIM(MaterialFileName)        
        CALL INFO(FunctionName,Message,level=3)
        InitialNumberOfRockRecords = NumberOfRockRecords
      END IF
      !------------------------------------------------------------------------------
      ! Allocate and read stuff
      !------------------------------------------------------------------------------
      !M = Model % Mesh % NumberOfNodes
      IF (AllocationsDone) THEN
        DEALLOCATE(&
             GlobalRockMaterial % ks0th,&
             GlobalRockMaterial % e1,&
             GlobalRockMaterial % bs,&
             GlobalRockMaterial % rhos0,&
             GlobalRockMaterial % cs0,&
             GlobalRockMaterial % Xi0,&
             GlobalRockMaterial % eta0,&
             GlobalRockMaterial % etak,&
             GlobalRockMaterial % hs0,&
             GlobalRockMaterial % Kgwh0, &
             GlobalRockMaterial % qexp, &
             GlobalRockMaterial % alphaL, &
             GlobalRockMaterial % alphaT, &
             GlobalRockMaterial % RadGen, &
             GlobalRockMaterial % acs, &
             GlobalRockMaterial % as0, &
             GlobalRockMaterial % aas, &
             GlobalRockMaterial % ks0, &
             GlobalRockMaterial % cks, &
             GlobalRockMaterial % Es0, &
             GlobalRockMaterial % nuS0, &
             GlobalRockMaterial % acsl, &
             GlobalRockMaterial % aasl, &
             GlobalRockMaterial % cksl, &
             GlobalRockMaterial % VariableBaseName)
      END IF
      ALLOCATE(&
           GlobalRockMaterial % ks0th(NumberOfRockRecords),&
           GlobalRockMaterial % e1(NumberOfRockRecords),&
           GlobalRockMaterial % bs(NumberOfRockRecords),&
           GlobalRockMaterial % rhos0(NumberOfRockRecords),&
           GlobalRockMaterial % cs0(NumberOfRockRecords),&
           GlobalRockMaterial % Xi0(NumberOfRockRecords),&
           GlobalRockMaterial % eta0(NumberOfRockRecords),&
           GlobalRockMaterial % etak(NumberOfRockRecords),&
           GlobalRockMaterial % hs0(NumberOfRockRecords),&
           GlobalRockMaterial % Kgwh0(3,3,NumberOfRockRecords),&
           GlobalRockMaterial % qexp(NumberOfRockRecords), &
           GlobalRockMaterial % alphaL(NumberOfRockRecords), &
           GlobalRockMaterial % alphaT(NumberOfRockRecords), &
           GlobalRockMaterial % RadGen(NumberOfRockRecords), &
           GlobalRockMaterial % acs(0:5,NumberOfRockRecords), &
           GlobalRockMaterial % as0(NumberOfRockRecords), &
           GlobalRockMaterial % aas(0:5,NumberOfRockRecords), &
           GlobalRockMaterial % ks0(NumberOfRockRecords), &
           GlobalRockMaterial % cks(0:5,NumberOfRockRecords), &
           GlobalRockMaterial % Es0(NumberOfRockRecords), &
           GlobalRockMaterial % nuS0(NumberOfRockRecords), &
           GlobalRockMaterial % acsl(NumberOfRockRecords), &     
           GlobalRockMaterial % aasl(NumberOfRockRecords), &
           GlobalRockMaterial % cksl(NumberOfRockRecords), &
           GlobalRockMaterial % VariableBaseName(NumberOfRockRecords),&
           STAT=OK)
      AllocationsDone = .TRUE.
      DataRead = .TRUE.
      IF (OK /= 0) THEN
        CLOSE(io)
        CALL FATAL(FunctionName, 'Allocation Error of input data array')
      END IF
      
      DO RockMaterialID=1,NumberOfRockRecords
        READ (io, *, END=30, IOSTAT=OK, ERR=40) GlobalRockMaterial % VariableBaseName(RockMaterialID), EntryNumber
        IF (EntryNumber /= RockMaterialID) THEN
          WRITE(Message,'(A,I3,A,I3)') &
               "Entry number", EntryNumber, "does not match expected number ",I
          CLOSE(io)
          CALL FATAL(FunctionName,Message)
        ELSE
          WRITE(Message,'(A,A,A,I3,A)')&
               "Material ", TRIM(GlobalRockMaterial % VariableBaseName(RockMaterialID)),&
               " entry number ", EntryNumber, " will be read in"
          CALL INFO(FunctionName,Message,Level=3)
        END IF
        WRITE(Message,'(A,I2,A,A)') "Input for Variable No.",RockMaterialID,": ",&
             GlobalRockMaterial % VariableBaseName(RockMaterialID)
        CALL INFO(FunctionName,Message,Level=9)
        READ (io, *, END=30, IOSTAT=OK, ERR=40) GlobalRockMaterial % Xi0(RockMaterialID), Comment
        READ (io, *, END=30, IOSTAT=OK, ERR=40) GlobalRockMaterial % eta0(RockMaterialID), Comment
        READ (io, *, END=30, IOSTAT=OK, ERR=40) GlobalRockMaterial % etak(RockMaterialID), Comment
        READ (io, *, END=30, IOSTAT=OK, ERR=40) GlobalRockMaterial % ks0th(RockMaterialID), Comment
        READ (io, *, END=30, IOSTAT=OK, ERR=40) GlobalRockMaterial % e1(RockMaterialID), Comment
        READ (io, *, END=30, IOSTAT=OK, ERR=40) GlobalRockMaterial % bs(RockMaterialID), Comment
        READ (io, *, END=30, IOSTAT=OK, ERR=40) GlobalRockMaterial % rhos0(RockMaterialID), Comment
        READ (io, *, END=30, IOSTAT=OK, ERR=40) GlobalRockMaterial % cs0(RockMaterialID), Comment
        READ (io, *, END=30, IOSTAT=OK, ERR=40) GlobalRockMaterial % hs0(RockMaterialID), Comment
        DO I=1,3
          DO J=1,3
            READ (io, *, END=30, IOSTAT=OK, ERR=40) GlobalRockMaterial % Kgwh0(I,J,RockMaterialID), Comment
          END DO
        END DO
        IF ((DIM == 2) .AND. dimensionalTensorReduction) THEN
          WRITE(Message,'(A,I2,A,A)') "Reducing tensor Kgwh0 for material no.",RockMaterialID
          CALL INFO(FunctionName,Message,Level=9)
          GlobalRockMaterial % Kgwh0(2,2,RockMaterialID) = GlobalRockMaterial % Kgwh0(3,3,RockMaterialID)
          GlobalRockMaterial % Kgwh0(1,2,RockMaterialID) = GlobalRockMaterial % Kgwh0(1,3,RockMaterialID)
          GlobalRockMaterial % Kgwh0(2,1,RockMaterialID) = GlobalRockMaterial % Kgwh0(3,1,RockMaterialID)
          GlobalRockMaterial % Kgwh0(3,1:3,RockMaterialID) = 0.0_dp
          GlobalRockMaterial % Kgwh0(1:3,3,RockMaterialID) = 0.0_dp
        END IF
        READ (io, *, END=30, IOSTAT=OK, ERR=40) GlobalRockMaterial % qexp(RockMaterialID), Comment
        READ (io, *, END=30, IOSTAT=OK, ERR=40) GlobalRockMaterial % alphaL(RockMaterialID), Comment
        READ (io, *, END=30, IOSTAT=OK, ERR=40) GlobalRockMaterial % alphaT(RockMaterialID), Comment
        READ (io, *, END=30, IOSTAT=OK, ERR=40) GlobalRockMaterial % RadGen(RockMaterialID), Comment
        READ (io, *, END=30, IOSTAT=OK, ERR=40) GlobalRockMaterial % acs(0:5,RockMaterialID),  Comment
        READ (io, *, END=30, IOSTAT=OK, ERR=40) GlobalRockMaterial % as0(RockMaterialID),  Comment
        READ (io, *, END=30, IOSTAT=OK, ERR=40) GlobalRockMaterial % aas(0:5,RockMaterialID),  Comment
        READ (io, *, END=30, IOSTAT=OK, ERR=40) GlobalRockMaterial % ks0(RockMaterialID),  Comment
        !--------------------
        READ (io, *, END=30, IOSTAT=OK, ERR=40) GlobalRockMaterial % cks(0:5,RockMaterialID),  Comment
        READ (io, *, END=30, IOSTAT=OK, ERR=40) GlobalRockMaterial % Es0(RockMaterialID), Comment
        READ (io, *, END=30, IOSTAT=OK, ERR=40) GlobalRockMaterial % nuS0(RockMaterialID), Comment
        READ (io, *, END=30, IOSTAT=OK, ERR=40) GlobalRockMaterial % acsl(RockMaterialID),  Comment
        READ (io, *, END=30, IOSTAT=OK, ERR=40) GlobalRockMaterial % aasl(RockMaterialID),  Comment
        READ (io, *, END=30, IOSTAT=OK, ERR=40) GlobalRockMaterial % cksl(RockMaterialID),  Comment
      END DO
      !CALL INFO(FunctionName,Message,Level=1)
30    CLOSE(io)
      IF (RockMaterialID < NumberOfRockRecords) THEN
        WRITE(Message,'(I3,A,I3)') RockMaterialID,"records read, which is smaller than given number ", NumberOfRockRecords
        CALL FATAL(FunctionName,Message)
      ELSE
        WRITE(Message,'(A,I2,A,A)') "Read ",NumberOfRockRecords," rock material records from file ", TRIM(MaterialFileName)
        CALL INFO(FunctionName,Message,Level=1)
      END IF
      RETURN
    END IF

40  CALL WARN(FunctionName,"I/O error! Last successfully read variable:")
    CALL WARN(FunctionName,Comment)
    CALL FATAL(FunctionName,"Stopping simulation")    
  END FUNCTION ReadPermafrostRockMaterial
  
  !---------------------------------------------------------------------------------------------  
  FUNCTION ReadPermafrostElementRockMaterial(MaterialFileName,Solver,DIM,SkipInit) RESULT(NumberOfRockRecords)
    IMPLICIT NONE
    CHARACTER(LEN=MAX_NAME_LEN), INTENT(IN) :: MaterialFileName
    TYPE(Solver_t) :: Solver
    INTEGER :: NumberOfRockRecords,DIM
    LOGICAL, OPTIONAL :: SkipInit
    !-----------------------------------------------------------
    CHARACTER(LEN=MAX_NAME_LEN) :: SubroutineName="ReadPermafrostElementRockMaterial"
    LOGICAL :: FirstTime=.TRUE., Parallel=.FALSE.
    TYPE(RockMaterial_t), TARGET :: LocalRockMaterial
    TYPE(Element_t), POINTER :: CurrentElement
    INTEGER, ALLOCATABLE :: GlobalToLocalPerm(:)
    INTEGER :: OK, CurrentNo, I, J, NoElements, LocalNoElements, &
         minglobalelementnumber, maxglobalelementnumber, mmaxglobalelementnumber, ierr
    INTEGER, PARAMETER :: io=21
    REAL(KIND=dp) :: ReceivingArray(50)    

    SAVE LocalRockMaterial, FirstTime, Parallel, minglobalelementnumber, maxglobalelementnumber,&
         GlobalToLocalPerm, NoElements, mmaxglobalelementnumber

    IF (PRESENT(SkipInit) .AND. FirstTime) CALL FATAL(SubroutineName,'Initialization error')
    
    IF (FirstTime) THEN
      NoElements = Solver % NumberOfActiveElements ! active elements in partition/serial mesh
      Parallel = ( ParEnv % PEs > 1 )
      IF ( Parallel ) THEN        
        DO I=1,NoElements
          CurrentElement => Solver % Mesh % Elements(I)
          IF (FirstTime) THEN
            minglobalelementnumber = CurrentElement % GElementIndex
            maxglobalelementnumber = minglobalelementnumber
            FirstTime = .FALSE.
          ELSE
            minglobalelementnumber = MIN((CurrentElement % GElementIndex),minglobalelementnumber)
            maxglobalelementnumber = MAX((CurrentElement % GElementIndex),maxglobalelementnumber)
          END IF
        END DO
        CALL MPI_ALLREDUCE(maxglobalelementnumber,mmaxglobalelementnumber,1,&
            MPI_INTEGER,MPI_MAX,ELMER_COMM_WORLD,ierr)
        
        !IF (ParEnv % myPE == 0) &
        !     PRINT *,"ReadPermafrostElementRockMaterial:",Parenv % myPE, "ming/maxg",&
        !     minglobalelementnumber,maxglobalelementnumber
        IF ((maxglobalelementnumber - minglobalelementnumber) < 1) &
             CALL FATAL("ReadPermafrostElementRockMaterial","Failed to create global to local permutation")
        ALLOCATE(GlobalToLocalPerm(maxglobalelementnumber - minglobalelementnumber + 1), STAT=OK)
        IF (OK /= 0) CALL FATAL("ReadPermafrostElementRockMaterial","Allocation error of GlobalToLocalPerm")
        GlobalToLocalPerm = 0
        DO I=1,NoElements
          CurrentElement => Solver % Mesh % Elements(I)
          GlobalToLocalPerm((CurrentElement % GElementIndex) - minglobalelementnumber + 1) = I
          !IF (ParEnv % myPE == 0) &
          !     PRINT *,"ReadPermafrostElementRockMaterial:",Parenv % myPE, &
          !     "GlobalToLocalPerm(",(CurrentElement % GElementIndex)," - ",minglobalelementnumber," + 1)=",I
        END DO
      ELSE
        minglobalelementnumber = 1
        maxglobalelementnumber = NoElements
        mmaxglobalelementnumber = NoElements
      END IF      
      ALLOCATE(&           
           GlobalRockMaterial % ks0th(NoElements),&
           GlobalRockMaterial % e1(NoElements),&
           GlobalRockMaterial % bs(NoElements),&
           GlobalRockMaterial % rhos0(NoElements),&
           GlobalRockMaterial % Xi0(NoElements),&
           GlobalRockMaterial % eta0(NoElements),&
           GlobalRockMaterial % etak(NoElements),&
           GlobalRockMaterial % hs0(NoElements),&
           GlobalRockMaterial % Kgwh0(3,3,NoElements),&
           GlobalRockMaterial % qexp(NoElements), &
           GlobalRockMaterial % alphaL(NoElements), &
           GlobalRockMaterial % alphaT(NoElements), &
           GlobalRockMaterial % RadGen(NoElements), &
           GlobalRockMaterial % cs0(NoElements),&
           GlobalRockMaterial % acs(0:5,NoElements), &
           GlobalRockMaterial % as0(NoElements), &
           GlobalRockMaterial % aas(0:5,NoElements), &
           GlobalRockMaterial % ks0(NoElements), &
           GlobalRockMaterial % cks(0:5,NoElements), &
           GlobalRockMaterial % Es0(NoElements),&
           GlobalRockMaterial % nus0(NoElements),&
           GlobalRockMaterial % acsl(NoElements), &
           GlobalRockMaterial % aasl(NoElements), &
           GlobalRockMaterial % cksl(NoElements), &
           GlobalRockMaterial % VariableBaseName(NoElements),&
           STAT=OK)
      OPEN(unit = io, file = TRIM(MaterialFileName), status = 'old',iostat = ok)
      IF (ok /= 0) THEN
        WRITE(Message,'(A,A)') 'Unable to open file ',TRIM(MaterialFileName)
        CALL FATAL(Trim(SubroutineName),Trim(message))        
      ELSE        
        !------------------------------------------------------------------------------
        ! Read in the number of records in file (first line integer)
        ! MIND: all receiving array numbers are shifted by -1 in index with resepect
        !       to J. Hartikainen's instructions in input_data_forsmark_2d_example.pdf!
        !------------------------------------------------------------------------------
        WRITE (Message,*) "Attempting to read ",mmaxglobalelementnumber,&
             " records from data file ",TRIM(MaterialFileName)
        CALL INFO(SubroutineName,Message,level=3)
        LocalNoElements = 0
        DO J=1,maxglobalelementnumber                             
          READ (io, *, END=50, ERR=60, IOSTAT=OK) CurrentNo, ReceivingArray(1:50)
          IF ( Parallel ) THEN
            IF (J < minglobalelementnumber) CYCLE
            I = GlobalToLocalPerm(J - minglobalelementnumber +1)
            !IF (I> 0) &
            !     PRINT *,"ReadPermafrostElementRockMaterial:", Parenv % myPE, &
            !     "GlobalToLocalPerm(",J," -", minglobalelementnumber," +1) =", &
            !     GlobalToLocalPerm(J - minglobalelementnumber +1)
            IF (I == 0) CYCLE
            LocalNoElements = LocalNoElements + 1
          ELSE
            I=J
            LocalNoElements = LocalNoElements + 1
          END IF
          CurrentElement => Solver % Mesh % Elements(I)
          !! IMPORTANT: Mind that all ReceivingArray numbers ar N-1 with respect to the document (input_data_forsmark_2d)
          GlobalRockMaterial % ks0th(I) = ReceivingArray(11) ! shall be changed to tensor
          !-----------------------------
          GlobalRockMaterial % e1(I) = ReceivingArray(33) ! e1 (mail from Juha 11.10.)
          !IF (GlobalRockMaterial % e1(I) > 0.01) PRINT *,"e1:", ReceivingArray(34)
          !IF (GlobalRockMaterial % e1(I) < 0.0) PRINT *,"e1:", ReceivingArray(34)
          GlobalRockMaterial % bs(I) = ReceivingArray(23) ! b11,1 (mail from Juha 11.10.)
          GlobalRockMaterial % rhos0(I) = ReceivingArray(1)
          GlobalRockMaterial % Xi0(I) = ReceivingArray(32)
          !-----------------------------
          GlobalRockMaterial % eta0(I) = ReceivingArray(30) ! eta_t (mail from Juha 11.10.)
          GlobalRockMaterial % etak(I) = ReceivingArray(31)
          GlobalRockMaterial % hs0(I) = 0.0_dp! will be removed
          !----------------------------- Hydrol. Conductivity
          !GlobalRockMaterial % Kgwh0 = 0.0_dp

          IF(DIM==2) THEN
            GlobalRockMaterial % Kgwh0(1,1,I) = ReceivingArray(35)
            GlobalRockMaterial % Kgwh0(2,2,I) = ReceivingArray(37)
            GlobalRockMaterial % Kgwh0(1,2,I) = ReceivingArray(39)
            GlobalRockMaterial % Kgwh0(2,1,I) = ReceivingArray(39)
            GlobalRockMaterial % Kgwh0(3,1:3,I) = 0.0_dp
            GlobalRockMaterial % Kgwh0(1:3,3,I) = 0.0_dp
          ELSE
            GlobalRockMaterial % Kgwh0(1,1,I) = ReceivingArray(35)
            GlobalRockMaterial % Kgwh0(2,2,I) = ReceivingArray(36)
            GlobalRockMaterial % Kgwh0(3,3,I) = ReceivingArray(37)
            GlobalRockMaterial % Kgwh0(1,2,I) = ReceivingArray(38)
            GlobalRockMaterial % Kgwh0(1,3,I) = ReceivingArray(39)
            GlobalRockMaterial % Kgwh0(2,3,I) = ReceivingArray(40)
            GlobalRockMaterial % Kgwh0(2,1,I) = ReceivingArray(38)
            GlobalRockMaterial % Kgwh0(3,1,I) = ReceivingArray(39)
            GlobalRockMaterial % Kgwh0(3,2,I) = ReceivingArray(40)
          END IF
          !-----------------------------
          GlobalRockMaterial % qexp(I) = ReceivingArray(41) !?????????????????????????????????????????????
          GlobalRockMaterial % alphaL(I) = ReceivingArray(47)
          GlobalRockMaterial % alphaT(I) = ReceivingArray(48)
          GlobalRockMaterial % RadGen(I) = ReceivingArray(29)
          !-----------------------------
          GlobalRockMaterial % cs0(I) = ReceivingArray(8)
          GlobalRockMaterial % acs(0,I) =  ReceivingArray(9)
          GlobalRockMaterial % acs(1,I) =  ReceivingArray(10)
          GlobalRockMaterial % acs(2:5,I) = 0.0_dp
          GlobalRockMaterial % acsl(I)= 1
          !-----------------------------
          GlobalRockMaterial % as0(I)= ReceivingArray(2)
          GlobalRockMaterial % aas(0,I) =  ReceivingArray(3)
          GlobalRockMaterial % aas(1,I) =  ReceivingArray(4)
          GlobalRockMaterial % aas(2:5,I) = 0.0_dp
          GlobalRockMaterial % aasl(I)= 1
          !-----------------------------
          GlobalRockMaterial % ks0(I)= ReceivingArray(5)
          GlobalRockMaterial % cks(0,I) = ReceivingArray(6)
          GlobalRockMaterial % cks(1,I) = ReceivingArray(7)
          GlobalRockMaterial % cks(2:5,I)= 0.0_dp
          GlobalRockMaterial % cksl(I)= 1
          !-----------------------------
          GlobalRockMaterial % Es0(I) = ReceivingArray(49)
          GlobalRockMaterial % nus0(I) = ReceivingArray(50)
          WRITE(Message,*) 'Element',I
          GlobalRockMaterial % VariableBaseName(I) = TRIM(Message)
        END DO
        GlobalRockMaterial % NumberOfRockRecords = NoElements
        NumberOfRockRecords = NoElements
50      CLOSE(io)
        IF (LocalNoElements < NoElements) THEN
          IF (Parallel) THEN
            WRITE (Message,*) 'Parallel proc.', ParEnv % myPe, '. Found ONLY ',&
                 LocalNoElements,' entries in file ',TRIM(MaterialFileName),&
                 ' for ', NoElements, ' elements in mesh partition.'
          ELSE
            WRITE (Message,*) 'Found ONLY ',LocalNoElements,' entries in file ',&
                 TRIM(MaterialFileName),&
                 ' for ', NoElements, ' elements in mesh'
          END IF
          CALL FATAL(TRIM(SubroutineName),Message)
        ELSE
          IF (Parallel) THEN
            PRINT *, TRIM(SubroutineName), ': Parallel proc.', ParEnv % myPe, '. Read ', NoElements, 'entries from file'
          ELSE
            WRITE (Message,*) 'Read ',LocalNoElements,' entries from file ',TRIM(MaterialFileName)
            CALL INFO(TRIM(SubroutineName),Message,Level=3)
          END IF
        END IF
      END IF
      IF (Parallel) DEALLOCATE(GlobalToLocalPerm)
      CALL MPI_BARRIER(ELMER_COMM_WORLD,ierr)
      FirstTime = .FALSE.
    ELSE
            NumberOfRockRecords = NoElements
    END IF
    RETURN
60  WRITE (Message,*) 'I/O error at entry ',CurrentNo,' of file ',TRIM(MaterialFileName)
    CALL FATAL(TRIM(SubroutineName),Message)
  END FUNCTION ReadPermafrostElementRockMaterial
  
  !---------------------------------------------------------------------------------------------
  FUNCTION ReadPermafrostConstants(Model, FunctionName,&
       DIM, GasConstant, N0, DeltaT, T0, p0, eps, Gravity) RESULT(Constantsread)
    !------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    CHARACTER(LEN=MAX_NAME_LEN) :: FunctionName
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    INTEGER :: DIM
    REAL(KIND=dp) :: GasConstant, N0, DeltaT, T0, p0,eps, Gravity(3) ! constants read only once 
    LOGICAL :: Constantsread
    !------------------------------------------------------------------------------
    REAL(KIND=dp), POINTER :: gWork(:,:)
    LOGICAL :: Found
    INTEGER :: I
    !------------------------------------------------------------------------------
    DIM = CoordinateSystemDimension()
    gWork => ListGetConstRealArray( Model % Constants,"Gravity",Found)
    IF (.NOT.Found) THEN
      Gravity = 0.0
      CALL WARN(FunctionName,'Gravity not found in Constants section. Setting to zero')
    ELSE
      Gravity = gWork(1:3,1)*gWork(4,1)
    END IF
    !------------------------------------------------------------------------------
    ! Constants
    ! GasConstant, N0, T0, p0, DeltaT, eps
    !------------------------------------------------------------------------------
    !IF (.NOT.ASSOCIATED(Model % Constants)) STOP
    GasConstant = GetConstReal(Model % Constants, "Gas Constant", Found)
    IF (.NOT.Found) THEN
      GasConstant = 8.3145_dp
      CALL INFO(FunctionName, ' "Gas Constant" not found in Constants and set to default value 8.3145',Level=3)
    END IF
    N0 = GetConstReal(Model % Constants, "Avogadro Number", Found)
    IF (.NOT.Found) THEN
      N0 = 6.022140857d23
      CALL INFO(FunctionName, ' "Avogadro Number" not found in Constants and set to default value 6.022140857E23',Level=3)
    END IF
    T0 = GetConstReal(Model % Constants, 'Reference Temperature', Found)
    IF (.NOT.Found) THEN
      T0 = 273.15_dp
      CALL INFO(FunctionName, ' "Reference Temperature" not found in Constants and set to default value T0=273.15',Level=3)
    END IF
    p0 = GetConstReal(Model % Constants, 'Reference Pressure', Found)
    IF (.NOT.Found) THEN
      p0 = 100132.0_dp
      CALL INFO(FunctionName, ' "Reference Pressure not found in Constants and set to default value p0=100132.0',Level=3)
    END IF
    DeltaT = GetConstReal(Model % Constants,"Permafrost DeltaT",Found)
    IF (.NOT.Found) THEN
      DeltaT = 1.0_dp
      CALL INFO(FunctionName, ' "Permafrost DeltaT" not found in Constants and set to default value DeltaT=1.0',Level=3)
    END IF
    Eps = GetConstReal(Model % Constants,"Permafrost eps",Found)
    IF (.NOT.Found) THEN
      eps = 0.99_dp
      CALL INFO(FunctionName, ' "Permafrost eps" not found in Constants and set to default value eps=0.99',Level=3)
    END IF
    ConstantsRead = .TRUE.
    CALL INFO(FunctionName,"-----------------------------------------------------------------",Level=9)
    CALL INFO(FunctionName,"Model Constants:", Level=9)
    WRITE(Message,*) "GasConstant, T0, p0, DeltaT, eps:"
    CALL INFO(FunctionName,Message, Level=9)
    WRITE(Message,*) GasConstant, T0, p0, DeltaT, eps
    CALL INFO(FunctionName,Message, Level=9)
    CALL INFO(FunctionName,"-----------------------------------------------------------------",Level=9)
  END FUNCTION ReadPermafrostConstants
  !---------------------------------------------------------------------------------------------
  ! assign single nodal variable
  !---------------------------------------------------------------------------------------------
  SUBROUTINE AssignSingleVar(Solver,Model,NodalVariable,VariableVar,VariablePerm,Variable,&
       VariableName,VariableDOFS,VariableExists,PrevNodalVariable, PrevVariable)
    IMPLICIT NONE
    
    TYPE(Solver_t) :: Solver
    TYPE(Model_t) :: Model
    REAL(KIND=dp),POINTER :: NodalVariable(:), Variable(:)
    LOGICAL :: VariableExists
    CHARACTER(LEN=MAX_NAME_LEN) :: VariableName
    TYPE(Variable_t), POINTER :: VariableVar
    INTEGER, POINTER :: VariablePerm(:)
    INTEGER :: VariableDOFS
    REAL(KIND=dp),POINTER,OPTIONAL :: PrevNodalVariable(:), PrevVariable(:)

    ! ----
    INTEGER :: N, istat
    CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='AssignSingleVar'
    
    
    !IF (.NOT.ASSOCIATED(VariableVar)) &
    VariableVar => VariableGet(Solver % Mesh % Variables,VariableName)
    IF (.NOT.ASSOCIATED(VariableVar)) THEN
      VariableExists = .FALSE.
      WRITE (Message,*) 'Variable ',TRIM(VariableName),' not found'
      CALL WARN(SolverName,Message) 
      RETURN
    ELSE
      VariableDOFS = VariableVar % DOFs
      VariablePerm => VariableVar % Perm
      Variable => VariableVar % Values
      IF (.NOT.ASSOCIATED(VariablePerm) .OR. .NOT.ASSOCIATED(Variable)) &
           CALL FATAL(SolverName, ' Error in assignments of variable pointers')
      !PRINT *, PRESENT(PrevVariable)
      IF (PRESENT(PrevVariable)) THEN
        PrevVariable => VariableVar % PrevValues(:,1)
        CALL INFO(SolverName,"Assigned previous time values variable pointer",Level=1)
      END IF
    END IF
    
    IF ((.NOT.VariableExists) .OR. (Model % Mesh % Changed)) THEN
      N = MAX( Solver % Mesh % MaxElementDOFs, Solver % Mesh % MaxElementNodes )
      IF (VariableExists) THEN
        DEALLOCATE(NodalVariable)
        IF (PRESENT(PrevNodalVariable)) DEALLOCATE(PrevNodalVariable)
      END IF
      ALLOCATE(NodalVariable(N*VariableDOFS),STAT=istat )      
      IF (PRESENT(PrevNodalVariable) .AND. (istat == 0))  ALLOCATE(PrevNodalVariable(N*VariableDOFS),STAT=istat )
      IF ( istat /= 0 ) THEN
        CALL FATAL(SolverName,"Allocation error")
      ELSE
        VariableExists = .TRUE.
	WRITE(Message,*) "Allocations done for nodal variable of ",TRIM(VariableName)
        CALL INFO(SolverName,Message,Level=1)
      END IF
    END IF    

  END SUBROUTINE AssignSingleVar
  !---------------------------------------------------------------------------------------------
  ! assign single nodal variable time derivative
  !---------------------------------------------------------------------------------------------
  SUBROUTINE AssignSingleVarTimeDer(Solver,Model,Element,NodalVariableTimeDer,&
       VariableVar,VariableTimeDerExists,dt)
    IMPLICIT NONE
    
    TYPE(Solver_t) :: Solver
    TYPE(Model_t) :: Model
    TYPE(Element_t) :: Element
    REAL(KIND=dp),POINTER :: NodalVariableTimeDer(:)
    REAL(KIND=dp) :: dt
    LOGICAL :: VariableTimeDerExists
    TYPE(Variable_t), POINTER :: VariableVar
    ! ----
    INTEGER :: VariableDOFS
    !LOGICAL :: AllocationsDone=.FALSE.
    INTEGER :: MaxNodes, I, J, istat,CurrentvariableNodeIndex
    INTEGER, POINTER :: VariablePerm(:)
    CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='AssignSingleVarTimeDer'
    REAL(KIND=dp), POINTER :: Variable(:),VariablePrev(:,:)
    
    SAVE MaxNodes

    IF (dt <= 0.0_dp) CALL FATAL(SolverName, "Negative or zero timestep")
    
    IF ((.NOT.VariableTimeDerExists) .OR. (Model % Mesh % Changed)) THEN
      MaxNodes = MAX( Solver % Mesh % MaxElementDOFs, Solver % Mesh % MaxElementNodes )
      IF (VariableTimeDerExists) THEN
        CALL INFO(SolverName,"Deallocation of nodal time derivtive")
        DEALLOCATE(NodalVariableTimeDer)
      END IF

      VariableDOFS = VariableVar % DOFs
      ALLOCATE(NodalVariableTimeDer(MaxNodes*VariableDOFS),STAT=istat )
      IF ( istat /= 0 ) THEN
        CALL FATAL(SolverName,"Allocation error")
      ELSE
        VariableTimeDerExists = .TRUE.
        CALL INFO(SolverName,"Allocations Done",Level=1)
      END IF
    END IF
    
    VariableDOFS = VariableVar % DOFs
    VariablePrev => VariableVar % PrevValues
    NodalVariableTimeDer(1:maxNodes*VariableDOFS) = 0.0_dp
    IF ( gettimestep() == 1 ) RETURN ! use zero value in 1st timestep
    
    IF (.NOT.ASSOCIATED(VariablePrev)) THEN
      VariableTimeDerExists = .FALSE.
     ELSE
      VariablePerm => VariableVar % Perm
      VariableDOFS = VariableVar % DOFs
      Variable => VariableVar % Values
      VariableTimeDerExists = .TRUE.
      IF (MaxNodes < GetElementNOFNodes(Element)) CALL FATAL(SolverName,"Number of Nodes exceeds allocation")
      DO I=1,GetElementNOFNodes(Element)
        CurrentVariableNodeIndex = VariablePerm(Element % NodeIndexes(I))
        DO J=1,VariableDOFS
          NodalVariableTimeDer((I-1)*VariableDOFS + J) = &
               (Variable((CurrentVariableNodeIndex - 1) * VariableDOFS + J) &
               - VariablePrev((CurrentVariableNodeIndex - 1) * VariableDOFS + J,1))/dt          
        END DO
        !PRINT *,"AssignSingleVarTimeDer:", Variable((CurrentVariableNodeIndex - 1) * VariableDOFS + J),&
        !     VariablePrev((CurrentVariableNodeIndex - 1) * VariableDOFS + J,1)
      END DO
    END IF

  END SUBROUTINE AssignSingleVarTimeDer

  ! compute element-wise single nodal variable
  SUBROUTINE ReadSingleVar(N,Element,VariablePerm,NodalVariable,Variable,VariableDOFs)
    IMPLICIT NONE
    
    INTEGER :: N,VariableDOFs
    INTEGER, POINTER :: VariablePerm(:)
    TYPE(Element_t) :: Element
    REAL(KIND=dp),POINTER :: NodalVariable(:),Variable(:)
    !-----------------------
    INTEGER :: I,J

    DO I=1,N
      DO J=1,VariableDOFs
        NodalVariable((VariableDOFs*I - 1) + J) = &
             Variable(VariableDOFs*(VariablePerm(Element % NodeIndexes(I))-1) + J)
      END DO
    END DO
  END SUBROUTINE ReadSingleVar

  !---------------------------------------------------------------------------------------------
  !---------------------------------------------------------------------------------------------
  !---------------------------------------------------------------------------------------------
  !---------------------------
  ! model parameters functions
  !---------------------------
  !---------------------------------------------------------------------------------------------
  !---------------------------------------------------------------------------------------------
  !
  !---------------------------------------------------------------------------------------------
  ! general functions 
  !---------------------------------------------------------------------------------------------
  REAL(KIND=dp) FUNCTION GeneralPolynomial(Variable,ReferenceValue,Normation,coeff,pdeg)
    IMPLICIT NONE
    !-------
    REAL(KIND=dp), INTENT(IN) :: Variable,ReferenceValue,Normation,coeff(0:5)
    INTEGER, INTENT(IN) :: pdeg
    REAL(KIND=dp) outval
    ! ------
    REAL(KIND=dp) currpot
    INTEGER :: i

    outval = 0.0_dp
    currpot = 1.0_dp
    DO i=0,pdeg
      outval = outval + coeff(i) * currpot
      currpot = currpot * (Variable - ReferenceValue)/Normation
    END DO
    GeneralPolynomial = outval
  END FUNCTION GeneralPolynomial
  !---------------------------------------------------------------------------------------------
  REAL(KIND=dp) FUNCTION GeneralIntegral(Variable,ReferenceValue,Normation,coeff0,coeff,pdeg)
    IMPLICIT NONE
    !-------
    REAL(KIND=dp), INTENT(IN) :: Variable,ReferenceValue,Normation,coeff0,coeff(0:5)
    INTEGER, INTENT(IN) :: pdeg
    REAL(KIND=dp) prefactor, summation
    ! ------
    REAL(KIND=dp) currpot
    INTEGER :: currdeg
    prefactor = coeff0*(Variable - ReferenceValue)
    summation = 0.0_dp
    DO currdeg=0,pdeg
      summation = summation +&
           ( coeff(currdeg) * &
           ( (Variable - ReferenceValue)/Normation )**(DBLE(currdeg)) )/(DBLE(currdeg)  + 1.0_dp)
    END DO
    !currpot = 1.0_dp
    !DO currdeg=0,pdeg
    !  outval = outval * coeff(currdeg) * currpot/(DBLE(currdeg) + 1.0_dp)
    !  currpot = currpot * (Variable - ReferenceValue)/Normation
    !END DO
    GeneralIntegral = prefactor * summation
  END FUNCTION GeneralIntegral
  !---------------------------------------------------------------------------------------------
  ! functions specific to heat transfer and phase change
  !---------------------------------------------------------------------------------------------
  FUNCTION GetXiAnderson(A,B,Beta,rhow,rhos0,T0,Temperature,Pressure,Porosity) RESULT(XiAnderson)
    REAL(KIND=dp), INTENT(IN) :: A,B,Beta,rhow,rhos0,T0,Temperature,Pressure,Porosity
    REAL(KIND=dp) :: Tstar, XiAnderson
    IF (Porosity <= 0.0) &
         CALL FATAL("PermafrostMaterials(GetXiAnderson)","Zero or negative porosity detected")
    Tstar = T0 - Beta * Pressure - Temperature
    XiAnderson =  MAX(MIN((rhos0/rhow)*(A*(Tstar**B)/Porosity),1.0_dp),0.0_dp)
  END FUNCTION GetXiAnderson
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION XiAndersonT(Xi,A,B,Beta,rhow,rhos0,T0,Temperature,Pressure,Porosity)
    REAL(KIND=dp), INTENT(IN) :: Xi,A,B,Beta,rhow,rhos0,T0,Temperature,Pressure,Porosity
    REAL(KIND=dp) :: Tstar
    IF (Porosity <= 0.0) &
         CALL FATAL("PermafrostMaterials(GetXiAndersonT)","Zero or negative porosity detected")
    Tstar = T0 - Beta * Pressure - Temperature
    IF (Xi == 1.0_dp .OR. Xi == 0.0_dp) THEN
      XiAndersonT = 0.0_dp
    ELSE
      XiAndersonT = -(rhos0/rhow)*(A*B*(Tstar**(B - 1.0_dp)))/Porosity
    END IF
  END FUNCTION XiAndersonT
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION XiAndersonP(Xi,A,B,Beta,rhow,rhos0,T0,Temperature,Pressure,Porosity)
    REAL(KIND=dp), INTENT(IN) :: Xi,A,B,Beta,rhow,rhos0,T0,Temperature,Pressure,Porosity
    REAL(KIND=dp) :: Tstar
    IF (Porosity <= 0.0_dp) &
         CALL FATAL("PermafrostMaterials(GetXiAndersonP)","Zero or negative porosity detected")
    Tstar = T0 - Beta * Pressure - Temperature
    IF (Xi == 1_dp .OR. Xi == 0.0_dp) THEN
      XiAndersonP = 0.0_dp
    ELSE
      XiAndersonP = -Beta*(rhos0/rhow)*(A*B*(Tstar**(B - 1.0_dp)))/Porosity
    END IF
  END FUNCTION XiAndersonP
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION XiAndersonEta(Xi,A,B,Beta,rhow,rhos0,T0,Temperature,Pressure,Porosity) 
    REAL(KIND=dp), INTENT(IN) :: Xi,A,B,Beta,rhow,rhos0,T0,Temperature,Pressure,Porosity
    REAL(KIND=dp) :: Tstar
    IF (Porosity <= 0.0) &
         CALL FATAL("PermafrostMaterials(GetXiAndersonEta)","Zero or negative porosity detected")
    Tstar = T0 - Beta * Pressure - Temperature
    IF (Xi == 1.0_dp .OR. Xi == 0.0_dp) THEN
      XiAndersonEta = 0.0_dp
    ELSE
      XiAndersonEta =  -(rhos0/rhow)*(A*(Tstar**B))/(Porosity*Porosity)
    END IF
  END FUNCTION XiAndersonEta
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION delta(CurrentSolventMaterial,&
       eps,DeltaT,T0,GasConstant)
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp), INTENT(IN) :: eps,DeltaT,T0,GasConstant
    REAL(KIND=dp) :: aux,Mw,hi0,cw0,ci0
    LOGICAL :: FirstTime=.TRUE.
    SAVE FirstTime,Mw,hi0,ci0,cw0
    IF(FirstTime) THEN
      !assign needed variables
      Mw = CurrentSolventMaterial % Mw
      hi0 = CurrentSolventMaterial % hi0
      ci0 = CurrentSolventMaterial % ci0
      cw0 = CurrentSolventMaterial % cw0
      FirstTime = .FALSE.
    END IF
    aux = -0.5_dp*hi0*DeltaT/T0 &
         + (cw0 - ci0)*((T0 + 0.5_dp*DeltaT)*LOG(1.0_dp + 0.5_dp*DeltaT/T0) - 0.5_dp*DeltaT)
    delta = aux*(eps*(1.0_dp - eps)/(2.0_dp*eps - 1.0_dp))* Mw/(GasConstant*(T0 + 0.5_dp*DeltaT))
    !IF (delta < 1.0d-10) PRINT *, "delta=", delta, "(aux,Mw,T0,DeltaT,eps)",hi0, Mw, T0,DeltaT,eps
  END FUNCTION delta
  !---------------------------------------------------------------------------------------------
  FUNCTION GetAcAlphatilde(CurrentSolventMaterial,ComputeIce) RESULT(acAlphatilde)
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    LOGICAL, INTENT(IN) :: ComputeIce
    !---------------------------------
    REAL(KIND=dp) :: acAlphaTilde(0:5)
    REAL(KIND=dp) :: acAlpha(0:5), sumation
    INTEGER :: acAlphal,I
    !assign needed properties
    acAlphaTilde(0:5) = 0.0_dp
    IF (ComputeIce) THEN
      acAlphal = CurrentSolventMaterial % acil
      acAlpha(0:5) = CurrentSolventMaterial % aci(0:5)
    ELSE
      acAlphal = CurrentSolventMaterial % acwl
      acAlpha(0:5) = CurrentSolventMaterial % acw(0:5)
    END IF
    ! acAlphal-entries 
    sumation = 0.0_dp
    DO I=acAlphal,1,-1
      sumation = acAlpha(I)- sumation
      acAlphaTilde(I) = ( (1.0_dp/(DBLE(I) + 1.0_dp)) - 1.0/DBLE(I) ) *  sumation
    END DO
    ! zero-entry only for acAlphaTilde(0)
    acAlphaTilde(0)= acAlpha(0)- sumation
  END FUNCTION GetAcAlphatilde
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION gwa(CurrentSolventMaterial,&
       p0,T0,rhow,Temperature,Pressure)
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp), INTENT(IN) :: p0,T0,rhow,Temperature,Pressure
    REAL(KIND=dp) :: cw0,kw0,bcw(0:5)
    INTEGER :: I,bcwl
    REAL(KIND=dp) :: acwtilde(0:5),aux
    LOGICAL :: FirstTime=.TRUE.

    SAVE FirstTime,acwtilde

    IF (FirstTime) THEN
      acwtilde = GetAcAlphatilde(CurrentSolventMaterial,.FALSE.)
      FirstTime = .FALSE.
    END IF

    aux = -(CurrentSolventMaterial % cw0)*(acwtilde(0) * Temperature * LOG(Temperature/T0) &
         - (Temperature - T0) &
         * GeneralPolynomial(Temperature,T0,T0,acwtilde,CurrentSolventMaterial % acwl))
    gwa = aux + ((Pressure - p0)*(1.0_dp + 0.5_dp*(CurrentSolventMaterial % kw0)*(Pressure - p0))/rhow)

    IF (gwa .NE. gwa) THEN
      PRINT *, "gwa:", gwa
      PRINT *, GeneralPolynomial(Temperature,T0,T0,acwtilde,CurrentSolventMaterial % acwl)
      PRINT *, acwtilde(0), Temperature, LOG(Temperature/T0)
      PRINT *, ((Pressure - p0)*(1.0_dp + 0.5_dp*(CurrentSolventMaterial % kw0)*(Pressure - p0))/rhow)
      PRINT *,rhow, Pressure
      CALL FATAL("PermafrostMaterials(gwa)","Error in gwa")
    END IF
  END FUNCTION gwa
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION gia(CurrentSolventMaterial,&
       p0,T0,rhoi,Temperature,Pressure)
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp), INTENT(IN) :: p0,T0,rhoi,Temperature,Pressure
    !---------------------------------
    INTEGER :: acil,I
    REAL(KIND=dp) ::acitilde(0:5),aux,aux1
    LOGICAL :: FirstTime=.TRUE.
    !-----------------------------------
    SAVE FirstTime,acitilde

    IF (FirstTime) THEN
      acitilde = GetAcAlphatilde(CurrentSolventMaterial,.TRUE.)
      FirstTime = .FALSE.
    END IF
    aux = -(CurrentSolventMaterial % hi0)*((Temperature - T0)/T0)&
         - (CurrentSolventMaterial % ci0) *(acitilde(0) * Temperature * LOG(Temperature/T0) &
         - (Temperature - T0) * GeneralPolynomial(Temperature,T0,T0,acitilde,CurrentSolventMaterial % acil))
    gia = aux + ((Pressure - p0)*(1.0_dp + 0.5_dp*(CurrentSolventMaterial % ki0)*(Pressure - p0))/rhoi)
  END FUNCTION gia
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION gwaT(CurrentSolventMaterial,&
       p0,T0,rhow,Temperature)
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp), INTENT(IN) :: p0,T0,rhow,Temperature
    !----------------------------
    INTEGER :: I
    REAL(KIND=dp) :: acwtilde(0:5),aux, currpot
    LOGICAL :: FirstTime=.TRUE.

    SAVE FirstTime, acwtilde

    IF (FirstTime) THEN
      acwtilde = GetAcAlphatilde(CurrentSolventMaterial,.FALSE.)
      FirstTime = .FALSE.
    END IF

    ! get the derivative
    aux = 0.0_dp
    currpot = 1.0_dp
    DO i=0,CurrentSolventMaterial % acwl
      !IF (pdeg == 0) PRINT*,i,"of",pdeg, currpot, Variable, ReferenceValue, Normation
      aux = aux + (i + 1.0_dp)*acwtilde(i) * currpot
      currpot = currpot * (Temperature - T0)/T0
    END DO
    gwaT = -(CurrentSolventMaterial % cw0)*(acwtilde(0) * (1.0_dp + LOG(Temperature/T0)) - aux)
    
    ! neglected term
    !gwaT = aux &
    !* GeneralPolynomial(watercont,1.0_dp,1.0_dp,&
    !     CurrentSolventMaterial % bcw(0:5),&
    !     CurrentSolventMaterial % bcwl)

    IF (gwaT .NE. gwaT) THEN
      PRINT *, "gwaT"
      CALL FATAL("PermafrostMaterials(gwaT)","Error in gwaT")
    END IF
  END FUNCTION gwaT
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION giaT(CurrentSolventMaterial,&
       p0,T0,rhoi,Temperature)
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp), INTENT(IN) :: p0,T0,rhoi,Temperature
    INTEGER :: I
    REAL(KIND=dp) :: acitilde(0:5),aux, currpot
    LOGICAL :: FirstTime=.TRUE.

    SAVE FirstTime, acitilde

    IF (FirstTime) THEN
      acitilde = GetAcAlphatilde(CurrentSolventMaterial,.TRUE.)
      !FirstTime = .FALSE.
    END IF

    ! get the derivative
    aux = 0.0_dp
    currpot = 1.0_dp
    DO i=0,CurrentSolventMaterial % acil
      aux = aux + (i + 1.0_dp)*acitilde(i) * currpot
      currpot = currpot * (Temperature - T0)/T0
    END DO
    aux = (CurrentSolventMaterial % ci0)*(acitilde(0) * (1.0_dp + LOG(Temperature/T0)) - aux)
    giaT = -(CurrentSolventMaterial % hi0)/T0 - aux
    IF (giaT .NE. giaT) THEN
      PRINT *, "giaT"
      CALL FATAL("PermafrostMaterials(gwa)","Error in giaT")
    END IF
    FirstTime = .FALSE.
  END FUNCTION giaT
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION deltaG(gwa,gia)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: gwa,gia
    deltaG = gwa - gia
  END FUNCTION deltaG
  !---------------------------------------------------------------------------------------------
  FUNCTION GetBi(CurrentSoluteMaterial,RockMaterialID,&
       Xi0Tilde,Salinity,Update) RESULT(bi)
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    REAL(KIND=dp), INTENT(IN) :: Xi0Tilde,Salinity
    INTEGER,  INTENT(IN) :: RockMaterialID
    REAL(KIND=dp):: bi(4)
    LOGICAL :: Update
    !----------
    REAL(KIND=dp)::  aux,d1,d2,e1

    IF (Update) THEN
      e1 = GlobalRockMaterial % e1(RockMaterialID)
      bi(3) = (1.0_dp - Xi0Tilde)*e1
      bi(4) = Xi0Tilde*e1
    ELSE
      d1 = CurrentSoluteMaterial % d1
      d2 = CurrentSoluteMaterial % d2
      aux = Salinity/(1.0_dp - Salinity)
      bi(1) = aux*(d1 + 0.5_dp*d2*aux)
      bi(2) = aux*(d1 + d2*aux)/(1.0_dp - Salinity)
      bi(3) = 0.0_dp
      bi(4) = 0.0_dp
    END IF
  END FUNCTION GetBi
  !---------------------------------------------------------------------------------------------
  FUNCTION GetBiYc(CurrentSoluteMaterial,Salinity) RESULT(biYc)
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    REAL(KIND=dp), INTENT(IN) :: Salinity
    REAL(KIND=dp):: biYc(2)
    !----------
    REAL(KIND=dp)::  aux,d1,d2

    d1 = CurrentSoluteMaterial % d1
    d2 = CurrentSoluteMaterial % d2

    aux = 1.0_dp/(1.0_dp - Salinity)
    biYc(1) = (d1 + d2*Salinity*aux)*aux*aux
    biYc(2) = (d1*(1.0_dp + Salinity) + d2*Salinity*(2.0_dp + Salinity))*aux**3.0_dp
  END FUNCTION GetBiYc
  !---------------------------------------------------------------------------------------------
  FUNCTION GetB(RockMaterialID,CurrentSolventMaterial,&
       Xi0tilde,delta,deltaG,GasConstant,bi,Temperature) RESULT(B)
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp), INTENT(IN) :: Xi0tilde,delta,deltaG,GasConstant,bi(4),Temperature
    INTEGER, INTENT(IN) :: RockMaterialID
    REAL(KIND=dp) :: B
    REAL(KIND=dp) :: e1,Mw
    Mw = CurrentSolventMaterial % Mw
    e1 = GlobalRockMaterial % e1(RockMaterialID)

    B =(Mw*deltaG/(GasConstant*Temperature) - bi(1) + bi(3))/(delta + bi(2) + bi(4)) 
        
    IF (B .NE. B) THEN
      PRINT *, "B:", Mw, deltaG,Temperature,bi(1),e1,delta,bi(2),bi(4)
      CALL FATAL("PermafrostMaterials(getB)","Error in getB")
    END IF
  END FUNCTION GetB
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION D(RockMaterialID,delta,bi)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: RockMaterialID
    REAL(KIND=dp), INTENT(IN) :: delta,bi(4)
    ! local
    D = delta/(delta + bi(2) + bi(4))
    !IF (D .NE. D) THEN
    !  PRINT *, "D"
    !  CALL FATAL("PermafrostMaterials(getD)","Error in getD")
    !END IF
  END FUNCTION D
  !---------------------------------------------------------------------------------------------
  FUNCTION GetXi0Tilde(RockMaterialID,Porosity) RESULT(Xi0tilde)
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    INTEGER, INTENT(IN) :: RockMaterialID    
    REAL(KIND=dp), INTENT(IN) :: Porosity
    REAL(KIND=dp) :: Xi0tilde
    REAL(KIND=dp) :: Xi0,eta0
    LOGICAL :: FirstTime = .TRUE.
    SAVE FirstTime

    Xi0 = GlobalRockMaterial % Xi0(RockMaterialID)
    eta0 = GlobalRockMaterial % eta0(RockMaterialID)
    IF (Porosity <= 0.0_dp) THEN
      IF (Xi0 == 0.0_dp) THEN
        Xi0tilde = 1.0_dp
      ELSE
        CALL FATAL("PermafrostMaterials(GetXiTilde)","Zero or negative porosity detected")
      END IF
    ELSE
      !Xi0tilde = MIN(Xi0 * (eta0/Porosity) * (1.0_dp - Porosity)/(1.0_dp - eta0),1.0_dp)
      Xi0tilde = Xi0 * (eta0/Porosity) * (1.0_dp - Porosity)/(1.0_dp - eta0)
    END IF
  END FUNCTION GetXi0Tilde
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION fw(RockMaterialID,CurrentSolventMaterial,&
       Xi0tilde,rhow,Xi,GasConstant,Temperature)
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    INTEGER, INTENT(IN) :: RockMaterialID
    REAL(KIND=dp), INTENT(IN) :: Xi0tilde,rhow,Xi,GasConstant,Temperature
    !-----------------
    REAL(KIND=dp) :: Mw, e1
    !----------------
    !   IF (Xi > Xi0tilde) THEN
    !     fw = 0.0_dp
    !   ELSE
    !     e1 = GlobalRockMaterial % e1(RockMaterialID)
    !     Mw = CurrentSolventMaterial % Mw
    !     fw = rhow*GasConstant*Temperature*e1*Xi0tilde/(Mw*Xi)
    !   END IF
    fw = 0.0_dp !! CHANGE BACK WHEN JUHA TELLS US TO DO SO !!
  END FUNCTION fw
  !---------------------------------------------------------------------------------------------
  FUNCTION GetXi(B,D) RESULT(Xi)
    REAL(KIND=dp), INTENT(IN) :: B,D
    REAL(KIND=dp) :: Xi
    Xi= 1.0_dp/(1.0_dp + 0.5_dp*B + SQRT(0.25_dp*B*B + D))
  END FUNCTION GetXi
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION XiT(CurrentSolventMaterial,&
       B,D,Xi,bi,p0,delta,deltaG,T0,gwa,gia,gwaT,giaT,GasConstant,Temperature)
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp), INTENT(IN) :: B,D,Xi,bi(4),p0,delta,deltaG,&
         T0,gwa,gia,gwaT,giaT,GasConstant,Temperature
    !local
    REAL(KIND=dp) :: aux1, aux2, aux3, Mw,e1, hi0,hw0,rhow0,rhoi0,cw0,ci0
    LOGICAL :: FirstTime=.TRUE.

    SAVE Mw,hi0,hw0,rhow0,rhoi0,cw0,ci0, FirstTime

    IF (FirstTime) THEN
      Mw = CurrentSolventMaterial % Mw
      hi0   = CurrentSolventMaterial % hi0  
      hw0   = CurrentSolventMaterial % hw0  
      rhow0 = CurrentSolventMaterial % rhow0
      rhoi0 = CurrentSolventMaterial % rhoi0
      cw0   = CurrentSolventMaterial % cw0  
      ci0   = CurrentSolventMaterial % ci0  
      FirstTime=.FALSE.
    END IF
    aux1 = 1.0_dp/(delta + bi(2) + bi(4))
    aux2 = (1.0_dp + B/SQRT(B*B + 4.0_dp*D))
    aux3 = ((gwa - gia)/Temperature - (gwaT - giaT))
    XiT = 0.5_dp*(Mw/(GasConstant*Temperature))*aux1*aux2*aux3*Xi*Xi
  END FUNCTION XiT
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION XiP(CurrentSolventMaterial,&
       B,D,bi,Xi,gwap,giap,delta,GasConstant,Temperature)
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp), INTENT(IN) :: B,D,bi(4),Xi,gwap,giap,delta,GasConstant,Temperature
    !local
    REAL(KIND=dp) :: aux1, aux2, rhow0,rhoi0, Mw
    LOGICAL :: FirstTime=.TRUE.

    SAVE Mw,rhow0,rhoi0,FirstTime
    IF (FirstTime) THEN
      Mw = CurrentSolventMaterial % Mw
      rhow0 = CurrentSolventMaterial % rhow0
      rhoi0 = CurrentSolventMaterial % rhoi0
      FirstTime=.FALSE.
    END IF
    IF (Temperature <= 0.0_dp) CALL FATAL("Permafrost (XiP)","(sub-)Zero Temperature detected")
    aux1 = 1.0_dp/(delta + bi(2) + bi(4))
    aux2 = (1.0_dp + B/SQRT(B*B + 4.0_dp*D))
    XiP = 0.5_dp * aux1 * aux2 *(giap - gwap)* Mw/(GasConstant*Temperature)*Xi*Xi
  END FUNCTION XiP
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION XiYc(B,D,bi,biYc,Xi,delta)
    IMPLICIT NONE
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    REAL(KIND=dp), INTENT(IN) :: B,D,bi(4),biYc(2),Xi,delta
    !local
    REAL(KIND=dp) :: aux1, aux2, aux3, aux_sqrt

    aux_sqrt = B*B + 4.0_dp*D
    aux1 = 1.0_dp/(delta + bi(2) + bi(4))
    aux2 = ( 1.0_dp + B/SQRT(aux_sqrt) )*(biYc(1) + B*biYc(2)) &
         + biYc(2)/(SQRT(aux_sqrt))
    aux3 = 2.0_dp*D*biYc(2)/(SQRT(aux_sqrt))
    XiYc = 0.5_dp*aux1*(aux2 + aux3)*Xi*Xi      
    !IF (XiYc .NE. XiYc) THEN
    !  PRINT *, "XiYc:", aux1, aux2
    !  PRINT *, B, D, biYc(1), biYc(2),delta
    !  STOP
    !END IF
  END FUNCTION XiYc
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION XiEta(RockMaterialID,&
       B,D,bi,biYc,Xi,delta,Porosity)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: RockMaterialID    
    REAL(KIND=dp), INTENT(IN) :: B,D,bi(4),biYc(2),Xi,delta,Porosity
    !local
    REAL(KIND=dp) :: aux1, aux2, aux3, aux_sqrt,Xi0,eta0

    Xi0 = GlobalRockMaterial % Xi0(RockMaterialID)
    eta0 = GlobalRockMaterial % eta0(RockMaterialID)
    
    aux_sqrt = B*B + 4.0_dp*D
    
    IF (Porosity >= 0.0_dp) THEN
      aux1 = 1.0_dp/(delta + bi(2) + bi(4))
      aux2 = ( 1.0_dp + B/SQRT(aux_sqrt) )*(1.0_dp + B)
      aux3 = 2.0_dp*D*biYc(2)/SQRT(aux_sqrt)
      XiEta = 0.5_dp*aux1*(aux2 + aux3) * (Xi0*eta0/(1.0_dp - eta0))&
           *(1.0_dp/(Porosity**2.0_dp))*Xi*Xi
    ELSE
      CALL WARN("PermafrostMaterials(XiEta)","Porosity out of physical range - returning zero")
      XiEta = 0.0_dp
    END IF
  END FUNCTION XiEta
  !----------------------------------------------------------------------
  SUBROUTINE GetXiHartikainen (RockMaterialID,&
       CurrentSoluteMaterial,CurrentSolventMaterial,&
       TemperatureAtIP,PressureAtIP,SalinityAtIP,PorosityAtIP,&
       Xi0tilde,deltaInElement,rhowAtIP,rhoiAtIP,&
       GasConstant,p0,T0,&
       XiAtIP,XiTAtIP,XiYcAtIP,XiPAtIP,XiEtaAtIP,&
       ComputeXi,ComputeXiT, ComputeXiYc, ComputeXiP, ComputeXiEta)

    IMPLICIT NONE

    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    INTEGER :: RockMaterialID
    REAL(KIND=dp), INTENT(IN) :: Xi0tilde,deltaInElement,rhowAtIP,rhoiAtIP
    REAL(KIND=dp), INTENT(IN) :: GasConstant,p0,T0
    REAL(KIND=dp), INTENT(IN) :: TemperatureAtIP,PressureAtIP,SalinityAtIP,PorosityAtIP
    REAL(KIND=dp), INTENT(OUT) :: XiAtIP,XiTAtIP,XiYcAtIP,XiPAtIP,XiEtaAtIP
    LOGICAL, INTENT(IN) :: ComputeXi,ComputeXiT, ComputeXiYc, ComputeXiP, ComputeXiEta
    !---------------------------
    REAL(KIND=dp) :: biAtIP(4),biYcAtIP(2),gwaAtIP,gwaTAtIP,gwapAtIP,&
         giaAtIP,giaTAtIP,giapAtIP,deltaGAtIP,DAtIP,BAtIP
    !---------------------------
    IF (ComputeXi .OR. (ComputeXiT .OR. ComputeXiYC .OR. ComputeXiP)) THEN
      biAtIP = GetBi(CurrentSoluteMaterial,RockMaterialID,&
           Xi0Tilde,SalinityAtIP,.FALSE.) 
      gwaAtIP = gwa(CurrentSolventMaterial,&
           p0,T0,rhowAtIP,TemperatureAtIP,PressureAtIP)     
      giaAtIP = gia(CurrentSolventMaterial,&
           p0,T0,rhoiAtIP,TemperatureAtIP,PressureAtIP)
      deltaGAtIP = deltaG(gwaAtIP,giaAtIP)
      DAtIP= D(RockMaterialID,deltaInElement,biAtIP)
      BAtIP = GetB(RockMaterialID,CurrentSolventMaterial,&
           Xi0tilde,deltaInElement,deltaGAtIP,GasConstant,biAtIP,TemperatureAtIP)
    ELSE
      CALL WARN("GetXiHartikainen","Nothing to be done - why did you call this routine?")
    END IF
    IF (ComputeXi)  THEN
      XiAtIP = GetXi(BAtIP,DAtIP)
      IF (XiAtIP < 0.0_dp) THEN
        WRITE (Message,*) 'Dedected invalid value for Xi=', XiAtIP
        XiAtIP = 0.0_dp
        CALL WARN("GetXiHartikainen",Message)
      ELSE IF (XiAtIP > 1.0_dp) THEN
        WRITE (Message,*) 'Dedected invalid value for Xi=', XiAtIP
        XiAtIP = 0.99_dp
        CALL WARN("GetXiHartikainen",Message)
      END IF
    END IF
    
    ! updates of derivatives
    IF (XiAtIP < Xi0tilde)  THEN
      biAtIP = GetBi(CurrentSoluteMaterial,RockMaterialID,&
           Xi0Tilde,SalinityAtIP,.TRUE.)
      XiAtIP = GetXi(BAtIP,DAtIP)
      IF (XiAtIP <= 0.0_dp) THEN
        WRITE(Message,*) 'Xi=',XiAtIP,' reset to 0.001'
        CALL WARN("GetXiHartikainen",Message)
        XiAtIP = 0.001_dp
      END IF
    END IF
    !----------------------------------------------------
    XiTAtIP = 0.0_dp
    XiYcAtIP = 0.0_dp
    XiPAtIP = 0.0_dp
    IF (ComputeXiT) THEN
      giaTAtIP = giaT(CurrentSolventMaterial,&
           p0,T0,rhoiAtIP,TemperatureAtIP)
      gwaTAtIP =  gwaT(CurrentSolventMaterial,&
           p0,T0,rhowAtIP,TemperatureAtIP)!   
      XiTAtIP= XiT(CurrentSolventMaterial,&
           BAtIP,DAtIP,XiAtIP,biAtIP,p0,&
           deltaInElement,deltaGAtIP,T0,gwaAtIP,giaAtIP,gwaTAtIP,giaTAtIP,GasConstant,TemperatureAtIP)
    END IF
    IF (ComputeXiYC) THEN
      biYcAtIP = GetBiYc(CurrentSoluteMaterial,SalinityAtIP)
      XiYcAtIP = XiYc(BAtIP,DAtIP,biAtIP,biYcAtIP,XiAtIP,deltaInElement)
    END IF
    IF (ComputeXiP) THEN
      giapAtIP = 1.0_dp/rhoiAtIP
      gwapAtIP = 1.0_dp/rhowAtIP
      XiPAtIP = XiP(CurrentSolventMaterial,&
           BAtIP,DAtIP,biAtIP,gwapAtIP,giapAtIP,XiAtIP,&
           deltaInElement,GasConstant,TemperatureAtIP)
    END IF
  END SUBROUTINE GetXiHartikainen
  !---------------------------------------------------------------------------------------------
  ! Densities and their derivatives, thermal expansion, isothermal chemical compaction and
  !     compressibility coefficients
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION rhos(RockMaterialID,&
       T0,p0,Temperature,Pressure,ConstVal)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: RockMaterialID 
    REAL(KIND=dp), INTENT(IN) :: T0,p0,Temperature,Pressure
    LOGICAL :: ConstVal
    !----------------------
    REAL(KIND=dp) :: aux1,aux2
    !----------------------
    rhos = GlobalRockMaterial % rhos0(RockMaterialID)
    IF (.NOT.ConstVal) THEN
      !aux1 = GeneralIntegral(Pressure,p0,p0,&
      !     GlobalRockMaterial % ks0(RockMaterialID),&
      !     GlobalRockMaterial % cks(0:5,RockMaterialID),&
      !     GlobalRockMaterial % cksl(RockMaterialID))
      ! a shortcut, as only cks(0) = 1
      aux1 = (GlobalRockMaterial % ks0(RockMaterialID)) * (Pressure - p0)
      aux2 = GeneralIntegral(Temperature,T0,T0,&
           GlobalRockMaterial % as0(RockMaterialID),&
           GlobalRockMaterial % aas(0:5,RockMaterialID),&
           GlobalRockMaterial % aasl(RockMaterialID))
      rhos = rhos * EXP(aux1 - aux2)
    END IF
  END FUNCTION rhos
  !---------------------------------------------------------------------------------------------
  REAL(KIND=dp) FUNCTION rhosT(RockMaterialID,rhos,T0,Temperature)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: rhos,T0,Temperature
    INTEGER, INTENT(IN) :: RockMaterialID
    REAL(KIND=dp) :: alphaS

    alphaS = GlobalRockMaterial % as0(RockMaterialID) *&
	 GeneralPolynomial(Temperature,T0,T0,&
         GlobalRockMaterial % aas(0:5,RockMaterialID),&
         GlobalRockMaterial % aasl(RockMaterialID))
    rhosT = rhos * alphaS
  END FUNCTION rhosT
!---------------------------------------------------------------------------------------------
  REAL(KIND=dp) FUNCTION rhosp(RockMaterialID,rhos,p0,Pressure)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: RockMaterialID
    REAL(KIND=dp), INTENT(IN) :: rhos,p0,Pressure
    !--------------------
    REAL(KIND=dp) ::  kappas
    !--------------------
    kappas = ( GlobalRockMaterial % ks0(RockMaterialID))
    rhosP = rhos * kappas
  END FUNCTION rhosp
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION rhow(CurrentSolventMaterial,T0,p0,Temperature,Pressure,ConstVal)
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp), INTENT(IN) :: T0,p0,Temperature,Pressure
    LOGICAL :: ConstVal
    !----------------------
    REAL(KIND=dp) :: aux1, aux11, aux2, aux22, aux3, aux33, watercont
    !----------------------
    rhow = CurrentSolventMaterial % rhow0
    IF (.NOT.ConstVal) THEN
      !aux1 = GeneralIntegral(Pressure,p0,p0,&
      !     CurrentSolventMaterial % kw0,&
      !     CurrentSolventMaterial % ckw(0:5),&
      !     CurrentSolventMaterial % ckwl)
      ! a shortcut, as only ckw(0) = 1
      aux1 = (CurrentSolventMaterial % kw0) * (Pressure - p0)
      aux2 = GeneralIntegral(Temperature,T0,T0,&
           CurrentSolventMaterial % aw0,&
           CurrentSolventMaterial % aaw(0:5),&
           CurrentSolventMaterial % aawl)
      rhow = rhow * EXP(aux1 - aux2)
      IF (rhow < 800.0) THEN
        PRINT *, "rhow:",  rhow,CurrentSolventMaterial % rhow0,aux1, aux2,Pressure,Temperature
      END IF
      IF (rhow .NE. rhow) THEN
        PRINT *, "rhow:", rhow,CurrentSolventMaterial % rhow0,aux1, aux2,Pressure,Temperature
        CALL FATAL("PermafrostMaterials(rhow)","Error in rhow")
      END IF
    END IF
  END FUNCTION rhow
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION rhowupdate(CurrentSolventMaterial,&
       previousrhow,Xi,Salinity,ConstVal)
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp), INTENT(IN) :: previousrhow,Xi,Salinity
    LOGICAL :: ConstVal
    !----------------------
    REAL(KIND=dp) :: aux3, aux33, watercont
    !----------------------
    rhowupdate = previousrhow   
    IF (.NOT.ConstVal) THEN
      watercont = MAX(0.0_dp, 1.0_dp - Salinity/Xi)
      aux3 = GeneralIntegral(watercont,1.0_dp,1.0_dp,&
           CurrentSolventMaterial % zw0,&
           CurrentSolventMaterial % bzw(0:5),&
           CurrentSolventMaterial % bzwl)
      !aux33 = (watercont - 1.0_dp)* (CurrentSolventMaterial % zw0) *&
      !     ( (CurrentSolventMaterial % bzw(0))&
      !     + 0.5_dp*((CurrentSolventMaterial % bzw(1)) * (watercont - 1.0_dp))) 
      rhowupdate = previousrhow * EXP(aux3)
      !IF (aux3 .NE. aux33) THEN
      ! PRINT *, "rhowupdate:", previousrhow, EXP(aux3), aux3, aux33, watercont,Salinity,Xi
      ! PRINT *, "zw0",CurrentSolventMaterial % zw0, "bzw",CurrentSolventMaterial % bzw(0:CurrentSolventMaterial % bzwl)
      ! PRINT *, CurrentSolventMaterial % bzwl, CurrentSolventMaterial % zw0, "---", CurrentSolventMaterial % bzw(0:5)
      !END IF
      !IF (rhowupdate .NE. rhowupdate) THEN
      !  PRINT *, "rhowupdate:"
      !  STOP
      !END IF
    END IF
  END FUNCTION rhowupdate
  !---------------------------------------------------------------------------------------------
  REAL(KIND=dp) FUNCTION rhowT(CurrentSolventMaterial,rhow,T0,Temperature)
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp), INTENT(IN) :: rhow,T0,Temperature
    !--------------------
    REAL(KIND=dp) :: alphaW
    !--------------------
    alphaW = (CurrentSolventMaterial % aw0) *&
         GeneralPolynomial(Temperature,T0,T0,&
         CurrentSolventMaterial % aaw(0:5),&
         CurrentSolventMaterial % aawl)
    rhowT = rhow * alphaW
  END FUNCTION rhowT
  !---------------------------------------------------------------------------------------------
  REAL(KIND=dp) FUNCTION rhowP(CurrentSolventMaterial,rhow,p0,Pressure)
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp), INTENT(IN) :: rhow,p0,Pressure
    !--------------------
    REAL(KIND=dp) ::  kappaW
    !--------------------
    !kappaW = (CurrentSolventMaterial % kw0) *&
    !     GeneralPolynomial(Pressure,p0,p0,&
    !     CurrentSolventMaterial % ckw(0:5),&
    !     CurrentSolventMaterial % ckwl)
    ! a shortcut, as only ckw(0) = 1.0
    kappaW = (CurrentSolventMaterial % kw0)
    rhowP = rhow * kappaW
  END FUNCTION rhowP
  !---------------------------------------------------------------------------------------------
  REAL(KIND=dp) FUNCTION rhowYc(CurrentSolventMaterial,rhow,Xi,Salinity)
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp), INTENT(IN) :: rhow,Xi,Salinity
    !--------------------
    REAL(KIND=dp) ::  zetaW, xw    
    !--------------------
    xw = 1.0_dp - (Salinity/Xi)
    zetaW = (CurrentSolventMaterial % zw0) *&
         GeneralPolynomial(xw,1.0_dp,1.0_dp,&
         CurrentSolventMaterial % bzw(0:5),&
         CurrentSolventMaterial % bzwl)
    rhowYc = -zetaW *rhow/Xi
  END FUNCTION rhowYc
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION rhoi(CurrentSolventMaterial,T0,p0,Temperature,Pressure,ConstVal)
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp), INTENT(IN) :: T0,p0,Temperature,Pressure
    LOGICAL :: ConstVal
    !----------------------
    REAL(KIND=dp) :: aux1, aux2
    !----------------------
    rhoi = CurrentSolventMaterial % rhoi0
    IF (.NOT.ConstVal) THEN
      !aux1 = GeneralIntegral(Pressure,p0,p0,&
      !     CurrentSolventMaterial % ki0 ,&
      !     CurrentSolventMaterial % cki(0:5),&
      !     CurrentSolventMaterial % ckil)
      ! a shortcut, as only cki(0) = 1
      aux1 = (CurrentSolventMaterial % ki0 ) * (Pressure - p0)
      aux2 = GeneralIntegral(Temperature,T0,T0,&
           CurrentSolventMaterial % ai0,&
           CurrentSolventMaterial % aai(0:5),&
           CurrentSolventMaterial % aail)
      rhoi = rhoi * EXP(aux1 - aux2)
    END IF
  END FUNCTION rhoi
  !---------------------------------------------------------------------------------------------
  REAL(KIND=dp) FUNCTION rhoiT(CurrentSolventMaterial,rhoi,T0,Temperature)
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp), INTENT(IN) :: rhoi,T0,Temperature
    !------------------------
    REAL(KIND=dp) :: alphaI
    !------------------------
    alphaI = (CurrentSolventMaterial % ai0) *&
         GeneralPolynomial(Temperature,T0,T0,&
         CurrentSolventMaterial % aai(0:5),&
         CurrentSolventMaterial % aail) 
    rhoiT = rhoi * alphaI
  END FUNCTION rhoiT
  !---------------------------------------------------------------------------------------------
  REAL(KIND=dp) FUNCTION rhoiP(CurrentSolventMaterial,rhoi,p0,Pressure)
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp), INTENT(IN) :: rhoi,p0,Pressure
    !------------------------
    REAL(KIND=dp):: kappaI
    !------------------------
    !kappaI = (CurrentSolventMaterial % ki0) *&
    !     GeneralPolynomial(Pressure,p0,p0,&
    !     CurrentSolventMaterial % cki(0:5),&
    !     CurrentSolventMaterial % ckil)
    ! a shortcut, as only cki(0) = 1
    kappaI = CurrentSolventMaterial % ki0
    rhoiP = rhoi * kappaI
  END FUNCTION rhoiP
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION rhoc(CurrentSoluteMaterial,T0,p0,Xi,Temperature,Pressure,Salinity,ConstVal)
    IMPLICIT NONE
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    REAL(KIND=dp), INTENT(IN) :: T0,p0,Xi,Temperature,Pressure,Salinity
    LOGICAL :: ConstVal
    !----------------------
    REAL(KIND=dp) :: aux1, aux2, aux3, xc
    !----------------------
    rhoc = CurrentSoluteMaterial % rhoc0
    IF (.NOT.ConstVal) THEN
      xc = Salinity/Xi
      aux1 = GeneralIntegral(Pressure,p0,p0,&
           CurrentSoluteMaterial % kc0,&
           CurrentSoluteMaterial % ckc(0:5),&
           CurrentSoluteMaterial % ckcl)
      aux2 = GeneralIntegral(Temperature,T0,T0,&
           CurrentSoluteMaterial % ac0,&
           CurrentSoluteMaterial % aac(0:5),&
           CurrentSoluteMaterial % aacl)
      aux3 = GeneralIntegral(xc,0.0_dp,1.0_dp,&
           CurrentSoluteMaterial % zc0,&
           CurrentSoluteMaterial % bzc(0:5),&
           CurrentSoluteMaterial % bzcl)
      rhoc = rhoc * EXP(aux1 - aux2 + aux3)
    END IF
  END FUNCTION rhoc
  !---------------------------------------------------------------------------------------------
  REAL(KIND=dp) FUNCTION rhocT(CurrentSoluteMaterial,rhoc,T0,Temperature,ConstVal)
    IMPLICIT NONE
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    REAL(KIND=dp), INTENT(IN) :: rhoc,T0,Temperature
    LOGICAL :: ConstVal
    !-------------------------
    REAL(KIND=dp):: alphaC
    !-------------------------
    IF (ConstVal) THEN
      rhocT = 0.0_dp
    ELSE      
      alphaC = (CurrentSoluteMaterial % ac0) * &
           GeneralPolynomial(Temperature,T0,T0,&
           CurrentSoluteMaterial % aac(0:5),&
           CurrentSoluteMaterial % aacl) 
      rhocT = rhoc * alphaC
    END IF
  END FUNCTION rhocT
  !---------------------------------------------------------------------------------------------
  REAL(KIND=dp) FUNCTION rhocP(CurrentSoluteMaterial,rhoc,ConstVal)
    IMPLICIT NONE
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    REAL(KIND=dp), INTENT(IN) :: rhoc
    LOGICAL :: ConstVal
    !---------------
    IF (ConstVal) THEN
      rhocP = 0.0_dp
    ELSE
      rhocP = rhoc * (CurrentSoluteMaterial % kc0)
    END IF
  END FUNCTION rhocP
  !---------------------------------------------------------------------------------------------
  REAL(KIND=dp) FUNCTION rhocYc(CurrentSoluteMaterial,rhoc,Xi,Salinity,ConstVal)
    IMPLICIT NONE
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    REAL(KIND=dp), INTENT(IN) :: rhoc, Xi, Salinity
    LOGICAL :: ConstVal
    !---------------
    REAL(KIND=dp):: xc, zc
    !---------------
    IF (ConstVal) THEN
      rhocYc = 0.0_dp
    ELSE
      xc = Salinity/Xi
      zc = (CurrentSoluteMaterial % zc0) * &
           GeneralPolynomial(xc,0.0_dp,1.0_dp,&
           CurrentSoluteMaterial % bzc(0:5),&
           CurrentSoluteMaterial % bzcl)
      rhocYc = rhoc * zc /Xi
    END IF
  END FUNCTION rhocYc
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION rhogw(rhow,rhoc,Xi,Salinity)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: rhow,rhoc,Xi,Salinity
    !------------
    REAL(KIND=dp) :: xc
    !------------
    xc = MAX(Salinity/Xi,0.0_dp)
    rhogw = rhow + xc*(rhoc - rhow)
  END FUNCTION rhogw
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION rhogw_driesner(rhow,dummyrhoc,Xi,Salinity)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: rhow,dummyrhoc,Xi,Salinity
    !------------
    REAL(KIND=dp) :: xc   
    !------------
    xc = MAX(Salinity/Xi,0.0_dp)
    rhogw_driesner = rhow*(18.01528_dp * (1.0_dp - xc) + 58.443_dp * xc)/18.01528_dp
  END FUNCTION rhogw_driesner
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION rhogwP(rhowp,rhocp,Xi,Salinity)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: rhowp,rhocp,Xi,Salinity
    !------------
    REAL(KIND=dp) :: xc   
    !------------
    xc = Salinity/Xi    
    rhogwP = (1.0_dp - xc)*rhowP + xc*rhocP
  END FUNCTION rhogwP
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION rhogwT(rhowT,rhocT,Xi,Salinity)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: rhowT,rhocT,Xi,Salinity
    !------------
    REAL(KIND=dp) :: xc   
    !------------    
    xc = Salinity/Xi    
    rhogwT = (1.0_dp - xc)*rhowT + xc*rhocT
  END FUNCTION rhogwT
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION rhogwYc(rhow, rhoc, rhowYc,rhocYc,Xi,Salinity)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: rhow, rhoc, rhowYc,rhocYc,Xi,Salinity
    !------------
    REAL(KIND=dp) :: xc   
    !------------    
    xc = Salinity/Xi    
    rhogwYc = ((1.0_dp - xc)*rhowYc + xc*rhocYc +  rhow + rhoc)/Xi
  END FUNCTION rhogwYc
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION cs(RockMaterialID,T0,Temperature,ConstVal)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: RockMaterialID 
    REAL(KIND=dp), INTENT(IN) :: T0,Temperature
    LOGICAL :: ConstVal
    !----------------------
    REAL(KIND=dp) :: aux
    !----------------------
    cs = GlobalRockMaterial % cs0(RockMaterialID)
    !PRINT *,"cs:", Temperature,T0,GlobalRockMaterial % cs0(RockMaterialID),&
    !     GlobalRockMaterial % acs(0:5,RockMaterialID),GlobalRockMaterial % acsl(RockMaterialID)
    IF (.NOT.ConstVal) &
         cs = cs * GeneralPolynomial(Temperature,T0,T0,&
         GlobalRockMaterial % acs(0:5,RockMaterialID),&
         GlobalRockMaterial % acsl(RockMaterialID))
  END FUNCTION cs
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION cw(CurrentSolventMaterial,T0,Xi,Temperature,Salinity,ConstVal)
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp), INTENT(IN) :: T0,Xi,Temperature,Salinity
    LOGICAL :: ConstVal
    !----------------------
    REAL(KIND=dp) :: aux1, aux2, watercont
    !----------------------
    cw = CurrentSolventMaterial % cw0
    IF (.NOT.ConstVal) THEN
      watercont = MAX(1.0_dp - Salinity/Xi,0.0_dp)
      aux1 = GeneralPolynomial(Temperature,T0,T0,&
           CurrentSolventMaterial % acw(0:5),&
           CurrentSolventMaterial % acwl)
      aux2 = GeneralPolynomial(watercont,1.0_dp,1.0_dp,&
           CurrentSolventMaterial % bcw(0:5),&
           CurrentSolventMaterial % bcwl)
      cw = cw * aux1 * aux2
    END IF
  END FUNCTION cw
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION ci(CurrentSolventMaterial,&
       T0,Temperature,ConstVal)
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp), INTENT(IN) :: T0,Temperature
    REAL(KIND=dp) :: ci0
    REAL(KIND=dp), DIMENSION(0:5) :: aci
    INTEGER :: acil
    LOGICAL :: ConstVal
    !----------------------
    ci = CurrentSolventMaterial % ci0     
    IF (.NOT.ConstVal) THEN
      !PRINT *, "ci:",ci0,aci(0:5),acil
      ci = ci *&
           GeneralPolynomial(Temperature,T0,T0,&
           CurrentSolventMaterial % aci(0:5),&
           CurrentSolventMaterial % acil)
    END IF
  END FUNCTION ci
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION cc(CurrentSoluteMaterial,&
       T0,Temperature,Salinity,ConstVal)
    IMPLICIT NONE
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    REAL(KIND=dp), INTENT(IN) :: T0,Temperature,Salinity
    LOGICAL :: ConstVal
    !----------------------
    REAL(KIND=dp) :: aux1, aux2
    !----------------------
    cc = CurrentSoluteMaterial % cc0
    IF (.NOT.ConstVal) THEN
      aux1 = GeneralPolynomial(Temperature,T0,T0,&
           CurrentSoluteMaterial % acc(0:5),&
           CurrentSoluteMaterial % accl)
      aux2 = GeneralPolynomial(Salinity,0.0_dp,1.0_dp,&
           CurrentSoluteMaterial % bcc(0:5),&
           CurrentSoluteMaterial % bccl)
      cc = cc*aux1*aux2
    END IF
  END FUNCTION cc
  !---------------------------------------------------------------------------------------------
  ! latent heat of water  
  REAL (KIND=dp) FUNCTION hw(CurrentSolventMaterial,&
       T0,Xi,Temperature,Salinity,ConstVal)
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp), INTENT(IN) :: T0,Xi,Temperature,Salinity
    LOGICAL :: ConstVal
    !----------------------
    REAL(KIND=dp) :: aux1, aux2, watercont
    !----------------------
    hw = CurrentSolventMaterial % hw0 
    IF (.NOT.ConstVal) THEN
      watercont = MAX(1.0_dp - Salinity/Xi,0.0_dp)
      aux1 = GeneralPolynomial(watercont,1.0_dp,1.0_dp,&
           CurrentSolventMaterial % bcw(0:5),&
           CurrentSolventMaterial % bcwl)
      aux2 = GeneralIntegral(Temperature,T0,T0,&
           CurrentSolventMaterial % cw0,&
           CurrentSolventMaterial % acw(0:5),&
           CurrentSolventMaterial % acwl)
      hw = hw + aux1 * aux2
    END IF
  END FUNCTION hw
  !---------------------------------------------------------------------------------------------
  ! latent heat of ice  
  REAL (KIND=dp) FUNCTION hi(CurrentSolventMaterial,&
       T0,Temperature,ConstVal)
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp), INTENT(IN) :: T0,Temperature
    LOGICAL :: ConstVal
    !----------------------
    REAL(KIND=DP) :: hi0,ci0
    REAL(KIND=DP), DIMENSION(0:5) :: aci
    INTEGER :: acil
    LOGICAL :: FirstTime = .TRUE.
    SAVE FirstTime,hi0,ci0,aci,acil
    !----------------------
    hi0 = CurrentSolventMaterial % hi0
    IF (ConstVal) THEN
      hi = hi0
    ELSE
      IF (FirstTime) THEN
        ci0 = CurrentSolventMaterial % ci0 
        aci(0:5) = CurrentSolventMaterial % aci(0:5)
        acil = CurrentSolventMaterial % acil
        FirstTime = .FALSE.
      END IF
      hi = hi0 + GeneralIntegral(Temperature,T0,T0,ci0,aci,acil)
    END IF
  END FUNCTION hi
  !---------------------------------------------------------------------------------------------
  ! General constituent thermal conductivity: kalpha0th and balpha have to be directly transferred
  FUNCTION GetKAlphaTh(kalpha0th,balpha,T0,Temperature)RESULT(kalphath)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: kalpha0th,balpha,T0,Temperature
    REAL(KIND=dp) :: kalphath
    !-------------------------
    kalphath = kalpha0th/( 1.0_dp + balpha*(Temperature - T0)/T0)
    !kalphath = kalpha0th
  END FUNCTION GetKAlphaTh
  !---------------------------------------------------------------------------------------------
  FUNCTION GetCGTT(Xi,XiT,rhos,rhow,rhoi,rhoc,cw,ci,cs,cc,hi,hw,&
       Porosity,Salinity)RESULT(CGTT)! All state variables or derived values
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: Xi,XiT,rhos,rhow,rhoi,rhoc,cw,ci,cs,cc,&
         hi,hw,Porosity,Salinity
    REAL(KIND=dp) :: CGTT
    !-------------------------
    REAL(KIND=dp) :: xc
    !-------------------------
    xc = Salinity/Xi
    CGTT = (1.0_dp - Porosity)*rhos*cs &
         + (Xi - Salinity) * Porosity * rhow * cw & ! mind xc * Xi = Salinity
         + Salinity * Porosity * rhoc * cc & ! mind xc * Xi = Salinity
         + (1.0_dp - Xi)*Porosity*rhoi*ci &
         + rhoi*(hw - hi)*Porosity*XiT
  END FUNCTION GetCGTT
  !---------------------------------------------------------------------------------------------
  FUNCTION GetCGTp(rhoi,hi,hw,XiP,Porosity)RESULT(CGTp)! All state variables or derived values
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: rhoi,hi,hw,XiP,Porosity
    REAL(KIND=dp) :: CGTp
    !-------------------------
    CGTp = Porosity*rhoi*(hw - hi)*XiP
  END FUNCTION GetCGTp
  !---------------------------------------------------------------------------------------------
  FUNCTION GetCGTyc(rhoi,hi,hw,XiYc,Porosity)RESULT(CGTyc)! All state variables or derived values
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: rhoi,hi,hw,XiYc,Porosity
    REAL(KIND=dp) :: CGTyc
    !-------------------------
    CGTyc = Porosity*rhoi*(hw - hi)*XiYc
  END FUNCTION GetCGTyc
  !---------------------------------------------------------------------------------------------
  ! functions specific to groundwater flow
  !---------------------------------------------------------------------------------------------
  FUNCTION GetJgwD(Kgwpp,KgwpT,Kgw,gradp,gradT,Gravity,rhogw,DIM,CryogenicSuction) RESULT(JgwD)
    IMPLICIT NONE
    REAL (KIND=dp), INTENT(IN) :: Kgwpp(3,3),KgwpT(3,3),Kgw(3,3),gradp(3),gradT(3),Gravity(3),&
         rhogw
    REAL (KIND=dp)  :: JgwD(3)
    LOGICAL, INTENT(IN):: CryogenicSuction
    !-------------------------
    INTEGER, INTENT(IN) :: DIM
    INTEGER :: i
    REAL (KIND=dp) :: fluxp(3),fluxT(3),fluxg(3)
    !-------------------------
    fluxT(1:DIM) = 0.0_dp
    JgwD = 0.0_dp
    DO i=1,DIM
      fluxp(i) = -1.0_dp * SUM(Kgwpp(i,1:DIM)*gradp(1:DIM))
      IF (CryogenicSuction) &
           fluxT(i) = -1.0_dp * SUM(KgwpT(i,1:DIM)*gradT(1:DIM))
      fluxg(i) =   rhogw * SUM(Kgw(i,1:DIM)*Gravity(1:DIM))
    END DO
    JgwD(1:DIM) = fluxp(1:DIM) + fluxT(1:DIM) + fluxg(1:DIM)
  END FUNCTION GetJgwD
  !---------------------------------------------------------------------------------------------
  FUNCTION GetKGTT(ksth,kwth,kith,kcth,Xi,&
       Salinity,Porosity,meanfactor)RESULT(KGTT) ! All state variables or derived values
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: ksth,kwth,kith,kcth,Xi,&
         Salinity,Porosity,meanfactor
    REAL(KIND=dp) :: KGTT(3,3)
    !-------------------------
    REAL(KIND=dp) :: KGaTT, KghTT, unittensor(3,3),xc
    !-------------------------
    xc = Salinity/Xi
    unittensor=RESHAPE([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0], SHAPE(unittensor))
    KGhTT = 1.0_dp/((1.0_dp - Porosity)/ksth + (1.0_dp - xc)*Xi*Porosity/kwth &
         + xc*Porosity/kcth + (1.0_dp - Xi)*Porosity/kith)
    KGaTT = (1.0_dp - Porosity)*ksth + (1.0_dp - xc)*Xi*Porosity*kwth &
         + xc*Porosity*kcth + (1.0_dp - Xi)*Porosity*kith
    KGTT = unittensor*((1.0_dp - meanfactor)*KGhTT + meanfactor * KGaTT)
  END FUNCTION GetKGTT
  !---------------------------------------------------------------------------------------------
  FUNCTION  GetDtd(RockMaterialID,Xi,Porosity,JgwD)RESULT(Dtd)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: Xi,Porosity,JgwD(3)
    REAL(KIND=dp) :: Dtd(3,3)
    INTEGER, INTENT(IN) :: RockMaterialID
    !-------------------------
    REAL(KIND=dp) :: unittensor(3,3),absJgwD,alphaL,alphaT
    INTEGER :: I,J
    !-------------------------
    unittensor=RESHAPE([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0], SHAPE(unittensor))
    absJgwD = SQRT(SUM(JgwD(1:3)*JgwD(1:3)))
    IF(absJgwD > 0.0_dp) THEN
      alphaL = GlobalRockMaterial % alphaL(RockMaterialID)
      alphaT = GlobalRockMaterial % alphaT(RockMaterialID)
      DO I=1,3
        DO J=1,3
          Dtd(I,J) = alphaT*absJgwD*unittensor(I,J) + (alphaL - alphaT)*JgwD(I)*JgwD(J)/absJgwD
        END DO
      END DO
    ELSE
      Dtd = 0.0_dp
    END IF
  END FUNCTION GetDtd
  !---------------------------------------------------------------------------------------------
  FUNCTION GetCgwTT(rhow,rhoc,cw,cc,Xi,Salinity)RESULT(CgwTT)! All state variables or derived values
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: rhow,rhoc,cw,cc,Xi,Salinity
    REAL(KIND=dp) :: CgwTT
    !-------------------------
    REAL(KIND=dp) :: xc
    !-------------------------
    xc = Salinity/Xi
    CgwTT = (1.0_dp - xc)*rhow*cw + xc*rhoc*cc
  END FUNCTION GetCgwTT
  !---------------------------------------------------------------------------------------------
  FUNCTION GetCgwpp(rhogw,rhoi,rhos,rhogwp,rhoip,rhosp,&
       kappaG,Xi,Xip,&
       RockMaterialID,Porosity)RESULT(Cgwpp)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: rhogw,rhoi,rhos,rhogwp,rhoip,rhosp,kappaG,Xi,Xip,Porosity
    INTEGER, INTENT(IN) :: RockMaterialID
    REAL(KIND=dp) :: Cgwpp
    !-------------------------
    Cgwpp = Porosity * ( (rhogw - rhoi) * Xip  + Xi * rhogwp + (1.0_dp - Xi)*rhoip ) &
         + (Xi * rhogw + (1.0_dp - Xi)*rhoi) * ( (1.0_dp - Porosity) * rhosp/rhos  + kappaG )
  END FUNCTION GetCgwpp
  !---------------------------------------------------------------------------------------------
  FUNCTION GetCgwpT(rhogw,rhoi,rhos,rhogwT,rhoiT,rhosT,Xi,XiT,Porosity)RESULT(CgwpT)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: rhogw,rhoi,rhos,rhogwT,rhoiT,rhosT,Xi,XiT,Porosity
    REAL(KIND=dp) :: CgwpT
    !-------------------------
    CgwpT = Porosity * ( (rhogw - rhoi) * XiT  + Xi * rhogwT + (1.0_dp - Xi)*rhoiT ) &
         + (Xi * rhogw + (1.0_dp - Xi)*rhoi) *(1.0_dp - Porosity) * rhosT/rhos
  END FUNCTION GetCgwpT
  !---------------------------------------------------------------------------------------------
  FUNCTION GetCgwpYc(rhogw,rhoi,rhogwYc,Xi,XiYc,Porosity)RESULT(CgwpYc)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: rhogw,rhoi,rhogwYc,Xi,XiYc,Porosity
    REAL(KIND=dp) :: CgwpYc
    !-------------------------
    CgwpYc = Porosity * ( (rhogw - rhoi) * XiYc  + Xi * rhogwYc )
  END FUNCTION GetCgwpYc
  !---------------------------------------------------------------------------------------------
  FUNCTION GetCgwpI1(rhogw,rhoi,Xi,kappaG,RockMaterialID)RESULT(CgwpI1)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: rhogw,rhoi,Xi,kappaG
    INTEGER, INTENT(IN) :: RockMaterialID
    REAL(KIND=dp) :: CgwpI1
    !-------------------------
    REAL(KIND=dp) :: kappas
    !-------------------------
    kappas = GlobalRockMaterial % ks0(RockMaterialID)
    CgwpI1 = (Xi * rhogw + (1.0_dp - Xi) * rhoi)*(kappaG)/3.0_dp
  END FUNCTION GetCgwpI1
  !---------------------------------------------------------------------------------------------
  REAL (KIND=dp) FUNCTION mugw(CurrentSolventMaterial,CurrentSoluteMaterial,&
       Xi,T0,Salinity,Temperature,ConstVal)
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    REAL(KIND=dp), INTENT(IN) :: Xi,T0,Salinity,Temperature
    LOGICAL :: ConstVal
    !-------------------------
    REAL(KIND=dp) :: nu1, nu2, xc
    !-------------------------
    IF (.NOT.ConstVal) THEN
      xc = MAX(Salinity/Xi,0.0_dp)
      xc = MIN(xc,0.2_dp)
      nu1 = (CurrentSolventMaterial % nu10) *&
           GeneralPolynomial(Temperature,T0,T0,&
           CurrentSolventMaterial % anw(0:5),&
           CurrentSolventMaterial % anwl)
      !      PRINT *,"mugw:anw,anwl,nu1", anw(0:5),anwl, nu1
      !      PRINT *,"mugw:bnc", bnc(0:5),bncl
      !      PRINT *,"bnc=",bnc(0:5),"bncl=",bncl,"xc=",xc
      !      STOP
      nu2 = (CurrentSoluteMaterial % nu20) *&
           GeneralPolynomial(xc,0.0_dp,1.0_dp,&
           CurrentSoluteMaterial % bnc(0:5),&
           CurrentSoluteMaterial % bncl)
      mugw = CurrentSolventMaterial % muw0 * &
           EXP(nu1 * (Temperature - T0) + nu2 * (xc - 0.0))
      IF ((mugw .NE. mugw) .OR. (mugw > HUGE(mugw))) THEN
        WRITE (Message,*) 'invalid value:' , mugw,&
             'Input values: muw0/nu1/nu2/T/T0/xc/Salinity/Xi: ',&
             CurrentSolventMaterial % muw0, nu1,nu2,Temperature, T0, xc, Salinity, Xi
       ! PRINT *, nu1 * (Temperature - T0), nu2 * (xc - 0.0)
        CALL FATAL("PermafrostMaterials(mugw)",Message)
      END IF
    ELSE
      mugw = CurrentSolventMaterial % muw0
    END IF
  END FUNCTION mugw
  !---------------------------------------------------------------------------------------------
  FUNCTION GetKGpe( RockMaterialID,CurrentSolventMaterial,Xi)RESULT(KGpe)
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    INTEGER, INTENT(IN) :: RockMaterialID 
    REAL(KIND=dp), INTENT(IN) :: Xi
    REAL(KIND=dp) :: KGpe(3,3)
!--------------------------
    REAL(KIND=dp) :: muw0,rhow0,qexp,Kgwh0(3,3),factor
    REAL(KIND=dp), PARAMETER :: gval=9.81_dp !hard coded, so match Kgwh0 with this value
    INTEGER :: I, J
    !-------------------------
    muw0 = CurrentSolventMaterial % muw0
    rhow0 = CurrentSolventMaterial % rhow0
    qexp = GlobalRockMaterial % qexp(RockMaterialID)
    Kgwh0(1:3,1:3) = GlobalRockMaterial % Kgwh0(1:3,1:3,RockMaterialID) ! hydro-conductivity
    ! transformation factor from hydr. conductivity to permeability hydr. conductivity tensor
    factor = (muw0)*(Xi**qexp)/(rhow0*gval)
    DO I=1,3
      DO J=1,3
        KGpe(i,j) = Kgwh0(i,j)*factor
      END DO
    END DO
  END FUNCTION GetKGpe
  !---------------------------------------------------------------------------------------------
  FUNCTION GetXikG0hy(RockMaterialID,Xi)RESULT(XikG0hy)
     IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    INTEGER, INTENT(IN) :: RockMaterialID 
    REAL(KIND=dp), INTENT(IN) :: Xi
    REAL(KIND=dp) :: XikG0hy(3,3)
    !--------------------------
    REAL(KIND=dp) :: Kgwh0(3,3), qexp
    Kgwh0(1:3,1:3) = GlobalRockMaterial % Kgwh0(1:3,1:3,RockMaterialID) ! hydro-conductivity
    qexp = GlobalRockMaterial % qexp(RockMaterialID)
    XikG0hy = (Xi**qexp) * Kgwh0
    !PRINT *,  "GetXikG0hy", XikG0hy, Xi,qexp,Kgwh0(1:3,1:3)
  END FUNCTION GetXikG0hy
  !---------------------------------------------------------------------------------------------
  FUNCTION GetKgw(RockMaterialID,CurrentSolventMaterial,mugw,Xi,MinKgw) RESULT(Kgw)
    
    IMPLICIT NONE
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    INTEGER, INTENT(IN) :: RockMaterialID
    REAL(KIND=dp), INTENT(IN) :: Xi,MinKgw,mugw
    REAL(KIND=dp) :: Kgw(3,3)
 
    !--------------------------
    REAL(KIND=dp) :: muw0,rhow0,qexp,Kgwh0(3,3),factor
    REAL(KIND=dp), PARAMETER :: gval=9.81_dp !hard coded, so match Kgwh0 with this value
    INTEGER :: I, J
    !-------------------------
    IF (mugw <= 0.0_dp) &
         CALL FATAL("PermafrostMaterials(GetKgw)","Unphysical viscosity detected")
    muw0 = CurrentSolventMaterial % muw0
    rhow0 = CurrentSolventMaterial % rhow0
    qexp = GlobalRockMaterial % qexp(RockMaterialID)
        ! this is OK
    !PRINT *,"Kgw:",muw0,mugw,rhow0,Kgwh0,Xi,factor

    Kgwh0(1:3,1:3) = GlobalRockMaterial % Kgwh0(1:3,1:3,RockMaterialID) ! hydro-conductivity
    ! transformation factor from hydr. conductivity to permeability hydr. conductivity tensor
    factor = (muw0/mugw)*(Xi**qexp)/(rhow0*gval)

    Kgw = 0.0_dp
    DO I=1,3
      DO J=1,3
        Kgw(i,j) = Kgwh0(i,j)*factor
      END DO
    END DO
    DO I=1,3
      Kgw(i,i) = MAX(Kgw(i,i),MinKgw)
    END DO
  END FUNCTION GetKgw
  !---------------------------------------------------------------------------------------------
  FUNCTION GetKgwpT(fw,XiT,Kgw)RESULT(KgwpT) ! All state variables or derived values
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: fw,XiT,Kgw(3,3)
    REAL(KIND=dp) :: KgwpT(3,3)
    !-------------------------
    KgwpT(1:3,1:3) = fw*XiT*Kgw(1:3,1:3)
  END FUNCTION GetKgwpT
  !---------------------------------------------------------------------------------------------
  FUNCTION GetKgwpp(fw,XiP,Kgw)RESULT(Kgwpp)! All state variables or derived values
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: fw,XiP,Kgw(3,3)
    REAL(KIND=dp) :: Kgwpp(3,3)
    !-------------------------
    Kgwpp(1:3,1:3) = (1.0_dp + fw*XiP)*Kgw(1:3,1:3)
  END FUNCTION GetKgwpp
  !---------------------------------------------------------------------------------------------
  ! functions specific to solute transport
  !---------------------------------------------------------------------------------------------
  FUNCTION GetKc(RockMaterialID,Dm,Xi,JgwD,Porosity)RESULT(Kc) 
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: Dm,Xi,JgwD(3),Porosity
    INTEGER, INTENT(IN) :: RockMaterialID
    REAL(KIND=dp) :: alphaL,alphaT,Kc(3,3), unittensor(3,3), aux, eL(3),absJgwD
    INTEGER :: I,J
    !-------------------------
    unittensor=RESHAPE([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0], SHAPE(unittensor))
    IF (Porosity <= 0.0_dp) &
         CALL FATAL("GetKc","Negative/Zero Porosity detected")
    IF (Xi <= 0.0_dp) &
         CALL FATAL("GetKc","Negative/Zero water content detected")
    Kc =  Dm * unittensor
    absJgwD = SQRT(SUM(JgwD(1:3) * JgwD(1:3)))
    IF (absJgwD > 0.0_dp) THEN
      alphaL = GlobalRockMaterial % alphaL(RockMaterialID)
      alphaT = GlobalRockMaterial % alphaT(RockMaterialID)
      eL = JgwD/absJgwD
      aux = absJgwD/(Porosity * Xi)   
      DO I=1,3
        DO J=1,3
          Kc(I,J) = Kc(I,J)  &
               + aux*((alphaL - alphaT)*eL(I)*eL(J)  + alphaT * unittensor(I,J))
        END DO
      END DO
    END IF
  END FUNCTION GetKc
  !---------------------------------------------------------------------------------------------
  FUNCTION GetConstKc(DispersionCoefficient)RESULT(Kc)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: DispersionCoefficient
    REAL(KIND=dp) :: Kc(3,3)   
    !-------------------------
    REAL(KIND=dp) :: unittensor(3,3)
    unittensor=RESHAPE([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0], SHAPE(unittensor))
    Kc = DispersionCoefficient  * unittensor
  END FUNCTION GetConstKc
  !---------------------------------------------------------------------------------------------
  REAL(KIND=dp) FUNCTION Dm(CurrentSoluteMaterial,N0,GasConstant,rhoc,mugw,Temperature)
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    REAL(KIND=dp), INTENT(IN) :: N0,GasConstant,rhoc,mugw,Temperature
    !-------------------------
    REAL(KIND=dp), PARAMETER :: bconst = 3.0 * PI
    REAL(KIND=dp) :: Mc, lcbar
    !-------------------------
    Mc = CurrentSoluteMaterial % Mc
    lcbar = (Mc/(rhoc * N0))**(1.0_dp/3.0_dp)
    Dm = GasConstant * Temperature / (bconst * mugw * lcbar * N0)
  END FUNCTION Dm
  !---------------------------------------------------------------------------------------------
  FUNCTION GetR(CurrentSoluteMaterial,CurrentSolventMaterial,GasConstant,rhow,rhoc,Xi,Temperature,Salinity) RESULT(r12)
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp), INTENT(IN) :: GasConstant,rhow,rhoc,Xi,Salinity,Temperature
    REAL(KIND=dp) :: r12(2)
    REAL(KIND=dp) :: d1, d2, Mc, Mw, aux, epsilonc
    !-------------------------
    d1 = CurrentSoluteMaterial % d1
    d2 = CurrentSoluteMaterial % d2
    Mc = CurrentSoluteMaterial % Mc
    Mw = CurrentSolventMaterial % Mw
    epsilonc = (Mc/Mw)*(rhow/rhoc)
    aux = Salinity/(Xi - Salinity)
    r12(1) = (1/epsilonc)*Mc*(1.0_dp - Salinity/Xi)/(rhoc * GasConstant * Temperature)
    r12(2) = epsilonc * ( d1 + (d1 + d2)*aux + d2*aux*aux )
  END FUNCTION GetR
  !---------------------------------------------------------------------------------------------
  FUNCTION  GetKcYcYc(Kc,r12) RESULT(KcYcYc)! All state variables or derived values
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: Kc(3,3),r12(2)
    REAL(KIND=dp) :: KcYcYc(3,3)
    !-------------------------
    KcYcYc(1:3,1:3) = r12(2) * Kc(1:3,1:3) 
  END FUNCTION GetKcYcYc
  !---------------------------------------------------------------------------------------------
  FUNCTION GetFc(rhoc,rhow,Gravity,r12,XiT,XiP,Xi,gradP,gradT) RESULT(fc)! All state variables or derived values
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: rhoc,rhow,Gravity(3),r12(2),XiT,XiP,Xi,gradP(3),gradT(3)
    REAL(KIND=dp) :: fc(3)
    !-------------------------
    fc(1:3) = r12(1)*(rhoc - rhow)*Gravity(1:3) + r12(2)*(XiT*gradT(1:3) + XiP*gradP(1:3))/Xi
  END FUNCTION GetFc
  !---------------------------------------------------------------------------------------------
  FUNCTION GetJcF(KcYcYc,Kc,fc,GradSalinity,Salinity) RESULT(JcF)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: KcYcYc(3,3),Kc(3,3),fc(3), GradSalinity(3),Salinity
    REAL(KIND=dp) :: JcF(3)
    !-----------
    !REAL(KIND=dp) :: diffFlux(3), extForce(3)
    INTEGER :: i
    DO i=1,3
      !diffFlux(i) = SUM(KcYcYc(i,1:3) * GradSalinity(1:3))
      !extForce(i) = SUM(Kc(i,1:3) * fc(1:3)) * Salinity
      JcF(i) = -SUM(KcYcYc(i,1:3)*GradSalinity(1:3)) +  SUM(Kc(i,1:3) * fc(1:3)) * Salinity
      !SUM(KgwAtIP(i,1:DIM)*Gravity(1:DIM))
    END DO
  END FUNCTION GetJcF
  !---------------------------------------------------------------------------------------------
  REAL(KIND=dp) FUNCTION CcYcT(rhocT,Porosity,Salinity)! All state variables or derived values
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: rhocT,Porosity, Salinity
    CcYcT = Porosity*Salinity*rhocT
  END FUNCTION CcYcT
  !---------------------------------------------------------------------------------------------
  REAL(KIND=dp) FUNCTION CcYcP(rhocP,Porosity, Salinity)! All state variables or derived values
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: rhocP,Porosity, Salinity
    !-------------------------
    CcYcP = Porosity*Salinity*rhocp
  END FUNCTION CcYcP
  !---------------------------------------------------------------------------------------------
  REAL(KIND=dp) FUNCTION CcYcYc(rhoc,rhocYc,Porosity, Salinity)! All state variables or derived values
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: rhoc,rhocYc,Porosity, Salinity
    !-------------------------
    CcYcYc = Porosity*(rhoc + Salinity*rhocYc)
  END FUNCTION CcYcYc
  !---------------------------------------------------------------------------------------------
  REAL(Kind=dp) FUNCTION RadiogenicHeatProduction(RockMaterialID,Depth,RefDepth)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: Depth,RefDepth
      INTEGER, INTENT(IN) :: RockMaterialID
    !---------
    RadiogenicHeatProduction = GlobalRockMaterial % RadGen(RockMaterialID) &
         * EXP(-Depth/RefDepth)    
  END FUNCTION RadiogenicHeatProduction
  !---------------------------------------------------------------------------------------------
  ! functions specific to ground deformation
  !---------------------------------------------------------------------------------------------
  REAL(Kind=dp) FUNCTION EG(CurrentSolventMaterial,RockMaterialID,Xi,Porosity)    
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: Xi,Porosity
     TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    INTEGER, INTENT(IN) :: RockMaterialID
    !PRINT *,"EG", Porosity, Xi, RockMaterialID, GlobalRockMaterial % Es0(RockMaterialID)
    EG = (1.0_dp - Porosity)*(GlobalRockMaterial % Es0(RockMaterialID))&
         /(1.0_dp - (GlobalRockMaterial % eta0(RockMaterialID))) &
         + Porosity * (1.0_dp - Xi) * (CurrentSolventMaterial % Ei0)
  END FUNCTION EG
  !---------------------------------------------------------------------------------------------
  REAL(Kind=dp) FUNCTION nuG(CurrentSolventMaterial,RockMaterialID,Xi,Porosity)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: Xi,Porosity
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    INTEGER, INTENT(IN) :: RockMaterialID
    !---------
    nuG = (1.0_dp - Porosity)*(GlobalRockMaterial % nuS0(RockMaterialID))&
         +  Porosity * (1.0_dp - Xi) * (CurrentSolventMaterial % nui0)
  END FUNCTION nuG
  !---------------------------------------------------------------------------------------------
  REAL(Kind=dp) FUNCTION betaG(CurrentSolventMaterial,RockMaterialID,Xi,Porosity)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: Xi,Porosity
     TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    INTEGER, INTENT(IN) :: RockMaterialID
    !---------
    betaG = (1.0_dp - Porosity)*(GlobalRockMaterial % betas(RockMaterialID)&
         +  Porosity * (1.0_dp - Xi) * (CurrentSolventMaterial % betai))
  END FUNCTION BetaG
  !---------------------------------------------------------------------------------------------
  REAL(Kind=dp) FUNCTION rhoG(rhos,rhogw,rhoi,Porosity,Salinity,Xi)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: rhos,rhogw,rhoi,Porosity,Salinity,Xi
    !---------
    rhoG = (1.0_dp - Porosity)*rhos + Porosity*Xi*(1.0_dp - Salinity)*rhogw &
         + Porosity*(1.0_dp - Xi)*rhoi
  END FUNCTION rhoG

  !---------------------------------------------------------------------------------------------
  FUNCTION KGuu(EG,nuG,DIM) RESULT(OutKGuu)
    REAL(KIND=dp), INTENT(IN) :: EG,nuG
    REAL(KIND=dp) OutKGuu(6,6)
    INTEGER, INTENT(IN) :: DIM
    !----------
    INTEGER :: I,J
    !---------
    OutKGuu = 0.0_dp
    DO I=1,DIM
      OutKGuu(I,I) = 1.0_dp - nuG
      OutKGuu(DIM+I,DIM+I) = 0.5_dp - nuG
      DO J=1,DIM
        IF (J /= I) OutKGuu(I,J) = nuG
      END DO
    END DO
  END FUNCTION KGuu
  !---------------------------------------------------------------------------------------------
  REAL(Kind=dp) FUNCTION kappaG(EG,nuG) ! needed directly in Darcy Model
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: EG,nuG
    !---------
    kappaG = (3.0_dp*(1.0_dp - 2.0_dp * nuG))/EG
  END FUNCTION KappaG
  !---------------------------------------------------------------------------------------------
END MODULE PermafrostMaterials
