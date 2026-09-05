!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This library is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU Lesser General Public
! *  License as published by the Free Software Foundation; either
! *  version 2.1 of the License, or (at your option) any later version.
! *
! *  This library is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! *  Lesser General Public License for more details.
! * 
! *  You should have received a copy of the GNU Lesser General Public
! *  License along with this library (in file ../LGPL-2.1); if not, write 
! *  to the Free Software Foundation, Inc., 51 Franklin Street, 
! *  Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Authors: fabien Gillet-Chaulet
! *  Email:   fabien.gillet-chaulet@univ-grenoble-alpes.fr
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 3 May 2022
! *
! ******************************************************************************/
!--------------------------------------------------------------------------------
!>  Module containing utility subroutines for geographical transformations
!>  utility routines for fwd (LonLat => xy) and inv. (xy => LonLat) projections
!>  Currently supported projections:
!     > polar stereographic projections north and south
!     > generic from proj definition (requires proj library)
!--------------------------------------------------------------------------------
      MODULE ProjUtils
      USE DefUtils

#ifdef HAVE_PROJ
      USE proj6_interface
#endif

      IMPLICIT NONE

      INTERFACE proj_inv
        MODULE PROCEDURE projinv_proj6,projinv_stereo
      END INTERFACE
      INTERFACE proj_fwd
        MODULE PROCEDURE projfwd_proj6,projfwd_stereo
      END INTERFACE

      LOGICAL :: PjInitialized=.FALSE.
      CHARACTER(LEN=MAX_NAME_LEN) :: proj_type,proj_string
      REAL(KIND=dp) :: RefLon,RefLat
      REAL(KIND=dp) :: MinLon=0._dp,MinLat=0._dp,MaxLon=90._dp,MaxLat=90._dp
      REAL(KIND=dp) :: xmax,ymax,xmin,ymin

#ifdef HAVE_PROJ
      TYPE(pj_t) :: pj
      TYPE(pj_area_t) :: pj_area
      CHARACTER(LEN=*),PARAMETER :: proj_stringf = 'EPSG:4326'
#endif

      !WGS84 ellipsoid parameters : radius and flattening 
      real(kind=dp),parameter :: a=6378137.0_dp,f=1._dp/298.257223563
      real(kind=dp),parameter :: e2=2._dp*f-f*f,e=sqrt(e2)

      real(kind=dp),parameter :: rad2deg=180._dp/Pi
      real(kind=dp),parameter :: deg2rad=Pi/180._dp
      CONTAINS

!------------------------------------------------------------------------------
! generic xy2LonLat
!------------------------------------------------------------------------------
      SUBROUTINE xy2LonLat(x,y,Lon,Lat)
        REAL(KIND=dp),INTENT(IN) :: x,y
        REAL(KIND=dp),INTENT(OUT) :: lon,lat

        IF (.NOT.PjInitialized) CALL ProjINIT

        SELECT CASE(proj_type)
          CASE('polar stereographic north')
             CALL proj_inv(x,y,Lon,Lat,RefLon,RefLat)
             Lon=Lon*rad2deg
             Lat=Lat*rad2deg
          CASE('polar stereographic south')
             CALL proj_inv(-x,-y,Lon,Lat,-RefLon,-RefLat)
             Lon=-Lon*rad2deg
             Lat=-Lat*rad2deg

          CASE('regular')
             Lon=MinLon+(x-xmin)*(MaxLon-MinLon)/(xmax-xmin)
             Lat=MinLat+(y-ymin)*(MaxLat-MinLat)/(ymax-ymin)

#ifdef HAVE_PROJ
          CASE('proj')
             CALL proj_inv(x,y,Lon,Lat)
#endif

          CASE DEFAULT
            CALL FATAL('xy2LonLat','unsuported projection type: '//TRIM(proj_type))
        END SELECT
      END SUBROUTINE xy2LonLat

!------------------------------------------------------------------------------
! generic LonLat (degrees) => xy (m)  
!------------------------------------------------------------------------------
      SUBROUTINE LonLat2xy(Lon,Lat,x,y)
        REAL(KIND=dp),INTENT(IN) :: lon,lat
        REAL(KIND=dp),INTENT(OUT) :: x,y

        IF (.NOT.PjInitialized) CALL ProjINIT

        SELECT CASE(proj_type)
          CASE('polar stereographic north')
             CALL proj_fwd(Lon*deg2rad,Lat*deg2rad,x,y,RefLon,RefLat)
          CASE('polar stereographic south')
             CALL proj_fwd(-Lon*deg2rad,-Lat*deg2rad,x,y,-RefLon,-RefLat)
             x=-x
             y=-y
#ifdef HAVE_PROJ
          CASE('proj')
             CALL proj_fwd(Lon,Lat,x,y)
#endif

          CASE DEFAULT
            CALL FATAL('LonLat2xy','unsuported projection type: '//TRIM(proj_type))
        END SELECT
      END SUBROUTINE LonLat2xy

!------------------------------------------------------------------------------
! Initialise proj variables
!------------------------------------------------------------------------------
      SUBROUTINE ProjINIT
         TYPE(Mesh_t), POINTER :: Mesh
         LOGICAL :: Parallel
         LOGICAL :: GotIt
         REAL(KIND=dp) :: val

         PjInitialized=.TRUE.

         proj_type=ListGetString(GetSimulation(),'projection type',UnFoundFatal=.True.)

         SELECT CASE(proj_type)
            CASE('polar stereographic north')
               RefLon = ListGetConstReal(GetSimulation(),'central_meridian',UnFoundFatal=.True.)
               RefLon = RefLon*deg2rad
               RefLat = ListGetConstReal(GetSimulation(),'latitude_of_origin',UnFoundFatal=.True.)
               RefLat = RefLat*deg2rad
            CASE('polar stereographic south')
               RefLon = ListGetConstReal(GetSimulation(),'central_meridian',UnFoundFatal=.True.)
               RefLon=RefLon*deg2rad
               RefLat = ListGetConstReal(GetSimulation(),'latitude_of_origin',UnFoundFatal=.True.)
               RefLat = RefLat*deg2rad

            CASE('regular')
               IF(ParEnv % PEs > 1) Parallel = .TRUE.
               Mesh => GetMesh()     
               xmax=MAXVAL(Mesh % Nodes % x )
               ymax=MAXVAL(Mesh % Nodes % y )
               xmin=MINVAL(Mesh % Nodes % x )
               ymin=MINVAL(Mesh % Nodes % y )

               IF (Parallel) THEN
                  xmax=ParallelReduction(xmax,2)
                  ymax=ParallelReduction(ymax,2)
                  xmin=ParallelReduction(xmin,1)
                  ymin=ParallelReduction(ymin,1)
              END IF

              val = ListGetConstReal(GetSimulation(),'Min Latitude',GotIt)
              IF (GotIt) MinLat=val
              val = ListGetConstReal(GetSimulation(),'Min Longitude',GotIt)
              IF (GotIt) MinLon=val
              val = ListGetConstReal(GetSimulation(),'Max Latitude',GotIt)
              IF (GotIt) MaxLat=val
              val = ListGetConstReal(GetSimulation(),'Max Longitude',GotIt)
              IF (GotIt) MaxLon=val

              IF ((MaxLat-MinLat).LT.0._dp) &
                 CALL FATAL("ProjINIT","(MaxLat-MinLat) <= 0")
              IF ((MaxLon-MinLon).LT.0._dp) &
                 CALL FATAL("ProjINIT","(MaxLon-MinLon) <= 0")

#ifdef HAVE_PROJ
           CASE('proj')
               !> Get the mesh projection CRS     
               proj_string=ListGetString(GetSimulation(),'CRS code',UnFoundFatal=.True.)

               ! > detroy pj pointer in case it has already been initialised
               IF (proj_associated(pj)) pj=proj_destroy(pj)

               ! > Create projection from EPSG:4326 (proj_stringf) to mesh CRS
               ! > use deafult_context for now
               ! > do not use pj_area  for now
               pj = proj_create_crs_to_crs(pj_default_ctx,TRIM(proj_stringf)//CHAR(0), &
                      TRIM(proj_string)//CHAR(0), pj_area_null)

               !> Raise error and stop is something wrong
               IF (.NOT.proj_associated(pj)) THEN
                     CALL print_proj_error(pj_default_ctx)
                     CALL FATAL('ProjINIT','proj not associated')
               ENDIF
               ! > Print some projection informations
               CALL print_proj_info(pj)
#endif
          CASE('proj4')
             CALL FATAL('ProjINIT',& 
                   'projection type <proj4> deprecated - use proj - see elmerice/Utils/Documentation/ProjUtils.md')

            CASE DEFAULT
               CALL FATAL('ProjINIT','unsuported projection type: '//TRIM(proj_type))

         END SELECT
      END SUBROUTINE ProjINIT

!------------------------------------------------------------------------------
! proj4 : Inverse projection : x,y (m) => Lon,Lat (degrees)
!------------------------------------------------------------------------------
      SUBROUTINE projinv_proj6(x,y,lon,lat)
        REAL(KIND=dp),INTENT(IN) :: x,y
        REAL(KIND=dp),INTENT(OUT) :: lon,lat
#ifdef HAVE_PROJ
        TYPE(pj_coord_t) :: coordp,coordg

        coordp%x = x
        coordp%y = y
        coordg = proj_trans(pj, pj_inv, coordp)
        lon = coordg % y
        lat = coordg % x
#else
        CALL FATAL('projinv_proj6','proj not supported')
#endif
      END SUBROUTINE projinv_proj6

!------------------------------------------------------------------------------
! proj4 : fwd projection  Lon,Lat (degrees) => x,y (m)
!------------------------------------------------------------------------------
      SUBROUTINE projfwd_proj6(lon,lat,x,y)
        REAL(KIND=dp),INTENT(IN) :: lon,lat
        REAL(KIND=dp),INTENT(OUT) :: x,y
#ifdef HAVE_PROJ
        TYPE(pj_coord_t) :: coordp,coordg

        coordg%x=lat
        coordg%y=lon
        coordp = proj_trans(pj, pj_fwd, coordg)
        x = coordp % x
        y = coordp % y
#else
        CALL FATAL('proj_proj','proj not supported')
#endif
      END SUBROUTINE projfwd_proj6


!------------------------------------------------------------------------------
! north polar stereographic: Inverse projection : x,y (m) => Lon,Lat (rad)
!  From J. P. Snyder, Map-projections - A working Manual, 1987  
!------------------------------------------------------------------------------
      SUBROUTINE projinv_stereo(x,y,lon,lat,rLon,rLat)
        REAL(KIND=dp),INTENT(IN) :: x,y
        REAL(KIND=dp),INTENT(IN) :: rLon,rLat
        REAL(KIND=dp),INTENT(OUT) :: lon,lat
        REAL(KIND=dp) :: phi_c,phi1,phi,lambda,tc,mc,rho,t,change
        REAL(KIND=dp),parameter :: tol=1.0d-9

        phi_c=rlat

        ! Eq. 20-16
        lambda=rlon+atan2(x,-y)

        tc=tan(pi/4._dp-phi_c/2._dp)/((1._dp-e*sin(phi_c))/(1._dp+e*sin(phi_c)))**(e/2._dp)
        mc=cos(phi_c)/sqrt(1._dp-e2*sin(phi_c)*sin(phi_c))
        rho=sqrt(x*x+y*y)
        !Eq. 21-40
        t=rho*tc/(a*mc)

        phi1=pi/2._dp-2._dp*atan(t)

        !Eq. 7-9
        phi=pi/2._dp-2._dp*atan(t*((1._dp-e*sin(phi1))/(1._dp+e*sin(phi1)))**(e/2._dp))
        change=2._dp*tol

        do while (change.GT.tol)
          phi1=phi
          phi=pi/2._dp-2._dp*atan(t*((1._dp-e*sin(phi1))/(1._dp+e*sin(phi1)))**(e/2._dp))
          change=abs(phi-phi1)
        enddo
          lat=phi
          lon=lambda

      END SUBROUTINE projinv_stereo

!------------------------------------------------------------------------------
! north polar stereographic: fwd projection  Lon,Lat (rad) => x,y (m)
!  From J. P. Snyder, Map-projections - A working Manual, 1987  
!------------------------------------------------------------------------------
      SUBROUTINE projfwd_stereo(lon,lat,x,y,rLon,rlat)
        REAL(KIND=dp),INTENT(IN) :: lon,lat
        REAL(KIND=dp),INTENT(IN) :: rLon,rLat
        REAL(KIND=dp),INTENT(OUT) :: x,y
        REAL(KIND=dp) :: mc,tc,t,r

        mc=cos(rlat)/sqrt(1._dp-e2*sin(rlat)*sin(rlat))
        tc=sqrt(((1._dp-sin(rlat))/(1._dp+sin(rlat)))*((1._dp+e*sin(rlat))/(1._dp-e*sin(rlat)))**e)

        ! Eq. 15-9a
        t=sqrt(((1._dp-sin(lat))/(1._dp+sin(lat)))*((1._dp+e*sin(lat))/(1._dp-e*sin(lat)))**e)

        ! Eq. 21-34
        r=a*mc*t/tc

        ! Eqs. 21-30; 21-31
        x=r*sin(lon-rlon)
        y=-r*cos(lon-rlon)

      END SUBROUTINE projfwd_stereo

#ifdef HAVE_PROJ

      SUBROUTINE print_proj_info(pj)
       TYPE(pj_t), INTENT(in) :: pj
       TYPE(pj_info_t)   :: prinfo
       TYPE(pj_proj_info_t)   :: pjinfo
       CHARACTER(*), PARAMETER :: Caller="proj_info"
       INTEGER, PARAMETER :: olevel=3

       CALL INFO(Caller,"==== Projection Initialisation ====", level=olevel)

       !> Get info about proj version:
       prinfo = proj_info()
       CALL INFO(Caller," > ==== Proj library ====", level=olevel)
       CALL INFO(Caller,"  release: "//TRIM(cstrtof(prinfo%release)), level=olevel)

       !> Get info about the transformation:
       pjinfo=proj_pj_info(pj)
       CALL INFO(Caller," > ==== Requested transformation ====", level=olevel)
       CALL INFO(Caller,'  id: '//TRIM(cstrtof(pjinfo%id)), level=olevel)
       CALL INFO(Caller,'  desc: '//TRIM(cstrtof(pjinfo%description)), level=olevel)
       CALL INFO(Caller,'  def: '//TRIM(cstrtof(pjinfo%definition)), level=olevel)
       CALL INFO(Caller,'  has inv: '//I2S(pjinfo%has_inverse), level=olevel)
       write(Message,*) ' accuracy: ',pjinfo%accuracy
       CALL INFO(Caller,Message, level=olevel)
       CALL INFO(Caller,'============================',level=olevel)

      END SUBROUTINE print_proj_info

      SUBROUTINE print_proj_error(ctx)
        TYPE(pj_context_t),INTENT(in) :: ctx
        CHARACTER(*), PARAMETER :: Caller="proj_error"
        INTEGER :: errno

        errno = proj_context_errno(ctx)
        IF (errno /= 0) THEN
         CALL INFO(Caller, &
                 TRIM(cstrtof(proj_context_errno_string(ctx,errno))), &
                 level=1)
        ENDIF

      END SUBROUTINE print_proj_error
#endif


      END MODULE ProjUtils
      
