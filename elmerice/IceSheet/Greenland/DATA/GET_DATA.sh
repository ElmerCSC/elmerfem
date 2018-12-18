#### GET MINIMAL DATA REQUIRED TO RUN THE GREENLAND APPLICATION
function usage {
 echo "This code will download and process standard data files to run GrIS applications"
 echo "  1- TOPOGRAPHY: BedMachineGreenland-2017-09-20.nc"
 echo "  2- OBS. VELOCITY: MEaSUREs Multi-year Greenland Ice Sheet Velocity Mosaic, Version 1"
 echo "  3- SLIP COEFFICIENT: greenland_elmerice_C1_v0.nc"
 echo "  4- VISCOSITY:  greenland_elmerice_MuMean_v0.nc"
 echo "  5- SMB : MAR v5.3 "
 echo 
 echo "requirement: "
 echo "  1- An account to download data from nsidc (https://nsidc.org/nsidc/register) "
 echo "  2- curl and wget to download the data files"
 echo "      (optionnally you will get the link to download teh files from your webbrowser)"
 echo "  3- http://www.gdal.org/   to convert from .tif to netcdf and reproject"
 echo "  4- nco to manipulate netcdf files (http://nco.sourceforge.net/)"
 echo
 exit 1
}


main() {
####################################################################
## GET TOPOGRAPHY
# IceBridge BedMachine Greenland, Version 3
#https://nsidc.org/data/idbmg4#
####################################################################
TOPOGRAPHY=BedMachineGreenland-2017-09-20.nc
if [ ! -f "$TOPOGRAPHY" ]; then
   echo "****************************************************************"
   echo "   $TOPOGRAPHY file not found downloading from nsidc using curl"
   echo "****************************************************************"

   sh scripts/GET_FROM_NSIDC.sh https://daacdata.apps.nsidc.org/pub/DATASETS/IDBMG4_BedMachineGr/BedMachineGreenland-2017-09-20.nc
fi
FILE_IS_HERE $TOPOGRAPHY

####################################################################
## GET VELOCITY DATA
#MEaSUREs Multi-year Greenland Ice Sheet Velocity Mosaic, Version 1
##http://nsidc.org/data/NSIDC-0670/versions/1#
####################################################################
VELOCITY_X=greenland_vel_mosaic250_vx_v1.tif
if [ ! -f "$VELOCITY_X" ]; then
   echo "****************************************************************"
   echo "$VELOCITY_X file not found downlowding from nsidc using curl"
   echo "****************************************************************"
   sh scripts/GET_FROM_NSIDC.sh https://n5eil01u.ecs.nsidc.org/MEASURES/NSIDC-0670.001/1995.12.01/greenland_vel_mosaic250_vx_v1.tif
fi
FILE_IS_HERE $VELOCITY_X

VELOCITY_Y=greenland_vel_mosaic250_vy_v1.tif
if [ ! -f "$VELOCITY_Y" ]; then
   echo "****************************************************************"
   echo "$VELOCITY_Y file not found downlowding from nsidc using curl"
   echo "****************************************************************"
   sh scripts/GET_FROM_NSIDC.sh https://n5eil01u.ecs.nsidc.org/MEASURES/NSIDC-0670.001/1995.12.01/greenland_vel_mosaic250_vy_v1.tif
fi
FILE_IS_HERE $VELOCITY_Y

#####################################
# Convert to netcdf
VELOCITY=greenland_vel_mosaic250_v1.nc
if [ ! -f "$VELOCITY" ]; then
   echo "****************************************************************"
   echo "$VELOCITY not found"
   echo "****************************************************************"
   sh scripts/VELOCITY_TO_NETCDF.sh $VELOCITY_X $VELOCITY_Y $VELOCITY
fi
FILE_IS_HERE $VELOCITY


####################################################################
### GET Slip Coefficient from elmerice
##  
####################################################################
Slip=greenland_elmerice_C1_v0.nc
if [ ! -f "$Slip" ]; then
   echo "****************************************************************"
   echo "$Slip not found"
   echo "****************************************************************"
   if [ ! -f "greenland_elmerice_C1_v0.nc.tgz" ]; then
     echo "****************************************************************"
     echo " Downloading greenland_elmerice_C1_v0.tgz "
     echo " http://elmerfem.org/elmerice/wiki/lib/exe/fetch.php?media=eis:greenland:present:greenland_elmerice_C1_v0.nc.tgz"
     echo "****************************************************************"
     wget http://elmerfem.org/elmerice/wiki/lib/exe/fetch.php?media=eis:greenland:present:greenland_elmerice_C1_v0.nc.tgz -O greenland_elmerice_C1_v0.nc.tgz
   fi
   tar xzf greenland_elmerice_C1_v0.nc.tgz
fi
FILE_IS_HERE $Slip

####################################################################
### GET Viscosity from elmerice
##  
####################################################################
Mu=greenland_elmerice_MuMean_v0.nc
if [ ! -f "$Mu" ]; then
   echo "****************************************************************"
   echo "$Mu not found"
   echo "****************************************************************"
   if [ ! -f "greenland_elmerice_mumean_v0.tgz" ]; then
     echo "****************************************************************"
     echo " Downloading greenland_elmerice_mumean_v0.tgz "
     echo " http://elmerfem.org/elmerice/wiki/lib/exe/fetch.php?media=eis:greenland:present:greenland_elmerice_mumean_v0.tgz"
     echo "****************************************************************"
     wget http://elmerfem.org/elmerice/wiki/lib/exe/fetch.php?media=eis:greenland:present:greenland_elmerice_mumean_v0.tgz -O greenland_elmerice_mumean_v0.tgz
   fi
   tar xzf greenland_elmerice_mumean_v0.tgz
fi
FILE_IS_HERE $Mu

####################################################################
### GET SMB FORCING
## GET SMB 21st century FORCING FROM MAR:
## get e.g. MARv3.5-yearly-CanESM2-rcp85-2006-2100.nc
####################################################################
SMB=MARv3.5-yearly-CanESM2-rcp85-2006-2100.nc 
link=ftp://ftp.climato.be/fettweis/MARv3.5.2/Greenland/CanESM2-rcp85_2006-2100_25km/monthly_outputs_interpolated_at_5km/MARv3.5-yearly-CanESM2-rcp85-2006-2100.nc
SMB_PROJ=MARv3.5-yearly-CanESM2-rcp85-2006-2100_EPSG3413.nc

if [ ! -f "$SMB_PROJ" ]; then

   if [ ! -f "$SMB" ] ; then
      echo "****************************************************************"
      echo " $SMB not found"
      echo "****************************************************************"
      echo "downlading data file from $link "
      wget $link
   fi 
   FILE_IS_HERE $SMB

   echo
   echo
   echo "****************************************************************"
   echo "Project SMB forcing to EPSG3413:"
   sh scripts/MAR_Bamber5km_to_EPSG3413.sh  $SMB
   echo 
   echo "****************************************************************"
fi 
FILE_IS_HERE $SMB_PROJ


## GET SMB FORCING FROM MAR:
## get e.g. MARv3.5.2-10km-yearly-ERA-Interim-1979-2014.nc
if [ ! -f "MARv3.5.2-10km-yearly-ERA-Interim-1979-2014.nc" ] ; then
   echo "****************************************************************"
   echo "MARv3.5.2-10km-yearly-ERA-Interim-1979-2014.nc not found"
   echo "****************************************************************"
   echo "downlading data file from ftp://ftp.climato.be/fettweis/"
   wget ftp://ftp.climato.be/fettweis/MARv3.5.2/Greenland/ERA-int_1979-2014_10km/monthly_outputs_interpolated_at_5km/MARv3.5.2-10km-yearly-ERA-Interim-1979-2014.nc
   echo
fi
FILE_IS_HERE MARv3.5.2-10km-yearly-ERA-Interim-1979-2014.nc

if [ ! -f "MARv3.5.2-10km-yearly-ERA-Interim-1979-2014_EPSG3413.nc" ]; then
   echo
   echo "****************************************************************"
   echo "Project SMB forcing to EPSG3413:"
   sh scripts/MAR_Bamber5km_to_EPSG3413.sh  MARv3.5.2-10km-yearly-ERA-Interim-1979-2014.nc
   echo
   echo "****************************************************************"
fi
FILE_IS_HERE MARv3.5.2-10km-yearly-ERA-Interim-1979-2014_EPSG3413.nc

if [ ! -f "MARv3.5.2-10km-yearly-ERA-Interim-1979-1999_Mean_EPSG3413.nc" ] ; then
   echo 
   echo "****************************************************************"
   echo "compute the 1979-1999 mean:"
      ncra -F -d TIME,1,21,1 MARv3.5.2-10km-yearly-ERA-Interim-1979-2014_EPSG3413.nc MARv3.5.2-10km-yearly-ERA-Interim-1979-1999_Mean_EPSG3413.nc
   echo "****************************************************************"
fi 
FILE_IS_HERE MARv3.5.2-10km-yearly-ERA-Interim-1979-1999_Mean_EPSG3413.nc

echo
echo
echo "SUCCESS! - Now you can run Greenland simulations - enjoy!"
echo
echo

}

##########################
##########################
function FILE_IS_HERE {
  if [ ! -f "$1" ]; then
   echo
   echo "File $1 still not found. Aborting !!"
   echo " see http://elmerice.elmerfem.org/wiki/doku.php?id=eis:greenland#present"
   echo " for details on how to get the dataset"
   echo
   exit
  else
   echo "File $1 is present"
  fi
}

##########################
##########################
if [ ! $# -eq 0 ]
  then
    usage
fi



main

