#!/bin/bash

function usage {
 echo 
 echo "sh VELOCITY_TO_NETCDF.sh FILE_Ux.tif FILE_Ux.tif FILE.nc"
 echo 
 echo "   convert and concatenate .tif files of the ux and uy component of the obs. velocity "
 echo "   to netcdf FILE.nc"
 echo
 exit 1
}

if [ ! $# -eq 3 ]
  then
    echo "not enough arguments supplied -- Aborting !!"
    echo 
    usage
fi


echo "convert vx ($1) and vy ($2)  files to netcdf ($3)"

## convert ux to netcdf
gdal_translate  -of netcdf -a_nodata -2000000000 $1 $3
ncrename -h -v Band1,vx $3

## convert uy to netcdf
gdal_translate  -of netcdf -a_nodata -2000000000 $2 tmp.nc
ncrename -h -v Band1,vy tmp.nc

## append data files 
ncks -h -A tmp.nc $3

## compute velocity norm
ncap2 -h  -A -s"vnorm=sqrt(vx*vx+vy*vy)" $3
ncatted -h -a _FillValue,vnorm,o,f,0.0 $3

rm -rf tmp.nc


