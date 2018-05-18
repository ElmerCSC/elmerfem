#!/bin/bash
Ice_DENSITY=910.0 # kg/m^3

#####################################
#####################################
function usage {
 echo "This code will reproject MARv3.5 forcing"
 echo "from the Bamber 5km grid to EPSG 3413 (Standard Nort Pole stereopolar projection for ISMIP6)"
 echo "convert SMBCORR variable from mm (water eq.) to m (ice eq.) using an ice density of $Ice_DENSITY kg/m^3" 
 echo 
 echo "Files can be downloaded from: ftp://ftp.climato.be/fettweis/"
 echo "use yearly files interpolated on the 5km grid"
 echo "usage : "
 echo "   sh MAR_Bamber5km_to_EPSG3413.sh <FILE_TO_PROCESS>"
 echo
 exit 1
}

#####################################
#####################################
main() {


FILE_IN=${1}

# check file dexist
if [ ! -f "$FILE_IN" ]; then
   echo "File $FILE does not exist. --ABORT-- "
   exit 1
fi

# number of tile steps
nTime=$(ncdmnsz TIME $FILE_IN)
# test that this is not empty
if [ -z "$nTime" ]; then
   echo "There is no TIME dimension?. --ABORT-- "
   exit 1
fi

echo "There is $nTime timesteps in input file $FILE_IN"

filename=$(basename "$FILE_IN")
filename="${filename%.*}"
FILE_OUT=$filename"_EPSG3413.nc"
echo "projecting all timesteps in $FILE_OUT"

if [ -f "$FILE_OUT" ]; then
   echo "output file already exist.  --ABORT--"
   exit 1
fi


#### Output Grid definition
## Grid resolution
res=5000
## Grid cell centers
xmin=-720000
xmax=960000
ymin=-3450000
ymax=-570000
#######################
# get corner defintion rq. this will perform only integer operation.
xlc=$(($xmin-$res/2))
xrc=$(($xmax+$res/2))
ylc=$(($ymin-$res/2))
yuc=$(($ymax+$res/2))
#echo 'grid corners :' $xlc $xrc $ylc $yuc
###############

for (( c=1; c<=nTime; c++ ))
do  
   # remove temporary files
   rm -rf tmp[1-4].nc
   echo "Reproject Time $c"

   ## extract current_time
   ncks  -O -F -v SMBCORR -d TIME,$c,$c,1 $FILE_IN tmp1.nc

   #Fix the corrdinates to the Bamber 5km grid definition
   ncap2 -O -v -s "smb=SMBCORR" -s "X=1.0*(X-5)-800.0" -s "Y=1.0*(Y-5)-3400.0" tmp1.nc tmp2.nc

   ### reprojection using gdalwrap
   gdalwarp -overwrite -of netcdf \
   -s_srs "+proj=stere +lat_0=90 +lat_ts=71 +lon_0=-39 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs" \
   -t_srs EPSG:3413  -tr $res $res -te $xlc $ylc $xrc $yuc \
   NETCDF:"tmp2.nc":smb tmp3.nc
   
   # make a time record dimension
   ncecat -O -u TIME -v Band1 tmp3.nc tmp4.nc

   # Append current time step to output File
   ncrcat -h --rec_apn tmp4.nc $FILE_OUT

done
# remove temporary files
rm -rf tmp[1-4].nc

## copy time variable dimension
ncks -h -A -v TIME $FILE_IN $FILE_OUT

## fix variable names and attributes
ncrename -h -v Band1,SMBCORR $FILE_OUT
ncap2 -h -A -s "smb_ice=SMBCORR/$Ice_DENSITY" $FILE_OUT
ncatted -h -a long_name,smb_ice,d,, -a long_name,SMBCORR,o,c,"Surface Mass Balance (with corrections)" -a units,smb_ice,o,c,"mIE/yr (Ice density=$Ice_DENSITY kg m^-3)" -a comment,global,a,c,"This is File $FILE_IN re-projected to EPSG:3413\n" $FILE_OUT

}


# ncdmnsz $dmn_nm $fl_nm : What is dimension size?
function ncdmnsz { ncks --trd -m -M ${2} | grep -E -i ": ${1}, size =" | cut -f 7 -d ' ' | uniq ; }
##

if [ $# -eq 0 ]
  then
    echo "No arguments supplied -- ABORT --"
    echo 
    usage
fi

main "$@"
