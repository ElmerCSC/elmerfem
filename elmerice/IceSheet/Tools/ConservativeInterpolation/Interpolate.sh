#! *****************************************************************************/
#! ******************************************************************************
#! *
#! *  Authors: F. Gillet-Chaulet (IGE-France)
#! *  Web:     http://elmerice.elmerfem.org
#! *  Original Date: 04/2019
#! * 
#! *****************************************************************************
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#  Example script for conservative projection of SMB and aSMB using CDO
#
## Parameters
DATA_DIR=DATA
PROJ_DIR=PROJ_DATA
mkdir -p $PROJ_DIR

aSMB_BASE_NAME=aSMB_MARv3.9-yearly-MIROC5-rcp85

STRUCT_GRID=grid_ISMIP6_GrIS_01000m.nc
ELMER_GRID=MESH_CDOGrid.txt
##
##

## First generate projection weigths
echo '***************************************************************************'
echo '***************************************************************************'
echo '***  generate conservative projection weights '
echo '     from structured grid : ' $STRUCT_GRID
echo '     to unstructured : ' $ELMER_GRID
echo '***************************************************************************'
echo '***************************************************************************'
cdo genycon,$ELMER_GRID -selname,SMB -setgrid,$STRUCT_GRID $DATA_DIR/MARv3.9-ERA-Interim-1980-1999.nc PROJ_WEIGHTS.nc

## Now remap
InputFile=$DATA_DIR/MARv3.9-ERA-Interim-1980-1999.nc
if [[ ! -f "$InputFile" ]]; then
    echo "$InputFile does not exist; aborting"
    exit
fi
OutPutFile=$PROJ_DIR/MARv3.9-ERA-Interim-1980-1999-UNST.nc
echo '***************************************************************************'
echo '***************************************************************************'
echo '***  remap file : ' $InputFile
echo '***************************************************************************'
echo '***************************************************************************'
cdo remap,$ELMER_GRID,PROJ_WEIGHTS.nc -selname,SMB -setgrid,$STRUCT_GRID $InputFile $OutPutFile

# Remap aSMB Files
for i in {2015..2100}
do
   InputFile=$DATA_DIR/"$aSMB_BASE_NAME"-"$i".nc
   if [[ ! -f "$InputFile" ]]; then
    echo "$InputFile does not exist; aborting"
    exit
   fi
   echo '***************************************************************************'
   echo '***************************************************************************'
   echo '***  remap file : ' $InputFile
   echo '***************************************************************************'
   echo '***************************************************************************'
   OutPutFile=$PROJ_DIR/"$aSMB_BASE_NAME"-UNST-"$i".nc

   cdo remap,$ELMER_GRID,PROJ_WEIGHTS.nc -selname,aSMB -setgrid,$STRUCT_GRID $InputFile $OutPutFile
done
# Concatenate yearly files
ncrcat PROJ_DATA/aSMB_MARv3.9-yearly-MIROC5-rcp85-UNST-2* PROJ_DATA/aSMB_MARv3.9-yearly-MIROC5-rcp85-UNST-2015-2100.nc
