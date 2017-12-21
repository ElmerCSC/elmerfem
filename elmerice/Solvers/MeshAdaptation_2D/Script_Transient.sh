#!/bin/bash 
#####################################################################################
#####################################################################################
#   TRANSIENT FIXED POINT MESH ADAPATATION ALGORITHM:
#      Alauzet, F., Frey, P.J., George, P.L., Mohammadi, B., 2007. 3D transient fixed point mesh adaptation for time-dependent problems: Application to CFD simulations. Journal of Computational Physics 222, 592–623. https://doi.org/10.1016/j.jcp.2006.08.012
### 
### IMPLEMENTATION BY:
###   F. GILLET-CHAULET, IGE, Grenoble; Nov. 2017
###   Doc: http://elmerice.elmerfem.org/wiki/doku.php?id=mesh:meshadaptation
#####################################################################################
display_usage() { 
  echo "TRANSIENT FIXED POINT MESH ADAPATATION ALGORITHM:"
  echo "  Provide/update template files:"
  echo "   - RUN_INIT.sif: Initialisation run"
  echo "   - RUN_I_J.sif : Configuration file for the physical simulation" 
  echo "   - MESH_OPTIM_I_J.sif : Configuration file for the mesh adaptation" 
  echo 
} 
if [[ ( "$@" == "--help") ||  "$@" == "-h" ]] 
then 
  display_usage
  exit 
fi 
#####################################################################################
#####################################################################################
ELMERSOLVER=ElmerSolver
##
TEMPLATES=SIFs
##
## TEST THAT TEMPLATE FILES EXIST
if [ ! -d "$TEMPLATES" ]; then
   echo "################################"
   echo "$TEMPLATES directory not found" 
   echo "################################"
   echo
   display_usage
   exit
else
  if [ ! -f "$TEMPLATES/RUN_INIT.sif" ]||[ ! -f "$TEMPLATES/RUN_I_J.sif" ]||[ ! -f "$TEMPLATES/MESH_OPTIM_I_J.sif" ]; then
   echo "################################"
   echo "Template files not found" 
   echo "################################"
   echo
   display_usage
   exit
  fi
fi
#####################################################################################
#####################################################################################
### I- SIMULATIONS PARAMETERS
dt=0.0005 # time step size
time_intervals=10 # number of time steps
output_interval=5 # output intervals 

### II- ALGORITHM PARAMETERS
##   1- PHYSICAL LOOP
imin=0
imax=0
##   2- REMESHING LOOP
jmax=1 # maxium number of remeshing loops

### III- CONVERGENCE TEST (Here compute the relative node number difference between 2 remeshing loops) 
converged (){
tol=0.05
if [ $2 -eq 0 ]
then
  return 1
else
  mesh1=Mesh_I"$1"_J"$2"
  jj=$(($2-1))
  mesh2=Mesh_I"$1"_J"$jj"
  result=$(echo $(head -n 1 $mesh1/mesh.header | awk '{print $1}')  $(head -n 1 $mesh2/mesh.header | awk '{print $1}') | awk '{print 2*sqrt(($1-$2)*($1-$2))/($1+$2)}')
  var=$(awk 'BEGIN{print "'$result'"<"'$tol'"?1:0}')
  if [ "$var" -eq 1 ];then
     return 0
  else
     return 1
  fi
fi
}
#####################################################################################
#####################################################################################
## COMPILATION
#make
##

### I- PHYSICAL LOOP:  i
for ((i=$imin ; i<=$imax ; i++))
do
  if [ $i -eq 0 ]
  then
  ## INITIALISATION RUN
    cp $TEMPLATES/RUN_INIT.sif .
    echo RUN_INIT.sif > ELMERSOLVER_STARTINFO
    $ELMERSOLVER 
  fi

  ### II- MESH ADAPTATION LOOP: j
  for ((j=0 ; j<=$jmax ; j++))
  do
   ##########################
   #### PHYSICAL SIMULATION
   ##########################
       echo 
       echo "----------------------------------------------------"
       echo "----------------------------------------------------"
       echo " Mesh Adapation steps (i,j):"
       echo $i $j
       echo "----------------------------------------------------"
       echo "----------------------------------------------------"
       echo 
     ################################  
     ### PARAMETERS IN PHYSICAL SIF
       #!!PARAMETERS:
       #!! <I>
       #!! <J>

       #!! <O_I>
       #!! <T_I>
       #!! <dT>
     sed  "s/<I>/$i/g;s/<J>/$j/g;s/<O_I>/$output_interval/g;s/<T_I>/$time_intervals/g;s/<dT>/$dt/g" $TEMPLATES/RUN_I_J.sif  > Run_"$i"_"$j".sif
     ################################  
     echo Run_"$i"_"$j".sif > ELMERSOLVER_STARTINFO
     $ELMERSOLVER

   #################################
   #### MESH ADAPTATION SIMULATION
   ################################
     ################################  
     ### PARAMETERS IN PHYSICAL SIF
       #!!PARAMETERS:
       #!! <INTERVALS>
       #!! <ExecRELOAD>
       #!! <SAVE>
     if ( [ $j -eq $jmax ] || converged $i $j )
     # MESH ADAPTATION HAS FINISHED
     then
        # PHYSICAL SIMULATION HAS FINISHED
        if [ $i -eq $imax ];then exit;fi 
	# i=i+1; only read last time step for mesh adaptation
        INTERVALS=1
        Exec="never"
        SAVE=I"$((i+1))"_J0
     else
     # j=j+1; read all time steps for the mesh adaptation
         INTERVALS=$((time_intervals/output_interval+1))
         Exec="always"
         SAVE=I"$i"_J"$((j+1))"
     fi
     ################################  
     sed  "s/<I>/$i/g;s/<J>/$j/g;s/<INTERVALS>/$INTERVALS/g;s/<ExecRELOAD>/$Exec/g;s/<SAVE>/$SAVE/g" $TEMPLATES/MESH_OPTIM_I_J.sif > MESH_OPTIM_"$i"_"$j".sif
     echo MESH_OPTIM_"$i"_"$j".sif > ELMERSOLVER_STARTINFO
     $ELMERSOLVER

  ### END j-LOOP
  done

## END i-LOOP
done
