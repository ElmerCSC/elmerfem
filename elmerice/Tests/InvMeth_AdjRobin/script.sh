#!/bin/sh
#

#Post-procesing on files to obtain file:  results.dat
`tail -5 Cost_Robin_Beta.dat | awk '{print$1}' > foo1.dat`
#Recuperation de la norme de la vitesse: (nrm iter 2 de chaque time step)
`grep "ComputeChange: NS (ITER=2) (NRM,RELC): (" OutputSIF_InvMeth_AdjRobin.log | grep "navier-stokes" | awk '{print$6}'  > NRMvel.dat`
`paste foo1.dat NRMvel.dat> results.dat`
`rm foo1.dat NRMvel.dat`

#Creation of file:  results.dat.names
`echo   1: value: time scalar variable > results.dat.names`
`echo   2: NRM: velocity >> results.dat.names`
