#!/bin/sh
#

#Post-procesing on files to obtain file:  results.dat
`tail -5 Cost_Robin_Beta.dat | awk '{print$0}' > foo1.dat`
`tail -5 gradientnormadjoint_robin_beta.dat | awk '{print$2}' > foo2.dat`
`paste foo1.dat foo2.dat > results.dat`
`rm foo1.dat foo2.dat`

#Creation of file:  results.dat.names
`echo   1: value: time scalar variable > results.dat.names`
`echo   2: >> results.dat.names`
`echo   3: >> results.dat.names`
`echo   4: >> results.dat.names`
`echo   5: >> results.dat.names`