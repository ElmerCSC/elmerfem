Solve linearized 2-D shell equations of Reissner-Naghdi type in the case of 
a straight cylindrical shell with membrane-dominated asymptotic behaviour.
The problem has been described in Pitkäranta et al. Shell deformation states 
and the finite element method: a benchmark study of cylindrical shells. 
Computer Methods in Applied Mechanics and Engineering 1995. 128:81-121. Here
the shell equations are written in lines of curvature coordinates (y1 = 
angular direction, y2 = axial direction) and a special solver code is
employed to handle the problem using the exact surface parametrization.

This test provides a case for verification of the p-version of FEM.
The energy norm of the error for different p-element definitions (over the
same 2 X 2 mesh) is found to be as follows (the shell thickness d=0.1): 

p:1       Relative energy error =   0.88896260271790639
p:2       Relative energy error =   0.29814339193540801 
p:2 b:1   Relative energy error =   0.23711572317485702
p:2 b:15  Relative energy error =   0.10031149151079155
p:3       Relative energy error =   0.17701104808582918
p:3 b:15  Relative energy error =    1.9351692322981600E-002
p:4       Relative energy error =    3.4526616304186068E-002
p:4 b:15  Relative energy error =    7.0424540514543998E-003
p:5       Relative energy error =    7.8498375202090825E-003
p:5 b:15  Relative energy error =    2.0809479830700603E-003
p:6       Relative energy error =    2.4502183740983921E-003
p:6 b:15  Relative energy error =    6.2443241030678861E-004
p:7       Relative energy error =    6.0846479859093326E-004
p:7 b:15  Relative energy error =    1.6535890420287850E-004
p:8       Relative energy error =    1.3328077843970739E-004
p:8 b:15  Relative energy error =    1.3328077843970739E-004  
p:9       Relative energy error =    6.7817426835655489E-005

Note that, except for the case p=2, the current version of Elmer doesn't 
offer a way to obtain an approximation based on the space Q_{p} of polynomials
of degree at most p in each variable separately, since Elmer doesn't employ
the bubbles corresponding to Q_{p}. For example the p-approximation space 
S_{p=8} of Elmer has 15 bubble functions and includes only Q_{4}. 
A deterioration of accuracy may thus happen over non-affine finite elements. 
The results above also show the error when the bubble augmentation 
corresponding S_{p=8} is used in connection with p-elements of lower order. 

For a very high p errors from floating-point arithmetic may limit the 
obtainable accuracy.
