NOTE: The capability to handle elementwise changing p-basis is undergoing
      evolution and the ability to treat "hanging DOFs" on element interfaces
      may not be perfect yet.

Here the principal purpose is to test the ability to use p-basis functions
of varying degree within a single model. The example case is a structural
shell model. This case has been derived from the test 

               ..tests/Classical2DShell/

Linearized 2-D shell equations of Reissner-Naghdi type are solved in the case
of a straight cylindrical shell with membrane-dominated asymptotic behaviour.
The problem has been described in Pitkäranta et al. Shell deformation states 
and the finite element method: a benchmark study of cylindrical shells. 
Computer Methods in Applied Mechanics and Engineering 1995. 128:81-121. Here
the shell equations are written in lines of curvature coordinates (y1 = 
angular direction, y2 = axial direction) and a special solver code is
employed to handle the problem using the exact surface parametrization. The
solver code outputs the energy norm of the error.

The part of the computational domain where a boundary layer has an important
role is associated with an additional body identifier so that a different
degree of approximation can be used within the boundary layer.

The degree of approximation is defined by using a special keyword construct

              Element = "p:%a_matc_function"

where the part "a_matc_function" is the name of a MATC function of a 4-tuple.
The first argument is the identifier of the body while the remaining components
are the coordinates of the element mid-point. At the time of writing, this
seems to be the only functional way to get varying p within a single model.

The energy norm of the error for different p-element definitions (over the
same 2 X 2 mesh, with p being fixed as p=8 in the boundary layer and p varying
elsewhere) is found to be as follows (the shell thickness d=0.1): 

p=1  Relative energy error =   0.85396972358657408
p=2  Relative energy error =   0.23800207437997356
p=3  Relative energy error =   0.14558431573734487
p=4  Relative energy error =   2.8912002775000799E-002
p=5  Relative energy error =   7.0313590547632018E-003
p=6  Relative energy error =   2.3734217622063724E-003
p=7  Relative energy error =   5.8687216090359958E-004
p=8  Relative energy error =   1.3339512982768768E-004
