# User Functions
This Section presents the various user functions developed to solve glaciological applications using the finite element code Elmer. User functions are used to define boundary conditions, initial conditions or any material parameter which depends on other variables.

User functions in the Elmer/Ice package are automatically compiled and can be called using the object file name ElmerIceUSF:
```
Temperature = Variable Coordinate 3 
 Real Procedure  "ElmerIceUSF" "NameUSF"
```
where *NameUSF* is the name of the user function you want to execute. The source code of the user functions of the Elmer/Ice package can be obtained from the Elmer svn in /trunk/elmerice/UserFunctions/.

A non-distributed user function has to be link to Elmer:

`elmerf90 MyUSF.f90 -o MyUSF`
In the SIF file, the call to the user function is done by pointing on the object file MyUSF, like this

```
Temperature = Variable Coordinate 3 
 Real Procedure "MyUSF" "NameUSF"
```
where *NameUSF* is the name of the fortran function in the file MyUSF.f90.

Instead of using a fortran user function, one can write directly in the SIF file a MATC function, as explained here. However, MATC functions present some limitations and cannot be used in all cases. Moreover, MATC functions are less computationally efficient, and should be avoided for very large problems.
