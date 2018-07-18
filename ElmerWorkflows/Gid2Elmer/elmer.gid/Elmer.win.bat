rem set basename = %1   
rem set directory = %2  
rem set ProblemDirectory = %3

cd %2

%3\gid2elmer.exe < %1.dat

%3\findparents.exe

del %2\mesh.boundary

ren %2\mesh.boundary.corrected mesh.boundary
