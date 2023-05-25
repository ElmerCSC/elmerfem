% Shear force

f = [0.2:0.2:4];

% Reference results from the paper Sze KY, Liu XH, Lo SH. Popular benchmark problems for 
% geometric nonlinear analysis of shells. Finite Elements in Analysis and Design 
% 2004, 40(11):1551-1569. The first column gives the shear force/4.

refres = [0.05 0.026 0.663
0.1 0.103 1.309
0.15 0.224 1.922
0.2 0.381 2.493
0.25 0.563 3.015
0.3 0.763 3.488
0.35 0.971 3.912
0.4 1.184 4.292
0.45 1.396 4.631
0.5 1.604 4.933
0.55 1.807 5.202
0.6 2.002 5.444
0.65 2.190 5.660
0.7 2.370 5.855
0.75 2.541 6.031
0.8 2.705 6.190
0.85 2.861 6.335
0.9 3.010 6.467
0.95 3.151 6.588
1.0 3.286 6.698];

load cantilever.dat;
plot(-cantilever(:,1),f(:),'r');
hold on;
plot(cantilever(:,3),f(:),'b');
plot(refres(:,2),refres(:,1)*4,'ro')
plot(refres(:,3),refres(:,1)*4,'bo')
xlabel('Displacement')
ylabel('Force')
title('Cantilever Benchmark')
