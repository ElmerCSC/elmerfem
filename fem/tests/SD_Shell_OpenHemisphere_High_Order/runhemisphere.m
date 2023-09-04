% The point load of the sif file, 
% multiply with 2 to get the load for the problem without symmetry conditions 
f_elmer = [10:10:200];

% The index must be chosen according to the number of the field variables
% (ind=11 for DOFs==6 and ind=14 for DOFs==9) 
ind = 14; 

load hemisphere.dat;
plot(hemisphere(:,1),2*f_elmer(:)/100,'r');
hold on;
plot(-hemisphere(:,ind),2*f_elmer(:)/100,'b');

xlabel('Displacement')
ylabel('Force/100')
title('Hemisphere Benchmark')


% Reference results from the paper Sze KY, Liu XH, Lo SH. Popular benchmark problems for 
% geometric nonlinear analysis of shells. Finite Elements in Analysis and Design 
% 2004, 40(11):1551-1569. 

refres = [0.05 0.855 0.955
0.1 1.499 1.840
0.15 1.969 2.604
0.2 2.321 3.261
0.25 2.596 3.833
0.30 2.819 4.339
0.35 3.002 4.790
0.40 3.158 5.196
0.45 3.291 5.565
0.50 3.406 5.902
0.55 3.508 6.212
0.60 3.598 6.497
0.65 3.678 6.761
0.70 3.750 7.006
0.75 3.816 7.234
0.80 3.875 7.448
0.85 3.929 7.647
0.90 3.979 7.835
0.95 4.025 8.011 
1.00 4.067 8.178];

plot(refres(:,2),refres(:,1)*4,'rx')
plot(refres(:,3),refres(:,1)*4,'bx')
