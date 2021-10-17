% The point load of the sif file, 
% multiply with 2 to get the load for the problem without symmetry conditions 
f_elmer = [10:10:200];

load hemisphere.dat;
plot(hemisphere(:,1),2*f_elmer(:)/100,'r');
hold on;
plot(-hemisphere(:,11),2*f_elmer(:)/100,'b');

xlabel('Displacement')
ylabel('Force/100')
title('Hemisphere Benchmark')



% Reference results from the paper Sze KY, Liu XH, Lo SH. Popular benchmark problems for 
% geometric nonlinear analysis of shells. Finite Elements in Analysis and Design 
% 2004, 40(11):1551-1569. 

%refres = [0.05 0.855 0.955
%0.1 1.499 1.840
%0.15 1.969 2.604
%0.2 2.321 3.261];

%plot(refres(:,2),refres(:,1)*4,'ro')
%hold on
%plot(refres(:,3),refres(:,1)*4,'bx')
