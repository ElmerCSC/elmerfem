% The end moment / the maximum moment M (M = 50*pi/3)

f = [0.05:0.05:1];

% Reference results from the paper Sze KY, Liu XH, Lo SH. Popular benchmark problems for 
% geometric nonlinear analysis of shells. Finite Elements in Analysis and Design 
% 2004, 40(11):1551-1569.

refres = [0.05 0.196 1.870
0.1 0.774 3.648
0.15 1.699 5.248
0.2 2.918 6.598
0.25 4.361 7.639
0.3 5.945 8.333
0.35 7.585 8.664 
0.4 9.194 8.637
0.45 10.688 8.281
0.5 12.000 7.639
0.55 13.073 6.775
0.6 13.871 5.758
0.65 14.377 4.665 
0.7 14.595 3.571
0.75 14.546 2.546
0.8 14.270 1.650
0.85 13.818 0.926
0.9 13.247 0.405
0.95 12.621 0.098
1.0 12.000 0.0];

load cantilever.dat;
n = size(cantilever,1);
plot(-cantilever(1:n,1),f(1:n),'ro');
hold on;
plot(cantilever(1:n,3),f(1:n),'bo');
plot(refres(:,2),refres(:,1),'r')
plot(refres(:,3),refres(:,1),'b')
xlabel('Displacement')
ylabel('Moment/M_{max}')
title('Cantilever Benchmark')
