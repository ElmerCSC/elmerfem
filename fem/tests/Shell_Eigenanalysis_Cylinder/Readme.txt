This case gives a FE model for a fan blade problem considered in

Olson M, Lindberg G: Vibration analysis of cantilevered curved plates using
a new cylindrical shell finite element. A report of Air Force Flight Dynamics
Laboratory, AFFDL-TR-68-150, 1968.

The report lists the first natural frequencies that were measured 
experimentally. The table below shows the first seven vibration frequencies 
obtained computationally with this FE model and the experimental values of
the report. Note that Elmer outputs true eigenvalues, so to get the natural 
frequency one needs to compute f = sqrt(lambda)/(2 pi), with lambda being 
the eigenvalue. 

Mode  FE model    Experimental  
1     95.3        86.6
2     145.4       135.5
3     251.7       258.9
4     362.5       350.6
5     404.7       395.2
6     535.6       531.1
7     752.6       743.2
