A benchmark problem where an open hemispherical shell is subject to point
loads. The reference results can be found in Sze KY, Liu XH, Lo SH. Popular 
benchmark problems for geometric nonlinear analysis of shells. Finite 
Elements in Analysis and Design 2004, 40(11):1551-1569.

Here the problem is solved by combining a high-order discretization and improved
surface reconstruction done over a standard surface mesh embedded in 3-D space.
Having a p-element definition automatically switches to a special formulation 
(sets "Cartesian Formulation = Logical True"). 

The following deflections at the two points where the loads are applied are
obtained for different element definitions:

    2.523397056307E-003  p:1
    1.822842204650E-001  p:2
    1.823444430702E-001  p:2 b:1
    1.848458106450E-001  p:3
    1.848616269838E-001  p:3 b:1
    1.849267786277E-001  p:3 b:6
    1.849362922672E-001  p:4

   -2.524376505111E-003  p:1
   -1.863177044741E-001  p:2
   -1.863777270085E-001  p:2 b:1
   -1.889294940331E-001  p:3
   -1.889448030200E-001  p:3 b:1
   -1.890120195621E-001  p:3 b:6
   -1.890220598746E-001  p:4

Note: With "p:3 b:6" the basis of Q_3 is included
