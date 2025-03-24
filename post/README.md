# ElmerPost
ElmerPost has been moved to a separate repository.
https://github.com/ElmerCSC/ElmerPost


Below is a copy-paste from there:
--

Repository for the old preprocessor of Elmer FEM suita, ElmerPost

ElmerPost is the original finite element postprocessor for Elmer FEM suite that was
written starting from 1995 by Juha Ruokolainen. It was used as the main postprocessing tool
until the improved support of external postprocessors gradually made it loose its foothold.
Currently Paraview is the most popular postprocessing tool for Elmer FEM suite. 

Bcause the code is more or less obsolite it was moved to a separate repository. If you have an old compilation system the subdirectory "post" here is one-to-one copy from the main repository. Hence you can make a symlink to that and everything in your compilation script should work as before.

ElmerPost is written in C with Tcl/Tk graphical user interface toolkit.
Also ElmerPost is integarated with MATC, a mathematical scripting utility, also written by Juha Ruokolainen.
In total ElmerPost has around 100,000 lines of code.