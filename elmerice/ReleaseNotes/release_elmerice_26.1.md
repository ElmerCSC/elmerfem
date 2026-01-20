Elmer/Ice Release Notes for version 26.1
========================================

Previous release: **9.0**  
Period covered: **Nov 11, 2020 - Jan 19, 2025**  

Trying to exclude the for Elmer/Ice relevant commits is difficult. The list obtained with
```bash
git log --since="2020-11-11" --stat -- elmerice
```
revealed 371 commits since Nov. 11th, 2020.

Apart from the core Elmer team at CSC (Juhani K., Mika M., Juha R., Peter R., Thomas Z.) git log shows contributions from Markus Mützel, Fabien G-C, Iain W.,  Rupert G., Julien B., Samuel C.,  Benjamin R., Mondher C., Olivier G., Joe T.,  Lucas B. to Elmer/Ice related contributions in this release.

### New Versioning Scheme

From this version onward we migrate for Elmer (and hence also Elmer/Ice therein) to Calendar Versioning such that
- First number (major) is the year of the 21st century, e.g. 26
- Second number (minor) is an ordinal number of releases in that year
- Third number (micro) is a growing number which for releases is always 0 and may be omitted.  


I. Overview of changes
---------------------
Just in time to the previous release the Elmer/Ice community introduced the concept of documentation inside this repository, namely,
- `elmerice/Solvers/Documentation` for Solvers
-`elmerice/UserFunctions/Documentation`for user-functions
Since then, these directories have been filled with the documents.



II. New Solver Modules
----------------------

III. Enhancement of existing solvers
------------------------------------

