# Mass Conservation

Validation and test case for the mass conservation method. 

The thickness solver is used to compute a steady-state thickness field using a given 2D velocity field and a mass balance.
The adjoint code of the thickness solver is used to optimise the input fields to match discrete thickness observations.

The example is based on an analytical solution for an ice shelf ramp in 2D, cf *e.g.* Greve, *Dynamics of Ice Sheets and Glaciers*, lectures note 2010 (ch. 5.4):

	- The ice thickness decreases linearly form Hgl in x=0 with a slope dhdx; i.e. H= Hg + dhdx * x

	- We have an inflow velocity in x=0, u=(Vgl,0)

	- there is a calving front in x=Lx

	- there is no friction and no velocicty accross the side boundaries in y=0 and y=Ly, i.e. uy=0

As there is no lateral variation along y, the velocity and SMB should be identical to the 1D analytical solution

## Content

	- **src** : source codes and inputs to run the tests

	- **DirectValidation** : The SSA solver and thickness solver are validated againt the analytical solution

	- **GradientValidation** : Validation of the adjoint based gradient of the steady-state thickness with respect to the velocity field

	- **Optimisation** : Optimise the input velocity field to match manufactured thickness observations. 
