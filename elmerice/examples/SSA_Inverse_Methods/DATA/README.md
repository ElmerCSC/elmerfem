# SSA inverse methods test cases

This directory contains data required to run the examples:

- MacAyeal_VELOCITIES.txt: x,y,u,v

   - contains position and velocity for the Mac Ayeal test case that have been produce with the model and no noise added.

- AddNoise.py: Python script to add random gaussian noise to MacAyeal_VELOCITIES.txt with mean 0 and standard deviation STD

   - Edit and chande STD value.
   - result stored in MacAyeal_VELOCITIES_NOISE.txt

- PlotHist.py: python script to plot histograms of velocity differences between MacAyeal_VELOCITIES.txt and MacAyeal_VELOCITIES_NOISE.txt


- To run the Ronne_Filchner you will need to get velocity observations:

   - You can download the [MEaSUREs InSAR-Based Antarctica Ice Velocity Map](https://nsidc.org/data/NSIDC-0484/versions/2)
   - You can subset the data file using (i.e. taking data within the domain and taking only fourth point) : 

```     
ncks -d x,-1515000.0,-535000.0,4 -d y,135000.0,1040000.0,4 antarctica_ice_velocity_450m_v2.nc RonneFilchner.nc
```	
