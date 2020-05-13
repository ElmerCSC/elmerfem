import numpy as np
from pylab import genfromtxt;

STD=1.0

fname="MacAyeal_VELOCITIES.txt"
mat = genfromtxt(fname);

npoints=np.size(mat,0)

noiseU = np.random.normal(0, STD, npoints)
noiseV = np.random.normal(0, STD, npoints)

rms=np.sqrt(np.sum((noiseU*noiseU)+(noiseV*noiseV))/npoints)
print(rms)

matn=mat
matn[:,2]=mat[:,2]+noiseU
matn[:,3]=mat[:,3]+noiseV

np.savetxt('MacAyeal_VELOCITIES_NOISE.txt',matn, delimiter=' ') 
