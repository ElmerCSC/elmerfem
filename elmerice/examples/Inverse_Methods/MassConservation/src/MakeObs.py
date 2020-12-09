###############################################################
## Generate synthetic thickness observations along lines at given 
#    x positions
###############################################################
import numpy as np

## parameters H = Hgl + dhdx * x
Hgl=400.0
dhdx=-1.0e-3

## mu and sigma to generate noise from a normal distribution
mu=0.0
sigma=20.0

## the ice true thickness
def H(x):
    z = Hgl + dhdx*x
    return z

## make the obs
def MakeObs(x,y):

    nn=len(x)*len(y)
    r = np.random.normal(mu, sigma, nn)

    with open('Hflight.txt', 'w') as f:
        for i in range(len(x)):
            for j in range(len(y)):
                h=H(x[i]) + r[i*len(y)+j]       
                f.write("{} {} {}\n".format(x[i],y[j],h))
    f.close()

if __name__ == "__main__":

    # generate 3 lines at x=
    x=np.linspace(50e3,150e3,3)

    # resolution along y
    y=np.linspace(0.0,50e3,51)

    MakeObs(x,y)

