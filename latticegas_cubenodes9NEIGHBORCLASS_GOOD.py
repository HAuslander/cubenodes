#!/cygdrive/c/Python27/python


global FreshProb
global MaxNodes

print("""

# cubnodes



The goal is to create graphs with no definite embeddingthat have a quasi-Sierpinski structure in n dimensions. 
As the graph "grows" by dividing it's interior cubes into an ever finer patchwork mosaic, we arrive at what is errectively and expanding universe (from the perspective of anyone residing within these graphs).

On these graphs, densities of random walks will be shown to satisfy Fick's laws of diffusion and
therefore the heat equation, and therefore the wave equation (by way of Brownian-Huygens propagation),
which lead to a Euclidean embedding arising naturally.

I.e. this is about how to create graphs where the Euclideanness of the 
underlying space is inferred.

This is easily done for 1-dimensional case, since any method of increasing nodes results

In this version, we start with cubes of length 8 in either dimension (toroidally connected)
We will typically divide a cube with TWO cuts allong each axis instead of just bisecting it with one
so as to preserve the parity of each node

ALSO, the main difference between this and an earlier (failed) version is we will allow adjacent 
(i.e. sharing a face or edge) cubes to be divided, which means
matching up previously created mid-nodes (so a cube can only hava fractality order that is
1 greater than surrounding cubes). In the previous version once a cube
was divided, no adjacent cubes could be divided (except the "diagonal" cubes that share a single point).

To distinguish cubes that have a fractality differing by 1 we initally allow each neighbor to have a previous neighbor feature (ie it's a little like a 2nd order
equation that requires two time steps to update)

A subsequent version allows us to create these graphs with no notion of previous neighbors.

STEPS:
        Divide a cube:
        1) pick a node and check that it has 2D neighbors
        2) pick D other neighbors (connected edges of the same color/magnification) that are not collinear and that have 2D neighbors each;
           check to see if they are part of a cube whose edges have the same color/magnification and whose nodes all have 2D neighbors.
        3) remember the slices can vary with  color/magnification; the important thing is that no cube gets left behind or put into a weird
        configuration that disallows it from being divided.

        4) the fact that you know the "intervening" nodes across an edge are actually subdivisions (because they have a higher color/magnification)
        and the fact that all the divided edges (even the cross-hatch ones inside) connect to the larger cube edges means that you can match
        up the nodes


    In this version, we will have 2 probabilities; a really small one and a higher one; the higher one will
    only act on nodes where two or more magnifications intersect (which only happens if one of its axes is
    an edge of a division, so that the nodes neighbors have less than 2D neighbors of their own)

""")

import sys
import os
import optparse
from copy import copy

import numpy as np
import pandas as pd

import os.path
import datetime, time

import glob

import pickle 
#import pandas as pn    

import matplotlib.pyplot as plt      
from matplotlib.backends.backend_pdf import PdfPages

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import random as rn
#from bisect import bisect_left
import bisect
import scipy.stats as stats
from scipy.stats import linregress

#import statgrid2dclassical # don't need this cython routine after all
 

def MonteCarlo(D, R, nruns=15, nstep=5):
    nsuccess = 0
    p = np.ones((D,))/float(D)
    
    sumvar = 0
    NSucc = 0
    for i in range(nruns):
        x = np.zeros((D,))
        pvec = np.ones((D,))/float(D)
        
        if D > 1:
            randvec = rn.multinomial(nstep, p) # rn.normal(size=(D,nstep))/nstep
        else:
            randvec = [nstep]
        for ii,ir in enumerate(randvec):
            myres = rn.multinomial(ir, [0.5, 0.5])
            x[ii] = myres[0] - myres[1]
            #print(x, randvec, myres)
        sumsq = np.sum(x * x)
        sumvar += sumsq
        if sumsq <= R * R:
            NSucc += 1

    print("%f %f %f" % (NSucc/float(nruns), sumvar/float(nruns), np.sqrt(sumvar/float(nruns))))

    return NSucc / float(nruns)
# MonteCarlo(2,1.666,1,1)

#
#D=1; k=1; R=1; t=1; NSteps=100
# PDEGauss(D,k,R,t,NSteps= 100)

from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt


def OLDPDEGauss(D,k,R,t,NSteps= 100):
    # remember sigma^2 = 2kt

    if D == 1:
        grid = np.zeros((NSteps*2 + 10))
        midpt = grid.shape[0]//2
        grid[midpt] = 1.0


        
        for istep in range(1,NSteps+1):
            chunk = copy(grid[(midpt-istep):(midpt+istep)])
            grid[(midpt-1-istep):(midpt-1+istep)] = grid[(midpt-1-istep):(midpt-1+istep)] + 0.5 * chunk
            grid[(midpt+1-istep):(midpt+1+istep)] = grid[(midpt+1-istep):(midpt+1+istep)] + 0.5 * chunk
            
            grid[istep % 2::2] = 0

        sumvar = 0
        for i in range(grid.shape[0]):
            sumvar += (i-midpt) * (i-midpt) * grid[i]
        print("var %f sd  %f" % (sumvar, np.sqrt(sumvar)))
        
        Radj = np.round(R*np.sqrt(NSteps))
        lo = int(np.max([0, midpt-Radj]))
        hi = int(np.min([grid.shape[0], midpt+Radj]))



        return np.sum( grid[lo : hi ]), np.sum(grid)
    

    if D == 2:
        grid = np.zeros((NSteps*2 + 10, NSteps*2 + 10))
        midpt = grid.shape[0]//2
        grid[midpt, midpt] = 1.0

        for istep in range(1,NSteps+1):
            gridcopy = copy(grid)
            for i in range(midpt-istep,midpt+istep+1):
                for j in range(midpt-istep,midpt+istep+1):
                    if (istep+i+j) % 2 == 1:
                        grid[i, j] = 0
                    else:
                        grid[i, j] = grid[i, j] + 0.25 *(gridcopy[i-1,j] + gridcopy[i+1,j] + gridcopy[i,j-1] + gridcopy[i,j+1]) 



            #chunk = copy(grid[(midpt-istep):(midpt+istep+1),(midpt-istep):(midpt+istep+1)])
            #grid[(midpt-1-istep):(midpt-1+istep+1),(midpt-istep):(midpt+istep+1)] = grid[(midpt-1-istep):(midpt-1+istep+1),(midpt-istep):(midpt+istep+1)] + 0.25 * chunk
            #grid[(midpt+1-istep):(midpt+1+istep+1),(midpt-istep):(midpt+istep+1)] = grid[(midpt+1-istep):(midpt+1+istep+1),(midpt-istep):(midpt+istep+1)]+ 0.25 * chunk
            #grid[(midpt-istep):(midpt+istep+1),(midpt-1-istep):(midpt-1+istep+1)] = grid[(midpt-istep):(midpt+istep+1),(midpt-1-istep):(midpt-1+istep+1)] + 0.25 * chunk
            #grid[(midpt-istep):(midpt+istep+1),(midpt+1-istep):(midpt+1+istep+1)] = grid[(midpt-istep):(midpt+istep+1),(midpt-1-istep):(midpt-1+istep+1)] + 0.25 * chunk
            
            #grid[istep % 2::2, istep % 2::2] = 0
            #grid[1+istep % 2::2, 1+istep % 2::2] = 0
            
        Radj = np.round(D * R * np.sqrt(NSteps))
        lo = int(np.max([0, midpt-Radj]))
        hi = int(np.min([grid.shape[0], midpt+Radj]))


        sumvar = 0
        for i in range(grid.shape[0]):
            for j in range(grid.shape[1]):
                sumvar += ((i-midpt) * (i-midpt) + (j-midpt) * (j-midpt)) * grid[i,j]
        print("var %f sd  %f %f" % (sumvar, np.sqrt(sumvar), np.sum(np.sum(grid))))
        bGraph= False
        if bGraph:
            x = np.arange(grid.shape[0]) # np.outer(np.linspace(-2, 2, 10), np.ones(10))
            y = np.arange(grid.shape[0]) # x.copy().T
            z = grid[x,y] # np.cos(x ** 2 + y ** 3)
            
            fig = plt.figure()
            
            # syntax for 3-D plotting
            ax = plt.axes(projection='3d')
            
            # syntax for plotting
            ax.plot_surface(x, y, z, cmap='viridis',\
                            edgecolor='green')
            ax.set_title('Surface plot geeks for geeks')
            plt.show()
        return np.sum( grid[lo:hi, lo:hi] )


# cutoff is a list of increasing
def OLDPDEGauss2(D, NSteps = 10):
    # remember sigma^2 = 2kt
    orighist = []
    if D == 1:
        grid = np.zeros((NSteps*2 + 10))
        midpt = grid.shape[0]//2
        grid[midpt] = 1.0


        
        for istep in range(1,NSteps+1):
            chunk = copy(grid[(midpt-istep):(midpt+istep)])
            
            grid[(midpt-1-istep):(midpt-1+istep)] = grid[(midpt-1-istep):(midpt-1+istep)] + 0.5 * chunk
            grid[(midpt+1-istep):(midpt+1+istep)] = grid[(midpt+1-istep):(midpt+1+istep)] + 0.5 * chunk
            
            grid[istep % 2::2] = 0

        FrontBndry = midpt
        BackBndry = midpt




        return np.sum( grid[lo : hi ]), np.sum(grid)
    

    if D == 2:
        grid = np.zeros((NSteps*2 + 10, NSteps*2 + 10))
        midpt = grid.shape[0]//2
        grid[midpt, midpt] = 1.0

        for istep in range(1,NSteps+1):
            gridcopy = copy(grid)
            for i in range(midpt-istep,midpt+istep+1):
                for j in range(midpt-istep,midpt+istep+1):
                    if (istep+i+j) % 2 == 1:
                        grid[i, j] = 0
                    else:
                        grid[i, j] = grid[i, j] + 0.25 *(gridcopy[i-1,j] + gridcopy[i+1,j] + gridcopy[i,j-1] + gridcopy[i,j+1]) 



            #chunk = copy(grid[(midpt-istep):(midpt+istep+1),(midpt-istep):(midpt+istep+1)])
            #grid[(midpt-1-istep):(midpt-1+istep+1),(midpt-istep):(midpt+istep+1)] = grid[(midpt-1-istep):(midpt-1+istep+1),(midpt-istep):(midpt+istep+1)] + 0.25 * chunk
            #grid[(midpt+1-istep):(midpt+1+istep+1),(midpt-istep):(midpt+istep+1)] = grid[(midpt+1-istep):(midpt+1+istep+1),(midpt-istep):(midpt+istep+1)]+ 0.25 * chunk
            #grid[(midpt-istep):(midpt+istep+1),(midpt-1-istep):(midpt-1+istep+1)] = grid[(midpt-istep):(midpt+istep+1),(midpt-1-istep):(midpt-1+istep+1)] + 0.25 * chunk
            #grid[(midpt-istep):(midpt+istep+1),(midpt+1-istep):(midpt+1+istep+1)] = grid[(midpt-istep):(midpt+istep+1),(midpt-1-istep):(midpt-1+istep+1)] + 0.25 * chunk
            
            #grid[istep % 2::2, istep % 2::2] = 0
            #grid[1+istep % 2::2, 1+istep % 2::2] = 0
            
        Radj = np.round(D * R * np.sqrt(NSteps))
        lo = int(np.max([0, midpt-Radj]))
        hi = int(np.min([grid.shape[0], midpt+Radj]))


        sumvar = 0
        for i in range(grid.shape[0]):
            for j in range(grid.shape[1]):
                sumvar += ((i-midpt) * (i-midpt) + (j-midpt) * (j-midpt)) * grid[i,j]
        print("var %f sd  %f %f" % (sumvar, np.sqrt(sumvar), np.sum(np.sum(grid))))
        bGraph= False
        if bGraph:
            x = np.arange(grid.shape[0]) # np.outer(np.linspace(-2, 2, 10), np.ones(10))
            y = np.arange(grid.shape[0]) # x.copy().T
            z = grid[x,y] # np.cos(x ** 2 + y ** 3)
            
            fig = plt.figure()
            
            # syntax for 3-D plotting
            ax = plt.axes(projection='3d')
            
            # syntax for plotting
            ax.plot_surface(x, y, z, cmap='viridis',\
                            edgecolor='green')
            ax.set_title('Surface plot geeks for geeks')
            plt.show()
        return np.sum( grid[lo:hi, lo:hi] )


def PDEGauss(D,R,k=1,t=1,NSteps=100):
    # remember sigma^2 = 2kt

    orighist = []
    
    if D == 1:
        grid = np.zeros((NSteps*2 + 10))
        midpt = grid.shape[0]//2
        grid[midpt] = 1.0


        
        for istep in range(1,NSteps+1):
            chunk = copy(grid[(midpt-istep):(midpt+istep)])
            grid[(midpt-1-istep):(midpt-1+istep)] = grid[(midpt-1-istep):(midpt-1+istep)] + 0.5 * chunk
            grid[(midpt+1-istep):(midpt+1+istep)] = grid[(midpt+1-istep):(midpt+1+istep)] + 0.5 * chunk
            
            grid[istep % 2::2] = 0
            if (istep + midpt) % 2 != 0:
                orighist.append(grid[midpt])

        sumvar = 0
        for i in range(grid.shape[0]):
            sumvar += (i-midpt) * (i-midpt) * grid[i]
        print("var %f sd  %f" % (sumvar, np.sqrt(sumvar)))
        
        Radj = np.round(R*np.sqrt(NSteps)) / np.sqrt(2*k*t)
        lo = int(np.max([0, midpt-Radj]))
        hi = int(np.min([grid.shape[0], midpt+Radj]))

        return orighist, np.sum( grid[lo : hi ]), np.sum(grid)

    if D == 2:
        grid = np.zeros((NSteps*2 + 6, NSteps*2 + 6))
        midpt = grid.shape[0]//2
        grid[midpt, midpt] = 1.0

        for istep in range(1,NSteps+1):
            gridcopy = copy(grid)
            for i in range(midpt-istep,midpt+istep+1):
                for j in range(midpt-istep,midpt+istep+1):
                    if (istep+i+j) % 2 == 1:
                        grid[i, j] = 0
                    else:
                        grid[i, j] = grid[i, j] + 1/4.0*(gridcopy[i-1,j] + gridcopy[i+1,j] + gridcopy[i,j-1] + gridcopy[i,j+1])

            if (istep + D*midpt) % 2 == 0:
                orighist.append(grid[midpt,midpt])



            
        Radj = np.round(R*np.sqrt(NSteps)) / np.sqrt(2*k*t)
        lo = int(np.max([0, midpt-Radj]))
        hi = int(np.min([grid.shape[0], midpt+Radj]))


        sumvar = 0
        for i in range(grid.shape[0]):
            for j in range(grid.shape[1]):
                sumvar += ((i-midpt) * (i-midpt) + (j-midpt) * (j-midpt)) * grid[i,j]

        retval = 0
        for i in range(grid.shape[0]):
            for j in range(grid.shape[1]):
                if ((i-midpt) * (i-midpt) + (j-midpt) * (j-midpt)) <= R * R:
                    retval += grid[i,j]
        print("var %f sd  %f %f" % (sumvar, np.sqrt(sumvar), np.sum(np.sum(grid))))
        return orighist, np.sum( grid[lo:hi, lo:hi] )

    

    if D == 3:
        grid = np.zeros((NSteps*2 + 6, NSteps*2 + 6, NSteps*2 + 6))
        midpt = grid.shape[0]//2
        grid[midpt, midpt, midpt] = 1.0

        for istep in range(1,NSteps+1):
            gridcopy = copy(grid)
            for i in range(midpt-istep,midpt+istep+1):
                for j in range(midpt-istep,midpt+istep+1):
                    for k in range(midpt-istep,midpt+istep+1):
                        if (istep+i+j+k) % 2 == 1:
                            grid[i, j, k] = 0
                        else:
                            grid[i, j, k] = grid[i, j, k] + 1/6.0*(gridcopy[i-1,j,k] + gridcopy[i+1,j,k] + gridcopy[i,j-1,k] + gridcopy[i,j+1,k] + gridcopy[i,j,k-1] + gridcopy[i,j,k+1])

            if (istep + D*midpt) % 2 == 0:
                orighist.append(grid[midpt,midpt,midpt])
                print(grid[midpt,midpt,midpt])


            #chunk = copy(grid[(midpt-istep):(midpt+istep+1),(midpt-istep):(midpt+istep+1)])
            #grid[(midpt-1-istep):(midpt-1+istep+1),(midpt-istep):(midpt+istep+1)] = grid[(midpt-1-istep):(midpt-1+istep+1),(midpt-istep):(midpt+istep+1)] + 0.25 * chunk
            #grid[(midpt+1-istep):(midpt+1+istep+1),(midpt-istep):(midpt+istep+1)] = grid[(midpt+1-istep):(midpt+1+istep+1),(midpt-istep):(midpt+istep+1)]+ 0.25 * chunk
            #grid[(midpt-istep):(midpt+istep+1),(midpt-1-istep):(midpt-1+istep+1)] = grid[(midpt-istep):(midpt+istep+1),(midpt-1-istep):(midpt-1+istep+1)] + 0.25 * chunk
            #grid[(midpt-istep):(midpt+istep+1),(midpt+1-istep):(midpt+1+istep+1)] = grid[(midpt-istep):(midpt+istep+1),(midpt-1-istep):(midpt-1+istep+1)] + 0.25 * chunk
            
            #grid[istep % 2::2, istep % 2::2] = 0
            #grid[1+istep % 2::2, 1+istep % 2::2] = 0
            
        Radj = np.round(R*np.sqrt(NSteps)) / np.sqrt(2*k*t)
        lo = int(np.max([0, midpt-Radj]))
        hi = int(np.min([grid.shape[0], midpt+Radj]))


        sumvar = 0
        for i in range(grid.shape[0]):
            for j in range(grid.shape[1]):
                for k in range(grid.shape[2]):
                    sumvar += ((i-midpt) * (i-midpt) + (j-midpt) * (j-midpt) + (k-midpt) * (k-midpt)) * grid[i,j,k]

        retval = 0
        for i in range(grid.shape[0]):
            for j in range(grid.shape[1]):
                for k in range(grid.shape[2]):
                    if ((i-midpt) * (i-midpt) + (j-midpt) * (j-midpt) + (k-midpt) * (k-midpt)) <= R * R:
                        retval += grid[i,j,k]
        print("var %f sd  %f %f" % (sumvar, np.sqrt(sumvar), np.sum(np.sum(np.sum(grid)))))
        return orighist, np.sum( grid[lo:hi, lo:hi, lo:hi] )


    if D == 4:
        grid = np.zeros((NSteps*2 + 6, NSteps*2 + 6, NSteps*2 + 6, NSteps*2 + 6))
        midpt = grid.shape[0]//2
        grid[midpt, midpt, midpt, midpt] = 1.0
        orighist = []

        for istep in range(1,NSteps+1):
            gridcopy = copy(grid)
            for i in range(midpt-istep,midpt+istep+1):
                for j in range(midpt-istep,midpt+istep+1):
                    for k in range(midpt-istep,midpt+istep+1):
                        for m in range(midpt-istep,midpt+istep+1):
                            if (istep+i+j+k+m) % 2 == 1:
                                grid[i, j, k, m] = 0
                            else:
                                grid[i, j, k, m] = grid[i, j, k, m] + 1/8.0*(gridcopy[i-1,j,k,m] + gridcopy[i+1,j,k,m] + gridcopy[i,j-1,k,m] + gridcopy[i,j+1,k,m] + gridcopy[i,j,k-1,m] + gridcopy[i,j,k+1,m] + gridcopy[i,j,k,m-1] + gridcopy[i,j,k,m+1])

            if (istep + D*midpt) % 2 == 0:
                orighist.append(grid[midpt,midpt,midpt,midpt])
                print(grid[midpt,midpt,midpt,midpt])


            #chunk = copy(grid[(midpt-istep):(midpt+istep+1),(midpt-istep):(midpt+istep+1)])
            #grid[(midpt-1-istep):(midpt-1+istep+1),(midpt-istep):(midpt+istep+1)] = grid[(midpt-1-istep):(midpt-1+istep+1),(midpt-istep):(midpt+istep+1)] + 0.25 * chunk
            #grid[(midpt+1-istep):(midpt+1+istep+1),(midpt-istep):(midpt+istep+1)] = grid[(midpt+1-istep):(midpt+1+istep+1),(midpt-istep):(midpt+istep+1)]+ 0.25 * chunk
            #grid[(midpt-istep):(midpt+istep+1),(midpt-1-istep):(midpt-1+istep+1)] = grid[(midpt-istep):(midpt+istep+1),(midpt-1-istep):(midpt-1+istep+1)] + 0.25 * chunk
            #grid[(midpt-istep):(midpt+istep+1),(midpt+1-istep):(midpt+1+istep+1)] = grid[(midpt-istep):(midpt+istep+1),(midpt-1-istep):(midpt-1+istep+1)] + 0.25 * chunk
            
            #grid[istep % 2::2, istep % 2::2] = 0
            #grid[1+istep % 2::2, 1+istep % 2::2] = 0
            
        Radj = np.round(R*np.sqrt(NSteps)) / np.sqrt(2*k*t)
        lo = int(np.max([0, midpt-Radj]))
        hi = int(np.min([grid.shape[0], midpt+Radj]))


        sumvar = 0
        for i in range(grid.shape[0]):
            for j in range(grid.shape[1]):
                for k in range(grid.shape[2]):
                    for m in range(grid.shape[3]):
                        sumvar += ((i-midpt) * (i-midpt) + (j-midpt) * (j-midpt) + (k-midpt) * (k-midpt) + (m-midpt) * (m-midpt)) * grid[i,j,k,m]

        retval = 0
        for i in range(grid.shape[0]):
            for j in range(grid.shape[1]):
                for k in range(grid.shape[2]):
                    for m in range(grid.shape[3]):
                        if ((i-midpt) * (i-midpt) + (j-midpt) * (j-midpt) + (k-midpt) * (k-midpt) + (m-midpt) * (m-midpt)) <= R * R:
                            retval += grid[i,j,k,m]
        print("var %f sd  %f %f" % (sumvar, np.sqrt(sumvar), np.sum(np.sum(np.sum(np.sum(grid))))))
        return orighist, retval, np.sum( grid[lo:hi, lo:hi, lo:hi, lo:hi] )






def PDETouch(D,R,k=1,t=1,NSteps=100):
    # remember sigma^2 = 2kt

    orighist = []
    
    if D == 1:
        grid = np.zeros((NSteps*2 + 10))
        midpt = grid.shape[0]//2
        grid[midpt] = 1.0


        ngrid = grid.shape[0]
        for istep in range(1,NSteps+1):
            for i in range(ngrid):
                if grid[(i+1) % ngrid] != 0 or  grid[i-1] != 0:
                    grid[i] = 1.0
            print(np.sum(grid))

        
        Radj = np.round(R*np.sqrt(NSteps)) / np.sqrt(2*k*t)
        lo = int(np.max([0, midpt-Radj]))
        hi = int(np.min([grid.shape[0], midpt+Radj]))


        return orighist, np.sum( grid[lo : hi ]), np.sum(grid)

    if D == 2:
        grid = np.zeros((NSteps*2 + 6, NSteps*2 + 6))
        midpt = grid.shape[0]//2
        grid[midpt, midpt] = 1.0

        for istep in range(1,NSteps+1):
            gridcopy = copy(grid)
            for i in range(midpt-istep,midpt+istep+1):
                for j in range(midpt-istep,midpt+istep+1):
                    if (istep+i+j) % 2 == 1:
                        grid[i, j] = 0
                    else:
                        if grid[i, j] + 1/4.0*(gridcopy[i-1,j] + gridcopy[i+1,j] + gridcopy[i,j-1] + gridcopy[i,j+1]) > 0:
                            grid[i, j] = 1.0
            print(np.sum(np.sum(grid)))
            
        Radj = np.round(R*np.sqrt(NSteps)) / np.sqrt(2*k*t)
        lo = int(np.max([0, midpt-Radj]))
        hi = int(np.min([grid.shape[0], midpt+Radj]))

        return orighist, np.sum( grid[lo:hi, lo:hi] )

    

    if D == 3:
        grid = np.zeros((NSteps*2 + 6, NSteps*2 + 6, NSteps*2 + 6))
        midpt = grid.shape[0]//2
        grid[midpt, midpt, midpt] = 1.0

        for istep in range(1,NSteps+1):
            gridcopy = copy(grid)
            for i in range(midpt-istep,midpt+istep+1):
                for j in range(midpt-istep,midpt+istep+1):
                    for k in range(midpt-istep,midpt+istep+1):
                        if (istep+i+j+k) % 2 == 1:
                            grid[i, j, k] = 0
                        else:
                            if grid[i, j, k] + 1/6.0*(gridcopy[i-1,j,k] + gridcopy[i+1,j,k] + gridcopy[i,j-1,k] + gridcopy[i,j+1,k] + gridcopy[i,j,k-1] + gridcopy[i,j,k+1]) > 0:
                                grid[i, j, k] = 1.0
            print(np.sum(np.sum(np.sum(grid))))

            
        Radj = np.round(R*np.sqrt(NSteps)) / np.sqrt(2*k*t)
        lo = int(np.max([0, midpt-Radj]))
        hi = int(np.min([grid.shape[0], midpt+Radj]))

        return orighist, np.sum( grid[lo:hi, lo:hi, lo:hi] )


    if D == 4:
        grid = np.zeros((NSteps*2 + 6, NSteps*2 + 6, NSteps*2 + 6, NSteps*2 + 6))
        midpt = grid.shape[0]//2
        grid[midpt, midpt, midpt, midpt] = 1.0
        orighist = []

        for istep in range(1,NSteps+1):
            gridcopy = copy(grid)
            for i in range(midpt-istep,midpt+istep+1):
                for j in range(midpt-istep,midpt+istep+1):
                    for k in range(midpt-istep,midpt+istep+1):
                        for m in range(midpt-istep,midpt+istep+1):
                            if (istep+i+j+k+m) % 2 == 1:
                                grid[i, j, k, m] = 0
                            else:
                                if grid[i, j, k, m] + 1/8.0*(gridcopy[i-1,j,k,m] + gridcopy[i+1,j,k,m] + gridcopy[i,j-1,k,m] + gridcopy[i,j+1,k,m] + gridcopy[i,j,k-1,m] + gridcopy[i,j,k+1,m] + gridcopy[i,j,k,m-1] + gridcopy[i,j,k,m+1]) > 0:
                                    grid[i, j, k, m] = 1.0
            print(np.sum(np.sum(np.sum(np.sum(grid)))))
        #print("var %f sd  %f %f" % (sumvar, np.sqrt(sumvar), np.sum(np.sum(np.sum(np.sum(grid))))))
        return orighist, retval, np.sum( grid[lo:hi, lo:hi, lo:hi, lo:hi] )





# x = PDEGauss(2,1,1,1,200)
# from scipy.stats import linregress
# x = x[0]
# rng = np.log(1.0+np.arange(len(x)))
# rng = np.log(2.0+2*np.arange(len(x)))
# linregress(rng[25:],np.log(x[25:])) 
# linregress(rng[5:],np.log(x[5:])) 


# for D==1 x = PDEGauss(1,1,1,1,1000)
# linregress(rng[5:],np.log(x[5:]))
# LinregressResult(slope=-0.4980745578031345, intercept=-0.5836424158694404, rvalue=-0.999996081510863, pvalue=0.0, stderr=6.279809885717687e-05, intercept_stderr=0.00033524625441152444)

# D=2: x = PDEGauss(2,1,1,1,200)
# linregress(rng[25:],np.log(x[25:]))
# LinregressResult(slope=-0.9951416059364524, intercept=-1.1691063151479835, rvalue=-0.9999996102857909, pvalue=1.0438090750384467e-224, stderr=0.00010282807344726438, intercept_stderr=0.00042097590415575496)
# linregress(rng[5:],np.log(x[5:]))
# LinregressResult(slope=-0.990191588294568, intercept=-0.5030823198805514, rvalue=-0.999991729731274, pvalue=3.793222773666564e-224, stderr=0.00041759497174586636, intercept_stderr=0.0018897155786101634)


# for D=3 and 60 steps:
# linregress(rng[15:],np.log(x[15:]))
# LinregressResult(slope=-1.5533282636116321, intercept=-1.2375911395201697, rvalue=-0.9999952328299576, pvalue=5.035110036788766e-34, stderr=0.001330265845095701, intercept_stderr=0.004154705760723693)
# linregress(rng[5:],np.log(x[5:]))
# LinregressResult(slope=-1.5902953860036997, intercept=-0.02118388065452237, rvalue=-0.999905580429158, pvalue=2.4615982662054983e-44, stderr=0.004557121900941132, intercept_stderr=0.01603284128963471)


# for D=4 and 45 steps:
# linregress(rng[5:],np.log(x[5:]))
# LinregressResult(slope=-1.9589960151550359, intercept=-0.38391952292197473, rvalue=-0.9999940919860749, pvalue=2.2397031491090677e-38, stderr=0.0017387037884019172, intercept_stderr=0.00571264195428811)




# for the SimplexGRID:
# for D == 1, 1536 nodes (expansion=9) and nonzero nodes == 501 (or 62% of the available (i.e. half of the) nodes that aren't zer just due to parity)
# LinregressResult(slope=-0.49709033171431016, intercept=-0.5877235071398355, rvalue=-0.9999942921289714, pvalue=0.0, stderr=0.00010774226774210002, intercept_stderr=0.0005043209734374606)


# by using the "touch() " routine to determine dimensionality, the
# 2-d case is (suprisingly) obtained when the dimension is THREE
#python latticegas_nodes.py -p 1.0 --expand 7 --dim 3
# 4
# 10
# 21
# 29
# 47
#...
# 33023
# 33053
# 33125
# 33197"""
#python latticegas_nodes.py -p 1.0 --expand 7 --dim 3
#
# LinregressResult(slope=2.0486173118957374, intercept=1.5696524351598118, rvalue=0.9991357582229897, pvalue=3.4832669505870735e-137, stderr=0.008609185588291437, intercept_stderr=0.03230831696888312)
#
# curiously, this also gives a good dimensionalito of 2 when we use the HeatEq
# LinregressResult(slope=-1.04339647064362, intercept=-0.7319608184893913, rvalue=-0.9785299129199432, pvalue=9.500359862682954e-21, stderr=0.04153215390505293, intercept_stderr=0.10903517409533228)
# (remember, the dimensionality is DOUBLE the negative of the slope)

# we can get a 3-d behavior by both
# python latticegas_nodes.py -p 1.0 --expand 5 --dim 7
#  LinregressResult(slope=2.9481044675708596, intercept=2.1850549155592276, rvalue=0.997656348048403, pvalue=3.0182990882004047e-21, stderr=0.049039313962914166, intercept_stderr=0.10860580783300199)
# (the Heat Eq approach is not particularly indicative of 3d)
# LinregressResult(slope=-1.1775219195732618, intercept=-0.9459094882681156, rvalue=-0.9771779348467556, pvalue=4.625669742657659e-21, stderr=0.04753316113564067, intercept_stderr=0.12621254597737303)

# and
# python latticegas_nodes.py -p 1.0 --expand 4 --dim 8
# LinregressResult(slope=3.0247052866749646, intercept=2.288381149734917, rvalue=0.9978055136176972, pvalue=1.6048934353109006e-15, stderr=0.05794148100877665, intercept_stderr=0.11286218565115406)
# (not sure which is closer to 3.0, since they both level off as the graph
# is saturated
# the Heat Eq approach doesn't give a very linear exponent (it's more kinked) -- maybe due to saturation of the graph?
# LinregressResult(slope=-1.3096881078585298, intercept=-0.782732078002887, rvalue=-0.9754348119707248, pvalue=2.8302231300729043e-13, stderr=0.06971476670420912, intercept_stderr=0.15756521196097967)




def Krnl(t,r,k,D):
    multip = 1.0/ (4 * math.pi * k * t)**(D/2.0) 
    return multip * np.exp(-r*r/(4*k*t))

def S_nminus1(nminus1):
    if nminus1 == 0:
        return 2 # r**0
    if nminus1 == 1:
        return 2 * math.pi # * r**1
    if nminus1 == 2:
        return 4 * math.pi # * r ** 2

    n = nminus1n + 1.0
    return 2 * (math.pi ** (n/2)) / math.gamma(n/2)


def Integ1(R,k,t,delta=0.001):
    D = 1
    multip = 1.0/ (4 * math.pi * k * t)**(D/2.0) 
    integg = 0
    for i in np.arange(-R, R+delta, delta):
        integg += np.exp(-i * i / 4.0 / t)
    
    return multip * integg * delta


def Integ2(R,k,t,delta=0.01):
    D = 2
    multip = 1.0/ (4 * math.pi * k * t)**(D/2.0) 
    integg = 0
    for i in np.arange(-R, R+delta, delta):
        for j in np.arange(-R, R+delta, delta):
            if i*i + j*j <= R * R:
                integg += np.exp(-(i * i + j * j) / 4.0 / t)
    
    return multip * integg * delta * delta


def IntegKrnl(D, R, k, t):
    multip = 1.0/ (4 * math.pi * k * t)**(D/2) / S_nminus1(D-1) 
    a = 4 * k * t
    volelementfac = S_nminus1(D-1)
    if D == 1:
        return np.sqrt(a * math.pi) * 0.5 * math.erfc(R/np.sqrt(a))
    if D == 2:
        return -volelementfac * a/2.0 * (np.exp(-R*R/a) - 1.0)


    

# in D dimensions the surface element is S_nminus1(D-1) r**(D-1) * dr
# so u = r**(D-2) while dv = -r * Krn(t,r,k,D) / 2 / k / t  (which implies that v = Krnl)
# i.e.  -2kt * integ( u * dv ) = -2kt * ( u(R)v(R)  - 0   - ingteg( (D-2) * r**(D-3) )  )


def ReadParams():

  p = optparse.OptionParser()
    
  p.add_option("-d", "--dimension", "--dim", default=2,
                  action="store", type="int", dest='dimension',
                  help="")

  p.add_option("--picklefile", default="",
                  action="store", dest='picklefile', type='string',
                  help="pickle file, to which will be added the random seed ")

  p.add_option("--seed", default=84848484,
                  action="store", dest='seed', type='int',
                  help="rand number seed ")

  p.add_option("--prob", "--probability", "-p", default=0.001,
                  action="store", dest='prob', type='float',
                  help="each node, at each time steop, will be 'expanded' with this probability")

  p.add_option("--xprob", "--xprobability", default=0.000,
                  action="store", dest='xprob', type='float',
                  help="this only applies to a division of a cube surrounded by undivided cubes -- that should happen so rarely (as opposed to cubes situated next a divided cube -- that it prevents a Pareto distribution from forming)")

  p.add_option("--expand", "--expansion", "-t", default=10,
                  action="store", dest='expansionruns', type='int',
                  help="how many time steps do we expand, using the --prob arg on all the nodes available (understanding that the available nodes may increase with time)")


  p.add_option("--length", default=-1,
                  action="store", dest='length', type='int',
                  help="if > 0 the target number of nodes will be at or slightly greater than this parameter**Dim")



  p.add_option("--maxnodes", "--maxnode", default=-1000000,
                  action="store", dest='maxnodes', type='int',
                  help="the division loop is exited when maxnodes exceed")


  
                
  return p.parse_args()                


"""
Direction will be + or - N where N is the dimension (1-initialized, not zero-initialized), e.g.

+2 = (0, +1...)
-3 = (0, 0, -1,...)

so positive x direction is 1, negative is -1, positive y is +2, negative z is -3, etc...

"""



class cubehistoryvalues_t():
    # cubes (and their parents) will be keyed/indexed by its m-tuple nodes where m = 2**D
    # the
    def __init__(self, ):
        self.Children = []
        # a) the new nodes created when the cube was subdivided b) the decimal equilent center and c) 
        self.DecimalEquivalentCenter = -1
        self.DecimalEquivalentCubeLength = -1

# https://stackoverflow.com/questions/18500541/how-to-flatten-a-tuple-in-python
def flatten(data):
    """ this is a useful function, but it isn't used  here
    """
    if isinstance(data, tuple):
        if len(data) == 0:
            return ()
        else:
            return flatten(data[0]) + flatten(data[1:])
    else:
        return (data,)



class node_t():
    def __init__(self, id):
        self.Id = id
        self.Neighbors = [] # here, every element of the list is of type neighbor_t
        self.Amplitude = 0
        self.PrevAmplitude = 0
        self.bDefunct = False
        self.Parity = 0 # this is only used as a temporary scratch pad, though if nslices is even, the parity of a point during the expansion will remain constant
        self.Coords = None # not used for anything but post analysis and forensic analysis on nodes where the dimension parameters are not close to self.NDim
        self.LastDivide = 0 # not strictly needed, introduced in order to "even out" the divisions so that the resultant lattice does
        # not have large gaps

    def NeighborIds(self,):
        thislist = [nbr.Id for nbr in self.Neighbors]
        thislist.sort()
        return thislist

    

    def GetNeighbor(self, Id):
        for inbr in self.Neighbors:
            if inbr.Id == Id:
                return inbr
        return None
    
    def NeighborInfo(self, NeighborId = -1):
        if not(NeighborId in self.NeighborIds()) and NeighborId >= 0:
            return {}
        myret = {}           
        for inbr in self.Neighbors:
            if inbr.Id == NeighborId:
                myret["Id"] = inbr.Id
            else:
                myret[inbr.Id] = {}
        return myret


class neighbor_t():
    def __init__(self, targetId):
        self.Id = targetId # this is the node of the neighbor
        self.PrevNeighbors = (-1, -1) # if this was split from a larger edge, it keeps the (immediately preceding) larger ecdge
        # and if it is one of the edges created in a slice, the PrevNeighbors are the far nodes (i.e. extending to edge of cube) of the edge
        # in the initial starting grid, the PrevNeighbors are the same as the
        """
        to check if a cube is divisible, make sure that all the desired nodes are all traversible by paths whose
        component edges are all the same magnification (i.e. color); it's OK if the edges have been spliced to ONE higher
        magnficiation, but not more than one -- i.e. it can enconter one grain that is one-level finer; of course, this 
        all assumes the cube has not already been divided, in which case, its "inside" nodes have arleady been created


        Divide a cube:
        1) check that it has 2D neighbors
        2) pick D other neighbors (connected edges of the same color/magnification) that are not collinear and that have 2D neighbors each;
           check to see if they are part of a cube whose edges have the same color/magnification and whose nodes all have 2D neighbors.
        3) remember the slices can vary with  color/magnification; the important thing is that no cube gets left behind or put into a weird
        configuration that disallows it from being divided.


        """
        


        
        

class cubenodeset_t():
    def __init__(self, NDim): #, NMaxNeighbors=-1):
        """
        The default construction will be an NDim+1 dimensional cube (so as to have a genus equivalent to a ball in NDim+1 dimensions)

        Once it's created it can then be beefed up.
        As a start, we'll just divide any avalable cube into 2^NDim smaller cubes. Any smaller cube can be further subdivided, etc.
        
        Maybe, when every edgelength is set to 1, this will have the appropriately dimensioned random lattice (both in terms of how it grows with "radius", i.e.
        it should grow as R^NDim, and aso w.r.t to the implied dimensionality of any heat equation we generate there).

        But maybe that won't work and will instead produce a weird Hausdorff dimension, or some arrangement where the radius-growth dimension and the
        heat-dquation-at-origin dimension are the same.

        In that case, we may have to randomly reduce the edges in this structure, so as to produce a more sensible dimension (though I'm not sure how).
        Maybe you just tune things? But if that's the case, that implies a 3d dimensional world is a delicate and finely-balanced creation which goes
        against the overall dumb-things-are-good-enogh-to-make-a-universe ethos.
        """


        self.NDim = NDim # this is used in the initialization to make a simplex of size NDim+2
        self.NDimPlus1 = NDim + 1 # since we connect the NDim-dimensional structure into a sphere of one higher dimension, this is the effective dimension
        #self.NodeDict = {}
        self.NodeVec = []
        self.NNode = 0 

        self.NodeHistory = {} # all cubes -- the values will be of class

        self.TorLen = 8 # 32
        self.Cubes = [] # soecify the cubes with a 2-tuple consisting of an axis (self.Dim neighboring nodes) and the opposing node; we could also
                        # specify the fractality, but that can be accessed by way of the nieghbors
                        # THESE ARE ONLY THE CUBES THAT -- AS LONG AS ALL THE NODES OF THE CUBE HAVE 2*DIM NEIGHBORS, CAN BE DIVIDED
                        # Note we store this so as to avoid having to compute each of the subcubes by taking neighbors of neighbors...etc.
                        # We could also compute everything over and over again, which would involve a lot of redundant computing given that every
                        # subcube has 2**NDim edges,  and each time we pick one of those nodes as a place to subdivide, we'd have to redo the
                        # same test, but whatever course we choose, we must know something about
                        # ALL the cubes that intersect at a point before we divide any of them; otherwise, we may wind up subdividing a given
                        # edge more than once, and make it very difficult to similarly subdivide any of the cubes that share that edge, whereas
                        # if we only subdivde an edge once, and then wait for all the other edges in nearby cubes to be divided once, we can match
                        # up  the subdivisions without too much trouble
                        # The final alternative is to not let that objection get in our way, and instead 
                        # allow whatever kind of repeated subdivision that happens, and if that means that the adjacent
                        # cubes can't be divided (because there are too man slices resulting from that earlier subdivision), then we simply won't do
                        # any dividing there. That will STILL produce a useable cube, but given how sparse the subdivided cubes are, you will need
                        # to look at really enormous collections before you get something that's smooth and free of any lumpiness as you zoom out,
                        # and you won't be able to do it with a desktop, and maybe not with any computer around these days unless the budget is huge
                        # and the amount of memory is prodigious

                        # NOTE: Nodes with less than 2*NDim neighbors CANNOT be subdivded (at some point, when the adjoining cubes are divided and
                        # filled in, so to speak) the nodes will then have the full 2*NDim neighbors, and THEN you can divide them.                        

        self.NodeCubes = {} # this is a dictionary of the 2**Dim cubes that touch a given node
        self.MagnificationRingLen = 100 # 4 EVENTUALLY, ONCE IT CHECKS OUT, SET IT TO 4
        
        self.NDivisions = 0

        histshape = tuple([self.TorLen] * self.NDim)
        self.Hist = np.zeros(histshape).astype("int")
        self.NodesLastHistogramUpdate = 0

        self.BoundaryNodes = [] # Boundary nodes are any nodes who have 2*D neighbors but
        # who have neighbers whose own neighbor count is incomplete (i.e. less than 2D)
        # these will form a surface in the D dimensional space (though it may be fractal so that
        # its dimension may not be D-1) and any cube division will be overwhelmingly more
        # likely to happen on such boundary nodes as opposed to happening in an area where
        # the local topology is uniformly rectangular. (Obviously, since we start with  a pur
        # rectangular space, the universe growth has to start with such a freakish event, and
        # once the divisions propagate along the boundary throughout the cube and thereby  restore
        # uniformity at some higher magnification, we'll need another freakish  event to kick things
        # off). The hope is that if the boundary nodes form a unformly fractal subspace, the
        # average number of nodes in any region of space is uniform throughout the space sp that
        # Fick's laws can be applied.


        #self.AllDirections = list(np.arange(1,self.NDim+1)) +  list(np.arange(-1,-self.NDim-1,-1))

    def SumAmp(self,bPrint=True):
        myret = np.sum([inode.Amplitude for inode in self.NodeVec])
        if bPrint:
            print("Sumamp ", myret)
        return myret

    def NonZero(self,):
        nnon = np.sum([int(inode.Amplitude != 0) for inode in self.NodeVec])
        denom = 1.0/self.NNode
        pct = denom * nnon
        twicepct = pct * 2
        #print("Nonzero ", nnon, pct, twicepct)
        return (nnon, pct)
    def UpdatePrev(self,):
        for inode in self.NodeVec:
            inode.PrevAmplitude = inode.Amplitude

    def Touch(self, NSteps):
        print("Touch")
        myretarr = []
        x = self.SumAmp(False)
        myretarr.append( x/float(self.NNode) )
        print(x/float(self.NNode))
        for istep in range(NSteps):
            for inode in self.NodeVec:
                thisamp = 0
                nbrs = [jnode.Id for jnode in inode.Neighbors]
                for inbr in nbrs:                    
                    if self.NodeVec[inbr].PrevAmplitude  != 0:
                        thisamp += 1
                    
                inode.Amplitude = np.min([1,thisamp])
            self.UpdatePrev()
            x = self.SumAmp(False)
            myretarr.append( x/float(self.NNode))
            print(x/float(self.NNode))
            #print("Touch ", istep)
        return myretarr


    def MakeParity(self, Id=0):
        def IsAllDefined():
            for i in self.NodeVec:
                if i.Parity == 0:
                    return False
            return True

        self.NodeVec[0].Parity = 1
        StepChunk = 100
        while not(IsAllDefined()):
            for istep in range(StepChunk):
                for inode in self.NodeVec:
                    #import pdb; pdb.set_trace()
                    if inode.Parity != 0:
                        nbrs = [jnode.Id for jnode in inode.Neighbors]
                        for inbr in nbrs:
                            if self.NodeVec[inbr].Parity == inode.Parity:
                                print("failed at ", inode.Id)
                                #import pdb; pdb.set_trace()
                                return
                            self.NodeVec[inbr].Parity = -inode.Parity

    def ClearParity(self, ):
        for inode in self.NodeVec:
            inode.Parity = 0



    def HeatEq(self, NSteps, iprint=0, bPrint=True):
        print("Heat")
        myretarr = [ self.NodeVec[iprint].Amplitude ]
        print( self.NodeVec[iprint].Amplitude )

        self.UpdatePrev()
        self.SumAmp(False)
        for istep in range(NSteps):
            for inode in self.NodeVec:
                thisamp = 0
                nbrs = [jnode.Id for jnode in inode.Neighbors]
                for inbr in nbrs:
                    denomin = 1.0/len(self.NodeVec[inbr].Neighbors)
                    #if self.NodeVec[inbr].PrevAmplitude * denomin != 0:
                    #    import pdb; pdb.set_trace()
                    thisamp += self.NodeVec[inbr].PrevAmplitude * denomin
                inode.Amplitude = thisamp

            self.UpdatePrev()
            if True: #istep % 2 == 1:
                #print(self.NodeVec[self.NNode//2-1].Amplitude)
                if bPrint:
                    print( self.NodeVec[iprint].Amplitude )
                myretarr.append( self.NodeVec[iprint].Amplitude )
            
        return myretarr



    def MonitoredHeatEq(self, NSteps, iprint=0, bPrint=True):
        """
        this stops the simulation when the value at the starting node (the max value)
        is less than, say, 100 times the "heat death value" of 1.0/NNodes
        It also monitors the total number of nodes greater than the

        """
        myretheatarr = [ self.NodeVec[iprint].Amplitude ]
        myrettoucharr = [ self.NodeVec[iprint].Amplitude ]

        heatdeathval = 1.0 / self.NNode
        threshold = 2 * heatdeathval
        myretarr = [ self.NodeVec[iprint].Amplitude ]
        self.UpdatePrev()
        self.SumAmp(False)
        for istep in range(NSteps):
            for inode in self.NodeVec:
                thisamp = 0
                nbrs = [jnode.Id for jnode in inode.Neighbors]
                for inbr in nbrs:
                    denomin = 1.0/len(self.NodeVec[inbr[0]].Neighbors)
                    #if self.NodeVec[inbr].PrevAmplitude * denomin != 0:
                    #    import pdb; pdb.set_trace()
                    thisamp += self.NodeVec[inbr].PrevAmplitude * denomin
                inode.Amplitude = thisamp


            self.UpdatePrev()
            if True: #istep % 2 == 1:
                #print(self.NodeVec[self.NNode//2-1].Amplitude)
                if bPrint:
                    print(self.NodeVec[iprint].Amplitude)
                myretheatarr.append( self.NodeVec[iprint].Amplitude )



            sumnonzero = 0
            for inode in self.NodeVec:
                if inode.Amplitude > 0:
                    sumnonzero += 1
            #myretthresharr.append( self.NodeVec[iprint].Amplitude )
            myrettoucharr.append( sumnonzero )

        print("Amp")
        #import pdb; pdb.set_trace()
        for iamp in myretheatarr:
            print(iamp)
        
        print("Touch")
        for iamp in myrettoucharr:
            print(iamp * heatdeathval)                

            self.UpdatePrev()
            if True: #istep % 2 == 1:
                #print(self.NodeVec[self.NNode//2-1].Amplitude)
                if bPrint:
                    print(self.NodeVec[iprint].Amplitude)
                myretarr.append( self.NodeVec[iprint].Amplitude )
            
        return myretarr

                    

    def Wipe(self, ):
        for inode in self.NodeVec:
            inode.PrevAmplitude = 0
            inode.Amplitude = 0


    

    #DO NOT USE
    def AddEdge(self, x): # x is a 2-tuple whose last index is greater than the first
        # https://stackoverflow.com/questions/8024571/insert-an-item-into-sorted-list-in-python

        if not(x in self.Edges):
            bisect.insort(self.Edges, x)

                

    def DivideCube(self, CoordNode, NodeCoord, nslices):

        #print("""get the fractality of any of the 2**D subcubes intersecting at any point;
        #    cannot increase the fractality of any divided (higherfractalit) subcube until
        #    all the others have been brought up to the existing level
        #""")



        def GetMag(CoordNode, NodeCoord):
            maglist = []
            for icoord,inode in CoordNode.items():
                for jcoord,jnode in CoordNode.items():
                    if inode == jnode:
                        continue
                    
                    if jnode in self.NodeVec[inode].NeighborIds():
                        maglist.append(self.NodeVec[inode].GetNeighbor(jnode).RingMagnification)
                    else:
                        
                        for inbr in self.NodeVec[inode].Neighbors:
                            if jnode in inbr.PrevNeighbors:
                                maglist.append( (inbr.RingMagnification + 1) % self.MagnificationRingLen  )
                                break
            
            return (self.RingMin(maglist) + 1) % self.MagnificationRingLen                 

        def NIntervening(vec, basis):
            retsum = 0
            for icomp in vec:
                if icomp != 0 and icomp != basis-1:
                    retsum += 1
            return retsum 

        def modall(vec, basis):
            for icomp in vec:
                if icomp % basis != 0:
                    return False
            return True

        def GetCoords(thistuple, node0, nslices, axisedges, CoordNode, NodeCoord):
            dimlen = float(nslices+1)
            thisaxis = np.array(self.NodeVec[node0].Coords) 
            coord0 = copy(thisaxis) # this is the anchor vector
            for i in range(self.NDim):
                coord_thisdim = np.array( self.NodeVec[axisedges[i]].Coords )
                try:
                    if np.max(np.abs(coord_thisdim - coord0)) > 1:
                        #import pdb; pdb.set_trace() 
                        for iicomp, icomp in enumerate(coord0):
                            if np.abs(coord_thisdim[iicomp] - coord0[iicomp]) > 1:
                                if coord_thisdim[iicomp] > coord0[iicomp]:
                                    coord0[iicomp] += self.TorLen
                                elif coord_thisdim[iicomp] < coord0[iicomp]:
                                    coord_thisdim[iicomp] += self.TorLen
                        
                    thisaxis = thisaxis + (coord_thisdim - coord0) * thistuple[i]/dimlen
                except:
                    print("BAD")
                    import pdb; pdb.set_trace()

            for iicomp, icomp in enumerate(thisaxis):
                if icomp < 0:
                    thisaxis[iicomp] += self.TorLen
                thisaxis[iicomp] %= self.TorLen
            return tuple(thisaxis)

        def rectify(x, modulus):
            if x >= modulus:
                return x % modulus
            if x < 0:
                return x + modulus
            return x


        def GetJthNeighbor(thisnodestrip, j, NodeB, NodeA, bReverse=False):
            if bReverse:
                earlierj = j + 1
            else:
                earlierj = j - 1
            for inbr in self.NodeVec[thisnodestrip[earlierj]].Neighbors:
                if inbr.PrevNeighbors == (NodeB, NodeA) and not(inbr.Id in thisnodestrip.values()):
                    thisnodestrip[j] = inbr.Id
                    return inbr.Id
            #import pdb; pdb.set_trace()
            if np.min([NodeB, NodeA]) > self.TorLen ** self.NDim: # i.e. if this is NOT one of the starting nodes in the original cubeset
                print("ERROR: how come no neighbor of the earlier element of the strip has the right PrevNeighbors?")
                import pdb; pdb.set_trace()
            return None


        #import pdb; pdb.set_trace() 
        #if Id == 885 and CoordNode[(1,1)] == 854:
        #    import pdb; pdb.set_trace()

        #if Id in (671, 670, 639, 638) and self.NNode > 3900:
        #    import pdb; pdb.set_trace()
        Id = CoordNode[ tuple([0 for i in range(self.NDim)])]

        cubelen = nslices + 1
        nodelen = cubelen + 1

        newcoordnode = {}
        newnodecoord = {}

        bHistory = True
        bAssignNewFractalityToOldNodes = False


        #import pdb; pdb.set_trace()
        for key,val in CoordNode.items():
            newkey = tuple([cubelen * icomp for icomp in key])
            newcoordnode[newkey] = val
            newnodecoord[val] = newkey
        

        
        if bHistory:
            # find coordinate axes
            #axes = []

            
            
            # DEPRECATED
            if bAssignNewFractalityToOldNodes:
                newfractality = self.NodeVec[Id].Fractality + 1
                print("new fractality", newfractality)
                for coord,id in CoordNode.items():
                    self.NodeVec[id].Fractality = newfractality

            axisedges = [] # the first node of the edge is just (0,0,...) -- i.e. Id -- and is understood
            for i in range(self.NDim):
                thisaxis = np.zeros((self.NDim,)).astype("int")
                thisaxis[i] = 1
                #axes.append(tuple(thisaxis))
                axisedges.append(CoordNode[tuple(thisaxis)])

        alreadydone = []
        # first we will see if any of the required have already been created due to divisions
        # in adjacent D-cubes

        # first find the nodes situated on a (previous) edge connecting two nodes of the cube
        
        allcoords_ranked = []        
        for j in range(nodelen**self.NDim):
            thisvec = self.BaseN(j, nodelen, self.NDim)
            sumintervening = NIntervening(thisvec, nodelen)
            allcoords_ranked.append((sumintervening, tuple(thisvec)))
            #sumintervening = NIntervening(thisvec, nodelen)
            #if sumintervening == 0:
            if modall(thisvec, cubelen): 
                oldcoord = [icomp//cubelen for icomp in thisvec]
                thisnode = CoordNode[tuple(oldcoord)]
                thiscoord = tuple(thisvec)                
                #newcoordnode[thiscoord] = thisnode
                #newnodecoord[thisnode] = thiscoord
                if not(thisvec in alreadydone):
                    alreadydone.append(tuple(thisvec))

            """
            if False and sumintervening == 1:
                retsum = 0
                whichax = -1
                for iicomp, icomp in enumerate(vec):
                    if icomp != 0 and icomp != basis:
                        whichax = iicomp
                        break
                nodeA = copy(thisvec)
                nodeA[iicomp] = 0
                nodeB = copy(thisvec)
                nodeB[iicomp] = nodelen
            """



        #for j in range(nodelen**self.NDim):
        #    thisvec = self.BaseN(j, nodelen, self.NDim)
        #    sumintervening = NIntervening(thisvec, nodelen)
        #    allcoords_ranked.append((sumintervening, tuple(thisvec)))
        
        if False:
        # now determine the magnification of t
            thisringmag = -1

            nbrs = self.NodeVec[coordnode[alreadydone[0]]].Neighbors
            tmpmag = -1
            inbr = nbrs[0]
            #import pdb; pdb.set_trace()
            if inbr.Id in NodeCoord.keys(): 
                thisringmag = (inbr.RingMagnification + 1)  % self.MagnificationRingLen
            else:
                thisringmag = inbr.RingMagnification
        #thisringmag = GetMag(CoordNode, NodeCoord)

        

        allcoords_ranked.sort()
        #founddimensions = [0]


        # now, try to see if there exist some nodes from adjacent and previously divided cubes (that we should therefore avoid trying to recreate)
        for isum,icoord in allcoords_ranked:
            if icoord in alreadydone:
                continue
            if isum == 0:
                print("Why is this not already in alreadydone? Get rid of this check if it never fails.")
                import pdb; pdb.set_trace()
                continue # if isum is zero, thie coord is one of the nodes
            else:       
                                
                thisvec = list(icoord)
                whichaxklist = []
                for iicomp, icomp in enumerate(thisvec):
                    if icomp != 0 and icomp != cubelen:
                        whichaxklist.append(iicomp)
                
                for whichax in whichaxklist:

                    coordA = copy(thisvec)                
                    coordA[whichax] = 0       
                    nodeA = CoordNode[tuple(np.array(coordA)//cubelen)]
                    
                    coordB = copy(thisvec)                
                    coordB[whichax] = cubelen
                    nodeB = CoordNode[tuple(np.array(coordB)//cubelen)]
                    
                    #prevnode = copy(nodeA)
                    #if Id == 1 and (icoord == (1, 0) or icoord == (2, 0)):
                    #    import pdb; pdb.set_trace()
                    
                    
                    
                    donepairs = []
                    for inbr in set(list(self.NodeVec[nodeA].Neighbors) + list(self.NodeVec[nodeB].Neighbors)):                        
                        # check if there exist some nodes created from some adjacent cube's division that match up in terms of having the same PrevNeighbors 

                        if inbr.PrevNeighbors == (nodeA, nodeB):

                            thispair = [nodeA, nodeB]
                            thispair.sort()
                            thispair = tuple(thispair)
                            if thispair in donepairs:
                                continue

                            thisnodestrip = {}
                            thisnodestrip[0] = nodeA
                            bReverse = False
                            for j in range(1, cubelen): 
                                #import pdb; pdb.set_trace()
                                thiscoord = copy(thisvec) #[i * cubelen for i in thisvec] 
                                thiscoord[whichax] = j
                                thiscoord = tuple(thiscoord)                                
                                thisnode = GetJthNeighbor(thisnodestrip, j, nodeA, nodeB, bReverse)
                                if thisnode is None: # this happens if we're trying to look for neighbors in between the starting nodes of the cubeset
                                    continue

                                if thisnode in newnodecoord.keys():
                                    continue

                                newcoordnode[thiscoord] = thisnode
                                newnodecoord[thisnode] = thiscoord
                                if not(thiscoord in alreadydone):
                                    alreadydone.append(tuple(thiscoord))

                                #if j == 1:
                                #    for inbr in self.NodeVec[nodeA]:
                                #        if inbr.Id == nodeB:
                                #            thismagnification = inbr.Magnification
                                #           thisringmag = thismagnification % self.MagnificationRingLen

                            donepairs.append(tuple(thispair))        
                            if cubelen-1 in thisnodestrip.keys():
                                if thisnodestrip[1] != nodeB:
                                    thisnode = GetJthNeighbor(thisnodestrip, cubelen, nodeA, nodeB, bReverse)
                                    if thisnode != nodeB:
                                        print("PROBLEM: we have more slices to try and match up than the cubelen-1 we expected")
                            elif np.max(list(thisnodestrip.keys())) > 1:
                                    import pdb; pdb.set_trace()
                                    print("PROBLEM: we didn't find all the intervening nodes we expected to find -- something wrong.")

                        elif inbr.PrevNeighbors == (nodeB, nodeA):


                            thispair = [nodeA, nodeB]
                            thispair.sort()
                            thispair = tuple(thispair)
                            if thispair in donepairs:
                                continue


                            thisnodestrip = {}
                            thisnodestrip[cubelen] = nodeB
                            bReverse = True
                            for j in range(cubelen-1, 0, -1): #range(1, cubelen): 
                                #import pdb; pdb.set_trace()
                                thiscoord = copy(thisvec)
                                thiscoord[whichax] = j # nodelen - 1 - j
                                thiscoord = tuple(thiscoord)                                
                                thisnode = GetJthNeighbor(thisnodestrip, j, nodeB, nodeA, bReverse)
                                if thisnode is None: # this happens if we're trying to look for neighbors in between the starting nodes of the cubeset
                                    continue

                                if thisnode in newnodecoord.keys():
                                    continue      
                                                 
                                newcoordnode[thiscoord] = thisnode
                                newnodecoord[thisnode] = thiscoord

                                if NIntervening(thiscoord, nodelen) == self.NDim:
                                    print("ERROR: found an existing node IN THE INTERIOR OF CUBE -- it was already divided; why was it allowed to happen again? ")
                                    import pdb; pdb.set_trace()
                                if not(thiscoord in alreadydone):
                                    alreadydone.append(tuple(thiscoord))

                                #if j == 1:
                                #    for inbr in self.NodeVec[nodeB]:
                                #        if inbr.Id == nodeA:
                                #            thismagnification = inbr.Magnification
                                #            thisringmag = thismagnification % self.MagnificationRingLen

                            donepairs.append(tuple(thispair))
                            bReverse = True
                            if 1 in thisnodestrip.keys():
                                if thisnodestrip[cubelen-1] != nodeA:
                                    thisnode = GetJthNeighbor(thisnodestrip, 0, nodeB, nodeA, bReverse)                             
                                    if thisnode != nodeA:                                    
                                        print("PROBLEM: we have more slices to try and match up than the cubelen-1 we expected")
                                        import pdb; pdb.set_trace()
                            elif np.min(list(thisnodestrip.keys())) < cubelen - 1:                                    
                                    print("PROBLEM: we didn't find all the intervening nodes we expected to find -- something wrong.")
                                    import pdb; pdb.set_trace()

                        #if not(isum in founddimensions):
                        #    founddimensions.append(isum)


        



        # do some checks -- delete section later


        #BEGIN DELETE
        #BEGIN DELETE
        #BEGIN DELETE
     
        alreadydone_ranked = []
        for j in alreadydone:
            sumintervening = NIntervening(j, nodelen)
            alreadydone_ranked.append((sumintervening, j))

        alreadydone_ranked.sort()

        #if Id == 1:
        #    import pdb; pdb.set_trace()
        if False and alreadydone_ranked[-1][0] == 0:
            # this means that no already-created interior nodes (aside from the cube coords themselves) were found
            nodeAvec = np.zeros((self.NDim,)).astype("int")
            nodeA = coordnode[tuple(nodeAvec)]
            nodeBvec = copy(nodeAvec)
            nodeBvec[0] = 1
            nodeB = coordnode[tuple(nodeBvec)]
            nbrnodes = self.NodeVec[nodeA].NeighborIds() 
            if not(nodeB in nbrnodes):
                print("How is the once-removed vector not a neighbor of the zero node?")
                import pdb; pdb.set_trace()
            for inbr in self.NodeVec[nodeA].Neighbors:
                if inbr.Id == nodeB:
                    newthisringmag = (inbr.RingMagnification + 1) % self.MagnificationRingLen
                    if newthisringmag != thisringmag:
                        print("ERROR: why is the mag calc'd here not the same as thismagnification")
                        import pdb; pdb.set_trace()

        #END DELETE
        #END DELETE
        #END DELETE



        
        # now, create the new nodes that don't exist, give them the right neighbors, and you'll be done  
        for j in range(nodelen**self.NDim):
            thisvec = self.BaseN(j, nodelen, self.NDim)
            if tuple(thisvec) in alreadydone:
                continue
            # DELETE THE NEXT 2 LINES -- REDUNDANT since all nodes that pass that test are already in alreadydone
            if modall(thisvec, cubelen):                   
                continue # if they're all modulo cubelen, that means they were in the starting set
            


            thisnode = self.CreateNode()
            #if thisnode.Id == 7162:
            #    import pdb; pdb.set_trace()
            thistuple = tuple(thisvec)

            

            newcoordnode[thistuple] = thisnode.Id
            newnodecoord[thisnode.Id] = thistuple
            if bHistory:
                #import pdb; pdb.set_trace()
                if thisnode.Coords is None:
                    thisnode.Coords = GetCoords(thistuple, Id, nslices, axisedges, CoordNode, NodeCoord)
                    #xy = thisnode.Coords
                    #if self.NNode > 7169 and xy[0] >= 19 and xy[0] <= 20 and xy[1] >= 30 and xy[1] <= 31:
                        #print("herree")
                        #import pdb; pdb.set_trace()
                else:
                    print("why are we here -- we should only be looping through new nodes that don't have coords")



        if Id > 0 and (thisnode == 19901) and self.NNode >= 48772:
                print("1033")
                import pdb; pdb.set_trace()

        
        # now delete whichever old neighbors are part of the original (unsliced) cube        
        for key, val in NodeCoord.items():

            thisnode = key
            thisvec = list(val)
            thesenbrs = [jnode for jnode in self.NodeVec[thisnode].Neighbors]
            thesenbrnodes = [jnode.Id for jnode in thesenbrs]

            for nbr,inbr in zip(thesenbrs, thesenbrnodes):
                if inbr in NodeCoord.keys():
                    #self.NodeVec[inbr].PrevNeighbors[key] = inbr
                    self.NodeVec[thisnode].Neighbors.remove(nbr)


            #self.NodeVec[thisnode].Neighbors.sort()
            self.NodeVec[thisnode].Neighbors.sort(key=lambda x: x.Id, reverse=False)


        # Note the number of neighbors depends on the number of components in the vector that are 0 or nslices;
        # if that number zero, we're in the interior of the crube and the number of neighbors is 2*D, if it's 1 then
        # we're (in 3d) in the middle of  face then there are only 2*D-1, if it's 2 then it's 2*D-2, etc.

        #import pdb; pdb.set_trace()
        for j in range(nodelen**self.NDim):
            thisvec = self.BaseN(j, nodelen, self.NDim)
            thisnode = newcoordnode[tuple(thisvec)]



            for iiv, iv in enumerate(thisvec):
                nbrvecUp = list(copy(thisvec))
                nbrvecDn = list(copy(thisvec))

                nbrvecUp[iiv] = iv + 1 # if greater than cubelen, then Up nbr already exists from before, and we leave this alone
                if iv + 1 <= cubelen:
                    nbrnodeUp = newcoordnode[tuple(nbrvecUp)]
                    if not(nbrnodeUp in self.NodeVec[thisnode].NeighborIds()):
                        # it's already in there from some adjacent cube-division
                        nbrnodeUp_prevA = copy(nbrvecUp)       
                        nbrnodeUp_prevB = copy(nbrvecUp)                
                        nbrnodeUp_prevA[iiv] = 0
                        nbrnodeUp_prevB[iiv] = cubelen
                        nodeId_prevA = newcoordnode[tuple(nbrnodeUp_prevA)]
                        nodeId_prevB = newcoordnode[tuple(nbrnodeUp_prevB)]

                        thisnbrUp = neighbor_t(nbrnodeUp)
                        prevnbrs = [nodeId_prevA, nodeId_prevB]
                        prevnbrs.sort() # note that by doing this, we remove any "directionality" prev neighbors 
                        # so that self.NodeVec[x].GetNeighbor(y).PrevNeighbors is equal to self.NodeVec[y].GetNeighbor(x).PrevNeighbors
                        thisnbrUp.PrevNeighbors = tuple(prevnbrs)

                        self.NodeVec[thisnode].Neighbors.append(thisnbrUp)

                        if len(self.NodeVec[thisnode].Neighbors) > 2*self.NDim:
                            print("Neighbor count exceeded -- more than 2 * Dim? -- DELETE THIS CHECK")
                            import pdb; pdb.set_trace()


                nbrvecDn[iiv] = iv - 1 # if less than 0, then dn nbr already exists from before, and we leave this alone
                if iv - 1 >= 0:                    
                    nbrnodeDn = newcoordnode[tuple(nbrvecDn)]
                    #try:
                    self.NodeVec[thisnode].NeighborIds()
                    #except:
                    #    print("hmm")
                    #    import pdb; pdb.set_trace()
                    if not(nbrnodeDn in self.NodeVec[thisnode].NeighborIds()):
                        #continue # it's already in there from some adjacent cube-division
                        nbrnodeDn_prevA = copy(nbrvecDn)     
                        nbrnodeDn_prevB = copy(nbrvecDn)
                        nbrnodeDn_prevA[iiv] = 0
                        nbrnodeDn_prevB[iiv] = cubelen
                        nodeId_prevA = newcoordnode[tuple(nbrnodeDn_prevA)]
                        nodeId_prevB = newcoordnode[tuple(nbrnodeDn_prevB)]

                        thisnbrDn = neighbor_t(nbrnodeDn) # thismagnification) 
                        prevnbrs = [nodeId_prevA, nodeId_prevB]
                        prevnbrs.sort()
                        thisnbrDn.PrevNeighbors = tuple(prevnbrs)
                        self.NodeVec[thisnode].Neighbors.append(thisnbrDn)

                        if len(self.NodeVec[thisnode].Neighbors) > 2 * self.NDim:
                            print("Neighbor count exceeded -- more than 2 * Dim? -- DELETE THIS CHECK")
                            import pdb; pdb.set_trace()

            self.NodeVec[thisnode].Neighbors.sort(key=lambda x: x.Id, reverse=False)

        for coord,node in newcoordnode.items():
            if len(self.NodeVec[node].Neighbors) == 2:
                print("too low")
                import pdb; pdb.set_trace()
            
            self.NodeVec[node].LastDivide = self.NDivisions
        self.NDivisions += 1
        #import pdb; pdb.set_trace()
        #print("done with divide")
                



    def DistinctNeighbor(self, Id, nbrX):
        #for inbr in nbrX:
        #    if inbr != Id:
        #        return inbr
        #return None
        myret = list(nbrX)
        if Id in myret:
            myret.remove(Id)

        return tuple(myret)



        
    def bIsCorner(self, nodeId):

        if nodeId < self.TorLen**self.NDim: # if it's in the "original set", it  is automatically a corner, even if we haven't figured out how to rig its prev neighbors
            return True
        thisnode = self.NodeVec[nodeId]
        if len(thisnode.Neighbors) != 2*self.NDim:
            return False
        for inbr in thisnode.Neighbors:
            if not(nodeId in inbr.PrevNeighbors):
                return False
        return True
    
    def FarNeighbors(self, nodeId):
        if not(self.bIsCorner(nodeId)):
            return []
        #ringmags = []
        #for inbr in self.NodeVec[nodeId].Neighbors:
        #    ringmags.append(inbr.RingMagnification)
        
        minmag = self.RingMin(ringmags)
        bGoLower = minmag == self.RingMax(ringmags)
        myretval = []
        for inbr in self.NodeVec[nodeId].Neighbors:
            if (inbr.RingMagnification == minmag and not(bGoLower)) or np.min(inbr.PrevNeighbors) < 0:
                myretval.append(inbr.Id)
            else:
                xlist = list(inbr.PrevNeighbors)
                xlist.remove(nodeId)
                myretval.append(xlist[0])
        return myretval


    def FartherNeighbor(self, nodeId, nbrId):

        prevnbrs = self.NodeVec[nodeId].GetNeighbor(nbrId).PrevNeighbors

        bCheck = False
        if bCheck:
            derivedneighbors = self.GetPreviousNeighbors(nodeId, nbrId)
            if derivedneighbors != prevnbrs:
                print("mismatch")
                import pdb; pdb.set_trace()


        return self.DistinctNeighbor(nodeId, prevnbrs)


    def CreateNode(self,):
        newnode = node_t(self.NNode)        
        #self.NodeDict[self.NNode] = newnode
        self.NodeVec.append(newnode)
        #import pdb; pdb.set_trace()
        if len(self.NodeVec) != self.NNode + 1:
            print("something wrong -- go back to using NodeDict?")
            import pdb; pdb.set_trace()

        self.NNode += 1
        return newnode

    def Distance(self, vec1, vec2):
        direc = None
        absdist = np.sum(np.abs(np.array(vec1)-np.array(vec2)).astype("int"))            
        #import pdb; pdb.set_trace()
        if absdist == 1:
            ii = 1
            for i,j in zip(vec1,vec2):
                if i + 1 == j:
                    return absdist, ii
                elif i - 1 == j:
                    return absdist, -ii
                ii += 1
        return absdist, None

    def CreateNCube(self, NDim=-1):
        global GuideDict
        global NodeCoord
        global CoordNode 

        """ this is the NDim "universe", 
            i.e. in 1 dim, this will be 4 nodes on a circle
            in 2 dim, this will be 16 nodes, with each "3" node connected (in each dimension) to an "0" node with a wrap-around edge
               each subsqare gets divided into a tic-tac-toe hash configuration of 9 smaller squares
            in 3 dim, this will be a cube composed of 64 nodes



        """

        def SortTuple(zarr):
            zarr = list(zarr)
            zarr.sort()
            return tuple(zarr)

        def UpDn(Up1DnN1, id, idim, torlen):
            vecid = self.BaseN(id, torlen, self.NDim)
            increment = 1 if Up1DnN1==1 else -1
            reldig = (vecid[idim-1] + increment) % torlen
            vecid[idim-1] = reldig
            return vecid

               
        
        if NDim == -1:
            NDim = self.NDim
        self.NDimPlus1 = NDim + 1
        self.NDim = NDim

        nodecoord = {}
        coordnode = {}
        
        torlen = self.TorLen # must be greater than 4, or the collinearity tests will be off
        if torlen < 6:
            print("Lenght of the 'universe' must be > 4. Otherwise, the bIsCollinear will fail, since 3 collinear nodes will wrap around the torus and create a bogus face and the routine will fail")
            print("The resultsnt 'cubes' will then include configurations that wrap around the torus in some way but are also collinear.")
        for ix in range(torlen**self.NDim):
            thisnode = self.CreateNode()
            thisvec = tuple(self.BaseN(ix, torlen, self.NDim))
            nodecoord[ix] = thisvec
            coordnode[thisvec] = ix          

       
        for idim in range(self.NDim):
            for id in range(torlen**self.NDim):
                revvec = UpDn(+1, id, idim, torlen)  
                thisupnbr = neighbor_t (coordnode[tuple(revvec)])
                prevnbrlist = [id, thisupnbr.Id]
                prevnbrlist.sort()
                thisupnbr.PrevNeighbors = tuple(prevnbrlist) #SortTuple((id, coordnode[tuple(revvec)])) # (-1,-1) #
                self.NodeVec[id].Neighbors.append(thisupnbr)
                revvec = UpDn(-1, id, idim, torlen)
                thisdnnbr = neighbor_t (coordnode[tuple(revvec)])
                prevnbrlist = [id, thisdnnbr.Id]
                prevnbrlist.sort()
                thisdnnbr.PrevNeighbors = tuple(prevnbrlist) #(-1,-1) #SortTuple((id, coordnode[tuple(revvec)]))

                self.NodeVec[id].Neighbors.append(thisdnnbr)
        
        
        for id in range(torlen**self.NDim):
            #self.NodeVec[id].Neighbors.sort()
            self.NodeVec[id].Neighbors.sort(key=lambda x: x.Id, reverse=False)
            
        
        idlistsofar = [xnode.Id for xnode in self.NodeVec]
        for inode in idlistsofar:
            """ we could simply add in the vectors as follows, and
                isolate and identify one cube per node,
                but to be consistent we'll suss out the surrounding
                cubes at each node, even though it's grossly redundant
        for ix in range(torlen**self.NDim):
            thesecoord = np.array(nodecoord[ix]).astype("int")
            for ix in range(2**self.NDim):
                thisvec = tuple(np.vec(self.BaseN(ix, 2, self.NDim)).astype("int") + thesecoord)
                thesecubenodes.append( CoordNode[ix] )
            thesecubenodes.sort()
            self.Cubes.append(tuple(thesecubenodes))
            
            """

            thesecubenodes = [  ]
            allaxes = self.ReturnAllAxes(inode)

            #import pdb; pdb.set_trace()
            for theseaxes in allaxes:
                
                bTryPrevNeighborsToo = False
                (loccoordnode, locnodecoord) = self.FindCubeCoord(inode, theseaxes)
                """
                thiscube = list(locnodecoord.keys())

                thiscube.sort()
                thiscube = tuple(thiscube)
                if not(thiscube in self.Cubes):
                    self.Cubes.append(thiscube)
                """


        """
        for inode in idlistsofar:
            thislist = []
            for icube in self.Cubes:
                if inode in icube:
                    thislist.append(icube)
            thislist.sort()
            #A s noted, self.NodeCubes is only used to make it easier to run through all the cubes converging on a point
            # as opposed to doing all that searching every time you consider which  cubes can be subdividied
            self.NodeCubes[inode] = thislist
        """
        
        NodeCoord = nodecoord
        CoordNode = coordnode 


        
        bHistory = True
        if bHistory:
            for id, coord in NodeCoord.items():
                self.NodeVec[id].Coords = coord
        
    
    def OLDbIsItAFace(self, Id, nodeA, nodeB):
        """
        Note that this routine will fail if the toroidal length of universe is 4 or less.
        In that case, even if Id and the neighboring nodes A and B are collinear (so that
        they do NOT form a face), the routine will fail since there is a far neighbor of A and B
        that wraps around the toroid.

        The other solution is, as previously attempted, to endow each edge with a property that can
        be likened to a direction, which would allow an easier way of determining whether two 
        neighbors are oppositely oriented from a given node.


        """
        #print("Make sure you replace Magnification with RingMag")

        #if Id == 33 and (nodeA in [1034, 1053]) and (nodeB in [1034, 1053]):
        #    print("ok")
        #    import pdb; pdb.set_trace()
            #self.bIsItAFace(Id, nodeA, nodeB)
            

        nbrs = self.NodeVec[Id].NeighborIds()
        nbrA = self.NodeVec[Id].GetNeighbor(nodeA)
        nbrB = self.NodeVec[Id].GetNeighbor(nodeB)
        minmag = self.RingMin([nbrA.RingMagnification, nbrB.RingMagnification])
        #maxmag = self.RingMax([nbrA.Magnification, nbrB.Magnification])


        if nbrA.Magnification == nbrB.Magnification:
            intersected = list(set(self.NodeVec[nodeA].NeighborIds()).intersection(set(self.NodeVec[nodeB].NeighborIds())))
            if Id in intersected:
                intersected.remove(Id)

            if len(intersected) == 1:
                return True
            else:
                
                if np.max([nodeA, nodeB]) < self.TorLen**self.NDim:
                    return False 
                #import pdb; pdb.set_trace()
                xlist = list(nbrA.PrevNeighbors)
                if Id in xlist:
                    xlist.remove(Id)
                #else:
                #    print("why?")
                #    import pdb; pdb.set_trace()
                farnbrsA = self.FarNeighbors(xlist[0])
                xlist = list(nbrB.PrevNeighbors)
                if Id in xlist:
                    xlist.remove(Id)
                farnbrsB = self.FarNeighbors(xlist[0])                
                intersected = list(set(farnbrsA).intersection(set(farnbrsB)))
                if Id in intersected:
                    intersected.remove(Id)
                return len(intersected) == 1 
            

        nbrLo = nbrA
        nbrHi = nbrB
        
        if self.RingLessThan(nbrB.RingMagnification, nbrA.RingMagnification):
            nbrLo, nbrHi = nbrHi, nbrLo

        farnbrsLo = self.FarNeighbors(nbrLo.Id)
        farnbrsHi = self.FarNeighbors(self.DistinctNeighbor(Id, nbrHi.PrevNeighbors))
        intersected = list(set(farnbrsLo).intersection(set(farnbrsHi)))
        if Id in intersected:
            intersected.remove(Id)
        return len(intersected) == 1 



    def bIsItAFace(self, Id, nodeA, nodeB):
        """
        Note that this routine will fail if the toroidal length of universe is 4 or less.
        In that case, even if Id and the neighboring nodes A and B are collinear (so that
        they do NOT form a face), the routine will fail since there is a far neighbor of A and B
        that wraps around the toroid.

        The other solution is, as previously attempted, to endow each edge with a property that can
        be likened to a direction, which would allow an easier way of determining whether two 
        neighbors are oppositely oriented from a given node.

        This has been rewritten so that it does not involve any use of Magnification, and relies
        instead on testing whether a node or its neighbors is along the edge (in which case its
        neighbor count is not 2D)

        """

        nbrs = self.NodeVec[Id].NeighborIds()
        nbrA = self.NodeVec[Id].GetNeighbor(nodeA)
        nbrB = self.NodeVec[Id].GetNeighbor(nodeB)

        intersectedprev =  list(set(nbrA.PrevNeighbors).intersection(set(nbrB.PrevNeighbors)))
        if Id in intersectedprev:
            intersectedprev.remove(Id)
        if len(intersectedprev) > 0:
            return False

        if len(nbrs) != 2 * self.NDim:
            print("Why is bIsItAFace being called for a node that does not have 2D neighbors?")
            #intersected = list(set(nbrnbrA).intersection(set(nbrnbrB)))
            #if Id in intersected:
            #    intersected.remove(Id)
            #else:
            #    print("why is Id not in the intersection, if these are both neighbors of ID")   
            #return len(intersected) == 1   
        nbrA = nbrA.Id
        nbrB = nbrB.Id

        if len(self.NodeVec[nbrA].Neighbors) != 2 * self.NDim:
            nbrA = self.FartherNeighbor(Id, nbrA)[0]
        if len(self.NodeVec[nbrB].Neighbors) != 2 * self.NDim:
            nbrB = self.FartherNeighbor(Id, nbrB)[0]

        nbrnbrA = [] 
        for inbr in self.NodeVec[nbrA].NeighborIds():
            if len(self.NodeVec[inbr].Neighbors) != 2 * self.NDim:
                nbrnbrA.append( self.FartherNeighbor(nbrA, inbr)[0] )
            else:
                nbrnbrA.append(inbr)
        nbrnbrB = [] 
        for inbr in self.NodeVec[nbrB].NeighborIds():
            if len(self.NodeVec[inbr].Neighbors) != 2 * self.NDim:
                nbrnbrB.append( self.FartherNeighbor(nbrB, inbr)[0] )
            else:
                nbrnbrB.append(inbr)
        intersected = list(set(nbrnbrA).intersection(set(nbrnbrB)))


        if Id in intersected:
            intersected.remove(Id) 
        return len(intersected) == 1



    def bNonCollinear(self, Id, nbrids):
        for iinbr, inbr in enumerate(nbrids):
            for jjnbr in range(iinbr+1, len(nbrids)):
                jnbr = nbrids[jjnbr]
                #if Id==1 and (inbr,jnbr)==(0,9):
                #    print("why")
                #    import pdb; pdb.set_trace()
                if not(self.bIsItAFace(Id, inbr, jnbr)):
                    return False
     
        return True

    def bNonCollinearOLD(self, Id, nbrids):
        blenA = len(self.NodeVec[Id].NeighborIds()) == 2 * self.NDim
        blen0 = len(self.NodeVec[nbrids[0]].NeighborIds()) == 2 * self.NDim
        blen1 = len(self.NodeVec[nbrids[1]].NeighborIds()) == 2 * self.NDim
        if blenA and blen0 and blen1:
            if nbrids[1] == self.NextLinearStep(nbrids[0], Id):
                return False
        if nbrids[1] == self.NextLinearStep(nbrids[0], Id):
            return False
        return True



    def RingMax(self, xarr):
        setxarr = set(xarr)
        if len(setxarr) > 2:
            print("RingMax() is used with the understanding that set(xarr) will have at most 2 values differing by at most 1--something is wrong")
            import pdb; pdb.set_trace()
            
        thismin = np.min(xarr)
        thismax = np.max(xarr)
        if thismax == self.MagnificationRingLen-1:
            if thismin != 0:
                print("if we cross the ring, then the minimum val must be 0 and the max must be FractalityRingLen-1")  
            return 0
        return thismax

    def RingMin(self, xarr):
        setxarr = set(xarr)
        if len(setxarr) > 2:
            print("RingMin() is used with the understanding that set(xarr) will have at most 2 values differing by at most 1--something is wrong")
            import pdb; pdb.set_trace()
            
        thismin = np.min(xarr)
        thismax = np.max(xarr)
        if thismax == self.MagnificationRingLen-1:
            if thismax != 0:
                print("if we cross the ring, then the minimum val must be 0 and the max must be FractalityRingLen-1")  
            return self.MagnificationRingLen-1
        return thismin
    
    def RingLessThan(self, ringmagA, ringmagB):
        return (ringmagA != ringmagB) and ringmagA == self.RingMin([ringmagA, ringmagB])
    def RingGreaterThan(self, ringmagA, ringmagB):
        return (ringmagA != ringmagB) and ringmagA == self.RingMax([ringmagA, ringmagB])
    def RingLessThanOrEqual(self, ringmagA, ringmagB):
        return ringmagA == self.RingMin([ringmagA, ringmagB])
    def RingGreaterThanOrEqual(self, ringmagA, ringmagB):
        return ringmagA == self.RingMax([ringmagA, ringmagB])

    def RingDecrement(self, ringmag):
        ringmag = ringmag - 1
        if ringmag < 0:
            ringmag += self.MagnificationRingLen
        return ringmag


    def RingIncrement(self, ringmag):
        return (ringmag + 1) % self.MagnificationRingLen

    def bAbuttingPeriphery(self, Id, nbrset):
        """
        Here's a bad situation: Id is within a divided cube, but the proposed axes nodes are on the edge.
        In this case, there is no sensible "far neighbor" to use in place of the axes nodes, because those
        axes nodes ARE the far neighbor

        """

        for inbr in nbrset:
            #prevnbrs = self.GetPreviousNeighbors(Id, inbr)
            prevnbrs = self.NodeVec[Id].GetNeighbor(inbr).PrevNeighbors
            try:
                if len(self.NodeVec[inbr].Neighbors) != 2 * self.NDim and not(Id in prevnbrs):
                    return True
            except:
                import pdb; pdb.set_trace()
        return False



    def ReturnAllAxes(self, Id):
        nbrs = self.NodeVec[Id].NeighborIds() 
        lennbrs = len(nbrs)
        if self.NDim == 1:
            return [ nbrs ]

        myret = []

        indexcombos = AllCombosNoReplacement(len(nbrs), self.NDim)
        combos = []
        for icomb in indexcombos:
            combos.append([nbrs[i] for i in icomb])
        


        #import pdb; pdb.set_trace()
        
        # iterate over all combos -- redundant, but easy to code


        
        dellist = []
        for iicombo, icombo in enumerate(combos):
            #if self.NNode >= 13600 and Id == 11759 and (tuple(icombo) == (11758,12722) or tuple(icombo) == (12722, 11758) ):
            #    print("boooab")
            #    import pdb; pdb.set_trace()
            #if Id == 9 and tuple(icombo) == (7028, 7030) and self.NNode >= 230244:
            #    import pdb; pdb.set_trace()
            #if Id == 0 and (tuple(icombo) == (4317, 12496) or tuple(icombo) == (12496, 4317)) and self.NNode == 106798:
            #    import pdb; pdb.set_trace()            #bnoncollinear = self.bNonCollinear(Id, icombo, True)
            #if Id == 0:
            #    import pdb; pdb.set_trace()
            if self.bAbuttingPeriphery(Id, icombo):
                 dellist.append(iicombo)
            if not(self.bNonCollinear(Id, icombo)):
                #if not(self.bNonCollinear(Id, icombo, True)):
                if not(iicombo in dellist):
                    dellist.append(iicombo)
                    
        
        dellist.sort()
        dellist.reverse()
        for idel in dellist:
            del combos[idel]

        # Note that if (x,y) is in combos, so is (y, x); they're different axes

        return combos


    # deprecate in favor of nextlinearstep
    def GoLinearly(self, FromId, ToId):
        nbrs = self.NodeVec[ToId].NeighborIds()
        for inbr in nbrs:
            if inbr == FromId:
                continue
            if not(self.bIsItAFace(ToId, FromId, inbr)):
                return inbr
        return None
                                 

    def NextLinearStep(self, FromId, ToId):
        FromNbrs = self.NodeVec[FromId].NeighborIds()
        ToNbrs = self.NodeVec[ToId].NeighborIds()

        for inbr in ToNbrs:
            inbrnbrs =  self.NodeVec[inbr].NeighborIds()
            intersected = list(set(FromNbrs).intersection(set(inbrnbrs)))
            if len(intersected) == 1:
                return inbr 
        return None
    

        


    def GetPreviousNeighbors(self, Id, nbrId):
        """
        Notes: You only ever need to know PrevNeighbors of a node with less than 2D neighbors
        The nodes of interest for a point with less than 2D neighbors, when it comes to locating
        prev neighbors, consists of other nodes with less than 2D neighbors. (Note that for D > 2,
        there may be multiple neighbors which lead to other nodes of less than 2D neigbors, in which
        case test for collinearity within that subspace; once you find a different number of neighbors 
        note the difference can be positive if you go from a cube edge to a corner, or else negative if
        you go from a face to a cube ednge)

        """

        def Loop(last, nextId):
            
            while len(self.NodeVec[nextId].Neighbors) != 2 * self.NDim and (len(self.NodeVec[nextId].Neighbors) >= len(self.NodeVec[last].Neighbors)):  
                #nextId = self.GoLinearly(prevlast, last)
                nbrsprevlast = copy(self.NodeVec[last].NeighborIds())
                for inbr in nbrsprevlast:
                    intersected = list(set(self.NodeVec[prevlast].NeighborIds()).intersection(set(self.NodeVec[last].NeighborIds())))
                    if prevlast in intersected:
                        intersected.remove(prevlast) 
                    else:
                        print("Why is prevlast not in the intersected set? Something is wrong.")
                        import pdb; pdb.set_trace()
                    if len(intersected) > 0:
                        continue
                    
                    last = copy(nextId)
                    nextId = copy(inbr)
            return(nextId)

        if len(self.NodeVec[Id].Neighbors) == 2 * self.NDim and len(self.NodeVec[nbrId].Neighbrs) == 2:
            print("Why do you need to get prev neighbors of an edge that has 2D neighbors at both nodes?")
            return None, None

        nextId = copy(nbrId)
        last = copy(Id)
        prevlast = copy(Id)        
        
        EndPtA = Loop(last, nextId) 

        nextId = copy(Id)
        last = copy(nbrId)
        prevlast = copy(nbrId)  

        EndPtB = Loop(last, nextId) 

        myret = [EndPtA, EndPtB]
        myret.sort()
        return myret


    def CloseNeighbor(self, Id, FarId):
        for thisnbr in self.NodeVec[Id].Neighbors:
            if FarId in thisnbr.PrevNeighbors:
                return thisnbr.Id
        return None
    
    def bCheckForInteriorPoints(self, Id, nbrsubset):
        bOnlyCheckForInterior = True
        coordnode,nodecoord = self.FindCubeCoord(Id, nbrsubset, bOnlyCheckForInterior)
        if coordnode is None:
            return False
        else:
            #print("Found interior node for ", Id, nbrsubset)
            #import pdb; pdb.set_trace()
            return True

    
    def bIsNonUniform(self, coordnode):
        """
        If the edges connecting the nodes have divisions (i.e. the neighbors are not the same as the corresponding nodes on the coordnode cube)
        then we can subdivide the cube. Otherwise, no
        """

        def OppEdges(icoord, coordnode):
            myretval = []
            for idim in range(self.NDim):
                thisvec = list(icoord)
                if thisvec[idim] == 0:
                    thisvec[idim] = 1
                elif thisvec[idim] == 1:
                    thisvec[idim] = 0
                else:
                    print("PROBLEM -- the icoord vector", icoord, "should only be composed of 0's and 1's, and instead there is a ", thisvec[idim])
                    import pdb; pdb.set_trace()
                thisvec = tuple(thisvec)
                if not thisvec in coordnode.keys():
                    print("PROBLEM -- the icoord vector", icoord, "should only be in coordnode keys ")
                    import pdb; pdb.set_trace()
                myretval.append(coordnode[thisvec])
            return myretval

        for icoord, inode in coordnode.items():
            inbrs = self.NodeVec[inode].NeighborIds()
            for jnode in OppEdges(icoord, coordnode):
                if not(jnode in inbrs):
                    return True
        return False

    def nbrsfrombin(self, nbrsubset, binvec, coordnode, nodecoord):
            """
            Used to get the neighbor nodes of the corner node
            by selecting the nodes corresponding to the 1's in binvec
            """

            def componentvec(iiv, coordnode, nodecoord):
                avec = np.zeros((self.NDim,)).astype("int")
                avec[iiv] = 1
                return coordnode[tuple(avec)]
                
            myret = []
            for iibv, bv in enumerate(binvec):
                if bv == 1:
                    #myret.append(nbrsubset[iibv])
                    myret.append(componentvec(iibv, coordnode, nodecoord))
            return myret


    # this is just a customized copy of FindCubeCoord()
    def FindOppositeNode(self, Id, nbrsubset):

        if len(self.NodeVec[Id].Neighbors) != 2 * self.NDim:
            #print("ERROR --  the corners of the cube do not have 2*D neighbors?")
            return None, None         
    
        Idcoord = tuple([0] * self.NDim)
        coordnode = {Idcoord:Id}
        nodecoord = {Id:Idcoord}

        for iic, inbrId in enumerate(nbrsubset):
            thisvec = np.zeros((self.NDim,)).astype("int")
            thisvec[iic] = 1
            thisvec = tuple(thisvec)
            coordnode[thisvec] = inbrId
            nodecoord[inbrId] = thisvec

        sortbyones = []        
        for j in range(2**self.NDim):
            thisnode = copy(Id) 
            binstr = bin(j)[2:]
            while len(binstr) < self.NDim:
                binstr = '0' + binstr
            lenbinstr = len(binstr)
            binvec = [int(chr) for chr in binstr]
            sortbyones.append((np.sum(binvec), tuple(binvec)))   
        
        sortbyones.sort()
        for sum1, binvec in sortbyones[1:]: # we already installed the first one
            if sum1 == 0:
                continue
            if sum1 == 1:
                if binvec in coordnode.keys():
                    continue
                thisnode = self.nbrsfrombin(nbrsubset, binvec)[0]
                coordnode[binvec] = thisnode
                nodecoord[thisnode] = binvec
                continue
            if sum1 >= 2:
                thesenodes = tuple(self.nbrsfrombin(nbrsubset, binvec, coordnode, nodecoord)) # note there are sum1 components of thesenodes
                nbrs = []
                #import pdb; pdb.set_trace()
                for itup in thesenodes:
                    subsetnbrs = []
                    thesenbrs = self.NodeVec[itup].NeighborIds()
                    for jnbrId in thesenbrs:
                        if len(self.NodeVec[jnbrId].Neighbors) == 2 * self.NDim:
                            subsetnbrs.append( jnbrId )
                        else:
                            farneighbor = self.FartherNeighbor(itup, jnbrId)[0]
                            #if len(self.NodeVec[farneighbor].Neighbors) != 2 * self.NDim and farneighbor >= self.TorLen ** self.NDim:
                            #    print("trouble")
                            #    #import pdb; pdb.set_trace()
                            subsetnbrs.append(farneighbor)
                    nbrs.append(subsetnbrs)
                    
                    

                intersected = list(set(nbrs[0]).intersection(*nbrs[1:]))

                dellist = []
                for ixtec in intersected:
                    if ixtec in  nodecoord.keys():
                        dellist.append(ixtec)
                for idel in dellist:                    
                    intersected.remove(idel) 
                    if (idel, ) != (Id,) and len(idel) != 0:
                            print("dellist should only contain 1 element! And it should be Id -- what's going on?")

                if len(intersected) != 1:
                    #import pdb; pdb.set_trace()
                    return None, None # rerun with minmag
                #import pdb; pdb.set_trace()
                if len(self.NodeVec[thisnode].Neighbors) != 2 * self.NDim:
                    return  None, None # rerun with minmag
                thisnode = intersected[0]
                coordnode[binvec] = thisnode
                nodecoord[thisnode] = binvec

                #if np.sum(binvec) >= self.NDim and bOnlyCheckForInterior:
                    #import pdb; pdb.set_trace()
                    #print("Found an interior point in the desired cube -- this cube has already been divided")
                #    return None, None
                continue

            nextnodes = []            
            for jsum1, jbinvec in sortbyones[1:]:
                if jsum1 == sum1-1:
                    thisnode = coordnode[jbinvec]
                    for jnbr in self.NodeVec[thisnode].NeighborIds():
                        if len(self.NodeVec[jnbr].Neighbors) == 2 * self.Ndim:
                            nextnodex.append(jnbr)
                        else:
                            nextnodex.append(self.FartherNeighbor(thisnode, jnbr)[0])
            
            #import pdb; pdb.set_trace()


            intersected = list(set(nextnodes[0]).intersection(*nextnodes))
            if len(intersected) != 1:
                return None, None

            thisnode = intersected[0]
            coordnode[binvec] = thisnode
            nodecoord[thisnode] = binvec
            if (len(thisnode).Neighbors()) != 2 * self.NDim:
                return None, None


        



        





        #    bIsUniform = True
        #    for inode in nodecoord.keys():
        #        bIsUniform = bIsUniform and len(self.NodeVec[inode].Neighbors) == 2 * self.NDim
        #        if not(bIsUniform):
        #            break
        #    if bIsUniform:
        #        #import pdb; pdb.set_trace()
        #        return None, None

        
        return coordnode, nodecoord


   
    def FindCubeCoord(self, Id, nbrsubset):
        """

        generate the hypercube coordinates,
        and sort them according to how many 1's they have,
        since the calc of those with a given number of ones (say M)
        are the intersection of the neighbors of the nodes
        generated by nodes correesponding to a number of ones that is M-1


        
        """
        def nbrsfrombin(nbrsubset, binvec, coordnode, nodecoord):
            """
            Used to get the neighbor nodes of the corner node
            by selecting the nodes corresponding to the 1's in binvec
            """

            def componentvec(iiv, coordnode, nodecoord):
                avec = np.zeros((self.NDim,)).astype("int")
                avec[iiv] = 1
                return coordnode[tuple(avec)]
                
            myret = []
            for iibv, bv in enumerate(binvec):
                if bv == 1:
                    #myret.append(nbrsubset[iibv])
                    myret.append(componentvec(iibv, coordnode, nodecoord))
            return myret


        #nbrids = self.NodeVec[Id].NeighborIds()

        #if self.NNode >= 13630 and Id == 11759 and (tuple(nbrsubset) == (11758,12722) or tuple(nbrsubset) == (12722, 11758) ):
        #    print("bab")
        #    import pdb; pdb.set_trace()


        #if self.NNode == 1036:
        #    import pdb; pdb.set_trace()

        # Check that the corners of the cube have 


        if len(self.NodeVec[Id].Neighbors) != 2 * self.NDim:
            #print("ERROR --  the corners of the cube do not have 2*D neighbors?")
            return None, None         
    




        Idcoord = tuple([0] * self.NDim)
        coordnode = {Idcoord:Id}
        nodecoord = {Id:Idcoord}

        bUsingFarNeighbors = False
        for iic, inbrId in enumerate(nbrsubset):
            thisvec = np.zeros((self.NDim,)).astype("int")

            relevantnode = inbrId
            if inbrId in self.NodeVec[Id].NeighborIds():
                #fails when Id is 885 and icombo is [2021, 2025]
                if len(self.NodeVec[inbrId].Neighbors) != 2 * self.NDim: # and not(bOnlyCheckForInterior):
                    relevantnode = self.FartherNeighbor(Id, inbrId)[0]
                    bUsingFarNeighbors = True
            
            if len(self.NodeVec[relevantnode].Neighbors) != 2 * self.NDim: # and not(bOnlyCheckForInterior):
                return None, None

            thisvec[iic] = 1
            thisvec = tuple(thisvec)
            coordnode[thisvec] = relevantnode
            nodecoord[relevantnode] = thisvec

        if bUsingFarNeighbors:         
            innercoordnode, innernodecoord = self.FindOppositeNode(Id, nbrsubset)
            if not(innercoordnode is None):
                # there are interior points in the cube already -- must abort
                return None, None

            
        #if self.NNode >= 7159 and Id == 4001: #and (tuple(nbrsubset) == (3998,4000) or tuple(nbrsubset) == (4000, 3998) ):
        #    print("bab")
        #    import pdb; pdb.set_trace()



        # first enter the 2**NDim coordinates of the "big cube" whose edges are of length cubelen
        sortbyones = []
        
        for j in range(2**self.NDim):
            #print(j)
            thisnode = copy(Id) # always start at the anchor node
            binstr = bin(j)[2:]
            while len(binstr) < self.NDim:
                binstr = '0' + binstr
            lenbinstr = len(binstr)
            binvec = [int(chr) for chr in binstr]
            sortbyones.append((np.sum(binvec), tuple(binvec)))   
        
        sortbyones.sort()
        for sum1, binvec in sortbyones[1:]: # we already installed the first one
            if sum1 == 0:
                continue
            if sum1 == 1:
                if binvec in coordnode.keys():
                    continue
                thisnode = self.nbrsfrombin(nbrsubset, binvec)[0]
                coordnode[binvec] = thisnode
                nodecoord[thisnode] = binvec
                continue
            if sum1 >= 2:


                thesenodes = tuple(self.nbrsfrombin(nbrsubset, binvec, coordnode, nodecoord)) # note there are sum1 components of thesenodes

                nbrs = []

                
                for itup in thesenodes:
                    subsetnbrs = []
                    thesenbrs = self.NodeVec[itup].NeighborIds()
                    for jnbrId in thesenbrs:
                        if len(self.NodeVec[jnbrId].Neighbors) == 2 * self.NDim:
                            subsetnbrs.append( jnbrId )
                        else:
                            farneighbor = self.FartherNeighbor(itup, jnbrId)[0]
                            #if len(self.NodeVec[farneighbor].Neighbors) != 2 * self.NDim and farneighbor >= self.TorLen ** self.NDim:
                            #    print("trouble")
                            #    #import pdb; pdb.set_trace()
                            subsetnbrs.append(farneighbor)
                    nbrs.append(subsetnbrs)
                    
                    
                
                intersected = list(set(nbrs[0]).intersection(*nbrs[1:]))

                dellist = []
                for ixtec in intersected:
                    if ixtec in  nodecoord.keys():
                        dellist.append(ixtec)
                
                for idel in dellist:                    
                    intersected.remove(idel) 
                    if (idel, ) != (Id,) and len(idel) != 0:
                            print("dellist should only contain 1 element! And it should be Id -- what's going on?")
                
                if len(intersected) != 1:
                    import pdb; pdb.set_trace()
                    return None, None # rerun with minmag
                #import pdb; pdb.set_trace()
                if len(self.NodeVec[thisnode].Neighbors) != 2 * self.NDim:
                    return  None, None # rerun with minmag
                thisnode = intersected[0]
                coordnode[binvec] = thisnode
                nodecoord[thisnode] = binvec

                #if np.sum(binvec) >= self.NDim and bOnlyCheckForInterior:
                    #import pdb; pdb.set_trace()
                    #print("Found an interior point in the desired cube -- this cube has already been divided")
                #    return None, None
                continue

            nextnodes = []            
            for jsum1, jbinvec in sortbyones[1:]:
                if jsum1 == sum1-1:
                    thisnode = coordnode[jbinvec]
                    for jnbr in self.NodeVec[thisnode].NeighborIds():
                        if len(self.NodeVec[jnbr].Neighbors) == 2 * self.Ndim:
                            nextnodex.append(jnbr)
                        else:
                            nextnodex.append(self.FartherNeighbor(thisnode, jnbr)[0])
            
            import pdb; pdb.set_trace()


            intersected = list(set(nextnodes[0]).intersection(*nextnodes))
            if len(intersected) != 1:
                return None, None

            thisnode = intersected[0]
            coordnode[binvec] = thisnode
            nodecoord[thisnode] = binvec
            if (len(thisnode).Neighbors()) != 2 * self.NDim:
                return None, None


        


        # begin check
        bCheck = False
        if bCheck:
            for coord,node in coordnode.items():
                if len(self.NodeVec[node].Neighbors) != 2 * self.NDim:
                    print("wrong1 -- failed coordinate creation")
                    #import pdb; pdb.set_trace()
            for coord,node in coordnode.items():
                nbrs = self.NodeVec[node].NeighborIds()
                for inbr in nbrs:
                    if inbr in nodecoord.keys():
                        continue
                    farnbr = self.FartherNeighbor(node, inbr)[0]
                    if not(farnbr in nodecoord.keys()):
                        if node >= self.TorLen ** self.NDim:
                            print("wrong2")
                            import pdb; pdb.set_trace()
        # end check

        # now, make sure there is some nonuniformity in the local neighborhood
        #if (Id == 34 or Id == 67 or Id == 69) and self.NNode > 1024: #               and tuple(nbrsubset) == (2325, 2327):
        #    print("findcoord")
            #import pdb; pdb.set_trace()  


       
        #if self.NNode > self.TorLen ** self.NDim + 1:
        #    if not(self.bIsNonUniform(coordnode)):
        #        return None, None

            #print("thisisgood")
            #import pdb; pdb.set_trace()


        if self.NNode > self.TorLen ** self.NDim + 1:
            if not(self.bIsNonUniform(coordnode)):
                if rn.random() > FreshProb:
                    return None, None
                else:
                    #print("Passed the uniformity test")
                    pass
            else:
                #print("mixed")
                pass

            #print("thisisgood", coordnode)
            #import pdb; pdb.set_trace()

        #    bIsUniform = True
        #    for inode in nodecoord.keys():
        #        bIsUniform = bIsUniform and len(self.NodeVec[inode].Neighbors) == 2 * self.NDim
        #        if not(bIsUniform):
        #            break
        #    if bIsUniform:
        #        #import pdb; pdb.set_trace()
        #        return None, None

        
        return coordnode, nodecoord

    def FindLargestUndividedCubes(self, Id, allaxes):
        # for every axis analyzed, the nodes will be sorted and added to
        # the alreadyanalyzed list, and checked against that
        alreadyanalyzed = []

        myret = {} # return the magnification of all indices -- only the axes having
        # the largest magnification (only 2 magnifications will be present in any possible
        # set of axes)

        for iax in allaxes:
            cpyiax = list(copy(iax))
            cpyiax.sort()
            cpyiax = tuple(iax)
            if cpyiax in alreadyanalyzed:
                myret[iax] = myret[cpyiax]
                continue



        if len(set(myret.values())) > 2:
            print("ONLY TWO MAG VALUES SHOULD BE PRESENT IN ANY allaxes -- WHAT WENT WRONG?")
            import pdb; pdb.set_trace()



    def CheckThatAllHave2DNeighbors(self, CoordNode, NodeCoord):
        for coord, node in CoordNode.items():
            if len(node.Neighbors) != 2 * self.NDim:
                return False
        return True

    def GetMagnification(self, CoordNode, NodeCoord):
        magnification = []
        bFoundMag = False
        #import pdb; pdb.set_trace()
        for inode, icoord in NodeCoord.items():
            for jnbr in self.NodeVec[inode].Neighbors:
                if jnbr.Id in NodeCoord.keys():
                    bFoundMag = True
                    magnification.append(jnbr.RingMagnification)
        if not bFoundMag:
            import pdb; pdb.set_trace("Hey, why no magnification/color?")
            return None
        if len(set(magnification)) > 2:
            import pdb; pdb.set_trace()
            print("Hey, how come there are more than two magnification/color levels? " + str(CoordNode))

            return None
        
        # always return the lower magnification, but since it's a ring, cannot just do min and max

        return self.RingMin(magnification)

        
    def PrintCoords(self, coordnodes, multiplier=1):
        for coord in sorted(coordnodes.keys()):
            id = coordnodes[coord]
            thiscoord = tuple([round(x*multiplier,4) for x in self.NodeVec[id].Coords])
            print(coord, thiscoord, coordnodes[coord])  


    def Base3(self, i, padlength=-1):
        ix = copy(i)
        retvec = []
        while ix > 0:
            thispart = ix % 3
            retvec.append(thispart) # we'll reverse at the end
            ix = ix // 3        
        if padlength > 0 and len(retvec) < padlength:
            padding = padlength - len(retvec)
            retvec = retvec + ([0] * padding)
        retvec.reverse()
        return retvec

    def BaseN(self, i, N=4, padlength=-1):
        ix = copy(i)
        retvec = []
        if i == 0:
            if padlength != -1:
                return [0] * padlength
            else:
                return [0]
        while ix > 0:
            thispart = ix % N
            retvec.append(thispart) # we'll reverse at the end
            ix = ix // N     
        if padlength > 0 and len(retvec) < padlength:
            padding = padlength - len(retvec)
            retvec = retvec + ([0] * padding)
        retvec.reverse()
        return retvec



    # this allows you to iterate across an ndimensional hypercube
    # not used
    def iterateN(x, niter, onlyinorder=False):
        # https://www.geeksforgeeks.org/python-check-if-list-is-sorted-or-not/
        def isascending(a):
            return all(a[i] <= a[i + 1] for i in range(len(a) - 1))
        def isdescending(a):
            return all(b[i] >= b[i + 1] for i in range(len(b) - 1))


        n = len(x)

        for i in range(n**niter):
            vec = BaseN(i, n, niter)
            if onlyinorder and not(isascending(vec)):
                continue

            # for now, just print it out to verify it's looping correctly
            print(i, vec, [x[i] for i in vec])



    def SwapWithPreviousNeighbors(self, Id, theseaxes):

        bFoundNew = False
        swappedaxis = []
        for iax in theseaxes:
            if iax in self.NodeVec[Id].PrevNeighbors.keys():
                swappedaxis.append(self.NodeVec[Id].PrevNeighbors[iax])
                bFoundNew = True
            else:
                swappedaxis.append(iax)
        return swappedaxis, bFoundNew




    def FindNearest(self, xarr):
        #import pdb; pdb.set_trace()
        xarr = np.array(xarr)
        best = None
        mindist = np.sqrt(self.NDim * self.TorLen * self.TorLen)
        for inode in self.NodeVec:
            thisdist = xarr - np.array(inode.Coords)
            thisdist = np.sqrt(np.sum(thisdist * thisdist))
            if thisdist < mindist:
                mindist = copy(thisdist)
                best = inode.Id
        return best
        

    def ClearHistogram(self,):
        self.Hist = np.zeros(self.Hist.shape).astype("int")
        self.NodesLastHistogramUpdate = 0

    def UpdateHistogram(self, ):
        for inode in range(self.NodesLastHistogramUpdate, len(self.NodeVec)):
            x,y = self.NodeVec[inode].Coords
            self.Hist[int(np.floor(x)), int(np.floor(y))] += 1
        self.NodesLastHistogramUpdate = copy(self.NNode)

    def HistogramPercentile(self,):
        def gethist(thisvec):
            return self.Hist[thisvec[0], thisvec[1]]

        gridsort = []
        for j in range(self.TorLen ** self.NDim):
            thisvec = self.BaseN(j, self.TorLen, self.NDim)
            gridsort.append( [gethist(thisvec), j, tuple(thisvec)] )
        
        gridsort = sorted(gridsort, key=lambda x: x[0])
        retpoints = {}
        for i in[0, 0.5, 0.8, 0.9, 0.95]:
            islot = int(i * len(gridsort))
            vec = np.array(gridsort[islot][2]) + 0.5
            retpoints[i] = (gridsort[islot][2], self.FindNearest(gridsort[islot][2]), self.FindNearest(vec), gridsort[islot][0])

        


        gridflat = self.Hist.flatten()
        gridflat.sort()
        cumgridflat = np.cumsum(gridflat)/float(np.sum(gridflat))
        maxvar = float(np.max(gridflat))
        N = len(gridflat)
        xarr = np.arange(len(cumgridflat))/float(N)
        #import pdb; pdb.set_trace()
        print("percentiles")
        for i in np.arange(0,100,0.5):
            pctile = np.round(np.percentile(gridflat, i),4)
            cumpct = np.interp(i/100.0, xarr, cumgridflat)
            print(i, pctile, cumpct)
        print("paretofrac")
        for i in np.arange(0,100,0.5):
            pctile = np.round(np.percentile(gridflat, i),4)/maxvar
            cumpct = np.interp(i/100.0, xarr, cumgridflat)
            revcumpct = 1.0 - np.interp(i/100.0, xarr, cumgridflat)
            print(i, pctile, cumpct, revcumpct)    
        
        mu = np.mean(gridflat)
        sig = np.std(gridflat)
        var = np.var(gridflat)
        alpha = mu * mu / var
        theta = var / mu
        print("GAMMA PARAMS: N %d mean %f std %f var %f alpha %f theta %f" % (N, mu, sig, var, alpha, theta))
        print(retpoints)
        return retpoints



    def ExpandManyRandomly(self, Prob, NRuns, ProbDiagonalize, listnslices):
        #print("dont forget to randomize nslices")

        #print("Assign sierpinski level to try and even out the ")
        global FreshProb
        global MaxNodes

        self.ClearHistogram()
        ndivide = 0

        nchunk = 1000
        nlastchunk = 0

        histmean = [0]
        histstd = [0]

        startnnode = copy(self.NNode)

        lastprintscreen = 0
        ii = 0
        for i in range(NRuns):
            if len(listnslices) > 1:
                nslices = rn.choice(listnslices)
            else:
                nslices = listnslices[0]

            print(i)
            #print("CHANGE THIS BACK")
            #thisnode = self.NodeVec[0] #rn.choice(self.NodeVec)
            #thisnode = rn.choice(self.NodeVec)
            
            #inode = thisnode.Id
            #if self.NNode >= 2324:
            #    print("Y", i, inode, self.NNode, self.NodeVec[2324].NeighborIds())
            #elif self.NNode >= 1182:
            #    print("Y", i, inode, self.NNode, self.NodeVec[1181].NeighborIds())
            #else:
            #    print("Y", i, inode, self.NNode)

            toprun = copy(self.NNode)
            bottomrun = 0
            if i == 0:
                bottomrun = int(rn.random()*(self.TorLen ** self.NDim))
            for inode in range(bottomrun, toprun):
                thisnode = self.NodeVec[inode]
                #inode = thisnode.Id


                #if i == 159:
                #    print("159")
                #    import pdb; pdb.set_trace()
                
                if len(thisnode.Neighbors) != 2 * self.NDim:
                    #print("node", inode, " has less than 2*D neighbors and therefore abuts a lower-mag cube. this attempt fails.")
                    continue
                #currentnodes = list(self.NodeDict.keys())
                
                #if thisnode.Parity != 1 and ((nslices % 2)==0):
                #    print("node %d has different parity than what we want, so we're skipping it" % (inode,))
                #    continue

                if rn.random() > Prob:
                    continue

                if len(thisnode.Neighbors) == 0:
                    continue
                if thisnode.bDefunct:
                    print("node %d is defunct" % (inode,))
                    continue


                #if thisnode.Id == 885: # and tuple(iax) == (2025, 2021):
                #        print("def bad")
                #        import pdb; pdb.set_trace()

                bUsePrevNeighborsToo = True





                allaxes = self.ReturnAllAxes(inode) # this doesn't test if any of the axes are usable, and only gets rid of combos that have collinear directions

                if len(allaxes) == 0:
                    #print("No good axes at point", inode)
                    continue

                """

                coords, noodes, mag = ggrid.PickADivisibleCube(thisnode, allaxes)
                ggrid.DivideCube(startnode01, coords, noodes, nslices)

                """

                #if i >= 1 and 66 in thisnode.NeighborIds():
                #    print("this is neighbor of 66", inode, i)
                #    import pdb; pdb.set_trace()



                coords, noodes, iax = self.PickADivisibleCube(inode, allaxes)






                if not(coords is None):

    
                    
                    tmpcrd = []
                    for crd,nd in sorted(coords.items()):
                        x = np.array(self.NodeVec[nd].Coords)*27
                        tmpcrd.append(x)
                        #print(nd, crd, x, self.NNode, iax)
                    
                    sumx = tmpcrd[0]
                    for ix in tmpcrd[1:]:
                        sumx = sumx + ix
                    meanx = sumx/len(tmpcrd)
                    
                    errx = []
                    
                    for ix in tmpcrd:
                        errx.append(  np.linalg.norm(ix-meanx)  )
                    
                    errx = np.array(errx)
                    relerrx = np.abs((errx - np.mean(errx))/errx)
                    
                    if np.max(relerrx) > 0.001:
                        #import pdb; pdb.set_trace()
                        print("The coords are getting very close", errx, relerrx)
                        #for icrd,inode in coords:
                        #    print(inode, self.NodeVec[inode].Coords) 

                    #print("NOW DIVIDE")
                    #import pdb; pdb.set_trace()   

                    self.DivideCube(inode, coords, noodes, nslices)
                    if self.NNode > MaxNodes:
                        print("reached maxnode upper limit", self.NNode)
                        break
                    if ii - lastprintscreen > 2048:
                        #self.ScatterPlot()
                        lastprintscreen = copy(ii)
                    ii  += 1


                    ndivide += 1
                    if self.NNode % 1000 == 0:


                        sum0 = 0
                        n = 0
                        for inode in self.NodeVec:
                            for inbr in inode.Neighbors:
                                n += 1
                                #if inbr.RingMagnification == 0:
                                #    sum0 += 1

                        if self.NNode >= self.TorLen ** self.NDim:
                            #print("Here's a division")
                            #import pdb; pdb.set_trace()
                            pass
                        print("Dividing", inode.Id, coords, self.NNode, sum0/2.0/n, sum0, ndivide, inode.Coords)
                        if ndivide >= nlastchunk + nchunk:
                            self.UpdateHistogram()
                            histmean.append( np.mean(self.Hist) )
                            histstd.append( np.std(self.Hist) )


                            nlastchunk = (ndivide // nchunk) * nchunk


                else:
                    #print("No good coords at point ", inode, " with axes", allaxes)
                    continue                


        if self.NNode - startnnode > (nslices + 1) ** self.NDim:
            print("HISTOGRAM")
            for i in range(len(histmean)):
                print(i, histmean[i], histstd[i])
            myretval = self.HistogramPercentile()

            return myretval
        else:
            return {}






    def ExpandTopCheckerboard(self, nslices):
        print("dont forget to randomize nslices")

        if self.NDim != 2:
            print("Only works for 2-dimensions")
            return

        
        packet = ([0, [1, 8]], 
                  [1, [2, 57]],
                  [2, [3, 10]],
                  [3, [4, 59]],
                  [4, [5, 12]],
                  [5, [6, 61]],
                  [6, [7, 14]],
                  [7, [0, 63]],
                  [16, [17, 24]],
                  [17, [ 9, 18]],
                  [18, [19, 26]],
                  [19, [11, 20]],
                  [20, [21, 28]],
                  [21, [13, 22]],
                  [20, [23, 30]],
                  [23, [15, 16]],
                  [32, [33, 40]],
                  [33, [25, 34]],
                  [34, [35, 42]],
                  [35, [27, 36]],
                  [36, [37, 44]],
                  [37, [29, 38]],
                  [38, [39, 46]],
                  [39, [31, 32]],
                  [48, [49, 56]],
                  [49, [41, 50]],
                  [50, [51, 58]],
                  [51, [43, 52]],
                  [52, [53, 60]],
                  [53, [45, 54]],
                  [54, [55, 62]],
                  [55, [47, 48]])





        #import pdb; pdb.set_trace()

        for inode, theseaxes in packet:

                (CoordNode, NodeCoord) = self.FindCubeCoord(inode, theseaxes)
                if CoordNode is None:
                    continue
                if False: #thisrand > Prob:
                    #import pdb; pdb.set_trace()
                    self.DiagonalizeCube(inode, CoordNode, NodeCoord)
                else:
                    #import pdb; pdb.set_trace()
                    print("div ", inode, CoordNode, NodeCoord, len(self.NodeVec))
                    self.DivideCube(inode, CoordNode, NodeCoord, nslices)

                    self.ClearParity()
                    self.MakeParity(0)
                    if len(self.NodeVec) > 600**self.NDim:
                        print("maxed out")
                        return
    def NNeighborHistogram(self,):

        hist = {}
        for i in range(2*self.NDim + 1):
            hist[i] = 0
        for inode in self.NodeVec:
            nnbr = len(inode.Neighbors)
            hist[nnbr] += 1
        
        for key in sorted(hist):
            val = hist[key]
            #print("%d neighbors: %d %f" % (key, val, val/float(len(self.NodeVec))))


    def ExpandTopCheckerboardGeneral(self, coordnode, nodecoord, nslices):
        print("dont forget to randomize nslices")

        def GTby1(vecB, vecA, cubelen):
            nbigger = 0
            for icompa, icompb in zip(vecA, vecB):
                if icompb == icompa + 1 or (icompb == icompa + 1 + cubelen):
                    nbigger += 1
                    continue
                elif  icompb == icompa:
                    continue
                else:
                    return False
            if nbigger == 1:
                return True
            return False


        
        Norig = copy(self.NNode)

        self.ClearParity()
        self.MakeParity(0)

        cubelen = int(np.round(Norig**(1.0/self.NDim)))







        #import pdb; pdb.set_trace()

        for iinode in sorted(nodecoord.keys()):

                #import pdb; pdb.set_trace()
      
                inode = self.NodeVec[iinode]

                if inode.Parity != 1:
                    continue
                maxneighbor = max(inode.Neighbors)
                #import pdb; pdb.set_trace()
                if maxneighbor >= len(coordnode):
                    #import pdb; pdb.set_trace()
                    continue

                nbrs = inode.Neighbors

                thiscoord = tuple(nodecoord[inode.Id])

                axes = []

                for inbr in nbrs:
                    #import pdb; pdb.set_trace()
                    try:
                        if GTby1(nodecoord[inbr], thiscoord, cubelen):
                            axes.append(inbr)
                    except:
                        import pdb; pdb.set_trace()


                
                (CoordNode, NodeCoord) = self.FindCubeCoord(inode.Id, axes)
                if CoordNode is None:
                    print("Why does this set of axes fail?")
                    import pdb; pdb.set_trace()
                    continue
                if False: #thisrand > Prob:
                    #import pdb; pdb.set_trace()
                    self.DiagonalizeCube(inode, CoordNode, NodeCoord)
                else:                    
                    print("OK", inode.Id, len(CoordNode))
                    #if inode.Id == 6:
                    #    import pdb; pdb.set_trace()
                    #print("div ", inode, CoordNode, NodeCoord, len(self.NodeVec))
                    self.DivideCube(inode, CoordNode, NodeCoord, nslices)
                    #import pdb; pdb.set_trace()





    # not good for anything
    def NodesWithinRadius(self, Id, Radius):
        # for any generated path, include the path length at Radis-2, Radius-4, Radius-6....
        NPath = self.NDimPlus1 ** Radius

        base = self.NDimPlus1
        for i in range(NPath):
            
            igen = copy(i)

            ipath = []
            for istep in range(Radius):
                ipath.append(igen % base)
                igen = igen // base

    # earlier badder version
    def HeatEqBAD(self, Id, T):
        # start with a 1 at point Id and then crunch the head equation for T steps
        NPath = self.NDimPlus1 ** Radius

        base = self.NDimPlus1
        for i in range(NPath):
            
            igen = copy(i)

            ipath = []
            for istep in range(Radius):
                ipath.append(igen % base)
                igen = igen // base


    def PickADivisibleCube(self, startnodeId, allaxes = None):
        """ 
        This returns the magnification of the cube that is to be divided;
        after the division, the smaller subcubes will have a larger magnification
        """


        if allaxes is None:
            allaxes = self.ReturnAllAxes(startnodeId)
        #import pdb; pdb.set_trace()
        resallax = []
        alreadyrejected = []
        alreadyaccepted = []
        for iax in allaxes:
            #import pdb; pdb.set_trace()

            #if startnodeId == 16905 and (tuple(iax) == (32581, 19901) or tuple(iax) == (19901, 32581)):
            #    print("hhh")
            #    import pdb; pdb.set_trace()


            testiax = copy(iax)
            testiax.sort()
            if tuple(testiax) in alreadyrejected:
                continue
            if tuple(testiax) in alreadyaccepted:
                # note that flipping order around permutes the axes of the coordinates, but since they're
                # not retained, and since all it does is permute the numbers assigned to the new nodes,
                # there's no point in retesting a set of axes if one of their permutations has already passed
                continue

            coords, noodes = self.FindCubeCoord(startnodeId, iax)
            #print (iax, coords, noodes, mag)
            #import pdb; pdb.set_trace()
            if not(coords is None): # and not(mag is None):
                resallax.append((coords, noodes, iax))
                thisiax = copy(iax)
                thisiax.sort()
                thisiax = tuple(iax)
                if not(thisiax) in alreadyaccepted:
                    alreadyaccepted.append(thisiax)
            else:
                reject_iax = copy(iax)
                reject_iax.sort()
                alreadyrejected.append(tuple(reject_iax))                
        
        if len(resallax) == 0:
            return None, None, None


        
        return( rn.choice(resallax) )

    
    def NodeFromCoords(self, coordtup, err=1.0e-5):
        x0 = np.array((coordtup))
        for i, ix in enumerage(x0):
            if ix < 0:
                x[ix] += self.TorLen


        myretval = []
        for inode in self.NodeVec:
            idist = np.linalg.norm(x0, np.array(inode.Coords))
            if idist <= err:
                myretval.append(inode.Id)
        return myretval
    
    def GetAxisNodes(self, Id):
        def DifferWhere(x,y):
            myretval = []
            for i in range(len(x)):
                if x[i] != y[i]:
                    myretval.append(i)
            return myretval
        nbrs = self.NodeVec[Id].Neighbors
        IdCoord = self.NodeVec[Id].Coords

        myretval = [None for i in range(self.NDim)]
        for inbr in nbrs:
            inbrcoord = self.NodeVec[inbr].Coords
            whichax = DifferWhere(inbrcoord, IdCoord)
            if len(whichax) == 1:
                nominaldist = inbrcoord[whichax[0]] - IdCoord[whichax[0]]
                if nominaldist > 1.0: # we know the max nbr dist is 1
                    nominaldist -= self.TorLen     
                if nominaldist > 0:
                    myretval[whichax[0]] = inbr
        return myretval



    def ScatterPlot(self, xminxmax=None, yminymax=None, bEcho=None):
        # even though this is designed for 2-d graphs, I will make the data-acquisition loop general

        dataarr = np.zeros((self.NDim, len(self.NodeVec)))

        iinput = 0
        for ii,inode in enumerate(self.NodeVec):
            thisx, thisy = inode.Coords
            if not(xminxmax is None):
                if thisx < xminxmax[0] or thisx > xminxmax[1]:
                    continue
            if not(yminymax is None):
                if thisy < yminymax[0] or thisy > yminymax[1]:
                    continue           
            if not(bEcho is None):
                print(inode.Id, inode.Coords, inode.Neighbors)     
            dataarr[:,iinput] = inode.Coords
            iinput += 1


        x = dataarr[0,:]
        y = dataarr[1,:]
        s = 4
        plt.scatter(x, y, s)
        if not(xminxmax is None):
            plt.axis([xminxmax[0], xminxmax[1], yminymax[0], yminymax[1]])
        plt.show()
    

    def ScatterPlotLabeled(self, xminxmax, yminymax, bEcho=None):
        # even though this is designed for 2-d graphs, I will make the data-acquisition loop general

        
        dataarr = []
        labels = []
        iinput = 0
        for ii,inode in enumerate(self.NodeVec):
            thisx, thisy = inode.Coords
            
            if xminxmax[0] <= 0 and thisx > self.TorLen//2:
                thisx -= self.TorLen
            if yminymax[0] <= 0 and thisy > self.TorLen//2:
                thisy -= self.TorLen
                
            if xminxmax[1] >= self.TorLen//2 and thisx < np.min([xminxmax[0], self.TorLen//2]) :
                thisx += self.TorLen
            if yminymax[1]  >= self.TorLen//2 and thisy < np.min([yminymax[0], self.TorLen//2]) :
                thisy += self.TorLen
                



            thislabel = str(inode.Id)

            if thisx < xminxmax[0] or thisx > xminxmax[1]:
                continue
            if thisy < yminymax[0] or thisy > yminymax[1]:
                continue           
            if not(bEcho is None):
                print(inode.Id, inode.Coords, inode.Neighbors)     
            dataarr.append( (thisx, thisy) )
            labels.append( thislabel )
            iinput += 1


        x = [icomp[0] for icomp in dataarr]
        y = [icomp[1] for icomp in dataarr]
        fig, ax = plt.subplots()
        s = 4
        ax.scatter(x, y, s)
        for i, txt in enumerate(labels):
            ax.annotate(txt, (x[i], y[i]))

        plt.show()
        

    def ScatterPlot3(self, xminxmax=None, yminymax=None, zminzmax=None, bEcho=None):
        # even though this is designed for 2-d graphs, I will make the data-acquisition loop general

        dataarr = np.zeros((self.NDim, len(self.NodeVec)))

        iinput = 0
        for ii,inode in enumerate(self.NodeVec):
            thisx, thisy, thisz = inode.Coords
            if not(xminxmax is None):
                if thisx < xminxmax[0] or thisx > xminxmax[1]:
                    continue
            if not(yminymax is None):
                if thisy < yminymax[0] or thisy > yminymax[1]:
                    continue           
            if not(zminzmax is None):
                if thisz < zminzmax[0] or thisz > zminzmax[1]:
                    continue   
            if not(bEcho is None):
                print(inode.Id, inode.Coords, inode.Neighbors)     
            dataarr[:,iinput] = inode.Coords
            iinput += 1


        x = dataarr[0,:]
        y = dataarr[1,:]
        z = dataarr[2,:]



        from mpl_toolkits.mplot3d import Axes3D


        # Create a figure and 3D axis
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111, projection='3d')

        # Create scatter plot
        ax.scatter3D(x, y, z, color='red', marker='o')

        # Labels
        ax.set_xlabel('X Axis')
        ax.set_ylabel('Y Axis')
        ax.set_zlabel('Z Axis')
        ax.set_title('Cubes x:%s, y:%s z:%s' % (str(xminxmax), str(yminymax), str(zminzmax)))


        plt.show()
        



    def ScatterPlotLabeled3(self, xminxmax, yminymax, zminzmax, bEcho=None):
        # even though this is designed for 2-d graphs, I will make the data-acquisition loop general
        # https://stackoverflow.com/questions/10374930/annotating-a-3d-scatter-plot
        
        dataarr = []
        labels = []
        iinput = 0
        for ii,inode in enumerate(self.NodeVec):
            thisx, thisy, thisz = inode.Coords
            
            if xminxmax[0] <= 0 and thisx > self.TorLen//2:
                thisx -= self.TorLen
            if yminymax[0] <= 0 and thisy > self.TorLen//2:
                thisy -= self.TorLen
            if zminzmax[0] <= 0 and thisz > self.TorLen//2:
                thisz -= self.TorLen
                
            if xminxmax[1] >= self.TorLen//2 and thisx < np.min([xminxmax[0], self.TorLen//2]) :
                thisx += self.TorLen
            if yminymax[1]  >= self.TorLen//2 and thisy < np.min([yminymax[0], self.TorLen//2]) :
                thisy += self.TorLen
            if zminzmax[1]  >= self.TorLen//2 and thisy < np.min([zminzmax[0], self.TorLen//2]) :
                thisz += self.TorLen
                



            thislabel = str(inode.Id)

            if thisx < xminxmax[0] or thisx > xminxmax[1]:
                continue
            if thisy < yminymax[0] or thisy > yminymax[1]:
                continue           
            if thisz < zminzmax[0] or thisz > zminzmax[1]:
                continue           
            if not(bEcho is None):
                print(inode.Id, inode.Coords, inode.Neighbors)     
            dataarr.append( (thisx, thisy, thisz) )
            labels.append( thislabel )
            iinput += 1


        x = [icomp[0] for icomp in dataarr]
        y = [icomp[1] for icomp in dataarr]
        z = [icomp[2] for icomp in dataarr]



        from mpl_toolkits.mplot3d import Axes3D


        # Create a figure and 3D axis
        s = 10
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')

        for i in range(len(x)): #plot each point + it's index as text above
            ax.scatter(x[i],y[i],z[i],color='b') 
            ax.text(x[i],y[i],z[i],  '%s' % (labels[i]), size=s, zorder=1,  
            color='k') 
        plt.show()



    def ScatterPlotLabeledColored3(self, xminxmax, yminymax, zminzmax, bEcho=None):
        # even though this is designed for 2-d graphs, I will make the data-acquisition loop general
        # https://stackoverflow.com/questions/10374930/annotating-a-3d-scatter-plot
        
        dataarr4 = []
        labels4 = []
        dataarr5 = []
        labels5 = []
        dataarr6 = []
        labels6 = []
        iinput = 0
        for ii,inode in enumerate(self.NodeVec):
            thisx, thisy, thisz = inode.Coords
            
            if xminxmax[0] <= 0 and thisx > self.TorLen//2:
                thisx -= self.TorLen
            if yminymax[0] <= 0 and thisy > self.TorLen//2:
                thisy -= self.TorLen
            if zminzmax[0] <= 0 and thisz > self.TorLen//2:
                thisz -= self.TorLen
                
            if xminxmax[1] >= self.TorLen//2 and thisx < np.min([xminxmax[0], self.TorLen//2]) :
                thisx += self.TorLen
            if yminymax[1]  >= self.TorLen//2 and thisy < np.min([yminymax[0], self.TorLen//2]) :
                thisy += self.TorLen
            if zminzmax[1]  >= self.TorLen//2 and thisy < np.min([zminzmax[0], self.TorLen//2]) :
                thisz += self.TorLen
                

            thislabel = str(inode.Id)

            if thisx < xminxmax[0] or thisx > xminxmax[1]:
                continue
            if thisy < yminymax[0] or thisy > yminymax[1]:
                continue           
            if thisz < zminzmax[0] or thisz > zminzmax[1]:
                continue           
            if not(bEcho is None):
                print(inode.Id, inode.Coords, inode.Neighbors)     
            
            nbrlen = len(inode.Neighbors)
            if nbrlen == 6:
                dataarr6.append( (thisx, thisy, thisz) )
                labels6.append( thislabel )
            elif nbrlen == 5:
                dataarr5.append( (thisx, thisy, thisz) )
                labels5.append( thislabel )
            elif nbrlen == 4:
                dataarr4.append( (thisx, thisy, thisz) )
                labels4.append( thislabel )
            iinput += 1


        x6 = [icomp[0] for icomp in dataarr6]
        y6 = [icomp[1] for icomp in dataarr6]
        z6 = [icomp[2] for icomp in dataarr6]


        x5 = [icomp[0] for icomp in dataarr5]
        y5 = [icomp[1] for icomp in dataarr5]
        z5 = [icomp[2] for icomp in dataarr5]


        x4 = [icomp[0] for icomp in dataarr4]
        y4 = [icomp[1] for icomp in dataarr4]
        z4 = [icomp[2] for icomp in dataarr4]



        from mpl_toolkits.mplot3d import Axes3D


        # Create a figure and 3D axis
        s = 8
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')

        for i in range(len(x6)): #plot each point + it's index as text above
            ax.scatter(x6[i],y6[i],z6[i],color='b') 
            ax.text(x6[i],y6[i],z6[i],  '%s' % (labels6[i]), size=s, zorder=1,  
            color='k') 
        for i in range(len(x5)): #plot each point + it's index as text above
            ax.scatter(x5[i],y5[i],z5[i],color='r') 
            ax.text(x5[i],y5[i],z5[i],  '%s' % (labels5[i]), size=s, zorder=1,  
            color='k') 
        for i in range(len(x4)): #plot each point + it's index as text above
            ax.scatter(x4[i],y4[i],z4[i],color='g') 
            ax.text(x4[i],y4[i],z4[i],  '%s' % (labels4[i]), size=s, zorder=1,  
            color='k') 


        plt.show()



    def ScatterPlotLabeledColoredPacked3(self, xminxmax, yminymax, zminzmax, bEcho=None):
        # even though this is designed for 2-d graphs, I will make the data-acquisition loop general
        # https://stackoverflow.com/questions/10374930/annotating-a-3d-scatter-plot
        
        dataarr4 = []
        labels4 = []
        dataarr5 = []
        labels5 = []
        dataarr6 = []
        labels6 = []
        iinput = 0
        for ii,inode in enumerate(self.NodeVec):
            if len(inode.Neighbors) != 2 * self.NDim:
                continue

            thisx, thisy, thisz = inode.Coords
            
            if xminxmax[0] <= 0 and thisx > self.TorLen//2:
                thisx -= self.TorLen
            if yminymax[0] <= 0 and thisy > self.TorLen//2:
                thisy -= self.TorLen
            if zminzmax[0] <= 0 and thisz > self.TorLen//2:
                thisz -= self.TorLen
                
            if xminxmax[1] >= self.TorLen//2 and thisx < np.min([xminxmax[0], self.TorLen//2]) :
                thisx += self.TorLen
            if yminymax[1]  >= self.TorLen//2 and thisy < np.min([yminymax[0], self.TorLen//2]) :
                thisy += self.TorLen
            if zminzmax[1]  >= self.TorLen//2 and thisy < np.min([zminzmax[0], self.TorLen//2]) :
                thisz += self.TorLen
                

            thislabel = str(inode.Id)

            if thisx < xminxmax[0] or thisx > xminxmax[1]:
                continue
            if thisy < yminymax[0] or thisy > yminymax[1]:
                continue           
            if thisz < zminzmax[0] or thisz > zminzmax[1]:
                continue           
            if not(bEcho is None):
                print(inode.Id, inode.Coords, inode.Neighbors)     
            
            nbrlen = len(inode.Neighbors)
            if nbrlen == 6:
                dataarr6.append( (thisx, thisy, thisz) )
                labels6.append( thislabel )
            elif nbrlen == 5:
                dataarr5.append( (thisx, thisy, thisz) )
                labels5.append( thislabel )
            elif nbrlen == 4:
                dataarr4.append( (thisx, thisy, thisz) )
                labels4.append( thislabel )
            iinput += 1


        x6 = [icomp[0] for icomp in dataarr6]
        y6 = [icomp[1] for icomp in dataarr6]
        z6 = [icomp[2] for icomp in dataarr6]


        x5 = [icomp[0] for icomp in dataarr5]
        y5 = [icomp[1] for icomp in dataarr5]
        z5 = [icomp[2] for icomp in dataarr5]


        x4 = [icomp[0] for icomp in dataarr4]
        y4 = [icomp[1] for icomp in dataarr4]
        z4 = [icomp[2] for icomp in dataarr4]



        from mpl_toolkits.mplot3d import Axes3D


        # Create a figure and 3D axis
        s = 8
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')

        for i in range(len(x6)): #plot each point + it's index as text above
            ax.scatter(x6[i],y6[i],z6[i],color='b') 
            ax.text(x6[i],y6[i],z6[i],  '%s' % (labels6[i]), size=s, zorder=1,  
            color='k') 
        for i in range(len(x5)): #plot each point + it's index as text above
            ax.scatter(x5[i],y5[i],z5[i],color='r') 
            ax.text(x5[i],y5[i],z5[i],  '%s' % (labels5[i]), size=s, zorder=1,  
            color='k') 
        for i in range(len(x4)): #plot each point + it's index as text above
            ax.scatter(x4[i],y4[i],z4[i],color='g') 
            ax.text(x4[i],y4[i],z4[i],  '%s' % (labels4[i]), size=s, zorder=1,  
            color='k') 


        plt.show()


    def ScatterPlotColored3(self, xminxmax, yminymax, zminzmax, bEcho=None):
        # even though this is designed for 2-d graphs, I will make the data-acquisition loop general
        # https://stackoverflow.com/questions/10374930/annotating-a-3d-scatter-plot
        
        dataarr4 = []
        labels4 = []
        dataarr5 = []
        dataarr6 = []
        iinput = 0
        for ii,inode in enumerate(self.NodeVec):
            thisx, thisy, thisz = inode.Coords
            
            if xminxmax[0] <= 0 and thisx > self.TorLen//2:
                thisx -= self.TorLen
            if yminymax[0] <= 0 and thisy > self.TorLen//2:
                thisy -= self.TorLen
            if zminzmax[0] <= 0 and thisz > self.TorLen//2:
                thisz -= self.TorLen
                
            if xminxmax[1] >= self.TorLen//2 and thisx < np.min([xminxmax[0], self.TorLen//2]) :
                thisx += self.TorLen
            if yminymax[1]  >= self.TorLen//2 and thisy < np.min([yminymax[0], self.TorLen//2]) :
                thisy += self.TorLen
            if zminzmax[1]  >= self.TorLen//2 and thisy < np.min([zminzmax[0], self.TorLen//2]) :
                thisz += self.TorLen
                

            thislabel = str(inode.Id)

            if thisx < xminxmax[0] or thisx > xminxmax[1]:
                continue
            if thisy < yminymax[0] or thisy > yminymax[1]:
                continue           
            if thisz < zminzmax[0] or thisz > zminzmax[1]:
                continue           
            if not(bEcho is None):
                print(inode.Id, inode.Coords, inode.Neighbors)     
            
            nbrlen = len(inode.Neighbors)
            if nbrlen == 6:
                dataarr6.append( (thisx, thisy, thisz) )

            elif nbrlen == 5:
                dataarr5.append( (thisx, thisy, thisz) )

            elif nbrlen == 4:
                dataarr4.append( (thisx, thisy, thisz) )

            iinput += 1


        x6 = [icomp[0] for icomp in dataarr6]
        y6 = [icomp[1] for icomp in dataarr6]
        z6 = [icomp[2] for icomp in dataarr6]


        x5 = [icomp[0] for icomp in dataarr5]
        y5 = [icomp[1] for icomp in dataarr5]
        z5 = [icomp[2] for icomp in dataarr5]


        x4 = [icomp[0] for icomp in dataarr4]
        y4 = [icomp[1] for icomp in dataarr4]
        z4 = [icomp[2] for icomp in dataarr4]



        from mpl_toolkits.mplot3d import Axes3D


        # Create a figure and 3D axis
        s = 8
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')

        for i in range(len(x6)): #plot each point + it's index as text above
            ax.scatter(x6[i],y6[i],z6[i],color='b') 
            
        for i in range(len(x5)): #plot each point + it's index as text above
            ax.scatter(x5[i],y5[i],z5[i],color='r') 
            
        for i in range(len(x4)): #plot each point + it's index as text above
            ax.scatter(x4[i],y4[i],z4[i],color='g') 
            


        plt.show()
       
    def Scat(self, Id, dist, wLabels=True):
        x = np.array(self.NodeVec[Id].Coords)        
        if self.NDim == 2:
            xx = (x[0] - 0.5*dist, x[0] + 0.5*dist)
            yy = (x[1] - 0.5*dist, x[1] + 0.5*dist)
            if wLabels:
                self.ScatterPlotLabeled(xx, yy)
            else:
                self.ScatterPlot(xx, yy)
        if self.NDim == 3:
            xx = (x[0] - 0.5*dist, x[0] + 0.5*dist)
            yy = (x[1] - 0.5*dist, x[1] + 0.5*dist)
            zz = (x[2] - 0.5*dist, x[2] + 0.5*dist)            
            if wLabels:
                self.ScatterPlotLabeled3(xx, yy, zz)
            else:
                self.ScatterPlot3(xx, yy, zz)

        
    def ScatColored(self, Id, dist, wLabels=True):

        x = np.array(self.NodeVec[Id].Coords)

        if self.NDim == 2:
            xx = (x[0] - 0.5*dist, x[0] + 0.5*dist)
            yy = (x[1] - 0.5*dist, x[1] + 0.5*dist)
            if wLabels:
                self.ScatterPlotLabeled(xx, yy)
            else:
                self.ScatterPlot(xx, yy)
        if self.NDim == 3:
            xx = (x[0] - 0.5*dist, x[0] + 0.5*dist)
            yy = (x[1] - 0.5*dist, x[1] + 0.5*dist)
            zz = (x[2] - 0.5*dist, x[2] + 0.5*dist)
            if wLabels:
                self.ScatterPlotLabeledColored3(xx, yy, zz)
            else:
                self.ScatterPlotColored3(xx, yy, zz)

        
    def ScatColoredPacked(self, Id, dist):

        x = np.array(self.NodeVec[Id].Coords)

        if self.NDim == 2:
            xx = (x[0] - 0.5*dist, x[0] + 0.5*dist)
            yy = (x[1] - 0.5*dist, x[1] + 0.5*dist)
            if wLabels:
                self.ScatterPlotLabeled(xx, yy)
            else:
                self.ScatterPlot(xx, yy)
        if self.NDim == 3:
            xx = (x[0] - 0.5*dist, x[0] + 0.5*dist)
            yy = (x[1] - 0.5*dist, x[1] + 0.5*dist)
            zz = (x[2] - 0.5*dist, x[2] + 0.5*dist)
            self.ScatterPlotLabeledColoredPacked3(xx, yy, zz)

        

    def Distributions(self,inode):
        for inod in ggrid.NodeVec:
            inod.Amplitude = 0
            inod.PrevAmplitude = 0

        #import pdb; pdb.set_trace()

        self.NodeVec[inode].PrevAmplitude = 1.0
        self.NodeVec[inode].Amplitude = 1.0

        print("Heat", ggrid.NNode)            
        self.HeatEq(250, inode, True)
        #ggrid.MonitoredHeatEq(60, inode, bPrint=True)

    
        for inod in ggrid.NodeVec:
            inod.Amplitude = 0
            inod.PrevAmplitude = 0

        #inode = 0
        self.NodeVec[inode].PrevAmplitude = 1.0
        self.NodeVec[inode].Amplitude = 1.0

        print("Touch", self.NNode)
        self.Touch(120)



def AllCombosNoReplacement(N, D):
    """
    all draws of 0...(N-1) elements of length D, no replacement
    """
    ncombos = 1
    for i in range(D):
        ncombos *= N-i

    myret = []
    #import pdb; pdb.set_trace()
    for i in range(ncombos):
        thisnumber = copy(i)
        vec = []
        thisrng = list(np.arange(N).astype("int"))
        for j in range(D):     
                   
            divsor = N - j
            thisdig = i % divsor
            vec.append(thisrng[thisdig])
            thisrng.pop(thisdig)            
            divsor = i // divsor
        #print(i, vec)
        myret.append(tuple(vec))
    return myret

def MoreStuff(ggrid, i, nneighbors, nsteps=60):
    igroup = [i]
    nloop = 0
    while nloop < 50:
        newgroup = []
        nloop  +=1
        for i in igroup:
            if len(ggrid.NodeVec[i].Neighbors) == nneighbors:
                for inod in ggrid.NodeVec:
                    inod.Amplitude = 0
                    inod.PrevAmplitude = 0


                ggrid.NodeVec[i].PrevAmplitude = 1.0
                ggrid.NodeVec[i].Amplitude = 1.0
                    

                ggrid.HeatEq(nsteps, i, True)
                return
            newgroup.extend(ggrid.NodeVec[i].Neighbors)
        igroup = copy(newgroup)
    
    print("Failed to find neighbor of len", nneighbors, "within 50 steps")
    return




"""


def NodesWithinRadius(Radius):
    id = 0
    NDimPlus1 = 2
    # for any generated path, include the path length at Radis-2, Radius-4, Radius-6....
    NPath = NDimPlus1 ** 2


    base = NDimPlus1
    for i in range(NPath):
        
        igen = copy(i)

        ipath = []
        for istep in range(Radius):
            ipath.append(igen % base)
            igen = igen // base
        print(i, ipath)



"""





def ProofNeighbors(ggrid):
    iloop = -1
    for id in range(ggrid.NNode):
        iloop += 1
        if iloop % 100 == 0:
            print(iloop)
        Nnbrs = 0
        oppneighbors = []
        for jnode in ggrid.NodeVec:
            for jnbr in jnode.Neighbors:
                if id == jnbr:
                    Nnbrs += 1
                    oppneighbors.append(jnode.Id)
        if Nnbrs != len(ggrid.NodeVec[id].Neighbors):
            print("WHOAH")
            import pdb; pdb.set_trace()
        oppneighbors.sort()
        if oppneighbors != ggrid.NodeVec[id].Neighbors:
            print("WHOAH2")
            import pdb; pdb.set_trace()


def getweightedfractality(ggrid):
    wtdsum = 0
    wtdsum2 = 0
    wtdamp = 0
    for inode in ggrid.NodeVec:
        wtdamp += inode.Amplitude
        #wtdnbr = np.array([jnode.Magnification for jnode in inode.Neighbors]).mean()
        #wtdnbr2 = np.array([jnode.Magnification * jnode.Magnification for jnode in inode.Neighbors]).mean()
        wtdsum += inode.Amplitude * wtdnbr
        wtdsum2 += inode.Amplitude * wtdnbr2

    wtdamp = 1.0 if wtdamp <= 0 else wtdamp
    return wtdsum/wtdamp, wtdsum2/wtdamp

if __name__ == '__main__':   
    
    opts, args = ReadParams()
    rn.seed(opts.seed)

    FreshProb = opts.xprob
    MaxNodes = opts.maxnodes

    TargetLength = opts.length

    GuideDict = {}
    NodeCoord = {}
    CoordNode = {}

    Prob = opts.prob
    NRuns = opts.expansionruns

    ggrid = cubenodeset_t(opts.dimension)
    
    ggrid.CreateNCube()

    """
    maxnbr = -1
    minnbr = 2 * ggrid.NDim + 1
    for inv in ggrid.NodeVec:
        if len(inv.Neighbors) < minnbr:
            minnbr = len(inv.Neighbors)
        if len(inv.Neighbors) > maxnbr:
            maxnbr = len(inv.Neighbors)

    import pdb; pdb.set_trace()
    """

    startnode00 = 0
    nslices = 2
    #nodecoord = {0:(0,0), 1:(0,1), 32:(1,0), 33:(1,1)}
    #coordnode = {}
    #for key,val in nodecoord.items():
    #    coordnode[val] = key

    if False:
        allaxes = ggrid.ReturnAllAxes(startnode00)

        
        targetax = allaxes[4] # should be [32,1]

        (coordnode, nodecoord) = ggrid.FindCubeCoord(0, targetax)
        
        ggrid.DivideCube(startnode00, coordnode, nodecoord, nslices)


        allaxes00 = ggrid.ReturnAllAxes(startnode00)
        
        startnode01 = 1 

        allaxes01 = ggrid.ReturnAllAxes(startnode01)

        
        coords, noodes, iax = ggrid.PickADivisibleCube(startnode01, allaxes01)
        if not(coords is None):
            ggrid.DivideCube(startnode01, coords, noodes, nslices)




        for iax in allaxes01:
            coords, noodes = ggrid.FindCubeCoord(startnode01, iax)
            print (iax, coords, noodes)

        nodecoord = {32:(1,0), 33:(1,1), 64:(2,0), 65:(2,1)}
        coordnode = {}
        for key,val in nodecoord.items():
            coordnode[val] = key



    bEnclosedCube = True

    if bEnclosedCube:


        GuideDict = {}
        NodeCoord = {}
        CoordNode = {}

        Prob = opts.prob
        NRuns = opts.expansionruns

        if opts.dimension == 2:

            ggrid = cubenodeset_t(opts.dimension)
            
            ggrid.CreateNCube()

            
            #targetax = allaxes[4] # should be [32,1]
            #import pdb; pdb.set_trace()

            i = 0
            print("start", ggrid.NNode); i += 1
            import pdb; pdb.set_trace()
            (coordnode, nodecoord) = ggrid.FindCubeCoord(32, [33,64])
            ggrid.DivideCube(coordnode, nodecoord, nslices)
            print(i, ggrid.NNode); i += 1
            (coordnode, nodecoord) = ggrid.FindCubeCoord(34, [35,66])
            ggrid.DivideCube(coordnode, nodecoord, nslices)
            print(i, ggrid.NNode); i += 1
            (coordnode, nodecoord) = ggrid.FindCubeCoord(1, [2,33])
            ggrid.DivideCube(coordnode, nodecoord, nslices)
            print(i, ggrid.NNode); i += 1
            (coordnode, nodecoord) = ggrid.FindCubeCoord(65, [66,97])
            ggrid.DivideCube(coordnode, nodecoord, nslices)
            print(i, ggrid.NNode); i += 1



            # this next section subdivides each of the 4 divided "petals" 
            if True:
                (coordnode, nodecoord) = ggrid.FindCubeCoord(1027, [1031,1028])
                ggrid.DivideCube(coordnode, nodecoord, nslices)

                (coordnode, nodecoord) = ggrid.FindCubeCoord(1051, [1055,1052])
                ggrid.DivideCube(coordnode, nodecoord, nslices)

                (coordnode, nodecoord) = ggrid.FindCubeCoord(1039, [1043,1040])
                ggrid.DivideCube(coordnode, nodecoord, nslices)

                (coordnode, nodecoord) = ggrid.FindCubeCoord(1063, [1067,1064])
                ggrid.DivideCube(coordnode, nodecoord, nslices)




            (coordnode, nodecoord) = ggrid.FindCubeCoord(33, [1053,1034])
            ggrid.DivideCube(coordnode, nodecoord, nslices)

            (coordnode, nodecoord) = ggrid.FindCubeCoord(33, [1049,1030])
            ggrid.DivideCube(coordnode, nodecoord, nslices)

            (coordnode, nodecoord) = ggrid.FindCubeCoord(1052, [1089,1053])
            ggrid.DivideCube(coordnode, nodecoord, nslices)
            (coordnode, nodecoord) = ggrid.FindCubeCoord(1053, [1120,1057])
            ggrid.DivideCube(coordnode, nodecoord, nslices)
            (coordnode, nodecoord) = ggrid.FindCubeCoord(1120, [1122,1121])
            ggrid.DivideCube(coordnode, nodecoord, nslices)
            (coordnode, nodecoord) = ggrid.FindCubeCoord(1034, [1120,1035])
            ggrid.DivideCube(coordnode, nodecoord, nslices)
            (coordnode, nodecoord) = ggrid.FindCubeCoord(1031, [1034,1082])
            ggrid.DivideCube(coordnode, nodecoord, nslices)
            print("same check")
            import pdb; pdb.set_trace()
            (coordnode, nodecoord) = ggrid.FindCubeCoord(1036, [1037,1039])
            ggrid.DivideCube(coordnode, nodecoord, nslices)

            # now do the central cube -- only one new node should be created
            startnode33 = 33

            allaxes33 = ggrid.ReturnAllAxes(startnode33)


            (coordnode, nodecoord) = ggrid.FindCubeCoord(startnode33, [34,65])
            ggrid.DivideCube(coordnode, nodecoord, nslices)
            print(ggrid.NNode)
            



            
            (coordnode, nodecoord) = ggrid.FindCubeCoord(32, [1026, 0])
            ggrid.DivideCube(coordnode, nodecoord, nslices)
            
            
            (coordnode, nodecoord) = ggrid.FindCubeCoord(startnode00, [1080, 1077])
            
            
            coordnode, nodecoord, mag = ggrid.PickADivisibleCube(startnode00)
            if not(coordnode is None):
                ggrid.DivideCube(coordnode, nodecoord, nslices)
        elif opts.dimension == 3:
            import pdb; pdb.set_trace()
            ggrid = cubenodeset_t(opts.dimension)
            ggrid.CreateNCube()

    #import pdb; pdb.set_trace()
    BigHist = copy(ggrid.Hist)
    for ibig in range(200):
        ggrid = cubenodeset_t(opts.dimension)
    
        ggrid.CreateNCube()
        myretval = ggrid.ExpandManyRandomly(Prob, NRuns, TargetLength, [1])
        BigHist = BigHist + ggrid.Hist
    
        print("just finished run", ibig)
        
        import pickle
        with open('blah2.pkl', 'wb') as fp:
            pickle.dump(BigHist, fp)
    np.savetxt('blah2.csv', BigHist)

    """
with open('blah2.pkl', 'rb') as fp:
    hist = pickle.load(fp)
hist.shape
hist[:3,:3]
np.percentile(hist, 90)
for i in range(1,100):
    print(np.percentile(hist,i))
x = np.arange(1,100)
y  = []
for i in range(1,100):
    y.append(np.percentile(hist,i))
len(y), len(x)
plt.plot(x,y)
plt.show()
z = [0]
for i in range(2,100):
    z.append(y[i]-y[i-1])
for i in range(1,len(y)):
    z.append(y[i]-y[i-1])
len(z)
z = [0]
for i in range(1,len(y)):
    z.append(y[i]-y[i-1])
len(z)
plt.plot(x,y)


z20 = np.array(z)*20
plt.plot(x,z20)
plt.show()
with open('blah2.pkl', 'rb') as fp:
    hist = pickle.load(fp)

    """

    bDebuggingEarly = True
    if bDebuggingEarly:    



        nsliceslist = [ nslices ]

        bGoForward = True
        bSkip = False
        if bGoForward:

            if False:
                ggrid.ExpandTopCheckerboardGeneral(CoordNode, NodeCoord, nsliceslist[0])
                print("The above creates a pretty uniform criss-cross of sierpinski-ness; if nslices is 1,")
                print(" then there are few nodes, and they  power law saturates quickly (so you have to )")
                print("ignore everything but the first few steps. For n=5 or more, the desired deminsionality is fairly evident so that")
                print(" these are pretty  close to 2 or 3 dimensions; I'm guessing they're not perfect due to roundoff error and some edge effects,")
                print("which hopefully diminish for very large surfaces or very high slice couts")


            ggrid = cubenodeset_t(opts.dimension)
            ggrid.CreateNCube()

            
            #import pdb; pdb.set_trace()
            myretval = ggrid.ExpandManyRandomly(Prob, NRuns, TargetLength, nsliceslist)
            #import pdb; pdb.set_trace()
            nodes = []
            for key,val in sorted(myretval.items()):
                nodes.append([key, val[0], val[1], val[2], val[3]])
            
            ggrid.ScatterPlot()
            import pdb; pdb.set_trace()
            for pct, icoord, inodeA, inodeB, nnodes in nodes:
                for inode in [inodeA, inodeB]:
                    for inod in ggrid.NodeVec:
                        inod.Amplitude = 0
                        inod.PrevAmplitude = 0

                    #import pdb; pdb.set_trace()

                    ggrid.NodeVec[inode].PrevAmplitude = 1.0
                    ggrid.NodeVec[inode].Amplitude = 1.0

                    print("for pct ", pct, "icoord", icoord, "inode", inode, "nnodesincube",nnodes)
                    print("Heat", ggrid.NNode)            
                    ggrid.HeatEq(250, inode, True)
                    #ggrid.MonitoredHeatEq(60, inode, bPrint=True)

                
                    for inod in ggrid.NodeVec:
                        inod.Amplitude = 0
                        inod.PrevAmplitude = 0

                    #inode = 0
                    ggrid.NodeVec[inode].PrevAmplitude = 1.0
                    ggrid.NodeVec[inode].Amplitude = 1.0

                    print("for pct ", pct, "icoord", icoord, "inode", inode, "nnodesincube",nnodes)
                    print("Touch", ggrid.NNode)
                    ggrid.Touch(120)

            
            for inod in ggrid.NodeVec:
                inod.Amplitude = 0
                inod.PrevAmplitude = 0

            inode = 0
            ggrid.NodeVec[inode].PrevAmplitude = 1.0
            ggrid.NodeVec[inode].Amplitude = 1.0


            import pdb; pdb.set_trace()

            #for i in [0, len(ggrid.NodeVec)//2, len(ggrid.NodeVec)//2+1, len(ggrid.NodeVec)-1]:
            heatresdict = {}
            touchresdict = {}

            for islice in  [(nslices,)]: #  [(2,),(4,),(2,4,6)]:         




                ggrid = cubenodeset_t(opts.dimension)
                ggrid.CreateNCube()

                #ggrid.ClearParity()

                #import pdb; pdb.set_trace()
                #ggrid.MakeParity(0)

                NodeCoordOrig = copy(NodeCoord)
                CoordNodeOrig = copy(CoordNode)


                #import pdb; pdb.set_trace()
                ggrid.ExpandManyRandomly(Prob, NRuns, TargetLength, islice)

                if False:
                    #import pdb; pdb.set_trace()
                    startnode01 = 41 # startnode10 = 32
                    import pdb; pdb.set_trace()
                    allaxes01 = ggrid.ReturnAllAxes(startnode01)   
                    coords,nodes = ggrid.FindCubeCoord(41, [7033,8829])    
                    coords, noodes, iax = ggrid.PickADivisibleCube(startnode01, allaxes01)
                    if not(coords is None):
                        ggrid.DivideCube(startnode01, coords, noodes, nslices)
                    import pdb; pdb.set_trace()

                    startnode01 = 32 # startnode10 = 32
                    import pdb; pdb.set_trace()
                    allaxes01 = ggrid.ReturnAllAxes(startnode01)       
                    coords, noodes, iax = ggrid.PickADivisibleCube(startnode01, allaxes01)
                    if not(coords is None):
                        ggrid.DivideCube(startnode01, coords, noodes, nslices)
                    import pdb; pdb.set_trace()
                    #import pdb; pdb.set_trace()

                    startnode01 = 63 # startnode10 = 32
                    import pdb; pdb.set_trace()
                    allaxes01 = ggrid.ReturnAllAxes(startnode01)       
                    coords, noodes, iax = ggrid.PickADivisibleCube(startnode01, allaxes01)
                    if not(coords is None):
                        ggrid.DivideCube(startnode01, coords, noodes, nslices)
                    import pdb; pdb.set_trace()

                    startnode01 = 32 # startnode10 = 32
                    import pdb; pdb.set_trace()
                    allaxes01 = ggrid.ReturnAllAxes(startnode01)       
                    coords, noodes, iax = ggrid.PickADivisibleCube(startnode01, allaxes01)
                    if not(coords is None):
                        ggrid.DivideCube(startnode01, coords, noodes, nslices)
                    import pdb; pdb.set_trace()

                    startnode01 = 7 # startnode10 = 32
                    #import pdb; pdb.set_trace()
                    allaxes01 = ggrid.ReturnAllAxes(startnode01)
                    coords, noodes, iax = ggrid.PickADivisibleCube(startnode01, allaxes01)
                    if not(coords is None):
                        ggrid.DivideCube(startnode01, coords, noodes, nslices)
                    import pdb; pdb.set_trace()


                    startnode01 = 64 # startnode10 = 32
                    #import pdb; pdb.set_trace()
                    allaxes01 = ggrid.ReturnAllAxes(startnode01)
                    coords, noodes, iax = ggrid.PickADivisibleCube(startnode01, allaxes01)
                    if not(coords is None):
                        ggrid.DivideCube(startnode01, coords, noodes, nslices)
                    #import pdb; pdb.set_trace()




                ggrid.ScatterPlot()
                import pdb; pdb.set_trace()
                #ggrid.ExpandMany(opts.prob, opts.expansionruns, 0, nslices) # opts.expansionruns)
                # ggrid.SubSimplices


                bCheckGridNodes = False # Takes a long time for big conglomerations!
                if bCheckGridNodes:
                    ProofNeighbors(ggrid)

                nnodes = len(ggrid.NodeVec)
                nnodesperdim = nnodes ** (1.0/ggrid.NDim)
                print("Generated grid has %d nodes, or roughly %.3f nodes in each of the %d genated dimensions" % ( nnodes, nnodesperdim, ggrid.NDim))

                ggrid.ClearParity()
                ggrid.MakeParity(0)


                nnodes = len(ggrid.NodeVec)
                
                nnodesperdim = nnodes ** (1.0/ggrid.NDim)
                #print("Generated grid has %d nodes, or roughly %.3f nodes in each of the %d genated dimensions" % ( nnodes, nnodesperdim, ggrid.NDim))
                
                print("defunct ", np.sum(np.array([ggrid.NodeVec[i].bDefunct for i in range(len(ggrid.NodeVec))])))
                #import pdb; pdb.set_trace()      

                # now analyze implies powers (and p-vals) over a bunch of nodes
                
                nsamples = 5
                nsamples = np.min((ggrid.NNode, nsamples//2))

                
                indsamples = rn.sample([i for i in range(ggrid.NNode)], nsamples)
                rn.shuffle(indsamples)
                




                for i in indsamples:







                    #initreg = ggrid.NodeVec[0].Neighbors
                    inode = ggrid.NodeVec[i].Id
                    leninit = len(ggrid.NodeVec[i].Neighbors)
                    print("target node", i, leninit, "nslices", islice)

                    
                    for inod in ggrid.NodeVec:
                        inod.Amplitude = 0
                        inod.PrevAmplitude = 0
                    ggrid.NodeVec[inode].PrevAmplitude = 1.0
                    ggrid.NodeVec[inode].Amplitude = 1.0
                        

                    #import pdb; pdb.set_trace()
                    if True:
                        print("Heat", ggrid.NNode)
                        arrdecay = ggrid.HeatEq(250, i, True)

                        #wtdheatfractality, wtdheatfractality2 = getweightedfractality(ggrid)
                        #heatstd = np.sqrt(wtdheatfractality2 - wtdheatfractality*wtdheatfractality)


                        
                        arr2 = arrdecay[::2]
                        rng = np.log(1.0+np.arange(len(arr2)))
                        rng2 = rng*rng
                        rng3 = rng*rng*rng
                        front = 0

                        slopedict = {}                
                        if len(rng) - front > 5:

                            for i in range(len(rng)):


                                if len(rng) - i - front < 5:
                                    break
                                #plt.plot(rng[front:], np.log(x2)[front:])
                                #plt.show()
                                nrng = len(rng)

                                heatres = linregress(rng[front:nrng-i],np.log(arr2[front:nrng-i]))
                                heatresparab = linregress(rng2[front:nrng-i],np.log(arr2[front:nrng-i]))
                                heatrestert = linregress(rng3[front:nrng-i],np.log(arr2[front:nrng-i]))
                                

                                #for i in arrdecay:
                                #    print(i)

                                ndata = len(rng[front:nrng-i])

                                print("slopeheat", heatres.slope, len(rng) - front)
                                print("RESHEAT: N: %d node %d nslice %s slope: %f pval %f  slopequadr: %f pvalquadr %f slopetert: %f pvaltert %f " % (ndata, i, str(islice), heatres.slope, heatres.pvalue, heatresparab.slope, heatresparab.pvalue, heatrestert.slope, heatrestert.pvalue))
                                slopedict[len(rng) - front] = {"HEATSLOPE":heatres.slope, "HEATPVAL":heatres.pvalue,  "TOUCHSLOPE":heatres.slope, "TOUCHPVAL":heatres.pvalue} #, "SLOPEDICT":slopedict}
                                #import pdb; pdb.set_trace()
                                heatresdict[(tuple(islice), len(rng)-i)] = {"HEATSLOPE":heatres.slope, "HEATPVAL":heatres.pvalue, }
                        

                    for inod in ggrid.NodeVec:
                        inod.Amplitude = 0
                        inod.PrevAmplitude = 0

                    ggrid.NodeVec[i].Amplitude = 1.0
                    ggrid.NodeVec[i].PrevAmplitude = 1.0

                    
                    
                    if True: 
                        print("Touch", ggrid.NNode)
                        touchgrow = ggrid.Touch(120)

                        #wtdtouchfractality, wtdtouchfractality2 = getweightedfractality(ggrid)
                        #touchstd = np.sqrt(wtdtouchfractality2 - wtdtouchfractality*wtdtouchfractality)


                        
                        arr2 = copy(touchgrow)
                        rng = np.log(1.0+np.arange(len(arr2)))
                        rng2 = rng*rng
                        rng3 = rng*rng*rng
                        front = 0

                        
                        #plt.plot(rng[front:], np.log(x2)[front:])
                        #plt.show()

                        #import pdb; pdb.set_trace()

                        slopedict = {}
                        if len(rng) - front > 5:                    
                            for i in range(len(rng)):
                                nrng = len(rng)
                                touchres = linregress(rng[front:nrng-i],np.log(arr2[front:nrng-i]))
                                touchresparab = linregress(rng2[front:nrng-i],np.log(arr2[front:nrng-i]))
                                touchrestert = linregress(rng3[front:nrng-i],np.log(arr2[front:nrng-i]))

                                #for i in arrdecay:
                                #    print(i)
                                ndata = len(rng[front:nrng-i])

                                print("touchslope", touchres.slope, len(rng) - front)

                                print("RESTOUCH: N: %d node %d nslice %s slope: %f pval %f  slopequadr: %f pvalquadr %f slopetert: %f pvaltert %f" % (ndata, i, str(islice), touchres.slope, touchres.pvalue, touchresparab.slope, touchresparab.pvalue, touchrestert.slope, touchrestert.pvalue))
                                #slopedict[len(rng) - front] = {"TOUCHSLOPE":heatres.slope, "TOUCHPVAL":heatres.pvalue, "TOUCHFRACAVG":wtdheatfractality, "TOUCHFRACSTD":heatstd, "TOUCHSLOPE":heatres.slope, "TOUCHPVAL":heatres.pvalue} #, "SLOPEDICT":slopedict}
                                touchresdict[(tuple(islice), len(rng)-i)] = {"TOUCHSLOPE":touchres.slope, "TOUCHPVAL":touchres.pvalue}
                            
                                #import pdb; pdb.set_trace()

                    

                    
                    for inod in ggrid.NodeVec:
                        inod.Amplitude = 0
                        inod.PrevAmplitude = 0
                    #ggrid.NonZero()
                    #import pdb; pdb.set_trace()
                    ggrid.NNeighborHistogram()
                    #import pdb; pdb.set_trace()
                    #resdict[(i, tuple(islice))] = {"HEATARR": arrdecay, "TOUCHARR": touchgrow, }

            heatdf = pd.DataFrame.from_dict(heatresdict,orient="index")
            touchdf = pd.DataFrame.from_dict(touchresdict,orient="index")
            #import pdb; pdb.set_trace()





"""




python regress_on_etf_spy_lag_correllimited_float_wOFFSET_lessnoetfnospy.py -f /Volumes/PNT3/polygon_data/%/qt_@_%.bz2 --tk F --nonspy XLY --startdt 2022-01-01 --enddt 2023-09-01 --regfile trash_per4_smallnospyfixed_regfileY_@_@@_220101to230901_%sec_lagord6_err_%%_4wndw.txt --freq 5 --lag 5 --threshold -0.15 > /dev/null
./BR_vwapREV.exe --threshold 0.00015 --submitthreshold 0.00015 --TEMPregressionfile /Users/hrvojehrgovcic/quant/polygon/trash_per4_smallnospyfixed_regfileY_F_XLY_220101to230901_5sec_lagord6_err_n0p15_4wndw.txt   --output /Users/hrvojehrgovcic/quant/polygon/trash/temp.txt -y 2023 --startdt 2023-01-04 --enddt 2023-05-26 --exclude /Users/hrvojehrgovcic/quant/data/HolidayExcludeDates.csv --symbol F  --etf XLY --file /Volumes/PNT3/polygon_data/qtfmt_%_@.csv --ticksize 1 --dnin 50 --upin 50 --upout 30450 --probability 2.0 --latency 0.005 --starttimehistory 01:30 --starttime 09:35 --endtime 15:58 --endtimedealinitiation 15:40 --endtimehistory 16:00  --meanreversionfilter +1 --allowgap --dlysum --maxlimitorder 1 --seed 500 --fee 0.14 --rebate 0.4 --append --specificlookback 0 --tracefreq 24hr --fillonlyontrade --mindeals 5 --insidedistance 0 --fastema 100 --slowema 50 --nocancelbecausequantities --vwapinterval 2 > blah_F_2sec_vwap_0p00015_per4fixed2.txt
python simplesr.py blah_F_5sec_vwap_0p00015_per4fixed2.txt
pmset sleepnow

./BR_vwapREV.exe --threshold 0.00015 --submitthreshold 0.00015 --TEMPregressionfile /Users/hrvojehrgovcic/quant/polygon/TRcoeff_nogapremovebadcodes_per4_smallnospyfixed_regfileY_F_XLY_220101to230901_2sec_lagord6_err_n0p15_4wndw.txt --output /Users/hrvojehrgovcic/quant/polygon/trash/temp.txt -y 2023 --startdt 2023-01-04 --enddt 2023-05-26 --exclude /Users/hrvojehrgovcic/quant/data/HolidayExcludeDates.csv --symbol F  --etf XLY --file /Volumes/PNT3/polygon_data/qtfmt_%_@.csv --ticksize 1 --dnin 50 --upin 50 --upout 30450 --probability 2.0 --latency 0.005 --starttimehistory 01:30 --starttime 09:35 --endtime 15:58 --endtimedealinitiation 15:40 --endtimehistory 16:00  --meanreversionfilter +1 --allowgap --dlysum --maxlimitorder 1 --seed 500 --fee 0.14 --rebate 0.4 --append --specificlookback 0 --tracefreq 24hr --fillonlyontrade --mindeals 5 --insidedistance 0 --fastema 100 --slowema 50 --nocancelbecausequantities --vwapinterval 2 > blah_F_2sec_vwap_0p00015_per4filteredcodes.txt
python simplesr.py  blah_F_2sec_vwap_0p00015_per4filteredcodes.txt
pmset sleepnow


def BaseN(i, N=4, padlength=-1):
    ix = copy(i)
    retvec = []
    if i == 0:
        if padlength != -1:
            return [0] * padlength
        else:
            return [0]
    while ix > 0:
        thispart = ix % N
        retvec.append(thispart) # we'll reverse at the end
        ix = ix // N     
    if padlength > 0 and len(retvec) < padlength:
        padding = padlength - len(retvec)
        retvec = retvec + ([0] * padding)
    retvec.reverse()
    return retvec


def laplac(grid, thisvec, NDim, NLEN):
    sumlap = 0
    for i in range(NDim):
        ivec = copy(thisvec)
        #print(ivec)
        ivec[i] = (ivec[i] + 1) % NLEN
        if NDim == 3:
            sumlap += grid[ivec[0], ivec[1], ivec[2]]/2.0/NDim
        elif NDim == 2:
            sumlap += grid[ivec[0], ivec[1]]/2.0/NDim


        ivec = copy(thisvec)
        ivec[i] = ivec[i] - 1
        if NDim == 3:
            sumlap += grid[ivec[0], ivec[1], ivec[2]]/2.0/NDim
        elif NDim == 2:
            sumlap += grid[ivec[0], ivec[1]]/2.0/NDim
    if NDim == 3:
        grid[thisvec[0], thisvec[1], thisvec[2]] += sumlap
    elif NDim == 2:
        grid[thisvec[0], thisvec[1]] += sumlap




print("Don't use this routine -- instead, use PDEGauss() above)

NLEN =  6

NDim = 3
grid = np.zeros(tuple([NLEN]*NDim))
NSteps = 100

if NDim == 3:
    grid[0,0,0] = 1
elif NDim == 2:
    grid[0,0] = 1

for t in range(NSteps):
    if NDim == 3:
        print(grid[0,0,0]) # , np.sum(grid))
    elif NDim == 2:
        print(grid[0,0])
    #import pdb; pdb.set_trace()
    for j in range(NLEN**NDim):
        thisvec = BaseN(j, NLEN, NDim)
        if (np.sum(thisvec) + t) % 2 == 0:
            continue
        #if tuple(thisvec) == (5,0,0) or tuple(thisvec) == (0,5,0) or tuple(thisvec) == (0,0,5):
        #    import pdb; pdb.set_trace()
        #import pdb; pdb.set_trace()
        laplac(grid, thisvec, NDim, NLEN)
    
    # cleanup
    
    for j in range(NLEN**NDim):
        thisvec = BaseN(j, NLEN, NDim)
        if (np.sum(thisvec) + t) % 2 == 1:
            continue
        if NDim == 3:
            grid[thisvec[0], thisvec[1], thisvec[2]] = 0
        elif NDim == 2:
            grid[thisvec[0], thisvec[1]] = 0










# python  ~/quant/latticegas_nodes.py  -t 10 -p 1.0

# this gives us a 3-d "grid" with about 100 * 100 * 100 nodes in the space
python  ~/quant/latticegas_nodes.py  -t 48 -p 0.1  -d 3 --seed 24

# this gives us a 3-d "grid" with about 20 * 20 * 20 nodes in the space
python  ~/quant/latticegas_nodes.py  -t 30 -p 0.1  -d 3 --seed 24

HOW TO INTEGRATE D-DIMENAIONAL GAUSSIAN RADIALLY

https://math.stackexchange.com/questions/3880777/n-dimensional-integral-of-radial-function
https://math.stackexchange.com/questions/1482747/integral-in-n-dimensional-spherical-coordinates?rq=1
https://en.wikipedia.org/wiki/List_of_integrals_of_Gaussian_functions#Indefinite_integrals
https://upload.wikimedia.org/wikipedia/commons/a/a2/Cumulative_function_n_dimensional_Gaussians_12.2013.pdf
https://deepblue.lib.umich.edu/bitstream/handle/2027.42/117496/jdw140914.pdf?sequence=1#:~:text=Exploiting%20the%20rotational%20invariance%2C%20and%20using%20identities,on%20the%20number%20of%20dimensions%20are%20discussed.
https://www.osti.gov/servlets/purl/1991086#:~:text=4%20Relation%20of%20Gaussian%20to,i=1  [NOT GOOD -- NO DEFINITE]

https://upload.wikimedia.org/wikipedia/commons/a/a2/Cumulative_function_n_dimensional_Gaussians_12.2013.pdf

See https://davidegerosa.com/nsphere/  or https://math.stackexchange.com/questions/3880777/n-dimensional-integral-of-radial-function
f(r) = exp(-r*r/2)

IntegFrom0toC of f(r) = Integ(0,C) f(r) S(r) dr

where S = Integ(S_N-a(r)) dA  
S is the surface of an N-1 dimensional sphere of radius r and therefore explicitly depends on r. 
The simplest closed form is 2 * pi**(n/2) / Gamma(n/2) * r**(n-1)
So in 2 dimensions, it is just the circumference of a circle, or 2*pi*r and we get the old volumes of revolution result: https://math.libretexts.org/Bookshelves/Calculus/Calculus_(OpenStax)/06%3A_Applications_of_Integration/6.03%3A_Volumes_of_Revolution_-_Cylindrical_Shells
and in 3d it is just the area of a sphere, or 4*pi*r*r

So by collecting one of the r's in r**(n-1) to sit with the Gaussian
you are integrating exp(-4*4/2)*r*dr * r**(n-2) 

which can be solved (via integration by parts) semi-explicitly (i.e. in terms of an erfc(r) fumction) plus a bunch of
definite integrals from 0 to C. In the literature, there is usually a change of variables so that C is 1 but the end result is the same.

You can write a python routine to do the integration-by-parts for any N



https://web.stanford.edu/class/math220b/handouts/heateqn.pdf  (page 24)



%run  /Users/hrvojehrgovcic/quant/latticegas_cubenodes2.py  -t 20  -p 1.0 --dim 2



%run  /Users/hrvojehrgovcic/quant/latticegas_cubenodes2.py  -t 8  -p 1.0 --dim 3

%run  /Users/hrvojehrgovcic/quant/latticegas_cubenodes2.py  -t 2000  -p 1.0 --dim 3

%run  /Users/hrvojehrgovcic/quant/latticegas_cubenodes3.py  -t 20  -p 1.0 --dim 3 --length 100 --expand 20

%run  /Users/hrvojehrgovcic/quant/latticegas_cubenodes4.py  -t 20  -p 1.0 --dim 2 --length 100 --expand 20

%run  /Users/hrvojehrgovcic/quant/latticegas_cubenodes6.py  -t 1000  -p 1.0 --dim 2 --length 200000 --expand 20000
"""