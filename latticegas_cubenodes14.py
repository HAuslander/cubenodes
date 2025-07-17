#!/cygdrive/c/Python27/python


global FreshProb
global MaxNodes


print("""

The goal is to create some weird quasi-Sierpinski structure in n dimensions
so that we're left with a family of random N-dimensional graphs for which 
the demand that every edge can be traversed in one time step
plus maybe dynamics (or maybe even just the radius) automatically embeds
the graph in a Euclidean space.

I.e. this is about how to create graphs where the Euclideanness of the 
underlying space is inferred.

Originally, had sought to create a suitably random N-dimensional grid as follows:

1) start with a Dim+1 simplex (i.e. you're creating a "simply connecte" genus-0 system of dimensionality Dim)
2) randomly choose to replace each node in this simplex with  a simplex (so in 2-d, you take a node with 3 neighbors and replace it
   with a triangle of nodes, with the new nodes connected to each other (and each one is connected to a different neighbor of the old node)

3) by re-using nodes (i.e. one of the new nodes that replaces an old node is just the old node, but with a different set of neighbors)
   we can avoid having to tabulate all the nodes in a dictionary, and instead use a list. Maybe that will help with speed and storage in a large system

The above approach failed miserably. The "large expanses" formed by the faces of the original simplex never got filled in, and
the dimensionality was something lower than that of the simplex.

In the updated version, we start with cubes of length 8 in either dimension (toroidally connected)
We will generally divide each dimension with TWO cuts instead of just bisecting with one
so as to preserve the parity of each node, but the number of cuts can possibly be increased further.
(It may be that for larger slice-per-cube divisions the resultants paths can be more problematic, but if so,
restricting ourselves to cuts of 1 or 2 is not particularly problematic.)

ALSO, the main difference between this and the earlier version is we will allow adjacent (i.e. sharing a face or edge) cubes to be divided, which means
matching up previously created mid-nodes (so a cube can only hava fractality 1 greater than surrounding cubes). In the previous version once a cube
was divided, no adjacent cubes could be divided (except the "diagonal" cubes that share a single point).


In general, it is best to prefer dividing cubes that are adjacent to some cube with higher fractality, than just
any cube; otherwise, the finer cubes will eventually predominate the space and the lower fractality lingerers will
keep getting left behind.

    Therefore, the are 2 probabilities to consider; a really small one and a higher one; the higher one will
    only act on nodes where two or more magnifications intersect (which only happens if one of its axes is
    an edge of a division, so that the nodes neighbors have less than 2D neighbors of their own)


    What to do next:
    1) preferentially divide cubes whose coord_paths have some "long" paths,
        as opposed to those where the paths are all of length two (so that
        all edges of the cube connect nearest neighbors)
    2) create a "sphere" of a given number of steps (or steps of length nslice+2 whose
        intervening nodes are an edge of a cube) and whenever a node has no
        available cubes for division (and maybe also exclude the ones whose cubes have
        corners connected by nearest neighbors) find the nearest node that does have
        some mixed fractality (i.e. cubes have some edges that have divided paths) and
        do that one -- maybe that will eventualy lead to the stuck node with no possible
        cube divisions to get unstuck



   Jul-2025; have encoutered situations in 3-d that require an overhaul of what it meanns for a node to be "packed".
   In 2-d, the number of neighbors is sufficient to distinguish all 'edge' nodes. In 3-d, it's true that edges of newly
   divided cubes can have 4 neighbors, and faces of newly divided cubes can have 5 neighbors, but consider 2 rubik's cubes
   arranged as 2 steps on a staircase, and therefore having one edge in common. The nodes along that edge all have 6 neighbors,
   even though they are both "newly divided". 

   So the updated criterion is "new-growth" (or newly divided) vs. old-growth. A node is part of an edge until all 2**D adjoining
   cubes (or regions) have an opposite cube. Eventually, all nodes will have 2**D "opposing nodes" (i.e. nodes that are the "opposite"
   node of some hypercube), and until that time they are "interfaces" or generalized "edges". Hopefully, some combination of NOppositeNode
   AND ALSO NNeighbors will be enough to specify. For now, each division will introduce a new set of opposing nodes.

   (Compare the staircase Rubik's cube, with 2 edges of a cube coinciding, with the more typical case of a node created along a cube's face. 
   Both those types of nodes have 4 opposing nodes per newly-created-node, but in the former case, there are 6 neighbors per node, while 
   in the latter there are 5; also, the grouping of the new nodes is different under each scenario, with two seperate additions of 2 
   opposing nodes in the former case, and one addition of 4 nodes in the latter case. We'll try and keep track of both, since in higher 
   dimensions that may be relevant as well, but in any number of dimensions, the notion that a newly created node is "new-growth" until
   such a time that it attains 2**D opposite-node regions is a notion that will generalize across dimensions)






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
 

 # https://stackoverflow.com/questions/49139727/creating-a-rational-class-and-operator-overloading-in-python
class RationalNumber: # to be used in coordinates
    numerator=1
    denominator=1
    lcm=1
    def __init__(self,numerator,denominator=1):
        if(denominator<=0):
            print('Denominator can not be <=0')
        else:
            gcd=1
            divisorsOfNum1 = []
            divisorsOfNum2 = []
            gotlcm=False
            for i in range(1, int(numerator)+1):
                if numerator % i == 0:
                    divisorsOfNum1.append(i)
            for i in range(1, int(denominator)+1):
                if denominator % i == 0:
                    divisorsOfNum2.append(i)
            for i in range(0, len(divisorsOfNum1)):
                for j in range(0, len(divisorsOfNum2)):
                    if divisorsOfNum1[i] == divisorsOfNum2[j]:
                        if(gotlcm==False):
                            gotlcm=True
                            self.lcm=divisorsOfNum1[i]
                        gcd = divisorsOfNum1[i]
                        continue
            self.numerator=numerator/gcd
            self.denominator=denominator/gcd
            #print(self.numerator,self.denominator)
    def __add__(self,other):
       numeratr=self.numerator*other.denominator + other.numerator*self.denominator
       denominatr=self.denominator*other.denominator
       return RationalNumber(numeratr,denominatr)
    def __sub__(self,other):
       numeratr=self.numerator*other.denominator - other.numerator*self.denominator
       denominatr=self.denominator*other.denominator
       return RationalNumber(numeratr,denominatr)
    def __truediv__(self,other):
       numeratr=self.numerator*other.denominator  
       denominatr=other.numerator*self.denominator
       return RationalNumber(numeratr,denominatr)
    def __mul__(self,other):
        numeratr=self.numerator*other.numerator  
        denominatr=other.denominator*self.denominator
        return RationalNumber(numeratr,denominatr)
    def __ne__(self,other):
        if self.numerator!=other.numerator  and other.denominator!=self.denominator:
            return True
        else:  
            return False
    def __eq__(self,other):
        if self.numerator==other.numerator  and other.denominator==self.denominator:
            return True
        else:  
            return False
    def __gt__(self,other):
        if (self.numerator/self.denominator)>(other.numerator/other.denominator):
            return True
        else:  
            return False
    def __lt__(self,other):
        if (self.numerator/self.denominator)<(other.numerator/other.denominator):
            return True
        else:  
            return False
    def __ge__(self,other):
        if (self.numerator/self.denominator)>=(other.numerator/other.denominator):
            return True
        else:  
            return False
    def __le__(self,other):
        if (self.numerator/self.denominator)<=(other.numerator/other.denominator):
            return True
        else:  
            return False
"""
r1=RationalNumber(48,36)
r2=RationalNumber(48,36)
print(r1+r2)
"""



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

  p.add_option("--expand", "--expansion", "-t", default=100,
                  action="store", dest='expansionruns', type='int',
                  help="how many time steps do we expand, using the --prob arg on all the nodes available (understanding that the available nodes may increase with time)")


  p.add_option("--length", default=-1,
                  action="store", dest='length', type='int',
                  help="if > 0 the target number of nodes will be at or slightly greater than this parameter**Dim")



  p.add_option("--maxnodes", "--maxnode", default=-1000000,
                  action="store", dest='maxnodes', type='int',
                  help="the division loop is exited when maxnodes exceed")


  p.add_option("--post",  default=False,
                  action="store_true", dest='PostAnalysis', 
                  help="the division loop is exited when maxnodes exceed")


  
                
  return p.parse_args()                



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



class listiterator_t:
    """ Deprecated in favor of dictionaryiterator_t() which for in the present file makes the jumble of list
    indices harder to parse than if we keep  keys and list_indices a separate concept
    """
    def __init__(self, baselist):
        if np.min(baselist) < 1:
            print("error in listiterator_t() -- components of base list must be > 0")
            import pdb; pdb.set_trace()
        self.BaseList = baselist
        self.N = len(baselist)
        self.Period = np.prod(baselist)
        self.Current = [ 0 for i in range(len(baselist))]
        self.CurrentN = 0
    
    def Increment(self,):
        bDone = False
        icol = self.N - 1
        while not(bDone):
            if self.Current[icol] == (self.BaseList[icol] - 1):
                self.Current[icol] = 0
                icol -= 1
                if icol < 0:
                    bDone = True
                
            else:
                bDone = True
                self.Current[icol] += 1

        self.CurrentN = self.GetCurrentN()

    def GetCurrentN(self,):
        imult = 1
        xsum = 0
        for icol in range(len(self.BaseList)-1,-1,-1):
            xsum += self.Current[icol] * imult
            imult *= self.BaseList[icol]
        return xsum

    

class dictionaryiterator_t:
    """ iterates over a dictionary of lists, so that if keys A, B, and C have lists of 3, 4, 5 elemennts each,
    you can iterate over every distinct combination
    """

    def __init__(self, basedict):
        if np.min(list(basedict.values())) < 1:
            print("error in dictionaryiterator_t() -- components of base list must be > 0")
            import pdb; pdb.set_trace()
        self.BaseDict = basedict
        self.KeyList = list(basedict.keys())
        self.KeyList.sort()
        self.N = len(basedict)
        self.Period = np.prod(list(basedict.values()))
        self.Current = copy(basedict)
        for key,val in self.Current.items():
            self.Current[key] = 0
        self.CurrentN = 0
    
    def Increment(self,):
        bDone = False
        ikey = self.N - 1
        while not(bDone):
            key = self.KeyList[ikey]
            if self.Current[key] == (self.BaseDict[key] - 1):
                self.Current[key] = 0
                ikey -= 1
                if ikey < 0:
                    bDone = True
                
            else:
                bDone = True
                self.Current[key] += 1

        self.CurrentN = self.GetCurrentN()

    def GetCurrentN(self,):
        imult = 1
        xsum = 0
        for ikey in range(len(self.BaseDict)-1,-1,-1):
            key = self.KeyList[ikey]
            xsum += self.Current[key] * imult
            imult *= self.BaseDict[key]
        return xsum

    


"""

lit = listiterator_t([5,7,4])
#lit = listiterator_t([1,1,1])
print("Period", lit.Period)

print(lit.BaseList, lit.Period)
import pdb; pdb.set_trace()
for i  in range(lit.Period):
    print(lit.Current, lit.CurrentN)
    lit.Increment()

"""

class binomialcoefficienttree_t():
    """
    a list of dictionaries whose keys and entries share similarities (as far as their multiplicity) with the D-th row of the pascal's 
    triangles (and respectively correspond to lines, faces, cubes... and other subcubes)
    Each key will be the binary representation of a number in [0...2**D]
    An entry whose key is a binary vector that has m ones will have a value that is a sub-dictionary with m keys, 
    corresponding to the m vectors (each of which correspond
    to entries in the mail directory that have already been filled out) that are obtained from the key by
    zeroing out one coefficient. 

    If m=2, then the implied structure (barring abnormalities) a face (i.e. an intersection of 2 line-segments), 
    if m=3, then it is a cube (i.e. an intersection of 3 faces), if m=4 then it is a tesserect (an intersection of 4 cubes).
    The values of these m subdirectory entries will be the m path indices from pathdump whose endpoints are the same node.
    If there are MULTIPLE entries (i.e. m paths converge in point A while m others converge in point B) then for now,
    it means the class needs to be branched into two seperate instances (which will be contained in a list). The Pascal's
    triangle reference means, in, say, the 2-dim case, there is one origin, with 2 nodes creating one face, and a far node.
    In the 3-dim case, there is one origin, having 3 axes, yielding 3 faces, yielding one far node.
    In the 4-d case, we'd have one origin, 4 axes yielding 6 faces whose intersections result in 4 cubes, yielding one far node.
    
    Alternatively, we will modify the notion of "well-formed cube" in such a way that the very existence of such cases 
    means the cubes that we are attempting to construct are flawed and the attempt should be aborted,
    in the same way that two distinct edges having the same endpoints is an indication that both should be removed. (This makes
    sense, since the "intersection paths" will eventually be edges of the created cubes, so that if we impose a non-muliplicity restriction
    on the starting axes that are themselves prospective edges of the cube, the same criterion should be mandated for all other edges.)
    If, on the other hand, there is NO m-fold intersection, that definitely means this particular dictionary is a dead-end (and means that
    one of the edges of the prospective cube cannot be found) and needs to be aborted.

    NOTE: The above structure is not used (instead the values of the directory will be a 2-tuple consisiting of a subdirectory and an endpoint,
    but we will leave the class here in order to explain how the relevant object will 



    """

    def __init__(self, pathdump):
        self.d = {}
        self.Pathdump = pathdump

    def EndPt(self, binvec, pathdump):
        if np.sum(binvec) == 0:
            return None # for reasons of completeness only
        else:
            if len(self.d[binvec]) == 0:
                return None
            for key, val in self.d[binvec].items():
                return self.Pathdump[val][-1]

class node_t():
    def __init__(self, id):
        self.Id = id
        self.Neighbors = [] # here, every element of the list is of type neighbor_t
        self.Amplitude = 0
        self.PrevAmplitude = 0
        self.bDefunct = False
        self.Parity = 0 # this is only used as a temporary scratch pad, though if nslices is even, the parity of a point during the expansion will remain constant
        self.Coords = None # not used for anything but post analysis and forensic analysis on nodes where the dimension parameters are not close to self.NDim
        self.RatCoords = None # rational-number versions of coords that should allow for more accurace; as with Coords, it is only used for visualization and debugging (for now. program continues to use plain floating-point coords); since it slows the program to a crawl for large numerators/denominators, it has been disabled
        self.LastDivide = 0 # not needed; it may prove useful in efforts "even out" the divisions so that the resultant lattice does not have large gaps

        self.NRegions = 0
        #self.RegionGrowth = [] # hopefully this will not need to be used, at least for 3D

        # region count gives the number of "colors" converging on a node. At the start of the cube creation, every node has 2**D cubes.
        # In, say, 3-d, When one of those cubes gets subdivided, the "middle" dube edges have 4 neighbors and 5 regions (i.e. the 4 regions that
        # adjoined before there was a node there, but one of those 4 regions is subdivided so 3 "long regions, and 2 divided ones"
        # the faces, have 5 neighbors and also 5 regions (one "planar" region on one side of face, plus 4 divided ones on the other)
        # in two dimensions, neighbor count can distinguish between the relevant sub-cubes (i.e. edges), but in 3 dimensions, we can
        # have 2 cubes share only a single edge (so that the cross section is two squares on a diagonal), and in that case there are
        # 6 neighbrs, and also 6 regions (2 "long ones", and two sets of two short ones). That is distinct from a fully-divided cube
        # with has 6 neighbors and 8 regions. If we only used neighbor count, then this type of intervening node would be considered
        # a "fully packed" node, so that if we were doing a search for paths along that edge, they would not span the entire edge, but
        # only go to the nearest neighbor. But if we can recognize that this is indeed an edge, we won't make that mistake.

        # It remains to be seen whether region count plus neighbor count is enough to uniquely specify the different possible "intervening"
        # nodes in four and higher dimeensions, and we may need additional specifiers, but the overall itent is clear: a well-formed cube
        # has edges consisting either of nearest neighbors (connectiong two points with maximal 2**D inntersecting regions), or else
        # an edge consisisting of 2+nslices nodes having two endpoints with maximal 2**D intersecting regions, and nslices intervening
        # nodes that have the exact same type of edge.



        

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
        self.NodeVec = [] # decided to use a vector instead of dictionary to hold the nodes
        self.NNode = 0 

        self.NodeHistory = {} # all cubes -- the values will be of class

        self.TorLen = 8 # 32
 

        # NOTE: cubes of a given fractality can only abut other smaller/larger cubes whose fractality differes by one; otherwise, trying to match up adjacent edges/faces/etc. becomes intractable                      

        self.NodeCubes = {} # this is a dictionary of the 2**Dim cubes that touch a given node
        self.MagnificationRingLen = 100 # 4 EVENTUALLY, ONCE IT CHECKS OUT, SET IT TO 4
        
        self.NDivisions = 0

        self.HistScale = 1 # can be some higher integer -- see below
        histshape = tuple([self.TorLen * self.HistScale] * self.NDim)
        self.Hist = np.zeros(histshape).astype("int")
        self.NodesLastHistogramUpdate = 0

        self.OnesCoord = tuple([1 for i in range(self.NDim)]) # this gets used frequently, so only calc it once

        self.NSlices = 2 # can be any positive integer; even integers will retain even/odd
        # parity of ndoes

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

        self.NeighborDistanceCount = {}

        self.bStop = False

        #self.AllDirections = list(np.arange(1,self.NDim+1)) +  list(np.arange(-1,-self.NDim-1,-1))
        sortbyones = []  
        for j in range(2 ** self.NDim): # note that here len(nbrs) is 2*D
            #thisnode = copy(Id) 
            binstr = bin(j)[2:]
            while len(binstr) < self.NDim:
                binstr = '0' + binstr
            lenbinstr = len(binstr)
            binvec = [int(chr) for chr in binstr]
            sortbyones.append((np.sum(binvec), tuple(binvec)))           
        sortbyones.sort()
        self.SortByOnes = sortbyones

    def MaxDist(self, inode):
        d = 0
        for inbr in self.NodeVec[inode].Neighbors:
            nbrd = self.l2node(inode, inbr)
            d = np.max([d, nbrd])
        return d

    def UpdateNeighborDistanceCount(self):
        #if self.NNode > 1000:
        #    import pdb; pdb.set_trace()

        self.NeighborDistanceCount = {}
        for iinode, inode in enumerate(self.NodeVec):
            x = np.array(inode.Coords)            
            for inbr in inode.Neighbors:
                thisdist = np.linalg.norm((np.array(self.NodeVec[inbr].Coords)) - x)
                if thisdist > 4:
                    thisdist -= 8
                    thisdist = np.abs(thisdist)
                thisinvdist = np.round(1.0/thisdist)
                if not(thisinvdist in self.NeighborDistanceCount):
                    self.NeighborDistanceCount[thisinvdist] = 1
                else:
                    self.NeighborDistanceCount[thisinvdist] += 1

        print("Hist update for nnode count", self.NNode)
        for key, val in sorted(self.NeighborDistanceCount.items()):
            print("d", key, val)
            if key == 1.0 and val  == 30 and self.NNode > 25000:
                self.bStop = True
    
    def FindNeighborDistances(self, InvDist):
        myretval = []
        for iinode, inode in enumerate(self.NodeVec):
            x = np.array(inode.Coords)            
            for inbr in inode.Neighbors:
                thisdist = np.linalg.norm((np.array(self.NodeVec[inbr].Coords)) - x)
                if thisdist > 4:
                    thisdist -= 8
                    thisdist = np.abs(thisdist)
                thisinvdist = np.round(1.0/thisdist)
                if thisinvdist == InvDist:
                    tup = [inode.Id, self.NodeVec[inbr].Id]
                    tup.sort()
                    tup = tuple(tup)
                    if not(tup in myretval):
                        myretval.append(tup)
        return myretval




    def SumAmp(self,bPrint=True):
        myret = np.sum([inode.Amplitude for inode in self.NodeVec])
        if bPrint:
            print("Sumamp ", myret)
        return myret


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
                nbrs = inode.Neighbors
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
                        nbrs = [jnode for jnode in inode.Neighbors]
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
                nbrs =copy(inode.Neighbors)
                for inbr in nbrs:
                    denomin = 1.0/len(self.NodeVec[inbr].Neighbors)
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
                nbrs = copy(inode.Neighbors)
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

    def bCheckCoordLinearity(self, A, B):
        """
        This is used for debugging only, to ensure that two points are related by a translation along one of the cubes axes
        """
        bSomeAxisIsSame = False
        coordsA = self.NodeVec[A].Coords
        coordsB = self.NodeVec[B].Coords
        for i in range(self.NDim):
            if coordsA[i] == coordsB[i]:
                bSomeAxisIsSame = True
                break
        
        return bSomeAxisIsSame
            

                    




    def MismatchedNeighbors(self, Id, DesiredCount, dontinclude=None):
        myretlist = []
        for nbrid in self.NodeVec[Id].Neighbors:
            if len(self.NodeVec[nbrid].Neighbors) != DesiredCount:
                if nbrid == dontinclude:
                    continue
                myretlist.append(nbrid)
        return myretlist

    def AcrossADividedEdge(self, Id, nbrId):
        # if we go along an edge, then we search for two end nodes with 2D neighbors; 
        # if we run THROUGH  a cube, then the starting and  endpoints are on an edge or face 
        # and  they  are the ones with lses than 2D neighbors, while the  intermediate steps have more

        # NOTE: when we're in a channel between a fractality  on one side and a lower one
        # on the other, the routine can return a weird knight's move, which is only a 
        # problem when the division would cause that channel to disappear (and therefore
        # lead to two  adjoining regions where the  ractality differs by two which is
        # verboten), so  it is best to weed out such scenarios before we need to make
        # use  of this  routine
        curr = copy(nbrId)
        prev = copy(Id)
        iloop = 0

        myretlist = [prev, curr]

        EndPointNeighborCount = len(self.NodeVec[Id].Neighbors)

        while len(self.NodeVec[curr].Neighbors) != EndPointNeighborCount:
            theseunpackedneighbors = self.MismatchedNeighbors(curr, EndPointNeighborCount, prev)

            dellist = []
            for iunp in theseunpackedneighbors:
                if self.bSimpleIsItAFace(curr, prev, iunp):
                    dellist.append(iunp)
            for idel in  dellist:
                theseunpackedneighbors.remove(idel)
                    
            if len(theseunpackedneighbors) > 1:
                print("more unpacked neighbors than expected in AcrossADividedEdge")
                import pdb; pdb.set_trace()
            if len(theseunpackedneighbors) == 0:
                break
            prev = copy(curr)
            curr = theseunpackedneighbors[0]
            myretlist.append(curr)
            iloop += 1
            if iloop > 10:
                break

        prevnbrs = self.NodeVec[prev].Neighbors
        currnbrs = self.NodeVec[curr].Neighbors

        allowablenext = []
        for icurr in currnbrs:
            if icurr == prev:
                continue
            Anbrs = prevnbrs
            Bnbrs = self.NodeVec[icurr].Neighbors
            intersected =  list(set(Anbrs).intersection(Bnbrs))  
            intersected.remove(curr)
            if len(intersected) == 0:
                allowablenext.append(icurr)
        if len(allowablenext) != 1:
            return myretlist
            #print("what is going on here in AcrossADividedEdge -- there should only be one noncollinear neighbor")
            #import pdb; pdb.set_trace()

        #if (Id, nbrId) == (1034, 1031):
        #    print("Y")
        #    import pdb; pdb.set_trace()


        return myretlist + [allowablenext[0]]


        print("AcrossADividedEdge has failed to go around the edge")
        import pdb; pdb.set_trace()









                    
    
    def FarNeighbor(self, Id, nbrId, PossibleTerminalNodes=None): 
        # PossibleTerminalNodes is just NodeCoord.keys() --  i.e. the outer ndoes of the cube being divided

     

        """
        Expected usage will be the nodes of a cube that is to be divided
        Will return Far Neighbor (i.e. other node along the edge of some cute)
        and all the intervening nodes (if any) from some previous split

        NOTE: This fails (i.e. finds a weird kinght's move  or other kinked far-neighbor path) for nodes of divided cubes that are one  step away from a cube where the fractality is two orders lower

        NOTE THAT if the a step of this path reaches a node with less than 2d neighbors,
        we have to do a custom  is-ot-a-face routine between that node and the two preceding ones

        If we do not do this, we have situations in which a divided region (say its fractality is order 1) has a further-divided region (fractalith 2) that is ALMOST next to the boundary
        (so that it is as close as can be to a region where the cubes are even larger, of fractality 0). If the first step hits that boundary edge between fractality 1 and 0,
        then the far-neighbor routine will 
        """



        if len(self.NodeVec[Id].Neighbors) == len(self.NodeVec[nbrId].Neighbors):
            return [Id, nbrId], True

        #if self.NNode == 1392 and Id == 1044 and nbrId == 1047:
        #    import pdb; pdb.set_trace()

        #if self.NNode == 1392 and Id == 1044 and nbrId == 1047:
        #    import pdb; pdb.set_trace()

        #if (Id, nbrId) == (33, 1030) and self.NNode == 1072:
        #    pass
            #import pdb; pdb.set_trace()
            #A = AcrossADividedEdge(Id, nbrId)


        return self.AcrossADividedEdge(Id, nbrId), True

        NLoop = 0


        while NLoop < 100: # this is only to prevent some infinite loop

            current = myretval[-1]
            prev = myretval[-2]

            prevneighbors = self.NodeVec[prev].Neighbors
            prevfarnbrs =  []
            for inbr in  prevneighbors:
                if not(self.bIsCorner(inbr)):
                    prevfarnbrs.append(self.AcrossADividedEdge(prev, inbr))
            
            currneighbors = self.NodeVec[current].Neighbors
            currfarnbrs =  []
            for inbr in  currneighbors:
                if not(self.bIsCorner(inbr)):
                    currfarnbrs.append(self.AcrossADividedEdge(current, inbr))
            

            
            PossibleNextSteps = []
            best_nbrs_inbr = []
            best_nbrs_farnbrs = []

            for inbr in self.NodeVec[current].Neighbors:
                if inbr == prev:
                    continue # we want the neigbor in the opposite direction of prev

                
                nbrs_inbr = self.NodeVec[inbr].Neighbors
                nbrs_farnbrs = []
                for jnbr  in nbrs_inbr:
                    if not(self.bIsCorner(inbr)):
                        nbrs_farnbrs.append(self.AcrossADividedEdge(inbr, jnbr))       


                Anbrs = prevneighbors + prevfarnbrs
                Bnbrs = nbrs_inbr + nbrs_farnbrs
                intersected =  list(set(Anbrs).intersection(Bnbrs))        

                if current in intersected:
                    intersected.remove(current)     
                
                if len(intersected) == 0:
                    PossibleNextSteps.append(inbr)
                    best_nbrs_inbr = copy(nbrs_inbr)
                    best_nbrs_farnbrs = copy(nbrs_farnbrs)
                
                else:
                    return myretval, True



   
    def bAreTwoAxesBadlyOriented(self, Id, nbr2subset): 
        # PossibleTerminalNodes is just NodeCoord.keys() --  i.e. the outer ndoes of the cube being divided
            
            if not(self.bIsCorner(Id)):
                print("why is bAreTwoAxesBadlyOriented being called for a node that does not have 2D neighbors?")
                import pdb; pdb.set_trace()
                return True
            
            relevantA = nbr2subset[0]
            relevantB = nbr2subset[1]

            if self.NodeVec[relevantA].NRegions != 2 ** self.NDim:
                relevantA = self.AcrossADividedEdge(Id, relevantA)[-1]
            if self.NodeVec[relevantB].NRegions != 2 ** self.NDim:
                relevantB = self.AcrossADividedEdge(Id, relevantB)[-1]

            relnbrsA = []
            for inbr in self.NodeVec[relevantA].Neighbors:
                if len(self.NodeVec[inbr].Neighbors) == 2 * self.NDim:
                    relnbrsA.append(inbr)
                else:                
                    relnbrsA.append( self.AcrossADividedEdge(relevantA, inbr)[-1] )

            relnbrsB = []
            for inbr in self.NodeVec[relevantB].Neighbors:
                if len(self.NodeVec[inbr].Neighbors) == 2 * self.NDim:
                    relnbrsB.append(inbr)
                else:
                    relnbrsB.append( self.AcrossADividedEdge(relevantB, inbr)[-1] )

            intersected =  list(set(relnbrsA).intersection(relnbrsB))  
            if not(Id in intersected):
                import pdb; pdb.set_trace()
            intersected.remove(Id)
            return len(intersected) == 0


    def BruteDivideCube(self, CoordNode, NodeCoord, coord_path):


        """
        sort the desired new coords by 
        A) number of intervening components (i.e. components not equal to 0 or cubelen)
        B) index of intervening node

        Having done this, any desired matchup between about-to-be-created and previously-created
        can be found by looking for the intersection of nodes already matched up.

        E.g. any node immediately next to two edges (e.g. the (1,1‚  node) is the intesection of
        two edge nodes, and then the rest (i.e. (2,1), (3,1)... and (1,2), (1,3), (1,4)...) can be 
        found the same way -- just decrement two comonents; there's gotta be one already created (if
        the adjacent node was divided)
        """
          
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
            #try:
            #    thisaxis_rational = list( self.NodeVec[node0].RatCoords )
            #except:
            #    import pdb; pdb.set_trace()
            coord0 = copy(thisaxis) # this is the anchor vector
            for i in range(self.NDim):
                coord_thisdim = np.array( self.NodeVec[axisedges[i]].Coords )
                #coord_thisdim_rational = [ self.NodeVec[axisedges[i]].RatCoords[j] for j 
                try:
                    if np.max(np.abs(coord_thisdim - coord0)) > 1:
                        for iicomp, icomp in enumerate(coord0):
                            if np.abs(coord_thisdim[iicomp] - coord0[iicomp]) > 1:
                                if coord_thisdim[iicomp] > coord0[iicomp]:
                                    coord0[iicomp] += self.TorLen
                                elif coord_thisdim[iicomp] < coord0[iicomp]:
                                    coord_thisdim[iicomp] += self.TorLen
                        
                    thisaxis = thisaxis + (coord_thisdim - coord0) * thistuple[i]/dimlen
                    #for j in range(len(thisaxis_rational)):
                    #    thisaxis_rational[j] = thisaxis_rational[j] + RationalNumber(((coord_thisdim[j] - coord0[j]) *  thistuple[i]), dimlen)
                except:
                    print("BAD")
                    import pdb; pdb.set_trace()

            for iicomp, icomp in enumerate(thisaxis):
                if icomp < 0:
                    thisaxis[iicomp] += self.TorLen
                    #thisaxis_rational[iicomp] = thisaxis_rational[iicomp] + self.TorLen
                thisaxis[iicomp] %= self.TorLen

                #num = thisaxis_rational[iicomp].numerator
                #den = thisaxis_rational[iicomp].denominator
                #num %=  (self.TorLen * den)
                #thisaxis_rational[iicomp] = RationalNumber(num, den)
            #return (tuple(thisaxis), tuple(thisaxis_rational))
            return (tuple(thisaxis), None)



        def allintervening(somevec, cubelen):
            myretval = []
            for i in list(somevec):
                if i != 0 and i != cubelen:
                    myretval.append(i)
            return tuple(myretval)
        
        def CoordRun(scaledsubcoord,scaledcoord, cubelen):
            retcoords = []
            for i in range(len(scaledsubcoord)):
                if scaledsubcoord[i] == 0 and scaledcoord[i] == cubelen:
                    for j in range(cubelen+1):
                        thiscoord = list(scaledsubcoord)
                        thiscoord[i] = j
                        thiscoord = tuple(thiscoord)
                        retcoords.append(thiscoord)
            return retcoords

        def ToBeIntersected(icoord,  cubelen):
            """
            let ninter be the number of intervening components of icoord

            return ninter icoords (each  of which has differes by one decremented component from icoord whose inersection is the node corresponding to icoord.
            For example, if icoord is (1,1), then the returned vectors are (1,0) and (0,1), and the non-trivial intersection of their neighbors will be the node corresponding
            to (1,1)

            If coord is (1,1,1), then returned values are (1,1,0),(1,0,1) and (0,1,1)
            If coord is (1,2,1), then returned values are (1,2,0),(1,1,1) and (0,2,1) -- note that because of the sorting of the nodes to be matched, each of the rturned
            nodes will have already been matched by the time we need to d the matching
            
            """

            myret = []

            for ii,i in enumerate(icoord):
                if i != 0 and i != cubelen:
                    thiscoord = list(icoord)
                    thiscoord[ii] -= 1
                    myret.append(tuple(thiscoord))
            return myret
            
        nslices = self.NSlices

        IdCoord = tuple([0] * self.NDim)

        if CoordNode is None:
            return 
        Id = CoordNode[IdCoord]

        cubelen = nslices + 1
        nodelen = cubelen + 1

        newcoordnode = {}
        newnodecoord = {}

        bHistory = True



        
        for key,val in CoordNode.items():
            newkey = tuple([cubelen * icomp for icomp in key])
            newcoordnode[newkey] = val
            newnodecoord[val] = newkey
        

        


        axisedges = [] # the first node of the edge is just (0,0,...) -- i.e. Id -- and is understood
        for i in range(self.NDim):
            thisaxis = np.zeros((self.NDim,)).astype("int")
            thisaxis[i] = 1
            #axes.append(tuple(thisaxis))
            axisedges.append(CoordNode[tuple(thisaxis)])

        alreadydone = []


        allcoords_ranked = []        
        for j in range(nodelen**self.NDim):
            thisvec = self.BaseN(j, nodelen, self.NDim)
            allinterven = allintervening(thisvec, cubelen)
            nintervening =len(allinterven)
            
            allcoords_ranked.append((nintervening, allinterven, tuple(thisvec)))
            if modall(thisvec, cubelen): 
                oldcoord = [icomp//cubelen for icomp in thisvec]
                thisnode = CoordNode[tuple(oldcoord)]
                thiscoord = tuple(thisvec)                
                if not(thisvec in alreadydone):
                    alreadydone.append(tuple(thisvec))


        
        # now, get rig of the 2nd and 3rd components of allcoords_ranked, since we don't explicitly need them after they've been used in the sort
        allcoords_ranked.sort()
        newcoords = []
        for x,y,z in allcoords_ranked:
            newcoords.append((x, z))
        allcoords_ranked = newcoords
        # now, allcoords is ranked not only by sum or nonintervening, but also by the lowest intervening component

        # Next we will see if any of the required have already been created due to divisions
        # in adjacent D-cubes

        # first explicitly search for the nodes situated on a (previous) edge connecting two nodes of the cube;
        # note any path in coord_path that has more than 2 nodes, those extra interior nodes must correspond to pre-existing edge nodes

        for icoord, dictpaths in coord_path.items():
            scaledcoord = tuple([cubelen*i for i in icoord])
            for subcoord, edgepath in dictpaths.items():
                scaledsubcoord = tuple([cubelen*i for i in subcoord])
                coordpath = CoordRun(scaledsubcoord, scaledcoord, cubelen)
                if len(edgepath) == 2:
                    continue            
                for jcoord, jnode in zip(coordpath[1:-1], edgepath[1:-1]):
                    newcoordnode[jcoord] = jnode
                    newnodecoord[jnode] = jcoord 
                alreadydone.extend(coordpath) 
        #  now that you've found the pre-existing edges (i.e. 1-dim structures) from adjacent and previously divided cubes do all the higher dimensinal pre-existing structures


        
        previsum = 0
        for isum,icoord in allcoords_ranked:
            if icoord in alreadydone:
                continue
            if isum <= 1:
                continue # already did these      
            else:   
                if previsum != isum:
                    vecintervening = [NIntervening(x, nodelen) for x in alreadydone]
                    if not(isum-1 in vecintervening):
                        # if we have not found any matchings for some previous (i.e. lower) isum, there is no point in looking forther
                        # i.e. if no edges have been found to have intervening already-created nodes, there won't be any faces or anything with higher dimensions
                        #import pdb; pdb.set_trace()
                        break
                if isum >= self.NDim:
                    #print("IF YOU ARE HERE, IT MEANS THERE ARE ALREADY INTERIOR POINTS IN THIS to-be-DIVIDED REGION -- i.e. the check-for-interior-nodes test failed. Fix it")
                    #import pdb; pdb.set_trace()
                    # in 2d there can only be pre-existing edges (1-d objects), 
                    # in 2d there can only be pre-existing edges and faces (1 and 2-d objects)
                    # we will allow the isum==NDim case as an extra sanity check
                    
                    # you're done since there has already been a check for interior points
                    break                

                #find the two existing nodes which we will need to find intersection of
                relevantcoords = ToBeIntersected(icoord, cubelen)
                relevantneighbornodes = []
                bFailed = False
                for jcoord in relevantcoords:
                    if not(jcoord in newcoordnode):
                        bFailed = True
                        break
                    relevantnode = newcoordnode[jcoord]
                    relevantneighbornodes.append( self.NodeVec[relevantnode].Neighbors )
                
                if bFailed:                     
                    # for each value of isum  the corresponding nodes are a package deal -- either they all exist,
                    #  or none exist, so the first time it fails, there is no point in going further
                    continue



                intersected = list(set(relevantneighbornodes[0]).intersection(*relevantneighbornodes[1:]))

                copyintersected = copy(intersected)
                for intnode in copyintersected:
                    if intnode in newnodecoord.keys():
                        intersected.remove(intnode)
                


                if len(intersected) != 1:
                    continue # if you cannot find anything this far, no sense in going higher



                inode = intersected[0]
                newcoordnode[icoord] = inode
                newnodecoord[inode] = icoord
                alreadydone.append(tuple(icoord))
                previsum = copy(isum)






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
            nbrnodes = self.NodeVec[nodeA].Neighbors 
            if not(nodeB in nbrnodes):
                print("How is the once-removed vector not a neighbor of the zero node?")
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

            thistuple = tuple(thisvec)

            

            newcoordnode[thistuple] = thisnode.Id
            newnodecoord[thisnode.Id] = thistuple
            if bHistory:
                if thisnode.Coords is None:
                    crds, ratcrds = GetCoords(thistuple, Id, nslices, axisedges, CoordNode, NodeCoord)
                    thisnode.Coords = crds
                    #thisnode.RatCoords = ratcrds
                    for ind in self.NodeVec:
                        if ind.Id == thisnode.Id:
                            continue
                        if np.sum(np.abs(np.array(ind.Coords) - np.array(thisnode.Coords))) < 10e-9:                            
                            print("Possibly identical coords? ", thisnode.Id, ind.Id, "GET RID OF THIS TEST AFTER DEBUGGING; IT MUST NOT REMAIN")
                            import pdb; pdb.set_trace()
                else:
                    print("why are we here -- we should only be looping through new nodes that don't have coords")
       
        # now delete whichever old neighbors are part of the original (unsliced) cube        
        for key, val in NodeCoord.items():
            thisnode = key
            thisvec = list(val)
            thesenbrs = copy(self.NodeVec[thisnode].Neighbors)

            for nbr in thesenbrs:
                if nbr in NodeCoord.keys():
                    self.NodeVec[thisnode].Neighbors.remove(nbr)

            self.NodeVec[thisnode].Neighbors.sort()

        # Note the number of neighbors depends on the number of components in the vector that are 0 or nslices;
        # if that number zero, we're in the interior of the crube and the number of neighbors is 2*D, if it's 1 then
        # we're (in 3d) in the middle of  face then there are only 2*D-1, if it's 2 then it's 2*D-2, etc.

        for j in range(nodelen**self.NDim):
            thisvec = self.BaseN(j, nodelen, self.NDim)
            thisnode = newcoordnode[tuple(thisvec)]


            for iiv, iv in enumerate(thisvec):
                nbrvecUp = list(copy(thisvec))
                nbrvecDn = list(copy(thisvec))

                nbrvecUp[iiv] = iv + 1 # if greater than cubelen, then Up nbr already exists from before, and we leave this alone
                if iv + 1 <= cubelen:
                    nbrnodeUp = newcoordnode[tuple(nbrvecUp)]
                    if not(nbrnodeUp in self.NodeVec[thisnode].Neighbors):
                        # it's already in there from some adjacent cube-division
                        nbrnodeUp_prevA = copy(nbrvecUp)       
                        nbrnodeUp_prevB = copy(nbrvecUp)                
                        nbrnodeUp_prevA[iiv] = 0
                        nbrnodeUp_prevB[iiv] = cubelen
                        nodeId_prevA = newcoordnode[tuple(nbrnodeUp_prevA)]
                        nodeId_prevB = newcoordnode[tuple(nbrnodeUp_prevB)]

                        self.NodeVec[thisnode].Neighbors.append(nbrnodeUp)

                        if len(self.NodeVec[thisnode].Neighbors) > 2*self.NDim:
                            print("Neighbor count exceeded -- more than 2 * Dim? -- DELETE THIS CHECK")
                            import pdb; pdb.set_trace()


                nbrvecDn[iiv] = iv - 1 # if less than 0, then dn nbr already exists from before, and we leave this alone
                if iv - 1 >= 0:                    
                    nbrnodeDn = newcoordnode[tuple(nbrvecDn)]
                    self.NodeVec[thisnode].Neighbors
                    if not(nbrnodeDn in self.NodeVec[thisnode].Neighbors):
                        #continue # it's already in there from some adjacent cube-division
                        nbrnodeDn_prevA = copy(nbrvecDn)     
                        nbrnodeDn_prevB = copy(nbrvecDn)
                        nbrnodeDn_prevA[iiv] = 0
                        nbrnodeDn_prevB[iiv] = cubelen
                        nodeId_prevA = newcoordnode[tuple(nbrnodeDn_prevA)]
                        nodeId_prevB = newcoordnode[tuple(nbrnodeDn_prevB)]

                        prevnbrs = [nodeId_prevA, nodeId_prevB]
                        prevnbrs.sort()
                        self.NodeVec[thisnode].Neighbors.append(nbrnodeDn)

                        if len(self.NodeVec[thisnode].Neighbors) > 2 * self.NDim:
                            print("Neighbor count exceeded -- more than 2 * Dim? -- DELETE THIS CHECK")
                            import pdb; pdb.set_trace()

            self.NodeVec[thisnode].Neighbors.sort()

        for coord,node in newcoordnode.items():
            if len(self.NodeVec[node].Neighbors) == 2:
                print("too low")
                import pdb; pdb.set_trace()
            
            self.NodeVec[node].LastDivide = self.NDivisions
        self.NDivisions += 1
        #print("done with divide")

        # finally, update the region count of all the participating nodes
        
        
        for inode in newnodecoord.keys():
            if self.NodeVec[inode].NRegions >= 2 ** self.NDim:
                continue # already reached terminal saturation/RegionCount
            nbrsinnew = 0
            for inbr in self.NodeVec[inode].Neighbors:
                if inbr in newnodecoord.keys():
                    nbrsinnew += 1
            
            thisregionincrement = NodeCountToHypercubeCount(nbrsinnew, self.NDim)
            self.NodeVec[inode].NRegions += thisregionincrement
            if self.NodeVec[inode].NRegions > 2 ** self.NDim:
                print("How did this happen,  since NRegions can be at most 2**D")
                import pdb; pdb.set_trace()
            #else:
            #    self.NodeVec[inode].RegionGrowth.append(thisregionincrement)
            #    self.NodeVec[inode].RegionGrowth.sort()
        
                






        
    def bIsCorner(self, nodeId):
       return len(self.NodeVec[nodeId].Neighbors) == 2 * self.NDim
    



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
                thisupnbr = coordnode[tuple(revvec)]
                prevnbrlist = [id, thisupnbr]
                prevnbrlist.sort()
                self.NodeVec[id].Neighbors.append(thisupnbr)
                revvec = UpDn(-1, id, idim, torlen)
                thisdnnbr = coordnode[tuple(revvec)]
                prevnbrlist = [id, thisdnnbr]
                prevnbrlist.sort()

                self.NodeVec[id].Neighbors.append(thisdnnbr)
        
        
        for id in range(torlen**self.NDim):
            self.NodeVec[id].Neighbors.sort()
            
                
        idlistsofar = [xnode.Id for xnode in self.NodeVec]
        for inode in idlistsofar:
            thesecubenodes = [  ]
            allaxes = self.ReturnAllAxes(inode)

            for theseaxes in allaxes:
                
                bTryPrevNeighborsToo = False       
                (loccoordnode, locnodecoord) = self.FindCubeCoord(inode, theseaxes)


        NodeCoord = nodecoord
        CoordNode = coordnode 
    

        for Id, coord in NodeCoord.items():
            self.NodeVec[Id].Coords = coord
            self.NodeVec[Id].RatCoords = tuple([RationalNumber(xint) for xint in coord])
            self.NodeVec[Id].NRegions =  2 ** self.NDim 
            #self.NodeVec[Id].RegionGrowth = [2 ** self.NDim]

    
    
    def bSimpleIsItAFace(self, Id, nodeA, nodeB):
        Anbrs = self.NodeVec[nodeA].Neighbors
        Bnbrs = self.NodeVec[nodeB].Neighbors 
        intersected =  list(set(Anbrs).intersection(Bnbrs))
        if Id in intersected:
            intersected.remove(Id)
        else:
            print(" bSimpleIsItAFace should only be called for a node and two neigbors of that node:", Id, nodeA, nodeB)
            import pdb; pdb.set_trace()
        if len(intersected) == 1:
            return True
        return False


    # Note: this is actually a misnomer in dimensions higer than 2, where even nodes with 2D neighbors can still be unpacked, so to speak
    def NearestPackedNeighbor(self, Id, nbr):
        if len(self.NodeVec[nbr].Neighbors) == 2 * self.NDim:
            return nbr
        else:
            thisrow, bIsOK = self.FarNeighbor(Id, nbr)
            #if len(thisrow) <= 2:
            #    print("This should return nslices+1 numbers -- what is wrong?")
            #    import pdb; pdb.set_trace()
            return thisrow[-1]     
               






    def BruteNeighborSearch(self, iPath, PathLen):
        if PathLen <= 2:
            return [ iPath ]
        
        # the following works in 2-d but gets you in trouble in higher dimensions
        #TgtNbrCount = len(self.NodeVec[iPath[-1]].Neighbors) 
        #if TgtNbrCount == 2 * self.NDim:
        #    return [ iPath ]

        TgtRegionCount = self.NodeVec[iPath[-1]].NRegions
        #TgtRegionGrowth = tuple(self.NodeVec[iPath[-1]].RegionGrowth)
        if TgtRegionCount == 2 ** self.NDim:
            return [ iPath ]        

        buildablepaths = [iPath]
        finishedpaths = []


        #import pdb; pdb.set_trace()
        for jextend in range(0, PathLen-2):
            dellist = []
            addedpaths = []
            for jjpath, jpath in enumerate(buildablepaths):
                if len(jpath) >= PathLen:
                    continue
                
                jnbrs = copy(self.NodeVec[jpath[-1]].Neighbors)
                for knbr in jnbrs:
                    if knbr in jpath:
                        continue
                    #kNbrLen = len(self.NodeVec[knbr].Neighbors)
                    #if kNbrLen == 2 * self.NDim:
                    kRegCnt = self.NodeVec[knbr].NRegions
                    #kRegGrowth = tuple(self.NodeVec[knbr].RegionGrowth)
                    if kRegCnt == 2 ** self.NDim:
                        if len(jpath) + 1 == PathLen:
                            cappedpath = jpath + [knbr]
                            finishedpaths.append(cappedpath)
                               
                        # in this case we have reached the end of the path
                        addedpaths.append( jpath + [knbr] )
                    elif kRegCnt == TgtRegionCount:
                        addedpaths.append( jpath + [knbr] )
                        # in this case, we continue the path along the same edge (as specified earlier in TgtRegionCount)


            if len(addedpaths) == 0:
                break
            buildablepaths = copy(addedpaths)
        return finishedpaths


    def SphereSweep(self, Id, NMax=100000, Radius=None):
        if Radius is None:
            Radius = (self.NSlices + 2) * self.NDim
        Sphere = [Id]

        myretval = []
        
        for irad in range(Radius):
            copysphere = copy(Sphere)
            for inode in copysphere:
                nbrs = self.NodeVec[inode].Neighbors
                for inbr in nbrs:
                    if inbr in Sphere:
                        continue
                    Sphere.append(inbr)
                    if self.NodeVec[inbr].NRegions != 2 ** self.NDim:
                        continue
                        
                    test_coordnode_nodecoord_coordpath_list = self.ReturnAllWellFormedCubes(inbr)
                    for jscenario in range(len(test_coordnode_nodecoord_coordpath_list)):
                        test_icoordnode, test_inodecoord, test_icoord_path = test_coordnode_nodecoord_coordpath_list[jscenario]
                        bisthismixed = self.bIsMixedFractality(test_icoord_path)

                        if bisthismixed:
                            if not(self.bIsNonUniform(test_icoordnode)):
                                print("why are bIsMixedFractality and bIsNonUniform not coinciding?")
                                import pdb; pdb.set_trace()
                            myretval.append((test_icoordnode, test_inodecoord, test_icoord_path))


            if len(myretval) > 0:
                print("SphereSweep found", str(len(myretval)), "workable mixed scenarios within sweep level", irad)
                return myretval
            if len(Sphere) > NMax:
                break
        
        return None
                    

            


        
        # the following works in 2-d but gets you in trouble in higher dimensions
        #TgtNbrCount = len(self.NodeVec[iPath[-1]].Neighbors) 
        #if TgtNbrCount == 2 * self.NDim:
        #    return [ iPath ]

        TgtRegionCount = self.NodeVec[iPath[-1]].NRegions
        #TgtRegionGrowth = tuple(self.NodeVec[iPath[-1]].RegionGrowth)
        if TgtRegionCount == 2 ** self.NDim:
            return [ iPath ]        

        buildablepaths = [iPath]
        finishedpaths = []


        #import pdb; pdb.set_trace()
        for jextend in range(0, PathLen-2):
            dellist = []
            addedpaths = []
            for jjpath, jpath in enumerate(buildablepaths):
                if len(jpath) >= PathLen:
                    continue
                
                jnbrs = copy(self.NodeVec[jpath[-1]].Neighbors)
                for knbr in jnbrs:
                    if knbr in jpath:
                        continue
                    #kNbrLen = len(self.NodeVec[knbr].Neighbors)
                    #if kNbrLen == 2 * self.NDim:
                    kRegCnt = self.NodeVec[knbr].NRegions
                    #kRegGrowth = tuple(self.NodeVec[knbr].RegionGrowth)
                    if kRegCnt == 2 ** self.NDim:
                        if len(jpath) + 1 == PathLen:
                            cappedpath = jpath + [knbr]
                            finishedpaths.append(cappedpath)
                               
                        # in this case we have reached the end of the path
                        addedpaths.append( jpath + [knbr] )
                    elif kRegCnt == TgtRegionCount:
                        addedpaths.append( jpath + [knbr] )
                        # in this case, we continue the path along the same edge (as specified earlier in TgtRegionCount)


            if len(addedpaths) == 0:
                break
            buildablepaths = copy(addedpaths)
        return finishedpaths

    def AllHave2DNeighbors(self, somelist):
        for inode in somelist:
            if len(self.NodeVec[inode].Neighbors) != 2 * self.NDim:
                return False
        return True


    def AllHave2totheDRegions(self, somelist):
        for inode in somelist:
            if len(self.NodeVec[inode].Neighbors) != 2 ** self.NDim:
                return False
        return True





    def ReturnAllWellFormedCubes(self, Id, PathLen = None, ForceAxis=None):
        """
        All corners of the cube must be  connected by paths of length NSlices+2
        or else 2 (in the latter case, the respective nodes are nearest neighbors);
        moreover all those intervening nodes must have the same number of neighbers,
        with the number being less than 2D

        To avoid problems arising from a region of fractality+1 only a step away from region
        of fractality-1, we must also impose the following restriction. If two paths have the
        same endpoints, then any of paths that is of length nslices+2 is excluded (i.e. if there
        exists a path in that same-endpoints set that is of length two -- so that it consists of
        nearest neighbors -- it will be retained, whereas all longer paths in that set are excluded.)

        Also, there can be no interior nodes -- and note interior noteds are only
        possible if ALL the paths connecting corners of the cube are of len NSlices+2
        """

        def IsAxisContained(iaxx, ForceAxis):
            alreadydone = []
            for ifax in ForceAxis:
                for iiax, iax in enumerate(iaxx):
                    if ifax in iax:
                        alreadydone.append(iiax)
            return (len(set(alreadydone)) == len(alreadydone)) and (len(alreadydone) == len(ForceAxis))
                
            
                


        if PathLen is None:
            PathLen = self.NSlices + 2


        nbrs = copy(self.NodeVec[Id].Neighbors)
        lennbrs = len(nbrs)
        if self.NDim == 1:
            return [ [Id, nbrs[0]], [Id, nbrs[1]] ]
        #if len(self.NodeVec[Id].Neighbors) != 2 * self.NDim:
        #    return []
        if self.NodeVec[Id].NRegions != 2 ** self.NDim:
            return []

        myretval = []

        nodecoord = {}
        coordnode = {}
        coordorigin = tuple([0] * self.NDim)
        nodecoord[Id] = tuple([0] * self.NDim)
        coordnode[coordorigin] = coordorigin


        # later on, the early part of routine was packaged into a separate routein (Chrysanthemum())
        sortbyones = self.SortByOnes

        pathdump = []
        pathdumpgeneration = []
        #
        alreadydone = []
        for inbr in nbrs:
            thesepaths = self.BruteNeighborSearch([Id, inbr], PathLen)
            for ipath in thesepaths:
                pathdump.append(ipath)
                pathdumpgeneration.append(0)
                alreadydone.append( (tuple(ipath), 0) )

        

        
        alreadydone = []
        for idim in range(1, self.NDim):
            iidump = -1
            for idump, idumpgen in zip(pathdump, pathdumpgeneration):
                iidump += 1
                if idumpgen != idim-1:
                    continue

                jnbrs = self.NodeVec[idump[-1]].Neighbors
                for jnbr in jnbrs:
                    #thisarg = (tuple(idump), idim)
                    #if idim == 1 and idump[-1] == 72 and jnbr == 522:
                    #    print("this gets we")
                    #    import pdb; pdb.set_trace()
                    thesepaths = self.BruteNeighborSearch([idump[-1], jnbr], PathLen)                       
                    for jp in thesepaths:
                        thisarg = (tuple(jp), idim)
                        intersected = list(set(idump).intersection(set(jp[1:])))
                        if len(intersected) == 0 and not( thisarg in alreadydone  ) :
                            pathdump.append(jp)
                            pathdumpgeneration.append(idim)
                            alreadydone.append( thisarg )


                        #else:
                            #alreadydone[(iidump, len(pathdump)-1)] = False

        #if self.NNode >= 612:
        #    print("look for 72")
        #    import pdb; pdb.set_trace()

        # now, delete any "long" paths that have identical endpoints 
        # (any "short" paths -- i.e. nearest neighbor paths -- will of necessity be unique)
        # NOTE: CONSIDER removing block since this sometimes winds up killing legit "straight" paths too, though the very fact that alternate kinky
        # versions begin and end at same point indicates we're in a defective region when it comes to trying to divide.
        dellist = []
        for iipath, ipath in enumerate(pathdump):
            iendpts = [ipath[0], ipath[-1]]
            iendpts.sort()
            iendpts = tuple(iendpts)
            biDelete = False
            for jjpath in range(iipath + 1, len(pathdump)):
                jpath = pathdump[jjpath]
                jendpts = [jpath[0], jpath[-1]]
                jendpts.sort()
                jendpts = tuple(jendpts)
                if iendpts == jendpts and pathdumpgeneration[iipath] == pathdumpgeneration[jjpath]:
                    #if len(jpath) == 4 and len(ipath) == 4:
                    #    if (self.l2node(jpath[0], jpath[1]) * 3 ==  self.l2node(jpath[0], jpath[-1])) or (self.l2node(ipath[0], ipath[1]) * 3 ==  self.l2node(ipath[0], ipath[-1])):
                    #        print("one is straight?")
                     #       import pdb; pdb.set_trace()
                    if len(jpath) > 2:
                        dellist.append(jjpath)
                    if len(ipath) > 2:
                        dellist.append(iipath)
                        biDelete = True
                        break
            if biDelete:
                continue
        dellist = list(set(dellist))
        dellist.sort()
        dellist.reverse()

        for idel in dellist:
            del pathdump[idel]
            del pathdumpgeneration[idel]

          


        alreadydone = {}       



        pathdumpgeneration = np.array(pathdumpgeneration).astype("int")
        gen0len = np.sum(pathdumpgeneration == 0) # these can be the D edges connecting to the corner node
        
        allaxcombos = AllUniqueCombosNoReplacement(gen0len, self.NDim)
        
        """
        print("Here are the axes we're about to test for ", Id)
        for icombo in allaxcombos:
            iax = [pathdump[iicombo] for iicombo in icombo]

            # check that all their nodes (except of course the starting) are distinct
            iaxx = [pathdump[iicombo][1:] for iicombo in icombo]
            intersected = list(set(iaxx[0]).intersection(*iaxx[1:]))
            if len(intersected) != 0:
                continue
            print(Id, icombo, iaxx)
        print("now we test them")
        """
        TestForInteriorPoints = {}
        for icombo in allaxcombos:

            iax = [pathdump[iicombo] for iicombo in icombo]

            # check that all their nodes (except of course the starting) are distinct
            iaxx = [pathdump[iicombo][1:] for iicombo in icombo]
            if not(ForceAxis is None):
                if not(IsAxisContained(iaxx, ForceAxis)):
                    continue

            intersected = list(set(iaxx[0]).intersection(*iaxx[1:]))
            if len(intersected) != 0:
                continue
            # the axis is like the "independent" paths and unambiguously start at the origin node
            # -- the rest are all "dependent" paths whose points  of origin depend on which independent path was chosen
            #if Id == 1087:
            #    import pdb; pdb.set_trace()
            #if Id == 1087 and tuple(iaxx[0]) == (1088) and tuple(iaxx[1]) == (1091):
            #    import pdb; pdb.set_trace()

            #if Id == 2055 and icombo == (0,1,2) and self.NNode == 7816:
            #    print("aaa")
            #    import pdb; pdb.set_trace()

            
            coordnode, nodecoord, coord_path = self.BruteFindCubeCoord(Id, iax, icombo, pathdump, pathdumpgeneration, gen0len, TestForInteriorPoints, sortbyones)
            if not(coordnode is None):
                myretval.append((coordnode, nodecoord, coord_path))  

            
            
        print("Well-formed cubes wth count", len(myretval), "at", Id, len(myretval), [i[2] for i in myretval])
        return myretval



    def ReturnAllAxesNoRestrictions(self, Id, bAllowPermutations=True):

        nbrs = copy(self.NodeVec[Id].Neighbors)
        lennbrs = len(nbrs)
        if self.NDim == 1:
            return [ nbrs ]

        # if corner of node does not have 2D neighbors, we cannot divide the cube
        #if len(self.NodeVec[Id].Neighbors) != 2 * self.NDim:
        if self.NodeVec[Id].NRegions != 2 ** self.NDim:
            return []
        # Likewise, check that the axes either have 2D neighbors or that their far neighbors (relative to Id)
        # have 2D neighbors; if neither is true, we can't form a good cube with it as an edge, 
        # so remove that node from mbrs  

        myret = []
        
        indexcombos = AllCombosNoReplacement(len(nbrs), self.NDim)
        combos = []
        for icomb in indexcombos:
            combos.append([nbrs[i] for i in icomb])
        
        if not(bAllowPermutations):
            alreadydone = []
            dellist = []
            for iicomb, icomb in enumerate(combos):
                ilist = list(icomb)
                ilist.sort()
                if tuple(ilist) in alreadydone:
                    dellist.append(iicomb)
                    continue
                alreadydone.append(tuple(ilist))
            dellist.reverse()
            for idel in dellist:
                del combos[idel]

        return combos
        

    def ReturnAllAxes(self, Id):

        nbrs = copy(self.NodeVec[Id].Neighbors)
        lennbrs = len(nbrs)
        if self.NDim == 1:
            return [ nbrs ]

        # if corner of node does not have 2D neighbors, we cannot divide the cube
        #if len(self.NodeVec[Id].Neighbors) != 2 * self.NDim:
        if self.NodeVec[Id].NRegions != 2 ** self.NDim:
            return []
        # Likewise, check that the axes either have 2D neighbors or that their far neighbors (relative to Id)
        # have 2D neighbors; if neither is true, we can't form a good cube with it as an edge, 
        # so remove that node from mbrs  


       
        dellist = []
        for inbr in nbrs:
            #if len(self.NodeVec[inbr].Neighbors) != 2 * self.NDim:
            if self.NodeVec[inbr].NRegions != 2 ** self.NDim:
                inbr_farnbr_strip = self.AcrossADividedEdge(Id, inbr)
                import pdb; pdb.set_trace()
                if len(inbr_farnbr_strip) != 2 + self.NSlices: #  2 * self.NDim:
                    dellist.append(inbr)
                    continue
                if self.AcrossADividedEdge(inbr_farnbr_strip[-1], inbr_farnbr_strip[-2])[-1] != Id:
                #if len(self.NodeVec[inbr_farnbr].Neighbors) != 2 * self.NDim:
                    dellist.append(inbr)
        for idel in dellist:
            nbrs.remove(idel)


        myret = []

        indexcombos = AllCombosNoReplacement(len(nbrs), self.NDim)
        combos = []
        for icomb in indexcombos:
            combos.append([nbrs[i] for i in icomb])
        

        
        dellist = []

        collinearitytest = {}
       
        dellist = []
        alreadydone = {}

        index2combos =  AllUniqueCombosNoReplacement(self.NDim, 2)
        for iicombo, icombo in enumerate(combos):
            #if self.NNode >= 13600 and Id == 11759 and (tuple(icombo) == (11758,12722) or tuple(icombo) == (12722, 11758) ):
            #    print("boooab")
            #    import pdb; pdb.set_trace()


            sortcombo = list(icombo)
            sortcombo.sort()
            sortcombo = tuple(icombo)
            if sortcombo in alreadydone:
                if not(alreadydone[sortcombo]):
                    dellist.append(iicombo)
                    continue

            bCollinearAxes = False
            for i2combo in index2combos:
                nodeA = icombo[i2combo[0]]
                nodeB = icombo[i2combo[1]]
                listAB = [nodeA, nodeB]
                listAB.sort()
                tlistAB = tuple(listAB)
                if tlistAB in collinearitytest:
                    ibAreCollinear = collinearitytest[tlistAB]
                else:
                    # axes are badly oriented if they're collinear or
                    # if the cube is in a midway-fractality channel between higher and lower fractality on either side

                    bAreCollinear = self.bAreTwoAxesBadlyOriented(Id, [nodeA, nodeB])
                    collinearitytest[tlistAB] = bAreCollinear
                if bAreCollinear:
                    bCollinearAxes = True
                    break
            if bCollinearAxes:
                dellist.append(iicombo)
                alreadydone[sortcombo] = False    
                continue
                            
            # first, check for interior nodes by trying to find opposite node without any farneighbor extensions
            coordnode, nodecoord = self.FindInterior(Id, icombo)
            if not(coordnode is None):
                if self.AllHave2DNeighbors(nodecoord.keys()) and self.AllHave2totheDRegions():
                    alreadydone[sortcombo] = True
                    continue # an acceptable set of coords
                else:
                    dellist.append(iicombo)
                    alreadydone[sortcombo] = False
                    continue

            # if we make it here then no interior point was found, and we do a more thoroguh (i.e. using far-neighbor extensions) effort to  construct a cbe
            #if 1049 in icombo and 1030 in icombo and  self.NNode == 1072:
            #    import pdb; pdb.set_trace()

            coordnode, nodecoord = self.FindCubeCoord(Id, icombo)

            if coordnode is None:
                dellist.append(iicombo)
                alreadydone[sortcombo] = False
            else:
                alreadydone[sortcombo] = True
            
                    
        
        dellist.sort()
        dellist.reverse()
        for idel in dellist:
            del combos[idel]

        # Note that if (x,y) is in combos, so is (y, x); they're different axes

        return combos


      


        

    
    
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
            inbrs = self.NodeVec[inode].Neighbors
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


    # this is just a customized copy of FindCubeCoord() -- SEEMS TO BE UNUSUED
    def FindOppositeNode(self, Id, nbrsubset, bAllowFarNeighbors=True, subdim=None):
        """
        This is only ever used (so far) with bAllowFarNeighbors set to False, but it may eventually have other uses.
        In 2d this is very similar to bIstItAFace
        """
        if subdim is None:
            subdim = self.NDim

        if self.NodeVec[Id].NRegions != 2 ** self.NDim:
            return None, None         
    
        Idcoord = tuple([0] * self.NDim)
        coordnode = {Idcoord:Id}
        nodecoord = {Id:Idcoord}

        for iic, inbrId in enumerate(nbrsubset):
            thisvec = np.zeros((subdim,)).astype("int")
            thisvec[iic] = 1
            thisvec = tuple(thisvec)
            coordnode[thisvec] = inbrId
            nodecoord[inbrId] = thisvec

        sortbyones = []        
        for j in range(2**subdim):
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
                
                if bAllowFarNeighbors:
                    print("REPLACE WITH ACROSSADIVIDEDEDGE")
                    thisnode = self.NearestPackedNeighbor(Id, self.nbrsfrombin(nbrsubset, binvec, coordnode, nodecoord)[0])
                else:
                    thisnode = coordnode[binvec]

                coordnode[binvec] = thisnode
                nodecoord[thisnode] = binvec
                continue
            if sum1 >= 2:
                thesenodes = tuple(self.nbrsfrombin(nbrsubset, binvec, coordnode, nodecoord)) # note there D elements are  just the sum1 components of thesenodes
                nbrs = []

                for itup in thesenodes:
                    subsetnbrs = []
                    thesenbrs = self.NodeVec[itup].Neighbors
                    for jnbrId in thesenbrs:
                        if bAllowFarNeighbors:
                            subsetnbrs.append( self.NearestPackedNeighbor(itup, jnbrId) )
                        else:
                            subsetnbrs.append( jnbrId )

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
                    return None, None # rerun with minmag

                if self.NodeVec[thisnode].NRegions != 2 ** self.NDim:
                    return  None, None # rerun with minmag
                thisnode = intersected[0]
                coordnode[binvec] = thisnode
                nodecoord[thisnode] = binvec

                continue

            nextnodes = []            
            for jsum1, jbinvec in sortbyones[1:]:
                if jsum1 == sum1-1:
                    thisnode = coordnode[jbinvec]
                    for jnbr in self.NodeVec[thisnode].Neighbors:
                        if len(self.NodeVec[jnbr].Neighbors) == 2 * self.Ndim:
                            nextnodex.append(jnbr)
                        else:
                            nextnodex.append(self.FartherNeighbor(thisnode, jnbr)[0])
            

            intersected = list(set(nextnodes[0]).intersection(*nextnodes))
            if len(intersected) != 1:
                return None, None

            thisnode = intersected[0]
            coordnode[binvec] = thisnode
            nodecoord[thisnode] = binvec

            if self.NodeVec[thisnode].NRegions != 2 ** self.NDim:
                return None, None


        
        return coordnode, nodecoord



    # this is just a customized copy of FindCubeCoord()
    def FindInterior(self, Id, nbrsubset):
        """
        This is just a stripped-down version of find opposite node
        Has no far-neighbor search; has no checks to see that the returned nodes have 2D niehgbors;
        that can all easily be added in the calling routine
        """

     
    
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
                
                thisnode = coordnode[binvec]

                coordnode[binvec] = thisnode
                nodecoord[thisnode] = binvec
                continue
            if sum1 >= 2:
                thesenodes = tuple(self.nbrsfrombin(nbrsubset, binvec, coordnode, nodecoord)) # note there D elements are  just the sum1 components of thesenodes
                nbrs = []
                #import pdb; pdb.set_trace()
                for itup in thesenodes:
                    subsetnbrs = []
                    thesenbrs = self.NodeVec[itup].Neighbors
                    for jnbrId in thesenbrs:
                        subsetnbrs.append( jnbrId )

                    nbrs.append(subsetnbrs)
                    
                    
                
                intersected = list(set(nbrs[0]).intersection(*nbrs[1:]))

                dellist = []
                for ixtec in intersected:
                    if ixtec in nodecoord.keys():
                        dellist.append(ixtec)
                for idel in dellist:                    
                    intersected.remove(idel) 
                    if (idel, ) != (Id,) and len(idel) != 0:
                            print("dellist should only contain 1 element! And it should be Id -- what's going on?")

                if len(intersected) != 1:
                    return None, None 
        
                thisnode = intersected[0]
                coordnode[binvec] = thisnode
                nodecoord[thisnode] = binvec

                continue

            nextnodes = []            
            for jsum1, jbinvec in sortbyones[1:]:
                if jsum1 == sum1-1:
                    thisnode = coordnode[jbinvec]
                    for jnbr in self.NodeVec[thisnode].Neighbors:
                        if len(self.NodeVec[jnbr].Neighbors) == 2 * self.Ndim:
                            nextnodex.append(jnbr)
                        else:
                            nextnodex.append(self.FartherNeighbor(thisnode, jnbr)[0])
            

            intersected = list(set(nextnodes[0]).intersection(*nextnodes))
            if len(intersected) != 1:
                return None, None

            thisnode = intersected[0]
            coordnode[binvec] = thisnode
            nodecoord[thisnode] = binvec
        
        return coordnode, nodecoord


   
    def FindCubeCoord(self, Id, nbrsubset):
        global FreshProb
        """

        generate the hypercube coordinates,
        and sort them according to how many 1's they have,
        since the calc of those with a given number of ones (say M)
        are the intersection of the neighbors of the nodes
        generated by nodes correesponding to a number of ones that is M-1

        NOTE: This has been deprecated in favor of ReturnWellFormedCubes, and without
        a neighbor_t class and PrevNeighbors, sometimes fails. It  has been retained
        in order to allow for easier comparison between this version of the code and
        The one using the aforementioned neighbor_t class. In cases where it fails,
        run FindWellFormedCube and use that and BruteDivid to shape the desired node
        configuration.


        
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

        #if self.NNode == 1036:
        #    import pdb; pdb.set_trace()

        if self.NodeVec[Id].NRegions != 2 ** self.NDim:
            #print("ERROR --  the corners of the cube do not have 2*D neighbors?")
            return None, None         
    




        Idcoord = tuple([0] * self.NDim)
        coordnode = {Idcoord:Id}
        nodecoord = {Id:Idcoord}



        
        bUsingFarNeighbor = False
        for iic, inbrId in enumerate(nbrsubset):
            thisvec = np.zeros((self.NDim,)).astype("int")

            relevantnode = inbrId
            if inbrId in self.NodeVec[Id].Neighbors:
                #fails when Id is 885 and icombo is [2021, 2025]
                if self.NodeVec[inbrId].NRegions != 2 ** self.NDim: # and not(bOnlyCheckForInterior):
                    #relevantnode = self.FartherNeighbor(Id, inbrId)[0]
                    thisrow, bIsOK  = self.FarNeighbor(Id, inbrId)
                    if not(bIsOK):
                        return None, None
                    relevantnode = thisrow[-1]   
                    bUsingFarNeighbor = True
            try:
                if self.NodeVec[relevantnode].NRegions != 2 ** self.NDim: # and not(bOnlyCheckForInterior):
                    return None, None
            except:
                import pdb; pdb.set_trace()



            thisvec[iic] = 1
            thisvec = tuple(thisvec)
            coordnode[thisvec] = relevantnode
            nodecoord[relevantnode] = thisvec
            if self.NodeVec[relevantnode].NRegions != 2 ** self.NDim:
                #print("ERROR --  the corners of the cube do not have 2*D neighbors?")
                return None, None 


        sortbyones = [] # sortying by number of "middle" or intervening components (the "ones" is a misnomer except in case where nslices=1)
        
        for j in range(2**self.NDim):
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
                continue 
                

            if sum1 >= 2:


                #if binvec == (1,1,1):
                #    import pdb; pdb.set_trace()
                thesenodes = tuple(self.nbrsfrombin(nbrsubset, binvec, coordnode, nodecoord)) # note there are sum1 components of thesenodes

                nbrs = []

                #import pdb; pdb.set_trace()
                for itup in thesenodes:




                    subsetnbrs = []
                    thesenbrs = self.NodeVec[itup].Neighbors
                    for jnbrId in thesenbrs:
                        if len(self.NodeVec[jnbrId].Neighbors) == 2 * self.NDim:
                            subsetnbrs.append( jnbrId )
                        else:
                            thisrow, bIsOK = self.FarNeighbor(itup, jnbrId)
                            if not(bIsOK):
                                return None, None
                            farneighbor = thisrow[-1]
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
                    return None, None # rerun with minmag

                if self.NodeVec[thisnode].NRegions != 2 ** self.NDim:
                    return  None, None # rerun with minmag
                thisnode = intersected[0]
                coordnode[binvec] = thisnode
                nodecoord[thisnode] = binvec
                if self.NodeVec[thisnode].NRegions != 2 ** self.NDim:
                    #print("ERROR --  the corners of the cube do not have 2*D neighbors?")
                    return None, None 

                continue
                

        


        # begin check
        bCheck = False
        if bCheck:
            for coord,node in coordnode.items():
                if self.NodeVec[node].NRegions != 2 ** self.NDim:
                    print("wrong1 -- failed coordinate creation")
                    #import pdb; pdb.set_trace()
            for coord,node in coordnode.items():
                nbrs = self.NodeVec[node].Neighbors
                for inbr in nbrs:
                    if inbr in nodecoord.keys():
                        continue
                    farnbr = self.FartherNeighbor(node, inbr)[0]
                    if not(farnbr in nodecoord.keys()):
                        if node >= self.TorLen ** self.NDim:
                            print("wrong2")
                            import pdb; pdb.set_trace()
        # end check



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

        return coordnode, nodecoord

    
    def DecomposeBinaryIntoSmallerDimensionalBinaries(self, binvec):
        myretval = []
        for iic, ic in enumerate(binvec):
            if ic == 1:
                thisvec = list(binvec)
                if thisvec[iic] == 1:
                    thisvec[iic] = 0
                    myretval.append(tuple(thisvec))
        return myretval




    def AddSegment(self, prevsegment, pathdump, pathdumpgeneration):

        def RangeGeneration(whichgeneration, pathdumpgeneration):
            if whichgeneration == 0:
                n = len([i==0 for i in pathdumpgeneration])
                return (0, n)
            nprev = np.sum([i<=(whichgeneration-1) for i in pathdumpgeneration])
            n = nprev + np.sum([ i == whichgeneration for i in pathdumpgeneration[nprev:]])
            return (nprev, n)


        for key, ipath in prevsegment[0].items():
            whichgen = np.sum(key) + 1 # because the argument is PREsegment, the generation is actually one plus that.
            segmentA = pathdump[ipath]
            break

        myretval = []            
        endpt = segmentA[-1]            
        ifrom, ito = RangeGeneration(whichgen, pathdumpgeneration)
        iip = ifrom - 1

        
        for ip in pathdump[ifrom:ito]:
            iip += 1
            if ip[0] == endpt:
                #thisseq = list(prevsegment[0]) + [iip]
                #thisseqx = iip
                segmentB = pathdump[iip][1:]
                intersected =  list(set(segmentA).intersection(segmentB))  
                if len(intersected) == 0:
                    myretval.append((iip, ip[-1]))
        
        return myretval  



    def SimpleIntersection(self, PathSegments):
        myretval = []

        countbykeys = copy(PathSegments)
        for key in PathSegments.keys():
            countbykeys[key] = len(PathSegments[key])

        dictit = dictionaryiterator_t(countbykeys)

        for iic in range(dictit.Period):
            thisdict = copy(dictit.Current)
            thisendpt = -1
            thisentry = {}

            bMatchForEveryKey = True
            for key in thisdict.keys():
                ival = PathSegments[key][dictit.Current[key]]
                if thisendpt < 0:
                    thisendpt = ival[-1]
                    thisentry[key] = ival[0]
                else:
                    if ival[-1] != thisendpt:
                        bMatchForEveryKey = False
                        break
                    else:
                        thisentry[key] = ival[0]
                
                    
            if bMatchForEveryKey:
                # check for redundant entries
                myretval.append((thisentry, thisendpt))
            dictit.Increment()    
                
        return myretval

    def bBruteFindInterior(self, Id, nbrs, sortbyones=None):

        PathLen = 2
        pathdump = []
        pathdumpgeneration = []
        alreadydone = []
        for inbr in nbrs:
            thesepaths = self.BruteNeighborSearch([Id, inbr], PathLen)
            for ipath in thesepaths:
                pathdump.append(ipath)
                pathdumpgeneration.append(0)
                alreadydone.append( ((Id, inbr), 0) )

        alreadydone = []
        for idim in range(1, self.NDim):
            iidump = -1
            for idump, idumpgen in zip(pathdump, pathdumpgeneration):
                iidump += 1
                if idumpgen != idim-1:
                    continue

                jnbrs = self.NodeVec[idump[-1]].Neighbors
                for jnbr in jnbrs:
                    thisarg = ((idump[-1], jnbr), idim)
                    thesepaths = self.BruteNeighborSearch([idump[-1], jnbr], PathLen)                       
                    for jp in thesepaths:
                        """
                        intersected = list(set(idump).intersection(set(jp[1:])))
                        if len(intersected) == 0 and not( thisarg in alreadydone  ) :
                        """
                        if not( thisarg in alreadydone  ) :
                            pathdump.append(jp)
                            pathdumpgeneration.append(idim)
                            alreadydone.append( thisarg )

        pathdumpgeneration = np.array(pathdumpgeneration).astype("int")
        gen0len = np.sum(pathdumpgeneration == 0) # these can be the D edges connecting to the corner node



        zerocoord = tuple([0 for i in range(self.NDim)])

        for iic in range(self.NDim):        
            thiscoord = sortbyones[1 + iic][1]
            coord_path[thiscoord] = [({zerocoord:iic}, nbrs[iic])]   #we will retain the more elaborate format used in other routines even though things are simpler here since all segments must have 2 elements
    
            # one element for coords whose sum1==1, but higher ones may have more elements than one so we will mk it a list
            # in this case, too

        for sum1, binvec in sortbyones[1:]: # we already installed the first one
            if sum1 == 0:
                continue
            if sum1 == 1:
                continue     
                # we have ve already added the appropriate far neigbors to  coordnode and nodecoord so we don't need the following block                

            if sum1 >= 2:


                #if (Id == 4944 and self.NNode == 6192):
                #    import pdb; pdb.set_trace()

                tobeintersected = {}

                thesecoords = self.DecomposeBinaryIntoSmallerDimensionalBinaries(binvec)

                #import pdb; pdb.set_trace()
                bItIsEmpty  = False
                for itup in thesecoords:
                    if len(coord_path[itup]) == 0:
                        coord_path[tuple(binvec)] = []
                        bItIsEmpty = True
                        break
                    

                    for ip in coord_path[itup]:
                        
                        ipextrax = self.AddSegment(ip, pathdump, pathdumpgeneration)
                        if not(ipextrax is None):
                            tobeintersected[itup] = ipextrax

                if bItIsEmpty:
                    continue
            

                if np.min([ len(i)  for i in tobeintersected.values()]) == 0:
                    return False

                #import pdb; pdb.set_trace()
                coord_path[binvec] = self.SimpleIntersection(tobeintersected)

                if len(coord_path[binvec]) == 0:
                    return False          

                if len(coord_path[binvec]) > 1:
                    print("ERROR in bBruteFindInterior: More than one possible even when all the paths are two-tuple nearest neighbors?", nbrsubset)
                    #import pdb; pdb.set_trace()


        
        contenders = []
        NContenders = len(coord_path[self.OnesCoord])
        if NContenders > 1:
            print("why is bBruteFindInterior finding multiple interior nodes even when all the edges are nearest neighbors?")
            import pdb; pdb.set_trace()
        return NContenders >= 1

  


    def BruteFindCubeCoord(self, Id, nbrsubset, pathdump_index_subset, pathdump, pathdumpgeneration, gen0len, TestForInteriorPoints, sortbyones=None):
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


        def GetEdgeDict(sub_coord_path):
            # since this occurs after pruning, we assume that all lists in sub_coord_path have only one element
            edge_dict = {}
            for ikey, ilist in sub_coord_path.items():
                subdict = ilist[0][0]
                for isubkey, ipathind in subdict.items():
                    newlist = [ikey, isubkey]
                    newlist.sort()
                    newtup = tuple(newlist)
                    edge_dict[newtup] = ipathind
            return edge_dict

        def ChangePathIndicesToNodeLists(sub_coord_path, pathdump):
            retdict = copy(sub_coord_path)
            for key, list_tupdict_endnode in sub_coord_path.items():
                thisdict = list_tupdict_endnode[0][0]
                newthisdict = copy(thisdict)
                for subkey, pathindex in thisdict.items():
                    thispath = pathdump[pathindex]
                    newthisdict[subkey] = thispath
                retdict[key] = newthisdict
            return retdict

        
        # used to weed out deformed/degenerate cube contenders
        def bSomeEdgesAreSameJustReversed(edge_dict, pathdump):
            alreadydone = []
            sortkeys = list(edge_dict.keys())
            sortkeys.sort()
        
            for i in range(len(sortkeys)):
                itup = sortkeys[i]
                ipathind = edge_dict[itup]
                ipath = pathdump[ipathind]
                copyipath = copy(ipath)
                copyipath.sort()
                if tuple(copyipath) in alreadydone:
                    return True
                alreadydone.append(tuple(copyipath))
            return False

        # used to weed out deformed/degenerate cube contenders
        def bSomeInterveningPathNodesAreNeighbors(edge_dict, pathdump):
            allintervening = []
            sortkeys = list(edge_dict.keys())
            sortkeys.sort()
        
            for i in range(len(sortkeys)):
                itup = sortkeys[i]
                ipathind = edge_dict[itup]
                ipath = pathdump[ipathind]

                if len(ipath) == 2:
                    continue
                iinterven = ipath[1:-1]
                iinterven_nbrs = []
                for inode in iinterven:
                    iinterven_nbrs.extend(self.NodeVec[inode].Neighbors)

                for j in range(len(sortkeys)):
                    if j == i:
                        continue
                    jtup = sortkeys[j]
                    jpathind = edge_dict[jtup]
                    jpath = pathdump[jpathind]
                    if len(jpath) == 2:
                        continue
                    jjinterven = jpath[1:-1]
                    intersected = list(set(iinterven_nbrs).intersection(jjinterven))
                    if len(intersected) > 0:
                        return True
            return False

        # used to weed out deformed/degenerate cube contenders
        def bNeighborExclusivity(edge_dict, pathdump):
            """
            No node along any edge can have a neighbor that is a node along some other edge
            (note: this does not catch identical edges for supposedly distinct corners if the two edge-paths are reverses of
            one another, but those situations are filtered out elsewhere)
            """
            def AllButThisEdge(edge_dict, excludedkey, pathdump):
                myretval = []
                for key, pathind in edge_dict.items():
                    if key == excludedkey:
                        continue
                    myretval.extend(pathdump[pathind])
                return myretval

            for key, pathind in edge_dict.items():
                thispath = pathdump[pathind]
                if len(thispath) == 2:
                    continue
                
                othernodes = AllButThisEdge(edge_dict, key, pathdump)
                for inode in thispath[1:-1]:
                    nbrs = self.NodeVec[inode].Neighbors
                    for inbrnode in nbrs:
                        if inbrnode in othernodes and not(inbrnode in thispath):
                            #import pdb; pdb.set_trace()
                            return False
            return True


            allintervening = []
            sortkeys = list(edge_dict.keys())
            sortkeys.sort()
        
            for i in range(len(sortkeys)):
                itup = sortkeys[i]
                ipathind = edge_dict[itup]
                ipath = pathdump[ipathind]

                if len(ipath) == 2:
                    continue
                iinterven = ipath[1:-1]
                iinterven_nbrs = []
                for inode in iinterven:
                    iinterven_nbrs.extend(self.NodeVec[inode].Neighbors)

                for j in range(len(sortkeys)):
                    if j == i:
                        continue
                    jtup = sortkeys[j]
                    jpathind = edge_dict[jtup]
                    jpath = pathdump[jpathind]
                    if len(jpath) == 2:
                        continue
                    jjinterven = jpath[1:-1]
                    intersected = list(set(iinterven_nbrs).intersection(jjinterven))
                    if len(intersected) > 0:
                        return True
            return False





        # used to prune away potential edges as the cube is being consolidated
        def PruneASingleContender(SubCoord_path, pathdump):
            # prune the segments with low sum1 that don't get used in paths associated with higher sum1
            keys = list(SubCoord_path.keys())
            ones_keys =[]
            for ikey in keys:
                ones_keys.append((np.sum(ikey), ikey))


            ones_keys.sort()
            maxkey = ones_keys[-1][1]
            ones_keys.reverse()

            for sum1, key in ones_keys:
                if sum1 == 1:
                    continue # we've already checked that the axis nodes are distinct.
                subkeys = SubCoord_path[key][0][0].keys()

                for isubkey in subkeys:
                    dellist = []

                    try:
                        for iscenario in range(len(SubCoord_path[isubkey])):
                            if len(SubCoord_path[isubkey]) <= 1:
                                continue                        
                            thisendpt = SubCoord_path[isubkey][iscenario][-1]

                            pathind = SubCoord_path[key][0][0][isubkey] # the key directory is higher up than subkey, and so will already be pruned to a single member by the time the routine gets here
                            if thisendpt !=  pathdump[pathind][0]:
                                dellist.append(iscenario)  
                    except:
                        import pdb; pdb.set_trace()                   
                    
                    if len(dellist) > 0:
                        if len(dellist) == len(SubCoord_path[isubkey]):
                            print("none left")
                            import pdb; pdb.set_trace()
                        dellist.reverse()
                        for idel in dellist:
                            del SubCoord_path[isubkey][idel]

            return SubCoord_path


        def bGenerationCheck(sub_coord_path, pathdumpgeneration):
            for key, listtuples in sub_coord_path.items():
                subdict = listtuples[0][0]
                sum1 = np.sum(key) - 1
                for pathind in subdict.values():
                    if pathdumpgeneration[pathind] != sum1:
                        return False
            return True

        # used to weed out deformed/degenerate cube contenders
        def bAllPathsDistinct(SubCoord_path):

            # now prune the segments with low sum1 that don't get used in paths associated with higher sum1
            keys = list(SubCoord_path.keys())
            keys.sort()

            edgepaths = []
            

            for ikey in keys:
                xdict = SubCoord_path[ikey][0][0]
                for ival in xdict.values():
                    edgepaths.append(ival)

            # check that all path indices are distinct
            if len(edgepaths) != len(set(edgepaths)):
                return False

            # check that all "inside" edge nodes belong only to one edge
            intervening_nodes = []
            for ipathindex in edgepaths:
                ipath = pathdump[ipathindex]
                if len(ipath) > 0:
                    intervening_nodes.append(ipath[1:-1]) # note this ignores 2-node paths
            if len(intervening_nodes) > 1:
                intersected = list(set(intervening_nodes[0]).intersection(*intervening_nodes[1:]))
            
                return len(intersected) == 0
            
            return True




        def GetCorners(sub_coord_path, newcoordnode, newnodecoord, pathdump, pathdumpgeneration):
            """
            Each corer of  cube connects to D other corners.
            Return corner nodes along with the 2nd (nearest  neighbor) node
            of the path connecting that corner node to those D other corners

            """




            def GetNeighborCorners(inode, sub_coord_path, newnodecoord, pathdump, pathdumpgeneration, coordbasis):
                """
                return the nodes of the D adjoining corners, their coords,
                the paths connecting them to the nde

                """
                icoord = newnodecoord[inode]
                icoordnumberofones = np.sum([j==1 for j in icoord])

                adjcoordlist = []
                adjnodelist = []
                for i in range(self.NDim):
                    adjcoord = np.array(icoord).astype("int") + coordbasis[i]
                    adjcoord = tuple([ j % 2 for j in adjcoord])
                    adjcoordlist.append(adjcoord) 
                    adjnodelist.append(newcoordnode[adjcoord])
                
                adjnodepaths = {}
                for iadj, adjcornernode in enumerate(adjnodelist):
                    thiscoord = adjcoordlist[iadj]
                    thisnumberofones = np.sum([j==1 for j in thiscoord])
                    adjnodepaths[adjcornernode] = []
                    for key in sub_coord_path.keys():
                        for icoord,ipath in sub_coord_path[key][0][0].items():

                            if pathdump[ipath][0] == inode and pathdump[ipath][-1] == adjnodelist[iadj]:
                                if not(pathdump[ipath][1] in adjnodepaths[adjcornernode]):
                                    adjnodepaths[adjcornernode].append(pathdump[ipath][1])
                                    if pathdumpgeneration[ipath] != thisnumberofones-1:
                                        print("Wrong generation - this coord", icoord, "has", thisnumberofones, "components equal to one, which indicates how many steps it is away from origin.") 
                                        import pdb; pdb.set_trace()
                            # could also be connected in the opposite orientation
                            elif pathdump[ipath][-1] == inode and pathdump[ipath][0] == adjnodelist[iadj]:
                                if not(pathdump[ipath][-2] in adjnodepaths[adjcornernode]):
                                    adjnodepaths[adjcornernode].append(pathdump[ipath][-2])
                                    if pathdumpgeneration[ipath] != icoordnumberofones-1:
                                        print("Wrong generation - this coord", icoord, "has", icoordnumberofones, "components equal to one, which indicates how many steps it is away from origin.") 
                                        import pdb; pdb.set_trace()



                    if len(adjnodepaths[adjcornernode]) != 1:
                        print("PROBLEM -- there should be one and only one path in sub_coord_path[111...] connecting", inode, "to ", adjcornernode, "whereas the number here is", len(adjnodepaths[adjcornernode]), "this cube is ill-formed.", newcoordnode)
                        import pdb; pdb.set_trace()
                        return None
                
                myretval = []
                for adjcornernode in adjnodelist:
                    myretval.append( adjnodepaths[adjcornernode][0])
                return myretval




            coordbasis = []
            for i in range(self.NDim):
                icoord = np.zeros((self.NDim,)).astype("int")
                icoord[i] = 1    
                coordbasis.append(icoord)

            cornerpaths = {}
            for inode in newnodecoord.keys():
                cornerpaths[inode] = GetNeighborCorners(inode, sub_coord_path, newnodecoord, pathdump, pathdumpgeneration, coordbasis)
                if cornerpaths[inode] is None:
                    return None
                else:
                    cornerpaths[inode] = tuple(cornerpaths[inode])



            # now check that the intervening nodes of corner paths are nowhere repeated (later we will also check that the 
            # intervening nodes are additionally not even neighbors of other intervening nodes)
            #if self.NNode >= 7816:
            #    print("If this works, get rid of this line and the preceding one, and tab the rest in the black to the left.")
            interveningnodes = []
            for inode in cornerpaths:
                if len(cornerpaths[inode]) > 2:
                    interveningnodes.extend(cornerpaths[inode][1:-1])
            
            if len(interveningnodes) != len(set(interveningnodes)):
                return None
            
            return cornerpaths





        if self.NodeVec[Id].NRegions != 2 ** self.NDim:
            #print("ERROR --  the corners of the cube do not have 2*D neighbors?")
            return None, None, None         
    
        bUsingFarNeighbor = False # if any of the 2**D edges are farnbrs as opposed to nearest neighbors, do not need to test for interior nodes

        # get the first D coords along the D-dimensional axis
        # startnode, pathdumpindex, nearestneighbor, endnode 


        coord_path = {}
        # note that for any key (i.e. binvec from 0 to 2**D) coord_path can both fork into another coord_path  (if the number of list entries for that key increases)
        # or else it can annihilate  if the number of list entries is found to be zero

        # values are list of n tuples where n is the sumbyones number for that coord:
        # for coords whose node is order 1 (i.e. axis), each 2-tuple contais the pathdump index and endpoint
        # for coords whose node is order 2 (i.e. faces connected to origin), each 3-tuple has two 2-tuples (each corresponding
        # to two pathdump indices, and having the same endpoint/startpoint resp) and a final endpoint -- so reallly 4 + 1 endpoints
        # for coords wose node is order 3 (i.e. cubes connected to origin), each 4-tuple has thee 3-tuple pathdump indices and an endpoint
        # etc.

        # algorithm proceeds forward according to sumbyones. Unless all the "faces" are doable, we abort and consider the axes
        # a failure, thereby winnowing down the possibliities without any complicataed iteration.
        # at the end, we have to delete all the failed tries, thereby leaving us (hopefully) with a single node at each of
        # the 2**D possibilities.

        # start with the axes, i.e. the coords for which sumbyones == 1

        if sortbyones is None:
            sortbyones = [] # sortying by number of "middle" or intervening components (the "ones" is a misnomer except in case where nslices=1)
            for j in range(2**self.NDim):
                thisnode = copy(Id) # always start at the anchor node
                binstr = bin(j)[2:]
                while len(binstr) < self.NDim:
                    binstr = '0' + binstr
                lenbinstr = len(binstr)
                binvec = [int(chr) for chr in binstr]
                sortbyones.append((np.sum(binvec), tuple(binvec)))   
            sortbyones.sort()

        zerocoord = tuple([0 for i in range(self.NDim)])
        for iic, iipath in enumerate(pathdump_index_subset):        
            thiscoord = sortbyones[1 + iic][1]
            coord_path[thiscoord] = [({zerocoord:iipath}, pathdump[iipath][-1])]
            # one element for coords whose sum1==1, but higher ones may have more elements than one so we will mk it a list
            # in this case, too

        for sum1, binvec in sortbyones[1:]: # we already installed the first one
            if sum1 == 0:
                continue
            if sum1 == 1:
                continue     
                # we have ve already added the appropriate far neigbors to  coordnode and nodecoord so we don't need the following block                

            if sum1 >= 2:

                tobeintersected = {}

                thesecoords = self.DecomposeBinaryIntoSmallerDimensionalBinaries(binvec)



                bItIsEmpty  = False
                for itup in thesecoords:
                    #import pdb; pdb.set_trace( )
                    if len(coord_path[itup]) == 0:
                        coord_path[tuple(binvec)] = []
                        bItIsEmpty = True
                        break
                    
                    # new version    
                    for ip in coord_path[itup]:
                        ipextrax = self.AddSegment(ip, pathdump, pathdumpgeneration)
                        if not(ipextrax is None):
                            tobeintersected[itup] = ipextrax





                if bItIsEmpty:
                    continue


                if np.min([ len(i)  for i in tobeintersected.values()]) == 0:
                    return None, None, None

                coord_path[binvec] = self.SimpleIntersection(tobeintersected)

                if len(coord_path[binvec]) == 0:
                    return None, None, None
              
                if len(coord_path[binvec]) > 1:
                    # we allow this for now, and weed out any malformed cubes later
                    #print("More than one possible?", nbrsubset)
                    #import pdb; pdb.set_trace()
                    pass


            
        # in general, there will still be several contenders in coord_path, one for each member of the list at "opposite ndoe"
        # whose binvec is tuple([1 for i in range(self.NDim)])
        # the presence of multiple contenders may be in and of itself problematic enogh to warrant excluding this particular
        # choice of potential hypercube, but I prefer, if possible, to winnow the acceptable number down to 1 (or 0) on other grounds

        
        contenders = []
        NContenders = len(coord_path[self.OnesCoord])
        good_coordnode_nodecoord_pairs = []
        node0 = tuple([0 for i in range(self.NDim)])

        





        

        for icontender in range(NContenders):
            bStllAGoodContender = True
            sub_coord_path = copy(coord_path)
            sub_coord_path = copy(coord_path)
            # each subconented will have only one possible set of coords and nodes, because we only allow one element at oppnode
            #sub_coord_path[self.OnesCoord] = [coord_path[self.OnesCoord][icontender]]
            sub_coord_path[self.OnesCoord] = [coord_path[self.OnesCoord][icontender]]



            if not(bGenerationCheck(sub_coord_path, pathdumpgeneration)):
                print("Failed generations? If this never hapens, get rid of the check")
                import pdb; pdb.set_trace()
                bGenerationCheck(sub_coord_path)
            #import pdb; pdb.set_trace()
            # NOTE: this routine gets rid of unused lower-order possible scenarios (e.g. multiple [100/010] intersections in the case of 
            # searching for a [110] node) and then proceeds; it may turn out to be the case that the very existence of such intersections
            # is aalready a sign that we should abort the attempt to find a well-formed cube with the axes used.
 
            sub_coord_path = PruneASingleContender(sub_coord_path, pathdump)   
            
            edge_dict = GetEdgeDict(sub_coord_path)

            # if same path serves as more than one edge, abort
            if len(edge_dict.values()) != len(set(edge_dict.values())):
                continue


            if bSomeEdgesAreSameJustReversed(edge_dict, pathdump):
                continue

            if bSomeInterveningPathNodesAreNeighbors(edge_dict, pathdump):
                continue
                # Why do we need this check at all? In 3 dimensions, you can have an empty cube surrounded by adjacent divided cubes on all
                # but one side (with the remaining side abutting a region that has same lower fractality as the hole within), and the returned
                # cube can be a "sandwich" cube that is shorter along one axis than the others. To forefend that,
                # we check that the intervening nodes (i.e. that are not corners of the cube) are not nearest neighbors
                # of any other such intervening node (this cannot happen in 2d because a sandwich would have some corners
                # with less than 2D neighbors)

            if not(bNeighborExclusivity(edge_dict, pathdump)):
                # this ensures that only the corners of the cube have neighbors that are in more than one edge
                # this ensures the edges "expand". so to speak. away from each other as they extend from any corner
                # as opposed to collapsing into what will form some weird degenerate cute
                continue

            if not(bAllPathsDistinct(sub_coord_path)):
                print("Not all distinct:", coord_path)
                import pdb; pdb.set_trace()
                continue

            nremaining = []
            for key,val in sub_coord_path.items():
                nremaining.append(len(val))
            if np.min(nremaining) == 0:
                import pdb; pdb.set_trace()
                continue

            newcoordnode = {node0:Id}
            newnodecoord = {Id:node0}
            
            for icoord in coord_path.keys():
                if len(sub_coord_path[icoord]) > 1:
                    print("Why was this not pruned?")
                    import pdb; pdb.set_trace()
                    x = PruneASingleContender(sub_coord_path, pathdump)  
                thisendpt = sub_coord_path[icoord][0][-1]

                newcoordnode[icoord] = thisendpt
                newnodecoord[thisendpt] = icoord
            
            if len(newnodecoord) != 2**self.NDim: # this happens (rarely) when the corner nodes of the cube are nondistinct, so that one or more entries get overwritten resulting in less than 2**D entries
                    continue

            # now, check that all the nodes connecting edges of the cube are distinct          
            # remember, all keys have a list of 1 element at this point  
            allpaths = []     
            
            for icoord, ilist in sub_coord_path.items():
                allpaths.extend(list(ilist[0][0].values()))
            bAllAreNear = True
            allnodes = []
            for ip in list(set(allpaths)):
                thispath = pathdump[ip]
                if len(thispath) != 2:
                    bAllAreNear = False
                allnodes.append(thispath)


            if len(allpaths) != len(set(allpaths)):
                import pdb; pdb.set_trace()
                bStllAGoodContender = False

            if bStllAGoodContender:
                # check that all the edges have no intervening nodes that are shared.
                for iipath, ipath in enumerate(allnodes):
                    for jjpath in range(iipath+1, len(allnodes)):
                        jpath = allnodes[jjpath]            
                        intersected = list(set(ipath).intersection(jpath))
                        # note we do not have to worry about the nodes of the cube being an in-between node of any path since cube nodes
                        # must have 2D elements and all in-between nodes have less than 2D neighbors

                        for i in intersected:
                            if i in newnodecoord.keys():
                                intersected.remove(i)
                        
                        if len(intersected) != 0:
                            #return None, None
                            bStllAGoodContender = False

            # test for interior nodes and make sure the space we wish to divide is empty
            # note, for some "crazy" nonlinear edges that have not yet been excluded, 
            # the region can have one or more short edges (i.e. with nearest neighbors)
            # and still have interior nodes, so we must always 
            # to see if there's an interior node and thereby exclude such crazy nonlinear edges
            # the only exception is if all edges are nearest-neighor-related in which case the
            # interior node is just the one corresponding to self.OnesCoord (i.e. the test fails)


         
            if bStllAGoodContender and not(bAllAreNear):
                   
                corners = GetCorners(sub_coord_path, newcoordnode, newnodecoord, pathdump, pathdumpgeneration)
                if corners is None:
                    return None, None, None
                # if same node serves as more than one corner, abort
                if len(list(corners.values())) != len(set(corners.values())):
                    return None, None, None
                # Either all the sourrounding cubes are divided, or else the cube itself is divided, so test if there is an interior node         
                if not(corners is None):
                    for inode, tupaxissubset in corners.items():
                        if not((inode, tupaxissubset) in TestForInteriorPoints):
                            TestForInteriorPoints[(inode, tupaxissubset)] = self.bBruteFindInterior(inode, tupaxissubset, sortbyones)                    
                        bIsThereInteriorPoint = TestForInteriorPoints[(inode, tupaxissubset)]
                        if bIsThereInteriorPoint:
                            bStllAGoodContender = False 
                            break
                            

            if bStllAGoodContender:
                exportable_sub_coord_path = ChangePathIndicesToNodeLists(sub_coord_path, pathdump)
                good_coordnode_nodecoord_pairs.append((newcoordnode, newnodecoord, exportable_sub_coord_path))
        
        if len(good_coordnode_nodecoord_pairs) > 1:
            print("multiple possibilities for coordnode and nodecoord? Will return None, None but find out why this is happening. ")
            import pdb; pdb.set_trace()
            return None, None, None
        elif len(good_coordnode_nodecoord_pairs) == 1:
            (newcoordnode, newnodecoord, exportable_sub_coord_path) = good_coordnode_nodecoord_pairs[0]
            return newcoordnode, newnodecoord, exportable_sub_coord_path
        else:
            #print("no good cubes available at", Id, self.NodeVec[Id].Coords, nbrsubset)
            #import pdb; pdb.set_trace()
            return None, None, None
                
        






        
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
        # only doing this for 2 and 3-dimensional case
        if self.NDim == 2:
            
        
            for inode in range(self.NodesLastHistogramUpdate, len(self.NodeVec)):
                
                x,y = self.NodeVec[inode].Coords 
                self.Hist[int(np.floor(x * self.HistScale)), int(np.floor(y * self.HistScale))] += 1

            self.NodesLastHistogramUpdate = copy(self.NNode)

    def HistNeighborCorrel(self, len=1, leny=0):
        # only doing this for 2 and 3-dimensional case
        if self.NDim != 2:

            lenhist = self.Histogram.shape[0]
            xvec = []
            yvec = []
            # 2-d
            for ix in range(self.Hist.shape[0]):
                for iy  in range(self.Hist.shape[1]):
                    
                    xdisplaced = (ix + len) % lenhist # can also displace in the negative direction, or in the
                    ydisplaced = (iy + leny) % lenhist
                    xvec.append(self.Histogram[ix,iy])
                    yvec.append(self.Histogram[ix+xdisplaced, iy+ydisplaced])
            
            return np.corrcoef(np.array(xvec), np.array(yvec))

            

    def HistogramPercentile(self,):
        if self.NDim != 2:
            return
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
        bPrintStuff = False

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

    def l2(self, x,y):
        # allows one to do a norm() across the toroid boundaries
        dvec = x - y
        for i in range(self.NDim):
            thisdist = x[i] - y[i]
            if thisdist >= self.TorLen - 1:
                dvec[i] -= self.TorLen
            if thisdist <= -self.TorLen + 1:
                dvec[i] += self.TorLen

        return np.sqrt(np.sum(dvec * dvec))

    def l2node(self, i, j):
        x = np.array(self.NodeVec[i].Coords)
        y = np.array(self.NodeVec[j].Coords)
        return self.l2(x, y)


    def bCubeIntegrityCheck(self, coordnode):
        """
        Simple check to make sure all the corners of cube are separated by appropriate distances
        """

        origincoord = tuple([0] * self.NDim)
        originnode = coordnode[origincoord]
        originvec = np.array(self.NodeVec[originnode].Coords)

        basedist = 0
        epsilon = 1.0e-6
        onemepsilon = 1.0 - epsilon
        onepepsilon = 1.0 + epsilon
        for i in range(self.NDim):
            thiscoord = list(origincoord)
            thiscoord[i] = 1
            thisvec = self.NodeVec[coordnode[tuple(thiscoord)]].Coords
            thisdist = self.l2(np.array(thisvec), originvec)

            if basedist == 0:
                basedist = copy(thisdist)
            else:
                if basedist < onemepsilon * thisdist or basedist > onepepsilon:
                    return False     

            
        
        for coord,node in coordnode.items():
            sumones = np.sum(coord)
            if sumones <= 1:
                continue # already did this
            thisvec = self.NodeVec[coordnode[coord]].Coords
            thisdist = self.l2(thisvec, originvec)
            target = np.sqrt(sumones) * basedist
            if thisdist < onemepsilon * target or thisdist > onepepsilon * target:
                return False
        
        return True


        origincoord = coordnode(tuple([0] * self.NDim))
    
    def bIsMixedFractality(self, coord_path):
        for coord, list_tups in coord_path.items():
            for isubkey, ipath in list_tups.items():
                if len(ipath) > 2:
                    return True
        return False


    # used for debugging and inspection -- this is just the loop from  expandmanyrandomly
    def ExpandOneNode(self, inode):
        
        if self.MaxDist(inode) == 1.0 and self.bStop:
            print("why not updating")
            import pdb; pdb.set_trace()
        import pdb; pdb.set_trace()
        coordnode_nodecoord_coordpath_list = self.ReturnAllWellFormedCubes(inode)

        if self.NodeVec[inode].NRegions != 2 ** self.NDim:
            return


        IsMixedFractalityArr = []
        MixedFractalityScenarios = []
        for iscenario in range(len(coordnode_nodecoord_coordpath_list)):
            icoordnode, inodecoord, icoord_path = coordnode_nodecoord_coordpath_list[iscenario]
            thisismixed = self.bIsMixedFractality(icoord_path) 
            IsMixedFractalityArr.append(thisismixed)
            if thisismixed:
                MixedFractalityScenarios.append((icoordnode, inodecoord, icoord_path))
        
        # if there are any mixed-fractality cubes, we'll always preferentially divide them rather than  
        if np.sum(IsMixedFractalityArr) > 0:
            coordnode_nodecoord_coordpath_list = copy(MixedFractalityScenarios)
            IsMixedFractalityArr = [ True ] * len(MixedFractalityScenarios)

        if len(coordnode_nodecoord_coordpath_list) == 0:
            # if this happens, we will look in the vicinity, and if we find a 
            # mixed scenario we will use that instead. We shall search in
            # stages goint out NDim steps (totally ad hoc choice)

            coordnode_nodecoord_coordpath_list = self.SphereSweep(inode)
            if coordnode_nodecoord_coordpath_list is None:
                print("Could find nothing in the vicinity?")
                import pdb; pdb.set_trace()
            else:
                print("Found alternate node with allowable divisions.")



        rn.shuffle(coordnode_nodecoord_coordpath_list)
        dellist = []
        for iielement, ielement in enumerate(coordnode_nodecoord_coordpath_list):
            subcoordnode, subnodecoord, subcoordpath = ielement
            #if not(self.bIsNonUniform(ielement[0])):
            if not(self.bIsNonUniform(subcoordnode)):
                #import pdb; pdb.set_trace()
                if rn.random() > FreshProb:
                    dellist.append(iielement) 
                else:
                    #print("Passed the uniformity test")
                    pass
            else:
                #print("mixed")
                pass

        if len(dellist) > 0:
            dellist.reverse()
            for idel in dellist:
                del coordnode_nodecoord_coordpath_list[i]
            
        if len(coordnode_nodecoord_coordpath_list) == 0:
            #print("No good coords at point ", inode, " for any possible set of axes there")
            return 


            # it has already been shuffled, so you can just return the last one.
        coords, noodes, coordpath = coordnode_nodecoord_coordpath_list[0]

            
            
            



        print("NDIV", ndivisions, inode)



        #coords, noodes, iax = self.PickADivisibleCube(inode, allaxes)


        #test for closeness of coords (just a debugging/monoitoring thing)
        bTestForSuspiciouslyCloseNodes = True

        if bTestForSuspiciouslyCloseNodes:

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


        if not(self.bCubeIntegrityCheck(coords)):
            print("BAD CUBE BEING DIVIDED -- what's going on?")
            import pdb; pdb.set_trace()
            self.bCubeIntegrityCheck(coords)


        
        self.BruteDivideCube(coords, noodes, coordpath)



    # old routine no longer used
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

        for iinode in sorted(nodecoord.keys()):

                #import pdb; pdb.set_trace()
      
                inode = self.NodeVec[iinode]

                if inode.Parity != 1:
                    continue
                maxneighbor = max(inode.Neighbors)
                if maxneighbor >= len(coordnode):
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
                 
                print("OK", inode.Id, len(CoordNode))
                self.DivideCube(CoordNode, NodeCoord)


    # deprecated -- no longer used
    def ExpandManyRandomly(self, Prob, NRuns):
        #print("dont forget to randomize nslices")

        #print("Assign sierpinski level to try and even out the ")
        global FreshProb
        global MaxNodes

        self.ClearHistogram()
        ndivide = 0

        nchunk = 1
        nlastchunk = 0

        histmean = [0]
        histstd = [0]

        startnnode = copy(self.NNode)

        lastprintscreen = 0
        ii = 0
        ndivisions = 0
        bDone = False
        for i in range(NRuns):
            if bDone:
                break

            if self.NNode < self.TorLen ** self.NDim:
                continue
                    

            if i % 10000 == 0:
                print(i)


            toprun = copy(self.NNode)
            bottomrun = 0

            ishuf = list(np.arange(0,self.NNode).astype("int"))
            rn.shuffle(ishuf)
            #    bottomrun = int(rn.random()*(self.TorLen ** self.NDim))


            for inode in ishuf: #range(bottomrun, toprun):
                
                if self.MaxDist(inode) == 1.0 and self.bStop:
                    print("why not updating")
                    import pdb; pdb.set_trace()
                coordnode_nodecoord_coordpath_list = self.ReturnAllWellFormedCubes(inode)

                if self.NodeVec[inode].NRegions != 2 ** self.NDim:
                    continue


                IsMixedFractalityArr = []
                MixedFractalityScenarios = []
                for iscenario in range(len(coordnode_nodecoord_coordpath_list)):
                    icoordnode, inodecoord, icoord_path = coordnode_nodecoord_coordpath_list[iscenario]
                    thisismixed = self.bIsMixedFractality(icoord_path) 
                    IsMixedFractalityArr.append(thisismixed)
                    if thisismixed:
                        MixedFractalityScenarios.append((icoordnode, inodecoord, icoord_path))
                
                # if there are any mixed-fractality cubes, we'll always preferentially divide them rather than  
                if np.sum(IsMixedFractalityArr) > 0:
                    coordnode_nodecoord_coordpath_list = copy(MixedFractalityScenarios)
                    IsMixedFractalityArr = [ True ] * len(MixedFractalityScenarios)


                                            
                if len(coordnode_nodecoord_coordpath_list) == 0:
                    # if this happens, we will look in the vicinity, and if we find a 
                    # mixed scenario we will use that instead. We shall search in
                    # stages goint out NDim steps (totally ad hoc choice)

                    coordnode_nodecoord_coordpath_list = self.SphereSweep(inode)
                    if coordnode_nodecoord_coordpath_list is None:
                        print("Could find nothing in the vicinity?")
                        import pdb; pdb.set_trace()
                    else:
                        #print("Found alternate node with allowable divisions.")
                        print("Found alternative choice -- no longe using ", inode)



                rn.shuffle(coordnode_nodecoord_coordpath_list)
                dellist = []
                for iielement, ielement in enumerate(coordnode_nodecoord_coordpath_list):
                    subcoordnode, subnodecoord, subcoordpath = ielement
                    #if not(self.bIsNonUniform(ielement[0])):
                    if not(self.bIsNonUniform(subcoordnode)):
                        #import pdb; pdb.set_trace()
                        if rn.random() > FreshProb:
                            dellist.append(iielement) 
                        else:
                            #print("Passed the uniformity test")
                            pass
                    else:
                        #print("mixed")
                        pass

                if len(dellist) > 0:
                    dellist.reverse()
                    for idel in dellist:
                        del coordnode_nodecoord_coordpath_list[i]
                    
                if len(coordnode_nodecoord_coordpath_list) == 0:
                    #print("No good coords at point ", inode, " for any possible set of axes there")
                    continue


                 # it has already been shuffled, so you can just return the last one.
                coords, noodes, coordpath = coordnode_nodecoord_coordpath_list[0]

                    
                 
                   



                print("NDIV", ndivisions, inode)
 


                #coords, noodes, iax = self.PickADivisibleCube(inode, allaxes)


                #test for closeness of coords (just a debugging/monoitoring thing)
                bTestForSuspiciouslyCloseNodes = True

                if bTestForSuspiciouslyCloseNodes:

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


                if not(self.bCubeIntegrityCheck(coords)):
                    print("BAD CUBE BEING DIVIDED -- what's going on?")
                    import pdb; pdb.set_trace()
                    self.bCubeIntegrityCheck(coords)


                
                self.BruteDivideCube(coords, noodes, coordpath)

                ndivisions += 1
                if self.NNode > MaxNodes:
                    print("reached maxnode upper limit", self.NNode)
                    bDone = True
                    break
                if ii - lastprintscreen > 2048:
                    #self.ScatterPlot()
                    lastprintscreen = copy(ii)
                ii  += 1


                ndivide += 1
                if ndivisions % 1 == 0:


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
                        self.UpdateNeighborDistanceCount()
                        
                        if ndivide % 4 == 0:
                            self.ExportPickle(self,"filedump.pkl")

                        histmean.append( np.mean(self.Hist) )
                        histstd.append( np.std(self.Hist) )


                        nlastchunk = (ndivide // nchunk) * nchunk






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


    # deprecated in favor of ReturnAllWellFormed... 
    def PickADivisibleCube(self, startnodeId, allaxes = None):

        if allaxes is None:
            allaxes = self.ReturnAllAxes(startnodeId)

        resallax = []
        alreadyrejected = []
        alreadyaccepted = []
        for iax in allaxes:

            testiax = copy(iax)
            testiax.sort()
            if tuple(testiax) in alreadyrejected:
                continue
            if tuple(testiax) in alreadyaccepted:
                # not retained, and since all it does is permute the numbers assigned to the new nodes,
                # there's no point in retesting a set of axes if one of their permutations has already passed
                continue
                
            coords, noodes = self.FindCubeCoord(startnodeId, iax)

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


    ### Next few routines are for plotting and locationm which is useful for debugging and visualization
    
    def NodeFromCoords(self, coordtup, err=1.0e-5):
        x0 = np.array((coordtup))
        for i, ix in enumerate(x0):
            if ix < 0:
                x0[i] += self.TorLen


        myretval = []
        for inode in self.NodeVec:
            idist = np.linalg.norm(x0 - np.array(inode.Coords))
            if idist <= err:
                myretval.append(inode.Id)
        return myretval
    
    # return the axis nodes obtained by moving in the POSITIVE direction away from Id
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
            if inode.NRegions != 2 ** self.NDim:
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

    
    # INCOMPLETE
    def PlotCoordPath(self, coordpath, dist = 0.1,  wLabels=True):
        def GetAllNodes(coordpath, wLabels):
            myretval = []
            for upkey, subdict in coordpath.items():
                for subkey, listnode in subdict.items():
                    myretval.extend(listnode)
            if wLabels:
                return myretval, [str(i) for i in myretval]
            else:
                return myretval
        if wLabels:
            listnodes, labels = GetAllNodes(coordpath, True)
        else:
            listnodes = GetAllNodes(coordpath, False)
       
        if self.NDim == 3:
            """
            xmin = np.min([self.NodeVec[i].Coords[0] for i in listnodes])
            xmax = np.max([self.NodeVec[i].Coords[0] for i in listnodes])
            ymin = np.min([self.NodeVec[i].Coords[1] for i in listnodes])
            ymax = np.max([self.NodeVec[i].Coords[1] for i in listnodes])
            zmin = np.min([self.NodeVec[i].Coords[2] for i in listnodes])
            zmax = np.max([self.NodeVec[i].Coords[2] for i in listnodes])
            """
            x = []
            y = []
            z = []
            print(listnodes)
            for inode in listnodes:
                crd = self.NodeVec[inode].Coords
                x.append(crd[0])
                y.append(crd[1])
                z.append(crd[2])

            import pdb; pdb.set_trace()
            from mpl_toolkits.mplot3d import Axes3D
            s = 10
            fig = plt.figure()
            ax = fig.add_subplot(projection='3d')

            for i in range(len(x)): #plot each point + it's index as text above
                ax.scatter(x[i],y[i],z[i],color='b') 
                if wLabels:
                    ax.text(x[i],y[i],z[i],  '%s' % (labels[i]), size=s, zorder=1, color='k') 

        return 
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
    
    def PlotSegment(self, A, B, ax, c='b', linewidth=3):
        if self.NDim == 2:
            x = [A[0], B[0]]
            y = [A[1], B[1]]   
            ax.plot(x,y)     
        if self.NDim == 3:
            x = [A[0], B[0]]
            y = [A[1], B[1]]   
            z = [A[2], B[2]]   

            from mpl_toolkits.mplot3d import Axes3D


            # Create a figure and 3D axis
            s = 10
            fig = plt.figure()
            ax = fig.add_subplot(projection='3d')
            ax.plot(x, y, z)    
            

    def FindNodes(x, err = 1.0e-9):
        myretval = []
        for i in range(len(self.NodeVec)):
            dist = self.l2(np.array(self.NodeVec[i].Coords), x)
            if dist <= err:
                myretval.append(i)
        return myretval


    def ExportPickle(self, ggrid, fname):
        import pickle
        fp = open(fname, "wb")
        pickle.dump(ggrid, fp)


    
    def Distributions(self,inode):
        for inod in ggrid.NodeVec:
            inod.Amplitude = 0
            inod.PrevAmplitude = 0

        #import pdb; pdb.set_trace()

        self.NodeVec[inode].PrevAmplitude = 1.0
        self.NodeVec[inode].Amplitude = 1.0

        print("Heat", ggrid.NNode)            
        self.HeatEq(250, inode, True)

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
        tvec = tuple(vec)
        if not tvec in myret:
            myret.append(tvec)
    return myret

def AllUniqueCombosNoReplacement(N, D):
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

        vec.sort()
        tvec = tuple(vec)
        if not(tvec) in  myret:
            myret.append(tvec)
    return myret


def NodeCountToHypercubeCount(HowManyNewNeighbors, NDim):
    """
    For any node in a newly divided grid, the number of neighbors THAT ARE PART OF THAT NEW GRID
    uniquely specifies the number of cubes (which can range, in powers of 2, from 1 to D**2) that new
    grid has taken up (out of the 2**D available).

    Eventually, as a node progresses from newly created entity, to "maturity -- i.e. as subdivisions keep happening, 
    it will attain 2**D distinct cubes adoining it, with
    each cube being part of a different subdivided cube, but those intermediary
    occupied regions will specify the kind of edge the node is on. (Note that nodes that are completely interior to
    the newly created cube will have 2**D separate hypercubes from the beginning, whereas those on an edge, or face -- or some
    higher hypercube -- will have successively larger numbers of hypercubes at the start)

    Note that even if you use two numbers -- i.e. NRegions and NNeighbors -- you will not be able to distinguish all edge configurations.
    For examples consider 3 cubes arranged like an L, compared with a single cube atop a plane. Buth will comprise six out of 8 available
    regions (in 3D) and both will have maxed out NNeighbor numbers of 6. In this case, the region growth will distinguish the two, but
    there can be two different cube-plus-plane scenarios (for which the remainder region is situated in some other quadraant) for which
    even the region growth number won't work. Hopefully, such scenarios will be disqualified by the fact that we cannot have multiple
    paths with same begin/end nodes. Also, even though the available tests can't distinguish the types of edges they are, the NRegions count
    is by itself sufficient to show that each of those scenarios is a kind of edge (which is more than can be said of NNeighbor), 
    and that may be good enough.
    """

    return 2 ** (HowManyNewNeighbors - NDim)

def HypercubeCountToNodeCount(HowManyHyperCubes, NDim):
    """ Inverse of the previous function """
    return NDim + int(np.log2(HowManyHyperCubes))



# a side-routine used for debugging; may be ignored
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
    import pickle


    opts, args = ReadParams()

    if opts.PostAnalysis:
        #import pdb; pdb.set_trace()
        f = open('config3d_75984.pkl', 'rb')
        ggrid = pickle.load(f)
        f.close()

       
        inode = 15603
        coordnode_nodecoord_coordpath_list = ggrid.SphereSweep(inode)
        import pdb; pdb.set_trace()
        


        if False:
            f = open('blah3.pkl', 'rb')
            coords, nodes, coordpath = pickle.load(f)
            f.close()
            import pdb; pdb.set_trace()
            ggrid.PlotCoordPath(coordpath)
    else:
        

        if opts.dimension == 3:
            from mpl_toolkits.mplot3d import Axes3D
        rn.seed(opts.seed)

        #import pdb; pdb.set_trace()
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
        # after cubenode is created, divide the cubes at the following nodes (always choose axes in the pos. direction):

        100, 110, 011 (you can actually omit the first in that list) and then, try to divide the cube at 111

        node_100 = NodeFromCoords((1,0,0))
        node_110 = NodeFromCoords((1,1,0))
        node_011 = NodeFromCoords((0,1,1))

        axisnodes_100 = GetAxisNodes(node_100)
        axisnodes_110 = GetAxisNodes(node_110)
        axisnodes_011 = GetAxisNodes(node_011)

            
        cubes_100 = self.ReturnAllWellFormedCubes(node_100)
        cubes_110 = self.ReturnAllWellFormedCubes(node_110)
        cubes_011 = self.ReturnAllWellFormedCubes(node_011)
        """
        bCreateProblemScenario = False
        if bCreateProblemScenario:
            
            node_100 = ggrid.NodeFromCoords((1,0,0))[0]
            node_110 = ggrid.NodeFromCoords((1,1,0))[0]
            node_010 = ggrid.NodeFromCoords((0,1,0))[0]

            axisnodes_100 = ggrid.GetAxisNodes(node_100)
            axisnodes_110 = ggrid.GetAxisNodes(node_110)
            axisnodes_010 = ggrid.GetAxisNodes(node_010)
 
            #import pdb; pdb.set_trace()
            cubes_100 = ggrid.ReturnAllWellFormedCubes(node_100, None, axisnodes_100)
            coordnode, nodecoord, coord_path = cubes_100[0]
            ggrid.BruteDivideCube(coordnode, nodecoord, coord_path)

            # this next one is not necessary
            cubes_110 = ggrid.ReturnAllWellFormedCubes(node_110, None, axisnodes_110)
            coordnode, nodecoord, coord_path = cubes_110[0]
            ggrid.BruteDivideCube(coordnode, nodecoord, coord_path)

            
            cubes_010 = ggrid.ReturnAllWellFormedCubes(node_010, None, axisnodes_010)
            coordnode, nodecoord, coord_path = cubes_010[0]
            ggrid.BruteDivideCube(coordnode, nodecoord, coord_path)

            import pdb; pdb.set_trace()
            node_000 = ggrid.NodeFromCoords((0,0,0))[0]
            axisnodes_000 = ggrid.GetAxisNodes(node_000)
            cubes_000 = ggrid.ReturnAllWellFormedCubes(node_000, None, axisnodes_000)
            coordnode, nodecoord, coord_path = cubes_000[0]
            ggrid.BruteDivideCube(coordnode, nodecoord, coord_path)

        if False:        
            list_of_coordnode_nodecoord_tuples_0 = ggrid.ReturnAllWellFormedCubes(0)
            (coordnode, nodecoord, coord_path) = list_of_coordnode_nodecoord_tuples_0[1]
            ggrid.BruteDivideCube(coordnode, nodecoord, coord_path)
            #import pdb; pdb.set_trace()
            #ggrid.ScatterPlotLabeled3((-0.1,3.1),(-0.1,3.1),(-0.1,3.1))

        



            startnode00 = 0
            # this uses "old" routines, but they are still useful for construction
            allaxes = ggrid.ReturnAllAxes(startnode00)
            
            targetax = allaxes[4] # should be [32,1]

            (coordnode, nodecoord) = ggrid.FindCubeCoord(0, targetax)
            
            ggrid.DivideCube(coordnode, nodecoord)


            allaxes00 = ggrid.ReturnAllAxes(startnode00)
            
            startnode01 = 1 

            allaxes01 = ggrid.ReturnAllAxes(startnode01)

            
            coords, noodes, iax = ggrid.PickADivisibleCube(startnode01, allaxes01)

            bCoordCheck = True # just for debugging

            if bCoordCheck:
                veclist = []
                sumvec = None
                for iinode, inode in enumerate(noodes.keys()):
                    thisvec = np.array(self.NodeVec[inode].Coords)
                    veclist.append(thisvec)
                    if iinode == 0:
                        sumvec = copy(thisvec)
                        
                    else:
                        sumvec += thisvec
                bSomeAreSame = False
                for j in range(1,ggrid.NDim):
                    if bSomeAreSame:
                        break
                    for i in range(ggrid.NDim):
                        if veclist[j][idim] == veclist[j-1][idim]:
                            bSomeAreSame = True
                            break
                if not(bSomeAreSame):
                    print("The nodes are not aligned along the axes")
                    import pdb; pdb.set_trace()

                meanvec = sumvec / flot(len(veclist))
                distlist = []
                for ivec in veclist:
                    thisdist = np.sqrt(np.sum((ivec - meanvec) * (ivec - meanvec)))
                    distlist.append(thisdist)
                
                if np.max(distlist) - np.min(distlist) > 1.04-3 *  np.mean(distlist):
                    print("these are not squares")




            if not(coords is None):
                ggrid.DivideCube(coords, noodes)




            for iax in allaxes01:
                coords, noodes = ggrid.FindCubeCoord(startnode01, iax)
                #print (iax, coords, noodes)

            nodecoord = {32:(1,0), 33:(1,1), 64:(2,0), 65:(2,1)}
            coordnode = {}
            for key,val in nodecoord.items():
                coordnode[val] = key


        bEnclosedCube = ggrid.NDim == 2

        if bEnclosedCube:


            GuideDict = {}
            NodeCoord = {}
            CoordNode = {}

            Prob = opts.prob
            NRuns = opts.expansionruns

            if True:
                
                ggrid = cubenodeset_t(opts.dimension)            
                
                ggrid.CreateNCube()

                
                FreshProb = 1.0
                i = 0
                print("start", ggrid.NNode); i += 1
                if False:
                
                    (coordnode, nodecoord) = ggrid.FindCubeCoord(32, [33,64])
                    ggrid.BruteDivideCube(coordnode, nodecoord)
                    print(i, ggrid.NNode); i += 1
                    (coordnode, nodecoord) = ggrid.FindCubeCoord(34, [35,66])
                    ggrid.BruteDivideCube(coordnode, nodecoord)
                    print(i, ggrid.NNode); i += 1
                    (coordnode, nodecoord) = ggrid.FindCubeCoord(1, [2,33])
                    ggrid.BruteDivideCube(coordnode, nodecoord)
                    print(i, ggrid.NNode); i += 1
                    (coordnode, nodecoord) = ggrid.FindCubeCoord(65, [66,97])
                    ggrid.BruteDivideCube(coordnode, nodecoord)
                    print(i, ggrid.NNode); i += 1

                    # this next section subdivides each of the 4 divided "petals" 

                    (coordnode, nodecoord) = ggrid.FindCubeCoord(1027, [1031,1028])
                    ggrid.BruteDivideCube(coordnode, nodecoord)

                    (coordnode, nodecoord) = ggrid.FindCubeCoord(1051, [1055,1052])
                    ggrid.BruteDivideCube(coordnode, nodecoord)

                    (coordnode, nodecoord) = ggrid.FindCubeCoord(1039, [1043,1040])
                    ggrid.BruteDivideCube(coordnode, nodecoord)

                    (coordnode, nodecoord) = ggrid.FindCubeCoord(1063, [1067,1064])
                    ggrid.BruteDivideCube(coordnode, nodecoord)



                    list_of_coordnode_nodecoord_tuples_33 = ggrid.ReturnAllWellFormedCubes(33)

                    import pdb; pdb.set_trace()
                    list_of_coordnode_nodecoord_tuples_1087 = ggrid.ReturnAllWellFormedCubes(1087)



                    coordnode, nodecoord  = list_of_coordnode_nodecoord_tuples_33[1]


                    ggrid.BruteDivideCube(coordnode, nodecoord)
                    coordnode, nodecoord  = list_of_coordnode_nodecoord_tuples_33[0]
                    ggrid.BruteDivideCube(coordnode, nodecoord)

                    import pdb; pdb.set_trace()

                    
                    (coordnode, nodecoord) = ggrid.FindCubeCoord(1052, [1089,1053])
                    ggrid.BruteDivideCube(coordnode, nodecoord)
                    (coordnode, nodecoord) = ggrid.FindCubeCoord(1053, [1120,1057])
                    ggrid.BruteDivideCube(coordnode, nodecoord)

                    (coordnode, nodecoord) = ggrid.FindCubeCoord(1120, [1122,1121])                
                    ggrid.BruteDivideCube(coordnode, nodecoord)                
                    (coordnode, nodecoord) = ggrid.FindCubeCoord(1034, [1120,1035])
                    ggrid.BruteDivideCube(coordnode, nodecoord)
                    (coordnode, nodecoord) = ggrid.FindCubeCoord(1031, [1034,1082])
                    ggrid.BruteDivideCube(coordnode, nodecoord)
                    print("same check")
                    import pdb; pdb.set_trace()
                    (coordnode, nodecoord) = ggrid.FindCubeCoord(1036, [1037,1039])
                    ggrid.BruteDivideCube(coordnode, nodecoord)















                    print("same check")
                    import pdb; pdb.set_trace()
                    print(ggrid.NNode)
                    # now do the central cube -- only one new node should be created
                    startnode33 = 33

                    list_of_coordnode_nodecoord_tuples_33 = ggrid.ReturnAllWellFormedCubes(startnode33)

                    import pdb; pdb.set_trace()
                    list_of_coordnode_nodecoord_tuples_1087 = ggrid.ReturnAllWellFormedCubes(1087)




                    coordnode, nodecoord  = list_of_coordnode_nodecoord_tuples_33[0]
                    ggrid.BruteDivideCube(coordnode, nodecoord) 

                    coordnode, nodecoord  = list_of_coordnode_nodecoord_tuples_33[1]
                    ggrid.BruteDivideCube(coordnode, nodecoord) 

                    import pdb; pdb.set_trace()

                    
                    (coordnode, nodecoord) = ggrid.FindCubeCoord(startnode00, [1080, 1077])
                    
                    coordnode, nodecoord, iax = ggrid.PickADivisibleCube(startnode00)
                    if not(coordnode is None):
                        ggrid.DivideCube(coordnode, nodecoord)

                    
                    (coordnode, nodecoord) = ggrid.FindCubeCoord(64, [1029, 96])
                    ggrid.DivideCube(coordnode, nodecoord)
                    (coordnode, nodecoord) = ggrid.FindCubeCoord(66, [1041, 1070])
                    ggrid.DivideCube(coordnode, nodecoord)
                    (coordnode, nodecoord) = ggrid.FindCubeCoord(2, [1058, 3])
                    ggrid.DivideCube(coordnode, nodecoord)

                    
                    (coordnode, nodecoord) = ggrid.FindCubeCoord(1102, [1105, 1103]); ggrid.DivideCube(coordnode, nodecoord, nslices)
                    
                    (coordnode, nodecoord) = ggrid.FindCubeCoord(1041, [1045, 1102]); ggrid.DivideCube(coordnode, nodecoord, nslices)
                    
                    (coordnode, nodecoord) = ggrid.FindCubeCoord(1070, [1102, 1071]); ggrid.DivideCube(coordnode, nodecoord, nslices)
                    (coordnode, nodecoord) = ggrid.FindCubeCoord(1067, [1070, 1068]); ggrid.DivideCube(coordnode, nodecoord, nslices)

                    #import pdb; pdb.set_trace()
                    
                    
                    (coordnode, nodecoord) = ggrid.FindCubeCoord(67, [1108, 68]); ggrid.DivideCube(coordnode, nodecoord, nslices)

                    #import pdb; pdb.set_trace()

                    # ggrid.ScatterPlotLabeled((0.9,4.1),(0.9,4.1))
                    list_of_coordnode_nodecoord_tuples = ggrid.ReturnAllWellFormedCubes(1103)
                    list_of_coordnode_nodecoord_tuples = ggrid.ReturnAllWellFormedCubes(1045)
                    list_of_coordnode_nodecoord_tuples = ggrid.ReturnAllWellFormedCubes(1102)
                




        if bEnclosedCube:


            GuideDict = {}
            NodeCoord = {}
            CoordNode = {}

            Prob = opts.prob
            NRuns = opts.expansionruns

            if False:
                
                ggrid = cubenodeset_t(opts.dimension)            
                
                ggrid.CreateNCube()

                
                FreshProb = 1.0
                i = 0
                print("start", ggrid.NNode); i += 1



                
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
                
                print(ggrid.NNode)


                
                
                # now do the central cube -- only one new node should be created
                startnode33 = 33
                allaxes33 = ggrid.ReturnAllAxes(startnode33)

                
                #(coordnode, nodecoord) = ggrid.FindCubeCoord(startnode33,  [1053,1034] ) # [34,65])
                #ggrid.DivideCube(coordnode, nodecoord, nslices)
                #print(ggrid.NNode)

                (coordnode, nodecoord) = ggrid.FindCubeCoord(1051,  [1055,1052] ) # [34,65])
                ggrid.DivideCube(coordnode, nodecoord, nslices)

                (coordnode, nodecoord) = ggrid.FindCubeCoord(1027,  [1031,1028] ) # [34,65])
                ggrid.DivideCube(coordnode, nodecoord, nslices)

                (coordnode, nodecoord) = ggrid.FindCubeCoord(1063,  [1067,1064] ) # [34,65])
                ggrid.DivideCube(coordnode, nodecoord, nslices)


                (coordnode, nodecoord) = ggrid.FindCubeCoord(1039,  [1043,1040] ) # [34,65])
                ggrid.DivideCube(coordnode, nodecoord, nslices)

                import pdb; pdb.set_trace()

                ggrid.ReturnAllWellFormedCubes(33)

        #import pdb; pdb.set_trace()
        FreshProb = opts.xprob
        NRuns = opts.expansionruns
        sshape = [ggrid.TorLen] * ggrid.NDim
        sshape = sshape + [NRuns]

        
        if ggrid.NDim == 2:
            BigHist = np.zeros(tuple(sshape)).astype("int")

        for ibig in range(NRuns):
            ggrid = cubenodeset_t(opts.dimension)
            #import pdb; pdb.set_trace()
            ggrid.CreateNCube()



            # put a deformation in, to give it some mixed-fractality stuff to grow around
            if ggrid.NDim == 2:
                listpossible = ggrid.ReturnAllWellFormedCubes(33)
                coordnode,nodecoord,coord_path = listpossible[0]

                ggrid.BruteDivideCube(coordnode,nodecoord,coord_path)
            else:
                
                try:
                    node_000 = ggrid.NodeFromCoords(tuple([1] * ggrid.NDim))[0]
                    axisnodes_000 = ggrid.GetAxisNodes(node_000) # this returns the axes in the "positive" direction
                    cubes_000 = ggrid.ReturnAllWellFormedCubes(node_000, None, axisnodes_000)
                    coordnode, nodecoord, coord_path = cubes_000[0]
                    ggrid.BruteDivideCube(coordnode, nodecoord, coord_path)
                except:
                    import pdb; pdb.set_trace()

 

            
            myretval = ggrid.ExpandManyRandomly(Prob, NRuns) 
            
            if ggrid.NDim == 2:
                BigHist[:,:,ibig] = ggrid.Hist

                if ibig == 0:
                    ggrid.Distributions(0)
                    ggrid.Distributions(ggrid.NNode//2)
        

            
                import pickle
                with open('blah5_slice2_expruns100_prob_0p01.pkl', 'wb') as fp:
                    pickle.dump(BigHist, fp)
            print("just finished run", ibig)
        #np.savetxt('blah2.csv', BigHist)



        bDebuggingEarly = True
        if bDebuggingEarly:    



            nsliceslist = [ ggrid.NSlices ]

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
                myretval = ggrid.ExpandManyRandomly(Prob, NRuns, nsliceslist)
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

                #for islice in  [(nslices,)]: #  [(2,),(4,),(2,4,6)]:      
                if True:  

                    ggrid = cubenodeset_t(opts.dimension)
                    ggrid.CreateNCube()

                    NodeCoordOrig = copy(NodeCoord)
                    CoordNodeOrig = copy(CoordNode)


                    #import pdb; pdb.set_trace()
                    ggrid.ExpandManyRandomly(Prob, NRuns, islice)

                    if False:
                        #import pdb; pdb.set_trace()
                        startnode01 = 41 # startnode10 = 32
                        import pdb; pdb.set_trace()
                        allaxes01 = ggrid.ReturnAllAxes(startnode01)   
                        coords,nodes = ggrid.FindCubeCoord(41, [7033,8829])    
                        coords, noodes, iax = ggrid.PickADivisibleCube(startnode01, allaxes01)
                        if not(coords is None):
                            ggrid.DivideCube(coords, noodes, nslices)
                        import pdb; pdb.set_trace()

                        startnode01 = 32 # startnode10 = 32
                        import pdb; pdb.set_trace()
                        allaxes01 = ggrid.ReturnAllAxes(startnode01)       
                        coords, noodes, iax = ggrid.PickADivisibleCube(startnode01, allaxes01)
                        if not(coords is None):
                            ggrid.DivideCube(coords, noodes, nslices)
                        import pdb; pdb.set_trace()
                        #import pdb; pdb.set_trace()

                        startnode01 = 63 # startnode10 = 32
                        import pdb; pdb.set_trace()
                        allaxes01 = ggrid.ReturnAllAxes(startnode01)       
                        coords, noodes, iax = ggrid.PickADivisibleCube(startnode01, allaxes01)
                        if not(coords is None):
                            ggrid.DivideCube(coords, noodes, nslices)
                        import pdb; pdb.set_trace()

                        startnode01 = 32 # startnode10 = 32
                        import pdb; pdb.set_trace()
                        allaxes01 = ggrid.ReturnAllAxes(startnode01)       
                        coords, noodes, iax = ggrid.PickADivisibleCube(startnode01, allaxes01)
                        if not(coords is None):
                            ggrid.DivideCube(coords, noodes, nslices)
                        import pdb; pdb.set_trace()

                        startnode01 = 7 # startnode10 = 32
                        #import pdb; pdb.set_trace()
                        allaxes01 = ggrid.ReturnAllAxes(startnode01)
                        coords, noodes, iax = ggrid.PickADivisibleCube(startnode01, allaxes01)
                        if not(coords is None):
                            ggrid.DivideCube(coords, noodes, nslices)
                        import pdb; pdb.set_trace()


                        startnode01 = 64 # startnode10 = 32
                        #import pdb; pdb.set_trace()
                        allaxes01 = ggrid.ReturnAllAxes(startnode01)
                        coords, noodes, iax = ggrid.PickADivisibleCube(startnode01, allaxes01)
                        if not(coords is None):
                            ggrid.DivideCube(coords, noodes, nslices)
                        #import pdb; pdb.set_trace()




                    ggrid.ScatterPlot()
                    import pdb; pdb.set_trace()

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
    

                    # now analyze implies powers (and p-vals) over a bunch of nodes
                    
                    nsamples = 5
                    nsamples = np.min((ggrid.NNode, nsamples//2))

                    
                    indsamples = rn.sample([i for i in range(ggrid.NNode)], nsamples)
                    rn.shuffle(indsamples)
                    




                    for i in indsamples:

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

                                    nrng = len(rng)

                                    heatres = linregress(rng[front:nrng-i],np.log(arr2[front:nrng-i]))
                                    heatresparab = linregress(rng2[front:nrng-i],np.log(arr2[front:nrng-i]))
                                    heatrestert = linregress(rng3[front:nrng-i],np.log(arr2[front:nrng-i]))
                                    

                                    ndata = len(rng[front:nrng-i])

                                    print("slopeheat", heatres.slope, len(rng) - front)
                                    print("RESHEAT: N: %d node %d nslice %s slope: %f pval %f  slopequadr: %f pvalquadr %f slopetert: %f pvaltert %f " % (ndata, i, str(islice), heatres.slope, heatres.pvalue, heatresparab.slope, heatresparab.pvalue, heatrestert.slope, heatrestert.pvalue))
                                    slopedict[len(rng) - front] = {"HEATSLOPE":heatres.slope, "HEATPVAL":heatres.pvalue,  "TOUCHSLOPE":heatres.slope, "TOUCHPVAL":heatres.pvalue} #, "SLOPEDICT":slopedict}
                            
                                    heatresdict[(tuple(islice), len(rng)-i)] = {"HEATSLOPE":heatres.slope, "HEATPVAL":heatres.pvalue, }
                            

                        for inod in ggrid.NodeVec:
                            inod.Amplitude = 0
                            inod.PrevAmplitude = 0

                        ggrid.NodeVec[i].Amplitude = 1.0
                        ggrid.NodeVec[i].PrevAmplitude = 1.0

                        
                        
                        if True: 
                            print("Touch", ggrid.NNode)
                            touchgrow = ggrid.Touch(120)
                            
                            arr2 = copy(touchgrow)
                            rng = np.log(1.0+np.arange(len(arr2)))
                            rng2 = rng*rng
                            rng3 = rng*rng*rng
                            front = 0


                            slopedict = {}
                            if len(rng) - front > 5:                    
                                for i in range(len(rng)):
                                    nrng = len(rng)
                                    touchres = linregress(rng[front:nrng-i],np.log(arr2[front:nrng-i]))
                                    touchresparab = linregress(rng2[front:nrng-i],np.log(arr2[front:nrng-i]))
                                    touchrestert = linregress(rng3[front:nrng-i],np.log(arr2[front:nrng-i]))

                                    ndata = len(rng[front:nrng-i])

                                    print("touchslope", touchres.slope, len(rng) - front)

                                    print("RESTOUCH: N: %d node %d nslice %s slope: %f pval %f  slopequadr: %f pvalquadr %f slopetert: %f pvaltert %f" % (ndata, i, str(islice), touchres.slope, touchres.pvalue, touchresparab.slope, touchresparab.pvalue, touchrestert.slope, touchrestert.pvalue))

                                    touchresdict[(tuple(islice), len(rng)-i)] = {"TOUCHSLOPE":touchres.slope, "TOUCHPVAL":touchres.pvalue}
                                
                                    #import pdb; pdb.set_trace()

                        

                        
                        for inod in ggrid.NodeVec:
                            inod.Amplitude = 0
                            inod.PrevAmplitude = 0

                        ggrid.NNeighborHistogram()

                heatdf = pd.DataFrame.from_dict(heatresdict,orient="index")
                touchdf = pd.DataFrame.from_dict(touchresdict,orient="index")
                #import pdb; pdb.set_trace()




"""


import pickle
f = open('blah.pkl', 'rb')
ggrid = pickle.load(f)
f.close()

f = open('blah3.pkl', 'rb')
coords,nodes,coordpath = pickle.load(f)
f.close()








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

NDIV 2540 11381
Dividing 29947 {(0, 0): 11381, (1, 0): 6038, (0, 1): 11380, (1, 1): 11383} 29948 0.0 0 2541 (4.333333333333333, 10.592592592592592)
NDIV 2541 3633
Dividing 29959 {(0, 0): 3633, (1, 0): 3629, (0, 1): 3634, (1, 1): 3630} 29960 0.0 0 2542 (15.481481481481481, 5.506172839506172)
NDIV 2542 24688
Dividing 29971 {(0, 0): 24688, (1, 0): 24684, (0, 1): 24687, (1, 1): 24683} 29972 0.0 0 2543 (4.493827160493828, 3.8518518518518516)
NDIV 2543 11416
Dividing 29983 {(0, 0): 11416, (1, 0): 11412, (0, 1): 11415, (1, 1): 11411} 29984 0.0 0 2544 (24.51851851851852, 10.506172839506172)
NDIV 2544 15620
Dividing 29993 {(0, 0): 15620, (1, 0): 15177, (0, 1): 15623, (1, 1): 15181} 29994 0.0 0 2545 (19.62962962962963, 13.555555555555555)
NDIV 2545 29146


%run  /Users/hrvojehrgovcic/quant/latticegas_cubenodes6.py  -t 1000  -p 1.0 --dim 2 --length 200000 --expand 20000

%run  /Users/hrvojehrgovcic/quant/latticegas_cubenodes12.py  -t 1000  --xprob 1.0 --maxnode 5000

%run  /Users/hrvojehrgovcic/quant/latticegas_cubenodes12.py  -t 10000  --xprob 1.0 --maxnode 8000 --dim 3
%run  /Users/hrvojehrgovcic/quant/latticegas_cubenodes12.py  -t 100000  --xprob 1.0 --maxnode 80000 --dim 3

%run  /Users/hrvojehrgovcic/quant/latticegas_cubenodes13.py  -t 100000  --xprob 1.0 --maxnode 80000 --dim 3 --post

%run  /Users/hrvojehrgovcic/quant/latticegas_cubenodes14.py  -t 100000  --xprob 1.0 --maxnode 80000 --dim 3
"""