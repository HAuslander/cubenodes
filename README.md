# cubenodes
Create expanding D-dimensional graphs on which dynamics (e.g. random walk, or Brownian-Huygens propagation) leads naturally to a Euclidean embedding



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
