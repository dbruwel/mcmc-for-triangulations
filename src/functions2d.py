### Generate a seed triangulation

import math, random
# import all functions etc. from regina
from regina import *

import sys
sys.path.append("./src/")

# Functions to perform MCMC in d=2
from functions2d import *

def surface(g,n):
    "This function produces a triangulation of the closed orientable surface of genus `g` with `n` triangles."
    if n%2 != 0:
        print("number of triangles must be even")
        return None
    if n < 4*g-2:
        print("a triangulation of a genus",g,"closed orientable surface must have at least",4*g-2,"triangles")
        return None
    
    # create empty surface triangulation
    s = Triangulation2()
    
    # Build fundamental domain disk with 4g-2 triangles and boundary a 4g-gon
    #
    # add first triangle
    s.newSimplex()
    for i in range(4*g-3):
        # add next triangle
        s.newSimplex()
        # take last two triangles and...
        a = s.simplex(s.size()-1)
        b = s.simplex(s.size()-2)
        # glue them together (gluing face iterates, face 2 gets never glued)
        if i%2 == 0:
            a.join(0,b,Perm3())
        else:
            a.join(1,b,Perm3())
    
    # special case if we have a surface
    if g==0:
        s.newSimplex()
        a = s.simplex(0)
        b = s.simplex(1)
        a.join(2,b,Perm3())
        a.join(1,b,Perm3())
        a.join(0,b,Perm3())
    # special case if we have a torus
    elif g==1:
        a = s.simplex(0)
        b = s.simplex(1)
        a.join(2,b,Perm3(2,0,1))
        a.join(1,b,Perm3(1,2,0))
    # glue antipodal edges to each other to obtain genus g surface
    else:
        for k in range(2*g-2):
            a = s.simplex(k)
            b = s.simplex(k+2*g)
            a.join(2,b,Perm3(1,0,2))
        a = s.simplex(2*g-2)
        b = s.simplex(4*g-3)
        a.join(2,b,Perm3(2,0,1))        
        a = s.simplex(2*g-1)
        b = s.simplex(0)
        a.join(2,b,Perm3(2,0,1))   
        
    # add triangles as needed
    while s.size() < n:
        s.pachner(s.simplex(0),True,True)
    return s
#######################################################


### Propose next move 

#######################################################
def alpha(n,gamma):
    "This function computes probability to go up depending on state size `n` and parameter `gamma`. A return value of i means an i-move is proposed next. "
    y= random.random()
    a = math.exp((-1)*gamma*n)
    # maybe to adjust later?
    b = 1
    if y < a:
        return 0
    elif y >= a and y < (b+a)/(b+1):
        return 1
    elif y >=(b+a)/(b+1):
        return 2
#######################################################

### Compute neighbouring triangulations

#######################################################
def neighbours(iso,f,g):
    "This function produces a dictionary of all `g`-neighbours (with isomorphism signatures as keys) of triangulation `iso` with f-vector `f`. This function uses non-standard isomorphism signatures and hence requires `regina` version 7.3 or newer."
    nbrs = {}
    # going up (0-2-moves at every triangle)
    if g == 0:
        for n in range(f[2]):
            # create copy of tri in standard iso sig labelling
            target = Triangulation2.fromIsoSig(iso)
            # test if move is possible and if so, perform it
            if target.pachner(target.triangle(n), True, False):
                target.pachner(target.triangle(n), False, True)
                # get isomorphism signature of result, add it to neighbours
                tiso = target.isoSig_RidgeDegrees()
                # add edge needed to flip to obtain this neighbour (in standard iso sig labelling)
                if not tiso in nbrs:
                    nbrs[tiso]=n
        return nbrs
    # horizontally (1-1-move or edge flip at every edge contained in two triangles)
    elif g == 1:
        for e in range(f[1]):
            # create copy of tri in standard iso sig labelling
            target = Triangulation2.fromIsoSig(iso)
            # test if move is possible and if so, perform it
            if target.pachner(target.edge(e), True, False):
                target.pachner(target.edge(e), False, True)
                # get isomorphism signature of result, add it to neighbours
                tiso = target.isoSig_RidgeDegrees()
                # add edge needed to flip to obtain this neighbour (in standard iso sig labelling)
                if not tiso in nbrs:
                    nbrs[tiso]=e
        return nbrs
    # going doen (2-0-move at every vertex of degree three in three distinct triangles)
    if g == 2:
        for v in range(f[0]):
            # create copy of tri in standard iso sig labelling
            target = Triangulation2.fromIsoSig(iso)
            # test if move is possible and if so, perform it
            if target.pachner(target.vertex(v), True, False):
                target.pachner(target.vertex(v), False, True)
                # get isomorphism signature of result, add it to neighbours
                tiso = target.isoSig_RidgeDegrees()
                # add edge needed to flip to obtain this neighbour (in standard iso sig labelling)
                if not tiso in nbrs:
                    nbrs[tiso]=v
        return nbrs
#######################################################

### Choose move from proposal


###############################################################################
def choosemove(iso, f, gamma):
    "This function takes a state triangulation given by isomorphism signature `iso` with f-vector `f` and paramter `gamma`. It computes a proposal, enumerates neighbours of `iso` and decides wether to perform the proposed move."
    # a==0 go up, a==1 horizontal, a==2 go down
    a = alpha(f[2],gamma)
    ngbrs = neighbours(iso,f,a)
    # get number of neighbours of this type
    num_ngbrs = len(ngbrs.keys())
    # random number for proposal to move or stay
    i = random.random()
    # setup done
    
    # go up
    if a == 0:
        # stay where you are
        if i > float(num_ngbrs)/float(f[2]):
            return iso, f
        # go up (very likely)
        else:
            return random.choice(list(ngbrs.keys())), [f[0]+1,f[1]+3,f[2]+2]
    # edge flip
    if a == 1:
        # stay where you are
        if i > float(num_ngbrs)/float(f[1]):
            return iso, f
        # move sideways (very likely)
        else:
            return random.choice(list(ngbrs.keys())), f        
    # go down
    elif a == 2:
        # stay where you are
        # avoid division by zero
        if f[2] == 2:
            return iso, f
        if i > float(num_ngbrs)/float(f[2]-2):
            return iso, f
        # go down (unlikely)
        else:
            return random.choice(list(ngbrs.keys())), [f[0]-1,f[1]-3,f[2]-2]
###############################################################################

### Main function



def randomise(iso, f, steps, gamma, interval, offset, name):
    "This is the main function taking in see triangulation `iso` with f-vector `f`. It performs a random walk in the Pachner graph of length `steps` with parameter `gamma`. Parameter `verbose` decides print behaviour, `name` is the filename for the output file."
    # initialise number of steps
    st = 0
    histogram = {}
    # save stuff to a file
    with open(name,"w") as fl:
        fl.write("")
    # open output file
    while st < steps + offset*interval-1:
        st += 1
        iso, f = choosemove(iso,f,gamma)
        if (interval != 0 and st % interval == 0 and st >= offset*interval):
            # open output file
            with open(name,"a") as fl:
                fl.write(iso+"\n")
            print("collecting triangulation",int((st-offset*interval)/interval+1),":",iso)
    return True
###############################################################################


def iterate(iso,gamma,steps=1):
    "Iterate the MCMC for 'steps' steps and return the iso of the final triangulation"
    st = 0
    s=Triangulation2.fromIsoSig(iso)
    f = s.fVector()
    while st < steps:
        st += 1
        iso,f = choosemove(iso,f,gamma)
    return iso

def mcmc2d(iso,gamma,samples=10,offset=0,interval=100,verbose=True,printToFile=False):
    "Collect 'samples' samples of triangulations starting at 'iso' with parameter 'gamma'. 'offset' is the number of triangulations to be burnt (discarded initially). 'interval' is the number of triangulations between successive samples"
    samp =0
    
    s=Triangulation2.fromIsoSig(iso)
    f=s.fVector()

    if printToFile!=False:
        name="outputs/mcmc2d_"+str(printToFile)+"_gamma_"+str(gamma)+"_samples_"+str(samples)+".txt"
    # burn
    if offset>0:
        iso = iterate(iso,gamma,steps=int(offset))
        s=Triangulation2.fromIsoSig(iso)
        f=s.fVector()
    
    while samp < samples:
        iso = iterate(iso,gamma,steps=int(interval))
        s=Triangulation2.fromIsoSig(iso)
        f=s.fVector()
        samp+=1
        if verbose:
            print("collecting triangulation",int(samp)," of ",int(samples)," :",iso)
        if printToFile:
            with open(name,"a") as fl:
                fl.write(iso+"\n")
            
