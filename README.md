# MCMC methods for sampling triangulations of manifolds

The is a fork of the original repository of the same name from https://github.com/jspreer.

The intention of this fork is to modify the sampling algorithm to include a "potential" function that encourages the sampler to explore areas of the space that we are interested in (for example tringulations with high degree Alexander polynomial, etc.)

**The original README is below.**

---

This repository contains supplementary material for the article 

<i>Sampling triangulations of manifolds using Monte Carlo methods</i> by <a href="https://www.maths.usyd.edu.au/u/ega/">Eduardo Altmann</a> and <a href="https://sites.google.com/view/jonathan-spreer/">Jonathan Spreer</a>, <a href="https://arxiv.org/abs/2310.07372"> arXiv:2310.07372 </a>.

Start at the notebooks <a href="https://github.com/edugalt/MCMCForTriangulations/blob/main/tutorial2d.ipynb" target=_blank> tutorial2d.ipynb </a> and <a href="https://github.com/edugalt/MCMCForTriangulations/blob/main/tutorial3d.ipynb" target=_blank> tutorial3d.ipynb </a> for minimal working examples in dimensions 2 and 3 (Algorithm 1 and 2), exploring the space of all triangulations of (generalised triangulations of) surfaces and/or (1-vertex triangulations of) 3-manifolds.

The source code of the MCMC method is in the folder src/

# Documentation

## Triangulations basics

A $d$-dimensional manifold $M$ is a topological space that locally looks like Euclidean $d$-space. A triangulation $T$ of $M$ is a subdivision of $M$ into $d$-simplices glued along their $(d-1)$-dimensional face such that the underlying space $|T|$ of $T$ (the $d$-simplices factored by their gluings) is homeomorphic to $M$.

The *face-vector* or *f-vector* $f(T) = (f_0, f_1, f_2, .... , f_d)$ of $T$ stores the number of $i$-dimensional simplices $f_i$ in $T$ in all dimensions $0 \leq i \leq d$. 

Example: The boundary of the tetrahedron $\Delta$ has face vector $f(\Delta) = (4,6,4)$.

## Regina basics

### 2d 
When $d=2$, $T$ is always a triangulation of a surface (or $2$-dimensional manifold) $S$, i.e., a set of triangles, identified along their edges. Data type is `regina.Triangulation2()`. Given a triagnulation `T`, type in `print T.detail()` to obtain an overview of what is stored in `T`, and print `T.fVector()` for its $f$-vector.

We can go from one triangulation `T` of $S$ to any other triangulation `T'` of $S$ by applying *Pachner moves* -- local modifications to `T` that change the isomorphism type of `T`, but not the topology of the underlying surface. Pachner moves in dimension two in regina are given by commands 

- `T.pachner(t)`: Stellar subdivision of triangle `t` into three triangles by placing a new vertex into the center of `t`. This move adds one vertex, three edges, and two triangles to `T` (it replaces `t` by three new triangles).
- `T.pachner(e)`: The *edge flip*. In the quadrilateral spanned by the two triangles of `T` containing edge `e`, remove `e` and replace it by the other diagonal `e'` in the quadrilateral. This move leaves the $f$-vector invariant.
- `T.pachner(v)`: Inverse of the stellar subdivision `T.pachner(t)`. This move removes vertex `v` to merge three triangles into one.



### 3d
When $d=3$, $T$ is always a triangulation of a $3$-dimensional manifold $M$, i.e., a set of tetrahedra, identified along their triangles. Data type is `regina.Triangulation3()`. Given a triagnulation `T`, type in `print T.detail()` to obtain an overview of what is stored in `T`, and print `T.fVector()` for its $f$-vector.

We can go from one triangulation `T` of $M$ to any other triangulation `T'` of $M$ by applying *Pachner moves* -- local modifications to `T` that change the isomorphism type of `T`, but not the topology of the underlying surface. Pachner moves in dimension three in regina are given by commands 

- `T.pachner(tet)`: Stellar subdivision of tetrahedron `tet` into four tetrahedra by placing a new vertex into the center of `tet`. This move adds one vertex, four edges, six triangles, and three tetrahedra to `T` (it replaces `tet` by four new tetrahedra). (This move is not used in Algorithm 2 below.)
- `T.pachner(t)`: Replaces two tetrahedra $\Delta_1$ and $\Delta_2$ joined along triangle `t` by three tetrahedra around a new edge with endpoints the two vertices of $\Delta_1$ and $\Delta_2$ opposite `t`. This move adds one edge, two triangles, and one tetrahedron.
- `T.pachner(e)`: This is the inverse of the previous move, taking three three tetrahedra joined around edge `e` and replacing it with two tetrahedra joined along the base triangle `t`. This move removes one edge, two triangles, and one tetrahedron.
- `T.pachner(v)`: Inverse of the stellar subdivision `T.pachner(tet)`. This move removes vertex `v` to merge four tetrahedra into one. (This move is not used in Algorithm 2 below.)

### Usage `c.pachner(face,check,perform)`

R3ginas functions to perform moves are all organised the same way: a move always returns `True` or `False` on whether or not the operation was allowed / successful. It changes the underlying triangulation (if applicable), and takes as argument the face where the move is performed (the face that is going to be removed from the triangulation), and two booleans as optional arguments. The first one determining whether to check that the move is possible, the second one to determine whether to perform the move.

I.e., the command `T.pachner(v,True,False)`  takes as input a triangulation `T` and checks whether vertex `v` is contained in exactly three distinct triangles, returns `True` or `False` accordingly, but does not alter `T`.
