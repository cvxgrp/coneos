coneos
============================================================ 
A C package for solving large-scale convex cone problems.

based on "Operator Splitting for Conic Optimization" by  
Brendan O’Donoghue, Eric Chu, Neal Parikh, and Stephen Boyd

This code provides a solver for convex cone problems. It is an
implementation of the algorithm described in [this
paper](http://www.stanford.edu/~boyd/). It provides both a direct and an
indirect solver in the form of a static library for inclusion in other
projects.

It simultaneously solves the primal cone program

	minimize     c'*x subject to   A*x + s == b s in K 
                 
and its dual

	maximize     -b'*y subject to   -A'*y == c y in K^* 

where `K` is a product cone of free cones, linear cones `{ x | x >= 0 }`, 
second-order cones `{ (t,x) | ||x||_2 <= t }`, and semi-definite cones `{ X | X psd }`
`K^*` is its dual cone.

Sparse and dense versions
---------- 
The sparse package uses the LDL and AMD packages written by Davis et. al.  
The dense package uses BLAS and LAPACK libraries. To use the dense version
you must point conOS.mk to the locaation of these installed libraries.

Solving SDPs
---------- 
In order to solve SDPs you must have BLAS and LAPACK installed (even for the sparse
version of the package). Point coneOS.mk to the location of these libraries. Without
these you can still solve SOCPs and LPs, however.

Installing 
---------- 
Typing `make` at the command line in either coneOSsparse or coneOSdense. It
will produce two libaries, `libconeosdir.a` and `libconeosindir.a` found under the
`lib` folder. As a byproduct, it will also produce two demo binaries under the
`bin` folder called `demo_direct` and `demo_indirect`.

One caveat: if you have a 32-bit version of Matlab and use the build process
below (for Matlab), then if you try to make the libraries (on a 64-bit
machine), you must `make purge` before `make` again.

File an issue with us if the build process fails for you.

### Compiling a Matlab mex file 
Running `make_coneos` in Matlab under the
`matlab` folder will produce two usable mex files

If `make_coneos` fails and complains about an incompatible architecture, edit the
`make_coneos` file according to the comments.

Remember to include the `matlab` directory in your Matlab path if you wish to
use the mex file in your Matlab code. The calling sequence is

	[x,y,status] = coneos_direct(data,cones,params)

where data contains A, b, c  
cones contains f (free/zero cone size), l (linear cone size), q (array of SOCs), s (array of SDCs)

Usage in C 
---------- 
If `make` completes successfully, it will produce two
static library files, `libconeosdir.a` and `libconeosindir.a` under the `lib`
folder. To include the libraries in your own source code, compile with the
linker option with `-L(PATH_TO_coneos)\lib` and `-lconeosdir` or `-lconeosindir` (as
needed).

These libraries (and `coneos.h`) expose only three API functions:

* Sol * coneos(Data \* d, Cone \* k, Info * info)
    
	This solves the problem specified in the `Data` and `Cone` structures,
    and returns the solution in the Sol struct and various information about the run in
    the Info struct.

* void freeData(Data \* d, Cone \* k)
    
	This frees the `Data` and `Cone` structures.
    
* void freeSol(Sol \* sol)

	This frees the `Sol` structure.
    
The four relevant data structures are:

    /* struct that containing standard problem data */
    typedef struct PROBLEM_DATA {
      int n, m; /* problem dimensions */
      /* problem data, A, b, c: */
      double * Ax; 
      int * Ai, * Ap; 
      int Anz;
      double * b, * c;
      int MAX_ITERS, CG_MAX_ITS;
      double EPS_ABS, ALPH, CG_TOL, UNDET_TOL, RHO_X;
      int VERBOSE, NORMALIZE;  // boolean
    } Data;

    typedef struct Cone_t {
        int f;          /* number of linear equality constraints */
        int l;          /* length of LP cone */
        int *q;         /* array of second-order cone constraints */
        int qsize;      /* length of SOC array */
        int *s;         /* array of SD constraints */
        int ssize;      /* length of SD array */
    } Cone;

    /* contains primal-dual solution vectors */
    typedef struct SOL_VARS {
        double * x, * y, *s; 
    } Sol;

    /* contains terminating information */
    typedef struct INFO {
        int iter;
        char status[16];
        double pobj;
        double dobj;
        double presid;
        double dresid;
        double gap;
        double time;
    } Info;

The data matrix `A` is specified in column-compressed format for coneOSsparse or dense
column major order for coneOSdense, and the vectors
`b` and `c` are specified as dense arrays. The solutions `x` (primal) and `y`
(dual) are returned as dense arrays. Cones are specified in terms of their
lengths; the only special one is the second-order cone, where the lengths are
provided as an array of second-order cone lengths (and a variable `qsize`
giving its length).


Scalability
----------- 
Note that this code is merely meant as an
implementation of the ideas in our paper. The actual code does not use more
than a single CPU. Nevertheless, for problems that fit in memory on a single
computer, this code will (attempt to) solve them.

To scale this solver, one must either provide a distributed solver for linear
systems or a distributed matrix-vector multiplication.
