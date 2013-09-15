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

Installing 
---------- 
Typing `make` at the command line should do the trick. It
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

	[x,y,status] = coneos_direct(A,b,c,cones,params)

Usage in C 
---------- 
If `make` completes successfully, it will produce two
static library files, `libconeosdir.a` and `libconeosindir.a` under the `lib`
folder. To include the libraries in your own source code, compile with the
linker option with `-L(PATH_TO_coneos)\lib` and `-lconeosdir` or `-lconeosindir` (as
needed).

These libraries (and `coneos.h`) expose only three API functions:

* Sol * coneos(Data \* d, Cone \* k)
    
	This solves the problem specified in the `Data` and `Cone` structures.  The
	solution is returned in a `Sol` structure.
    
* void freeData(Data \* d, Cone \* k)
    
	This frees the `Data` and `Cone` structures.
    
* void freeSol(Sol \* sol)

	This frees the `Sol` structure.
    
The three relevant data structures are:

    typedef struct PROBLEM_DATA {
      idxint n, m; /* problem dimensions */
      /* problem data, A, b, c: */
      double * Ax; 
      idxint * Ai, * Ap; 
      double * b, * c;
  
      Params * p;
    } Data;
        
    typedef struct PROBLEM_PARAMS {
      idxint MAX_ITERS, CG_MAX_ITS;
      double EPS_ABS, ALPHA, CG_TOL;
      idxint VERBOSE, NORMALIZE;  // boolean
    } Params;

    typedef struct SOL_VARS {
      idxint n, m; /* solution dimensions */
      double *x, *s, *y; 
      char status[16];
    } Sol;

    typedef struct CONE {
        idxint f;          /* number of linear equality constraints */
        idxint l;          /* length of LP cone */
        idxint *q;             /* array of second-order cone constraints */
        idxint qsize;      /* length of SOC array */
    } Cone;


The data matrix `A` is specified in column-compressed format and the vectors
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
