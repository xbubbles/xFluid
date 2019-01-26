# xFluid Source

A brief introduction of the algorithms that xFluid used.

## Framework
Routine of time step:  
1. Update directed distance field; update material.  
2. Reconstruct fluid surface mesh.  
3. Apply body forces.  
4. Pressure solve; apply pressure.  
5. Extrapolate velocity field.  
6. Advect particles and update velocity field.  

&emsp;
&emsp;

## Pressure solver
xFluid contains two different pressure solvers:   
1. Normal pressure sovler;  
2. Variational pressure solver.  

### Normal pressure solver
This pressure solver directly calculate pressure by solve the uncompressible equation.   
Discrete the equation in 3D grids, and we will get a positive definite symmetric linear system. We solve the system by **PCG** solver.  


### Variational pressure solver
This solver is an implementation of  
*"A fast variational framework for accurate solid-fluid coupling", Batty, SIGGRAPH'07*.  
https://cs.uwaterloo.ca/~c2batty/papers/Batty07.pdf  
In xFluid, we concentrate fluid only.  
In this varational framework, we solve pressure by minimizing the kinetic energy of the fluid system. We first write down the total kinetic energy of fluid, and then take derivate, after which we will get a positive definite symmetric linear system.

&emps;
&emps;

## Advection
The Advection framework xFluid used is FLIP&PIC.

&emps;
&emps;

## Meshing
The particles in FLIP will be used here to calculate a directed distance field. We construct a surface mesh of fluid using Marching Cubes method.  
&emps;
A link to an introduction of Marching Cubes method:  
http://paulbourke.net/geometry/polygonise/  

