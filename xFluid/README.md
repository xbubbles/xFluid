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


## Pressure solver
xFluid contains two different pressure solvers:   
1. Normal pressure sovler;  
2. Variational pressure solver.  

** Normal pressure solver **
This pressure solver directly calculate pressure by solve the uncompressible equation:  
$$ \nabla u = 0$$
