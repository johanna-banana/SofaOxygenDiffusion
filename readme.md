# SofaOxygenDiffusion
This is a project to wrap Sofa code in a C++ class for use directly in Matlab. The simulation executes a diffusion simulation for emulating oxygen diffusion and uptake through tissue.

# Use in Matlab

The required files (not including binary libraries from Sofa):
 - SofaOxygenDiffusion.mexw64
 - mex_interface.m
 - str2fun.m

The object can be created in matlab using either

`sofasim = mex_interface(str2fun('SofaOxygenDiffusion'))`

or

`sofasim = mex_interface(@SofaOxygenDiffusion)` <- This option does not require the "str2fun.m" file as mentioned above.

## Matlab testing
These lines of code should run when in the folder where the .mexw64 file is (after compiling).

```matlab
diffusion_fem = mex_interface(@SofaOxygenDiffusion);
diffusion_fem.init_simulation();
diffusion_fem.run_simulation();
```
