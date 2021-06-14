# plasma_physics_reactor
Code which generates equilibrium profiles of a plasma, given geometric and main factors.

For the first run, you need to install Interpos code on MATLAB. With Interpos directory on your workspace MATLAB path, just type "mex interposf90" on the terminal.

After that, change all saving/loading paths in the code to a directory in your computer. The functions which you need to do this are 'iterate_pedestal_reactor', 'create_shape_reactor', 'get_fbt_params_reactor', 'build_greenfunc_reactor', 'create_em_model_reactor', 'gs_solver', 'find_xpoints'.

Before running the code, run 'build_greenfunc_reactor'. This function creates a table with values of the Green function of the system (sum of green functions of each coil). After that, run 'create_em_mode_reactor'. That will create a 'static' file. This function make take a while to run. Both Green Function table and static files are already in this repository, but I strongly suggest you to run this functions. **Every time you change a coil position and/or add/remove a coil, you need to run both functions in the order above**.


Every run will be saved in a directory called '_shotXXXXXX_'. The directory name is defined in the main function (iterate_pedestal_reactor) on the variable _shot_. Make sure that you do not run a shot over an important simulation. 

All the input parameters are given in among lines 4 to 23. These values ajust all the equilibria for a shot. All the variables ending with 0 are parameters on the center of plasma (e.g. Te0). All the variables ending with 'ped' are parameters on the top of pedestal structure. The pedestal top is in the position 'psi_n = 1 - d_e' for electrons or 'psi_n = 1 - d_i' for ions, where 'd_e, d_i' are the pedestal width (in psi_n space) for electrons and ions. The variables ending with 'sep' are parameters on the separatrix (psi_n = 1).

If there are any problems or bugs within the code, please send an e-mail to felipe.machado.salvador@usp.br describing what and where is the problem.
