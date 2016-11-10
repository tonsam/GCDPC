# GCDPC
An efficient and scalable variant of Density Peaks Clustering algorithm for large scale datasets

data.zip			            datasets to be used, please unzip it before running the scripts

fun_dim_halo.m			      improved algorithm that clusters data of any dimensions (sparse integral grids，with halo step line 141-176)

fun_imp_clean.m			      improved algorithm special for 2D data (2D integral image, with or without halo removing need to change part of code)

fun_imp_dim.m			        improved algorithm that clusters data of any dimensions (sparse integral grids，with halo step line 141-176)

fun_ori.m			            original DPC implement

integral.m			          high dimension integrid implement for dc evaluating

mex_dc_evaluate.mexw64		mex function to evaluate dc via 2D integral image

mex_permutohedral.mexw64	permutohedral lattice that used to speed up global Gaussian filtering

normalize.m			          function to normalized the data before each algorithm

DPCProcess.m              the common process function that the hierarchy alogrithm use

fun_imp_hierarchy.m       the main function of hierarchy algorithm, based on the input parameters to select which level to run.

parameter_CBD_trips_locations_0901.mat  the example input data of the second level clustering

mex_permutohedral.mexw64  the new version of permutohedral which normalize the output in the end
