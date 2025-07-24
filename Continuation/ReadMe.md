*These codes have either been written in 2021-2025 by Dan J. Hill or adapted from codes by David J.B. Lloyd (University of Surrey) and Daniele Avitabile, "Numerical computation of coherent structures in spatially-extended neural networks", Second International Conference on Mathematical Neuroscience, Antibes Juan-les-Pins, 2016.*

**For a more thorough tutorial in secant continuation and access to the original codes, see https://github.com/danieleavitabile/continuation-spatially-extended-systems-tutorial.**

Before accessing any of the available folders, you should initialise this code by running

> Init;

These codes are for the continuation of localised dihedral patterns in the 2-3 Swift-Hohenberg equation with a spatial heterogeneity. Here are your next steps:

If you want to find and continue localised dihedral patterns in the 2-3 SH equation, type:

> branch = Cont_Patch(p, Dir);

where p=[mu, gamma, kappa, k, R, m, n+1] contains the parameters of the system including the bifurcation parameter mu, the respective quadratic and cubic coefficients (gamma, kappa), the bifurcating angular mode k, the width of the heterogeneity R, the dihedral lattice order m, and the dimension (n+1) of the reduced ODE system. Dir indicates the direction and step size of the continuation routine, which must be 'pl' (plus), 'mn' (minus), 'sp' (small plus), or 'sm' (small minus). For example:

> branch = Cont_Patch([0.27, 0, -1, 1, 2.8, 1, 15],'pl');

If you want to plot the bifurcation curve of a localised solution, type:
> ExploreBifurcationDiagram('FolderName/branch.mat',idVar);

where FolderName is the data folder that is created by Cont_Patch, and idVar is the column of the branch measures to be plotted. For example, for a data folder called 'D2_Patch_pl', you would type:

>> ExploreBifurcationDiagram('D1_Bif_2.8_pl/branch.mat',5);
