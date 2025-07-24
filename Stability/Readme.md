Codes for producing linear stability curves and implicit bifurcation curves for the Swift--Hohenberg equation with a spatial heterogeneity.

Here are your next steps:

To produce figures like Figure 4(a), type:

Stability_Plot(R, k);

where R denotes the width of the spatial heterogeneity and k is the angular fourier index of the linear eigenfunction.

To produce figures like Figure 5, type: 

mult_match(kmax);

This code plots the first implicit curve of bifurcation points in (epsilon,R)-space, where kmax is the max angular Fourier index to be plotted (plots include k=0,1,...,kmax). 
