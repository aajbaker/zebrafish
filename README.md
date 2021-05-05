# zebrafish
Written in collaboration with the Fraser Lab at USC. Calculates lagged correlations between neural activity throughout the zebrafish brain while both awake and asleep in an attempt to identify distinct patterns of neural behavior in awake and sleeping zebrafish.

smaller_data_tif: 4D neuroimaging data taken from sleep/wake zebrafish

plot_lag.py: takes in 4D data (3D images x time) and outputs figures that describe the analysis.

Figure1: The correlations of each region with respect to the seed region overlaid onto a sample 2D plane

Figure2: The lagged time series of a sample region with respect to the time series of the seed

Figure3: The correlations of the sample region with respect to the seed at each lag step

This was part of a larger body of work conducted over the summer session.
