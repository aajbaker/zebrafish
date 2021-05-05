"""
Calculates lagged correlations between neural activity throughout the zebrafish brain 
while both awake and asleep in an attempt to identify distinct patterns of neural behavior in awake and sleeping zebrafish.
"""

__author__ = "Aaron Baker"
__contact__ = "aaron@aajbaker.com"
__credits__ = [""]
__date__ = "2018/09"
__version__ = "0.0.1"

import numpy as np
import matplotlib.pyplot as plt
import tifffile as tiff
import matplotlib.patches as patch
import matplotlib.colors as colors
import matplotlib.cm as cmx
import scipy.stats as sci
import math
import matplotlib

### PARAMETERS ###

data = tiff.imread('smaller_data.tif')
seed_x, seed_y, seed_z =  14, 12, 4
z_pos = 4
min_lag = -20
max_lag = 20
thresh = 0.7
p_thresh = .00001
sig_ts = []
ts_example = []

def corr_with_lag(array1, array2, min_lag, max_lag):
	""" Returns the lag number that gives the maximum correlation between two arrays. """
	assert max_lag < len(array1) and len(array1) == len(array2)
	correlations, lag_num = [], None
	#find it for the negative lags (the leads) if the minimum is less than 0
	for lag in range(-min_lag, 0, -1):
		corr = np.corrcoef(array2[:(len(array2) - lag)], array1[lag:])[1, 0]
		correlations.append(corr)
		if corr == max(correlations):
			max_corr = corr
			lag_num = -lag
	#find the lag with max correlation with lag from 0 to max_lag
	for lag in range(0, max_lag + 1):
		corr = np.corrcoef(array1[:(len(array1) - lag)], array2[lag:])[1, 0]
		correlations.append(corr)
		if corr == max(correlations):
			max_corr = corr
			lag_num = lag
	return max_corr, lag_num, correlations

def lag_with_thresh(array1, array2, thresh):
	corr, lag, correlations = corr_with_lag(array1, array2, min_lag, max_lag)
	if thresh == None or corr > thresh:
		return corr, lag, correlations
	else:
		return float('nan'), float('nan'), []

def zeros_array(data, dimentions):
	""" Make a matrix of 0s either 2D or 3D based on the sized of those dimentions in data. """
	assert dimentions == 2 or dimentions == 3
	z, y, x = len(data[0]), len(data[0][0]), len(data[0][0][0])
	if dimentions == 3:	
		return np.zeros((z, y, x))
	else:
		return np.zeros((y, x))

def time_series(x, y, z, data):
	""" Calculates the time series from the 3D coordinate across time t from original data. """
	series = np.array([])
	for t in data:
		series = np.append(series, t[z][y][x])
		f0 = min(series)
		norm = [(f - f0)/f0 for f in series]
	return norm

def fill_lag(matrix, z, data, seed):
	""" Takes in a 2D matrix, its z value and mutates the 2D matrix with 
	the time series across time t from the data for each point. """
	for y in range(len(matrix)):
		for x in range(len(matrix[y])):
			timeser = time_series(x, y, z, data)
			corr, lag, correlations = lag_with_thresh(seed, timeser, thresh)
			matrix[y][x] = lag
			if not(math.isnan(lag)) and timeser != seed:
				sig_ts.append((timeser, lag, corr, (x, y), correlations))
	return matrix

def most_correlated(n, ts):
	highest_ts = []
	for lag in range(min_lag, max_lag+1):
		constrained = [t[2] for t in ts if t[1] == lag]
		if constrained:
			max_corr = max(constrained)
			for t in ts:
				if t[2] == max_corr:
					highest_ts.append(t)
	return highest_ts

### BUILDING THE MATRICES FOR THE LAGS ###

seed = time_series(seed_x, seed_y, seed_z, data)
empty_lags = zeros_array(data, 2)
lags = fill_lag(empty_lags, z_pos, data, seed)
highest_ts = most_correlated(5, sig_ts)
example = highest_ts[int(len(highest_ts)//1.3)]
for i in [min_lag, min_lag//2, 0, max_lag//2, max_lag]:
	if i < 0:
		ts_example.append((example[0][:i], i))
	else:
		ts_example.append((example[0][i:], i))

### FIGURE BUILDING ###
matplotlib.rcParams.update({'font.size': 15})
#coloring the scale in accordance with the coloring of the above subplot; source: https://stackoverflow.com/questions/8931268/using-colormaps-to-set-color-of-line-in-matplotlib
rainbow = cma = plt.get_cmap('rainbow') 
cNorm  = colors.Normalize(vmin=min_lag, vmax=max_lag)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=rainbow)

#Figure 1: Overlay of correlations on pixels of visual 2D data
fig1 = plt.figure(1)
#make the subplot for the heat map of lags on the original image dimensions
a = fig1.add_subplot(111)
cmap = matplotlib.cm.rainbow
norm = colors.Normalize(vmin=min_lag, vmax=max_lag)
plt.imshow(data[0][z_pos], cmap='gray')
imgplot = plt.imshow(lags, cmap=cmap)
imgplot.set_norm(norm)
#adding a circle patch over the plot where the example is located
a.add_patch(patch.RegularPolygon((example[3][0], example[3][1]), 4, .7, fill=False, orientation=math.pi/4, color='yellow'))
#adding a circle patch over the plot where the seed is located
if z_pos == seed_z:	
	a.add_patch(patch.Circle((seed_x, seed_y), radius=.5, color='green'))
a.set_title('Lag Map (R > ' + str(thresh) + ')')
a.axis('off')
#adding the reference for the heat map
plt.colorbar(fraction=0.03, pad=0.04)

#Figure 2: time series of the significant lags
fig2, axs = plt.subplots(len(ts_example)+2, 1, sharex=True)
fig2.subplots_adjust(hspace=0)
#plot all the time series from a single pixel at different lags, colored by lag times
for i in range(len(ts_example)):
	lag_num = ts_example[i][1]
	y_vals = ts_example[i][0]
	if lag_num < 0:
		x_vals = range(abs(lag_num), abs(lag_num)+len(y_vals))
	else:
		x_vals = range(0, len(y_vals))
	colorVal = scalarMap.to_rgba(lag_num)
	axs[i].plot(x_vals, y_vals, color=colorVal)
	axs[i].get_yaxis().set_ticks([])
	axs[i].get_yaxis().set_label_text(str(lag_num))
#use the subplot between the example and seed as a buffer
axs[len(ts_example)].axis('off')
#plot the seed time series in black
axs[len(ts_example)+1].plot(seed, color='black')
axs[len(ts_example)+1].get_yaxis().set_ticks([])
axs[len(ts_example)+1].get_yaxis().set_label_text('Seed')
plt.suptitle('Time Series at Different Lags for X=' + str(example[3][0]) + ', Y=' + str(example[3][1]))
fig2.text(0.04, 0.5, 'Lag Number', va='center', rotation='vertical')
fig2.text(0.4, 0.02, 'Time in Frames', va='center')

#Figure 3: Scatter plot of correlations at specific lag times
fig3 = plt.figure(3)
plt.scatter(range(min_lag, max_lag+1), example[4], c=range(min_lag, max_lag+1), cmap='rainbow', s=190)
plt.ylabel('Correlation Coefficient')
plt.xlabel('Lag Number')
plt.title('Correlation With Seed at Each Lag')

plt.show()


