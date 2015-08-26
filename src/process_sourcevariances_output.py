####################################################################
# Author: Christopher Plumberg
# Date: August 24, 2015
####################################################################
# script intended for plotting/analysis of output files
# (containing interpolation grid points) from sourcevariances codes
####################################################################

from numpy import *
import sys
from os import path
#from glob import glob
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.cm as cm

n_resonances = 320
# contains particle_info array index of chosen resonance for analysis
chosen_resonance_to_plot = 1

####################################################################
def plotSVdata(givenSVdata, pxgrid, pygrid, filename):
	#fig, axs = plt.subplots(1, 1)
	plotfontsize = 12
	fig = plt.figure()
	#ax = fig.add_subplot(111, projection='3d')
	ax = fig.gca(projection='3d')
	localcm = plt.cm.get_cmap('rainbow')
	frac = 0.9
	signs = 0.5 * frac * (sign(givenSVdata) + 1.0) + 0.5 * frac
	ax.scatter(pxgrid, pygrid, log(abs(givenSVdata)), c=signs, edgecolor='', s=5, cmap=localcm)
	ax.view_init(elev=90., azim=0.0)
	# sc = ...; plt.colorbar(sc)
	fig.suptitle(filename, fontsize=plotfontsize)
	#plt.show()

	return plt

####################################################################
def reformatSVdata(givenSVdata, pxgrid, pygrid, filename):
	

	return plt


####################################################################
if __name__ == "__main__":
	# set the names of all input files from the command line
	#thermal_sourcevariance_grid_filenames = map(lambda x: path.basename(x), sys.argv[1:])
	thermal_sourcevariance_grid_filenames = sys.argv[1:]
	
	# get the number of source variances to analyze
	number_of_sourcevariances = len(thermal_sourcevariance_grid_filenames)

	# current working directory
	currentWorkingDirectory = path.dirname(thermal_sourcevariance_grid_filenames[0])

	# load grid points in pT - pphi space
	# (this part will have to be modified for more general pT and pphi files)
	pTpts, pTwts = loadtxt(currentWorkingDirectory + '/pT_pts.dat').T
	pphipts, pphiwts = loadtxt(currentWorkingDirectory + '/pphi_pts.dat').T
	pTgrid, pphigrid = meshgrid(pTpts, pphipts)
	npT = len(pTpts)
	npphi = len(pphipts)
	pxgrid = pTgrid * cos(pphigrid)
	pygrid = pTgrid * sin(pphigrid)
	
	useSVdata = False

	# for each source variance, read in appropriate file and plot results
	if useSVdata:
		for iSV in range(number_of_sourcevariances):
			inFile = thermal_sourcevariance_grid_filenames[iSV]
			SVdata = loadtxt(inFile).reshape([n_resonances, npT, npphi])
			#SVplot = plotSVdata(SVdata[chosen_resonance_to_plot], pxgrid, pygrid, path.basename(inFile))
			#outfile = path.splitext(inFile)[0] + '.pdf'
			#print outfile
			#SVplot.show()
			#SVplot.savefig(outfile, format='pdf')
			outfile = path.splitext(inFile)[0] + '_res' + str(chosen_resonance_to_plot) +'_reformatted.dat'
			print 'Saving to $CWD/' + path.basename(outfile)
			# this format is useful for gnuplot, which is a little easier to use for right now...
			finalResult = vstack((pxgrid.flatten(), pygrid.flatten(), SVdata[chosen_resonance_to_plot].flatten())).T
			savetxt(outfile, finalResult, fmt='%10.8f   '*3)
	else:
		for iSV in range(number_of_sourcevariances):
			inFile = thermal_sourcevariance_grid_filenames[iSV]
			data = loadtxt(inFile)
			outfile = path.splitext(inFile)[0] + '_reformatted.dat'
			print 'Saving to $CWD/' + path.basename(outfile)
			# this format is useful for gnuplot, which is a little easier to use for right now...
			finalResult = vstack((pxgrid.flatten(), pygrid.flatten(), data.flatten())).T
			savetxt(outfile, finalResult, fmt='%10.8f   '*3)

	print 'Finished all.'

# End of file