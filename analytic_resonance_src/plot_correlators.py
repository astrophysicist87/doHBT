#!/usr/bin/env python
from numpy import *
from pylab import *
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

eps = 1.e-10
eps2barvec = ['0', '0.05', '0.1']
v2barvec = ['0', '0.05', '0.1']
KTvec = ['0', '1', '2']
colors = ['red', 'blue', 'green']
#qaxes = ['$q_x$', '$q_y$']
#qaxes2 = ['qx', 'qy']
qaxes = ['$q_o$', '$q_s$']
qaxes2 = ['qo', 'qs']

##################################################################
##################################################################
##################################################################

def generate_plot(ie2b, iv2b, sliceaxis, slicevalue):
	#set-up
	plotfontsize = 18
	fig = plt.figure()
	ax = fig.add_axes([0.125,0.125,0.825,0.825])
	
	#set axes limits
	xlower, xupper = -0.75, 0.75
	ylower, yupper = 0.85, 2.05
	
	ax.axhline(1.0, color='black', linewidth=2)
	
	for iKT in xrange(len(KTvec)):
		#load data
		filename = 'results/azavg_correlfunct1D_os_Pion_p_KT_%(sKT)s_eps2bar_%(se2b)s_v2bar_%(sv2b)s.dat' % {"sKT": KTvec[iKT], "se2b": ie2b, "sv2b": iv2b}
		data = loadtxt(filename)
		
		#plot data
		plotdata = data[(where(abs(data[:,sliceaxis]-slicevalue) <= eps))]
		nonsliceaxis = 1 - sliceaxis
		ax.plot(plotdata[:,nonsliceaxis],plotdata[:,2], color=colors[iKT], linestyle='solid', linewidth=2, label=r'$K_T = %(sKT)s$ GeV' % {"sKT": KTvec[iKT]})
		ax.plot(plotdata[:,nonsliceaxis],plotdata[:,4], color=colors[iKT], linestyle='--', linewidth=2)

	#update axes limits and axes labels
	ax.axis([xlower, xupper, ylower, yupper])
	ax.set_xlabel(qaxes[nonsliceaxis] + ' (GeV)', {'fontsize': plotfontsize + 5})
	ax.set_ylabel(r'$\left< C \right>_{\Phi_K}$, $\left< \tilde{C} \right>_{\Phi_K}$', {'fontsize': plotfontsize + 5})
	
	ax.legend(loc='best',ncol=3,prop={'size': plotfontsize - 4})
	
	plt.title(r'$\bar{\epsilon}_2 = %(se2b)s,\,\bar{v}_2 = %(sv2b)s$' % {"se2b": ie2b, "sv2b": iv2b}, fontsize = plotfontsize)
	outfile = 'results/azavg_correlfunct1D_Pion_p_%(qa)s=0_eps2bar_%(se2b)s_v2bar_%(sv2b)s.pdf' % {"se2b": ie2b, "sv2b": iv2b, "qa": qaxes2[sliceaxis]}
	#plt.savefig(outfile, format='pdf')
	#print 'Saved to', outfile
	plt.show()




##################################################################
##################################################################
##################################################################


if __name__ == "__main__":
	slicevalue = 0.0
	for ie2b in eps2barvec:
		for iv2b in v2barvec:
			for isa in xrange(2):
				#isa = ith sliceaxis
				generate_plot(ie2b, iv2b, isa, slicevalue)

# End of file
