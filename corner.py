import os
import sys
import numpy as np
import pandas as pd
import seaborn as sns
sns.set_context('paper')
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.special import erf
from scipy.interpolate import interpn
from matplotlib.ticker import FuncFormatter

## -----------------------------------------------------------------------------------------
# Utility function to round a number to 1 significant digit
def round_to_1(x):
    return round(x, -int(np.floor(np.log10(x))))

# Utility function to round a number to match the order of magnitude of a reference value
def round_to_reference(x, y):
    return round(x, -int(np.floor(np.log10(y))))

def scientific_notation(x):
	expo = np.floor(np.log10(x))
	if expo==0: return f'${x:.0f}$'
	x = np.log10(x) - expo
	x = 10**x
	return f'${x:.0f}{{\\cdot}}10^{{{expo:.0f}}}$'


## -----------------------------------------------------------------------------------------
def plot_corner(data, columns, labels, filename, skip, nbins, logpars, 
				density=True, ms=2, cmap='viridis', histcolor='k', fs=6):

	Npar = len(columns) - skip # Number of parameters to plot (excluding seg,chi2,dof columns)

	# Set up a grid of subplots
	fig, ax = plt.subplots(Npar, Npar, sharex='col', dpi=200, figsize=plt.figaspect(1), 
		gridspec_kw={"hspace": 0, 'wspace': 0})

    # Loop to configure each subplot in the corner plot grid
	for k in range(Npar):    		## k == rows
		for q in range(Npar):		## q == columns
			if q>k: 
				ax[k][q].axis("off")
			
			elif q==k:

				ax[q][q].spines['top'].set_linewidth(0)
				ax[q][q].spines['right'].set_linewidth(0)
				ax[q][q].spines['left'].set_linewidth(0.5)			
				ax[q][q].spines['bottom'].set_linewidth(0.5)
				if q==0: ax[q][q].spines['left'].set_linewidth(0)
				
			else:
				for _,s in ax[k][q].spines.items(): s.set_linewidth(0.5)
		
			## Deactivate all ticks
			ax[k][q].tick_params(axis='both', which='both', length=0, width=0, 
			bottom=False, left=False,labelbottom=False, labelleft=False)


	for qq in range(Npar-1):
		## Activate left-ticks in first column (without histogram)
		ax[qq+1][0].tick_params(axis='y', which='major', length=3, width=0.25, 
		labelleft=True, left=True, labelsize=fs-1, pad=0.75)

		## Activate bottom-ticks in last row (without histogram)
		ax[-1][qq].tick_params(axis='x', which='major', length=3, width=0.25, 
		labelbottom=True, bottom=True, labelsize=fs-1, pad=0.75)



    # Plot data for diagonal histograms and off-diagonal scatter plots
	for q in range(Npar):				## COLUMNS
		xd = data[columns[q+skip]]  # Select data for the current parameter

		# Determine whether to use a logarithmic scale for the x-axis
		if columns[q+skip] in logpars:
			ax[q][q].set_xscale('log')
			xb = np.geomspace(min(xd), max(xd), nbins)
		else:
			xb = np.linspace(min(xd), max(xd), nbins)

		# Plot histogram for the diagonal
		ax[q][q].hist(xd, bins=xb, histtype='stepfilled', color=histcolor, 
		alpha=0.25, density=density)
		
		ax[q][q].hist(xd, bins=xb, histtype='step', lw=0.5, color='k', 
		zorder=10, density=density)
		
		# Label the parameter
		ax[q][q].set_title(labels[q+skip], loc='center', pad=4, fontsize=fs)  

		ax[q][q].set_xlim(0.99 * min(xd), 1.01 * max(xd))


		for k in range(q+1, Npar):			## ROW
			yd = data[columns[k+skip]]

			# Check for log parameters and create unique bins
			if columns[k+skip] in logpars:
				ax[k][q].set_yscale('log')
				yb = np.geomspace(min(yd[yd > 0]), max(yd), nbins)
			else:
				yb = np.linspace(min(yd), max(yd), nbins)

			xb = np.unique(xb)
			yb = np.unique(yb)

			# Remove NaN and invalid values from data
			valid_mask = ~(xd.isna() | yd.isna())
			xd = xd[valid_mask]
			yd = yd[valid_mask]

			dd, xe, ye = np.histogram2d(xd, yd, bins=(xb, yb), density=density)

			zd = interpn(
				(0.5 * (xe[1:] + xe[:-1]), 0.5 * (ye[1:] + ye[:-1])),
				dd, np.vstack([xd, yd]).T, method="splinef2d",
				bounds_error=False, fill_value=1e-10 )

			idx = zd.argsort()
			x, y, z = xd.iloc[idx], yd.iloc[idx], zd[idx]

			ax[k][q].scatter(x, y, c=z, cmap=cmap, s=ms, marker='.', 
			edgecolor='none', zorder=10)
			ax[k][q].set_xlim(0.99 * min(xd), 1.01 * max(xd))
			ax[k][q].set_ylim(0.99 * min(yd), 1.01 * max(yd))

			if q==0:
				if columns[k+skip] in logpars:
					yticks = np.geomspace(min(yd), max(yd), 6)
					yticks = [yticks[1],yticks[2],yticks[3],yticks[4]]
					ylabels = [scientific_notation(x) for x in yticks]
				else:
					yticks = np.linspace(min(yd), max(yd), 6)
					yticks = [yticks[1],yticks[2],yticks[3],yticks[4]]
					ylabels = [f'${x:.1f}$' for x in yticks]
			
			
				ax[k][0].set_yticks(yticks, labels=ylabels)
				ax[k][0].set_ylabel(labels[k+skip], labelpad=2, fontsize=fs)
			
		
			if k==Npar-1:
				if columns[q+skip] in logpars:
					xticks = np.geomspace(min(xd), max(xd), 6)
					xticks = [xticks[1],xticks[2],xticks[3],xticks[4]]
					xlabels = [scientific_notation(x) for x in xticks]
				else:
					xticks = np.linspace(min(xd), max(xd), 6)
					xticks = [xticks[1],xticks[2],xticks[3],xticks[4]]
					xlabels = [f'${x:.1f}$' for x in xticks]
			
				ax[-1][q].set_xticks(xticks, labels=xlabels, rotation=90)
				ax[-1][q].set_xlabel(labels[q+skip], labelpad=2, fontsize=fs)  

	fig.align_labels()

	# Save the plot as a PNG file
	plt.savefig(f'{filename}.png', bbox_inches='tight', dpi=300)
	plt.savefig(f'{filename}.pdf')



## -------------------------------------------------------------------------------------
# Define the model name and columns used in the data
texpo = sys.argv[1]
pindx = sys.argv[2]
sigma = 3

# Base columns applicable to all models
pcols = ["chisquare", "dof", "nhabs", "temp", "abund", "tau", "pshnorm", "plnorm", "cons"]

plabs = ['$\\chi^2$', 'dof', 'N$_{\\rm H}$', '$kT$', "$Z$", "$\\tau_{\\rm u}$",
		'$\\mathcal{N}_{\\rm psh}$', '$\\mathcal{N}_{\\rm pl}$', "$C_{\\rm AB}$"]

pdata = pd.read_csv(f"simpars_p{pindx}_t{texpo}_s{sigma}.asc", sep='\\s+', names=pcols, 
					usecols=range(len(pcols)))

# Specify parameters to be plotted on a logarithmic scale
logpars = ["tau", "pshnorm", "plnorm"]

## -----------------------------------------------------------------------------

# Generate the corner plot
cmap = 'plasma'
color = 'indigo'

plot_corner(pdata, pcols, plabs, f'corner_p{pindx}_t{texpo}_s{sigma}', 2, 20, 
logpars, False, 5, cmap, color, 7)



