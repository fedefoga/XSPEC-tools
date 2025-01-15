from xspec import *
import matplotlib.pyplot as plt
import numpy as np
import os, sys, subprocess
import seaborn as sns
sns.set_context('paper')

# Prevent HEADAS prompts
os.environ["HEADASNOQUERY"] = "True"

texpo = int(sys.argv[1])			# Exposure time in ksecs
pindx = sys.argv[2]      	# Cyclotron width
sigma = sys.argv[3]
sname = "IRAS $18293-0914$" 

# XSPEC settings
Xset.chatter = 10
Xset.parallel.leven = 8 
Fit.query = 'yes'

# Model files
model1 = f"pshock_p{pindx}"  # Alternative hypothesis

ftest, fprob, pnorm = np.loadtxt(f'simftest_p{pindx}_t{texpo}_s{sigma}.asc', unpack=True)


## --------------------------------------------------------------
# Load hypothesis model

Xset.restore(model1+'.xcm')

m1 = AllModels(1)
m1.setPars({10:f"{np.mean(pnorm)}"})

## --------------------------------------------------------------------------------

# Fake spectra settings
fsA = FakeitSettings(background="nubkg.pha",
                     response="nustar.rmf",
                     arf="nustar.arf",
                     exposure=float(texpo)*1000.0,
                     backExposure=float(texpo)*1000.0,
                     correction=1.0,
                     fileName='simp_A.fak')

fsB = FakeitSettings(fsA)
fsB.fileName = 'simp_B.fak'

# Generate fake spectra
AllData.fakeit(nSpectra=2, settings={1: fsA, 2: fsB}, filePrefix='simp_')

# Group spectra
for cam in ['A', 'B']:
    cmd = (
        f"ftgrouppha clobber=YES "
        f"infile=simp_{cam}.fak "
        f"outfile=simp_{cam}_grp.fak "
        f"backfile=simp_{cam}_bkg.fak "
        f"respfile=nustar.rmf "
        f"grouptype=optmin "
        f"groupscale=30 "
    )
    subprocess.call(cmd, shell=True)

# Load grouped spectra
AllData(f"1:1 simp_A_grp.fak 2:2 simp_B_grp.fak")
AllData.ignore("bad")
AllData.ignore("*:**-3.0 79.-**")

s1 = AllData(1)
s2 = AllData(2)
rate1 = s1.rate
rate2 = s2.rate
print(f'Expected count rates FPMA: {rate1[0]:.3f} +/- {rate1[1]:.3f}')
print(f'Expected count rates FPMB: {rate2[0]:.3f} +/- {rate2[1]:.3f}')

AllModels.setEnergies('0.1 100 100 log')

m1 = AllModels(1)
m2 = AllModels(2)

m2(1).untie()
m2(1).frozen = False
Fit.perform()

# -----------------------------------------------------------------

Plot.xAxis = "keV"
Plot.add = True
Plot.background = True

Plot("eedata,model")
ncomps = Plot.nAddComps(1)

xd1 = np.asarray(Plot.x(1))
yd1 = np.asarray(Plot.y(1))
xe1 = np.asarray(Plot.xErr(1))
ye1 = np.asarray(Plot.yErr(1))
ym1 = np.asarray(Plot.model(1))
yb1 = np.asarray(Plot.backgroundVals(1))

xd2 = np.asarray(Plot.x(2))
yd2 = np.asarray(Plot.y(2))
xe2 = np.asarray(Plot.xErr(2))
ye2 = np.asarray(Plot.yErr(2))
ym2 = np.asarray(Plot.model(2))
yb2 = np.asarray(Plot.backgroundVals(2))

yc1 = []
[yc1.append(Plot.addComp(addCompNum=q+1, plotGroup=1)) for q in range(ncomps)]

yc2 = []
[yc2.append(Plot.addComp(addCompNum=q+1, plotGroup=2)) for q in range(ncomps)]

# -------------------------------------------------------------------

nrows = 3

fig,(ax) = plt.subplots(nrows,1, dpi=200, gridspec_kw={'hspace':0, 'height_ratios':[1,0.5,0.5]}, 
                        sharex=True, figsize=plt.figaspect(1))
fig.align_labels()
                        
for k in range(nrows):
	ax[k].tick_params(direction='in', which='both', top=True, right=True, labelsize='medium')
	ax[k].tick_params(which='major', length=3.5, width=0.65)
	ax[k].tick_params(which='minor', length=2, width=0.65)
	for pos in ['top','left','right','bottom']:
		ax[k].spines[pos].set_linewidth(0.75)

ax[0].set_title(sname + f"\t $\\Delta T = {texpo:d}$ ks" + f"\t $\\Gamma = {pindx}$", 
fontsize=10)

ax[0].set_ylabel('$E^2$ $F(E)$', fontsize='medium')
ax[1].set_ylabel('$\\Delta \\chi$', fontsize='medium')
ax[-1].set_ylabel('$\\Delta \\chi$', fontsize='medium')
ax[-1].set_xlabel('Energy (keV)', fontsize='medium')

ax[0].set_yscale('log')
ax[-1].set_xscale('log')

ax[1].axhline(0,lw=1,c='k',zorder=10)
ax[-1].axhline(0,lw=1,c='k',zorder=10)


ax[0].set_ylim(0.01*max(max(ym1),max(ym2)),1.75*max(max(ym1),max(ym2)))

ax[-1].set_xlim(2.8,80)
xticks = [3,6,10,20,40,70]
ax[-1].set_xticks(xticks)
ax[-1].set_xticklabels([f"${x}$" for x in xticks])


xticks = [-3,0,3]
ax[1].set_ylim(-4,4)
ax[1].set_yticks(xticks)
ax[1].set_yticklabels([f"${x}$" for x in xticks])

xticks = [0,3]
ax[-1].set_ylim(-3,5)
ax[-1].set_yticks(xticks)
ax[-1].set_yticklabels([f"${x}$" for x in xticks])


##----------------------------------------------------------------------

ax[0].errorbar(xd1, yd1, xerr=xe1, yerr=ye1, ls='', lw=0.75, c='r', zorder=1, label='FPMA')
ax[0].errorbar(xd1, ym1, lw=0.75, c='k', zorder=10)
[ax[0].errorbar(xd1, yc1[k], lw=1, c='r', ls=':', zorder=5) for k in range(ncomps)]

ax[0].errorbar(xd2, yd2, xerr=xe2, yerr=ye2, ls='', lw=0.75, c='b', zorder=1, label='FPMB')
ax[0].errorbar(xd2, ym2, lw=0.75, c='k', zorder=10)
[ax[0].errorbar(xd2, yc2[k], lw=1, c='b', ls=':', zorder=5) for k in range(ncomps)]

ax[0].legend(loc='upper right', ncol=1, fontsize='small', frameon=False, handletextpad=0.25)

ax[0].scatter(xd1, yb1, s=3, marker='d', c='r', zorder=1, alpha=0.75)
ax[0].scatter(xd2, yb2, s=3, marker='d', c='b', zorder=1, alpha=0.75)



##----------------------------------------------------------------------

#m1(6).values = 1e-10
Plot("ldata,model")
xd1 = np.asarray(Plot.x(1))
yd1 = np.asarray(Plot.y(1))
xe1 = np.asarray(Plot.xErr(1))
ye1 = np.asarray(Plot.yErr(1))
ym1 = np.asarray(Plot.model(1))
xd2 = np.asarray(Plot.x(2))
yd2 = np.asarray(Plot.y(2))
xe2 = np.asarray(Plot.xErr(2))
ye2 = np.asarray(Plot.yErr(2))
ym2 = np.asarray(Plot.model(2))

ax[1].errorbar(xd1, (yd1-ym1)/ye1, yerr=1.0, xerr=xe1, lw=0.75, ls='', c='r', zorder=0)
ax[1].errorbar(xd2, (yd2-ym2)/ye2, yerr=1.0, xerr=xe2, lw=0.75, ls='', c='b', zorder=0)



m1(10).values = 1e-10
#m1(10).frozen = True
#Fit.perform()

Plot("ldata,model")
xd1 = np.asarray(Plot.x(1))
yd1 = np.asarray(Plot.y(1))
xe1 = np.asarray(Plot.xErr(1))
ye1 = np.asarray(Plot.yErr(1))
ym1 = np.asarray(Plot.model(1))
xd2 = np.asarray(Plot.x(2))
yd2 = np.asarray(Plot.y(2))
xe2 = np.asarray(Plot.xErr(2))
ye2 = np.asarray(Plot.yErr(2))
ym2 = np.asarray(Plot.model(2))

ax[2].errorbar(xd1, (yd1-ym1)/ye1, yerr=1.0, xerr=xe1, lw=0.75, ls='', c='r', zorder=0)
ax[2].errorbar(xd2, (yd2-ym2)/ye2, yerr=1.0, xerr=xe2, lw=0.75, ls='', c='b', zorder=0)


plt.savefig(f"simftest_p{pindx}_t{texpo}_s{sigma}.png", bbox_inches='tight', dpi=130)
plt.savefig(f"simftest_p{pindx}_t{texpo}_s{sigma}.pdf", bbox_inches='tight')


