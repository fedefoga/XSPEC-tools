from xspec import *
import numpy as np
import os, sys, subprocess


# Prevent HEADAS prompts
os.environ["HEADASNOQUERY"] = "True"

# Input arguments
try:
	texpo = sys.argv[1]  		# Exposure time in ksecs
	pindx = sys.argv[2]  		# Cyclotron energy
	sigma = sys.argv[3]
except IndexError:
    print("Usage: fakeit.py <texpo> <pindx> <sigma>")
    sys.exit(1)

# XSPEC settings
Xset.chatter = 0
Xset.parallel.leven = 8 
Fit.query = 'yes'

# Model files
model0 = "pshock"  			# Null hypothesis
model1 = f"pshock_p{pindx}"  	# Alternative hypothesis

## ----------------------------------------------------------------------

# Load null hypothesis model
Xset.restore(model1 + '.xcm')
m1 = AllModels(1)
norm = m1(10).values[0] 
m1(10).values[0] = norm * (float(sigma)+1)

# Fake spectra settings
fsA = FakeitSettings(background="nubkg.pha", response="nustar.rmf", arf="nustar.arf",
                     exposure=float(texpo)*1000.0, backExposure=float(texpo)*1000.0,
                     correction=1.0, fileName='sim_A.fak')

fsB = FakeitSettings(fsA)
fsB.fileName = 'sim_B.fak'

# Generate fake spectra
AllData.fakeit(nSpectra=2, settings={1: fsA, 2: fsB}, filePrefix='sim_')

# Group spectra
for cam in ['A', 'B']:
    cmd = (
        f"ftgrouppha clobber=YES "
        f"infile=sim_{cam}.fak "
        f"outfile=sim_{cam}_grp.fak "
        f"backfile=sim_{cam}_bkg.fak "
        f"respfile=nustar.rmf "
        f"grouptype=optmin "
        f"groupscale=30 "
    )
    subprocess.call(cmd, shell=True)

# Load grouped spectra
AllData(f"1:1 sim_A_grp.fak 2:2 sim_B_grp.fak")
AllData.ignore("bad")
AllData.ignore("*:**-3.0 70.-**")
#AllModels.setEnergies('0.01 100 100 log')

m1 = AllModels(1)
m2 = AllModels(2)
m2(1).untie()
m2(1).frozen = False
m1(10).frozen = False  ## Thaw powernorm

Fit.perform()

chisq1, dof1 = Fit.statistic, Fit.dof

npars = m1.nParameters

line = f'{chisq1:e} {dof1:e} '
for nn in range(npars):
	if m1(nn+1).frozen : continue
	if m1(nn+1).link != "" : continue
	line += f"{m1(nn+1).values[0]:e} "

line += f"{m2(1).values[0]:e} "


# Save results
output_file = f'simpars_p{pindx}_t{texpo}_s{sigma}.asc'

with open(output_file, 'a+') as f:
    f.write(line+'\n')




