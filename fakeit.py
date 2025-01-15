from xspec import *
import matplotlib.pyplot as plt
import numpy as np
import os, sys, subprocess
from astropy.io import fits

from ftest import *

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
AllData.ignore("*:**-3.0 79.-**")
AllModels.setEnergies('0.01 100 1000 log')

m1 = AllModels(1)
m2 = AllModels(2)
m2(1).untie()
m2(1).frozen = False

m1(10).frozen = False  ## Thaw powernorm

Fit.perform()
chisq1, dof1 = Fit.statistic, Fit.dof

Plot('data,model')
r1 = (np.asarray(Plot.y(1))-np.asarray(Plot.model(1)))/np.asarray(Plot.yErr(1))
r2 = (np.asarray(Plot.y(2))-np.asarray(Plot.model(2)))/np.asarray(Plot.yErr(2))
rs1 = np.append(r1,r2)

norm = m1(10).values[0]

## ---------------------------------------------------------------------------------

Xset.restore(model0+'.xcm')  # Null hypothesis fitting
m1 = AllModels(1)
m2 = AllModels(2)
m2(1).untie()
m2(1).frozen = False

Fit.perform()
chisq0, dof0 = Fit.statistic, Fit.dof

Plot('data,model')
r1 = (np.asarray(Plot.y(1))-np.asarray(Plot.model(1)))/np.asarray(Plot.yErr(1))
r2 = (np.asarray(Plot.y(2))-np.asarray(Plot.model(2)))/np.asarray(Plot.yErr(2))
rs0 = np.append(r1,r2)

## -----------------------------------------------------------------------------------

ftest, fprob = ftest(rs0, rs1)


# Save results
output_file = f'simftest_p{pindx}_t{texpo}_s{sigma}.asc'
with open(output_file, 'a+') as f:
    f.write(f'{ftest:e} {fprob:e} {norm:e}\n')



## for kk in {1..1000}; do for tt in {60..240..30}; do for pp in 1.5 2 2.5; do python3 fakeit.py $tt $pp ; done ; done ; done



