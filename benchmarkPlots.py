import pandas as pa
import matplotlib as mpl
import matplotlib.pyplot as mp
import numpy as np
import sys

# Get max no. of processes used
maxProcs = int(sys.argv[1])
dirString = sys.argv[2]
nProcs = np.arange(maxProcs)
nProcs = nProcs+1
fullBoxMat = np.zeros((maxProcs,6)) # 4
atomMoveMat = np.zeros((maxProcs,6)) # 4

# Loop over number of processes, extracting all data from appropriate output file
for i in range (1,maxProcs+1):
  counter = i-1
  string = '-data.txt'
  data = pa.read_csv(str(i)+string,  delim_whitespace=True, header=None)
  data = data.values
  nData = len(data)
  fullBoxMat[counter,:] = data[0,:]
  atomMoveData = data[1:nData,:]
  atomMoveData = np.sum(atomMoveData,axis=0)
  atomMoveMat[counter,:] = atomMoveData

fullT0 = fullBoxMat[0,0]
moveT0 = atomMoveMat[0,0]
for j in range (0,maxProcs):
  fullBoxMat[j,0]=fullT0/fullBoxMat[j,0]
  atomMoveMat[j,0]=moveT0/atomMoveMat[j,0]

# Set font sizes for figure
SMALL_SIZE = 17
MEDIUM_SIZE = 19
BIGGER_SIZE = 21
mp.rc('font', size=SMALL_SIZE)          # controls default text sizes
mp.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
mp.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
mp.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
mp.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
mp.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize

# Plot the data
fig=mp.figure()
ax1=fig.add_subplot(111)
ax1.scatter(nProcs,fullBoxMat[:,0],s=75,c="g",marker="o",label="Full Box")
mp.grid()
mp.xlabel(r"$N_{p}}$")
mp.ylabel(r"t$_1$ / t$_{N_{p}}$")
#mp.legend(loc="upper left")
mp.savefig(dirString+'/'+'fullBoxSpeedUp.pdf',bbox_inches = "tight")

figg=mp.figure()
ax2=figg.add_subplot(111)
ax2.scatter(nProcs,atomMoveMat[:,0],s=75,c="g",marker="o",label="Atom Move")
mp.grid()
mp.xlabel(r"$N_{p}}$")
mp.ylabel(r"t$_1$ / t$_{N_{p}}$")
#mp.legend(loc="upper left")
mp.savefig(dirString+'/'+'atomMoveSpeedUp.pdf',bbox_inches = "tight")
