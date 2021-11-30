import pandas as pa
import matplotlib as mpl
import matplotlib.pyplot as mp
import numpy as np
import sys

# Get max no. of processes and nodes used
maxProcs = int(sys.argv[1])
maxNodes = int(sys.argv[2])
dirString = sys.argv[3]
iMax = int(maxProcs/maxNodes)
if (maxNodes > 1):
  nProcs = np.arange(0,iMax+1)
  nProcs = nProcs*maxNodes
  nProcs[0] = 1
else:
  nProcs = np.arange(0,iMax+1)
  nProcs = nProcs + 1
fullBoxMat = np.zeros((iMax+1,8)) # 4
atomMoveMat = np.zeros((iMax+1,8)) # 4

# Loop over number of processes, extracting all data from appropriate output file
counter = 0
for i in range (1,iMax+1):
  #counter = i-1
  num = i*maxNodes
  print (num)
  string = '-data.txt'
  if (i == 1):
    data = pa.read_csv(str(i)+string,  delim_whitespace=True, header=None)
    data = data.values
    nData = len(data)
    fullBoxMat[counter,:] = data[4,:]
    atomMoveData = data[1:nData,:]
    atomMoveData = np.sum(atomMoveData,axis=0)
    atomMoveMat[counter,:] = atomMoveData
    counter = counter + 1
    if (maxNodes > 1):
      data = pa.read_csv(str(num)+string,  delim_whitespace=True, header=None)
      data = data.values
      nData = len(data)
      fullBoxMat[counter,:] = data[4,:]
      atomMoveData = data[1:nData,:]
      atomMoveData = np.sum(atomMoveData,axis=0)
      atomMoveMat[counter,:] = atomMoveData
      counter = counter + 1
  else:
    data = pa.read_csv(str(num)+string,  delim_whitespace=True, header=None)
    data = data.values
    nData = len(data)
    fullBoxMat[counter,:] = data[4,:]
    atomMoveData = data[1:nData,:]
    atomMoveData = np.sum(atomMoveData,axis=0)
    atomMoveMat[counter,:] = atomMoveData
    counter = counter + 1

fullT0 = fullBoxMat[0,4]
moveT0 = atomMoveMat[0,4]
for j in range (0,iMax+1):
  fullBoxMat[j,4]=fullT0/fullBoxMat[j,4]
  atomMoveMat[j,4]=moveT0/atomMoveMat[j,4]

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
ax1.scatter(nProcs,fullBoxMat[:,4],s=75,c="g",marker="o",label="Full Box")
mp.grid()
mp.xlabel(r"$N_{p}}$")
mp.ylabel(r"t$_1$ / t$_{N_{p}}$")
#mp.legend(loc="upper left")
mp.savefig(dirString+'/'+'fullBoxOptSumSpeedUp.pdf',bbox_inches = "tight")

figg=mp.figure()
ax2=figg.add_subplot(111)
ax2.scatter(nProcs,atomMoveMat[:,4],s=75,c="g",marker="o",label="Atom Move")
mp.grid()
mp.xlabel(r"$N_{p}}$")
mp.ylabel(r"t$_1$ / t$_{N_{p}}$")
#mp.legend(loc="upper left")
mp.savefig(dirString+'/'+'atomMoveOptSumSpeedUp.pdf',bbox_inches = "tight")
