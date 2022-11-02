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
atomMoveTriplets = np.zeros(iMax+1)
fullBoxTriplets = np.zeros(iMax+1)
fullBoxSetUp = np.zeros(iMax+1)
fullBoxExp = np.zeros(iMax+1)
fullBoxRoot = np.zeros(iMax+1)
atomMoveSetUp = np.zeros(iMax+1)
atomMoveExp = np.zeros(iMax+1)
atomMoveShare = np.zeros(iMax+1)

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
    fullBoxMat[counter,:] = data[0,:]
    fullBoxSetUp[counter] = fullBoxMat[counter,1]+fullBoxMat[counter,7]
    fullBoxExp[counter] = fullBoxMat[counter,2]
    fullBoxRoot[counter] = fullBoxMat[counter,6]
    fullBoxTriplets[counter] = fullBoxMat[counter,3]+fullBoxMat[counter,4]+fullBoxMat[counter,5]
    atomMoveData = data[1:nData,:]
    atomMoveData = np.sum(atomMoveData,axis=0)
    atomMoveMat[counter,:] = atomMoveData
    atomMoveSetUp[counter] = atomMoveMat[counter,1]+atomMoveMat[counter,6]+atomMoveMat[counter,7]
    atomMoveExp[counter] = atomMoveMat[counter,2]
    atomMoveShare[counter] = atomMoveMat[counter,7]
    atomMoveTriplets[counter] = atomMoveMat[counter,3]+atomMoveMat[counter,4]+atomMoveMat[counter,5]
    counter = counter + 1
    if (maxNodes > 1):
      data = pa.read_csv(str(num)+string,  delim_whitespace=True, header=None)
      data = data.values
      nData = len(data)
      fullBoxMat[counter,:] = data[0,:]
      fullBoxSetUp[counter] = fullBoxMat[counter,1]+fullBoxMat[counter,7]
      fullBoxExp[counter] = fullBoxMat[counter,2]
      fullBoxRoot[counter] = fullBoxMat[counter,6]
      fullBoxTriplets[counter] = fullBoxMat[counter,3]+fullBoxMat[counter,4]+fullBoxMat[counter,5]
      atomMoveData = data[1:nData,:]
      atomMoveData = np.sum(atomMoveData,axis=0)
      atomMoveMat[counter,:] = atomMoveData
      atomMoveSetUp[counter] = atomMoveMat[counter,1]+atomMoveMat[counter,6]+atomMoveMat[counter,7]
      atomMoveExp[counter] = atomMoveMat[counter,2]
      atomMoveShare[counter] = atomMoveMat[counter,7]
      atomMoveTriplets[counter] = atomMoveMat[counter,3]+atomMoveMat[counter,4]+atomMoveMat[counter,5]
      counter = counter + 1
  else:
    data = pa.read_csv(str(num)+string,  delim_whitespace=True, header=None)
    data = data.values
    nData = len(data)
    fullBoxMat[counter,:] = data[0,:]
    fullBoxSetUp[counter] = fullBoxMat[counter,1]+fullBoxMat[counter,7]
    fullBoxExp[counter] = fullBoxMat[counter,2]
    fullBoxRoot[counter] = fullBoxMat[counter,6]
    fullBoxTriplets[counter] = fullBoxMat[counter,3]+fullBoxMat[counter,4]+fullBoxMat[counter,5]
    atomMoveData = data[1:nData,:]
    atomMoveData = np.sum(atomMoveData,axis=0)
    atomMoveMat[counter,:] = atomMoveData
    atomMoveSetUp[counter] = atomMoveMat[counter,1]+atomMoveMat[counter,6]+atomMoveMat[counter,7]
    atomMoveExp[counter] = atomMoveMat[counter,2]
    atomMoveShare[counter] = atomMoveMat[counter,7]
    atomMoveTriplets[counter] = atomMoveMat[counter,3]+atomMoveMat[counter,4]+atomMoveMat[counter,5]
    counter = counter + 1


fullT0 = fullBoxMat[0,0]
moveT0 = atomMoveMat[0,0]
tripT0 = atomMoveTriplets[0]
fTripT0 = fullBoxTriplets[0]

'''
for j in range (0,iMax+1):
  fullBoxMat[j,0]=fullT0/fullBoxMat[j,0]
  atomMoveMat[j,0]=moveT0/atomMoveMat[j,0]
  atomMoveTriplets[j] = tripT0/atomMoveTriplets[j]
  fullBoxTriplets[j] = fTripT0/fullBoxTriplets[j]
'''

fullT0 = int(fullT0)
moveT0 = int(moveT0)
tripT0 = int(tripT0)
fTripT0 = int(fTripT0)

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
ax1.scatter(nProcs,fullBoxMat[:,0],s=75,c="g",marker="o",label=r"t$_{\mathrm{total}}$")
ax1.plot(nProcs,fullBoxSetUp,c="r",label=r"t$_{\mathrm{set}}$")
ax1.plot(nProcs,fullBoxSetUp+fullBoxExp,c="b",label=r"t$_{\mathrm{set}}$+t$_{\mathrm{exp}}$")
#ax1.plot(nProcs,fullBoxTriplets+fullBoxSetUp+fullBoxExp,c="m",label=r"t$_{\mathrm{trip}}$+t$_{\mathrm{set}}$+t$_{\mathrm{exp}}$")
ax1.fill_between(x=nProcs,y1=fullBoxSetUp,y2=0,color="r",alpha=0.2)
ax1.fill_between(x=nProcs,y1=fullBoxSetUp+fullBoxExp,y2=fullBoxSetUp,color="b",alpha=0.2)
for i in range(0,iMax+1):
  ax1.plot([nProcs[i],nProcs[i]],[fullBoxMat[i,0]-fullBoxTriplets[i],fullBoxMat[i,0]],c="m")
mp.grid()
mp.xlim(0,121)
#mp.ylim(0,1.5)
mp.xlabel(r"$N_{p}}$")
mp.ylabel(r"t$_{N_{p}}$ / s")
mp.legend(loc="upper right")
mp.savefig(dirString+'/'+'fullBoxContribPlot-Np_120_Zoom_NormDens.pdf',bbox_inches = "tight")

figg=mp.figure()
ax2=figg.add_subplot(111)
ax2.scatter(nProcs,atomMoveMat[:,0],s=75,c="g",marker="o",label=r"t$_{\mathrm{total}}$")
ax2.plot(nProcs,atomMoveSetUp,c="r",label=r"t$_{\mathrm{set}}$")
ax2.plot(nProcs,atomMoveSetUp+atomMoveExp,c="b",label=r"t$_{\mathrm{set}}$+t$_{\mathrm{exp}}$")
#ax2.plot(nProcs,atomMoveTriplets+atomMoveSetUp+atomMoveExp,c="m",label=r"t$_{\mathrm{trip}}$+t$_{\mathrm{set}}$+t$_{\mathrm{exp}}$")
ax2.fill_between(x=nProcs,y1=atomMoveSetUp,y2=0,color="r",alpha=0.2)
ax2.fill_between(x=nProcs,y1=atomMoveSetUp+atomMoveExp,y2=atomMoveSetUp,color="b",alpha=0.2)
for i in range(0,iMax+1):
  ax2.plot([nProcs[i],nProcs[i]],[atomMoveMat[i,0]-atomMoveTriplets[i],atomMoveMat[i,0]],c="m")
mp.grid()
mp.xlim(0,121)
#mp.ylim(0,4.1)
mp.xlabel(r"$N_{p}}$")
mp.ylabel(r"t$_{N_{p}}$ / s")
mp.legend(loc="upper right")
mp.savefig(dirString+'/'+'atomMoveContribPlot-Np_120_Zoom_NormDens.pdf',bbox_inches = "tight")
