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
    fullBoxTriplets[counter] = fullBoxMat[counter,3]+fullBoxMat[counter,4]+fullBoxMat[counter,5]
    atomMoveData = data[1:nData,:]
    atomMoveData = np.sum(atomMoveData,axis=0)
    atomMoveMat[counter,:] = atomMoveData
    atomMoveTriplets[counter] = atomMoveMat[counter,3]+atomMoveMat[counter,4]+atomMoveMat[counter,5]
    counter = counter + 1
    if (maxNodes > 1):
      data = pa.read_csv(str(num)+string,  delim_whitespace=True, header=None)
      data = data.values
      nData = len(data)
      fullBoxMat[counter,:] = data[0,:]
      fullBoxTriplets[counter] = fullBoxMat[counter,3]+fullBoxMat[counter,4]+fullBoxMat[counter,5]
      atomMoveData = data[1:nData,:]
      atomMoveData = np.sum(atomMoveData,axis=0)
      atomMoveMat[counter,:] = atomMoveData
      atomMoveTriplets[counter] = atomMoveMat[counter,3]+atomMoveMat[counter,4]+atomMoveMat[counter,5]
      counter = counter + 1
  else:
    data = pa.read_csv(str(num)+string,  delim_whitespace=True, header=None)
    data = data.values
    nData = len(data)
    fullBoxMat[counter,:] = data[0,:]
    fullBoxTriplets[counter] = fullBoxMat[counter,3]+fullBoxMat[counter,4]+fullBoxMat[counter,5]
    atomMoveData = data[1:nData,:]
    atomMoveData = np.sum(atomMoveData,axis=0)
    atomMoveMat[counter,:] = atomMoveData
    atomMoveTriplets[counter] = atomMoveMat[counter,3]+atomMoveMat[counter,4]+atomMoveMat[counter,5]
    counter = counter + 1

fullT0 = fullBoxMat[0,0]
moveT0 = atomMoveMat[0,0]
tripT0 = atomMoveTriplets[0]
fTripT0 = fullBoxTriplets[0]
for j in range (0,iMax+1):
  fullBoxMat[j,0]=fullT0/fullBoxMat[j,0]
  atomMoveMat[j,0]=moveT0/atomMoveMat[j,0]
  atomMoveTriplets[j] = tripT0/atomMoveTriplets[j]
  fullBoxTriplets[j] = fTripT0/fullBoxTriplets[j]

'''
fullT0 = int(fullT0)
moveT0 = int(moveT0)
#tripT0 = 150*tripT0
tripT0 = int(tripT0)
fTripT0 = int(fTripT0)
'''
#print(fullT0)#+' '+moveT0+' '+tripT0+' '+fTripT0)

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
a,b = np.polyfit(nProcs[0:5],fullBoxMat[0:5,0],1)
print(a,b)
aStr=str(round(a,2))
bStr=str(round(b,2))
strFullT0=str(round(fullT0,2))
fig=mp.figure()
ax1=fig.add_subplot(111)
ax1.scatter(nProcs,fullBoxMat[:,0],s=75,c="g",marker="o",label=r"t$_1$="+strFullT0+" s")
ax1.plot(nProcs,a*nProcs+b,c='r',label=aStr+r"$x$"+r"+"+bStr)
mp.ylim(0,25)
mp.grid()
mp.xlabel(r"$N_{p}}$")
mp.ylabel(r"t$_1$ / t$_{N_{p}}$")
mp.legend(loc="upper left")
mp.savefig(dirString+'/'+'fullBoxSpeedUpConds_Np-40_normDens.pdf',bbox_inches = "tight")

a,b = np.polyfit(nProcs[0:5],atomMoveMat[0:5,0],1)
print(a,b)
aStr=str(round(a,2))
bStr=str(round(b,2))
strMoveT0=str(round(moveT0,2))
figg=mp.figure()
ax2=figg.add_subplot(111)
ax2.scatter(nProcs,atomMoveMat[:,0],s=75,c="g",marker="o",label=r"t$_1$="+strMoveT0+" s")
ax2.plot(nProcs,a*nProcs+b,c='r',label=aStr+r"$x$"+r"+"+bStr)
#mp.ylim(0,25)
mp.grid()
mp.xlabel(r"$N_{p}}$")
mp.ylabel(r"t$_1$ / t$_{N_{p}}$")
mp.legend(loc="upper left")
mp.savefig(dirString+'/'+'atomMoveSpeedUpConds_Np-40_normDens.pdf',bbox_inches = "tight")

a,b = np.polyfit(nProcs[0:5],atomMoveTriplets[0:5],1)
print(a,b)
aStr=str(round(a,2))
bStr=str(round(b,2))
strTripT0=str(round(tripT0,2))
figgg=mp.figure()
ax3=figgg.add_subplot(111)
ax3.scatter(nProcs,atomMoveTriplets,s=75,c="g",marker="o",label=r"t$_1$="+strTripT0+" s")
ax3.plot(nProcs,a*nProcs+b,c='r',label=aStr+r"$x$"+r"+"+bStr)
#mp.ylim(0,40)
mp.grid()
mp.xlabel(r"$N_\mathrm{p}}$")
mp.ylabel(r"t$_1$ / t$_{N_{\mathrm{p}}}$")
mp.legend(loc="upper left")
mp.savefig(dirString+'/'+'calc-barrier_moveSpeedUp_normDens.pdf',bbox_inches = "tight")

a,b = np.polyfit(nProcs[0:5],fullBoxTriplets[0:5],1)
print(a,b)
aStr=str(round(a,2))
bStr=str(round(b,2))
strFTripT0=str(round(fTripT0,2))
figggg=mp.figure()
ax4=figggg.add_subplot(111)
ax4.scatter(nProcs,fullBoxTriplets,s=75,c="g",marker="o",label=r"t$_1$="+strFTripT0+" s")
ax4.plot(nProcs,a*nProcs+b,c='r',label=aStr+r"$x$"+r"+"+bStr)
#mp.ylim(0,40)
mp.grid()
mp.xlabel(r"$N_\mathrm{p}}$")
mp.ylabel(r"t$_1$ / t$_{N_{\mathrm{p}}}$")
mp.legend(loc="upper left")
mp.savefig(dirString+'/'+'calc-barrier_fullSpeedUp_normDens.pdf',bbox_inches = "tight")
