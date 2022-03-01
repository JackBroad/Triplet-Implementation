import pandas as pa
import matplotlib as mpl
import matplotlib.pyplot as mp
import numpy as np
import sys

# Set up time vectors for full sim box
trueFullBoxVec = ([438.,72.2,39.1,29.7,27.1,22.,18.5,17.1,17.9])
dummyFullBoxVec = ([337.,63.1,42.6,29.5,19.8,20.3,15.4,13.8,12.6])

# Do the same for the atom move code and nProcs
trueMoveVec = ([2.17,.527,.234,.189,.156,.122,.0965,.0871,.0746])
dummyMoveVec = ([2.02,.376,.248,.168,.109,.0989,.0890,.0729,.0631])
nProcVec = ([1,5,10,15,20,25,30,35,40])

# Convert to speed-up data
trueFull0 = trueFullBoxVec[0]
dummyFull0 = dummyFullBoxVec[0]
trueMove0 = trueMoveVec[0]
dummyMove0 = dummyMoveVec[0]
for j in range (0,len(nProcVec)):
  trueFullBoxVec[j]=trueFull0/trueFullBoxVec[j]
  dummyFullBoxVec[j]=dummyFull0/dummyFullBoxVec[j]
  trueMoveVec[j] = trueMove0/trueMoveVec[j]
  dummyMoveVec[j] = dummyMove0/dummyMoveVec[j]

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
ax1.scatter(nProcVec,trueFullBoxVec,s=75,c="g",marker="o",label="True Sum")
ax1.scatter(nProcVec,dummyFullBoxVec,s=75,c="r",marker="o",label="Dummy Sum")
mp.grid()
mp.xlabel(r"$N_{p}}$")
mp.ylabel(r"t$_1$ / t$_{N_{p}}$")
mp.legend(loc="upper left")
mp.savefig('scatterByTriplets/realVsDummySum/fullBoxRealVsDummySum.pdf',bbox_inches = "tight")

figg=mp.figure()
ax2=figg.add_subplot(111)
ax2.scatter(nProcVec,trueMoveVec,s=75,c="g",marker="o",label="True Sum")
ax2.scatter(nProcVec,dummyMoveVec,s=75,c="r",marker="o",label="Dummy Sum")
mp.grid()
mp.xlabel(r"$N_{p}}$")
mp.ylabel(r"t$_1$ / t$_{N_{p}}$")
mp.legend(loc="upper left")
mp.savefig('scatterByTriplets/realVsDummySum/atomMoveRealVsDummySum.pdf',bbox_inches = "tight")
