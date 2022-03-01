import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as mp


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

####################Exp Time############################

timeExp = ([0.0193911,0.0096833,0.0060620,0.0033405,0.0026031,0.0025322,0.0020811,0.0017583,0.0014351,0.0017056,0.0019554,0.0016424,0.0017559,0.0013582,0.0015980,0.0014347,0.0015443,0.0015263,0.0019368,0.0009713,0.0013051])
tZeroExp = timeExp[0]
nProc = ([1,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64,68,72,76,80])
for i in range(0,len(timeExp)):
  timeExp[i] = tZeroExp/timeExp[i]


fig1=mp.figure()
ax1=fig1.add_subplot(111)
ax1.scatter(nProc,timeExp,s=75,c="g",marker="o",label=r"t$_1$ = "+str(tZeroExp)+" s")
mp.grid()
mp.xlabel(r"$N_{p}}$")
mp.ylabel(r"t$_1$ / t$_{N_{p}}$")
mp.legend(loc="upper left")
mp.savefig('expMoveSpeedUpSharedMem.pdf',bbox_inches = "tight")
