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

timeExp = ([19.07*10**(-3),9.974*10**(-3),5.101*10**(-3),3.564*10**(-3),2.545*10**(-3),2.364*10**(-3),1.899*10**(-3),1.694*10**(-3),1.511*10**(-3),1.707*10**(-3),1.835*10**(-3),1.527*10**(-3),1.712*10**(-3),1.153*10**(-3),1.592*10**(-3),1.650*10**(-3)])
tZeroExp = timeExp[0]
nProc = ([1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30])
#timeExp = timeExp[0:8] #Take first eight as plateau comes in thereafter
for i in range(0,len(timeExp)):
  timeExp[i] = tZeroExp/timeExp[i]
nProcExp = nProc[0:8]

m,c = np.polyfit(nProc,timeExp,1) # When interested in best fit params, uncomment timeExp line and change nProc to nProcExp here (do the same for accept data)

fig1=mp.figure()
ax1=fig1.add_subplot(111)
ax1.scatter(nProc,timeExp,s=75,c="g",marker="o",label=r"t$_1$ = "+str(tZeroExp)+" s")
mp.grid()
mp.xlabel(r"$N_{p}}$")
mp.ylabel(r"t$_1$ / t$_{N_{p}}$")
mp.legend(loc="upper left")
mp.savefig('expMoveSpeedUpSharedMem.pdf',bbox_inches = "tight")

print('exp',' ',m,' ',c)

####################Triplet Time############################

timeTrip = ([0.7720008,0.3919833,0.1928831,0.1282216,0.0969206,0.0775794,0.0652192,0.0545406,0.0489311,0.0420128,0.038765,0.0354619,0.0313368,0.0295397,0.02860998,0.027179])
tZeroTrip = timeTrip[0]
for i in range(0,len(timeTrip)):
  timeTrip[i] = tZeroTrip/timeTrip[i]

m,c = np.polyfit(nProc,timeTrip,1)

fig2=mp.figure()
ax2=fig2.add_subplot(111)
ax2.scatter(nProc,timeTrip,s=75,c="g",marker="o",label=r"t$_1$ = "+str(tZeroExp)+" s")
mp.grid()
mp.xlabel(r"$N_{p}}$")
mp.ylabel(r"t$_1$ / t$_{N_{p}}$")
mp.legend(loc="upper left")
mp.savefig('tripSumSpeedUpSharedMem.pdf',bbox_inches = "tight")

print('trip',' ',m,' ',c)

####################Partial Sum Time############################

timeSum = ([0.0853331,0.0428485,0.0213259,0.0142429,0.0107096,0.0085474,0.0071187,0.0060843,0.0053480,0.0047499,0.004288,0.0039144,0.0035850,0.0032802,0.00305217,0.002859])
tZeroSum = timeSum[0]
for i in range(0,len(timeTrip)):
  timeSum[i] = tZeroSum/timeSum[i]

m,c = np.polyfit(nProc,timeSum,1)

fig3=mp.figure()
ax3=fig3.add_subplot(111)
ax3.scatter(nProc,timeSum,s=75,c="g",marker="o",label=r"t$_1$ = "+str(tZeroExp)+" s")
mp.grid()
mp.xlabel(r"$N_{p}}$")
mp.ylabel(r"t$_1$ / t$_{N_{p}}$")
mp.legend(loc="upper left")
mp.savefig('partSumSpeedUpSharedMem.pdf',bbox_inches = "tight")

print('sum',' ',m,' ',c)

####################Accept Time############################

timeAcc = ([0.0006326,0.0003212,0.0001527,0.0001147,0.0000913,0.0000625,0.0000613,0.0000482,0.0000439,0.0000343,0.0000360,0.0000311,0.0000334,0.0000292,0.0000243,0.0000336])
tZeroAcc = timeAcc[0]
#timeAcc = timeAcc[0:15]
for i in range(0,len(timeAcc)):
  timeAcc[i] = tZeroAcc/timeAcc[i]
nProcAcc = nProc[0:15]

m,c = np.polyfit(nProc,timeAcc,1) # For true parameters, change nProc to nProcAcc and uncomment line above

fig4=mp.figure()
ax4=fig4.add_subplot(111)
ax4.scatter(nProc,timeAcc,s=75,c="g",marker="o",label=r"t$_1$ = "+str(tZeroExp)+" s")
mp.grid()
mp.xlabel(r"$N_{p}}$")
mp.ylabel(r"t$_1$ / t$_{N_{p}}$")
mp.legend(loc="upper left")
mp.savefig('acceptSpeedUpSharedMem.pdf',bbox_inches = "tight")

print('acc',' ',m,' ',c)

####################Reject Time############################

timeRej = ([0.0011379,0.0005794,0.0002726,0.0001775,0.0001468,0.0001118,0.0001080,0.0000786,0.0000815,0.0000689,0.0000653,0.0000605,0.0000559,0.0000517,0.0000514,0.0000467])
tZeroRej = timeRej[0]
for i in range(0,len(timeTrip)):
  timeRej[i] = tZeroRej/timeRej[i]

m,c = np.polyfit(nProc,timeRej,1)

fig5=mp.figure()
ax5=fig5.add_subplot(111)
ax5.scatter(nProc,timeRej,s=75,c="g",marker="o",label=r"t$_1$ = "+str(tZeroExp)+" s")
mp.grid()
mp.xlabel(r"$N_{p}}$")
mp.ylabel(r"t$_1$ / t$_{N_{p}}$")
mp.legend(loc="upper left")
mp.savefig('rejectSpeedUpSharedMem.pdf',bbox_inches = "tight")

print('rej',' ',m,' ',c)
