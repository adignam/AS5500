import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

#Load data for variables from text files
mhalo=np.loadtxt('mhalo.txt') /0.7                  #10^10 Mo/h
mstar=np.loadtxt('mstar.txt') /0.7                  #10^10 Mo/h
f_half=np.loadtxt('t50.txt')/1000.                  #f1/2
red_half=np.loadtxt('red_half.txt')                 #Lookback time
StellarMass=np.loadtxt('StellarMass.txt')
StellarMassHist=np.loadtxt('StellarMassHist.txt')
SfrDisk=np.loadtxt('SfrDisk.txt')                   #Mo / yr
SfrBulge=np.loadtxt('SfrBulge.txt')
SFR=SfrDisk+SfrBulge

#Initialise sSFR arrays
sSFR=np.zeros(19956)
sSFRDisk=np.zeros(19956)
sSFRBulge=np.zeros(19956)

#Create condition that stellar mass must be greater than 0.1 solar masses
class_ok = np.where((mstar > 0.1))

#Populate sSFR arrays
for i in range(0,(len(StellarMass[class_ok])-1)):
    if StellarMass[i] > 0:
        sSFR[i]=SFR[i]/StellarMass[class_ok][i]
        sSFRDisk[i]=SfrDisk[i]/StellarMass[class_ok][i]
        sSFRBulge[i]=SfrBulge[i]/StellarMass[class_ok][i]
    i=i+1
###############################################################################

#Log of halo mass
log_mhalo=np.log10(mhalo)+10.

###############################################################################

#Use data to create the mass-weighted average
file=np.loadtxt('Millennium-output-times.dati')
bins=file[:,0]
t=file[:,4]
mwa=np.zeros(19956)     
top=np.zeros((19956,64))
bottom=np.zeros((19956,64))
i=0

for i in range(0,(len(bins-1))):

    top[:,i]=(t[i]*StellarMassHist[:,i])
    bottom[:,i]=(StellarMassHist[:,i])

topsum=np.sum(top,axis=0)
botsum=np.sum(bottom,axis=0)

for i in range (0,(len(topsum)-1)):
    #if botsum[i]>0:
    mwa[i] = topsum[i]/botsum[i]

###############################################################################
###############################################################################

#Plot of Mass Weighted Average vs Mhalo
plt.scatter(log_mhalo,(mwa/1000))
plt.title("MWA vs Mhalo")
plt.xlabel("Mhalo (log10)")
plt.ylabel("MWA")
plt.ylim(-1,15)

###############################################################################
###############################################################################
'''
#Plot of Mstar/Mhalo vs Mhalo
plt.scatter(log_mhalo,((mstar/mhalo)) )
plt.title("Mstar/Mhalo vs Mhalo")
plt.xlabel("Mhalo (log10)")
plt.ylabel("Mstar/Mhalo")
plt.xlim((11.5,14.5))
'''
###############################################################################
'''
#Plot of Mstar/Mhalo vs f 1/2
plt.scatter(f_half,((mstar/mhalo)) )
plt.title("Mstar/Mhalo vs f 1/2")
plt.xlabel("f 1/2")
plt.ylabel("Mstar/Mhalo")
'''
###############################################################################
###############################################################################
'''
#Plotting SFR vs Mhalo(log)

#Hexbin plot of SFR
fig, axs = plt.subplots(ncols=3, sharey=True)
ax = axs[0]
#hb = ax.hexbin(log_mhalo, SFR, gridsize=40, bins='log', cmap='inferno')
#cb = fig.colorbar(hb, ax=ax)
#ax.set_title("SFR vs Mhalo(log) using hexbin")
xmin=min(log_mhalo)
xmax=max(log_mhalo)
ymin=min(SFR)
ymax=max(SFR)
ax.axis([xmin, xmax, ymin, ymax])

#Plot of SFR
ax=axs[0]
ax.scatter(log_mhalo, SFR )
ax.set_title("SFR vs Mhalo(log)")
ax.set_xlabel("Mhalo")
ax.set_ylabel("SFR")
ax.axis([xmin, xmax, ymin, ymax])

#Plot of SFR in Disk
ax=axs[1]
ax.scatter(log_mhalo, SfrDisk )
ax.set_title("SfrDisk vs Mhalo(log)")
ax.set_xlabel("Mhalo")
ax.set_ylabel("SfrDisk")
ymin=min(SfrDisk)
ymax=max(SfrDisk)
ax.axis([xmin, xmax, ymin, ymax])

#Plot of SFR in Bulge
ax=axs[2]
ax.scatter(log_mhalo, SfrBulge )
ax.set_title("SfrBulge vs Mhalo(log)")
ax.set_xlabel("Mhalo")
ax.set_ylabel("SfrBulge")
ax.axis([xmin, xmax, ymin, ymax])
ymin=min(SfrBulge)
ymax=max(SfrBulge)
'''
###############################################################################
###############################################################################
'''
#Plotting SFR versus F 1/2

#Set up subplots
fig, axs = plt.subplots(ncols=3, sharey=True)
ax = axs[0]
xmin=min(f_half)
xmax=max(f_half)
ymin=min(SFR)
ymax=max(SFR)

#Plot of SFR
ax=axs[0]
ax.scatter(f_half, SFR )
ax.set_title("SFR vs f 1/2")
ax.set_xlabel("f 1/2")
ax.set_ylabel("SFR")
ax.axis([xmin, xmax, ymin, ymax])

#Plot of SFR in Disk
ax=axs[1]
ax.scatter(f_half, SfrDisk )
ax.set_title("SfrDisk vs f 1/2")
ax.set_xlabel("f 1/2")
ax.set_ylabel("SfrDisk")
ymin=min(SfrDisk)
ymax=max(SfrDisk)
ax.axis([xmin, xmax, ymin, ymax])

#Plot of SFR in Bulge
ax=axs[2]
ax.scatter(f_half, SfrBulge )
ax.set_title("SfrBulge vs f 1/2")
ax.set_xlabel("f 1/2")
ax.set_ylabel("SfrBulge")
ymin=min(SfrBulge)
ymax=max(SfrBulge)
ax.axis([xmin, xmax, ymin, ymax])
'''
###############################################################################
###############################################################################
'''
#Plotting sSFR versus Mhalo

#Set up subplots
fig, axs = plt.subplots(ncols=3, sharey=True)
ax = axs[0]
xmin=min(log_mhalo)
xmax=max(log_mhalo)
ymin=min(sSFR)
ymax=max(sSFR)

#Plot of sSFR
ax=axs[0]
ax.scatter(log_mhalo, sSFR )
ax.set_title("sSFR vs Mhalo")
ax.set_xlabel("Mhalo")
ax.set_ylabel("sSFR")
ax.axis([xmin, xmax, ymin, ymax])

#Plot of sSFR in Disk
ax=axs[1]
ax.scatter(log_mhalo, sSFRDisk )
ax.set_title("sSfrDisk vs Mhalo")
ax.set_xlabel("Mhalo")
ax.set_ylabel("sSfrDisk")
ymin=min(sSFRDisk)
ymax=max(sSFRDisk)
ax.axis([xmin, xmax, ymin, ymax])

#Plot of sSFR in Bulge
ax=axs[2]
ax.scatter(log_mhalo, sSFRBulge )
ax.set_title("sSfrBulge vs Mhalo")
ax.set_xlabel("Mhalo")
ax.set_ylabel("sSfrBulge")
ymin=min(sSFRBulge)
ymax=max(sSFRBulge)
ax.axis([xmin, xmax, ymin, ymax])
'''
###############################################################################
###############################################################################
'''
#Plotting sSFR versus F 1/2

#Set up subplots
fig, axs = plt.subplots(ncols=3, sharey=True)
ax = axs[0]
xmin=min(f_half)
xmax=max(f_half)
ymin=min(sSFR)
ymax=max(sSFR)

#Plot of sSFR
ax=axs[0]
ax.scatter(f_half, sSFR )
ax.set_title("sSFR vs f 1/2")
ax.set_xlabel("f 1/2")
ax.set_ylabel("sSFR")
ax.axis([xmin, xmax, ymin, ymax])

#Plot of sSFR in Disk
ax=axs[1]
ax.scatter(f_half, sSFRDisk )
ax.set_title("sSfrDisk vs f 1/2")
ax.set_xlabel("f 1/2")
ax.set_ylabel("sSfrDisk")
ymin=min(sSFRDisk)
ymax=max(sSFRDisk)
ax.axis([xmin, xmax, ymin, ymax])

#Plot of sSFR in Bulge
ax=axs[2]
ax.scatter(f_half, sSFRBulge )
ax.set_title("sSfrBulge vs f 1/2")
ax.set_xlabel("f 1/2")
ax.set_ylabel("sSfrBulge")
ymin=min(sSFRBulge)
ymax=max(sSFRBulge)
ax.axis([xmin, xmax, ymin, ymax])
'''
###############################################################################

###############################################################################

#Spearman rank correlation coefficient
#rs=spearmanr(log_mhalo, (mstar/mhalo))  # rs=0.816150642461418
#rs=spearmanr(log_mhalo, SFR)            # rs=0.42876820663520371
#rs=spearmanr(log_mhalo, sSFR)           # rs=0.12956006509877613
#print(rs)
#Add in theta as a feature