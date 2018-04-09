import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from scipy.optimize import curve_fit


#Load data for variables from text files
mhalo=np.loadtxt('mhalo.txt') /0.7                  #10^10 Mo/h
mstar=np.loadtxt('mstar.txt') /0.7                  #10^10 Mo/h
'''
f_half=np.loadtxt('t50.txt')/1000.                  #f1/2
red_half=np.loadtxt('red_half.txt')                 #Lookback time
StellarMass=np.loadtxt('StellarMass.txt')
StellarMassHist=np.loadtxt('StellarMassHist.txt')
star_gas_angle=np.loadtxt('star_gas_angle')
SfrDisk=np.loadtxt('SfrDisk.txt')                   #Mo / yr
SfrBulge=np.loadtxt('SfrBulge.txt')
SFR=SfrDisk+SfrBulge

#Initialise sSFR arrays
sSFR=np.zeros(19956)
sSFRDisk=np.zeros(19956)
sSFRBulge=np.zeros(19956)
'''
#Create condition that stellar mass must be greater than 0.1 solar masses
class_ok = np.where((mstar > 0.1))
'''
#Populate sSFR arrays
for i in range(0,(len(StellarMass[class_ok])-1)):
    if StellarMass[i] > 0:
        sSFR[i]=SFR[i]/StellarMass[class_ok][i]
        sSFRDisk[i]=SfrDisk[i]/StellarMass[class_ok][i]
        sSFRBulge[i]=SfrBulge[i]/StellarMass[class_ok][i]
    i=i+1
'''
###############################################################################

ms_mh = ((mstar[class_ok]/mhalo[class_ok]))

#Log of halo mass
log_halo=np.log10(mhalo)+10.
log_mhalo=np.log10(mhalo[class_ok])+10.

bins=10
mass_bin = np.linspace(min(log_mhalo),max(log_mhalo), bins)
m1=[]
m2=[]
m3=[]
m4=[]
x=mass_bin[0:(len(mass_bin)-1)]
y1=[]
y2=[]
y3=[]
y4=[]
y5=[]
y6=[]
y7=[]
y8=[]
y9=[]
temp_y=np.zeros((1000))

for i in range (0,len(log_mhalo)):
    for j in range (0,len(mass_bin)):
        j=0
        if log_mhalo[i]>mass_bin[j] and log_mhalo[i]<mass_bin[j+1]:
            #m1.append(log_mhalo[i])
            y1.append(ms_mh[i])
            #temp_y[j]=ms_mh[i]
            #print(np.append(temp_y,ms_mh[i]))
        j=j+1
        if log_mhalo[i]>mass_bin[j] and log_mhalo[i]<mass_bin[j+1]:
            #m2.append(log_mhalo[i]) 
            y2.append(ms_mh[i])
        j=j+1
        if log_mhalo[i]>mass_bin[j] and log_mhalo[i]<mass_bin[j+1]:
            #m3.append(log_mhalo[i])    
            y3.append(ms_mh[i])
        j=j+1
        if log_mhalo[i]>mass_bin[j] and log_mhalo[i]<mass_bin[j+1]:
            #m4.append(log_mhalo[i])    
            y4.append(ms_mh[i])        
        j=j+1
        if log_mhalo[i]>mass_bin[j] and log_mhalo[i]<mass_bin[j+1]:
            #m4.append(log_mhalo[i])    
            y5.append(ms_mh[i])        
        j=j+1
        if log_mhalo[i]>mass_bin[j] and log_mhalo[i]<mass_bin[j+1]:
            #m4.append(log_mhalo[i])    
            y6.append(ms_mh[i])        
        j=j+1
        if log_mhalo[i]>mass_bin[j] and log_mhalo[i]<mass_bin[j+1]:
            #m4.append(log_mhalo[i])    
            y7.append(ms_mh[i])        
        j=j+1
        if log_mhalo[i]>mass_bin[j] and log_mhalo[i]<mass_bin[j+1]:
            #m4.append(log_mhalo[i])    
            y8.append(ms_mh[i])        
        j=j+1
        if log_mhalo[i]>mass_bin[j] and log_mhalo[i]<mass_bin[j+1]:
            #m4.append(log_mhalo[i])    
            y9.append(ms_mh[i])

y1mean=np.mean(y1)
y2mean=np.mean(y2)
y3mean=np.mean(y3)
y4mean=np.mean(y4)
y5mean=np.mean(y5)
y6mean=np.mean(y6)
y7mean=np.mean(y7)
y8mean=np.mean(y8)
y9mean=np.mean(y9)


y=[]
y.append(y1mean)
y.append(y2mean)
y.append(y3mean)
y.append(y4mean)
y.append(y5mean)
y.append(y6mean)
y.append(y7mean)
y.append(y8mean)
y.append(y9mean)

plt.plot(x,y,'r')
#Plot of Mstar/Mhalo vs Mhalo
plt.scatter(log_mhalo,ms_mh,s=1 , label='Total Mass')
plt.title("Mstar/Mhalo vs Halo Mass")
plt.xlabel("Halo Mass (log(M/$h^-1$ Msun))")
plt.ylabel("Mstar/Mhalo")
plt.legend()
###############################################################################
'''#f half
f_half=f_half[class_ok]

age_bin = np.linspace(min(f_half),max(f_half), 30)
mean=[]
i=0
while i < len(f_half):
    mean.append(np.mean(f_half[i:i+100]))
    i=i+100
'''
###############################################################################
'''
#Use data to create the mass-weighted age
file=np.loadtxt('Millennium-output-times.dati')
bins=file[:,0]
t=file[:,4]
mwa=np.zeros(7245)              #Using 7245 due to class_ok making array smaller 
top=np.zeros((7245,64))
bottom=np.zeros((7245,64))
i=0

for i in range(0,(len(bins-1))):

    top[:,i]=(t[i]*StellarMassHist[class_ok][:,i])
    bottom[:,i]=(StellarMassHist[class_ok][:,i])

topsum=np.sum(top,axis=1)
botsum=np.sum(bottom,axis=1)

mwa=(topsum/botsum)/1e3
'''
###############################################################################
'''
#T50: the time at which the central galaxy had formed 50% of its mass
t50=np.zeros(7245)
m50=StellarMass[class_ok]/2
while i<63: 
    half=StellarMassHist[class_ok][:,i] 
    for j in range (0,(len(half)-1)):     
        if half[j]==m50[j]:         
            t50[j]=red_half[j] 
            print("ok") 
    i=i+1 
'''
###############################################################################
'''
#Plot of star gas angle vs mhalo
plt.scatter(log_halo,star_gas_angle)
plt.title("Star gas angle vs Halo Mass")
plt.xlabel("Halo Mass (log(M/$h^-1$ Msun))")
plt.ylabel("Star gas angle")
'''
###############################################################################
###############################################################################
'''
#Plot of t50 vs Mhalo
plt.scatter(log_mhalo,StellarMassHist[class_ok][:,32])
plt.title("t50 vs Halo Mass")
plt.xlabel("Halo Mass (log(M/$h^-1$ Msun))")
plt.ylabel("t50")
'''
###############################################################################
###############################################################################
'''
#Plot of Mass Weighted Age vs Mhalo
plt.scatter(log_mhalo,mwa, s=1)
plt.title("MWA vs Halo Mass")
plt.xlabel("Halo Mass (log(M/$h^-1$ Msun))")
plt.ylabel("MWA (Gyr)")

#Plot of Mass Weighted Average vs f half
plt.scatter(f_half[class_ok],mwa, s=1)
plt.title("MWA vs f 1/2")
plt.xlabel("f 1/2 (Gyr)")
plt.ylabel("MWA (Gyr)")
plt.ylim(-1,10)
'''
###############################################################################
###############################################################################

###############################################################################
'''
#Plot of Mstar/Mhalo vs f 1/2
plt.scatter(f_half[0:(len(m1))],ms_mh[0:(len(m1))],s=1 )
plt.scatter(f_half[0:(len(m2))],ms_mh[0:(len(m2))],s=1 )
plt.scatter(f_half[0:(len(m3))],ms_mh[0:(len(m3))],s=1 )
plt.scatter(f_half[0:(len(m4))],ms_mh[0:(len(m4))],s=1 )
plt.title("Mstar/Mhalo vs f 1/2")
plt.xlabel("f 1/2 (Gyr)")
plt.ylabel("Mstar/Mhalo")
'''
###############################################################################
###############################################################################
'''
#Plotting SFR vs Mhalo(log)

#Hexbin plot of SFR
fig, axs = plt.subplots(ncols=3, sharey=True)
ax = axs[0]
hb = ax.hexbin(log_mhalo, SFR[class_ok], gridsize=[10,10], bins='log', cmap='inferno')
cb = fig.colorbar(hb, ax=ax)
ax.set_title("SFR vs Halo Mass")
xmin=min(log_mhalo)
xmax=max(log_mhalo)
ymin=min(SFR)
ymax=60#max(SFR)
ax.axis([xmin, xmax, ymin, ymax])

#Plot of SFR
ax=axs[0]
ax.scatter(log_mhalo, SFR[class_ok], s=1 )
ax.set_title("SFR vs Halo Mass")
ax.set_xlabel("Halo Mass (log(M/$h^-1$ Msun))")
ax.set_ylabel("SFR (Msun/yr)")
ax.axis([xmin, xmax, ymin, ymax])

#Plot of SFR in Disk
ax=axs[1]
ax.scatter(log_mhalo, SfrDisk[class_ok], s=1 )
ax.set_title("SFR Disk vs Halo Mass")
ax.set_xlabel("Halo Mass (log(M/$h^-1$ Msun))")
ax.set_ylabel("SFR Disk (Msun/yr)")
ymin=min(SfrDisk)
ymax=60#max(SfrDisk)
ax.axis([xmin, xmax, ymin, ymax])

#Plot of SFR in Bulge
ax=axs[2]
ax.scatter(log_mhalo, SfrBulge[class_ok], s=1 )
ax.set_title("SFR Bulge vs Halo Mass")
ax.set_xlabel("Halo Mass (log(M/$h^-1$ Msun))")
ax.set_ylabel("SFR Bulge (Msun/yr)")
ax.axis([xmin, xmax, ymin, ymax])
ymin=min(SfrBulge)
ymax=60#max(SfrBulge)
'''
###############################################################################
###############################################################################
'''
#Plotting SFR versus F 1/2

#Set up subplots
fig, axs = plt.subplots(ncols=2, sharey=True)
ax = axs[0]
hb = ax.hexbin(f_half[class_ok], SFR[class_ok], gridsize=40, bins='log', cmap='inferno')
ax = axs[0]
ax.set_title("SFR vs f 1/2")
xmin=min(f_half)
xmax=max(f_half)
ymin=min(SFR)
ymax=max(SFR)

#Plot of SFR
ax=axs[1]
ax.scatter(f_half[class_ok], SFR[class_ok], s=1 )
ax.set_title("SFR vs f 1/2")
ax.set_xlabel("f 1/2 (Gyr)")
ax.set_ylabel("SFR (Msun/yr)")
ax.axis([xmin, xmax, ymin, ymax])

#Plot of SFR in Disk
ax=axs[1]
ax.scatter(f_half, SfrDisk, s=1 )
ax.set_title("SfrDisk vs f 1/2")
ax.set_xlabel("f 1/2 (Gyr)")
ax.set_ylabel("SfrDisk (Msun/yr)")
ymin=min(SfrDisk)
ymax=max(SfrDisk)
ax.axis([xmin, xmax, ymin, ymax])

#Plot of SFR in Bulge
ax=axs[2]
ax.scatter(f_half, SfrBulge, s=1 )
ax.set_title("SfrBulge vs f 1/2")
ax.set_xlabel("f 1/2 (Gyr)")
ax.set_ylabel("SfrBulge (Msun/yr)")
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
hb = ax.hexbin(log_mhalo, sSFR[class_ok], gridsize=[15,40], bins='log',cmap='inferno')#bins='log', 
cb = fig.colorbar(hb, ax=ax)
ax.set_title("sSFR vs Halo Mass")
xmin=min(log_mhalo)
xmax=max(log_mhalo)
ymin=min(sSFR)
ymax=60#max(sSFR)

#Plot of sSFR
ax=axs[0]
ax.scatter(log_mhalo, sSFR[class_ok],s=1 )
ax.set_title("sSFR vs Halo Mass")
ax.set_xlabel("Halo Mass (log(M/$h^-1$ Msun))")
ax.set_ylabel("sSFR ($yr^-1$)")
ax.axis([xmin, xmax, ymin, ymax])

#Plot of sSFR in Disk
ax=axs[1]
ax.scatter(log_mhalo, sSFRDisk[class_ok],s=1 )
ax.set_title("sSFR Disk vs Halo Mass")
ax.set_xlabel("Halo Mass (log(M/$h^-1$ Msun))")
ax.set_ylabel("sSFR Disk ($yr^-1$)")
ymin=min(sSFRDisk)
ymax=60#max(sSFRDisk)
ax.axis([xmin, xmax, ymin, ymax])

#Plot of sSFR in Bulge
ax=axs[2]
ax.scatter(log_mhalo, sSFRBulge[class_ok], s=1 )
ax.set_title("sSFR Bulge vs Halo Mass")
ax.set_xlabel("Halo Mass (log(M/$h^-1$ Msun))")
ax.set_ylabel("sSFR Bulge ($yr^-1$)")
ymin=min(sSFRBulge)
ymax=60#max(sSFRBulge)
ax.axis([xmin, xmax, ymin, ymax])
'''
###############################################################################
###############################################################################
'''
#Plotting sSFR versus F 1/2

#Set up subplots
fig, axs = plt.subplots(ncols=3, sharey=True)
#ax = axs[0]
#hb = ax.hexbin(f_half[class_ok], sSFR[class_ok], gridsize=40, bins='log', cmap='inferno')
#ax.set_title("sSFR vs f 1/2")
xmin=min(f_half)
xmax=max(f_half)
ymin=min(sSFR)
ymax=max(sSFR)

#Plot of sSFR
ax=axs[0]
ax.scatter(f_half[class_ok], sSFR[class_ok], s=1 )
ax.set_title("sSFR vs f 1/2")
ax.set_xlabel("f 1/2 (Gyr)")
ax.set_ylabel("sSFR")
ax.axis([xmin, xmax, ymin, ymax])

#Plot of sSFR in Disk
ax=axs[1]
ax.scatter(f_half[class_ok], sSFRDisk[class_ok], s=1 )
ax.set_title("sSfrDisk vs f 1/2")
ax.set_xlabel("f 1/2 (Gyr)")
ax.set_ylabel("sSfrDisk")
ymin=min(sSFRDisk)
ymax=max(sSFRDisk)
ax.axis([xmin, xmax, ymin, ymax])

#Plot of sSFR in Bulge
ax=axs[2]
ax.scatter(f_half[class_ok], sSFRBulge[class_ok], s=1 )
ax.set_title("sSfrBulge vs f 1/2")
ax.set_xlabel("f 1/2 (Gyr)")
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
