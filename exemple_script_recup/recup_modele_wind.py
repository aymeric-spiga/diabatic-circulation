import pylab as py
import numpy as np
import ppcompute
import matplotlib.pyplot as plt
import numpy.matlib
from ppclass import pp
from numpy import loadtxt

### Petite fonction utilisee plus tard pour convertir les tableaux DYNAMICO dans la taille de ceux radiatif (Cherche l indice de la valeur lat la plus proche de celle prise)
def find_x(array,x):
    y = int(np.argmin(np.abs(np.ceil(array[None].T - x)),axis=0))
    return int(y)


###Initialisation du domaine, des profils vertical et latitudinal + caracteristiques planete

#Recuperation des coefficients pour calculer le champ de pression 2D
coeffbp = loadtxt("coeffbp.txt", delimiter=",")

fff = "../simulations_dynamiques/Xhistins_195.nc"
time_step=300#"761040,500000000" #309 #"420,500" "360,420"

### VARIABLES DYNAMICO
lsub = pp(file=fff,var="ls",x=0,y=0,verbose=True).getf()
ls=(lsub/np.pi)*180.
u,lon,lat,zmod,t = pp(file=fff,var="u",t=time_step,x="-180,180",verbose=True).getfd()
v,lon,lat,zmod,t = pp(file=fff,var="v",t=time_step,x="-180,180",verbose=True).getfd()
temp,lon,lat,zmod,t = pp(file=fff,var="temperature",t=time_step,x="-180,180",verbose=True).getfd()
umod=pp(file=fff,var="u",t=time_step,x="-180,180",verbose=True).get()
vmod=pp(file=fff,var="v",t=time_step,x="-180,180",verbose=True).get()
psmod=pp(file=fff,var="ps",t=time_step,x="-180,180",verbose=True).get()


fffrad = "../SIMUS/test2d_saturne/diagfi.nc"
time_steprad="3540,3600"

### Variables modele radiatif
tetarad,lonrad,latrad,zmodrad,trad = pp(file=fffrad,var="teta",t=time_steprad,x="-180,180",verbose=True,changetime="correctls").getfd()



# ,logy=1 ,useindex="1000" ,changetime="correctls"

#Pour verifier que le ls est bien celui mis dans les plot (fixe)
print 'Interval de LS =', min(ls), max(ls)


#########  IMPORTANT  #########
#### Selectionner ici la tranche d'altitude a surveiller :
#Varie de 0 a 64 
infz=16    #borne inf de Z ###0   ### 16=100mb (base stratosphere)
topz=32    #borne sup de Z ###32   ### 32=3mb (haut de Dynamico)


nn=1   # Indice pour retirer les poles sur les latitudes
latrad=latrad[nn:-nn]
tetarad=tetarad[infz:topz,nn:-nn]

### ICI lat ne contient pas les poles, pas besoin de corriger.
phi = lat*np.pi/180.
sinphi = np.sin(phi)
cosphi = np.cos(phi)   #Transformation des latitudes en phi

zp=zmod

##### NOTA BENE: Il n y a pas de champ 2D d altitude dans DYNAMICO, il faut le faire a la main avec un fichier .txt fourni par Spiga et multiplie par pression de surface
#z2D=zmod

coeffbp2D=np.matlib.repmat(coeffbp,np.shape(lat)[0],1).T
p2D=coeffbp2D*psmod.f

#########  IMPORTANT  #########
###Definitions de toute les constantes planetaires

a=58210000.            #Rayon de la planete en m
g=10.                  #g acceleration
T0=135.                #Temperature au "sol"
cp=11500.              #Capacite calorifique
Rm=8.3144621           #Constante des gaz parfaits
M=2.34                 #Masse molaire gaz
H=(Rm*T0)/(g*M)*1000.  #Calcul hauteur caracteristique en m
Rs=(Rm/M)*1000.        #Constante specifique des gaz parfaits
ro=p2D/(Rs*temp)      #Calcul du champ de masse volumique (density)
p0=300000.             #Pression a la base de l'atmosphere


###Valeurs de P depuis l'equation hypsometrique en moyennant T a chaque point du profil vertical si altitude disponible mais pas P dans le fichier nc
### ATTENTION, cette methode est ici inutile, cf plus haut: NOTA BENE
#p2D=np.empty(np.shape(z2D))
#p=np.empty(np.shape(zp))
#p2D[0,:]=p0
#p[0]=p0
#for i in range(1,np.shape(p2D)[0]):
#for i in range(1,np.shape(p)[0]):
#  p[i]=p[i-1]/np.exp((g/(Rs*(np.mean(temp[i,:]+temp[i-1,:])/2.)))*(zp[i]-zp[i-1])*1000.)
  #for j in range(1,np.shape(p2D)[1]):
    #p2D[i,j]=p2D[i-1,j]/np.exp((g/(Rs*((temp[i,j]+temp[i-1,j])/2.)))*(z2D[i,j]-z2D[i-1,j])*1000.)


p=np.mean(p2D,1)
print p
p=p[infz:topz]
u=u[infz:topz,:]
v=v[infz:topz,:]
temp=temp[infz:topz,:]
zp=zmod[infz:topz]

### Convertir les tableaux du modele DYNAMICO dans la taille de ceux radiatif
ind=np.zeros(np.shape(latrad)[0])
for i in range(np.shape(latrad)[0]):
  ind[i]=find_x(lat,latrad[i])

ind=ind.astype(int)
#ind2D=np.matlib.repmat(ind,np.shape(zp)[0],1).astype(int)
temptest=temp.T[ind]
temptest=temptest.T
print p


plt.figure()
plt.contourf(lat,p,u,levels=np.linspace(-np.max(abs(u)),np.max(abs(u)),100),cmap=plt.cm.seismic)
plt.title('U Ls=340-355')
plt.colorbar()
plt.yscale('log')
plt.ylabel('Pressure (Pa)')
plt.xlabel('Latitude')
plt.gca().invert_yaxis()
#plt.savefig('V_final')

plt.figure()
plt.contourf(lat,p,v,levels=np.linspace(-np.max(abs(v)),np.max(abs(v)),100),cmap=plt.cm.seismic)
plt.title('V Ls=340-355')
plt.colorbar()
plt.yscale('log')
plt.ylabel('Altitude km ')
plt.xlabel('Latitude')
plt.gca().invert_yaxis()
#plt.savefig('V_final')

plt.figure()
plt.contourf(lat,zp,temp,levels=np.linspace(np.min(temp),np.max(temp),100))
plt.colorbar()
CS=plt.contour(lat,zp,temp,levels=np.round(np.linspace(80,160,9)),colors='black')
plt.clabel(CS, inline=0.1, fontsize=10,fmt='% 3d')
plt.title('Temperature (K) Ls=340-355') 
plt.ylabel('Altitude km ')
plt.xlabel('Latitude')
#plt.yscale('log')

plt.figure()
plt.contourf(lat,p,temp,levels=np.linspace(np.min(temp),np.max(temp),100))
plt.colorbar()
CS=plt.contour(lat,p,temp,levels=np.round(np.linspace(80,160,9)),colors='black')
plt.clabel(CS, inline=0.1, fontsize=10,fmt='% 3d')
plt.title('Temperature (K) Ls=340-355' ) 
plt.ylabel('Pressure (Pa)')
plt.xlabel('Latitude')
plt.yscale('log')
plt.gca().invert_yaxis()

plt.figure()
plt.contourf(latrad,p,temptest,levels=np.linspace(np.min(temptest),np.max(temptest),100))
plt.colorbar()
CS=plt.contour(latrad,p,temptest,levels=np.round(np.linspace(80,160,9)),colors='black')
plt.clabel(CS, inline=0.1, fontsize=10,fmt='% 3d')
plt.title('Temperature (K) Ls=340-355' ) 
plt.ylabel('Pressure (Pa)')
plt.xlabel('Latitude')
plt.yscale('log')
plt.gca().invert_yaxis()



plt.show(block=False)
input("Hit it and quit it! (press enter to close all)")
plt.close()
