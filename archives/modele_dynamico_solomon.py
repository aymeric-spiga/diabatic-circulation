import pylab as py
import numpy as np
import ppcompute
import matplotlib.pyplot as plt
import numpy.matlib
from ppclass import pp
from numpy import loadtxt
from scipy import integrate

### Petite fonction utilisee plus tard pour convertir les tableaux DYNAMICO dans la taille de ceux radiatif (Cherche l indice de la valeur lat la plus proche de celle prise)
def find_x(array,x):
    y = int(np.argmin(np.abs(np.ceil(array[None].T - x)),axis=0))
    return int(y)


###Initialisation du domaine, des profils vertical et latitudinal + caracteristiques planete

#Recuperation des coefficients pour calculer le champ de pression 2D
coeffbp = loadtxt("../coeffbp.txt", delimiter=",")

fff = "../../simulations_dynamiques/Xhistins_195.nc"
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


fffrad = "../../simulations_nodyn/diagfi_saturn_rings.nc"
time_steprad="3570,3600"

### Variables modele radiatif
tetarad,lonrad,latrad,zmodrad,trad = pp(file=fffrad,var="teta",t=time_steprad,x="-180,180",verbose=True,changetime="correctls").getfd()
zdtsw= pp(file=fffrad,var="zdtsw",t=time_steprad,x="-180,180",verbose=True,changetime="correctls").get()
zdtlw= pp(file=fffrad,var="zdtlw",t=time_steprad,x="-180,180",verbose=True,changetime="correctls").get()


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

p2D=p2D[infz:topz,:]
temp=temp[infz:topz,:]
ro=ro[infz:topz,:]
p=np.mean(p2D,1)
print p
#p=p[infz:topz]
u=u[infz:topz,:]
v=v[infz:topz,:]
zp=zmod[infz:topz]
z=zp
Qrad=zdtsw+zdtlw
Q=Qrad.f[infz:topz,nn:-nn]
phi = latrad*np.pi/180.
sinphi = np.sin(phi)
cosphi = np.cos(phi)

### Convertir les tableaux du modele DYNAMICO dans la taille de ceux radiatif
ind=np.zeros(np.shape(latrad)[0])
for i in range(np.shape(latrad)[0]):
  ind[i]=find_x(lat,latrad[i])

ind=ind.astype(int)
#ind2D=np.matlib.repmat(ind,np.shape(zp)[0],1).astype(int)
tempind=temp.T[ind]
temp=tempind.T
p2Dind=p2D.T[ind]
p2D=p2Dind.T
roind=ro.T[ind]
ro=roind.T

(nz,nx)=np.shape(temp)
teta=(nz,nx)
teta=np.empty(teta)
teta=temp*((p0/p2D)**(Rs/cp)) 

print p

###Calcul des derivees de teta et temp en phi,z et p

dtetadz,dtetadphi=np.gradient(teta,z*1000.,phi, edge_order=2)
#dtetadphi,dtetadz = ppcompute.deriv2d(teta,phi,z*1000)
dtetadp,dtetadphi=np.gradient(teta,p,phi, edge_order=2)#deriv(teta,phi,p)
#dtetadphi,dtetadp = ppcompute.deriv2d(teta,phi,p)
dtempdz,dtempdphi=np.gradient(temp,z*1000.,phi, edge_order=2)
#dtempdphi,dtempdz = ppcompute.deriv2d(temp,phi,z*1000.)


###Initialisation de w(omega) premiere iteration avec v=0
wini=np.empty((nz,nx))
winisave=np.empty((nz,nx))

###3 equations qui calcul W avec 3 differentes methodes pour N2 (frequence de brunt vaisala)

#equations qui convergent le mieux
#Methode 1
wini=Q[:,:]/(((H/Rs)*g/teta*dtetadz)+dtempdz)

#Methode 2
#wini=Q[:,:]/((H/Rs)*((-ro*g**2)/teta)*dtetadp+dtempdz)

#1 equation equivalente qui diverge
#wini=Q[:,:]/((H/Rs)*((g/T0)*((-g/cp)+dtempdz))+dtempdz)

#Equation qui calcul omega
#mask = np.where(dtetadp != 0)
#wini[mask]=Q[mask]/((temp[mask]/teta[mask])*dtetadp[mask])


#Sauvegarde wini pour le plotter avant application du terme correctif + Initialisation de termes
winisave[:,:]=wini[:,:]
vini=np.empty((nz,nx))
psy=np.empty((nz,nx))
epsy=np.empty((nz,1))

###Initialisation de v depuis wini(v=0) : Correspond a la premiere iteration sachant que epsy= terme correctif sur l integral pour faire tendre la vitesse du pole oppose a 0

for i in range(nz):
  #Integrale avec methode simpson
  numer = integrate.simps(wini[i,:]*cosphi[0:nx+1],phi[0:nx+1],axis=0)
  denom = integrate.simps(cosphi[0:nx+1],phi[0:nx+1],axis=0)
  epsy[i] = -numer/denom
  #Integrale avec methode trapezoide (insuffisante, NaN sur certaines donnees)
  #numer = -np.trapz(wini[i,:]*cosphi[0:nx+1],phi[0:nx+1],axis=0)
  #denom = np.trapz(cosphi[0:nx+1],phi[0:nx+1],axis=0)
  for j in range(nx):   
    wini[i,j]=wini[i,j]+epsy[i]
    integrand = (wini[i,0:j+1])*cosphi[0:j+1]*a
    psy[i,j]=integrate.simps(integrand,phi[0:j+1],axis=0)
    #psy[i,j]=np.trapz(integrand,phi[0:j+1],axis=0)


#Calcul de V a partir de W
dropsydz,dropsydphi = np.gradient(ro*psy,z*1000.,phi, edge_order=2)
#dropsydphi,dropsydz = ppcompute.deriv2d(ro*psy,phi,z*1000.)
vini=-(1./cosphi)*dropsydz*1./ro


#Variables necessaires a la boucle + Initialisation des variables
h=1 #=nombre d iterations, on a deja "une" iteration en v=0
v_next = np.empty((nz,nx))
w_next = np.empty((nz,nx))
w_prev = np.empty((nz,nx))
eps = np.empty(nz)
v_next[:,:] = vini[:,:]
w_prev[:,:] = wini[:,:]
psyini=np.empty((nz,nx))
psyini[:,:]=psy[:,:]

# !!!!!!!!! Choix de la precision et du pourcent a comparer dans la boucle
percent = 1.
precision = 100000.

###Calcul V et W tant que cela converge suivant la precision et le percent

#####DEBUT Boucle While #####
#while precision >= percent/100.:
for h in range(2,4):  #Pour controler l iteration (decommentez h=h+1 dans la boucle)

  ###3 equations qui calcul W avec 3 differentes methodes pour N2 (frequence de brunt vaisala)
  
  #equations qui convergent le mieux
  #Methode 1
  w_next=(Q[:,:]-(v_next[:,:]*dtempdphi/a))/(((H/Rs)*g/teta*dtetadz)+dtempdz)
  #Methode 2
  #w_next=(Q[:,:]-(v_next[:,:]*dtempdphi/a))/((H/Rs)*((-ro*g**2)/teta)*dtetadp+dtempdz)

  #1 equation equivalente qui diverge
  #w_next=(Q[:,:]-(v_next[:,:]*dtempdphi/a))/((H/Rs)*((g/T0)*((-g/cp)+dtempdz))+dtempdz)

  #Equation qui calcul omega
  #w_next[mask]=(Q[mask]-(v_next[mask]*dtempdphi[mask]/a))/((temp[mask]/teta[mask])*dtetadp[mask])

  for i in range(nz):
    numer = -integrate.simps(w_next[i,:]*cosphi[0:nx+1],phi[0:nx+1],axis=0)
    denom = integrate.simps(cosphi[0:nx+1],phi[0:nx+1],axis=0)
    eps[i] = numer/denom
    for j in range(nx):   
      w_next[i,j]=w_next[i,j]+eps[i]
      integrand = (w_next[i,0:j+1])*cosphi[0:j+1]*a
      psy[i,j]=integrate.simps(integrand,phi[0:j+1],axis=0)
  w = np.where(np.abs(w_next[:,:])-np.abs(w_prev[:,:]) != 0)
  precision = np.mean(np.abs(w_next[w] - w_prev[w])/np.abs(0.5*(w_next[w]+w_prev[w])))
  print precision, np.max(w_next[w]), np.min(w_next[w])
  dropsydz,dropsydphi = np.gradient(ro*psy,z*1000.,phi, edge_order=2)
  #dropsydphi,dropsydz = ppcompute.deriv2d(ro*psy,phi,z*1000.)
  v_next=-(1./cosphi)*dropsydz*1./ro
  w_prev[:,:] = w_next[:,:]
  #h += 1

#####Fin Boucle While #####

#print precision*100.,percent
print('On a realise ',h,'iterations')


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
plt.contourf(latrad,p,temp,levels=np.linspace(np.min(temp),np.max(temp),100))
plt.colorbar()
CS=plt.contour(latrad,p,temp,levels=np.round(np.linspace(80,160,9)),colors='black')
plt.clabel(CS, inline=0.1, fontsize=10,fmt='% 3d')
plt.title('Temperature (K) Ls=340-355' ) 
plt.ylabel('Pressure (Pa)')
plt.xlabel('Latitude')
plt.yscale('log')
plt.gca().invert_yaxis()

plt.figure()
plt.contourf(latrad,p,teta,levels=np.linspace(np.min(teta),np.max(teta),100))
plt.colorbar()
plt.title('Teta (K) Ls=340-355' ) 
plt.ylabel('Pressure (Pa)')
plt.xlabel('Latitude')
plt.yscale('log')
plt.gca().invert_yaxis()

plt.figure()
plt.contourf(latrad,p,Q,levels=np.linspace(np.min(Q),np.max(Q),100))
plt.colorbar()
plt.title('Chauffage radiatif Ls=340-355' ) 
plt.ylabel('Pressure (Pa)')
plt.xlabel('Latitude')
plt.yscale('log')
plt.gca().invert_yaxis()

plt.figure()
plt.contourf(latrad,p,w_next,levels=np.linspace(-np.max(abs(w_next)),np.max(abs(w_next)),100),cmap=plt.cm.seismic)
plt.title('Final W (m.s-1) Ls=%s' %time_step) 
plt.colorbar()
plt.yscale('log')
plt.ylabel('Pressure (Pa)')
plt.xlabel('Latitude')
plt.gca().invert_yaxis()

plt.figure()
plt.contourf(latrad,p,v_next,levels=np.linspace(-np.max(abs(v_next)),np.max(abs(v_next)),100),cmap=plt.cm.seismic)
plt.title('Final V (m.s-1) Ls=%s' %time_step) 
plt.colorbar()
plt.yscale('log')
plt.ylabel('Pressure (Pa)')
plt.xlabel('Latitude')
plt.gca().invert_yaxis()
#plt.savefig('V_final')


plt.show(block=False)
input("Hit it and quit it! (press enter to close all)")
plt.close()
