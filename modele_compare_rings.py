import pylab as py
import numpy as np
#import ppcompute
import matplotlib.pyplot as plt
import numpy.matlib
from ppclass import pp
from scipy import integrate

###Emplacement du fichier, choix de la fenetre de temps, recuperation des variables puis choix de la grille (altitude, lat ..)

# -> Choix du fichier du modele radiatif/convectif
norings = "../simulations_nodyn/diagfi_saturn_norings.nc"
rings = "../simulations_nodyn/diagfi_saturn_rings.nc"
#temporaire.nc
#../simulations_nodyn/diagfi_saturn_rings.nc #rings
#../simulations_nodyn/diagfi_saturn_norings.nc #norings

# -> Choix du pas de temps
time_step="3540,3560"
#"3240,3600"=8e annee en entier (moyenne annuelle)
#"3360,3380"=134 EteHN 
#"3540,3560"=308 hiverHN
#"3570,3600"=histins195 (330-360) correspond a DYNAMICO


# -> Recuperation des variables dans les fichiers moyennes zonalement avec ls corrige

teta,lon,lat,zmod,t = pp(file=norings,var="teta",t=time_step,x="-180,180",verbose=True,changetime="correctls").getfd()
zdtsw= pp(file=norings,var="zdtsw",t=time_step,x="-180,180",verbose=True,changetime="correctls").get()
zdtlw= pp(file=norings,var="zdtlw",t=time_step,x="-180,180",verbose=True,changetime="correctls").get()
pmod= pp(file=norings,var="p",t=time_step,x="-180,180",verbose=True,changetime="correctls").get()
tempmod= pp(file=norings,var="temp",t=time_step,x="-180,180",verbose=True,changetime="correctls").get()
tetar,lon,lat,zmod,tr = pp(file=rings,var="teta",t=time_step,x="-180,180",verbose=True,changetime="correctls").getfd()
zdtswr= pp(file=rings,var="zdtsw",t=time_step,x="-180,180",verbose=True,changetime="correctls").get()
zdtlwr= pp(file=rings,var="zdtlw",t=time_step,x="-180,180",verbose=True,changetime="correctls").get()
pmodr= pp(file=rings,var="p",t=time_step,x="-180,180",verbose=True,changetime="correctls").get()
tempmodr= pp(file=rings,var="temp",t=time_step,x="-180,180",verbose=True,changetime="correctls").get()

#### Pour utiliser des LS: mettre ,changetime="correctls" et time_step=#Ls
#### Sinon ,useindex="1000" marche 



#########  IMPORTANT  #########
#### Selectionner ici la tranche d'altitude a surveiller :

##Pour plotter le maximum possible 16-56 <- ZONE CONSEILLEE
##Pour comparer a DYNAMICO 16-30
##Circulation bas strato 16-35
##Circulation haute strato et meso 16-56 ou 25-56

##10000Pa=100mbar=16 (base stratosphere)
##1000Pa=10mbar=25
##Sommet DYNAMICO a 30 (en realite 32 mais on neglige 2 points)
##100Pa=1mbar= 35
##10Pa=0.1mbar=45
##1Pa=0.01mbar=56  (haut stratosphere/mesosphere)

infz=16     #borne inf de Z ###0
topz=56     #borne sup de Z ###64

nn=1        #Indice pour retirer les poles sur les latitudes
lat=lat[nn:-nn]     #On retire les poles (sinon cosphi -> 0 et division par 0)
phi = lat*np.pi/180.
sinphi = np.sin(phi)
cosphi = np.cos(phi)   #Transformation des latitudes en phi

#Decoupe de l'altitude et creation du champ de taux radiatif
zp=zmod[infz:topz]
Qrad=zdtsw+zdtlw
Qradr=zdtswr+zdtlwr
Q=Qrad.f[infz:topz,nn:-nn]
Qr=Qradr.f[infz:topz,nn:-nn]


###Creation d'un champ de pression et temperature a partir de Z et teta
p2D=pmod.f[infz:topz,nn:-nn]
p2Dr=pmodr.f[infz:topz,nn:-nn]
p=p2D[:,0]
temp=tempmod.f[infz:topz,nn:-nn]
tempr=tempmodr.f[infz:topz,nn:-nn]
teta=teta[infz:topz,nn:-nn]
tetar=tetar[infz:topz,nn:-nn]


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
ro=p2D/(Rs*temp)       #Calcul du champ de masse volumique (density)
ror=p2Dr/(Rs*tempr)
p0=300000.             #Pression a la base de l'atmosphere

##### Si le fichier possede une sortie teta mais pas temp(ou inversement):
(nz,nx)=np.shape(teta)
#temp=(nz,nx)
#temp=np.empty(temp)
#temp=teta*((p2D/p0)**(Rs/cp)) 

#(nz,nx)=np.shape(temp)
#teta=(nz,nx)
#teta=np.empty(teta)
#teta=temp*((p0/p2D)**(Rs/cp)) 

###3 Methodes de calcul de z
#NOTA BENE -> La 3eme methode est plus precise mais elle pourrait etre amelioree
###1-Calcul de z en considerant Temperature constante
#z=-H*np.log(pmod.f[:,0]/p0)/1000.
###2-Valeurs de z dans le fichier du modele
#z=np.array([1.0,3.0,8.0,15.0,25.0,35.0,47.0,59.0,71.0,83.0,95.0,107.0,119.0,131.0,143.0,155.0,167.0,179.0,191.0,203.0,215.0,227.0,239.0,251.0,263.0,275.0,287.0,299.0,311.0,323.0,335.0,347.0,359.0,371.0,383.0,395.0,407.0,419.0,431.0,443.0,455.0,467.0,479.0,491.0,503.0,515.0,527.0,539.0,551.0,563.0,575.0,587.0,599.0,611.0,623.0,635.0,647.0,659.0,671.0,683.0,695.0,707.0,719.0,731.0])
###3-Valeurs de z depuis l'equation hypsometrique en moyennant T a chaque point du profil vertical
z=np.empty(np.shape(pmod.f)[0])
zr=np.empty(np.shape(pmodr.f)[0])
z[0]=-H*np.log(pmod.f[0,0]/p0)/1000.
zr[0]=-H*np.log(pmodr.f[0,0]/p0)/1000.
for i in range(1,np.shape(z)[0]):
  z[i]=z[i-1]+((Rs/g)*np.mean((tempmodr.f[i,:]+tempmodr.f[i-1,:])/2.)*np.log(pmodr.f[i-1,0]/pmodr.f[i,0]))/1000.
  zr[i]=zr[i-1]+((Rs/g)*np.mean((tempmodr.f[i,:]+tempmodr.f[i-1,:])/2.)*np.log(pmodr.f[i-1,0]/pmodr.f[i,0]))/1000.
z=z[infz:topz]
zr=zr[infz:topz]



###Valeurs de P depuis l'equation hypsometrique en moyennant T a chaque point du profil vertical si altitude disponible mais pas P dans le fichier nc
#p2D=np.empty(np.shape(z2D))
#p=np.empty(np.shape(z2D)[0])
#p2D[0,:]=p0
#for i in range(1,np.shape(p2D)[0]):
#  for j in range(1,np.shape(p2D)[1]):
#    p2D[i,j]=p2D[i-1,j]/np.exp((g/(Rs*((tempmod.f[i,j]+tempmod.f[i-1,j])/2.)))*(z2D[i,j]-z2D[i-1,j])*1000.)
#p=np.mean(p2D,1)
#p=p[infz:topz]



###Calcul des derivees de teta et temp en phi,z et p

dtetadz,dtetadphi=np.gradient(teta,z*1000.,phi, edge_order=2)
#dtetadphi,dtetadz = ppcompute.deriv2d(teta,phi,z*1000)
dtetadp,dtetadphi=np.gradient(teta,p,phi, edge_order=2)#deriv(teta,phi,p)
#dtetadphi,dtetadp = ppcompute.deriv2d(teta,phi,p)
dtempdz,dtempdphi=np.gradient(temp,z*1000.,phi, edge_order=2)
#dtempdphi,dtempdz = ppcompute.deriv2d(temp,phi,z*1000.)
dtetadzr,dtetadphir=np.gradient(tetar,zr*1000.,phi, edge_order=2)
dtetadpr,dtetadphir=np.gradient(tetar,p,phi, edge_order=2)
dtempdzr,dtempdphir=np.gradient(tempr,zr*1000.,phi, edge_order=2)



###Initialisation de w(omega) premiere iteration avec v=0
wini=np.empty((nz,nx))
winisave=np.empty((nz,nx))
winir=np.empty((nz,nx))
winisaver=np.empty((nz,nx))


###3 equations qui calcul W avec 3 differentes methodes pour N2 (frequence de brunt vaisala)


#equations qui convergent le mieux
#Methode 1 <- CONSEILLEE et utilisee pour le rapport
wini=Q[:,:]/(((H/Rs)*g/teta*dtetadz)+dtempdz)
winir=Qr[:,:]/(((H/Rs)*g/tetar*dtetadzr)+dtempdzr)
#Methode 2
#wini=Q[:,:]/((H/Rs)*((-ro*g**2)/teta)*dtetadp+dtempdz)
#winir=Qr[:,:]/((H/Rs)*((-ror*g**2)/tetar)*dtetadpr+dtempdzr)

#1 equation equivalente qui diverge
#wini=Q[:,:]/((H/Rs)*((g/T0)*((-g/cp)+dtempdz))+dtempdz)
#winir=Qr[:,:]/((H/Rs)*((g/T0)*((-g/cp)+dtempdzr))+dtempdzr)

#Equation qui calcul omega
#mask = np.where(dtetadp != 0)
#wini[mask]=Q[mask]/((temp[mask]/teta[mask])*dtetadp[mask])
#mask = np.where(dtetadpr != 0)
#winir[mask]=Qr[mask]/((tempr[mask]/tetar[mask])*dtetadpr[mask])


#Sauvegarde wini pour le plotter avant application du terme correctif + Initialisation de termes
winisave[:,:]=wini[:,:]
vini=np.empty((nz,nx))
psy=np.empty((nz,nx))
epsy=np.empty((nz,1))
winisaver[:,:]=winir[:,:]
vinir=np.empty((nz,nx))
psyr=np.empty((nz,nx))
epsyr=np.empty((nz,1))

###Initialisation de v depuis wini(v=0) : Correspond a la premiere iteration sachant que epsy= terme correctif sur l integral pour faire tendre la vitesse du pole oppose a 0

for i in range(nz):
  epsyr[i]=-(np.trapz(winir[i,:]*cosphi[0:nx+1],phi[0:nx+1],axis=0))/(np.trapz(cosphi[0:nx+1],phi[0:nx+1],axis=0))
  #Integrale avec methode simpson
  numer = integrate.simps(wini[i,:]*cosphi[0:nx+1],phi[0:nx+1],axis=0)
  denom = integrate.simps(cosphi[0:nx+1],phi[0:nx+1],axis=0)
  numerr = integrate.simps(winir[i,:]*cosphi[0:nx+1],phi[0:nx+1],axis=0)
  denomr = integrate.simps(cosphi[0:nx+1],phi[0:nx+1],axis=0)
  #Integrale avec methode trapezoide (insuffisante, NaN sur certaines donnees)
  #numer = np.trapz(wini[i,:]*cosphi[0:nx+1],phi[0:nx+1],axis=0)
  #denom = np.trapz(cosphi[0:nx+1],phi[0:nx+1],axis=0)
  #numerr = np.trapz(winir[i,:]*cosphi[0:nx+1],phi[0:nx+1],axis=0)
  #denomr = np.trapz(cosphi[0:nx+1],phi[0:nx+1],axis=0)
  epsy[i] = -numer/denom
  epsyr[i] = -numerr/denomr
  for j in range(nx):   
    wini[i,j]=wini[i,j]+epsy[i]
    winir[i,j]=winir[i,j]+epsyr[i]
    integrand = (wini[i,0:j+1])*cosphi[0:j+1]*a
    integrandr = (winir[i,0:j+1])*cosphi[0:j+1]*a
    psy[i,j]=integrate.simps(integrand,phi[0:j+1],axis=0)
    psyr[i,j]=integrate.simps(integrandr,phi[0:j+1],axis=0)
    #psy[i,j]=np.trapz(integrand,phi[0:j+1],axis=0)
    #psyr[i,j]=np.trapz(integrandr,phi[0:j+1],axis=0)


#Calcul de V a partir de W
dropsydz,dropsydphi = np.gradient(ro*psy,z*1000.,phi, edge_order=2)
dropsydzr,dropsydphir = np.gradient(ror*psyr,zr*1000.,phi, edge_order=2)
#dropsydphi,dropsydz = ppcompute.deriv2d(ro*psy,phi,z*1000.)
vini=-(1./cosphi)*dropsydz*1./ro
vinir=-(1./cosphi)*dropsydzr*1./ror


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
v_nextr = np.empty((nz,nx))
w_nextr = np.empty((nz,nx))
w_prevr = np.empty((nz,nx))
epsr = np.empty(nz)
v_nextr[:,:] = vinir[:,:]
w_prevr[:,:] = winir[:,:]
psyinir=np.empty((nz,nx))
psyinir[:,:]=psyr[:,:]


# Choix de la precision et du pourcent a comparer dans la boucle 
# uniquement necessaire si boucle while incluse
#percent = 1.
#precision = 100000.


###Calcul V et W tant que cela converge suivant la precision et le percent

#####DEBUT Boucle While #####
#while precision >= percent/100.:
for h in range(2,11):  #Pour controler l iteration (decommentez h=h+1 dans la boucle)

  ###3 equations qui calcul W avec 3 differentes methodes pour N2 (frequence de brunt vaisala)
  
  #equations qui convergent le mieux
  #Methode 1
  w_next=(Q[:,:]-(v_next[:,:]*dtempdphi/a))/(((H/Rs)*g/teta*dtetadz)+dtempdz)
  w_nextr=(Qr[:,:]-(v_nextr[:,:]*dtempdphir/a))/(((H/Rs)*g/tetar*dtetadzr)+dtempdzr)
  #Methode 2
  #w_next=(Q[:,:]-(v_next[:,:]*dtempdphi/a))/((H/Rs)*((-ro*g**2)/teta)*dtetadp+dtempdz)
  #w_nextr=(Qr[:,:]-(v_nextr[:,:]*dtempdphir/a))/((H/Rs)*((-ror*g**2)/tetar)*dtetadpr+dtempdzr)

  #1 equation equivalente qui diverge
  #w_next=(Q[:,:]-(v_next[:,:]*dtempdphi/a))/((H/Rs)*((g/T0)*((-g/cp)+dtempdz))+dtempdz)
  #w_nextr=(Qr[:,:]-(v_nextr[:,:]*dtempdphir/a))/((H/Rs)*((g/T0)*((-g/cp)+dtempdzr))+dtempdzr)

  #Equation qui calcul omega
  #w_next[mask]=(Q[mask]-(v_next[mask]*dtempdphi[mask]/a))/((temp[mask]/teta[mask])*dtetadp[mask])
  #w_nextr[mask]=(Qr[mask]-(v_nextr[mask]*dtempdphir[mask]/a))/((tempr[mask]/tetar[mask])*dtetadpr[mask])

  for i in range(nz):
    numer = integrate.simps(w_next[i,:]*cosphi[0:nx+1],phi[0:nx+1],axis=0)
    denom = integrate.simps(cosphi[0:nx+1],phi[0:nx+1],axis=0)
    numerr = integrate.simps(w_nextr[i,:]*cosphi[0:nx+1],phi[0:nx+1],axis=0)
    denomr = integrate.simps(cosphi[0:nx+1],phi[0:nx+1],axis=0)
    eps[i] = -numer/denom
    epsr[i] = -numerr/denomr
    for j in range(nx):   
      w_next[i,j]=w_next[i,j]+eps[i]
      w_nextr[i,j]=w_nextr[i,j]+epsr[i]
      integrand = (w_next[i,0:j+1])*cosphi[0:j+1]*a
      integrandr = (w_nextr[i,0:j+1])*cosphi[0:j+1]*a
      psy[i,j]=integrate.simps(integrand,phi[0:j+1],axis=0)
      psyr[i,j]=integrate.simps(integrandr,phi[0:j+1],axis=0)
  w = np.where(np.abs(w_next[:,:])-np.abs(w_prev[:,:]) != 0)
  wr = np.where(np.abs(w_nextr[:,:])-np.abs(w_prevr[:,:]) != 0)
  precision = np.mean(np.abs(w_next[w] - w_prev[w])/np.abs(0.5*(w_next[w]+w_prev[w])))
  precisionr = np.mean(np.abs(w_nextr[wr] - w_prevr[wr])/np.abs(0.5*(w_nextr[wr]+w_prevr[wr])))
  print precision, np.max(w_next[w]), np.min(w_next[w])
  print precisionr, np.max(w_nextr[wr]), np.min(w_nextr[wr])
  dropsydz,dropsydphi = np.gradient(ro*psy,z*1000.,phi, edge_order=2)
  dropsydzr,dropsydphir = np.gradient(ror*psyr,z*1000.,phi, edge_order=2)
  #dropsydphi,dropsydz = ppcompute.deriv2d(ro*psy,phi,z*1000.)
  v_next=-(1./cosphi)*dropsydz*1./ro
  v_nextr=-(1./cosphi)*dropsydzr*1./ror
  w_prev[:,:] = w_next[:,:]
  w_prevr[:,:] = w_nextr[:,:]
  #h += 1

#####Fin Boucle While #####

#print precision*100.,percent
print('On a realise ',h,'iterations')
#Tout les plots des resultats!

print np.mean((v_nextr-v_next)/(v_nextr))
print np.mean((w_nextr-w_next)/(w_nextr))

#####Tout les plots !#####
### IMPORTANT si Ls=integer -> mettre Ls=%i si Ls="string" -> Ls=%s

plt.figure()
plt.contourf(lat,zp,temp,levels=np.linspace(np.min(temp),np.max(temp),100),cmap=plt.cm.jet)
plt.colorbar()
CS=plt.contour(lat,zp,temp,levels=np.round(np.linspace(80,160,9)),colors='black')
plt.clabel(CS, inline=0.1, fontsize=10,fmt='% 3d')
plt.title('Temperature without rings (K) Ls=%s' %time_step) 
plt.yscale('log')
plt.ylabel('Pressure (Pa)')
plt.xlabel('Latitude')
plt.gca().invert_yaxis()

plt.figure()
plt.contourf(lat,zp,tempr,levels=np.linspace(np.min(tempr),np.max(tempr),100),cmap=plt.cm.jet)
plt.colorbar()
CS=plt.contour(lat,zp,tempr,levels=np.round(np.linspace(80,160,9)),colors='black')
plt.clabel(CS, inline=0.1, fontsize=10,fmt='% 3d')
plt.title('Temperature with rings (K) Ls=%s' %time_step) 
plt.yscale('log')
plt.ylabel('Pressure (Pa)')
plt.xlabel('Latitude')
plt.gca().invert_yaxis()

plt.figure()
plt.contourf(lat,zp,psy/cosphi,levels=np.arange(-np.max(abs(psy/cosphi)),np.max(abs(psy/cosphi)),step=np.max(abs(psy/cosphi))/100.),cmap=plt.cm.seismic)
plt.title('Final Streamfunction without rings Ls=%s' %time_step) 
plt.colorbar()
plt.yscale('log')
plt.gca().invert_yaxis()

plt.figure()
plt.contourf(lat,zp,psyr/cosphi,levels=np.arange(-np.max(abs(psyr/cosphi)),np.max(abs(psyr/cosphi)),step=np.max(abs(psyr/cosphi))/100.),cmap=plt.cm.seismic)
plt.title('Final Streamfunction with rings Ls=%s' %time_step) 
plt.colorbar()
plt.yscale('log')
plt.gca().invert_yaxis()

plt.figure()
plt.contourf(lat,zp,w_next,levels=np.arange(-np.max(abs(w_next)),np.max(abs(w_next)),step=np.max(abs(w_next))/50.),cmap=plt.cm.seismic)
plt.title('Final W without rings (m.s-1) Ls=%s' %time_step) 
plt.colorbar()
plt.yscale('log')
plt.ylabel('Pressure (Pa)')
plt.xlabel('Latitude')
plt.gca().invert_yaxis()

plt.figure()
plt.contourf(lat,zp,w_nextr,levels=np.arange(-np.max(abs(w_nextr)),np.max(abs(w_nextr)),step=np.max(abs(w_nextr))/50.),cmap=plt.cm.seismic)
plt.title('Final W with rings (m.s-1) Ls=%s' %time_step) 
plt.colorbar()
plt.yscale('log')
plt.ylabel('Pressure (Pa)')
plt.xlabel('Latitude')
plt.gca().invert_yaxis()

plt.figure()
plt.contourf(lat,zp,v_next,levels=np.arange(-np.max(abs(v_next)),np.max(abs(v_next)),step=np.max(abs(v_next))/100.),cmap=plt.cm.seismic)
plt.title('Final V without rings (m.s-1) Ls=%s' %time_step) 
plt.colorbar()
plt.yscale('log')
plt.ylabel('Pressure (Pa)')
plt.xlabel('Latitude')
plt.gca().invert_yaxis()

plt.figure()
plt.contourf(lat,zp,v_nextr,levels=np.arange(-np.max(abs(v_nextr)),np.max(abs(v_nextr)),step=np.max(abs(v_nextr))/100.),cmap=plt.cm.seismic)
plt.title('Final V with rings (m.s-1) Ls=%s' %time_step) 
plt.colorbar()
plt.yscale('log')
plt.ylabel('Pressure (Pa)')
plt.xlabel('Latitude')
plt.gca().invert_yaxis()

plt.figure(figsize=(12,5))
plt.subplot(1, 2, 1)
plt.contourf(lat,zp,w_nextr-w_next,levels=np.arange(-np.max(abs(w_nextr-w_next)),np.max(abs(w_nextr-w_next)),step=np.max(abs(w_nextr-w_next))/100.),cmap=plt.cm.seismic)
plt.title('W induit par les anneaux (m.s-1) Ls=%s' %time_step) 
plt.colorbar(format='%.0e')
plt.yscale('log')
plt.ylabel('Pressure (Pa)')
plt.xlabel('Latitude')
plt.gca().invert_yaxis()

#plt.figure()
plt.subplot(1, 2, 2)
plt.contourf(lat,zp,v_nextr-v_next,levels=np.arange(-np.max(abs(v_nextr-v_next)),np.max(abs(v_nextr-v_next)),step=np.max(abs(v_nextr-v_next))/100.),cmap=plt.cm.seismic)
plt.title('V induit par les anneaux (m.s-1) Ls=%s' %time_step) 
plt.colorbar(format='%.0e')
plt.yscale('log')
plt.ylabel('Pressure (Pa)')
plt.xlabel('Latitude')
plt.gca().invert_yaxis()

plt.figure(figsize=(7,9))
plt.subplot(2, 1, 1)
plt.contourf(lat,zp,w_nextr-w_next,levels=np.arange(-np.max(abs(w_nextr-w_next)),np.max(abs(w_nextr-w_next)),step=np.max(abs(w_nextr-w_next))/100.),cmap=plt.cm.seismic)
plt.title('W induit par les anneaux (m.s-1) Ls=%s' %time_step) 
plt.colorbar(format='%.0e')
plt.yscale('log')
plt.ylabel('Pressure (Pa)')
#plt.xlabel('Latitude')
plt.gca().invert_yaxis()

#plt.figure()
plt.subplot(2, 1, 2)
plt.contourf(lat,zp,v_nextr-v_next,levels=np.arange(-np.max(abs(v_nextr-v_next)),np.max(abs(v_nextr-v_next)),step=np.max(abs(v_nextr-v_next))/100.),cmap=plt.cm.seismic)
plt.title('V induit par les anneaux (m.s-1) Ls=%s' %time_step) 
plt.colorbar(format='%.0e')
plt.yscale('log')
plt.ylabel('Pressure (Pa)')
plt.xlabel('Latitude')
plt.gca().invert_yaxis()

plt.figure()
plt.contourf(lat,zp,psy*ro/cosphi,levels=np.arange(-np.max(abs(psy*ro/cosphi)),np.max(abs(psy*ro/cosphi)),step=np.max(abs(psy*ro/cosphi))/100.),cmap=plt.cm.seismic)
plt.title('Final Mass-Streamfunction without rings Ls=%s' %time_step) 
plt.colorbar()
plt.yscale('log')
plt.gca().invert_yaxis()

plt.figure()
plt.contourf(lat,zp,psyr*ror/cosphi,levels=np.arange(-np.max(abs(psyr*ror/cosphi)),np.max(abs(psyr*ror/cosphi)),step=np.max(abs(psyr*ror/cosphi))/100.),cmap=plt.cm.seismic)
plt.title('Final Mass-Streamfunction with rings Ls=%s' %time_step) 
plt.colorbar()
plt.yscale('log')
plt.gca().invert_yaxis()

plt.figure()
spe=np.sqrt(abs(v_nextr[:,1:-2])**2+abs(w_nextr[:,1:-2]**2))
plt.streamplot(lat[1:-2],z,v_nextr[:,1:-2],w_nextr[:,1:-2]*1000.,color=spe,cmap='coolwarm')
plt.title('Streamlines with factor vert/horiz=1000 with rings Ls=%s' %time_step) 
#plt.yscale('log')
plt.colorbar()
plt.ylabel('Altitude (km)')
plt.xlabel('Latitude (nearly 180000 km)')
#plt.gca().invert_yaxis()

plt.figure()
spe=np.sqrt(abs(v_next[:,1:-2])**2+abs(w_next[:,1:-2]**2))
plt.streamplot(lat[1:-2],z,v_next[:,1:-2],w_next[:,1:-2]*1000.,color=spe,cmap='coolwarm')
plt.title('Streamlines with factor vert/horiz=1000 without rings Ls=%s' %time_step) 
#plt.yscale('log')
plt.colorbar()
plt.ylabel('Altitude (km)')
plt.xlabel('Latitude (nearly 180000 km)')
#plt.gca().invert_yaxis()

plt.figure()
spe=np.sqrt(abs(v_nextr[:,1:-2]-v_next[:,1:-2])**2+abs((w_nextr[:,1:-2]-w_next[:,1:-2])**2))
plt.streamplot(lat[1:-2],z,v_nextr[:,1:-2]-v_next[:,1:-2],w_nextr[:,1:-2]*1000.-w_next[:,1:-2]*1000.)#,color=spe,cmap='autumn')
plt.title('Streamlines induite par les anneaux Ls=%s' %time_step) 
#plt.yscale('log')
#plt.colorbar()
plt.ylabel('Altitude (km)')
plt.xlabel('Latitude (nearly 180000 km)')
#plt.gca().invert_yaxis()

plt.figure()
spe=np.sqrt(abs(v_nextr[:,1:-2]-v_next[:,1:-2])**2+abs((w_nextr[:,1:-2]-w_next[:,1:-2])*1000.**2))
plt.streamplot(lat[1:-2],z,v_nextr[:,1:-2]-v_next[:,1:-2],w_nextr[:,1:-2]*1000.-w_next[:,1:-2]*1000.,color=spe,cmap='autumn')
plt.title('Streamlines with minus without rings Ls=%s' %time_step) 
#plt.yscale('log')
plt.colorbar()
plt.ylabel('Altitude (km)')
plt.xlabel('Latitude (nearly 180000 km)')
#plt.gca().invert_yaxis()

####COMMENTEZ ci dessous pour aller plus loin

plt.show(block=False)
input("Hit it and quit it! (press enter to close all)")
plt.close()
exit()

######ATTENTION les valeurs sont bonnes uniquement si le tableau d'altitude va bien de 16 a 56

plt.figure()
plt.plot(lat[21:106],w_next[19,21:106],'--b')#19 1mbar
plt.plot(lat[21:106],w_next[29,21:106],'--r')#29 [21:106] 0.1mbar
plt.plot(lat[21:106],w_nextr[19,21:106],'b')#19 1mbar
plt.plot(lat[21:106],w_nextr[29,21:106],'r')#29 [21:106] 0.1mbar
#plt.plot(lat[21:106],w_next[38,21:106],'g')#38 0.01mbar
plt.yticks(np.arange(-0.001,0.001,step=0.0002)) #pour valeur fixe
#plt.yticks(np.arange(-0.0002,0.0002,step=0.00002)) #pour moyenne annuelle
plt.xticks(np.arange(-50.,55.,step=10.))
plt.legend(('1mbar no rings','0.1mbar no rings','1mbar rings','0.1mbar rings'),loc='upper right')
plt.ylabel('W(m.s-1)')
plt.xlabel('Latitude')
plt.title('W final Ls=%s' %time_step) 
plt.grid(True)

plt.figure()
#plt.plot(lat[21:106],np.mean(v_next[3:13,21:106],0),'r--')#3:13 mean50-5mbar
plt.plot(lat[21:106],v_next[17,21:106],'g--')#17 2mbar
#plt.plot(lat[21:106],np.mean(v_nextr[3:13,21:106],0),'r')#3:13 mean50-5mbar
plt.plot(lat[21:106],v_nextr[17,21:106],'g')#17 2mbar
plt.plot(lat[21:106],v_next[19,21:106],'--b')#19 1mbar
plt.plot(lat[21:106],v_nextr[19,21:106],'b')#19 1mbar
plt.yticks(np.arange(-0.06,0.06,step=0.01)) #pour valeur fixe
#plt.yticks(np.arange(-0.01,0.01,step=0.001)) #pour moyenne annuelle
plt.xticks(np.arange(-50.,55.,step=10.))
plt.legend(('2mbar no rings','2mbar rings','1mbar no rings','1mbar rings'),loc='upper right')
plt.ylabel('W(m.s-1)')
plt.xlabel('Latitude')
plt.title('V final Ls=%s' %time_step) 
plt.grid(True)

plt.show(block=False)
input("Hit it and quit it! (press enter to close all)")
plt.close()
