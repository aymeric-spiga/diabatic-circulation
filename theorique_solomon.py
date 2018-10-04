import pylab as py
import numpy as np
import ppcompute
import matplotlib.pyplot as plt
import numpy.matlib
from scipy import integrate

###Initialisation du domaine, des profils vertical et latitudinal + caracteristiques planete

nz, nx = (90, 180)  #Domaine 9 valeurs de Z et 18 de latitudes 
z=np.linspace(15.,85.,nz)  #Profil vertical en km
#z=np.array([16.,21.,28.,34.,40.,50.,60.,70.,80.]) #Profil vertical theorique
lat=np.empty((1,nx))
lat=np.linspace(-85., 85., nx)   #Profil latitudinal
phi = lat*np.pi/180.
sinphi = np.sin(phi)
cosphi = np.cos(phi)   #Transformation des latitudes en phi


###Creation d un champ radiatif 2D, refroidit a gauche et rechauffe a droite: On cree un linspace au sommet de l'atmosphere et on le fait tendre vers 0 vers le bas
Q=np.empty((nz,nx))
Q[nz-1,:] = np.linspace(0.00015, -0.00030, nx) #Valeurs pour +-0.5K/heure (+-12K/jour)
for i in range(nx):
  Q[:,i] = py.flip(np.linspace(Q[nz-1,i], 0, nz),0)


####Chauffage code pour cas plus complexe
#yorgl = np.linspace(-0.00015, 0.00015, nx)*np.cos(np.pi*lat/150.)
#for i in range(nz):
#  Q[i,:] = yorgl*np.cos(np.pi*(z[i]+50.)/100.)


###Creation d un champ de temperature potentielle teta qui varie sur la latitude: On cree la premiere et derniere colonne et on lisse avec un linspace

teta=(nz,nx)
teta=np.empty(teta)
teta[:,0]=np.linspace(400.,410.,nz) #400.,1400. 450.,1400.
teta[:,nx-1]=np.linspace(450.,470.,nz)  #500.,1200. 400.,1200.
for i in range(nz):
  teta[i,:] = np.linspace(teta[i,0],teta[i,nx-1],nx) #Creation de la forme teta en 2D


###Definitions de toute les constantes planetaires + Creation d'un champ de pression et temperature a partir de Z et teta

a=6371000.             #Rayon de la planete en m
g=9.81                 #g acceleration
T0=288.                #Temperature au "sol"
cp=1004.               #Capacite calorifique
Rm=8.3144621           #Constante des gaz parfaits
M=29.0                 #Masse molaire gaz
H=(Rm*T0)/(g*M)*1000.  #Calcul hauteur caracteristique en m
Rs=(Rm/M)*1000.        #Constante specifique des gaz parfaits
p0=100000.             #Pression a la base de l'atmosphere

p=p0*np.exp((-z*1000.)/(H))  #Pression theorique sur terre en Pascal
p2D=np.rot90(np.matlib.repmat(p,nx,1),axes=(1,0))   #creation champ 2D
temp=(nz,nx)
temp=np.empty(temp)
temp=teta*(p2D/p0)**(Rs/cp)    #creation temperature theorique
ro=p2D/(Rs*temp)       #Calcul du champ de masse volumique (density)


###Calcul des derivees de teta et temp en phi,z et p

dtetadz,dtetadphi=np.gradient(teta,z*1000.,phi, edge_order=2)
dtetadp,dtetadphi=np.gradient(teta,p,phi, edge_order=2)
dtempdz,dtempdphi=np.gradient(temp,z*1000.,phi, edge_order=2)


###Initialisation de w,vitesse verticale, premiere iteration avec v=0

wini=np.empty((nz,nx))
winisave=np.empty((nz,nx))

###3 equations qui calcul W avec 3 differentes methodes pour N2 (frequence de brunt vaisala)

#equations qui convergent le mieux
#wini=Q[:,:]/((H/Rs)*((-ro*g**2)/teta)*dtetadp+dtempdz)
wini=Q[:,:]/(((H/Rs)*g/teta*dtetadz)+dtempdz)

#1 equation equivalente qui diverge
#wini=Q[:,:]/((H/Rs)*((g/T0)*((-g/cp)+dtempdz))+dtempdz)

#Equation qui calcul omega
#wini=Q[:,:]/((temp/teta)*dtetadp)


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
for h in range(2,11):  #Pour controler l iteration (decommentez h=h+1 dans la boucle)

  ###3 equations qui calcul W avec 3 differentes methodes pour N2 (frequence de brunt vaisala)

  #equations qui convergent le mieux
  #w_next=(Q[:,:]-(v_next[:,:]*dtempdphi/a))/((H/Rs)*((-ro*g**2)/teta)*dtetadp+dtempdz)
  w_next=(Q[:,:]-(v_next[:,:]*dtempdphi/a))/(((H/Rs)*g/teta*dtetadz)+dtempdz)

  #1 equation equivalente qui diverge
  #w_next=(Q[:,:]-(v_next[:,:]*dtempdphi/a))/((H/Rs)*((g/T0)*((-g/cp)+dtempdz))+dtempdz)

  #Equation qui calcul omega
  #w_next=(Q[:,:]-(v_next[:,:]*dtempdphi/a))/((temp/teta)*dtetadp)

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
#Tout les plots des resultats!

#####Tout les plots !#####

plt.figure()
plt.contourf(lat,z,winisave,levels=np.linspace(np.min(winisave),np.max(winisave),100))
plt.title('W initial avant correction')
plt.colorbar()

plt.figure(figsize=(12,5))
plt.subplot(1, 2, 1)
plt.contourf(lat,z,Q,levels=np.linspace(np.min(Q),np.max(Q),100),cmap=plt.cm.jet)
plt.title('Taux de chauffage radiatif (K.s-1)')
plt.ylabel('Altitude (km)')
plt.xlabel('Latitude')
plt.colorbar(format='%.0e')

#plt.figure()
plt.subplot(1, 2, 2)
plt.contourf(lat,z,temp,levels=np.linspace(np.min(temp),np.max(temp),100),cmap=plt.cm.jet)
plt.title('Temperature usuelle (K)')
plt.ylabel('Altitude (km)')
plt.xlabel('Latitude')
plt.colorbar(format='%.0f')

plt.figure()
plt.contourf(lat,z,teta,levels=np.linspace(np.min(teta),np.max(teta),100),cmap=plt.cm.jet)
plt.title('Temperature potentielle')
plt.ylabel('Altitude (km)')
plt.xlabel('Latitude')
plt.colorbar()

plt.figure()
plt.contourf(lat,z,p2D,levels=np.linspace(np.min(p2D),np.max(p2D),100),cmap=plt.cm.jet)
plt.title('p')
plt.ylabel('Altitude (km)')
plt.xlabel('Latitude')
plt.colorbar()

plt.figure()
plt.contourf(lat,z,wini,levels=np.linspace(np.min(wini),np.max(wini),100))
plt.title('W initial')
plt.ylabel('Altitude (km)')
plt.xlabel('Latitude')
plt.colorbar()

plt.figure()
plt.contourf(lat,z,psyini,levels=np.linspace(np.min(psyini),np.max(psyini),100))
plt.title('Psy Initial')
plt.ylabel('Altitude (km)')
plt.xlabel('Latitude')
plt.colorbar()

plt.figure()
plt.contourf(lat,z,vini,levels=np.linspace(np.min(vini),np.max(vini),100))
plt.title('V initial')
plt.ylabel('Altitude (km)')
plt.xlabel('Latitude')
plt.colorbar()

plt.figure()
plt.contourf(lat,z,psy,levels=np.linspace(np.min(psy),np.max(psy),100))
plt.title('Psy final')
plt.ylabel('Altitude (km)')
plt.xlabel('Latitude')
plt.colorbar()

plt.figure(figsize=(12,5))
plt.subplot(1, 2, 1)
plt.contourf(lat,z,w_next,levels=np.arange(-np.max(abs(w_next)),np.max(abs(w_next)),step=np.max(abs(w_next))/100.),cmap=plt.cm.seismic)
plt.title('W final (m.s-1)')
plt.ylabel('Altitude (km)')
plt.xlabel('Latitude')
plt.colorbar(format='%.2f')

#plt.figure()
plt.subplot(1, 2, 2)
plt.contourf(lat,z,v_next,levels=np.arange(-np.max(abs(v_next)),np.max(abs(v_next+1)),step=np.max(abs(v_next))/100.),cmap=plt.cm.seismic)
plt.title('V final (m.s-1)')
plt.ylabel('Altitude (km)')
plt.xlabel('Latitude')
plt.colorbar(format='%.2f')
#plt.yscale('log')
#plt.savefig('V_final')

plt.figure()
plt.contourf(lat,z,psy*ro/cosphi,levels=np.linspace(np.min(psy*ro/cosphi),np.max(psy*ro/cosphi),100))
plt.title('Psymass final')
plt.ylabel('Altitude (km)')
plt.xlabel('Latitude')
plt.colorbar()

plt.figure()
spe=np.sqrt(abs(v_next[:,1:-2])**2+abs(w_next[:,1:-2])**2)
plt.streamplot(lat[1:-2],z,v_next[:,1:-2],w_next[:,1:-2]*250.)#,density=spe)#,cmap='autumn')#,cmap=plt.cm.seismic)
plt.title('Streamlines with factor vert/horiz=250') 
#plt.yscale('log')
plt.ylabel('Altitude (km)')
plt.xlabel('Latitude (nearly 180000 km)')
#plt.gca().invert_yaxis()


plt.figure()
spe=np.sqrt(abs(v_next[:,1:-2])**2+abs(w_next[:,1:-2])*1000.**2)
plt.streamplot(lat[1:-2],z,v_next[:,1:-2],w_next[:,1:-2]*250.)#,density=spe)#,cmap='autumn')#,cmap=plt.cm.seismic)
plt.title('Streamlines with factor vert/horiz=250') 
#plt.yscale('log')
plt.ylabel('Altitude (km)')
plt.xlabel('Latitude (nearly 180000 km)')
#plt.gca().invert_yaxis()


plt.show(block=False)
input("Hit it and quit it! (press enter to close all)")
plt.close()
