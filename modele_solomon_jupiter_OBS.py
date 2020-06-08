import pylab as py
import numpy as np
#import ppcompute
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy.matlib
from ppclass import pp
from scipy import integrate

###Emplacement du fichier, choix de la fenetre de temps, recuperation des variables puis choix de la grille (altitude, lat ..)

# -> Choix du fichier contenant les températures observées ainsi que 
#    les taux de chauffage et redroidissement calculés à partir de températures observées également
#fff = "/planeto/sglmd/Jupiter/diagfi_aurora_Fint7-485_cloud_900mbar_suite2_avec_teta.nc"
fff = "/planeto/sglmd/Jupiter/diagfi_3D_construit_FLETCHER_taux_chauffage_aero_Zhang_BAS.nc"

# -> Choix du pas de temps
time_step=1
#time_step="360,720"
#"3240,3600"=8e annee en entier (moyenne annuelle)
#"3360,3380"=134 EteHN 
#"3540,3560"=308 hiverHN
#"3570,3600"=histins195 (330-360) correspond a DYNAMICO-195
#"3560,3590"=histins195 (320-350) correspond a DYNAMICO-170

infz=13     #borne inf de Z ###0  niveau 16 = 140 mbar = 76km
topz=48     #borne sup de Z ###64 niveau 56 = 0.01mbar = 316km

nn=3        #Indice pour retirer les poles sur les latitudes

# -> Recuperation des variables dans le fichier moyennes zonalement avec ls corrige

teta,lon,lat,zmod,t = pp(file=fff,var="teta",t=time_step,x="-180,180",verbose=True,changetime="correctls").getfd()
zdtsw= pp(file=fff,var="zdtsw",t=time_step,x="-180,180",verbose=True,changetime="correctls").getf()
zdtlw= pp(file=fff,var="zdtlw",t=time_step,x="-180,180",verbose=True,changetime="correctls").getf()
Q = zdtsw + zdtlw
Q = Q[infz:topz,nn:-nn]
pmod= pp(file=fff,var="p",t=time_step,x="-180,180",verbose=True,changetime="correctls").get()
tempmod= pp(file=fff,var="temp",t=time_step,x="-180,180",verbose=True,changetime="correctls").get()
psurf= pp(file=fff,var="ps",t=time_step,x="-180,180",verbose=True,changetime="correctls").get()

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

lat=lat[nn:-nn]     #On retire les poles (sinon cosphi -> 0 et division par 0)
phi = lat*np.pi/180.
sinphi = np.sin(phi)
cosphi = np.cos(phi)   #Transformation des latitudes en phi

#Decoupe de l'altitude et creation du champ de taux radiatif
zp=zmod[infz:topz]
#Qrad=zdtsw+zdtlw
#Q=Qrad.f[infz:topz,nn:-nn]


###Creation d'un champ de pression et temperature a partir de Z et teta
p2D=pmod.f[infz:topz,nn:-nn]
p=p2D[:,0]
temp=tempmod.f[infz:topz,nn:-nn]
teta=teta[infz:topz,nn:-nn]
print p


#########  IMPORTANT  #########
###Definitions de toute les constantes planetaires

#Saturn
#a=58210000.            #Rayon de la planete en m
#g=10.                  #g acceleration
#T0=135.                #Temperature au "sol"

#Jupiter:
a=69911000.            #Rayon de la planete en m
g=24.79                #g acceleration
T0=170.                #Temperature au "sol" = 1 bar?
cp=11500.              #Capacite calorifique
Rm=8.3144621           #Constante des gaz parfaits
M=2.34                 #Masse molaire gaz
H=(Rm*T0)/(g*M)*1000.  #Calcul hauteur caracteristique en m
Rs=(Rm/M)*1000.        #Constante specifique des gaz parfaits
ro=p2D/(Rs*temp)       #Calcul du champ de masse volumique (density)
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
z2D=np.zeros(np.shape(pmod.f))
z=np.zeros(np.shape(pmod.f)[0])
z2D[0,:]=0     #-H*np.log(pmod.f[0,:]/p0)/1000.
for i in range(1,np.shape(z2D)[0]):
  for j in range(1,np.shape(z2D)[1]):
    z2D[i,j]=z2D[i-1,j]+((Rs/g)*((tempmod.f[i,j]+tempmod.f[i-1,j])/2.)*np.log(pmod.f[i-1,j]/pmod.f[i,j]))/1000.
z=np.mean(z2D,1)
z=z[infz:topz]



###Valeurs de P depuis l'equation hypsometrique en moyennant T a chaque point du profil vertical si altitude disponible mais pas P dans le fichier nc
#p2D=np.empty(np.shape(z2D))
#p=np.empty(np.shape(z2D)[0])
#p2D[0,:]=psurf.f
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



###Initialisation de w(omega) premiere iteration avec v=0
wini=np.empty((nz,nx))
winisave=np.empty((nz,nx))


###3 equations qui calcul W avec 3 differentes methodes pour N2 (frequence de brunt vaisala)


#equations qui convergent le mieux
#Methode 1 <- CONSEILLEE et utilisee pour le rapport
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
#Tout les plots des resultats!
print np.mean(abs(w_next))
print np.mean(abs(v_next))
print np.mean(abs(Q))


#f = open("wstar_BAS.txt",'w') 
#for fff in w_next[12,:]:
#    f.write("%+8.5e\n" % fff)
#f.close
#f = open("lat.txt",'w')
#for lll in lat:
#    f.write("%+8.5e\n" % lll)
#f.close
#f = open("p.txt",'w')
#for ppp in p:
#    f.write("%+8.5e\n" % ppp)
#f.close


#####Tout les plots !#####
### IMPORTANT si Ls=integer -> mettre Ls=%i si Ls="string" -> Ls=%s
#colormap ,cmap=plt.cm.seismic

plt.figure(figsize=(12,5))
plt.subplot(1, 2, 1)
#plt.contourf(lat,zp,Q,levels=np.linspace(np.min(Q),np.max(Q),100),cmap=plt.cm.jet)
plt.contourf(lat,zp/100,Q,levels=np.arange(-3.e-6,3.e-6,step=1e-7),cmap=plt.cm.jet)
plt.title('Chauffage radiatif (K.s-1) Ls=%s' %time_step) 
plt.colorbar(format='%.1e')
plt.yscale('log')
plt.ylabel('Pressure (mbar)')
plt.xlabel('Latitude')
plt.gca().invert_yaxis()

#plt.figure()
plt.subplot(1, 2, 2)
#plt.contourf(lat,zp,temp,levels=np.linspace(np.min(temp),np.max(temp),100),cmap=plt.cm.jet)
plt.contourf(lat,zp/100,temp,levels=np.arange(100,200,step=1),cmap=plt.cm.jet)
plt.colorbar(format='%.0f')
CS=plt.contour(lat,zp/100,temp,levels=np.round(np.linspace(105,175,15)),colors='black')
plt.clabel(CS, inline=0.1, fontsize=10,fmt='% 3d')
plt.title('Temperature (K) Ls=%s' %time_step) 
plt.yscale('log')
plt.ylabel('Pressure (mbar)')
plt.xlabel('Latitude')
plt.gca().invert_yaxis()

plt.figure()
#plt.contourf(lat,zp,Q,levels=np.linspace(np.min(Q),np.max(Q),100),cmap=plt.cm.jet)
plt.contourf(lat,zp/100,Q,levels=np.arange(-3.5e-6,1.e-6,step=1e-7),cmap=plt.cm.jet)

#plt.title('Net heating rates (K.s-1) Ls=%s' %time_step) 
plt.colorbar(format='%.1e')
CS=plt.contour(lat,zp/100,Q,levels=np.arange(-1e-6,1.e-6,step=2e-7),colors='black')
plt.yscale('log')
plt.ylabel('Pressure (mbar)')
plt.xlabel('Latitude')
plt.gca().invert_yaxis()
plt.savefig('Net_heating_rates_OBS_aero_skewed_Zhang_BAS.eps', format='eps')



#plt.figure()
#plt.contourf(lat,zp,teta,levels=np.linspace(np.min(teta),np.max(teta),100),cmap=plt.cm.jet)
#plt.colorbar()
#plt.title('Teta (K) Ls=340-355' ) 
#plt.ylabel('Pressure (Pa)')
#plt.xlabel('Latitude')
#plt.yscale('log')
#plt.gca().invert_yaxis()

plt.figure()
plt.contourf(lat,zp,psy/cosphi,levels=np.arange(-np.max(abs(psy/cosphi)),np.max(abs(psy/cosphi)),step=np.max(abs(psy/cosphi))/100.),cmap=plt.cm.seismic)
plt.title('Final Streamfunction Ls=%s' %time_step) 
plt.colorbar()
plt.yscale('log')
plt.gca().invert_yaxis()

plt.figure(figsize=(14,5))
plt.subplot(1, 2, 1)
#plt.contourf(lat,zp/100,w_next,levels=np.arange(-np.max(abs(w_next)),np.max(abs(w_next)),step=np.max(abs(w_next))/100.),cmap=plt.cm.seismic)
plt.contourf(lat,zp/100,w_next,levels=np.arange(-0.5e-3,0.5e-3,step=0.00001),cmap=plt.cm.seismic)
plt.title('w* (m.s-1)') 
plt.colorbar(format='%.0e')
plt.yscale('log')
plt.ylabel('Pressure (mbar)')
plt.xlabel('Latitude')
plt.gca().invert_yaxis()

#plt.figure()
plt.subplot(1, 2, 2)
#plt.contourf(lat,zp,v_next,levels=np.arange(-np.max(abs(v_next)),np.max(abs(v_next)),step=np.max(abs(v_next))/50.),cmap=plt.cm.seismic)
plt.contourf(lat,zp/100,v_next,levels=np.arange(-0.5,0.5,step=0.01),cmap=plt.cm.seismic)
plt.title('v* (m.s-1)') 
plt.colorbar(format='%.0e')
plt.yscale('log')
plt.ylabel('Pressure (mbar)')
plt.xlabel('Latitude')
plt.gca().invert_yaxis()
plt.savefig('Vitesses_OBS_aero_skewed_Zhang_BAS.eps', format='eps')

#plt.figure(figsize=(7,9))
#plt.subplot(2, 1, 1)
#plt.contourf(lat,zp,w_next,levels=np.arange(-np.max(abs(w_next)),np.max(abs(w_next)),step=np.max(abs(w_next))/100.),cmap=plt.cm.seismic)
#plt.title('Final W (m.s-1) Ls=%s' %time_step) 
#plt.colorbar(format='%.0e')
#plt.yscale('log')
#plt.ylabel('Pressure (Pa)')
#plt.xlabel('Latitude')
#plt.gca().invert_yaxis()

#plt.figure()
#plt.subplot(2, 1, 2)
#plt.contourf(lat,zp,v_next,levels=np.arange(-np.max(abs(v_next)),np.max(abs(v_next)),step=np.max(abs(v_next))/100.),cmap=plt.cm.seismic)
#plt.title('Final V (m.s-1) Ls=%s' %time_step) 
#plt.colorbar(format='%.0e')
#plt.yscale('log')
#plt.ylabel('Pressure (Pa)')
#plt.xlabel('Latitude')
#plt.gca().invert_yaxis()

plt.figure()

##plt.contourf(lat,zp,psy*ro/cosphi,levels=np.arange(-np.max(abs(psy*ro/cosphi)),np.max(abs(psy*ro/cosphi)),step=np.max(abs(psy*ro/cosphi))/100.),cmap=plt.cm.seismic)
#plt.contourf(lat,zp/100,np.log10(psy*ro/cosphi)) #,levels=np.arange(-1.e2,1.e2,step=1) ,cmap=plt.cm.gist_stern) #seismic)
##plt.contourf(lat,zp/100,psy*ro/cosphi,locator=ticker.MaxNLocator(100), aspect='auto',origin='lower' ,cmap=plt.cm.seismic)
#CS=plt.contour(lat,zp/100,psy*ro/cosphi,levels=np.arange(-1e2,1.e2,step=10),colors='black')
##CS=plt.contour(lat,zp/100,psy*ro/cosphi,levels=np.round(np.logspace(0.1,100,10)),colors='black')






lev = [\
-1e3,\
-5e2,\
-1e2,\
-5e1,\
-1e1,\
-5e0,\
-1e0,\
-5e-1,\
-1e-1,\
+1e-1,\
+5e-1,\
+1e0,\
+5e0,\
+1e1,\
+5e1,\
+1e2,\
+5e2,\
+1e3]

lev = [\
-5e1,\
-1e1,\
-5e0,\
-1e0,\
-5e-1,\
+5e-1,\
+1e0,\
+5e0,\
+1e1,\
+5e1]

CS=plt.contour(lat,zp/100,psy*ro/cosphi,levels=lev,colors='black')
plt.clabel(CS, CS.levels, inline=False, fmt="%.1f", fontsize=10)

fifi = psy*ro/cosphi
fifi[fifi <= lev[0]] = lev[0]
fifi[fifi >= lev[-1]] = lev[-1]

plt.contourf(lat,zp/100,fifi,levels=lev,cmap=plt.cm.seismic)

plt.title('Mass streamfunction') 
#plt.colorbar()
plt.yscale('log')
plt.ylabel('Pressure (mbar)')
plt.xlabel('Latitude')
plt.gca().invert_yaxis()

#plt.figure()
#plt.quiver(lat[1:-2],z,v_next[:,1:-2],w_next[:,1:-2]*1000.)#,cmap=plt.cm.seismic)
#plt.title('Vectors Ls=%s' %time_step) 
##plt.yscale('log')
#plt.ylabel('Altitude (km)')
#plt.xlabel('Latitude')

plt.figure()
spe=np.sqrt(abs(v_next[:,1:-2])**2+abs(w_next[:,1:-2])**2)
plt.streamplot(lat[1:-2],z,v_next[:,1:-2],w_next[:,1:-2]*900.)#,color=spe,cmap='autumn')#,cmap=plt.cm.seismic)
#plt.title('Streamlines' %time_step) 
#plt.yscale('log')
#plt.colorbar()
plt.ylabel('Altitude (km)')
plt.xlabel('Latitude')
#plt.gca().invert_yaxis()
plt.savefig('Streamlines_OBS_aero_skewed_Zhang_BAS.eps', format='eps')

plt.show(block=False)
input("Hit it and quit it! (press enter to close all)")
plt.close()
exit()

#plt.figure()
#spe=np.sqrt(abs(v_next[:,1:-2])**2+abs(w_next[:,1:-2])**2)
#plt.streamplot(lat[1:-2],zp/100,v_next[:,1:-2],w_next[:,1:-2]*900.)#,color=spe,cmap='autumn')#,cmap=plt.cm.seismic)
#plt.yscale('log')
#plt.ylabel('Pressure (mbar)')
#plt.xlabel('Latitude')
#plt.gca().invert_yaxis()
#plt.savefig('Streamlines_p.eps', format='eps')


#plt.figure()
#spe=np.sqrt(abs(v_next[:,1:-2])**2+abs(w_next[:,1:-2])*1000.**2)
#plt.streamplot(lat[1:-2],z,v_next[:,1:-2],w_next[:,1:-2]*1000.,color=spe,cmap='autumn')#,cmap=plt.cm.seismic)
#plt.title('Streamlines with factor vert/horiz=1000 Ls=%s' %time_step) 
#plt.yscale('log')
#plt.colorbar()
#plt.ylabel('Altitude (km)')
#plt.xlabel('Latitude (nearly 180000 km)')
#plt.gca().invert_yaxis()

v_nextnan=np.empty(np.shape(v_next))
v_nextnan[:,:]=v_next[:,:]
v_nextnan[:,49:77]=np.NaN

w_nextnan=np.empty(np.shape(w_next))
w_nextnan[:,:]=w_next[:,:]
w_nextnan[:,49:77]=np.NaN


#[:,42:84] 30/-25
#[:,49:77] 20/-20


#plt.figure(figsize=(14,8))
#plt.subplot(2, 1, 1)
#plt.contourf(lat,zp,v_nextnan,levels=np.arange(-np.nanmax(abs(v_nextnan)),np.nanmax(abs(v_nextnan)),step=np.nanmax(abs(v_nextnan))/100.),cmap=plt.cm.seismic)
#plt.title('Final V (m.s-1) Ls=%s' %time_step) 
#plt.colorbar(format='%.0e')
#plt.yscale('log')
#plt.ylabel('Pressure (Pa)')
#plt.xlabel('Latitude')
#plt.gca().invert_yaxis()

#plt.figure(figsize=(10,5))
#plt.subplot(2, 1, 2)
#plt.contourf(lat,zp,w_nextnan,levels=np.arange(-np.nanmax(abs(w_nextnan)),np.nanmax(abs(w_nextnan)),step=np.nanmax(abs(w_nextnan))/100.),cmap=plt.cm.seismic)
#plt.title('Final W (m.s-1) Ls=%s' %time_step) 
#plt.colorbar(format='%.0e')
#plt.yscale('log')
#plt.ylabel('Pressure (Pa)')
#plt.xlabel('Latitude')
#plt.gca().invert_yaxis()

####COMMENTEZ ci dessous pour aller plus loin
plt.show(block=False)
input("Hit it and quit it! (press enter to close all)")
plt.close()
exit()

######ATTENTION les valeurs sont bonnes uniquement si le tableau d'altitude va bien de 16 a 58

plt.figure()
plt.plot(lat[21:106],w_next[19,21:106],'b--')#19 1mbar
plt.plot(lat[21:106],w_next[29,21:106],'r--')#29 [21:106] 0.1mbar
#plt.plot(lat[21:106],w_next[38,21:106],'g')#38 0.01mbar
plt.yticks(np.arange(-0.001,0.001,step=0.0002)) #pour valeur fixe
#plt.yticks(np.arange(-0.0002,0.0002,step=0.00002)) #pour moyenne annuelle
plt.xticks(np.arange(-50.,55.,step=10.))
plt.legend(('1mbar','0.1mbar'),loc='upper right')
plt.ylabel('W(m.s-1)')
plt.xlabel('Latitude')
plt.title('W final Ls=%s' %time_step) 
plt.grid(True)

plt.figure(figsize=(8,5))
#plt.plot(lat[21:106],np.mean(v_next[3:13,21:106],0),'r--')#3:13 mean50-5mbar
#plt.plot(lat[27:99],v_next[17,27:99],'g--')#17 2mbar
#plt.plot(lat[27:99],v_next[19,27:99],'b--')#19 1mbar
plt.plot(lat,v_next[17,:],'g--')#17 2mbar
plt.plot(lat,v_next[19,:],'b--')#19 1mbar
#plt.yticks(np.arange(-0.07,0.07,step=0.01)) #pour valeur fixe
plt.yticks(np.arange(-0.006,0.006,step=0.001)) #pour moyenne annuelle
#plt.xticks(np.arange(-50.,55.,step=10.))
plt.xticks(np.arange(-90.,95.,step=20.))
plt.legend(('2mbar','1mbar'),loc='upper right')
plt.ylabel('V (m.s-1)')
plt.xlabel('Latitude')
plt.title('V final (m.s-1) Ls=%s' %time_step) 
plt.grid(True)

plt.figure()
plt.plot(lat,w_next[39,:],'g--')#17 1Pa
#plt.plot(lat,w_next[48,:],'b--')#19 0.1Pa
#plt.plot(lat[21:106],w_next[19,21:106],'b--')#19 1mbar
#plt.plot(lat[21:106],w_next[29,21:106],'r--')#29 [21:106] 0.1mbar
#plt.plot(lat[21:106],w_next[38,21:106],'g')#38 0.01mbar
#plt.yticks(np.arange(-0.001,0.001,step=0.0002)) #pour valeur fixe
#plt.yticks(np.arange(-0.0002,0.0002,step=0.00002)) #pour moyenne annuelle
#plt.xticks(np.arange(-50.,55.,step=10.))
plt.legend(('1Pa','0.1mbar'),loc='upper right')
plt.ylabel('W(m.s-1)')
plt.xlabel('Latitude')
plt.title('W final Ls=%s' %time_step) 
plt.grid(True)

plt.figure()
#plt.plot(lat[21:106],np.mean(v_next[3:13,21:106],0),'r--')#3:13 mean50-5mbar
#plt.plot(lat[27:99],v_next[17,27:99],'g--')#17 2mbar
#plt.plot(lat[27:99],v_next[19,27:99],'b--')#19 1mbar
plt.plot(lat,v_next[39,:],'g--')#17 1Pa
#plt.plot(lat,v_next[48,:],'b--')#19 0.1Pa
#plt.yticks(np.arange(-0.07,0.07,step=0.01)) #pour valeur fixe
#plt.yticks(np.arange(-0.008,0.008,step=0.001)) #pour moyenne annuelle
#plt.xticks(np.arange(-50.,55.,step=10.))
#plt.xticks(np.arange(-90.,95.,step=20.))
plt.legend(('1Pa','1mbar'),loc='upper right')
plt.ylabel('W(m.s-1)')
plt.xlabel('Latitude')
plt.title('V final Ls=%s' %time_step) 
plt.grid(True)

print zp[39]

plt.show(block=False)
input("Hit it and quit it! (press enter to close all)")
plt.close()
