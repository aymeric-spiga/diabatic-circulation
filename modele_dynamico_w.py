import pylab as py
import numpy as np
#import ppcompute
import matplotlib.pyplot as plt
import numpy.matlib
from ppclass import pp
from numpy import loadtxt
from scipy import integrate

### Petite fonction utilisee plus tard pour convertir les tableaux DYNAMICO dans la taille des tableaux du modele radiatif (La fonction cherche l indice de la valeur de latitude la plus proche de celle prise)
### P.S: On a plus besoin de effectuer des calculs sur la meme grille mais je laisse la fonction ici au ou elle est utile pour des travaux futurs.
def find_x(array,x):
    y = int(np.argmin(np.abs(np.ceil(array[None].T - x)),axis=0))
    return int(y)




###Initialisation du domaine, des profils vertical et latitudinal + caracteristiques planete

#Recuperation des coefficients pour calculer le champ de pression 3D
coeffbp = loadtxt("coeffbp.txt", delimiter=",")

fff = "../simulations_dynamiques/Xhistins_170.nc" #"../../aslmd/DYNAMICO_SATURN/Xhistins_172.nc" #"../simulations_dynamiques/Xhistins_170.nc"
time_step="761040,6600000000" #309

### VARIABLES DYNAMICO: Toutes les variables sont ici en 3D

lsub = pp(file=fff,var="ls",x=0,y=0,verbose=True).getf()
ls=(lsub/np.pi)*180.
u,lon,lat,zmod,t = pp(file=fff,var="u",t=time_step,verbose=True).getfd()
v,lon,lat,zmod,t = pp(file=fff,var="v",t=time_step,verbose=True).getfd()
omega,lon,lat,zmod,t = pp(file=fff,var="omega",t=time_step,verbose=True).getfd()
temp,lon,lat,zmod,t = pp(file=fff,var="temperature",t=time_step,verbose=True).getfd()
psmod=pp(file=fff,var="ps",t=time_step,verbose=True).get()


#fffrad = "../SIMUS/test2d_saturne/diagfi.nc"
#time_steprad="3570,3600"

### Variables modele radiatif
#tetarad,lonrad,latrad,zmodrad,trad = pp(file=fffrad,var="teta",t=time_steprad,x="-180,180",verbose=True,changetime="correctls").getfd()
#zdtsw= pp(file=fffrad,var="zdtsw",t=time_steprad,x="-180,180",verbose=True,changetime="correctls").get()
#zdtlw= pp(file=fffrad,var="zdtlw",t=time_steprad,x="-180,180",verbose=True,changetime="correctls").get()

# ,logy=1 ,useindex="1000" ,changetime="correctls"


#Pour verifier que le ls est bien celui mis dans les plot (fixe)
print 'Interval de LS =', min(ls), max(ls)


#########  IMPORTANT  #########
#### Selectionner ici la tranche d'altitude a surveiller :
#Varie de 0 a 32 
infz=16    #borne inf de Z ###0   ### 16=100mb (base stratosphere) base conseillee
topz=30    #borne sup de Z ###32   ### 30=5mb (haut de Dynamico) plafond conseillee


nn=1   # Indice pour retirer les poles sur les latitudes 
#1=-89.25,89.25; 50=-64.25,64.25

lat=lat[nn:-nn,:]
#latrad=latrad[nn:-nn]
#tetarad=tetarad[infz:topz,nn:-nn]

phi = lat[:,0]*np.pi/180.
sinphi = np.sin(phi)
cosphi = np.cos(phi)   #Transformation des latitudes en phi

##### NOTA BENE: Il n y a pas de champ 2D d altitude dans DYNAMICO, il faut le faire a la main avec un fichier .txt fourni par Spiga et multiplie par pression de surface
#z2D=zmod

###Calcul du champ de pression 3D depuis les coefficient (plus precis que les autres estimations)
psmod.f=psmod.f[nn:-nn,:]
coeffbp2D=np.matlib.repmat(coeffbp,np.shape(lat)[0],1)
coeffbp3D=np.tile(coeffbp2D,(np.shape(lon)[1],1,1)).T
p2D=coeffbp3D*psmod.f
temp=temp[:,nn:-nn,:]

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
#p=np.empty(np.shape(z))
#p2D[0,:]=p0
#p[0]=p0
#for i in range(1,np.shape(p2D)[0]):
#for i in range(1,np.shape(p)[0]):
#  p[i]=p[i-1]/np.exp((g/(Rs*(np.mean(temp[i,:]+temp[i-1,:])/2.)))*(z[i]-z[i-1])*1000.)
  #for j in range(1,np.shape(p2D)[1]):
    #p2D[i,j]=p2D[i-1,j]/np.exp((g/(Rs*((temp[i,j]+temp[i-1,j])/2.)))*(z2D[i,j]-z2D[i-1,j])*1000.)


### ICI on decoupe tout les champs suivant l'altitude et la latitude + creation de vecteur 1D depuis 2D ou champd 2D depuis 3D
p2D=p2D[infz:topz,:,:]
temp=temp[infz:topz,:,:]
ro=ro[infz:topz,:,:]
p=np.mean(p2D,1)
#p=p[infz:topz]
u=u[infz:topz,nn:-nn]
v=v[infz:topz,nn:-nn]
omega=omega[infz:topz,nn:-nn]
z=zmod[infz:topz]
omega2D=np.mean(omega,2)
w=-omega/(ro*g)
u2D=np.mean(u,2)
v2D=np.mean(v,2)
w2D=np.mean(w,2)
temp2D=np.mean(temp,2) 
#Qrad=zdtsw+zdtlw
#Q=Qrad.f[infz:topz,nn:-nn]
#phirad = latrad*np.pi/180.
#sinphirad = np.sin(phirad)
#cosphirad = np.cos(phirad)


##### Convertir les tableaux du modele DYNAMICO dans la taille de ceux radiatif
# Si un jour besoin de comparer.
#ind=np.zeros(np.shape(latrad)[0])
#for i in range(np.shape(latrad)[0]):
#  ind[i]=find_x(lat,latrad[i])


#ind=ind.astype(int)
#ind2D=np.matlib.repmat(ind,np.shape(z)[0],1).astype(int)
#tempind=temp.T[ind]
#tempind=tempind.T
#p2Dind=p2D.T[ind]
#p2Dind=p2Dind.T
#roind=ro.T[ind]
#roind=roind.T

#creation de teta depuis temp
(nz,ny,nx)=np.shape(temp)
teta=(nz,ny,nx)
teta=np.empty(teta)
teta=temp*((p0/p2D)**(Rs/cp))
teta2D=np.mean(teta,2) 


#(nzi,nxi)=np.shape(tempind)
#tetaind=(nzi,nxi)
#tetaind=np.empty(tetaind)
#tetaind=tempind*((p0/p2Dind)**(Rs/cp)) 


####### ICI on cree les variables necessaire aux calculs des termes de forcage par tourbillons
#Pour moyenner zonalement et mettre la meme valeur sur chaque point en longitude , la fonction tile permet de cree un vecteur 3D depuis un champ 2D avec la meme valeur partout sur le 3eme axe.
#notation explique sur la ligne suivante
#b=barre,p=prime,e=etoile
vb2D=np.mean(v,2)
vb=np.tile(vb2D.T,(np.shape(lon)[1],1,1)).T
#vb=np.matlib.repmat(vb2D,np.shape(lat)[0],1).T

wb2D=np.mean(w,2)
wb=np.tile(wb2D.T,(np.shape(lon)[1],1,1)).T

tetab2D=np.mean(teta,2)
tetab=np.tile(tetab2D.T,(np.shape(lon)[1],1,1)).T
#tetab=np.matlib.repmat(tetab1D,np.shape(lat)[0],1).T

dtetabdz2D,dtetabdphi2D=np.gradient(tetab2D,z*1000.,phi, edge_order=2)
dtetabdz=np.tile(dtetabdz2D.T,(np.shape(lon)[1],1,1)).T
dtetabdphi=np.tile(dtetabdphi2D.T,(np.shape(lon)[1],1,1)).T

vp=v-vb
wp=w-wb
tetap=teta-tetab

vptetapb2D=np.mean(vp*tetap,2)
vptetapb=np.tile(vptetapb2D.T,(np.shape(lon)[1],1,1)).T
wptetapb2D=np.mean(wp*tetap,2)
wptetapb=np.tile(wptetapb2D.T,(np.shape(lon)[1],1,1)).T


#####A partir de la on calcul les termes en 3D
vterm=(ro*vptetapb)/dtetabdz
vterm2D=np.mean(vterm,2)
dvtermdz2D,dvtermdphi2D=np.gradient(vterm2D,z*1000.,phi, edge_order=2)
dvtermdz=np.tile(dvtermdz2D.T,(np.shape(lon)[1],1,1)).T
vterm=-(1./ro)*dvtermdz
vterm2D=np.mean(vterm,2)
vbe=vb+vterm
vbe2D=np.mean(vbe,2)

cosphi2D=np.matlib.repmat(cosphi,np.shape(z)[0],1)
cosphi3D=np.tile(cosphi2D.T,(np.shape(lon)[1],1,1)).T

wterm=(cosphi3D*vptetapb)/dtetabdz
wterm2D=np.mean(wterm,2)
dwtermdz2D,dwtermdphi2D=np.gradient(wterm2D,z*1000.,phi, edge_order=2)
dwtermdphi=np.tile(dwtermdphi2D.T,(np.shape(lon)[1],1,1)).T
wterm=(1./a*cosphi3D)*dwtermdphi
wterm2D=np.mean(wterm,2)
wbe=wb+wterm
wbe2D=np.mean(wbe,2)

###Rectified eddy forcing term

tetaterm=ro*(((vptetapb*dtetabdphi)/a*dtetabdz)+wptetapb)
tetaterm2D=np.mean(tetaterm,2)
dtetatermdz2D,dtetatermdphi2D=np.gradient(tetaterm2D,z*1000.,phi, edge_order=2)
dtetatermdz=np.tile(dtetatermdz2D.T,(np.shape(lon)[1],1,1)).T
tetaterm=-(1./ro)*dtetatermdz
tetaterm2D=np.mean(tetaterm,2)
reft=tetaterm2D


print min(p[:,0])
print max(p[:,0])

plt.figure()
plt.contourf(lat[:,0],p[:,0],u2D,levels=np.arange(-np.max(abs(u2D)),np.max(abs(u2D)),step=np.max(abs(u2D))/100.),cmap=plt.cm.seismic)
plt.title('U Ls=330-345')
plt.colorbar()
plt.yscale('log')
plt.ylabel('Pressure (Pa)')
plt.xlabel('Latitude')
plt.gca().invert_yaxis()
#plt.savefig('V_final')

plt.figure()
plt.contourf(lat[:,0],p[:,0],v2D,levels=np.arange(-np.max(abs(v2D)),np.max(abs(v2D)),step=np.max(abs(v2D))/100.),cmap=plt.cm.seismic)
plt.title('V Ls=330-345')
plt.colorbar()
plt.yscale('log')
plt.ylabel('Altitude km ')
plt.xlabel('Latitude')
plt.gca().invert_yaxis()
#plt.savefig('V_final')

plt.figure()
plt.contourf(lat[:,0],p[:,0],w2D,levels=np.arange(-np.max(abs(w2D)),np.max(abs(w2D)),step=np.max(abs(w2D))/100.),cmap=plt.cm.seismic)
plt.title('W Ls=330-345')
plt.colorbar()
plt.yscale('log')
plt.ylabel('Altitude km ')
plt.xlabel('Latitude')
plt.gca().invert_yaxis()
#plt.savefig('V_final')

plt.figure()
plt.contourf(lat[:,0],p[:,0],temp2D,levels=np.linspace(np.min(temp),np.max(temp),100),cmap=plt.cm.jet)
plt.colorbar()
CS=plt.contour(lat[:,0],p[:,0],temp2D,levels=np.round(np.linspace(80,160,9)),colors='black')
plt.clabel(CS, inline=0.1, fontsize=10,fmt='% 3d')
plt.title('Temperature (K) Ls=330-345' ) 
plt.ylabel('Pressure (Pa)')
plt.xlabel('Latitude')
plt.yscale('log')
plt.gca().invert_yaxis()

#plt.figure()
#plt.contourf(latrad,p,tempind,levels=np.linspace(np.min(tempind),np.max(tempind),100))
#plt.colorbar()
#CS=plt.contour(latrad,p,tempind,levels=np.round(np.linspace(80,160,9)),colors='black')
#plt.clabel(CS, inline=0.1, fontsize=10,fmt='% 3d')
#plt.title('Temperature (Indice Rad)(K) Ls=340-355' ) 
#plt.ylabel('Pressure (Pa)')
#plt.xlabel('Latitude')
#plt.yscale('log')
#plt.gca().invert_yaxis()

plt.figure()
plt.contourf(lat[:,0],p[:,0],teta2D,levels=np.linspace(np.min(teta),np.max(teta),100),cmap=plt.cm.jet)
plt.colorbar()
plt.title('Teta (K) Ls=330-345' ) 
plt.ylabel('Pressure (Pa)')
plt.xlabel('Latitude')
plt.yscale('log')
plt.gca().invert_yaxis()

#plt.figure()
#plt.contourf(latrad,p,Q,levels=np.linspace(np.min(Q),np.max(Q),100))
#plt.colorbar()
#plt.title('Chauffage radiatif Ls=340-355' ) 
#plt.ylabel('Pressure (Pa)')
#plt.xlabel('Latitude')
#plt.yscale('log')
#plt.gca().invert_yaxis()

plt.figure(figsize=(14,8))
plt.subplot(2, 2, 1)
plt.contourf(lat[:,0],p[:,0],wterm2D,levels=np.arange(-np.max(abs(wterm2D)),np.max(abs(wterm2D)),step=np.max(abs(wterm2D))/100.),cmap=plt.cm.seismic)
plt.colorbar(format='%.0e')
CS=plt.contour(lat[:,0],p[:,0],abs(wterm2D),levels=[0.000008],colors='black',linewidths=0.4,linestyles='dashed')
plt.clabel(CS, inline=0.1, fontsize=1)#,fmt='% 3d')
plt.title('Terme de forcage W Ls=330-345' ) 
plt.ylabel('Pressure (Pa)')
#plt.xlabel('Latitude')
plt.yscale('log')
plt.gca().invert_yaxis()

#plt.figure()
plt.subplot(2, 2, 2)
plt.contourf(lat[:,0],p[:,0],vterm2D,levels=np.arange(-np.max(abs(vterm2D)),np.max(abs(vterm2D)),step=np.max(abs(vterm2D))/100.),cmap=plt.cm.seismic)
plt.colorbar(format='%.0e')
CS=plt.contour(lat[:,0],p[:,0],abs(vterm2D),levels=[0.0055],colors='black',linewidths=0.4,linestyles='dashed')
plt.clabel(CS, inline=0.1, fontsize=1)#,fmt='% 3d')
plt.title('Terme de forcage V Ls=330-345' ) 
plt.ylabel('Pressure (Pa)')
#plt.xlabel('Latitude')
plt.yscale('log')
plt.gca().invert_yaxis()

#plt.figure()
plt.subplot(2, 1, 2)
plt.contourf(lat[:,0],p[:,0],reft,levels=np.arange(-np.max(abs(reft)),np.max(abs(reft)),step=np.max(abs(reft))/100.),cmap=plt.cm.seismic)
plt.colorbar(format='%.0e')
CS=plt.contour(lat[:,0],p[:,0],abs(vterm2D),levels=[0.00000004],colors='black',linewidths=0.4,linestyles='dashed')
plt.clabel(CS, inline=0.1, fontsize=1)#,fmt='% 3d')
plt.title('Terme de forcage thermodynamique (K.s-1) Ls=330-345' ) 
plt.ylabel('Pressure (Pa)')
plt.xlabel('Latitude')
plt.yscale('log')
plt.gca().invert_yaxis()

plt.figure()
plt.contourf(lat[:,0],p[:,0],vbe2D,levels=np.arange(-np.max(abs(vbe2D)),np.max(abs(vbe2D)),step=np.max(abs(vbe2D))/100.),cmap=plt.cm.seismic)
plt.colorbar()
plt.title('moyzonal(v)* Ls=330-345' ) 
plt.ylabel('Pressure (Pa)')
plt.xlabel('Latitude')
plt.yscale('log')
plt.gca().invert_yaxis()

plt.figure()
plt.contourf(lat[:,0],p[:,0],wbe2D,levels=np.arange(-np.max(abs(wbe2D)),np.max(abs(wbe2D)),step=np.max(abs(wbe2D))/100.),cmap=plt.cm.seismic)
plt.colorbar()
plt.title('moyzonal(w)* Ls=330-345' ) 
plt.ylabel('Pressure (Pa)')
plt.xlabel('Latitude')
plt.yscale('log')
plt.gca().invert_yaxis()

vbe2Dnan=np.empty(np.shape(vbe2D))
vbe2Dnan[:,:]=vbe2D[:,:]
vbe2Dnan[:,138:219]=np.NaN

wbe2Dnan=np.empty(np.shape(wbe2D))
wbe2Dnan[:,:]=wbe2D[:,:]
wbe2Dnan[:,138:219]=np.NaN

#[:,138:219] lat=20/-20
#[:,117:238]  lat=30/-25



plt.figure(figsize=(14,8))
plt.subplot(2, 1, 1)
plt.contourf(lat[:,0],p[:,0],vbe2Dnan,levels=np.arange(-np.nanmax(abs(vbe2Dnan)),np.nanmax(abs(vbe2Dnan)),step=np.nanmax(abs(vbe2Dnan))/100.),cmap=plt.cm.seismic)
plt.colorbar(format='%.0e')
plt.title('moyzonal(v)* (m.s-1) Ls=330-345' ) 
plt.ylabel('Pressure (Pa)')
#plt.xlabel('Latitude')
plt.yscale('log')
plt.gca().invert_yaxis()

#plt.figure(figsize=(10,5))
plt.subplot(2, 1, 2)
plt.contourf(lat[:,0],p[:,0],wbe2Dnan,levels=np.arange(-np.nanmax(abs(wbe2Dnan)),np.nanmax(abs(wbe2Dnan)),step=np.nanmax(abs(wbe2Dnan))/100.),cmap=plt.cm.seismic)
plt.colorbar(format='%.0e')
plt.title('moyzonal(w)* (m.s-1) Ls=330-345' ) 
plt.ylabel('Pressure (Pa)')
plt.xlabel('Latitude')
plt.yscale('log')
plt.gca().invert_yaxis()



plt.show(block=False)
input("Hit it and quit it! (press enter to close all)")
plt.close()
