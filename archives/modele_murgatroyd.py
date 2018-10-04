import pylab as py
import numpy as np
import ppcompute
import matplotlib.pyplot as plt
import numpy.matlib
from ppclass import pp

###Initialisation du domaine, des profils vertical et latitudinal + caracteristiques planete

fff = "../../simulations_nodyn/diagfi_saturn_rings.nc"
#../../simulations_nodyn/temporaire.nc
#../../simulations_nodyn/diagfi_saturn_rings.nc #rings
#../../simulations_nodyn/diagfi_saturn_norings.nc #norings

#fff = "../../simulations_nodyn/temporaire.nc"
teta,lon,lat,zmod,t = pp(file=fff,var="teta",t=500,x="-180,180",verbose=True,useindex="1000").getfd()
zdtsw= pp(file=fff,var="zdtsw",t=500,x="-180,180",verbose=True,useindex="1000").get()
zdtlw= pp(file=fff,var="zdtlw",t=500,x="-180,180",verbose=True,useindex="1000").get()
pmod= pp(file=fff,var="p",t=500,x="-180,180",verbose=True,useindex="1000").get()
tempmod= pp(file=fff,var="temp",t=500,x="-180,180",verbose=True,useindex="1000").get()
# ,logy=1 t="420,500" t=500
phi = lat*np.pi/180.
sinphi = np.sin(phi)
cosphi = np.cos(phi)   #Transformation des latitudes en phi
a=58210000. #Rayon de la planete en m

###Creation d un champ radiatif 2D, refroidit a gauche et rechauffe a droite: On cree un linspace au sommet de l'atmosphere et on le fait tendre vers 0 vers le bas
Qrad=zdtsw+zdtlw
Q=Qrad.f
z=zmod[::-1]

plt.figure()
plt.contourf(lat,z,Q,levels=np.linspace(np.min(Q),np.max(Q),100))
plt.title('Chauffage et refroidissement radiatif')
plt.colorbar()
plt.yscale('log')


###Creation d un champ de temperature potentielle teta qui varie sur la latitude: On cree la premiere et derniere colonne et on lisse avec un linspace

plt.figure()
plt.contourf(lat,z,teta,levels=np.linspace(np.min(teta),np.max(teta),100))
plt.title('Temperature potentielle')
plt.colorbar()
plt.yscale('log')


###Creation d'un champ de pression et temperature a partir de Z et teta

p2D=pmod.f
p=p2D[:,0]
temp=tempmod.f

plt.figure()
plt.contourf(lat,z,temp,levels=np.linspace(np.min(temp),np.max(temp),100))
plt.title('Temperature usuelle')
plt.colorbar()
plt.yscale('log')


###Calcul des derivees de teta et temp en phi,z et p

dtetadp,dtetadphi=np.gradient(teta,p,phi, edge_order=2)#deriv(teta,phi,p)
dtetadz,dtetadphi=np.gradient(teta,z*1000.,phi, edge_order=2)
dtempdz,dtempdphi=np.gradient(temp,z*1000.,phi, edge_order=2)

(nz,nx)=np.shape(teta)
dtetadp=dtetadp+0.000000000000001
###Initialisation de w(omega) premiere iteration avec v=0
wini=np.empty((nz,nx))
mask = np.where(dtetadp != 0)
wini[mask]=Q[mask]/((1./teta[mask])*dtetadp[mask]*temp[mask])
vini=np.ones((nz,nx))
epsy=np.ones((nz,1))
dw=np.ones((nz,nx))

dwinidp,dwinidphi = np.gradient(wini,p,phi, edge_order=2)

###Initialisation de v depuis wini(v=0) : Correspond a la premiere iteration sachant que epsy= terme correctif sur l integral pour faire tendre la vitesse du pole oppose a 0

for i in range(nz):
  epsy[i]=-(np.trapz(dwinidp[i,nx-1]*cosphi[0:nx+1],phi[0:nx+1],axis=0))/(np.trapz(cosphi[0:nx+1],phi[0:nx+1],axis=0))
  #epsy[i]=0.
  for j in range(nx):   
    integrand = (dwinidp[i,j]+epsy[i])*cosphi[0:j+1]
    vini[i,j]=np.trapz(integrand,phi[0:j+1],axis=0)
    vini[i,j]=-(a/cosphi[j])*vini[i,j]
    tab=np.gradient(p)
    dw[i,j]=epsy[i]*tab[i]


Q=(vini[:,:]/a*dtempdphi)+(wini[:,:]+dw[:,:])*(1/teta[:,:]*dtetadp*temp[:,:])
#Variables necessaires a la boucle + Initialisation des variables

h=1 #=nombre d iterations, on a deja "une" iteration en v=0
v_next = np.ones((nz,nx))
v_prev = np.ones((nz,nx))
eps = np.ones(nz)
v_prev[:,:] = vini[:,:]


# !!!!!!!!! Choix de la precision et du pourcent a comparer dans la boucle
percent = 1.
precision = 100000.

plt.figure()
plt.contourf(lat,z,dtetadp,levels=np.linspace(np.min(dtetadp),np.max(dtetadp),100))
plt.title('dtetadp')
plt.colorbar()
plt.yscale('log')

plt.figure()
plt.contourf(lat,z,p2D,levels=np.linspace(np.min(p2D),np.max(p2D),100))
plt.title('p')
plt.colorbar()
plt.yscale('log')

plt.figure()
plt.contourf(lat,z,wini,levels=np.linspace(np.min(wini),np.max(wini),100))
plt.title('W(omega) initial')
plt.colorbar()
plt.yscale('log')

plt.figure()
plt.contourf(lat,z,dwinidp,levels=np.linspace(np.min(dwinidp),np.max(dwinidp),100))
plt.title('dwinidp')
plt.colorbar()
plt.yscale('log')

plt.figure()
plt.contourf(lat,z,dwinidphi,levels=np.linspace(np.min(dwinidphi),np.max(dwinidphi),100))
plt.title('dwinidphi')
plt.colorbar()
plt.yscale('log')

plt.figure()
plt.contourf(lat,z,vini,levels=np.linspace(np.min(vini),np.max(vini),100))
plt.title('vini')
plt.colorbar()
plt.yscale('log')

plt.figure()
plt.contourf(lat,z,Q,levels=np.linspace(np.min(Q),np.max(Q),100))
plt.title('Q1')
plt.colorbar()
plt.yscale('log')

#plt.show(block=False)
#input("Hit it and quit it! (press enter to close all)")
#plt.close()
#exit()

#Boucle While, calcul V et W tant que cela converge suivant la precision et le percent
#while precision >= percent/100.:
for h in range(2,5):  #Pour controler l iteration (decommentez h=h+1 dans la boucle)
  w_next=(Q[:,:]-(v_prev[:,:]*dtempdphi/a))/((1/teta)*dtetadp*temp)
  dwdp,dwdphi=np.gradient(w_next,p,phi, edge_order=2)
  for i in range(nz):
    eps[i]=-(np.trapz(dwdp[i,nx-1]*cosphi[0:nx+1],phi[0:nx+1],axis=0))/(np.trapz(cosphi[0:nx+1],phi[0:nx+1],axis=0))
    #eps[i]=0.
    for j in range(nx):   
      v_next[i,j]=np.trapz((dwdp[i,j]+eps[i])*cosphi[0:j+1],phi[0:j+1],axis=0) 
      v_next[i,j]=-(a/cosphi[j])*v_next[i,j]
      tab=np.gradient(p)
      dw[i,j]=eps[i]*tab[i]
  w = np.where(np.abs(v_next[:,:])-np.abs(v_prev[:,:]) != 0)
  precision = np.mean(np.abs(v_next[w] - v_prev[w])/np.abs(0.5*(v_next[w]+v_prev[w])))
  print precision, np.max(v_next[w]), np.min(v_next[w])
  Q=(v_next[:,:]/a*dtempdphi)+(w_next[:,:]+dw[:,:])*(1/teta[:,:]*dtetadp*temp[:,:])
  v_prev[:,:] = v_next[:,:]
  #h += 1
#Fin de boucle

#print precision*100.,percent
print('On a realise ',h,'iterations')
#Tout les plots des resultats!

plt.figure()
plt.contourf(lat,z,w_next,levels=np.linspace(np.min(w_next),np.max(w_next),100))
plt.title('W(omega) final')
plt.colorbar()
plt.yscale('log')

plt.figure()
plt.contourf(lat,z,v_next,levels=np.linspace(np.min(v_next),np.max(v_next),100))
plt.title('V final')
plt.colorbar()
plt.yscale('log')
#plt.savefig('V_final')

plt.figure()
plt.contourf(lat,z,vini-v_next,levels=np.linspace(np.min(vini-v_next),np.max(vini-v_next),100))
plt.title('Diff Vini-Vfinal')
plt.colorbar()
plt.yscale('log')

plt.figure()
plt.contourf(lat,z,dwdp,levels=np.linspace(np.min(dwdp),np.max(dwdp),100))
plt.title('dwdp')
plt.colorbar()
plt.yscale('log')

plt.figure()
plt.contourf(lat,z,Q,levels=np.linspace(np.min(Q),np.max(Q),100))
plt.title('Qfin')
plt.colorbar()
plt.yscale('log')


plt.show(block=False)
input("Hit it and quit it! (press enter to close all)")
plt.close()
