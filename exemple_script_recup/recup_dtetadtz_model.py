from ppclass import pp
import ppplot
import ppcompute
import numpy as np

fff = "temporaire.nc"
teta,lon,lat,z,t = pp(file=fff,var="teta",t=500,x="-180,180",verbose=True,useindex="1000").getfd()
phi = lat*np.pi/180.
sinphi = np.sin(phi)
cosphi = np.cos(phi)
dtetadphi,dtetadz = ppcompute.deriv2d(teta,phi,z)
pl = ppplot.plot2d()
pl.verbose = True
pl.logy = 1
pl.f = dtetadz
pl.c = teta 
pl.x = lat
pl.y = z
pl.ylabel = 'altitude (km)'
pl.xlabel = 'latitude ($^{\circ}$N)'
pl.title = '$dteta/dz$'
pl.makeshow()

