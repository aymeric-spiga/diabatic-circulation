from ppclass import pp

teta = pp()
teta.verbose = True
teta.logy = 1
teta.file ="temporaire.nc"
teta.var ="teta"
teta.x ="-180,180"
#teta.t ="2"
teta.changetime="correctls"
teta.y="-90,90"
teta.get()
teta.printme()
teta.ylabel= "teta(K)"
teta.title = "Temperature potentielle"
teta.defineplot()
teta.makeplot()



#VARIABLES:  [u'Time', u'controle', u'latitude', u'longitude', u'altitude', u'aps', u'bps', u'ap', #u'bp', u'soildepth', u'aire', u'phisinit', u'Ls', u'tsurf', u'ps', u'temp', u'teta', u'p', #u'ISR', u'ASR', u'OLR', u'zdtsw', u'zdtlw']


