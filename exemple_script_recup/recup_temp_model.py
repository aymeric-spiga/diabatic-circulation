from ppclass import pp

temp = pp()
temp.verbose = True
#temp.logy = 1
temp.file = "../temporaire.nc"
temp.var = "temp"
temp.x ="-180,180"
#temp.t = "2"
temp.changetime="correctls"
temp.y="-90,90"
temp.get()
temp.printme()
temp.ylabel= "temperature(K)"
temp.title = "Temperature"
temp.defineplot()
temp.makeplot()



#VARIABLES:  [u'Time', u'controle', u'latitude', u'longitude', u'altitude', u'aps', u'bps', u'ap', #u'bp', u'soildepth', u'aire', u'phisinit', u'Ls', u'tsurf', u'ps', u'temp', u'teta', u'p', #u'ISR', u'ASR', u'OLR', u'zdtsw', u'zdtlw']


