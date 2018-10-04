from ppclass import pp

zdtsw = pp()
zdtsw.verbose = True
#zdtsw.logy = 1
zdtsw.file = "temporaire.nc"
zdtsw.var = "zdtsw"
zdtsw.x ="-180,180"
zdtsw.t = "416,500" #500
zdtsw.useindex = "1000"
#zdtsw.y="-90,90"
zdtsw.get()
zdtsw.printme()


print "test"
altitude = pp()
altitude.verbose = True
#zdtsw.logy = 1
altitude.file = "temporaire.nc"
altitude.var = "altitude"
altitude.get()
altitude.printme()

print altitude.f

zdtlw = pp()
zdtlw << zdtsw
zdtlw.var = "zdtlw"
zdtlw.get()
zdtlw.printme()

Qrad=zdtsw+zdtlw
Qrad.ylabel= "Log P"
Qrad.title = "Chauffage radiatif"
Qrad.defineplot()
Qrad.makeplot()



#VARIABLES:  [u'Time', u'controle', u'latitude', u'longitude', u'altitude', u'aps', u'bps', u'ap', #u'bp', u'soildepth', u'aire', u'phisinit', u'Ls', u'tsurf', u'ps', u'temp', u'teta', u'p', #u'ISR', u'ASR', u'OLR', u'zdtsw', u'zdtlw']


