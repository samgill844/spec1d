from astropy.table import Table 

def create_log():
    files = glob.glob('a*.fits')
    files.sort() 

    time = [] # HJD-OBS
    target = [] # OBJECT
    #ra = [] # TARG-RA
    #dec = [] # TARG-DEC
    arc_lamp = [] # ARC-LAMP
    exptype = [] # EXPTYPE
    for i in range(len(files)):
        header = fits.open(files[i])[0].header 
        time.append(float(header['HJD-OBS']))
        target.append(header['OBJECT'])
        arc_lamp.append(header['ARC-LAMP'])
        exptype.append(header['EXPTYPE'])

    t = Table()
    t.add_column(Column(files, name='FILE'))
    t.add_column(Column(np.array(time) + 2400000, name='HJD'))
    t.add_column(Column(target, name='TARGET'))
    t.add_column(Column(arc_lamp, name='ARC-LAMP'))
    t.add_column(Column(exptype, name='EXPTYPE'))

    t.write('log.fits', overwrite=True)

