"""
"""

import numpy as np
import astropy.io.fits as fits
import h5py

def GALFAHI_cubeinfo(fitsfile):
    from astropy.io import fits
    import numpy as np
   
    data = fits.getdata(fitsfile)
    header = fits.getheader(fitsfile)
    ra = header['CRVAL1'] + header['CDELT1'] * (
                np.arange(data.shape[2])+1 - header['CRPIX1'])             # degree
    dec = header['CRVAL2'] + header['CDELT2'] * (
                np.arange(data.shape[1])+1 - header['CRPIX2'])             # degree
    vlsr = (header['CRVAL3'] + header['CDELT3'] * (
                np.arange(data.shape[0])+1 - header['CRPIX3']))*(10**(-3)) # km/s
    delv = header['CDELT3']*(10**(-3))                                     # km/s
    return data, header, vlsr, dec, ra, delv


### first cube
this_cube_ra1 = ''
this_cube_dec1 = ''
ifile1 = 'GALFAdata/data/GALFA_HI_RA+DEC_%06.2f+%05.2f_W.fits' % (this_cube_ra1, this_cube_dec1)
cubeinfo1 = GALFAHI_cubeinfo(ifile1)

data1 = cubeinfo1[0]
vlsr1 = cubeinfo1[2]
dec1 = cubeinfo1[3]
ra1 = cubeinfo1[4]
indv = np.all([vlsr1>-500.00, vlsr1<=100], axis=0)
data1 = data1[indv]

### second cube
this_cube_ra2 = ''
this_cube_dec2 = ''
ifile2 = 'GALFAdata/data/GALFA_HI_RA+DEC_%06.2f+%05.2f_W.fits' % (this_cube_ra2, this_cube_dec2)
cubeinfo2 = GALFAHI_cubeinfo(ifile2)

data2 = cubeinfo2[0]
vlsr2 = cubeinfo2[2]
dec2 = cubeinfo2[3]
ra2 = cubeinfo2[4]
data2 = data2[indv]



ind = dec2>dec1.max()
dec2 = dec2[ind]
data2 = data2[:, ind, :] 

vlsr = vlsr1
ra = ra1
dec = np.concatenate([dec1, dec2])
data = np.concatenate([data1, data2], axis=1)
# indra = ra>=22.5

m33cube = data
m33ra = ra
m33dec = dec
m33v = vlsr[indv]
f = h5py.File('GALFAdata/specialData/M33_WrightC_DR2.h5', 'w') ###look at the names of these

f.create_dataset('m33cube', data=m33cube)
f.create_dataset('ra', data=m33ra)
f.create_dataset('dec', data=m33dec)
f.create_dataset('vlsr', data=m33v)

f.close()