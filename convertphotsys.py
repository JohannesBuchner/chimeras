#!/usr/bin/env python
# coding: utf-8
import numpy as np
from astropy.table import Table
import tqdm
import matplotlib.pyplot as plt
import sys

def ab_flux(ABmag):
    """convert AB magnitude to flux in uJy"""
    return (10**(23-(((ABmag)+48.6)/2.5))) * 1000000

def flux_ab(flux):
    return 2.5*(23-np.log10(flux/1000000)) - 48.6

dataC = Table.read(sys.argv[1])
offs = Table.read('GALOI15PP_SED.fits')
print('Model table has %d entries' % len(offs))

# new, corrected flux columns:
dataC['u_sdss'] = -999.9
dataC['r_sdss'] = -999.9
dataC['i_sdss'] = -999.9
dataC['z_sdss'] = -999.9
dataC['u_sdss_err'] = -999.9
dataC['r_sdss_err'] = -999.9
dataC['i_sdss_err'] = -999.9
dataC['z_sdss_err'] = -999.9
# no conversion for the following bands, they are quite similar and the errors are more important
dataC['J_2mass'] = dataC['WFCAM_J']
dataC['J_2mass_err'] = dataC['WFCAM_J_err']
dataC['H_2mass'] = dataC['CFHT_H']
dataC['H_2mass_err'] = dataC['CFHT_H_err']
dataC['Ks_2mass'] = dataC['CFHT_K']
dataC['Ks_2mass_err'] = dataC['CFHT_K_err']

print("looking for galaxy model for each object")
nomatch = []
for row in tqdm.tqdm(dataC):
    # in the model database, find the entry with the same model, rounded redshift, extinction, age:
    m = np.isclose(offs['Mod'], row['MODEL'])
    r = np.isclose(offs['redshift'], np.round(row['redshift'], 1))
    ex = np.logical_and(np.isclose(offs['E(B-V)'], row['EXTINCTION']), np.isclose(offs['Age'],row['AGE']))
    sel = np.logical_and(np.logical_and(m, r), ex)
    red = offs[sel]
    if len(red) > 0:
        # convert filter system from SUBARU to SDSS with model
        for band in 'z', 'i', 'r':
            model_mag_delta = red[band + '_sdss'] - red['SUBARU_' + band]
            row[band + '_sdss'] = ab_flux(flux_ab(row['SUBARU_' + band]) + model_mag_delta)
            ratio = row[band + '_sdss'] / row['SUBARU_' + band]
            row[band + '_sdss_err'] = row['SUBARU_%s_err' % band] * ratio

        # convert filter system from CFHT to SDSS with model
        delta_u = red['u_sdss'] - red['CFHT_u']
        row['u_sdss'] = ab_flux(flux_ab(row['CFHT_u']) + delta_u)
        u_ratio = row['u_sdss'] / row['CFHT_u']
        row['u_sdss_err'] = row['CFHT_u_err'] * u_ratio
        del u_ratio, delta_u
    else:
        nomatch.append(row['id'])



print("no match:", nomatch)


for band in ['u', 'r', 'i', 'z']:
    orig_band = 'SUBARU_' + band if band != 'u' else 'CFHT_u'
    dataC[band + '_offset'] = flux_ab(dataC[band+'_sdss']) - flux_ab(dataC[orig_band])
    offsets = dataC[np.isfinite(dataC[band+'_offset'])][band+'_offset']
    plt.figure()
    plt.hist(offsets, bins='auto', histtype='step')
    plt.xlabel("Magnitude difference new-old")
    plt.title('%s: (median: %.2f, 95%% are smaller than %.2f mag)' % (band, np.median(np.abs(offsets)), np.quantile(np.abs(offsets), 0.95)))
    plt.savefig('convertphotsys_band-%s.pdf' % band)
    plt.close()


dataC.write(sys.argv[2], format='fits', overwrite=True)

x = dataC['redshift']
y = flux_ab(dataC['i_sdss'])-flux_ab(dataC['z_sdss'])
plt.scatter(x,y, label = 'SUBARU')
plt.scatter(x,y, label = 'sdss')
plt.xlabel('Redshift')
plt.ylabel('i-z')
plt.savefig('convertphotsys_i-z.pdf')
plt.close()
