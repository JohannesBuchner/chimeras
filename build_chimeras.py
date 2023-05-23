#!/usr/bin/env python
# coding: utf-8

import numpy as np
import sys
from astropy.table import Table, vstack
import tqdm
import joblib

def ab_flux(ABmag):
    """convert AB magnitude to flux in uJy"""
    return (10**(23-(((ABmag)+48.6)/2.5))) * 1000000

def magerr2fluxerr(mag, magerr):
    """convert the lower (brighter) magnitude to flux, return the deviation in flux."""
    return ab_flux(mag - magerr) - ab_flux(mag)

def convert_flux_and_err_where_measured(mag, magerr):
    """convert safely where error is positive, return nan otherwise."""
    flux = np.where(magerr > 0, ab_flux(mag), np.nan)
    flux_err = np.where(magerr > 0, magerr2fluxerr(mag, magerr), np.nan)
    return flux, flux_err

def flux(c, q, w):
    """sum fluxes, with weight w"""
    return c + w * q

def fluxerr(ce, qe, w):
    """gaussian error propagation, with weight w"""
    return np.sqrt(ce**2 + (w*qe)**2)

# WEIGHTS
QSO_weights = [0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1]

bands = ['u_sdss', 'r_sdss', 'i_sdss', 'z_sdss', 'J_2mass', 'H_2mass', 'Ks_2mass', 'IRAC1', 'IRAC2']

def make_combinations(binss_orig, rowC, COSMOS_id, QSO_weight, QSO_IDs, colnames):
    """for the given QSO_weight compute summed fluxes and errors"""
    binss = binss_orig.copy()
    # assign unique ids for each weight
    binss["id"] = ['%s_%s_%s' % (j, COSMOS_id, QSO_weight) for j in QSO_IDs]
    binss["chimera_QSO_weight"] = QSO_weight
    # COMBINE FLUXES
    for band in bands:
        binss[band] = flux(rowC[band], binss[band], QSO_weight)
        band_err = band + "_err"
        binss[band_err] = fluxerr(rowC[band_err], binss[band_err], QSO_weight)

    return binss[colnames]

def make_chimeras_for(rowC, dataQ, GAL_columns_to_copy, colnames):
    """for the galaxy *rowC*, find matching quasars and create hybrids."""
    redshift = rowC['redshift']
    COSMOS_ID = rowC['id']
    closeby_candidates = np.abs(dataQ["redshift_QSO"] - redshift) < 0.01
    binss_orig = dataQ[closeby_candidates]
    # ADD COSMOS COLUMNS
    if len(binss_orig) == 0:
        print("\nno quasar match for COSMOS id %s at z=%.3f" % (rowC['id'], rowC['redshift']))
        return []
    else:
        binss_orig["redshift"] = redshift
        binss_orig["ID_COSMOS"] = COSMOS_ID
        for col in GAL_columns_to_copy:
            binss_orig[col + "_GAL"] = rowC[col]
        QSO_IDs = binss_orig["ID_SDSS"]

        #return [make_combinations(binss_orig, rowC, COSMOS_ID, QSO_weight, QSO_IDs, colnames) for QSO_weight in QSO_weights]
        #return joblib.Parallel(-1)(joblib.delayed(make_combinations)(binss_orig, rowC, COSMOS_ID, QSO_weight, QSO_IDs, colnames) for QSO_weight in QSO_weights)
        return [vstack(joblib.Parallel(2)(joblib.delayed(make_combinations)(binss_orig, rowC, COSMOS_ID, QSO_weight, QSO_IDs, colnames) for QSO_weight in QSO_weights))]

def main():
    # Import COSMOS galaxies
    dataC = Table.read(sys.argv[1])
    print("Galaxy input catalog:", len(dataC))

    for band in bands:
        # keep only sources where the uncertainty in each band is positive (detections)
        # this is needed to propagate the error bars, because upper limits would be more complicated to handle
        dataC = dataC[dataC[band + '_err'] > 0]
        del band

    print("after removing cases with non-positive uncertainties:", len(dataC))

    # Import QSO
    dataQ = Table.read(sys.argv[2])
    print("Quasar input catalog:", len(dataQ))
    dataQ.rename_column('z', "redshift_QSO")
    dataQ["ID_SDSS"] = dataQ["SDSS"].astype(str)
    for colname in dataQ.colnames:
        if colname.replace('_err', '') not in bands and not colname.endswith('_QSO') and not colname == 'ID_SDSS':
            # print("deleting column", colname)
            del dataQ[colname]
    
    GAL_columns_to_copy = []
    colnames = ["id", "redshift", "redshift_QSO", "ID_COSMOS", "ID_SDSS", "chimera_QSO_weight"]
    for band in bands:
        colnames += [band, band + '_err']
    
    print("Creating chimeras...")
    with open(sys.argv[3], 'w') as fout:
        for COSMOS_ID, t in zip(dataC['id'], joblib.Parallel(4)((joblib.delayed(make_chimeras_for)(rowC, dataQ, GAL_columns_to_copy, colnames) for rowC in tqdm.tqdm(dataC)))):
            if len(t) > 0:
                filename = 'chimeras/%s.fits' % COSMOS_ID
                del t[0]['redshift_QSO']
                t[0].write(filename, overwrite=True)
                fout.write(filename + "\n")
    

if __name__ == '__main__':
    main()
