#!/usr/bin/env python
# coding: utf-8
"""Download the SDSS DR7 spectra of the quasars paired with hosts that have a spectrum.

Spectra are fetched by plate-MJD-fiber from the SDSS science archive (DR7 plates,
run2d 26) and stored in spectra/qso/ in the same format as the galaxy spectra:
vacuum wavelength [A] and f_lambda with its error [erg/s/cm2/A].

usage: fetch_qso_spectra.py chimeras-fullinfo.fits galaxy-spectra.fits qso-spectra.fits
"""
import io
import os
import sys
import urllib.request

import numpy as np
from astropy.io import fits
from astropy.table import Table

OUTDIR = 'spectra/qso'
SAS = 'https://data.sdss.org/sas/dr17/sdss/spectro/redux/26/spectra/lite/%04d/spec-%04d-%05d-%04d.fits'
SPEC_RES = 2000


def fetch(plate, mjd, fiber, z, path):
    if not os.path.exists(path):
        with urllib.request.urlopen(SAS % (plate, plate, mjd, fiber), timeout=300) as response:
            data = fits.open(io.BytesIO(response.read()))[1].data
        ivar = np.asarray(data['ivar'], dtype=float)
        err = np.where(ivar > 0, 1 / np.sqrt(np.where(ivar > 0, ivar, 1)), np.inf) * 1e-17
        tab = Table(dict(wave_vac=10**np.asarray(data['loglam'], dtype=float),
                         flux=np.asarray(data['flux'], dtype=float) * 1e-17, flux_err=err))
        tab['wave_vac'].unit = 'Angstrom'
        tab['flux'].unit = tab['flux_err'].unit = 'erg / (s cm2 Angstrom)'
        tab.meta.update(SURVEY='SDSS-DR7', SPECID='%d-%d-%d' % (plate, mjd, fiber), SPEC_RES=SPEC_RES, Z=z)
        tab.write(path, overwrite=True)
    return path


data = Table.read(sys.argv[1])
galaxies = Table.read(sys.argv[2])
data = data[np.isin(data['ID_COSMOS'], galaxies['ID_COSMOS'])]
_, first = np.unique(np.asarray(data['ID_SDSS']), return_index=True)
qsos = data[first]
print('quasars paired with hosts that have a spectrum:', len(qsos))
os.makedirs(OUTDIR, exist_ok=True)

rows = []
for q in qsos:
    plate, mjd, fiber = int(q['Plate_QSO']), int(q['MJD_QSO']), int(q['Fiber_QSO'])
    path = '%s/spec-%04d-%05d-%04d.fits' % (OUTDIR, plate, mjd, fiber)
    try:
        fetch(plate, mjd, fiber, float(q['z_QSO']), path)
    except Exception as e:
        print('  failed %s (%d-%d-%d): %s' % (q['ID_SDSS'], plate, mjd, fiber, e))
        continue
    rows.append((q['ID_SDSS'], plate, mjd, fiber, float(q['z_QSO']), path))
out = Table(rows=rows, names=['ID_SDSS', 'plate', 'mjd', 'fiber', 'redshift_QSO', 'spectrum_file'])
print('quasar spectra:', len(out))
out.write(sys.argv[3], overwrite=True)
