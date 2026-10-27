#!/usr/bin/env python
# coding: utf-8
"""Find and download flux-calibrated optical spectra of the chimera host galaxies.

Sources (both absolutely flux calibrated):
  DESI DR1   (via NOIRLab SPARCL; vacuum wavelengths, R ~ 2000-5000)
  LEGA-C     (ESO Phase 3 archive; air wavelengths, R = 2500)

Each host (ID_COSMOS) is matched within MATCH_ARCSEC of its COSMOS position; DESI
spectra must also agree with the chimera redshift within MAX_DZ (LEGA-C has no
redshift in its archive metadata). Of all matches the spectrum with the highest
median S/N per pixel is kept. Spectra are stored in spectra/gal/ in a common format:
vacuum wavelength [A] and f_lambda with its error [erg/s/cm2/A].

usage: fetch_galaxy_spectra.py chimeras-fullinfo.fits galaxy-spectra.fits
"""
import io
import os
import sys
import time
import urllib.parse
import urllib.request

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
import astropy.units as u

MATCH_ARCSEC = 1.0
MAX_DZ = 0.01
OUTDIR = 'spectra/gal'
ESO_TAP = 'https://archive.eso.org/tap_obs/sync'
ESO_FILE = 'https://dataportal.eso.org/dataportal_new/file/'


def air_to_vacuum(wave):
    """Morton (2000) air-to-vacuum conversion [A]."""
    s2 = (1e4 / wave)**2
    return wave * (1 + 8.336624212083e-5 + 2.408926869968e-2 / (1.301065924522e2 - s2)
                   + 1.599740894897e-4 / (3.892568793293e1 - s2))


def median_snr(flux, err):
    ok = np.isfinite(flux) & np.isfinite(err) & (err > 0)
    return float(np.median(flux[ok] / err[ok])) if ok.any() else np.nan


def write_spectrum(path, wave, flux, err, **meta):
    tab = Table(dict(wave_vac=wave, flux=flux, flux_err=err))
    tab['wave_vac'].unit = 'Angstrom'
    tab['flux'].unit = tab['flux_err'].unit = 'erg / (s cm2 Angstrom)'
    tab.meta.update(meta)
    tab.write(path, overwrite=True)


def desi_candidates(hosts):
    """DESI DR1 spectra near the hosts, as (host index, sparcl_id, redshift)."""
    from sparcl.client import SparclClient
    client = SparclClient()
    pad = 2 * MATCH_ARCSEC / 3600
    constraints = {
        'data_release': ['DESI-DR1'],
        'ra': [float(hosts['ra'].min() - pad), float(hosts['ra'].max() + pad)],
        'dec': [float(hosts['dec'].min() - pad), float(hosts['dec'].max() + pad)],
    }
    found = client.find(outfields=['sparcl_id', 'ra', 'dec', 'redshift', 'redshift_warning'],
                        constraints=constraints, limit=1000000)
    recs = [r for r in found.records if r['redshift_warning'] == 0]
    print('DESI DR1 spectra in the field (no redshift warning):', len(recs))
    if not recs:
        return client, []
    coords = SkyCoord([r['ra'] for r in recs] * u.deg, [r['dec'] for r in recs] * u.deg)
    ih, ir, _, _ = coords.search_around_sky(SkyCoord(hosts['ra'] * u.deg, hosts['dec'] * u.deg), MATCH_ARCSEC * u.arcsec)
    out = [(i, recs[j]['sparcl_id'], recs[j]['redshift']) for i, j in zip(ih, ir)
           if abs(recs[j]['redshift'] - hosts['redshift'][i]) < MAX_DZ]
    return client, out


def fetch_desi(client, sparcl_ids, outdir, batch=50, tries=5):
    """Download DESI DR1 spectra not yet in outdir, in batches, retrying server errors."""
    todo = [s for s in sparcl_ids if not os.path.exists('%s/desi_%s.fits' % (outdir, s))]
    for start in range(0, len(todo), batch):
        for attempt in range(tries):
            try:
                recs = client.retrieve(uuid_list=todo[start:start + batch], dataset_list=['DESI-DR1'],
                                       include=['sparcl_id', 'wavelength', 'flux', 'ivar', 'redshift']).records
                break
            except Exception as e:
                print('  SPARCL retrieve failed (%s), retrying' % type(e).__name__)
                time.sleep(10 * (attempt + 1))
        else:
            raise RuntimeError('SPARCL retrieve kept failing')
        for rec in recs:
            ivar = np.asarray(rec['ivar'], dtype=float)
            err = np.where(ivar > 0, 1 / np.sqrt(np.where(ivar > 0, ivar, 1)), np.inf) * 1e-17
            write_spectrum('%s/desi_%s.fits' % (outdir, rec['sparcl_id']), np.asarray(rec['wavelength'], dtype=float),
                           np.asarray(rec['flux'], dtype=float) * 1e-17, err,
                           SURVEY='DESI-DR1', SPECID=str(rec['sparcl_id']), SPEC_RES=3000, Z=rec['redshift'])


def legac_candidates(hosts):
    """LEGA-C spectra near the hosts, as (host index, ESO dataset id, R)."""
    query = ("SELECT dp_id, s_ra, s_dec, em_res_power FROM ivoa.ObsCore "
             "WHERE obs_collection='LEGA-C' AND dataproduct_type='spectrum'")
    params = urllib.parse.urlencode(dict(REQUEST='doQuery', LANG='ADQL', FORMAT='csv', QUERY=query)).encode()
    with urllib.request.urlopen(ESO_TAP, params, timeout=300) as response:
        tab = Table.read(response.read().decode(), format='ascii.csv')
    print('LEGA-C spectra in the ESO archive:', len(tab))
    coords = SkyCoord(np.asarray(tab['s_ra']) * u.deg, np.asarray(tab['s_dec']) * u.deg)
    ih, ir, _, _ = coords.search_around_sky(SkyCoord(hosts['ra'] * u.deg, hosts['dec'] * u.deg), MATCH_ARCSEC * u.arcsec)
    return [(i, str(tab['dp_id'][j]), float(tab['em_res_power'][j])) for i, j in zip(ih, ir)]


def fetch_legac(dp_id, resolving_power, path):
    if not os.path.exists(path):
        with urllib.request.urlopen(ESO_FILE + dp_id, timeout=300) as response:
            data = fits.open(io.BytesIO(response.read()))[1].data
        wave = np.ravel(data['WAVE']).astype(float)
        flux = np.ravel(data['FLUX']).astype(float) * 1e-19
        err = np.ravel(data['ERR']).astype(float) * 1e-19
        bad = (np.ravel(data['QUAL']) != 0) | ~(err > 0)
        err[bad] = np.inf
        write_spectrum(path, air_to_vacuum(wave), flux, err, SURVEY='LEGA-C', SPECID=dp_id,
                       SPEC_RES=resolving_power)
    return Table.read(path)


data = Table.read(sys.argv[1])
_, first = np.unique(np.asarray(data['ID_COSMOS']), return_index=True)
hosts = Table(dict(ID_COSMOS=data['ID_COSMOS'][first], ra=data['ALPHA_J2000_GAL'][first],
                   dec=data['DELTA_J2000_GAL'][first], redshift=data['redshift'][first]))
print('chimera hosts:', len(hosts))
os.makedirs(OUTDIR, exist_ok=True)

client, desi = desi_candidates(hosts)
legac = legac_candidates(hosts)
print('candidate spectra: DESI DR1 %d, LEGA-C %d' % (len(desi), len(legac)))

best = {}
fetch_desi(client, sorted({specid for _, specid, _ in desi}), OUTDIR)
for i, specid, _ in desi:
    path = '%s/desi_%s.fits' % (OUTDIR, specid)
    spec = Table.read(path)
    best.setdefault(i, []).append(('DESI-DR1', path, median_snr(spec['flux'], spec['flux_err']), spec.meta['SPEC_RES']))
for i, dp_id, resolving_power in legac:
    path = '%s/legac_%s.fits' % (OUTDIR, dp_id.replace(':', '_'))
    spec = fetch_legac(dp_id, resolving_power, path)
    best.setdefault(i, []).append(('LEGA-C', path, median_snr(spec['flux'], spec['flux_err']), spec.meta['SPEC_RES']))

rows = []
for i, cands in sorted(best.items()):
    survey, path, snr, resolving_power = max(cands, key=lambda c: np.nan_to_num(c[2], nan=-1))
    rows.append((hosts['ID_COSMOS'][i], hosts['redshift'][i], survey, path, snr, resolving_power, len(cands)))
out = Table(rows=rows, names=['ID_COSMOS', 'redshift', 'survey', 'spectrum_file', 'snr_median', 'spec_res', 'n_candidates'])
print('hosts with a spectrum:', len(out), ' by survey:', dict(zip(*np.unique(out['survey'], return_counts=True))))
out.write(sys.argv[2], overwrite=True)
