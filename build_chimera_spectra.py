#!/usr/bin/env python
# coding: utf-8
"""Build chimera spectra: host galaxy spectrum + weighted quasar spectrum.

For each chimera whose host has a spectrum (galaxy-spectra.fits) and whose quasar
has an SDSS spectrum (qso-spectra.fits):

1. Each component is corrected for Milky Way extinction like its photometry
   (O'Donnell 1994 law, R_V = 3.1, via the extinction package): the host with the COSMOS E(B-V) used by
   gal/convert_gal.py, the quasar with the SFD E(B-V) behind the Shen et al. (2011)
   dereddened magnitudes, recovered as (rmag - rmag0) / 2.751. Then each component is
   scaled to its own photometry, so the spectrum agrees with the chimera photometry: the host to its converted SDSS photometry (*_sdss_GAL), the
   quasar to its SDSS photometry (*_sdss_QSO). The scale is the median ratio of
   photometry to synthetic top-hat photometry over the u/r/i/z bands the spectrum
   fully covers.
2. The quasar is moved to the chimera (host) redshift: wavelengths scale by
   (1+z)/(1+z_QSO), f_lambda by (1+z_QSO)/(1+z). Quasars are paired with hosts
   within |dz| < 0.01, so this changes fluxes by about 1%.
3. Both are convolved to a common resolving power R (Gaussian line-spread function,
   CHIMERA_SPECTRA_R, default 1500) and resampled onto a log-wavelength grid with
   2.5 pixels per resolution element, over the wavelength range covered by both.
4. composite = host + chimera_QSO_weight * quasar, with errors added in quadrature.

Each chimera is written to chimeras-spectra/<id>.fits with columns wave_obs (vacuum,
A), flux, flux_err, flux_gal and flux_qso (the weighted quasar), all in mJy; the
output table lists the chimeras with their provenance.

usage: build_chimera_spectra.py chimeras-fullinfo.fits gal/test_mass_onlygal_extragal.fits galaxy-spectra.fits qso-spectra.fits chimeras-spectra.fits
"""
import os
import sys

import extinction
import numpy as np
from astropy.table import Table

TARGET_R = float(os.environ.get('CHIMERA_SPECTRA_R', '1500'))
PIXELS_PER_RESEL = 2.5
OUTDIR = 'chimeras-spectra'
C_AA = 2.99792458e18
# SFD A_r / E(B-V) for SDSS r, used by Shen et al. (2011) for the dereddened magnitudes
SDSS_R_EXTINCTION = 2.751
trapezoid = getattr(np, 'trapezoid', None) or np.trapz  # numpy < 2 has only trapz
# top-hat approximations of the SDSS bands used by the chimera photometry [A]
BANDS = {'u_sdss': (3300, 3800), 'r_sdss': (5600, 6850), 'i_sdss': (6950, 8300), 'z_sdss': (8350, 9600)}


def deredden(wave, flam, err, ebv, rv=3.1):
    """Correct for Milky Way extinction with the O'Donnell (1994) law."""
    corr = 10**(0.4 * extinction.odonnell94(np.asarray(wave, dtype=float), rv * ebv, rv))
    return flam * corr, err * corr


def to_mjy(wave, flam):
    return flam * wave**2 / C_AA * 1e26


def synthetic_mjy(wave, flam, err, lo, hi):
    """Mean f_nu [mJy] in a top-hat band, or nan if the band is not fully covered."""
    ok = np.isfinite(flam) & np.isfinite(err)
    if not ok.any() or wave[ok].min() > lo or wave[ok].max() < hi:
        return np.nan
    sel = ok & (wave >= lo) & (wave <= hi)
    # photon-weighted mean f_nu of a top-hat response: int f_lambda lambda dlambda / int c/lambda dlambda
    return trapezoid(flam[sel] * wave[sel], wave[sel]) / trapezoid(C_AA / wave[sel], wave[sel]) * 1e26


def photometric_scale(wave, flam, err, row, suffix):
    """Median photometry / synthetic photometry over the covered bands, and the bands used."""
    ratios, used = [], []
    for band, (lo, hi) in BANDS.items():
        phot = float(row[band + suffix])
        synth = synthetic_mjy(wave, flam, err, lo, hi)
        if np.isfinite(phot) and phot > 0 and np.isfinite(synth) and synth > 0:
            ratios.append(phot / synth)
            used.append(band.split('_')[0])
    return (float(np.median(ratios)) if ratios else np.nan), ''.join(used)


def convolve_to_r(wave, flam, err, native_r):
    """Inverse-variance weighted Gaussian smoothing from native_r to TARGET_R on a log-wavelength grid."""
    good = np.isfinite(flam) & np.isfinite(err) & (err > 0)
    lnw = np.log(wave)
    step = np.median(np.diff(lnw))
    grid = np.arange(lnw[0], lnw[-1], step)
    w = np.interp(grid, lnw, np.where(good, 1 / np.where(good, err, 1)**2, 0.0))
    w[np.interp(grid, lnw, good.astype(float)) < 1] = 0.0  # pixels next to bad ones carry no weight
    f = np.interp(grid, lnw, np.where(good, flam, 0.0))
    sigma = np.sqrt(max(TARGET_R**-2 - native_r**-2, 0.0)) / 2.3548 / step
    if sigma > 0.3:
        x = np.arange(-int(np.ceil(4 * sigma)), int(np.ceil(4 * sigma)) + 1)
        kernel = np.exp(-0.5 * (x / sigma)**2)
        num, den = np.convolve(w * f, kernel, 'same'), np.convolve(w, kernel, 'same')
        var = np.convolve(w, kernel**2, 'same') / np.where(den > 0, den, 1)**2
    else:
        num, den, var = w * f, w, 1 / np.where(w > 0, w, 1)
    ok = den > 0.5 * np.max(den) * 1e-6
    return np.exp(grid), np.where(ok, num / np.where(ok, den, 1), np.nan), np.where(ok, np.sqrt(var), np.inf)


data = Table.read(sys.argv[1])
cosmos = Table.read(sys.argv[2])
galaxies = Table.read(sys.argv[3])
qsos = Table.read(sys.argv[4])
ebv_by_id = dict(zip(np.asarray(cosmos['NUMBER']).astype(int), np.asarray(cosmos['EBV'], dtype=float)))
gal_by_id = {int(g['ID_COSMOS']): g for g in galaxies}
qso_by_id = {str(q['ID_SDSS']).strip(): q for q in qsos}
data = data[np.isin(data['ID_COSMOS'], galaxies['ID_COSMOS'])]
print('chimeras whose host has a spectrum:', len(data))
os.makedirs(OUTDIR, exist_ok=True)

cache = {}
def load(path):
    if path not in cache:
        cache[path] = Table.read(path)
    return cache[path]

rows = []
for row in data:
    g, q = gal_by_id[int(row['ID_COSMOS'])], qso_by_id.get(str(row['ID_SDSS']).strip())
    if q is None:
        continue
    gspec, qspec = load(g['spectrum_file']), load(q['spectrum_file'])
    gw, gf, ge = (np.asarray(gspec[c], dtype=float) for c in ('wave_vac', 'flux', 'flux_err'))
    qw, qf, qe = (np.asarray(qspec[c], dtype=float) for c in ('wave_vac', 'flux', 'flux_err'))
    ebv_gal = ebv_by_id[int(row['ID_COSMOS'])]
    ebv_qso = (float(row['rmag_QSO']) - float(row['rmag0_QSO'])) / SDSS_R_EXTINCTION
    gf, ge = deredden(gw, gf, ge, ebv_gal)
    qf, qe = deredden(qw, qf, qe, ebv_qso)
    gscale, gbands = photometric_scale(gw, gf, ge, row, '_GAL')
    qscale, qbands = photometric_scale(qw, qf, qe, row, '_QSO')
    if not (np.isfinite(gscale) and np.isfinite(qscale)):
        continue
    z, zq, weight = float(row['redshift']), float(q['redshift_QSO']), float(row['chimera_QSO_weight'])
    wave_scale, flux_scale = (1 + z) / (1 + zq), (1 + zq) / (1 + z)
    gw, gf, ge = convolve_to_r(gw, gf * gscale, ge * gscale, float(g['spec_res']))
    qw, qf, qe = convolve_to_r(qw * wave_scale, qf * qscale * flux_scale, qe * qscale * flux_scale,
                               float(qspec.meta['SPEC_RES']))
    lo, hi = max(gw[0], qw[0]), min(gw[-1], qw[-1])
    if hi <= lo:
        continue
    wave = np.exp(np.arange(np.log(lo), np.log(hi), 1 / (PIXELS_PER_RESEL * TARGET_R)))
    gal = to_mjy(wave, np.interp(wave, gw, gf))
    qso = to_mjy(wave, weight * np.interp(wave, qw, qf))
    err = to_mjy(wave, np.hypot(np.interp(wave, gw, ge), weight * np.interp(wave, qw, qe)))
    err[~np.isfinite(gal + qso)] = np.inf
    out = Table(dict(wave_obs=wave, flux=gal + qso, flux_err=err, flux_gal=gal, flux_qso=qso))
    out['wave_obs'].unit = 'Angstrom'
    for c in ('flux', 'flux_err', 'flux_gal', 'flux_qso'):
        out[c].unit = 'mJy'
    out.meta.update(ID=str(row['id']), REDSHIFT=z, Z_QSO=zq, QSO_WT=weight, SPEC_RES=TARGET_R,
                    GAL_SURV=str(g['survey']), GAL_FILE=str(g['spectrum_file']), QSO_FILE=str(q['spectrum_file']),
                    EBV_GAL=ebv_gal, EBV_QSO=ebv_qso, GAL_SCAL=gscale, QSO_SCAL=qscale, GAL_BAND=gbands, QSO_BAND=qbands,
                    QSO_WSCL=wave_scale, QSO_FSCL=flux_scale)
    path = '%s/%s.fits' % (OUTDIR, row['id'])
    out.write(path, overwrite=True)
    rows.append((row['id'], z, row['ID_COSMOS'], row['ID_SDSS'], weight, str(g['survey']), path, len(wave),
                 wave[0], wave[-1], gscale, qscale, gbands, qbands, float(g['snr_median'])))

out = Table(rows=rows, names=['id', 'redshift', 'ID_COSMOS', 'ID_SDSS', 'chimera_QSO_weight', 'galaxy_survey',
                               'spectrum_file', 'n_pixels', 'wave_min', 'wave_max', 'galaxy_scale', 'qso_scale',
                               'galaxy_bands', 'qso_bands', 'galaxy_snr_median'])
out.meta['SPEC_RES'] = TARGET_R
print('chimera spectra written:', len(out), ' resolving power: %.0f' % TARGET_R)
if len(out):
    print('  host scale (photometry / spectrum) median %.2f, 16-84%% %.2f-%.2f' % tuple(np.percentile(out['galaxy_scale'], [50, 16, 84])))
    print('  quasar scale median %.2f, 16-84%% %.2f-%.2f' % tuple(np.percentile(out['qso_scale'], [50, 16, 84])))
out.write(sys.argv[5], overwrite=True)
