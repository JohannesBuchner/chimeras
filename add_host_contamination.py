#!/usr/bin/env python
# coding: utf-8
"""Estimate how much quasar host light each chimera carries.

The quasar photometry includes the quasar's own host galaxy, which is added to
the COSMOS galaxy with the same weight as the quasar. Following Shen et al. (2011),
the quasar host contributes a fraction host_frac_5100 of the (uncorrected) L5100
(computed in convert_QSOs.py). This script compares the weighted quasar host light
at rest-frame 5100A with the COSMOS galaxy's own rest-frame 5100A luminosity and
adds:

  logL5100_GAL                log10 lambda L_lambda(5100A) of the COSMOS galaxy [erg/s],
                              interpolated from its photometry
  QSO_host_to_GAL_5100        weight * host_frac_5100 * L5100_QSO / L5100_GAL
  QSO_host_contamination_dex  log10(1 + QSO_host_to_GAL_5100): approximate bias in the
                              galaxy stellar mass (and luminosity) if the quasar host is
                              attributed to the COSMOS galaxy, assuming similar M/L

Galaxy-property tests can cut on QSO_host_contamination_dex (e.g. < 0.05);
AGN-property tests can use all chimeras.
"""
import sys

import numpy as np
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table

# Same cosmology as Shen et al. (2011) luminosities.
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
# Galaxy bands (mJy) and approximate effective wavelengths [A].
BANDS = {
    'u_sdss': 3557, 'r_sdss': 6166, 'i_sdss': 7480, 'z_sdss': 8932,
    'J_2mass': 12350, 'H_2mass': 16620, 'Ks_2mass': 21590,
    'IRAC1': 35500, 'IRAC2': 44930,
}


def galaxy_log_l5100(data):
    """log lambda L_lambda at rest 5100A by log-log interpolation of the galaxy photometry."""
    wave = np.array(list(BANDS.values()), dtype=float)
    flux = np.vstack([np.asarray(data[band + '_GAL'], dtype=float) for band in BANDS]).T
    z = np.asarray(data['redshift'], dtype=float)
    target = np.log10(5100.0 * (1 + z))
    out = np.full(len(data), np.nan)
    for i in range(len(data)):
        ok = np.isfinite(flux[i]) & (flux[i] > 0)
        if ok.sum() < 2:
            continue
        lw, lf = np.log10(wave[ok]), np.log10(flux[i, ok])
        if not lw[0] <= target[i] <= lw[-1]:
            continue
        out[i] = np.interp(target[i], lw, lf)
    # nu L_nu (= lambda L_lambda) at rest 5100A from F_nu [mJy] at observed 5100(1+z)A.
    nu_obs = 2.99792458e18 / (5100.0 * (1 + z))
    d_l = cosmo.luminosity_distance(z).to('cm').value
    return out + np.log10(1e-26 * nu_obs * 4 * np.pi * d_l**2)


data = Table.read(sys.argv[1])
data['logL5100_GAL'] = galaxy_log_l5100(data)
qso_host = (np.asarray(data['chimera_QSO_weight'], dtype=float)
            * np.asarray(data['host_frac_5100_QSO'], dtype=float)
            * 10**np.asarray(data['logL5100_est_QSO'], dtype=float))
data['QSO_host_to_GAL_5100'] = qso_host / 10**data['logL5100_GAL']
data['QSO_host_contamination_dex'] = np.log10(1 + data['QSO_host_to_GAL_5100'])

contam = np.asarray(data['QSO_host_contamination_dex'])
print('chimeras:', len(data), ' with contamination estimate:', np.isfinite(contam).sum())
for limit in 0.01, 0.05, 0.1, 0.3:
    print('  QSO host contamination > %.2f dex: %d' % (limit, np.sum(contam > limit)))
data.write(sys.argv[2], overwrite=True)
