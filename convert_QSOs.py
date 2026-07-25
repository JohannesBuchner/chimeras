#!/usr/bin/env python
# coding: utf-8
import numpy as np
import sys
from astropy.table import Table

MIN_LOG_L5100 = 45.0

def ab_flux(ABmag):
    """convert AB magnitude to flux in mJy"""
    return (10**(23-(((ABmag)+48.6)/2.5))) * 1000

def magerr2fluxerr(mag, magerr):
    """convert the lower (brighter) magnitude to flux, return the deviation in flux."""
    return ab_flux(mag - magerr) - ab_flux(mag)

def convert_flux_and_err_where_measured(mag, magerr):
    """convert safely where error is positive, return nan otherwise."""
    flux = np.where(magerr > 0, ab_flux(mag), np.nan)
    flux_err = np.where(magerr > 0, magerr2fluxerr(mag, magerr), np.nan)
    return flux, flux_err

dataQ = Table.read(sys.argv[1])
print("Quasar input catalog:", len(dataQ))
# convert magnitudes to fluxes, including errors
dataQ["u_sdss"], dataQ["u_sdss_err"] = convert_flux_and_err_where_measured(dataQ['umag0'], dataQ['e_umag'])
dataQ["r_sdss"], dataQ["r_sdss_err"] = convert_flux_and_err_where_measured(dataQ['rmag0'], dataQ['e_rmag'])
dataQ["i_sdss"], dataQ["i_sdss_err"] = convert_flux_and_err_where_measured(dataQ['imag0'], dataQ['e_imag'])
dataQ["z_sdss"], dataQ["z_sdss_err"] = convert_flux_and_err_where_measured(dataQ['zmag0'], dataQ['e_zmag'])

dataQ["J_2mass"], dataQ["J_2mass_err"] = convert_flux_and_err_where_measured(dataQ['Jmag']+0.91, dataQ['e_Jmag'])
dataQ["H_2mass"], dataQ["H_2mass_err"] = convert_flux_and_err_where_measured(dataQ['Hmag']+1.39, dataQ['e_Hmag'])
dataQ["Ks_2mass"], dataQ["Ks_2mass_err"] = convert_flux_and_err_where_measured(dataQ['Ksmag']+1.85, dataQ['e_Ksmag'])
# convert with (very small) Assef conversions of quasars from WISE to IRAC
dataQ["IRAC1"], dataQ["IRAC1_err"] = convert_flux_and_err_where_measured(dataQ['W1mag'] + 2.699 - 0.01, dataQ['e_W1mag'])
dataQ["IRAC2"], dataQ["IRAC2_err"] = convert_flux_and_err_where_measured(dataQ['W2mag'] + 3.339 - 0.07, dataQ['e_W2mag'])

# keep only those with some data in NIR and MIR
has_NIR = np.logical_or(dataQ["J_2mass_err"]>0, np.logical_or(dataQ["H_2mass_err"]>0, dataQ["Ks_2mass_err"]>0))
has_MIR = np.logical_or(dataQ["IRAC1_err"]>0, dataQ["IRAC2_err"]>0)
log_l5100 = np.asarray(np.ma.filled(dataQ["logL5100"], np.nan), dtype=float)
# Keep unmeasured sources; apply the lower limit only where L5100 is available.
passes_l5100_cut = ~np.isfinite(log_l5100) | (log_l5100 > MIN_LOG_L5100)
print("Quasar with NIR:", has_NIR.sum())
print("Quasar with MIR:", has_MIR.sum())
print("Quasar with NIR+MIR:", np.logical_and(has_NIR, has_MIR).sum())
print(f"Quasar with missing logL5100 or logL5100 > {MIN_LOG_L5100:g}:", passes_l5100_cut.sum())

selection = has_NIR & has_MIR & passes_l5100_cut
print("Quasar with NIR+MIR and luminosity cut:", selection.sum())
dataQ[selection].write(sys.argv[2], overwrite=True)
