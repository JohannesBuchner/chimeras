#!/usr/bin/env python
# coding: utf-8

# Donley AGN selection criterion based on infrared colors from IRAC

import numpy as np
import astropy, astropy.table
import astropy.stats
import sys

galc = astropy.table.Table.read(sys.argv[1])

x = np.log10(galc['IRAC3']/galc['IRAC1'])
y = np.log10(galc['IRAC4']/galc['IRAC2'])

# for magnitudes
#x = (IRAC3-IRAC1)/(-2.5)
#y = (IRAC4-IRAC2)/(-2.5)

xy = np.logical_and(x>=0.08, y>=0.15)
yy = np.logical_and(y>=(1.21*x)-0.27, y<=(1.21*x)+0.27)
num = np.logical_and(x, y)

r1 = np.logical_and(galc['IRAC1']<galc['IRAC2'], galc['IRAC2']<galc['IRAC3'])
ref = np.logical_and(r1, galc['IRAC3']<galc['IRAC4'])

agn = np.logical_and(ref, num)

np.savetxt(sys.argv[2], galc[agn]['id'], fmt='%s')
