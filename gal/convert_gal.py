import numpy
import astropy.io.fits as pyfits
import tqdm
import sys
import os

"""
mag_AB=-2.58(log10(flux_muJ)+23.9
1muJ=10^-29 erg/s/cm^2/Hz

mag_B=-2.5log10(flux)+-23.9
and from microJ to the units needed, just multiply by 10^-29.

for IRac mag=18.1647
should be 196.84 microJ +- 0.9

returns muJy
"""
def ab_flux(ABmag):
    """convert AB magnitude to flux in uJy"""
    return (10**(23-(((ABmag)+48.6)/2.5))) * 1000000
def magerr2fluxerr(mag, magerr):
    """convert the lower (brighter) magnitude to flux, return the deviation in flux."""
    return ab_flux(mag - magerr) - ab_flux(mag)

print('reading input file...')
f = pyfits.open(sys.argv[1])
print('%d galaxies' % len(f[1].data))
ra = f[1].data['ALPHA_J2000']
dec = f[1].data['DELTA_J2000']
#mask = numpy.logical_and(
#	numpy.logical_and(ra > 149.5, ra < 150.7),
#	numpy.logical_and(dec > 1.65, dec < 2.7))

# filter based on SExtractor flags, X-ray detected, not a star, has redshift
# match each X-ray detection in file A
# with each object in file B

isolated_mask = f[1].data['Ks_FLAGS'] == 0
# defected photometry: dead pixels, nearby star, 

print('%d are clean of blending' % isolated_mask.sum())
neighbors = numpy.loadtxt('cosmos_acs_iphot_200709.coords')
#ra, dec = neighbors.transpose()
#mask = numpy.logical_and(
#	numpy.logical_and(ra > 149.5, ra < 150.7),
#	numpy.logical_and(dec > 1.65, dec < 2.7))
#neighbors = neighbors[mask,:]
nneighbors = []
for ra, dec in zip(tqdm.tqdm(f[1].data['ALPHA_J2000']), f[1].data['DELTA_J2000']):
	dist = ((ra - neighbors[:,0])**2 + (dec - neighbors[:,1])**2)**0.5
	#print dist.min()
	nneighbors.append((dist < 2. / 60. / 60).sum())

nneighbors = numpy.array(nneighbors)
print('within 2": 0:%d 1:%d 2:%d >2:%d' % ((nneighbors == 0).sum(), (nneighbors == 1).sum(), (nneighbors == 2).sum(), (nneighbors > 2).sum()))
isolated_mask = numpy.logical_and(isolated_mask, nneighbors <= 1)
print('%d are also without HST neighbors' % isolated_mask.sum())
if os.environ.get('UNBLENDED', '1') == '1':
	print('cutting others, so we have isolated sources only ...')
	f[1].data = f[1].data[isolated_mask]
else:
	print('cutting these away, so we have blended sources only ...')
	f[1].data = f[1].data[~isolated_mask]

print('%d: final sample' % len(f[1].data))

# look at flux columns (FLUX in name, APER3)
# ignore the ones with MAG, APER2
def errfield(c):
	c1 = c.replace('FLUX', 'FLUXERR')
	if c1 in allcolumns: 
		return c1
	c2 = c.replace('FLUX', 'FLUX_ERR')
	if c2 in allcolumns: 
		return c2
	assert False, (c, c1, c2, allcolumns)

ignore = ['FLUX_RADIUS', 
	'FLUXPEAK_90CM', 'FLUXPEAK_20CM', 'FLUXINT_90CM', 'FLUXINT_20CM', 
	'FLUX_XMM_0.2_2', 'FLUX_XMM_0.5_2', 'FLUX_XMM_2_10', 'FLUX_XMM_5_10', 'FLUX_XMM_0.5_10',
	'FLUX_CHANDRA_0.2_2', 'FLUX_CHANDRA_0.5_2', 'FLUX_CHANDRA_2_10', 'FLUX_CHANDRA_5_10', 'FLUX_CHANDRA_0.5_10',
	'FLUX_NUSTAR_3_24', 'FLUX_NUSTAR_3_8', 'FLUX_NUSTAR_8_24'
	]
allcolumns = [c.name.upper() for c in f[1].data.columns]
columns = [c for c in allcolumns if 'FLUX' in c and 
	'ERR' not in c and 'APER2' not in c and c not in ignore]

columns_errs = [errfield(c) for c in columns]
for c in columns_errs:
	 if c not in allcolumns:
	 	print(c, 'missing')

assert len(columns) == len(columns_errs), (len(columns), len(columns_errs))

filternames = []
#newcolumns = list(f[1].data.columns)
newcolumns = []

# rename columns:
for l in open('filternames.txt'):
	filtername, colname, EBVfactor = l.split()
	EBVfactor = float(EBVfactor)
	# filtername = filtername.replace('filters/', 'cosmos/')
	colname = colname.upper()
	guesses = [colname + '_FLUX_APER3', 'FLUX_' + colname, colname + '_FLUX_APER3', colname + '_FLUX']
	matches = [c for c in columns if (colname + '_') in c.replace('FLUX', '').replace('ERR', '').replace('APER3', '') + '_']
	for c in guesses + matches:
		if c.upper() in columns:
			break
	ce = errfield(c)

	#print 'processing', c, '-->', filtername
	i = allcolumns.index(c)
	bad = f[1].data[ce] < 0
	if colname + '_FLAGS' in allcolumns:
		if os.environ.get('KEEP_BLENDED_PHOTOMETRY', '0') != '1':
			blended = (f[1].data[colname + '_FLAGS'] != 0)
			bad = numpy.logical_or(bad, blended)
			if blended.any():
				print('   blending in %s for %d' % (colname, blended.sum()))
	
	flux = f[1].data[c]
	fluxerr = f[1].data[ce]
	
	if 'GALEX' in colname:
		# that value is actually a AB magnitude and we have to convert it
		# to flux
		mag = flux
		magerr = fluxerr
		flux = ab_flux(mag)
		fluxerr = magerr2fluxerr(mag, magerr)
	
	# mark zero flux as bad
	bad = numpy.logical_or(bad, flux == 0)
		
	EBV = f[1].data['EBV']
	# these are in uJy, convert to mJy
	# Apply E(B-V) from Milky Way 
	newvalues = flux * 0.001 * 10**(0.4 * EBVfactor * EBV)
	newvalues[bad] = -9999
	newcolumns.append(pyfits.Column(name=filtername,format='E',array=newvalues, unit='mJy'))
	#f[1].data[c] *= 0.001 * 10**(-0.4 * EBVfactor * EBV)
	#f[1].data[c][bad] = -9999
	#f[1].data.columns[i].unit = 'mJy'
	#f[1].data.columns[i].name = filtername
	filternames.append(filtername)
	#print 'processing', ce, '-->', filtername + '_err'
	i = allcolumns.index(ce)
	newvalues = fluxerr * 0.001 * 10**(0.4 * EBVfactor * EBV)
	newvalues[bad] = -9999
	newcolumns.append(pyfits.Column(name=filtername + '_err',format='E',array=newvalues, unit='mJy'))
	#f[1].data[ce] *= 0.001
	#f[1].data[ce][bad] = -9999
	#f[1].data.columns[i].unit = 'mJy'
	#f[1].data.columns[i].name = filtername + '_err'
	filternames.append(filtername + '_err')

info_columns = []
zcol = f[1].data['zs']

i = allcolumns.index('ID')
col = f[1].data.columns[i]
info_columns.append(pyfits.Column(name='id',format=col.format,array=col.array))
info_columns.append(pyfits.Column(name='redshift',format='E',array=zcol))

copy_columns = ['ALPHA_J2000', 'DELTA_J2000', 'MODEL', 'AGE', 'EXTINCTION', 'L_NU', 'L_R', 'L_K']
for param in 'MASS', 'SFR', 'SSFR':
    copy_columns += [param + '_BEST', param + '_MED', param + '_MED_MIN68', param + '_MED_MAX68']

for colname in copy_columns:
	i = allcolumns.index(colname)
	col = f[1].data.columns[i]
	info_columns.append(pyfits.Column(name=col.name,format=col.format,array=col.array))

print(', '.join(filternames))

tbhdu = pyfits.BinTableHDU.from_columns(pyfits.ColDefs(info_columns + newcolumns))
hdu = f[0]
import datetime, time
now = datetime.datetime.fromtimestamp(time.time())
nowstr = now.isoformat()
nowstr = nowstr[:nowstr.rfind('.')]
hdu.header['CREATOR'] = """Johannes Buchner <johannes.buchner.acad@gmx.com>"""
hdu.header['DATE'] = nowstr
hdulist = pyfits.HDUList([hdu, tbhdu])
print('writing to "%s"...' % sys.argv[2])
hdulist.writeto(sys.argv[2], overwrite=True)
print('done.')
