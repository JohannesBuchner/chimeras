from getdist import plots, MCSamples
import matplotlib.pyplot as plt
import numpy as np
import astropy.table
import corner
import sys
import tqdm

#inputfilename = 'input.fits'
#outputfilename = 'input-unique.fits'
inputfilename = sys.argv[1]
outputfilename = sys.argv[2]

data = astropy.table.Table.read(inputfilename)
for c in data.colnames:
    if '_err' in c:
        del data[c]
ids = data['id']
columns = ['redshift'] + [c for c in data.colnames[data.colnames.index('u_sdss'):-1] if '_err' not in c]
print("columns:", columns)
print("rows:", len(data))

all_fluxes = np.empty((len(data), len(columns)))
for i, c in enumerate(columns):
    all_fluxes[:,i] = np.log10(np.abs(np.array(data[c])))
    print(i, c, np.isfinite(all_fluxes[:,i]).mean())
for c in data.colnames:
    if c != 'id':
        del data[c]

print("plotting absolute flux grid...")
#all_fluxes[~np.isfinite(all_fluxes[:,1]), 1] = -6
#all_fluxes[~np.isfinite(all_fluxes)] = -5
corner.corner(np.nan_to_num(all_fluxes, nan=-5), labels=columns);
plt.savefig('resample-all.pdf')
plt.close()


print("plotting relative flux grid...")
rel_fluxes = all_fluxes.copy()
# normalise relative to r-band
rel_fluxes[:,1:] -= rel_fluxes[:,2:3]
rel_fluxes[:,2] = all_fluxes[:,2]

corner.corner(np.nan_to_num(rel_fluxes, nan=-4), labels=columns);
plt.savefig('resample-rel.pdf')
plt.close()

nrows, ncols = rel_fluxes.shape

print("finding indices...")
lo = np.nanmin(rel_fluxes, axis=0)
hi = np.nanmax(rel_fluxes, axis=0)

# for each column, compute index from 0 to 9 (last bin 9 is for invalid/unknown)
col_indices = np.nan_to_num(8.99999 * (rel_fluxes - lo.reshape((1, -1))) / (hi - lo).reshape((1, -1)), nan=9.).astype(int)

# for each row, sum up indices
row_indices = (col_indices * (10**np.arange(ncols)).reshape((1, -1))).sum(axis=1)
del rel_fluxes

_, unique_selection, unique_weight = np.unique(row_indices, return_index=True, return_counts=True)
print("unique:", len(unique_selection))
print("plotting selected distribution...")
settings = dict(names = columns, labels = columns, sampler='nested', 
                settings={'smooth_scale_2D':0.3, 'smooth_scale_1D':0.3})

samples  = MCSamples(samples=np.nan_to_num(all_fluxes, nan=-5), label='All',  **settings)
samples3 = MCSamples(samples=np.nan_to_num(all_fluxes[unique_selection,:], nan=-5), label='Unique', **settings)

g = plots.get_subplot_plotter()
g.settings.num_plot_contours = 5
g.triangle_plot([samples, samples3], filled=False)
del all_fluxes, samples, samples3
plt.savefig('resample-result.pdf')
plt.close()

print("writing out selection...")
unique_table = data[unique_selection]
del data
unique_table['resample_weight'] = unique_weight
unique_table.write(outputfilename, overwrite=True)
