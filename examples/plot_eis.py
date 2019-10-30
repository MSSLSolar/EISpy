"""
Loading and plotting EIS Level 0 files
======================================
"""

import os
import eispy.cube
import matplotlib.pyplot as plt
import matplotlib.colors as mcolor
import numpy as np

fname = 'eis_l0_20100723_143210.fits'
remote_dir = 'http://solar.ads.rl.ac.uk/MSSL-data/eis/level0/2010/07/23'

if not os.path.exists(fname) and not os.path.exists(f'{fname}.gz'):
    import urllib.request
    dl = f'{remote_dir}/{fname}.gz'
    urllib.request.urlretrieve(dl, f'{fname}.gz')

if not os.path.exists(fname):
    import gzip
    with gzip.open(f'{fname}.gz', 'rb') as f:
        with open(fname, 'wb') as g:
            g.write(f.read())


observation = eispy.cube.read(fname)
print(observation)
print(observation.wavelengths)

wlens = ['FE XIII 203.830', 'FE XV 284.160',
         'FE XVI 262.980', 'CA XIV 194.100',
         'CA XVI 208.500']
fig = plt.figure(figsize=(14, 8))
for i, wlen in enumerate(wlens):
    intensity = observation[wlen].total_intensity
    ax = fig.add_subplot(1, len(wlens), i + 1, projection=intensity.wcs)
    intensity.plot(axes=ax, norm=mcolor.LogNorm())
    ax.set_title(wlen)

    xlim = ax.get_xlim()
    ax.set_xlim(xlim[1], xlim[0])
    cdelt = intensity.wcs.wcs.cdelt
    ax.set_aspect(np.abs(cdelt[1] / cdelt[0]))

    ax.set_xlabel('')
    ax.set_ylabel('')


plt.show()
