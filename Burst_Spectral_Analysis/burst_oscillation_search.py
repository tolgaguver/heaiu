# %%
%matplotlib inline

# load auxiliary libraries
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

# import stingray
import stingray
# choose style of plots, `seaborn-talk` produce nice big figures
plt.style.use('seaborn-talk')
# %%
f = fits.open('/Users/tolga/ownCloud/burst_characterization_v4/4U_1608-522/burst1/ni0050070102_0mpu7_cl_bary_fpm_osc_b1.evt')
dt = f[1].header['TIMEDEL']
toa = f[1].data['Time']
dt = 0.00012207031
lc = stingray.Lightcurve.make_lightcurve(toa=toa, dt=dt)
# %%
# Create the lightcurve object
lc.plot(labels=['Time (s)', 'Counts / bin'])
plt.title('Lightcurve')
# %%
dynspec = stingray.DynamicalPowerspectrum(lc=lc, segment_size=1.0, norm='leahy')# %%
dynspec.dyn_ps
extent = min(dynspec.time), max(dynspec.time), max(dynspec.freq), min(dynspec.freq)
plt.imshow(dynspec.dyn_ps, origin="lower", aspect="auto", vmin=1.98, vmax=3.0,
           interpolation="none", extent=extent)
plt.colorbar()
plt.ylim(618,622)
plt.xlim(min(dynspec.time)+15,min(dynspec.time)+55)
# %%
