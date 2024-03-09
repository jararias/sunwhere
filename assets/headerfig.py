
import pylab as pl
import pandas as pd

import sunwhere


times = pd.date_range('2024-07-01', '2024-07-02', freq='1min')
site_lat, site_lon = 36.949, -3.822

sw = sunwhere.sites(times, site_lat, site_lon)
df = (sw.elevation.to_dataframe('elev')
      .reset_index().get(['time', 'elev'])
      .assign(
          tst=sw.true_solar_time.to_numpy(),
          sza=sw.zenith.isel(location=0).to_numpy(),
          saa=sw.azimuth.isel(location=0).to_numpy()))

fig, ax = pl.subplots(1, 1, figsize=(16, 3.5), layout='constrained')

ax.scatter(x='saa', y='elev', c='elev', cmap='autumn', s=2500,
           data=df.where(df.elev >= 0, float('nan')), clip_on=False)

ax.set(ylim=(-10, 85))
ax.tick_params(labelleft=False, labelbottom=False, left=False, bottom=False)

ax.text(0.5, 0.35, 'sunwhere', ha='center', va='center',
        transform=ax.transAxes, fontsize=140,
        family='', weight='medium')

for spine in ax.spines.values():
    spine.set_visible(False)

pl.savefig('headerfig.png')
