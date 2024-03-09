
import locale
from datetime import datetime

from tabulate import tabulate

import numpy as np
import pandas as pd

import sunwhere


locale.setlocale(locale.LC_TIME, locale.getlocale())


def calculate_position(time, lat, lon, algorithm, refraction, timezone):
    # TODO: decorate input & output!!

    if time == 'now':
        time = datetime.now()

    times = np.array(
        pd.to_datetime(time).tz_localize(None), ndmin=1, dtype='datetime64[ns]')
    times = pd.to_datetime(times).tz_localize(timezone)

    sw = sunwhere.sites(times, lat, lon, algorithm=algorithm, refraction=refraction)

    sr_local = pd.Timestamp(sw.sunrise(units="local").isel(time=0, location=0).item())
    ss_local = pd.Timestamp(sw.sunset(units="local").isel(time=0, location=0).item())

    s = pd.Series(
        {
            'Latitude': f'{lat:.4f}\u00b0',
            'Longitude': f'{lon:.4f}\u00b0',
            'Zenith': f'{sw.sza.isel(time=0, location=0):.4f}\u00b0',
            'Azimuth': f'{sw.saa.isel(time=0, location=0):.4f}\u00b0',
            'Declination': f'{sw.dec.isel(time=0):.4f}\u00b0',
            'Sunrise': f'{sr_local.tz_localize(sw.timezone):%H:%M:%S(%Z)}',
            'Sunset': f'{ss_local.tz_localize(sw.timezone):%H:%M:%S(%Z)}'
        }
    ).to_frame('ephemeris')

    ts = pd.to_datetime(time).tz_localize(timezone)
    print(' '.join(map(str.capitalize, ts.strftime('%A, %d %B %Y %H:%M:%S').split()))
          + ts.strftime('(%Z):'))

    print(tabulate(s, tablefmt='double_outline', colalign=['left', 'right']))
