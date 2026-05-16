
import locale
from datetime import datetime, UTC

from tabulate import tabulate

import pandas as pd

import sunwhere


locale.setlocale(locale.LC_TIME, locale.getlocale())


def calculate_position(time, lat, lon, algorithm, refraction):

    # Online service to verify the results of NREL's SPA
    # https://www.kso.ac.at/beobachtungen/ephem_api.php?date=20230715&time=12%3A00%3A00&lat=40.0000&lon=0.0000

    if time == 'now':
        time_local = datetime.now().astimezone()  # local hour
    else:
        time_local = pd.to_datetime(time).tz_localize(None).astimezone()  # local hour
    time_utc = datetime.now(UTC)  # UTC hour

    sw = sunwhere.sites(time_local, lat, lon, algorithm=algorithm, refraction=refraction)
    sr_utc = pd.Timestamp(sw.sunrise(units="utc").isel(time=0, site=0).item())
    ss_utc = pd.Timestamp(sw.sunset(units="utc").isel(time=0, site=0).item())
    sr_local = pd.Timestamp(sw.sunrise(units="local").isel(time=0, site=0).item())
    ss_local = pd.Timestamp(sw.sunset(units="local").isel(time=0, site=0).item())
    eot = sw.eot.isel(time=0).item()
    eot_min = int(eot)
    eot_sec = 60*(abs(eot) % 1)

    ts_local = pd.to_datetime(time_local)
    ts_utc = pd.to_datetime(time_utc)

    s = pd.Series(
        {
            'Day of week': ts_local.strftime("%A").capitalize(),
            'Date': ' '.join(map(str.capitalize, ts_local.strftime("%d %B %Y").split())),
            'Time': ts_local.strftime("%H:%M:%S(%Z)"),
            'Time (UTC)': ts_utc.strftime("%H:%M:%S"),
            'Latitude': f'{lat:.4f}\u00b0',
            'Longitude': f'{lon:.4f}\u00b0',
            'Declination': f'{sw.dec.isel(time=0):.4f}\u00b0',
            'Sun-Earth distance': f'{1/(sw.ecf.isel(time=0)**0.5):.4f} AU',
            'Zenith': f'{sw.sza.isel(time=0, site=0):.4f}\u00b0',
            'Elevation': f'{sw.elevation.isel(time=0, site=0):.4f}\u00b0',
            'Azimuth (0\u00b0 south)': f'{sw.saa.isel(time=0, site=0):.4f}\u00b0',
            'Sunrise(UTC)': f'{sr_utc.tz_localize(sw.timezone):%H:%M:%S}',
            'Sunset(UTC)': f'{ss_utc.tz_localize(sw.timezone):%H:%M:%S}',
            f'Sunrise({time_local.tzname()})': f'{sr_local.tz_localize(sw.timezone):%H:%M:%S}',
            f'Sunset({time_local.tzname()})': f'{ss_local.tz_localize(sw.timezone):%H:%M:%S}',
            'Equation of time': f'{eot_min}:{eot_sec:.1f} min'
        }
    ).to_frame('ephemeris')

    print(tabulate(s, tablefmt='double_outline', colalign=['left', 'right']))
