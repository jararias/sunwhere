
from pathlib import Path

import urllib3
import numpy as np
import pandas as pd
from loguru import logger
import dateutil.parser as dtparser


def get_jpl_ephemeris(out_filename, start_time, stop_time, site_lon,
                      site_lat, site_elev=0., step_size=1):
    # site_elev in km
    # step_size in minutes

    server = 'https://ssd.jpl.nasa.gov/api/horizons.api'
    format = 'format=text'
    _options = {
        'COMMAND': '10',                # the Sun
        'OBJ_DATA': 'YES',
        'MAKE_EPHEM': 'YES',
        'EPHEM_TYPE': 'OBSERVER',
        'CENTER': 'coord@399',
        # 'START_TIME': '2024-01-01',   # format: yyyy-mm-dd hh:mm
        # 'STOP_TIME': '2024-01-02',    # format: yyyy-mm-dd hh:mm
        # 'STEP_SIZE': '1 MINUTES',
        'QUANTITIES': '4,20',           # solar azimuth, elevation, and sun-earth distance (i.e., range)
        'COORD_TYPE': 'GEODETIC',       # to provide the location coordinates as lon and lat
        # 'SITE_COORD': '-3.5,37.5,0',  # location's lon, lat and altitude
        'REF_SYSTEM': 'ICRF',
        'CAL_FORMAT': 'BOTH',
        'CAL_TYPE': 'M',
        'TIME_DIGITS': 'SECONDS',
        'ANG_FORMAT': 'DEG',            # units for angles
        'APPARENT': 'AIRLESS',          # airless for NO atmospheric refraction
        'RANGE_UNITS': 'AU',            # units for the sun-earth distance
        'SUPPRESS_RANGE_RATE': 'YES',
        'SKIP_DAYLT': 'NO',
        'SOLAR_ELONG': '0,180',
        'EXTRA_PREC': 'NO',
        'R_T_S_ONLY': 'NO',
        'CSV_FORMAT': 'NO'
    }

    options = _options | {
        'START_TIME': pd.to_datetime(start_time).strftime('%Y-%m-%d %H:%M'),
        'STOP_TIME': pd.to_datetime(stop_time).strftime('%Y-%m-%d %H:%M'),
        'SITE_COORD': f'{site_lon:.6f},{site_lat:.6f},{site_elev:.3f}',
        'STEP_SIZE': f'{step_size:.0f} MINUTES'
    }

    max_lines = 90024  # max. allowed request lines !!!
    n_times = len(pd.date_range(options['START_TIME'], options['STOP_TIME'],
                                freq=f"{options['STEP_SIZE'].split()[0]}min"))
    if n_times > max_lines:
        raise ValueError(f'your request involves {n_times} data lines, which '
                         f'exceeds the max. allowed request lines of {max_lines}')

    url = server + '?' + format + '&' + '&'.join([f"{k}='{v}'" for k, v in options.items()])

    http = urllib3.PoolManager()
    response = http.request('GET', url)
    data = response.data.decode('utf-8')
    http.clear()

    path = Path(out_filename)
    if not path.parent.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
    path.open(mode='w').write(data)

    return path


def read_jpl_ephemeris_file(filename):

    def parse_line(line):
        dt64 = 'datetime64[ns]'
        dt = np.datetime64(dtparser.parse(line[1:21])).astype(dt64).astype('float64')
        jd = float(line[22:39])
        # sun_moon_presence = line[40:42]
        az = float(line[43:54])
        el = float(line[55:66])
        r = float(line[67:84])  # in AU
        return dt, jd, az, el, r

    with open(filename, 'r') as f:
        lines = f.readlines()
    first_data_line = lines.index('$$SOE\n')+1
    last_data_line = lines.index('$$EOE\n', first_data_line)
    a = np.array([parse_line(line) for line in lines[first_data_line: last_data_line]])
    return (pd.DataFrame(a)
            .set_index(pd.to_datetime(a[:, 0])).drop(columns=0)
            .set_axis(['julian_day', 'azimuth0-360', 'elevation', 'R'], axis=1)
            .rename_axis('time_utc', axis=0)
            .assign(
                azimuth=lambda df: df['azimuth0-360']-180,
                zenith=lambda df: 90.-df['elevation'],
                ecf=lambda df: 1/df['R']**2)
            .drop(columns=['julian_day', 'azimuth0-360', 'elevation', 'R'])
            .tz_localize('UTC'))


def load_jpl_ephemeris(year=2024, site_lon=-3.822, site_lat=36.949):
    out_dir = Path('data')
    if not out_dir.exists():
        out_dir.mkdir(parents=True, exist_ok=True)

    whole_data_filename = out_dir / f'jpl_ephemeris_highres_{year}.parquet'
    if not whole_data_filename.exists():
        data = []
        filenames = []
        for month in range(1, 13):
            out_filename = out_dir / f'jpl_ephemeris_highres_{year}-{month:02d}.txt'
            if not out_filename.exists():
                logger.info(f'downloading high-res ephemeris to file {out_filename}')
                start_time = pd.to_datetime(f'{year}-{month:02d}-01 00:00')
                stop_time = (start_time + pd.Timedelta('35days')).replace(day=1) - pd.Timedelta('1min')
                get_jpl_ephemeris(out_filename, start_time, stop_time, site_lon, site_lat, step_size=12)
                filenames.append(out_filename)
                logger.success('data downloaded')
            logger.info(f'loading high-res ephemeris from file {out_filename}')
            data.append(read_jpl_ephemeris_file(out_filename).drop(columns='julian_day'))
        df = pd.concat(data, axis=0)
        df.insert(0, 'site_lat', site_lat)
        df.insert(0, 'site_lon', site_lon)

        df.to_parquet(whole_data_filename)

        for filename in filenames:
            filename.unlink()

    return pd.read_parquet(whole_data_filename).tz_localize('UTC')
