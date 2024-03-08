
from datetime import datetime

import numpy as np
import pandas as pd

import sunwhere


def calculate_position(time, lat, lon, algorithm, refraction):
    # TODO: decorate input & output!!

    if time == 'now':
        time = datetime.now()

    times = np.array(pd.to_datetime(time, utc=True).to_numpy(),
                     ndmin=1, dtype='datetime64[ns]')
    sw = sunwhere.sites(times, lat, lon, algorithm=algorithm, refraction=refraction)
    print(
        time,
        sw.zenith.isel(time=0, location=0).to_numpy(),
        sw.azimuth.isel(time=0, location=0).to_numpy(),
        sw.dec.isel(time=0).to_numpy(),
        sw.sunrise(units='hour').isel(time=0, location=0).to_numpy(),
        sw.sunset(units='hour').isel(time=0, location=0).to_numpy()
    )
