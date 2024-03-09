
import numpy as np

from .._core import evaluate


def true_solar_time(times_utc, longitude):
    dt64 = 'datetime64[ns]'
    times_utc = np.array(times_utc, ndmin=1)
    longitude = np.array(longitude, ndmin=1)
    eot = evaluate(times_utc, 0., 0.).eot  # eot, minutes
    # for compat with numpy < 1.18.0, I use reshape instead of expand_dims
    new_shape = tuple(list(eot.shape) + [1]*longitude.ndim)
    expanded_eot = np.reshape(eot, new_shape)
    expanded_times_utc = np.reshape(
        np.array(times_utc, dtype=dt64).astype(np.float128),
        new_shape
    )
    expanded_longitude = np.expand_dims(longitude, axis=0)
    utc_to_tst = 60 * (4 * expanded_longitude + expanded_eot)
    times_tst = np.array(expanded_times_utc + 1e9 * utc_to_tst, dtype=dt64)
    return np.squeeze(times_tst)


def universal_time_coordinated(times_tst, longitude):
    dt64 = 'datetime64[ns]'
    times_tst = np.array(times_tst, ndmin=1)
    longitude = np.array(longitude, ndmin=1)
    eot = evaluate(times_tst, 0., 0.).eot  # eot, minutes
    # for compat with numpy < 1.18.0, I use reshape instead of expand_dims
    new_shape = tuple(list(eot.shape) + [1]*longitude.ndim)
    expanded_eot = np.reshape(eot, new_shape)
    expanded_times_tst = np.reshape(
        np.array(times_tst, dtype=dt64).astype(np.float128),
        new_shape
    )
    expanded_longitude = np.expand_dims(longitude, axis=0)
    tst_to_utc = -60 * (4 * expanded_longitude + expanded_eot)
    times_utc = np.array(expanded_times_tst + 1e9 * tst_to_utc, dtype=dt64)
    return np.squeeze(times_utc)
