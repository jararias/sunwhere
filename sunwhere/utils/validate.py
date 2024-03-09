
import numpy as np
import pandas as pd


def check_ndim(ndim):
    ndim = int(ndim)
    if ndim not in (0, 1, 2):
        raise ValueError(
            f'expected ndim in (0, 1, 2) but got {ndim}')
    return ndim


def check_timelatlon(times, latitude, longitude, ndim, dtype=None, dt64=None):
    """
    Check that times, latitude and longitude are 0-dim or 1-dim
    array-like instances, that latitude and longitude values are
    not out of allowed bounds, and that the shapes are consistent
    with the value of ndim
    """

    dtype = dtype or np.float64
    dt64 = dt64 or np.datetime64(1, 's')  # equivalent to 'datetime64[s]'

    # validate times...
    #   the argument `utc` of pd.to_datetime "localizes" timezone-naive
    #   inputs as UTC, while timezone-aware inputs are "converted" to UTC
    times = np.array(
        pd.to_datetime(times, utc=True).tz_localize(None),  # naive datetime
        ndmin=1, dtype=dt64)
    if times.ndim > 1:
        raise ValueError(
            'illegal shape: expected a 0-dim or 1-dim array of '
            f'times but got a {times.ndim}-dim array')

    # validate longitudes...
    longitude = np.array(longitude, ndmin=1, dtype=dtype)
    if longitude.ndim > 1:
        raise ValueError(
            'illegal shape: expected a 0-dim or 1-dim array of '
            f'longitudes but got a {longitude.ndim}-dim array')

    # the `where` argument is included from numpy 1.20 onwards, but
    # sun_geometry requires older versions
    # if not (np.all(-180 <= longitude, where=~np.isnan(longitude)) and
    #         np.all(longitude < +180., where=~np.isnan(longitude))):

    if not (np.all(-180 <= longitude[~np.isnan(longitude)]) and
            np.all(longitude[~np.isnan(longitude)] < +180.)):
        raise ValueError(
            'longitude out of bounds: expected -180 <= longitude < 180 '
            'but got values beyond')

    # validate latitudes...
    latitude = np.array(latitude, ndmin=1, dtype=dtype)
    if latitude.ndim > 1:
        raise ValueError(
            'illegal shape: expected a 0-dim or 1-dim array of '
            f'latitudes but but got {longitude.ndim}-dim array')

    # if not (np.all(-90 <= latitude, where=~np.isnan(latitude)) and
    #         np.all(latitude <= +90., where=~np.isnan(latitude))):

    if not (np.all(-90 <= latitude[~np.isnan(latitude)]) and
            np.all(latitude[~np.isnan(latitude)] <= +90.)):
        raise ValueError(
            'latitude out of bounds: expected -90 <= latitude <= 90 '
            'but got values beyond')

    # validate consistency against ndim...
    if (ndim == 0) and not times.shape == longitude.shape == latitude.shape:
        raise ValueError(
            'shape mismatch: expected equal shape for times, latitude '
            f'and longitude, but got {times.shape} for times, '
            f'{longitude.shape} for longitude and {latitude.shape} for '
            f'latitude')

    if (ndim == 1) and longitude.shape != latitude.shape:
        raise ValueError(
            'shape mismatch: expected equal shape for latitude and '
            f'longitude, but got {longitude.shape} for longitude '
            f'and {latitude.shape} for latitude')

    return times, latitude, longitude


# def is_lonlat_scalar(longitude, latitude):
#     try:
#         longitude.item()
#         latitude.item()
#         return True
#     except ValueError:
#         return False


def check_dimensions(times_dim, lats_dim, lons_dim, required_dims):
    req_times_dim, req_lats_dim, req_lons_dim = required_dims
    if not ((times_dim == req_times_dim) and
            (lats_dim == req_lats_dim) and
            (lons_dim == req_lons_dim)):
        raise ValueError(
            f'wrong dimensions: expected {req_times_dim}-dim times_utc, '
            f'got {times_dim}-dim; expected {req_lats_dim}-dim latitude, '
            f'got {lats_dim}-dim; expected {req_lons_dim}-dim longitude, '
            f'get {lons_dim}-dim')
