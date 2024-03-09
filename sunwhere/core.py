'''
Author: Jose A Ruiz-Arias
Email: jararias at uma.es
Date: June 2017
Updated: June 2022 (Updated coefficients in PSA per Blanco et al., 2020)
Updated: April 2022 (Added support for NaT values)
Updated: April 2022 (Added support for numpy and numexpr. Removed multiproccesing support)
...
'''

# pylint: disable=c-extension-no-member
# pylint: disable=protected-access
# pylint: disable=global-variable-undefined

import numpy as np
from loguru import logger

from .utils import safe_import, validate


__ALGORITHMS__ = {
    'psa': {
        'numpy': safe_import('.py_psa', package='sunwhere.libspa'),
        'numexpr': safe_import('.ne_psa', package='sunwhere.libspa')
    },
    'soltrack': {
        'numpy': safe_import('.py_soltrack', package='sunwhere.libspa'),
        'numexpr': safe_import('.ne_soltrack', package='sunwhere.libspa')
    },
    'iqbal': {
        'numpy': safe_import('.py_iqbal', package='sunwhere.libspa'),
        'numexpr': safe_import('.ne_iqbal', package='sunwhere.libspa')
    },
    'nrel': {
        'numexpr': safe_import('.ne_nrel', package='sunwhere.libspa')
    }
}

logger.disable(__name__)


def evaluate(times, latitude, longitude, algorithm='psa', ndim=1,
             refraction=True, engine='numexpr'):
    """
    Calculate solar position from the observers' latitudes and longitudes,
    at any given time instants. The solar position is described in terms
    of solar zenith angle and solar azimuth angle.

    Parameters
    ----------
    times: scalar, sequence, or rank-1 array-like with shape (K,) of:
          Python datetime instances, Pandas Timestamp instances, Numpy
          datetime64 instances, or Pandas DatetimeIndex
        Time instants at which computing the solar position. If they are
        not timezone-aware, UTC is assumed.
    latitude: scalar, sequence, or rank-1 array-like with shape (J,)
        Latitudes (degrees, positive northward) where the solar position
        is computed. Latitude must be such that -90 <= latitude <= 90.
    longitude: scalar, sequence, or rank-1 array-like with shape (I,)
        Longitudes (degrees, positive eastward) where the solar position
        is computed. Longitude must be such that -180 <= longitude < 180.
    ndim: integer, {0, 1, 2}
        Calculation mode. It determines how times, latitudes and longitudes
        should be treated in order to reduce computation time.
        =0: when every location has a different reference timestamp. For
            instance, for remote sensing observations gathered from a moving
            platform, such as a satellite, that results in non-regular
            latitude-longitude grids and location-dependent scanning times.
            `times_utc`, `latitude` and `longitude` must all have the same
            size, i.e., K=J=I, which is also that of the outputs.
        =1: when solar position is evaluated at arbitrary locations throughout
            the same times. `latitude` and `longitude` must have the same size,
            i.e., J=I. The shape of the output solar zenith and azimuth angles
            is (K, J), while the shape of the location-independent output
            variables is (K,).
        =2: when solar position is evaluated throughout a regular grid of
            latitudes and longitudes and common times. The output solar zenith
            and azimuth angles have shape (K, J, I), while the location-independent
            variables have shape (K,).
    algorithm: string, {`psa`, `nrel`, `iqbal`}
        Solar position algorithm. Default: `psa`.
        =nrel: NREL's Solar Position Algorithm [1] for the period from the year
            -2000 to 6000. Its expected uncertainty is +/- 0.0003 degrees.
        =psa: Plataforma Solar de Almería's (PSA) algorithm [2] with updated coefficients
            for the period 2020-2050 [3] that reduce the average error to 0.0024 degrees.
        =iqbal: algorithm described in Iqbal, M. [^4]. It is less accurate than
            `psa`, especially for solar azimuth angle, and only slightly faster.
        =soltrack: similar in performance to psa [5]
    refraction: bool
        Atmospheric refraction correction.
    engine: string, {`numpy`, `numexpr`}
        Baseline code implementation to perform the calculations. `numexpr` is faster.

    Returns
    -------
    A Sunpos instance

    References
    ----------
    .. [1]: Reda, I. and Andreas, A.,2003. Solar Position Algorithm for Solar Radiation
            Applications. 55 pp.; NREL Report No. TP-560-34302, Revised January 2008.
            URL: http://www.nrel.gov/docs/fy08osti/34302.pdf
    .. [2]: Blanco-Muriel et al., 2001. Computing the solar vector. Solar Energy, 70(5),
            431-441. doi: 10.1016/S0038-092X(00)00156-0.
    .. [3]: Blanco et al., 2020. Updating the PSA sun position algorithm. Solar Energy,
            212, 339-341. doi: 10.1016/j.solener.2020.10.084.
    .. [4]: Iqbal, M. An introduction to solar radiation. Academic Press, 1983.
    .. [5]: van der Sluys M and van Kan P, 2022. SolTrack: a free, fast and accurate routine
            to compute the position of the Sun doi: 10.48550/arXiv.2209.01557
            https://github.com/MarcvdSluys/SolTrack-Python
    """

    if algorithm not in __ALGORITHMS__:
        raise ValueError(f'missing algorithm `{algorithm}`. Valid algorithms '
                         f'are: {", ".join(__ALGORITHMS__.keys())}')

    if engine not in __ALGORITHMS__.get(algorithm):
        raise ValueError(f'missing engine `{engine}` for algorithm '
                         f'`{algorithm}`. Valid engines for this algorithm '
                         f'are: {", ".join(__ALGORITHMS__.get(algorithm))}')

    logger.debug(f'selected {algorithm} algorithm and {engine} engine')

    # validate inputs...
    ndim = validate.check_ndim(ndim)
    times_utc, latitude, longitude = validate.check_timelatlon(
        times, latitude, longitude, ndim, np.float64)

    # NOTE: times_utc is a naive datetime64[ns], in UTC

    # NOTE: hereinafter, times_utc, latitude and longitude are
    # 1-dim array such that:
    #   if ndim == 0, times.shape == latitude.shape == longitude.shape
    #   if ndim in [1, 2], latitude.shape == longitude.shape

    # support for NaT in times_utc...
    nonat = ~np.isnat(times_utc)
    times_nonat = times_utc[nonat]
    lats_nonat = latitude[nonat] if ndim == 0 else latitude
    lons_nonat = longitude[nonat] if ndim == 0 else longitude

    # dimensions...
    n_times = times_nonat.size
    n_lats = lats_nonat.size
    n_lons = lons_nonat.size
    out_shape = {
        0: (n_times,),                # n_times == n_lats == n_lons
        1: (n_times, n_lats),         # n_lats == n_lons
        2: (n_times, n_lats, n_lons)
    }[ndim]

    # load function to compute solar position...
    spa_func = __ALGORITHMS__.get(algorithm).get(engine).sunpos
    spa_args = [ndim, refraction]

    # ... and call it
    namespace = spa_func.__module__ + '.' + spa_func.__name__
    logger.debug(f'running simulation with funtion {namespace}')
    result = spa_func(times_nonat, lons_nonat, lats_nonat, *spa_args)

    # fill with NaN where times is NaT...
    nonat = ~np.isnat(times_utc)
    n_times = times_utc.size
    n_lats = latitude.size
    n_lons = longitude.size
    out_shape = {
        0: (n_times,),                # n_times == n_lats == n_lons
        1: (n_times, n_lats),         # n_lats == n_lons
        2: (n_times, n_lats, n_lons)
    }[ndim]

    def reshape(variable):
        if (values := result.get(variable, None)) is None:
            return None

        if variable in ('ecf', 'eot', 'declination'):
            output = np.full((n_times,), np.nan)
            output[nonat] = values
            return output

        output = np.full(out_shape, np.nan)
        output[nonat, ...] = values
        return output

    # here, times, latitude and longitude are 1-dim
    # ecf, eot and declination have shape=(n_times,)
    # zenith and azimuth are: dim-1, when ndim=0 (i.e.,
    # transect); dim-2, when ndim=1 (i.e., sites); or
    # dim-3, when ndim=2 (i.e., regular_grid)

    return {
        'times': times,
        'latitude': latitude,
        'longitude': longitude,
        'algorithm': algorithm,
        'engine': engine,
        'refraction': refraction,
        'ecf': reshape("ecf"),
        'eot': reshape("eot"),
        'declination': reshape("declination"),
        'zenith': reshape("zenith"),
        'azimuth': reshape("azimuth")
    }
