
# pylint: disable=protected-access

import numpy as np
import pandas as pd

from ._core import evaluate
from ._base import Sunpos
from .utils.validate import check_dimensions


# TYPES OF TIME GRIDS:
#  (SYNCHRONOUS) MEANS ALL GRID CELLS SHARE THE SAME TIMESTAMPS
#  ASYNC(HRONOUS) MEANS THE GRID CELLS HAVE POTENTIALLY DIFFERENT TIMESTAMPS

# TYPES OF SPATIAL GRIDS:
#  REGULAR MEANS SAME LATITUDES ROW-WISE AND SAME LONGITUDES COLUMN-WISE
#  OTHERWISE, THE GRID IS SAID NON-REGULAR
#  TIME-INVARIANT IF THE LOCATIONS ARE FIXED, TIME-VARIANT OTHERWISE


def sites(times, latitude, longitude, algorithm='psa',
          refraction=True, engine='numexpr'):
    """Solar position across an arbitrary number of locations throughout a common time grid.

    Latitude and longitude *must* have the same size.

    Args:
        times (1-dim array of datetime64, or similar, of size K): Time
            instants at which solar position is evaluated. Similar types
            are, for instance, datetime and pandas DatetimeIndex.
        latitude (scalar or 1-dim array like of floats of size J): Latitudes
            (degrees) where solar position is evaluated. It must be in the
            range [-90, 90].
        longitude (scalar or 1-dim array like of floats of size J): Longitudes
            (degrees) where solar position is evaluated. It must be in the
            range [-180, 180].
        algorithm ({_nrel_, _psa_, _soltrack_, _iqbal_}): Solar position algorithm.
            _nrel_ is for NREL's SPA [1], valid for the years -2000 to 6000.
            Its expected uncertainty is +/- 0.0003 degrees. _psa_ is for the
            Plataforma Solar de Almería's (PSA) algorithm [2] with updated
            coefficients for the period 2020-2050 [3]. Its expected average
            error is 0.002 degrees, but it is faster than _nrel_. _soltrack_
            [4] is similar in performance to _psa_. _iqbal_ is for the algorithm
            described in Iqbal, M. [5], which has lower accuracy than the former
            algorithms.
        refraction (bool): Whether atmospheric refraction must be considered.
        engine ({_numpy_, _numexpr_}): Baseline code implementation to perform
            the solar position calculations. _numexpr_ is expected to be faster
            than numpy, especially for big spatial and/or temporal grids.

    Returns:
        A Sunpos instance in which the shape of the variables that are not
        location-dependent (e.g., declination) is (K,) and that of the
        location-dependent variables (e.g., zenith) is (K, J).

    References:

        [1] Reda I and Andreas A, 2003. Solar Position Algorithm for Solar
        Radiation Applications. 55 pp.; NREL Report No. TP-560-34302, Revised
        January 2008 [pdf](http://www.nrel.gov/docs/fy08osti/34302.pdf)
        [url](https://midcdmz.nrel.gov/spa/).

        [2] Blanco-Muriel, M. et al. 2001. Computing the solar vector. Solar
        Energy, Vol. 70, pp. 431-441
        doi: [10.1016/S0038-092X(00)00156-0](https://doi.org/10.1016/S0038-092X(00)00156-0).

        [3] Blanco, M. et al. 2020. Updating the PSA sun position algorithm.
        Solar Energy, Vol. 212, pp. 339-341
        doi: [10.1016/j.solener.2020.10.084](https://doi.org/10.1016/j.solener.2020.10.084).

        [4] van der Sluys M and van Kan P, 2022. SolTrack. A free, fast and accurate routine
        to compute the position of the Sun
        doi: [10.48550/arXiv.2209.01557](https://doi.org/10.48550/arXiv.2209.01557)
        [code](https://github.com/MarcvdSluys/SolTrack-Python)

        [5] Iqbal, M., An introduction to solar radiation. Academic Press. 1983
        [url](https://www.sciencedirect.com/book/9780123737502/an-introduction-to-solar-radiation)
    """
    # I use ndmin=1 to allow scalar inputs
    #   the argument `utc` of pd.to_datetime "localizes" timezone-naive
    #   inputs as UTC, while timezone-aware inputs are "converted" to UTC
    times_utc = np.array(
        pd.to_datetime(times, utc=True).tz_localize(None),  # naive datetime
        ndmin=1, dtype='datetime64[s]')
    sites_lats = np.array(latitude, ndmin=1, dtype=np.float64)
    sites_lons = np.array(longitude, ndmin=1, dtype=np.float64)

    check_dimensions(times_utc.ndim, sites_lats.ndim, sites_lons.ndim, (1, 1, 1))

    if sites_lats.shape != sites_lons.shape:
        raise ValueError(
            'shape mismatch: expected equal shape for latitude and '
            f'longitude, but got {sites_lats.shape} for latitude and '
            f'{sites_lons.shape} for longitude')

    # NOTE: ndim=1: two dimensional space: (n_times, n_locations)
    # NOTE: times_utc is a naive datetime64[s]

    solpos = evaluate(times_utc, sites_lats, sites_lons, algorithm,
                      ndim=1, refraction=refraction, engine=engine)

    return Sunpos(
        times=times,
        latitude=sites_lats,
        longitude=sites_lons,
        algorithm=algorithm,
        engine=engine,
        refraction=refraction,
        usecase='sites',
        ecf=solpos.get('ecf'),
        eot=solpos.get('eot'),
        declination=solpos.get('declination'),
        zenith=solpos.get('zenith'),
        azimuth=solpos.get('azimuth')
    )


def regular_grid(times, latitude, longitude, algorithm='psa',
                 refraction=True, engine='numexpr'):
    """Solar position across a lon-lat regular grid throughout a common time grid.

    Args:
        times (1-dim array of datetime64, or similar, of size K): Time
            instants at which solar position is evaluated. Similar types
            are, for instance, datetime and pandas DatetimeIndex.
        latitude (scalar or 1-dim array like of floats of size J): Latitudes
            (degrees) where solar position is evaluated. It must be in the
            range [-90, 90].
        longitude (scalar or 1-dim array like of floats of size I): Longitudes
            (degrees) where solar position is evaluated. It must be in the
            range [-180, 180].
        algorithm ({_nrel_, _psa_, _soltrack_, _iqbal_}): Solar position algorithm.
            _nrel_ is for NREL's SPA [1], valid for the years -2000 to 6000.
            Its expected uncertainty is +/- 0.0003 degrees. _psa_ is for the
            Plataforma Solar de Almería's (PSA) algorithm [2] with updated
            coefficients for the period 2020-2050 [3]. Its expected average
            error is 0.002 degrees, but it is faster than _nrel_. _soltrack_
            [4] is similar in performance to _psa_. _iqbal_ is for the algorithm
            described in Iqbal, M. [5], which has lower accuracy than the former
            algorithms.
        refraction (bool): Whether atmospheric refraction must be considered.
        engine ({_numpy_, _numexpr_}): Baseline code implementation to perform
            the solar position calculations. _numexpr_ is expected to be faster
            than numpy, especially for big spatial and/or temporal grids.

    Returns:
        A Sunpos instance in which the shape of the variables that are not
        location-dependent (e.g., declination) is (K,) and that of the
        location-dependent variables (e.g., zenith) is (K, J, I).

    References:

        [1] Reda I and Andreas A, 2003. Solar Position Algorithm for Solar
        Radiation Applications. 55 pp.; NREL Report No. TP-560-34302, Revised
        January 2008 [pdf](http://www.nrel.gov/docs/fy08osti/34302.pdf)
        [url](https://midcdmz.nrel.gov/spa/).

        [2] Blanco-Muriel, M. et al. 2001. Computing the solar vector. Solar
        Energy, Vol. 70, pp. 431-441
        doi: [10.1016/S0038-092X(00)00156-0](https://doi.org/10.1016/S0038-092X(00)00156-0).

        [3] Blanco, M. et al. 2020. Updating the PSA sun position algorithm.
        Solar Energy, Vol. 212, pp. 339-341
        doi: [10.1016/j.solener.2020.10.084](https://doi.org/10.1016/j.solener.2020.10.084).

        [4] van der Sluys M and van Kan P, 2022. SolTrack. A free, fast and accurate routine
        to compute the position of the Sun
        doi: [10.48550/arXiv.2209.01557](https://doi.org/10.48550/arXiv.2209.01557)
        [code](https://github.com/MarcvdSluys/SolTrack-Python)

        [5] Iqbal, M., An introduction to solar radiation. Academic Press. 1983
        [url](https://www.sciencedirect.com/book/9780123737502/an-introduction-to-solar-radiation)
    """
    # I use ndmin=1 to allow scalar input times
    #   the argument `utc` of pd.to_datetime "localizes" timezone-naive
    #   inputs as UTC, while timezone-aware inputs are "converted" to UTC
    times_utc = np.array(
        pd.to_datetime(times, utc=True).tz_localize(None),  # naive datetime
        ndmin=1, dtype='datetime64[s]')
    grid_lats = np.array(latitude, ndmin=1, dtype=np.float64)
    grid_lons = np.array(longitude, ndmin=1, dtype=np.float64)

    check_dimensions(times_utc.ndim, grid_lats.ndim, grid_lons.ndim, (1, 1, 1))

    # ndim=2: three dimensional space: (n_times, n_lats, n_lons)

    solpos = evaluate(times, grid_lats, grid_lons, algorithm,
                      ndim=2, refraction=refraction, engine=engine)

    return Sunpos(
        times=times,
        latitude=grid_lats,
        longitude=grid_lons,
        algorithm=algorithm,
        engine=engine,
        refraction=refraction,
        usecase='regular_grid',
        ecf=solpos.get('ecf'),
        eot=solpos.get('eot'),
        declination=solpos.get('declination'),
        zenith=solpos.get('zenith'),
        azimuth=solpos.get('azimuth')
    )


def transect(times, latitude, longitude, algorithm='psa',
             refraction=True, engine='numexpr'):
    """Solar position across a transect.

    In a transect, the observer's position changes throughout time.

    Args:
        times (1-dim array of datetime64, or similar, of size K): Time
            instants at which solar position is evaluated. Similar types
            are, for instance, datetime and pandas DatetimeIndex.
        latitude (scalar or 1-dim array like of floats of size K): Latitudes
            (degrees) where solar position is evaluated. It must be in the
            range [-90, 90].
        longitude (scalar or 1-dim array like of floats of size K): Longitudes
            (degrees) where solar position is evaluated. It must be in the
            range [-180, 180].
        algorithm ({_nrel_, _psa_, _soltrack_, _iqbal_}): Solar position algorithm.
            _nrel_ is for NREL's SPA [1], valid for the years -2000 to 6000.
            Its expected uncertainty is +/- 0.0003 degrees. _psa_ is for the
            Plataforma Solar de Almería's (PSA) algorithm [2] with updated
            coefficients for the period 2020-2050 [3]. Its expected average
            error is 0.002 degrees, but it is faster than _nrel_. _soltrack_
            [4] is similar in performance to _psa_. _iqbal_ is for the algorithm
            described in Iqbal, M. [5], which has lower accuracy than the former
            algorithms.
        refraction (bool): Whether atmospheric refraction must be considered.
        engine ({_numpy_, _numexpr_}): Baseline code implementation to perform
            the solar position calculations. _numexpr_ is expected to be faster
            than numpy, especially for big spatial and/or temporal grids.

    Returns:
        A Sunpos instance in which the shape of all variables is (K,).

    References:

        [1] Reda I and Andreas A, 2003. Solar Position Algorithm for Solar
        Radiation Applications. 55 pp.; NREL Report No. TP-560-34302, Revised
        January 2008 [pdf](http://www.nrel.gov/docs/fy08osti/34302.pdf)
        [url](https://midcdmz.nrel.gov/spa/).

        [2] Blanco-Muriel, M. et al. 2001. Computing the solar vector. Solar
        Energy, Vol. 70, pp. 431-441
        doi: [10.1016/S0038-092X(00)00156-0](https://doi.org/10.1016/S0038-092X(00)00156-0).

        [3] Blanco, M. et al. 2020. Updating the PSA sun position algorithm.
        Solar Energy, Vol. 212, pp. 339-341
        doi: [10.1016/j.solener.2020.10.084](https://doi.org/10.1016/j.solener.2020.10.084).

        [4] van der Sluys M and van Kan P, 2022. SolTrack. A free, fast and accurate routine
        to compute the position of the Sun
        doi: [10.48550/arXiv.2209.01557](https://doi.org/10.48550/arXiv.2209.01557)
        [code](https://github.com/MarcvdSluys/SolTrack-Python)

        [5] Iqbal, M., An introduction to solar radiation. Academic Press. 1983
        [url](https://www.sciencedirect.com/book/9780123737502/an-introduction-to-solar-radiation)
    """
    # I use ndmin=1 to allow scalar inputs
    #   the argument `utc` of pd.to_datetime "localizes" timezone-naive
    #   inputs as UTC, while timezone-aware inputs are "converted" to UTC
    times_utc = np.array(
        pd.to_datetime(times, utc=True).tz_localize(None),
        ndmin=1, dtype='datetime64[s]')
    transect_lats = np.array(latitude, ndmin=1, dtype=np.float64)
    transect_lons = np.array(longitude, ndmin=1, dtype=np.float64)

    check_dimensions(times_utc.ndim, transect_lats.ndim, transect_lons.ndim, (1, 1, 1))

    if not (transect_lats.shape == transect_lons.shape == times_utc.shape):
        raise ValueError(
            'shape mismatch: expected equal shape for times, latitude and '
            f'longitude, but got {times_utc.shape} for times, {transect_lats.shape} '
            f'for latitude and {transect_lons.shape} for longitude')

    # ndim=0: one dimensional space: (n_times,) == (n_lats,) == (n_lons,)

    solpos = evaluate(times, transect_lats, transect_lons, algorithm,
                      ndim=0, refraction=refraction, engine=engine)

    return Sunpos(
        times=times,
        latitude=transect_lats,
        longitude=transect_lons,
        algorithm=algorithm,
        engine=engine,
        refraction=refraction,
        usecase='transect',
        ecf=solpos.get('ecf'),
        eot=solpos.get('eot'),
        declination=solpos.get('declination'),
        zenith=solpos.get('zenith'),
        azimuth=solpos.get('azimuth')
    )
