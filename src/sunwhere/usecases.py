
# pylint: disable=protected-access

from typing import Literal, Sequence

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


def sites(
    times: np.ndarray[tuple[int], np.datetime64] | pd.DatetimeIndex,
    latitude: np.ndarray[tuple[int]] | Sequence[float] | float,
    longitude: np.ndarray[tuple[int]] | Sequence[float] | float,
    algorithm: Literal['psa', 'nrel', 'iqbal'] = 'psa',
    refraction: bool = True,
    engine: Literal['numexpr', 'numpy'] = 'numexpr',
    site_names: Sequence[str] | None = None,
) -> Sunpos:
    """
    Compute solar position at arbitrarily dispersed locations over a common time series.

    This function calculates solar position (zenith and azimuth angles, declination,
    etc.) for an arbitrary number of fixed geographic locations throughout a common
    time grid. All locations share the same time instants.

    Parameters
    ----------
    times : array-like of datetime-like, shape (K,)
        Time instants at which solar position is evaluated. Can be datetime.datetime,
        pandas.DatetimeIndex, numpy.datetime64, or any type accepted by
        pandas.to_datetime(). Timezone-naive inputs are assumed to be UTC.
    latitude : float or array-like of float, shape (J,)
        Latitudes in degrees where solar position is evaluated.
        Valid range: [-90, 90]. Must have the same size as longitude.
    longitude : float or array-like of float, shape (J,)
        Longitudes in degrees where solar position is evaluated.
        Valid range: [-180, 180). Must have the same size as latitude.
    algorithm : {'psa', 'nrel', 'iqbal'}, optional
        Solar position algorithm to use:
        - 'psa': Plataforma Solar de Almería algorithm [1] with updated
          coefficients [2]. Fast and accurate (avg error ~0.002°). Default.
        - 'nrel': NREL's SPA algorithm [3]. Most accurate (±0.0003°) but slower.
          Valid for years -2000 to 6000.
        - 'iqbal': Iqbal's algorithm [4]. Lower accuracy, mainly for educational use.
    refraction : bool, optional
        If True (default), applies atmospheric refraction correction to zenith angle.
        Typically adds ~0.5° correction near the horizon.
    engine : {'numexpr', 'numpy'}, optional
        Computation engine. 'numexpr' (default) is faster for large arrays,
        'numpy' may be faster for small arrays.
    site_names : array-like of str, shape (J,), optional
        Custom names for each site. If provided, must have the same length as
        latitude and longitude. These names will be used as a coordinate for the
        'site' dimension, enabling labeled selection like ``result.sza.sel(site='Madrid')``.
        If None (default), sites are indexed numerically (0, 1, 2, ...).

    Returns
    -------
    sunpos : Sunpos
        Object containing solar position data as xarray DataArrays with dimensions:
        
        - Time-dependent only (K,): dec, eot, ecf
        - Time and location dependent (K, J): sza, saa, elevation, azimuth, zenith, cosz

    Notes
    -----
    - Latitude and longitude arrays must have the same length.
    - For a single location, use scalar values for latitude and longitude.
    - All times must be unique and monotonically increasing for best performance.
    - The refraction correction uses a simple model suitable for standard atmospheric
      conditions at sea level.

    Examples
    --------
    Calculate solar position for three cities over one day:

    >>> import pandas as pd
    >>> import sunwhere
    >>> 
    >>> # Define locations: Madrid, New York, Tokyo
    >>> cities = ['Madrid', 'New York', 'Tokyo']
    >>> lats = [40.4168, 40.7128, 35.6762]
    >>> lons = [-3.7038, -74.0060, 139.6503]
    >>> 
    >>> # Create hourly time series for 2024-06-21 (summer solstice)
    >>> times = pd.date_range('2024-06-21', periods=24, freq='h', tz='UTC')
    >>> 
    >>> # Compute solar position with site names
    >>> result = sunwhere.sites(times, lats, lons, site_names=cities)
    >>> 
    >>> # Access solar zenith angle with labeled coordinates
    >>> print(result.sza)  # shape: (24, 3) with site names as coordinate
    >>> 
    >>> # Select data for a specific city using name
    >>> madrid_sza = result.sza.sel(site='Madrid')
    >>> 
    >>> # Find solar noon (minimum zenith angle) for Madrid
    >>> solar_noon_idx = result.sza.sel(site='Madrid').argmin()
    >>> print(f"Solar noon in Madrid: {times[solar_noon_idx]}")

    Single location example:

    >>> # Single location (scalar coordinates)
    >>> times = pd.date_range('2024-01-01', periods=365, freq='D', tz='UTC')
    >>> result = sunwhere.sites(times, 40.0, -3.0, algorithm='nrel')
    >>> 
    >>> # Get solar azimuth angle
    >>> azimuth = result.saa  # shape: (365, 1)

    References
    ----------
    1. Blanco-Muriel, M. et al., 2001. Computing the solar vector.
       Solar Energy, 70(5), 431-441.
       https://doi.org/10.1016/S0038-092X(00)00156-0

    2. Blanco, M. et al., 2020. Updating the PSA sun position algorithm.
       Solar Energy, 212, 339-341.
       https://doi.org/10.1016/j.solener.2020.10.084

    3. Reda, I. and Andreas, A., 2003. Solar Position Algorithm for Solar
       Radiation Applications. NREL Report No. TP-560-34302.
       https://www.nrel.gov/docs/fy08osti/34302.pdf

    4. Iqbal, M., 1983. An introduction to solar radiation. Academic Press.

    See Also
    --------
    regular_grid : Solar position over a regular lat-lon grid
    transect : Solar position along a moving path
    """
    # I use ndmin=1 to allow scalar inputs
    #   the argument `utc` of pd.to_datetime "localizes" timezone-naive
    #   inputs as UTC, while timezone-aware inputs are "converted" to UTC
    times_utc = np.array(
        pd.to_datetime(times, utc=True).tz_localize(None),  # naive datetime, UTC
        ndmin=1, dtype='datetime64[ns]')
    sites_lats = np.array(latitude, ndmin=1, dtype=np.float64)
    sites_lons = np.array(longitude, ndmin=1, dtype=np.float64)

    check_dimensions(times_utc.ndim, sites_lats.ndim, sites_lons.ndim, (1, 1, 1))

    if sites_lats.shape != sites_lons.shape:
        raise ValueError(
            'shape mismatch: expected equal shape for latitude and '
            f'longitude, but got {sites_lats.shape} for latitude and '
            f'{sites_lons.shape} for longitude')

    # Validate site_names if provided
    if site_names is not None:
        site_names = np.array(site_names, ndmin=1, dtype=str)
        if site_names.ndim != 1:
            raise ValueError(
                f'expected 1-dim array for site_names, but got {site_names.ndim}-dim array')
        if site_names.shape[0] != sites_lats.shape[0]:
            raise ValueError(
                f'shape mismatch: site_names must have the same length as latitude and longitude, '
                f'but got {site_names.shape[0]} for site_names and {sites_lats.shape[0]} for latitude/longitude')

    # NOTE: ndim=1: two dimensional space: (n_times, n_locations)
    # NOTE: times_utc is a naive datetime64[s], UTC

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
        site_names=site_names,
        ecf=solpos.get('ecf'),
        eot=solpos.get('eot'),
        declination=solpos.get('declination'),
        zenith=solpos.get('zenith'),
        azimuth=solpos.get('azimuth')
    )


def regular_grid(
    times: np.ndarray[tuple[int], np.datetime64] | pd.DatetimeIndex,
    latitude: np.ndarray[tuple[int]] | Sequence[float] | float,
    longitude: np.ndarray[tuple[int]] | Sequence[float] | float,
    algorithm: Literal['psa', 'nrel', 'iqbal'] = 'psa',
    refraction: bool = True,
    engine: Literal['numexpr', 'numpy'] = 'numexpr',
) -> Sunpos:
    """
    Compute solar position over a regular latitude-longitude grid.

    This function calculates solar position for a regular 2D grid of latitudes
    and longitudes throughout a common time series. The output has shape
    (time, latitude, longitude), suitable for gridded climate data or maps.

    Parameters
    ----------
    times : array-like of datetime-like, shape (K,)
        Time instants at which solar position is evaluated. Can be datetime.datetime,
        pandas.DatetimeIndex, numpy.datetime64, or any type accepted by
        pandas.to_datetime(). Timezone-naive inputs are assumed to be UTC.
    latitude : float or array-like of float, shape (J,)
        1-D array of latitudes in degrees defining the grid rows.
        Valid range: [-90, 90]. Can be a scalar for a single latitude.
    longitude : float or array-like of float, shape (I,)
        1-D array of longitudes in degrees defining the grid columns.
        Valid range: [-180, 180). Can be a scalar for a single longitude.
    algorithm : {'psa', 'nrel', 'iqbal'}, optional
        Solar position algorithm to use:
        
        - 'psa': Plataforma Solar de Almería algorithm [2] with updated
          coefficients [3]. Fast and accurate (avg error ~0.002°). Default.
        - 'nrel': NREL's SPA algorithm [1]. Most accurate (±0.0003°) but slower.
          Valid for years -2000 to 6000.
        - 'iqbal': Iqbal's algorithm [4]. Lower accuracy, mainly for educational use.
    refraction : bool, optional
        If True (default), applies atmospheric refraction correction to zenith angle.
        Typically adds ~0.5° correction near the horizon.
    engine : {'numexpr', 'numpy'}, optional
        Computation engine. 'numexpr' (default) is significantly faster for large
        grids (>1000 points), 'numpy' may be faster for small grids.

    Returns
    -------
    sunpos : Sunpos
        Object containing solar position data as xarray DataArrays with dimensions:
        
        - Time-dependent only (K,): dec, eot, ecf
        - Time and space dependent (K, J, I): sza, saa, elevation, azimuth, zenith, cosz

    Notes
    -----
    - This function is optimized for regular grids where latitudes and longitudes
      form a rectangular mesh.
    - For irregular or scattered points, use `sites()` instead.
    - The grid is created using NumPy broadcasting: all combinations of lat×lon
      are computed for each time instant.
    - Memory usage scales as K × J × I × 8 bytes per output variable.

    Examples
    --------
    Create a global grid and compute solar zenith angle:

    >>> import numpy as np
    >>> import pandas as pd
    >>> import sunwhere
    >>> 
    >>> # Define a 5° resolution global grid
    >>> lats = np.arange(-90, 91, 5.0)
    >>> lons = np.arange(-180, 180, 5.0)
    >>> 
    >>> # One day at noon UTC
    >>> times = pd.date_range('2024-06-21 12:00', periods=1, freq='h', tz='UTC')
    >>> 
    >>> # Compute solar position
    >>> result = sunwhere.regular_grid(times, lats, lons)
    >>> 
    >>> # Solar zenith angle has shape (1, 37, 72)
    >>> print(result.sza.shape)
    (1, 37, 72)
    >>> 
    >>> # Create a map of solar zenith angle
    >>> import matplotlib.pyplot as plt
    >>> result.sza[0].plot(x='lon', y='lat')
    >>> plt.title('Solar Zenith Angle at Summer Solstice Noon UTC')
    >>> plt.show()

    High-resolution regional grid:

    >>> # European domain
    >>> lats = np.linspace(35, 70, 141)  # 0.25° resolution
    >>> lons = np.linspace(-10, 40, 201)
    >>> 
    >>> # Full year, daily resolution
    >>> times = pd.date_range('2024-01-01', '2024-12-31', freq='D', tz='UTC')
    >>> 
    >>> # Use NREL algorithm for highest accuracy
    >>> result = sunwhere.regular_grid(times, lats, lons, algorithm='nrel')
    >>> 
    >>> # Compute average zenith angle over the year
    >>> annual_mean_sza = result.sza.mean(dim='time')

    References
    ----------
    1. Reda, I. and Andreas, A., 2003. Solar Position Algorithm for Solar
       Radiation Applications. NREL Report No. TP-560-34302.
       https://www.nrel.gov/docs/fy08osti/34302.pdf
    2. Blanco-Muriel, M. et al., 2001. Computing the solar vector.
       Solar Energy, 70(5), 431-441.
       https://doi.org/10.1016/S0038-092X(00)00156-0
    3. Blanco, M. et al., 2020. Updating the PSA sun position algorithm.
       Solar Energy, 212, 339-341.
       https://doi.org/10.1016/j.solener.2020.10.084
    4. Iqbal, M., 1983. An introduction to solar radiation. Academic Press.

    See Also
    --------
    sites : Solar position at multiple fixed locations
    transect : Solar position along a moving path
    """
    # I use ndmin=1 to allow scalar input times
    #   the argument `utc` of pd.to_datetime "localizes" timezone-naive
    #   inputs as UTC, while timezone-aware inputs are "converted" to UTC
    times_utc = np.array(
        pd.to_datetime(times, utc=True).tz_localize(None),  # naive datetime, UTC
        ndmin=1, dtype='datetime64[s]')
    grid_lats = np.array(latitude, ndmin=1, dtype=np.float64)
    grid_lons = np.array(longitude, ndmin=1, dtype=np.float64)

    check_dimensions(times_utc.ndim, grid_lats.ndim, grid_lons.ndim, (1, 1, 1))

    # ndim=2: three dimensional space: (n_times, n_lats, n_lons)

    solpos = evaluate(times_utc, grid_lats, grid_lons, algorithm,
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


def transect(
    times: np.ndarray[tuple[int], np.datetime64] | pd.DatetimeIndex,
    latitude: np.ndarray[tuple[int]] | Sequence[float] | float,
    longitude: np.ndarray[tuple[int]] | Sequence[float] | float,
    algorithm: Literal['psa', 'nrel', 'iqbal'] = 'psa',
    refraction: bool = True,
    engine: Literal['numexpr', 'numpy'] = 'numexpr',
) -> Sunpos:
    """
    Compute solar position along a moving path (transect).

    This function calculates solar position for a trajectory where each time
    instant corresponds to a different geographic location. The output has shape
    (time,), suitable for moving observers like satellites, aircraft, or ships.

    Parameters
    ----------
    times : array-like of datetime-like, shape (K,)
        Time instants at which solar position is evaluated. Can be datetime.datetime,
        pandas.DatetimeIndex, numpy.datetime64, or any type accepted by
        pandas.to_datetime(). Timezone-naive inputs are assumed to be UTC.
    latitude : float or array-like of float, shape (K,)
        Latitude in degrees for each time instant. Valid range: [-90, 90].
        If scalar, the same latitude is used for all times (meridional transect).
    longitude : float or array-like of float, shape (K,)
        Longitude in degrees for each time instant. Valid range: [-180, 180).
        If scalar, the same longitude is used for all times (zonal transect).
    algorithm : {'psa', 'nrel', 'iqbal'}, optional
        Solar position algorithm to use:
        - 'psa': Plataforma Solar de Almería algorithm [2] with updated
          coefficients [3]. Fast and accurate (avg error ~0.002°). Default.
        - 'nrel': NREL's SPA algorithm [1]. Most accurate (±0.0003°) but slower.
          Valid for years -2000 to 6000.
        - 'iqbal': Iqbal's algorithm [4]. Lower accuracy, mainly for educational use.
    refraction : bool, optional
        If True (default), applies atmospheric refraction correction to zenith angle.
        Typically adds ~0.5° correction near the horizon.
    engine : {'numexpr', 'numpy'}, optional
        Computation engine. 'numexpr' (default) is typically faster for long
        transects (>100 points), 'numpy' may be faster for short transects.

    Returns
    -------
    sunpos : Sunpos
        Object containing solar position data as xarray DataArrays with dimension
        (K,) for all variables: dec, eot, ecf, sza, saa, elevation, azimuth, 
        zenith, cosz.

    Notes
    -----
    - This function is designed for scenarios where position changes with time,
      such as satellite ground tracks, flight paths, or ship routes.
    - If latitude and longitude are both scalars, all times share the same location
      (equivalent to a single site with varying times).
    - For multiple fixed locations, use `sites()` instead for better performance.
    - The time array and position arrays must have the same length if both are arrays.

    Examples
    --------
    Satellite ground track over one orbit:

    >>> import numpy as np
    >>> import pandas as pd
    >>> import sunwhere
    >>> 
    >>> # ISS-like orbit: ~90 minute period
    >>> times = pd.date_range('2024-01-15 00:00', periods=100, freq='54s', tz='UTC')
    >>> 
    >>> # Simplified sinusoidal ground track (not physically accurate)
    >>> lats = 51.6 * np.sin(2 * np.pi * np.arange(100) / 100)
    >>> lons = np.linspace(-180, 180, 100, endpoint=False)
    >>> 
    >>> # Compute solar position along the track
    >>> result = sunwhere.transect(times, lats, lons)
    >>> 
    >>> # Check when satellite is in sunlight (elevation > 0)
    >>> in_sunlight = result.elevation > 0
    >>> print(f"Sunlit for {in_sunlight.sum().item()} of {len(times)} points")

    Aircraft flight path from New York to Tokyo:

    >>> # Great circle approximation (10-hour flight)
    >>> times = pd.date_range('2024-03-20 10:00', periods=121, freq='5min', tz='UTC')
    >>> 
    >>> # Linear interpolation (simplified - not actual great circle)
    >>> lats = np.linspace(40.7, 35.7, 121)  # NYC to Tokyo
    >>> lons = np.linspace(-74.0, 139.7, 121)
    >>> 
    >>> # Compute solar position
    >>> result = sunwhere.transect(times, lats, lons, algorithm='nrel')
    >>> 
    >>> # Find solar elevation throughout flight
    >>> import matplotlib.pyplot as plt
    >>> result.elevation.plot()
    >>> plt.axhline(0, color='k', linestyle='--', label='Horizon')
    >>> plt.ylabel('Solar Elevation (degrees)')
    >>> plt.title('Sun Position During NYC-Tokyo Flight')
    >>> plt.legend()
    >>> plt.show()

    Meridional transect at fixed longitude:

    >>> # Moving north along Prime Meridian
    >>> times = pd.date_range('2024-06-21 12:00', periods=181, freq='h', tz='UTC')
    >>> lats = np.linspace(-90, 90, 181)  # South Pole to North Pole
    >>> lon = 0.0  # Prime Meridian (scalar)
    >>> 
    >>> result = sunwhere.transect(times, lats, lon)
    >>> # Solar zenith angle variation from pole to pole
    >>> result.sza.plot()

    References
    ----------
    1. Reda, I. and Andreas, A., 2003. Solar Position Algorithm for Solar
       Radiation Applications. NREL Report No. TP-560-34302.
       https://www.nrel.gov/docs/fy08osti/34302.pdf
    2. Blanco-Muriel, M. et al., 2001. Computing the solar vector.
       Solar Energy, 70(5), 431-441.
       https://doi.org/10.1016/S0038-092X(00)00156-0
    3. Blanco, M. et al., 2020. Updating the PSA sun position algorithm.
       Solar Energy, 212, 339-341.
       https://doi.org/10.1016/j.solener.2020.10.084
    4. Iqbal, M., 1983. An introduction to solar radiation. Academic Press.

    See Also
    --------
    sites : Solar position at multiple fixed locations
    regular_grid : Solar position over a latitude-longitude grid
    """
    # I use ndmin=1 to allow scalar inputs
    #   the argument `utc` of pd.to_datetime "localizes" timezone-naive
    #   inputs as UTC, while timezone-aware inputs are "converted" to UTC
    times_utc = np.array(
        pd.to_datetime(times, utc=True).tz_localize(None),
        ndmin=1, dtype='datetime64[s]')  # naive datetime, UTC
    transect_lats = np.array(latitude, ndmin=1, dtype=np.float64)
    transect_lons = np.array(longitude, ndmin=1, dtype=np.float64)

    check_dimensions(times_utc.ndim, transect_lats.ndim, transect_lons.ndim, (1, 1, 1))

    if not (transect_lats.shape == transect_lons.shape == times_utc.shape):
        raise ValueError(
            'shape mismatch: expected equal shape for times, latitude and '
            f'longitude, but got {times_utc.shape} for times, {transect_lats.shape} '
            f'for latitude and {transect_lons.shape} for longitude')

    # ndim=0: one dimensional space: (n_times,) == (n_lats,) == (n_lons,)

    solpos = evaluate(times_utc, transect_lats, transect_lons, algorithm,
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
