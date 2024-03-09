
# pylint: disable=broad-except

# https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html

import numpy as np
import pandas as pd
import xarray as xr


class Sunpos:
    """Container class for [sunwhere](https://github.com/jararias/sunwhere) outputs.

    Do not instantiate!!
    """

    def __init__(self, times, latitude, longitude,
                 algorithm, engine, refraction, usecase,
                 ecf, eot, declination, zenith, azimuth):
        """Creates a Sunpos' instance.

        Provides access to the solar geometry parameters, generally as
        [xarray's DataArrays](https://docs.xarray.dev/en/stable/generated/xarray.DataArray.html).

        Args:
            times (sequence of numpy datetime64, or convertible to it): times
              where solar geometry is evaluated.
            latitude (sequence of floats): latitude of the locations where
              solar geometry is evaluated. Must be in the range [-90, 90].
            longitude (sequence of floats): longitude of the locations where
              solar geometry is evaluated. Must be in the range [-180, 180).
            algorithm (str): solar position algorithm: nrel, psa or
              iqbal.
            engine (str): code implementation: numpy or numexpr
            refraction (bool): whether atmospheric refraction is considered
            usecase (str): sites, regular_grid or transect
            ecf (1-dim array-like of floats): sun-earth distance (eccentricity)
              correction factor.
            eot (1-dim array-like of floats): equation of time, minutes
            declination (1-dim array-like of floats): solar declination,
              degrees
            zenith (array of floats): solar zenith angle, degrees (1-dim for
              transect, 2-dim for sites and 3-dim for regular_grid)
            azimuth (array of floats): solar azimuth angle, degrees (1-dim for
              transect, 2-dim for sites and 3-dim for regular_grid)

        Raises:
            ValueError: if the inputs are not of the proper type or shape
        """

        # I need that self._times is timezone-aware and ndmin==1. To ensure
        # ndmin==1, I must convert to a numpy datetime64, but numpy complains
        # when dealing with timezone-aware times. Hence, first I use
        # pd.to_datetime to retain the timezone of the input times (if
        # timezone-aware) or assume UTC if the input times are naive. Then, I
        # convert to a numpy array of datetime64[ns] with ndmin=1, but de-localizing
        # the times to prevent numpy complains. Afterwards, I convert back to
        # pandas using again pd.to_datetime and localize
        _times = pd.to_datetime(times)

        if _times.tz is None:
            self._timezone = 'UTC'
        elif hasattr(_times.tz, 'zone'):
            self._timezone = _times.tz.zone
        else:
            self._timezone = 'UTC'

        self._times = pd.to_datetime(
            np.array(_times.tz_localize(None), ndmin=1, dtype='datetime64[ns]')
        ).tz_localize(self._timezone, ambiguous='infer')

        self._times_utc = pd.to_datetime(self._times).tz_convert('UTC')

        self._latitude = np.array(latitude, ndmin=1)
        self._longitude = np.array(longitude, ndmin=1)

        if self._latitude.ndim != self._longitude.ndim != 1:
            raise ValueError('expected 1-dim latitude and longitude, but '
                             f'got {self._latitude.ndim}-dim latitude and '
                             f'{self._longitude.ndim}-dim longitude')

        if usecase not in ('sites', 'regular_grid', 'transect'):
            raise ValueError(f'unknown use case `{usecase}`')
        self._usecase = usecase

        self._ecf = np.array(ecf)
        self._eot = np.array(eot)
        self._declination = np.array(declination)
        self._zenith = np.array(zenith)
        self._azimuth = np.array(azimuth)

        n_times = len(self._times)
        n_lats = self._latitude.size
        n_lons = self._longitude.size

        if self._usecase == 'sites':
            # sanity check
            assert n_lats == n_lons
            assert self._ecf.ndim == self._eot.ndim == self._declination.ndim == 1
            assert self._ecf.shape == self._eot.shape == self._declination.shape == (n_times,)
            assert self._zenith.ndim == self._azimuth.ndim == 2
            assert self._zenith.shape == self._azimuth.shape == (n_times, n_lats)
            # xarray coordinates
            coords = {
                'time': np.array(self._times.tz_localize(None), dtype='datetime64[ns]'),
                'location': range(n_lats),
                'latitude': ('location', self._latitude),
                'longitude': ('location', self._longitude)}
            # xarray DataArrays
            kwargs = {'coords': {k: coords[k] for k in ['location']}, 'dims': ['location']}
            self._latitude = xr.DataArray(self._latitude, name='latitude', **kwargs)
            self._longitude = xr.DataArray(self._longitude, name='longitude', **kwargs)
            kwargs = {'coords': {k: coords[k] for k in ['time']}, 'dims': ['time']}
            self._ecf = xr.DataArray(self._ecf, name='ecf', **kwargs)
            self._eot = xr.DataArray(self._eot, name='eot', **kwargs)
            self._declination = xr.DataArray(self._declination, name='declination', **kwargs)
            kwargs = {'coords': coords, 'dims': ['time', 'location']}
            self._zenith = xr.DataArray(self._zenith, name='zenith', **kwargs)
            self._azimuth = xr.DataArray(self._azimuth, name='azimuth', **kwargs)

        if self._usecase == 'regular_grid':
            # sanity check
            assert self._ecf.ndim == self._eot.ndim == self._declination.ndim == 1
            assert self._ecf.shape == self._eot.shape == self._declination.shape == (n_times,)
            assert self._zenith.ndim == self._azimuth.ndim == 3
            assert self._zenith.shape == self._azimuth.shape == (n_times, n_lats, n_lons)
            # xarray coordinates
            coords = {
                'time': np.array(self._times.tz_localize(None), dtype='datetime64[ns]'),
                'latitude': self._latitude, 'longitude': self._longitude}
            # xarray DataArrays
            kwargs = {'coords': {k: coords[k] for k in ['latitude']}, 'dims': ['latitude']}
            self._latitude = xr.DataArray(self._latitude, name='latitude', **kwargs)
            kwargs = {'coords': {k: coords[k] for k in ['longitude']}, 'dims': ['longitude']}
            self._longitude = xr.DataArray(self._longitude, name='longitude', **kwargs)
            kwargs = {'coords': {k: coords[k] for k in ['time']}, 'dims': ['time']}
            self._ecf = xr.DataArray(self._ecf, name='ecf', **kwargs)
            self._eot = xr.DataArray(self._eot, name='eot', **kwargs)
            self._declination = xr.DataArray(self._declination, name='declination', **kwargs)
            kwargs = {'coords': coords, 'dims': ['time', 'latitude', 'longitude']}
            self._zenith = xr.DataArray(self._zenith, name='zenith', **kwargs)
            self._azimuth = xr.DataArray(self._azimuth, name='azimuth', **kwargs)

        if self._usecase == 'transect':
            # sanity check
            assert n_times == n_lats == n_lons
            assert self._ecf.ndim == self._eot.ndim == self._declination.ndim == 1
            assert self._ecf.shape == self._eot.shape == self._declination.shape == (n_times,)
            assert self._zenith.ndim == self._azimuth.ndim == 1
            assert self._zenith.shape == self._azimuth.shape == (n_times,)
            # xarray coordinates
            coords = {
                'time': np.array(self._times.tz_localize(None), dtype='datetime64[ns]'),
                'latitude': ('time', self._latitude),
                'longitude': ('time', self._longitude)}
            # xarray DataArrays
            kwargs = {'coords': coords, 'dims': ['time']}
            self._latitude = xr.DataArray(self._latitude, name='latitude', **kwargs)
            self._longitude = xr.DataArray(self._longitude, name='longitude', **kwargs)
            self._ecf = xr.DataArray(self._ecf, name='ecf', **kwargs)
            self._eot = xr.DataArray(self._eot, name='eot', **kwargs)
            self._declination = xr.DataArray(self._declination, name='declination', **kwargs)
            self._zenith = xr.DataArray(self._zenith, name='zenith', **kwargs)
            self._azimuth = xr.DataArray(self._azimuth, name='azimuth', **kwargs)

        self._algorithm = algorithm
        self._engine = engine
        self._refraction = refraction

    @property
    def latitude(self):
        """DataArray: Latitudes where solar geometry is evaluated, degrees."""
        return self._latitude.rename('latitude').assign_attrs(
            units='degrees north')

    @property
    def longitude(self):
        """DataArray: Longitudes where solar geometry is evaluated, degrees."""
        return self._longitude.rename('longitude').assign_attrs(
            units='degrees east')

    @property
    def times(self):
        """Array of datetime64: Times where solar geometry is evaluated."""
        return self._times   # np.array(self._times, dtype='datetime64[ns]')

    @property
    def times_utc(self):
        """Array of datetime64: Universal coordinated times where solar geometry is evaluated."""
        return self._times_utc  # np.array(self._times_utc, dtype='datetime64[ns]')

    @property
    def algorithm(self):
        """str: Solar position algorithm."""
        return self._algorithm

    @property
    def engine(self):
        """str: calculation engine."""
        return self._engine

    @property
    def has_refraction(self):
        """bool: whether solar zenith angle consideres atmospheric refraction."""
        return self._refraction

    @property
    def usecase(self):
        """str: usage case."""
        return self._usecase

    @property
    def timezone(self):
        """str: times time zone."""
        return self._timezone

    @property
    def true_solar_time(self):
        """DataArray: True solar time, known also to as local apparent time."""
        dt64 = np.datetime64(1, 'ns')

        utc_f = xr.DataArray(
            np.array(self._times_utc.tz_localize(None), dtype=dt64).astype('float64'),
            coords={'time': np.array(self._times.tz_localize(None), dtype=dt64)},
            dims=['time'])

        tst_f = utc_f + (self.eot + 4*self.longitude) * 60 * 1e9  # nano-seconds
        return tst_f.astype(dt64).rename('true_solar_time').assign_attrs(
            description='True solar time (a.k.a. local apparent time, LAT)')

    @property
    def local_standard_time(self):
        """Array of datetime64: Times where solar geometry is evaluated."""
        return self.times

    @property
    def ecf(self):
        """DataArray: sun-earth orbit's eccentricity correction factor."""
        return self._ecf.rename('ecf').assign_attrs(
            units='-',
            description='sun-earth distance correction factor'
        )

    @property
    def eot(self):
        """DataArray: equation of time, in minutes."""
        return self._eot.rename('eot').assign_attrs(
            units='minutes',
            description='equation of time')

    @property
    def declination(self):
        """DataArray: solar declination, in degrees."""
        return np.degrees(self._declination).rename('declination').assign_attrs(
            units='degrees',
            description='solar declination')

    @property
    def dec(self):
        """DataArray: solar declination, in degrees. Alias for `declination`."""
        return self.declination

    @property
    def zenith(self):
        """DataArray: solar zenith angle, in degrees [0, 180]."""
        return np.degrees(self._zenith).rename('zenith').assign_attrs(
            units='degrees',
            range='[0, 180]',
            description='solar zenith angle')

    @property
    def sza(self):
        """DataArray: solar zenith angle, in degrees [0, 180]. Alias for `zenith`."""
        return self.zenith.rename('sza')

    @property
    def elevation(self):
        """DataArray: solar elevation angle, in degrees [-90, 90]."""
        return (90. - self.zenith).rename('elevation').assign_attrs(
            units='degrees',
            range='[-90, 90]',
            description='solar elevation angle')

    @property
    def azimuth(self):
        """DataArray: solar azimuth angle, in degrees [-180 am, 180 pm]."""
        return np.degrees(self._azimuth).rename('azimuth').assign_attrs(
            units='degrees',
            range='(-180 am, 180 pm)',
            description='solar azimuth angle'
        )

    @property
    def saa(self):
        """DataArray: solar azimuth angle, in degrees [-180 am, 180 pm]. Alias for `azimuth`."""
        return self.azimuth.rename('saa')

    @property
    def cosz(self):
        """DataArray: cosine of solar zenith angle."""
        return np.cos(self._zenith).rename('cosz').assign_attrs(
            units='-',
            range='[-1, 1]',
            description='solar zenith angle\'s cosine')

    def eth(self, ISC=1361.1, am_correction='none'):
        """Extraterrestrial horizontal solar irradiance.

        Calculates the extraterrestrial horizontal solar irradiance, in W m^-2,
        with optional corrections to account for solar eclipse obscuration
        and for the effect of high optical airmass if eth is used to evaluate
        clearness index

        Args:
            ISC (float): Solar constant, in W m<supersript>-2</superscript>
            am_correction (str): whether use airmass correction. Useful if eth is
              intended to evaluate the clearness index. Allowed values are `none`,
              `perez_1990` [1], `skarveith_1998` [2] and
              `gonzalez_and_calbo_1999` [3]

        Returns:
            A xarray's DataArray

        Raises:
            ValueError: if am_correction is unknown

        References:

            [1] Perez et al. (1990) Making full use of the clearness index for
            parameterizing hourly insolation conditions. Sol Energy, Vol. 45(2),
            pp 111-114. doi 10.1016/0038-092X(90)90036-C

            [2] Skartveit et al. (1998) An hourly diffuse fraction model with
            correction for variability and surface albedo. Sol Energy, Vol. 63(3),
            pp. 173-183. doi 10.1016/S0038-092X(98)00067-X

            [3] González and Calbó (1999) Influence of the global radiation
            variability on the hourly diffuse fraction correlations. Sol Energy,
            Vol. 65(2), pp. 119-131. doi 10.1016/S0038-092X(98)00121-2
        """

        eth = ISC*self.ecf*np.maximum(0., self.cosz)

        if am_correction == 'none':
            corr = 1.
        elif am_correction == 'perez_1990':
            am = self.airmass('kasten_and_young_1989')
            corr = 0.1 + 1.031*np.exp(-1.4 / (0.9 + (9.4 / am)))
        elif am_correction == 'skarveith_1998':
            corr = 0.83 - 0.56*np.exp(-0.06*(90.-self.sza))
        elif am_correction == 'gonzalez_and_calbo_1999':
            am = self.airmass('gueymard_2003')
            dCDA = 0.124 - 0.0285*np.log(am)
            corr = 0.229 + 0.957*np.exp(-1.74*am*dCDA)
        else:
            raise ValueError(f'unknown am_correction {am_correction}')

        return (eth*corr).rename('eth').assign_attrs(
            units='W m-2',
            description='extraterrestrial horizontal solar irradiance')

    def airmass(self, parameterization='gueymard_2003'):
        """Atmosphere relative optical air mass.

        Calculates the relative optical air mass of the atmosphere according
        to various parameterizations

        Args:
            parameterization (str): `kasten_1965` [1], `kasten_and_young_1989` [2],
            `gueymard_2001` [3] and `gueymard_2003` [4]

        Returns:
            A xarray's DataArray

        Raises:
            ValueError: if the parameterization is unknown

        References:

            [1] Table V in Kasten, F. (1965) A new table and approximation
            formula for the relative optical air mass. CRREL Tech. Report 136

            [2] Kasten, F. and Young, A.T. (1989) Revised optical air mass
            tables and approximation formula. Applied Optics. Vol. 28(22),
            pp. 4735-4738. doi 10.1364/AO.28.004735

            [3] Table A.1 in Gueymard, C.A. (2001) Parameterized transmittance
            model for direct bean and circumsolar spectral irradiance. Sol
            Energy, Vol. 71(5), pp. 325-346. doi 10.1016/S0038-092X(01)00054-8

            [4] Eq. B.8 in Gueymard, C.A. (2003) Direct solar transmittance
            and irradiance predictions with broadband models. Part I - Detailed
            theoretical performance assessment. Sol Energy, Vol. 74, pp. 355-379
            doi 10.1016/S0038-092X(03)00195-6
        """

        if parameterization == 'kasten_1965':
            Da = np.maximum(1e-4, 93.885 - self.sza)
            rec_am = self.cosz + 0.15*Da**(-1.253)

        elif parameterization == 'kasten_and_young_1989':
            Da = np.maximum(1e-4, 96.07995 - self.sza)
            rec_am = self.cosz + 0.50572*Da**(-1.6364)

        elif parameterization == 'gueymard_2001':
            Da = np.maximum(1e-4, 96.4836 - self.sza)
            rec_am = self.cosz + 0.45665*(self.sza**0.07)*(Da**(-1.697))

        elif parameterization == 'gueymard_2003':
            Da = np.maximum(1e-4, 96.741 - self.sza)
            rec_am = self.cosz + 0.48353*(self.sza**0.095846)*(Da**(-1.754))

        else:
            raise ValueError(
                f'unknown airmass parameterization {parameterization}')

        am = 1./rec_am.where(self.sza <= 90., other=np.nan)
        return am.clip(1., np.inf).rename('airmass').assign_attrs(
            units='-',
            description='relative optical air mass'
        )

    def sunrise(self, units='deg'):
        """Sunrise.

        Calculates the sunrise angle or sunrise time

        Args:
            units (str): `deg` for sunrise angle in degrees, `rad` for
            radians, `tst` for true solar time sunrise, `utc` for UTC and
            `local` for local time sunrise

        Returns:
            A xarray's DataArray

        Raises:
            ValueError: if units are unknown

        References:

            [1] Eq. (1.5.4) in Iqbal, M., An Introduction to Solar Radiation,
            Academic Press, 1983.
        """
        assert units in ('rad', 'deg', 'tst', 'utc', 'local')

        tanlat = np.tan(np.radians(self.latitude))
        tandec = np.tan(np.radians(self.declination))
        tantan = -tandec*tanlat

        domain = (tantan >= -1) & (tantan <= 1)
        wsr = np.arccos(tantan.where(domain, other=np.nan))  # radians

        if units == 'deg':
            wsr = np.degrees(wsr)

        if units in ('tst', 'utc', 'local'):
            wsr = 12. - (np.degrees(wsr)/15.)
            delta_ns = wsr.copy(
                data=(wsr*3600*1e9).astype('timedelta64[ns]'))
            tst_d = self.true_solar_time.copy(
                # if I convert directly from self.true_solar_time (which is
                # a DataArray) to `datetime64[D]` xarray complains. It raises
                # a UserWarning because I am downgrading the resolution from
                # ns to D. Hence, I use `to_numpy()` to convert to a numpy
                # array, then cast to `datetime64[D]`, then to `datetime64[ns]`
                # and the I create the DataArray tst_d. All this is only to
                # prevent that annoying UserWarning !!
                data=(self.true_solar_time.to_numpy()
                      .astype('datetime64[D]')
                      .astype('datetime64[ns]'))
            )
            wsr = delta_ns + tst_d

            if units in ('utc', 'local'):
                delta_hr = self.eot + 4*self.longitude
                delta_ns = delta_hr.copy(
                    data=(delta_hr*60*1e9).astype('timedelta64[ns]'))
                wsr = wsr - delta_ns

                if units == 'local':
                    utcoffset = xr.DataArray(
                        self._times.tz_localize(None) - self._times_utc.tz_localize(None),
                        coords={'time': wsr.coords['time']}, dims=('time',))
                    wsr = wsr + utcoffset

        return wsr.rename('wsr').assign_attrs(
            units={'deg': 'degrees', 'rad': 'radians', 'tst': 'TST',
                   'utc': 'UTC', 'local': self.timezone}.get(units, units),
            description='sunrise'
        )

    def sunset(self, units='deg'):
        """Sunset.

        Calculates the sunset angle or sunset time

        Args:
            units (str): `deg` for sunset angle in degrees, `rad` for
            radians, `tst` for true solar time sunset, `utc` for UTC and
            `local` for local time sunset

        Returns:
            A xarray's DataArray

        Raises:
            ValueError: if units are unknown

        References:

            [1] Eq. (1.5.4) in Iqbal, M., An Introduction to Solar Radiation,
            Academic Press, 1983.
        """
        assert units in ('rad', 'deg', 'tst', 'utc', 'local')

        if units in ('tst', 'utc', 'local'):
            wsr = self.sunrise(units)
            dl_hr = self.daylight_length()
            wss = wsr + dl_hr.copy(data=(dl_hr*3600*1e9).astype('timedelta64[ns]'))

        if units == 'deg':
            wss = (12. - wss)*15

        if units == 'rad':
            wss = np.radians((12. - wss)*15)

        return wss.rename('wss').assign_attrs(
            units={'deg': 'degrees', 'rad': 'radians', 'tst': 'TST',
                   'utc': 'UTC', 'local': self.timezone}.get(units, units),
            description='sunset'
        )

    def daylight_length(self):
        """Daylight length, in hours.

        References:

            [1] Eq. (1.5.5) in Iqbal, M., An Introduction to Solar Radiation,
            Academic Press, 1983.
        """
        return (2./15.)*self.sunrise(units='deg').rename('daylight_length').assign_attrs(
            units='hour',
            description='daylight length'
        )

    def incidence(self, sfc_slope, sfc_azimuth):
        """Cosine of the solar incidence angle.

        Calculates the cosine of the solar incidence angle.

        Args:
            sfc_slope (float): surface slope, in degrees
            sfc_azimuth (float): surface orientation, in degrees east positive

        References:

            [1] Eq. (1.6.5a) in Iqbal, M., An Introduction to Solar Radiation,
            Academic Press, 1983.
        """
        # hour angle...
        tst = self.true_solar_time.to_numpy()
        hour = self.true_solar_time.copy(
            data=tst - tst.astype('datetime64[D]')
        ).astype('float64') * 1e-9 / 3600
        ha = hour.copy(data=15*(12-hour))

        fcirc = lambda a: (np.cos(np.radians(a)), np.sin(np.radians(a)))  # noqa: E731

        cosbeta, sinbeta = fcirc(sfc_slope)
        cosgamma, singamma = fcirc(sfc_azimuth)
        coslat, sinlat = fcirc(self.latitude)
        cosdec, sindec = fcirc(self.dec)
        coshour, sinhour = fcirc(ha)
        coss = ((sinlat*cosbeta - coslat*sinbeta*cosgamma)*sindec +
                (coslat*cosbeta + sinlat*sinbeta*cosgamma)*cosdec*coshour +
                cosdec*sinhour*sinbeta*singamma)
        return coss.rename('incidence').assign_attrs(
            units='-', description='cosine of the angle of incidence')
