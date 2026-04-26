
# flake8: noqa: F841
# pylint: disable=unused-variable,unused-argument,too-many-locals

import numpy as np
import numexpr as ne


PI = np.pi
TWOPI = 2.*PI
HALFPI = 0.5*PI
RAD2DEG = np.degrees(1)
DEG2RAD = np.radians(1)

# PRESSURE = 1010.0  # atmospheric pressure, hPa
# TEMPERATURE = 10.  # air temperature, degC
# ATM_REFR = (PRESSURE/1010.0) * (283./(273+TEMPERATURE))


def sunpos(unixtime, longitude, latitude, ndim=1, with_refraction=True):

    # unixtime: seconds since 1970-01-01T00

    ndim = int(ndim)
    assert ndim in (0, 1, 2), 'unexpected ndim %d' % ndim
    # =0: scalar computation (all calculations are performed cell by cell)
    # =1: time-optimized calculation: the time grid is common to all cells
    #     but the space-related calculations are peformed cell by cell
    # =2: regular grid: latitude is constant row-wise and longitude is
    #     constant column-wise. The time grid may be, or may not, common
    #     to all grid cells

    # inputs must be 1d arrays
    unixtime = np.array(unixtime, ndmin=1, dtype=np.float64).ravel()
    longitude = np.array(longitude, ndmin=1, dtype=np.float32).ravel()
    latitude = np.array(latitude, ndmin=1, dtype=np.float32).ravel()

    if ndim == 0:

        # -- one dimensional space: (n_locations,) == (n_times,)

        condition = unixtime.size == longitude.size == latitude.size
        assert condition, 'unixtime, latitude and longitude must have same shape'

        sza, saa, ecf, decl, eot = _sunpos_one_dimensional(
            unixtime, longitude, latitude, with_refraction)

    if ndim == 1:

        # -- two dimensional space: (n_times, n_locations)

        condition = longitude.size == latitude.size
        assert condition, 'longitude and latitude must have same shape'

        sza, saa, ecf, decl, eot = _sunpos_two_dimensional(
            unixtime, longitude, latitude, with_refraction)

    if ndim == 2:

        # -- three dimensional space: (n_times, n_lats, n_lons)

        sza, saa, ecf, decl, eot = _sunpos_three_dimensional(
            unixtime, longitude, latitude, with_refraction)

    return {'zenith': sza, 'azimuth': saa, 'ecf': ecf,
            'eot': eot, 'declination': decl}


def time_dependent_calculations(unixtime):

    nex = ne.evaluate

    tJD = nex('unixtime / 86400 - 10957.5')  # time in julian days since 1 Jan 2000 UT
    tJC = tJD / 36525.0                      # time in julian centuries since 1 Jan 2000 UT
    tJC2 = tJC**2

    # compute the ecliptic longitude of the Sun and the obliquity of the ecliptic and nutation

    l0 = nex('4.895063168 + 628.331966786*tJC + 5.291838e-6*tJC2')    # Mean longitude, Eq 2
    ma = nex('6.240060141 + 628.301955152*tJC - 2.682571e-6*tJC2')    # Mean anomaly, Eq 3

    sec = nex(
        '(3.34161088e-2 - 8.40725e-5*tJC - 2.443e-7*tJC2)'
        '*sin(ma) + (3.489437e-4 - 1.76278e-6*tJC)*sin(2*ma)')        # Sun's equation of the centre, Eq 4

    omg = nex('2.1824390725 - 33.7570464271*tJC + 3.622256e-5*tJC2')  # Lon. of Moon's mean ascending node, Eq 6
    dpsi = nex('-8.338601e-5*sin(omg)')                               # Nutation in longitude, Eq 7
    ecc = nex('0.016708634 - 0.000042037*tJC - 0.0000001267*tJC2')    # Eccentricity of the Earth's orbit, Eq 8
    R = nex('1.0000010178*(1.0 - ecc**2)/(1.0 + ecc*cos(ma+sec))')    # Geocentric distance of the Sun in AU, Eq 10

    L = nex(                                                          # Apparent geocentric longitude, Eq 13
        '((l0 + sec) +'                                               # True longitude, Eq 5
        '+ (-9.93087e-5/R) +'                                         # Aberration, Eq 12
        ' dpsi) % TWOPI'
    )

    obliquity = nex(                                                  # True obliquity of the ecliptic, Eq 17
        '(0.409092804222 - 2.26965525e-4*tJC - 2.86e-9*tJC2) +'       # Mean obliquity of the ecliptic, Eq 15    
        '(4.4615e-5*cos(omg))'                                        # Nutation in obliquity, Eq 16
    )

    # convert ecliptic coordinates to geocentric equatorial coordinates (right ascension and declination)

    sin_L = nex('sin(L)')
    cos_L = nex('cos(L)')
    cos_o = nex('cos(obliquity)')
    sin_o = nex('sin(obliquity)')
    ra_uncorr = nex('arctan2(cos_o*sin_L, cos_L) % TWOPI')   # Uncorrected right ascension, Eq 18, 0 <= azimuth < 2pi
    declination_uncorr = nex('arcsin(sin_o*sin_L)')          # Uncorrected declination, Eq 19

    # equation of time, minutes (from NREL's SPA)

    x_deg = nex('(l0 - ra_uncorr + dpsi*cos_o)*RAD2DEG')
    eot = nex('((-0.0057183 + x_deg) % 360)*4')
    eot = nex('where(eot < -20., eot + 1440, eot)')
    eot = nex('where(eot > 20., eot - 1440, eot)')

    # sun-earth distance correction factor

    ecf = nex('1 / R**2')

    # Apparent GMST is needed to convert from equatorial to horizontal coordinates
    gmst = nex('4.89496121 + 6.300388098985*tJD + 6.77e-6*tJC2')      # Greenwich mean sidereal time, Eq 22
    agst = nex('gmst + dpsi*cos_o')                                   # Apparent Greenwich sidereal time, from Eqs 21 and 23

    return agst, ra_uncorr, declination_uncorr, ecf, eot


def _sunpos_one_dimensional(unixtime, longitude, latitude, with_refraction):

    nex = ne.evaluate

    agst, ra_uncorr, declination_uncorr, ecf, eot = time_dependent_calculations(unixtime)

    ha = nex('agst - ra_uncorr + longitude*DEG2RAD')
    sin_ha = nex('sin(ha)')  # local hour angle, rad, Eqs 21 and 20
    cos_ha = nex('cos(ha)')

    sin_dec = nex('sin(declination_uncorr)')
    cos_dec = nex('cos(declination_uncorr)')
    tan_dec = nex('tan(declination_uncorr)')

    sin_lat = nex('sin(latitude*DEG2RAD)')
    cos_lat = nex('cos(latitude*DEG2RAD)')

    sin_alt = nex('sin_lat*sin_dec + cos_lat*cos_dec*cos_ha')  # -pi/2 <= solar altitude <= pi/2, Eq 24
    azimuth = nex('arctan2(sin_ha, cos_ha*sin_lat - tan_dec*cos_lat)')  # 0 <= solar azimth < 2pi, Eq 25

    # correct for parallax
    alt_uncorr = nex('arcsin(sin_alt)')
    alt = nex('alt_uncorr - 4.2635e-5*cos(alt_uncorr)')      # horizontal parallax = 8.794" = 4.2635e-5 rad

    # correct for atmospheric refraction (ignore T and P), Eq 27
    if with_refraction is True:
        alt = nex('alt + (2.967e-4 / tan(alt + 3.1376e-3/(alt + 8.92e-2)))')

    zenith = HALFPI - alt

    # zenith, in radians; azimuth, in radians (0 south);
    # ecf, in AU; delta, in radians; eot, in minutes
    return zenith, azimuth, ecf, declination_uncorr, eot


def _sunpos_two_dimensional(unixtime, longitude, latitude, with_refraction):

    nex = ne.evaluate

    agst, ra_uncorr, declination_uncorr, ecf, eot = time_dependent_calculations(unixtime)

    ha = (agst - ra_uncorr)[:, None] + (longitude*DEG2RAD)[None, :]  # local hour angle, rad, Eqs 21 and 20
    sin_ha = nex('sin(ha)')
    cos_ha = nex('cos(ha)')

    sin_dec = nex('sin(declination_uncorr)')[:, None]
    cos_dec = nex('cos(declination_uncorr)')[:, None]
    tan_dec = sin_dec / cos_dec

    sin_lat = nex('sin(latitude*DEG2RAD)')[None, :]
    cos_lat = nex('cos(latitude*DEG2RAD)')[None, :]

    sin_alt = nex('sin_lat*sin_dec + cos_lat*cos_dec*cos_ha')   # -pi/2 <= solar altitude <= pi/2, Eq 24
    azimuth = nex('arctan2(sin_ha, cos_ha*sin_lat - tan_dec*cos_lat)')  # 0 <= solar azimth < 2pi, Eq 25

    # correct for parallax
    alt_uncorr = nex('arcsin(sin_alt)')
    alt = nex('alt_uncorr - 4.2635e-5*cos(alt_uncorr)')  # horizontal parallax = 8.794" = 4.2635e-5 rad

    # correct for atmospheric refraction (ignore T and P), Eq 27
    if with_refraction is True:
        alt = nex('alt + (2.967e-4 / tan(alt + 3.1376e-3/(alt + 8.92e-2)))')

    zenith = HALFPI - alt

    # zenith, in radians; azimuth, in radians (0 south);
    # ecf, in AU; delta, in radians; eot, in minutes
    return zenith, azimuth, ecf, declination_uncorr, eot


def _sunpos_three_dimensional(unixtime, longitude, latitude, with_refraction):

    nex = ne.evaluate

    agst, ra_uncorr, declination_uncorr, ecf, eot = time_dependent_calculations(unixtime)

    ha = ((agst - ra_uncorr)[:, None, None] +
          (longitude*DEG2RAD)[None, None, :])  # local hour angle, rad, Eqs 21 and 20
    sin_ha = nex('sin(ha)')
    cos_ha = nex('cos(ha)')

    sin_dec = nex('sin(declination_uncorr)')[:, None, None]
    cos_dec = nex('cos(declination_uncorr)')[:, None, None]
    tan_dec = nex('tan(declination_uncorr)')[:, None, None]

    sin_lat = nex('sin(latitude*DEG2RAD)')[None, :, None]
    cos_lat = nex('cos(latitude*DEG2RAD)')[None, :, None]

    sin_alt = nex('sin_lat*sin_dec + cos_lat*cos_dec*cos_ha')   # -pi/2 <= solar altitude <= pi/2, Eq 24
    azimuth = nex('arctan2(sin_ha, cos_ha*sin_lat - tan_dec*cos_lat)')  # 0 <= solar azimth < 2pi, Eq 25

    # correct for parallax
    alt_uncorr = nex('arcsin(sin_alt)')
    alt = nex('alt_uncorr - 4.2635e-5*cos(alt_uncorr)')  # horizontal parallax = 8.794" = 4.2635e-5 rad

    # correct for atmospheric refraction (ignore T and P), Eq 27
    if with_refraction is True:
        alt = nex('alt + (2.967e-4 / tan(alt + 3.1376e-3/(alt + 8.92e-2)))')

    zenith = HALFPI - alt

    # zenith, in radians; azimuth, in radians (0 south);
    # ecf, in AU; delta, in radians; eot, in minutes
    return zenith, azimuth, ecf, declination_uncorr, eot
