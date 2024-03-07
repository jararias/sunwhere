
# flake8: noqa: F841
# pylint: disable=unused-variable,unused-argument,too-many-locals

import numpy as np
import numexpr as ne


PI = np.pi
TWOPI = 2.*PI
HALFPI = 0.5*PI
EMRAU = 0.000042587565907513806  # Earth mean radius, AU
DEG2RAD = np.radians(1)
RAD2DEG = np.degrees(1)


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

    jde = nex('unixtime / 86400 - 10957.5')  # julian days since 1 Jan 2000 UT
    decimal_hour = nex('((unixtime / 86400) % 1)*24')
    gmst = nex('6.697096103 + 6.570984737e-2*jde + decimal_hour')

    # ecliptic coordinates..
    omega = nex('2.267127827 - 9.300339267e-4*jde')
    L = nex('4.895036035 + 1.720279602e-2*jde')  # mean longitude, radians
    g = nex('6.239468336 + 1.720200135e-2*jde')  # mean anomaly, radians
    l = nex(  # ecliptic longitude, radians
        'L - 1.544353226e-4 + 3.338320972e-2*sin(g)'
        '+ 3.497596876e-4*sin(2*g) - 8.689729360e-6*sin(omega)')
    # obliquity of the ecliptic, radians
    ep = nex('4.090904909e-1 - 6.213605399e-9*jde + 4.418094944e-5*cos(omega)')

    # celestial coordinates (right ascension and declination)
    ra = nex('arctan2(cos(ep)*sin(l), cos(l))')  # right ascension, radians
    ra = nex('where(ra < 0, ra + TWOPI, ra)')
    delta = nex('arcsin(sin(ep)*sin(l))')  # declination, radians

    # equation of time, in minutes (Michalsky, 1988)
    dXsun = np.int32(np.abs(L)/TWOPI)
    eot = nex('4*(L + TWOPI*where(L < 0., 1. + dXsun, -dXsun) - ra)*RAD2DEG')
    eot = nex('where(eot > 20., eot - 1440., eot)')
    eot = nex('where(eot < -20., eot + 1440., eot)')

    # sun-earth actual distance in AU
    # R = nex('1.00014-0.01671*cos(g)-0.00014*cos(2.*g)')  # Michalsky, 1988
    # From solartrack: https://github.com/MarcvdSluys/SolTrack-Python
    #   It is more accurate than Michalsky, 1988
    #   Sluys and Kan, 2022 (https://arxiv.org/pdf/2209.01557.pdf)
    jdc = nex('jde / 36525.0')
    ecc = nex('0.016708634 - 0.000042037*jdc - 0.0000001267*jdc*jdc')
    ma = nex('6.240060141 + 628.301955152*jdc - 2.682571e-6*jdc*jdc')
    sec = nex('((3.34161088e-2 - 8.40725e-5*jdc - 2.443e-7*jdc*jdc)*sin(ma) +'
              '(3.489437e-4 - 1.76278e-6*jdc)*sin(2*ma))')
    R = nex('1.0000010178*(1.0 - ecc**2)/(1.0 + ecc*cos((ma+sec)))')

    # sun-earth distance (eccentricity) correction factor
    ecf = nex('1. / R**2')

    return delta, ra, gmst, eot, ecf


def _sunpos_one_dimensional(unixtime, longitude, latitude, with_refraction):

    nex = ne.evaluate

    delta, ra, gmst, eot, ecf = time_dependent_calculations(unixtime)

    hour_angle = nex('(gmst*15 + longitude)*DEG2RAD - ra')  # hour angle, radians
    sinhour = nex('sin(hour_angle)')
    coshour = nex('cos(hour_angle)')

    sinlat = nex('sin(latitude*DEG2RAD)')
    coslat = nex('cos(latitude*DEG2RAD)')

    sindec = nex('sin(delta)')
    cosdec = nex('cos(delta)')
    tandec = nex('tan(delta)')

    zenith = nex('arccos(coslat*coshour*cosdec + sindec*sinlat)')
    zenith = nex('zenith + EMRAU*sin(zenith)')  # ..parallax Correction

    if with_refraction is True:
        # adapted from the NREL's SPA..
        sun_elev = HALFPI - zenith
        A = 1  #  assumed A=1: A = (pressure / 1010.) * (283 / (273+temperature))
        refr_corr = nex('A*1.02/(60*tan(sun_elev + 3.1376e-3 / (sun_elev + 8.9186e-2)))')
        zenith = nex('HALFPI - where(sun_elev >= 0.83337, sun_elev + refr_corr, sun_elev)')

    azimuth = nex('arctan2(-sinhour, tandec*coslat - sinlat*coshour)')
    azimuth = nex('where(azimuth < 0., azimuth + PI, azimuth - PI)')

    # zenith, in radians; azimuth, in radians (0 south);
    # ecf, in AU; delta, in radians; eot, in minutes
    return zenith, azimuth, ecf, delta, eot


def _sunpos_two_dimensional(unixtime, longitude, latitude, with_refraction):

    nex = ne.evaluate

    delta, ra, gmst, eot, ecf = time_dependent_calculations(unixtime)

    lat = latitude[None, :]
    lon = longitude[None, :]

    # local coordinates
    ra = ra[:, None]
    gmst = gmst[:, None]
    hour_angle = nex('(gmst*15 + lon)*DEG2RAD - ra')  # hour angle, radians

    sinlat = nex('sin(lat*DEG2RAD)')
    coslat = nex('cos(lat*DEG2RAD)')
    sinhour = nex('sin(hour_angle)')
    coshour = nex('cos(hour_angle)')

    decli = delta[:, None]
    sindec = nex('sin(decli)')
    cosdec = nex('cos(decli)')
    tandec = nex('tan(decli)')

    zenith = nex('arccos(coslat*coshour*cosdec + sindec*sinlat)')
    zenith = nex('zenith + EMRAU*sin(zenith)')  # ..parallax Correction

    if with_refraction is True:
        # adapted from the NREL's SPA..
        sun_elev = HALFPI - zenith
        A = 1  #  assumed A=1: A = (pressure / 1010.) * (283 / (273+temperature))
        refr_corr = nex('A*1.7802e-2/(60*tan(sun_elev + 3.1376e-3 / (sun_elev + 8.9186e-2)))')
        zenith = nex('HALFPI - where(sun_elev >= -1.4545e-2, sun_elev + refr_corr, sun_elev)')

    azimuth = nex('arctan2(-sinhour, tandec*coslat - sinlat*coshour)')
    azimuth = nex('where(azimuth < 0., azimuth + PI, azimuth - PI)')

    # zenith, in radians; azimuth, in radians (0 south);
    # ecf, in AU; delta, in radians; eot, in minutes
    return zenith, azimuth, ecf, delta, eot


def _sunpos_three_dimensional(unixtime, longitude, latitude, with_refraction):

    nex = ne.evaluate

    delta, ra, gmst, eot, ecf = time_dependent_calculations(unixtime)

    lat = latitude[None, :, None]
    lon = longitude[None, None, :]

    # local coordinates
    ra = ra[:, None, None]
    gmst = gmst[:, None, None]
    hour_angle = nex('(gmst*15 + lon)*DEG2RAD - ra')  # hour angle, radians

    sinlat = nex('sin(lat*DEG2RAD)')
    coslat = nex('cos(lat*DEG2RAD)')
    sinhour = nex('sin(hour_angle)')
    coshour = nex('cos(hour_angle)')

    decli = delta[:, None, None]
    sindec = nex('sin(decli)')
    cosdec = nex('cos(decli)')
    tandec = nex('tan(decli)')

    zenith = nex('arccos(coslat*coshour*cosdec + sindec*sinlat)')
    zenith = nex('zenith + EMRAU*sin(zenith)')  # ..parallax Correction

    if with_refraction is True:
        # adapted from the NREL's SPA..
        sun_elev = HALFPI - zenith
        A = 1  #  assumed A=1: A = (pressure / 1010.) * (283 / (273+temperature))
        refr_corr = nex('A*1.02/(60*tan(sun_elev + 3.1376e-3 / (sun_elev + 8.9186e-2)))')
        zenith = nex('HALFPI - where(sun_elev >= 0.83337, sun_elev + refr_corr, sun_elev)')

    azimuth = nex('arctan2(-sinhour, tandec*coslat - sinlat*coshour)')
    azimuth = nex('where(azimuth < 0., azimuth + PI, azimuth - PI)')

    # zenith, in radians; azimuth, in radians (0 south);
    # ecf, in AU; delta, in radians; eot, in minutes
    return zenith, azimuth, ecf, delta, eot
