
# pylint: disable=too-many-locals,too-many-statements

import numpy as np


PI = np.pi
TWOPI = 2.*PI
HALFPI = 0.5*PI
EMRAU = 0.000042587565907513806


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
    jde = unixtime / 86400 - 10957.5  # julian days since 1 Jan 2000 UT
    decimal_hour = ((unixtime / 86400) % 1)*24

    # ecliptic coordinates..
    omega = 2.267127827 - 9.300339267e-4*jde  #
    L = 4.895036035 + 1.720279602e-2*jde  # mean longitude, radians
    g = 6.239468336 + 1.720200135e-2*jde  # mean anomaly, radians
    eclon = (  # ecliptic longitude, radians
        L - 1.544353226e-4 + 3.338320972e-2*np.sin(g)
        + 3.497596876e-4*np.sin(2*g) - 8.689729360e-6*np.sin(omega))
    # obliquity of the ecliptic, radians
    ep = 4.090904909e-1 - 6.213605399e-9*jde + 4.418094944e-5*np.cos(omega)

    sinl = np.sin(eclon)

    # celestial coordinates (right ascension and declination)
    ra = np.arctan2(np.cos(ep)*sinl, np.cos(eclon))  # right ascension
    ra = np.where(ra < 0., ra + TWOPI, ra)
    delta = np.arcsin(np.sin(ep)*sinl)  # declination, radians

    gmst = 6.697096103 + 6.570984737e-2*jde + decimal_hour

    # sun-earth actual distance in AU
    # R = 1.00014-0.01671*np.cos(g)-0.00014*np.cos(2.*g)  # Michalsky, 1988
    # From solartrack: https://github.com/MarcvdSluys/SolTrack-Python
    #   It is more accurate than Michalsky, 1988
    #   Sluys and Kan, 2022 (https://arxiv.org/pdf/2209.01557.pdf)
    jdc = jde / 36525.0
    ecc = 0.016708634 - 0.000042037*jdc - 0.0000001267*jdc*jdc
    ma = 6.240060141 + 628.301955152*jdc - 2.682571e-6*jdc*jdc
    sec = ((3.34161088e-2 - 8.40725e-5*jdc - 2.443e-7*jdc*jdc)*np.sin(ma) +
           (3.489437e-4 - 1.76278e-6*jdc)*np.sin(2*ma))
    R = 1.0000010178*(1.0 - ecc**2)/(1.0 + ecc*np.cos((ma+sec)))

    ecf = 1. / R**2  # sun-earth eccentricity correction factor

    # equation of time, in minutes (Michalsky, 1988)
    dXsun = np.int32(np.abs(L)/TWOPI)
    eot = 4. * np.degrees(L + TWOPI*np.where(L < 0., 1. + dXsun, -dXsun) - ra)
    eot = np.where(eot > 20., eot - 1440., eot)
    eot = np.where(eot < -20., eot + 1440., eot)

    return gmst, ra, delta, ecf, eot


def _sunpos_one_dimensional(unixtime, longitude, latitude, with_refraction):

    gmst, ra, delta, ecf, eot = time_dependent_calculations(unixtime)

    Nr = len(unixtime)
    zenith = np.zeros((Nr,))  # solar zenith angle, radians
    azimuth = np.zeros((Nr,))  # solar azimuth angle, radians

    # local coordinates
    lmst = gmst*15 + longitude  # local mean sidereal time, degrees
    hour_angle = np.radians(lmst) - ra  # hour angle, radians

    latrad = np.radians(latitude)
    coslat = np.cos(latrad)
    sinlat = np.sin(latrad)

    sinhour = np.sin(hour_angle)
    coshour = np.cos(hour_angle)

    sindec = np.sin(delta)
    cosdec = np.cos(delta)
    tandec = np.tan(delta)

    acos_zenith = coslat*coshour*cosdec + sindec*sinlat
    zenith = np.arccos(np.round(acos_zenith, 6))  # radians
    zenith = zenith + EMRAU * np.sin(zenith)  # parallax correction

    if with_refraction is True:
        # adapted from the NREL's SPA..
        sun_elev = HALFPI - zenith
        A = 1  # assumed A=1: A = (pressure / 1010.) * (283 / (273+temperature))
        refr_corr = A * 1.7802e-2/(60*np.tan(sun_elev + 3.1376e-3 / (sun_elev + 8.9186e-2)))
        zenith = HALFPI - np.where(sun_elev >= -1.4545e-2, sun_elev + refr_corr, sun_elev)

    azimuth = np.arctan2(-sinhour, tandec*coslat - sinlat*coshour)
    azimuth = np.where(azimuth < 0., azimuth + PI, azimuth - PI)

    # zenith, in radians; azimuth, in radians (0 south);
    # ecf, in AU; delta, in radians; eot, in minutes
    return zenith, azimuth, ecf, delta, eot


def _sunpos_two_dimensional(unixtime, longitude, latitude, with_refraction):

    gmst, ra, delta, ecf, eot = time_dependent_calculations(unixtime)

    Nt = len(unixtime)
    Nr = len(latitude)

    zenith = np.zeros((Nt, Nr))  # solar zenith angle, radians
    azimuth = np.zeros((Nt, Nr))  # solar azimuth angle, radians

    # local coordinates
    lmst = gmst[:, None]*15 + longitude[None, :]  # local mean sidereal time, degrees
    hour_angle = np.radians(lmst) - ra[:, None]  # hour angle, radians

    latrad = np.radians(latitude[None, :])
    sinlat = np.sin(latrad)
    coslat = np.cos(latrad)

    sinhour = np.sin(hour_angle)
    coshour = np.cos(hour_angle)

    sindec = np.sin(delta[:, None])
    cosdec = np.cos(delta[:, None])
    tandec = np.tan(delta[:, None])

    acos_zenith = coslat*coshour*cosdec + sindec*sinlat
    zenith = np.arccos(np.round(acos_zenith, 6))  # radians
    zenith = zenith + EMRAU * np.sin(zenith)  # parallax correction

    if with_refraction is True:
        # adapted from the NREL's SPA..
        sun_elev = HALFPI - zenith
        A = 1  # assumed A=1: A = (pressure / 1010.) * (283 / (273+temperature))
        refr_corr = A * 1.7802e-2/(60*np.tan(sun_elev + 3.1376e-3 / (sun_elev + 8.9186e-2)))
        zenith = HALFPI - np.where(sun_elev >= -1.4545e-2, sun_elev + refr_corr, sun_elev)

    azimuth = np.arctan2(-sinhour, tandec*coslat - sinlat*coshour)
    azimuth = np.where(azimuth < 0., azimuth + PI, azimuth - PI)

    # zenith, in radians; azimuth, in radians (0 south);
    # ecf, in AU; delta, in radians; eot, in minutes
    return zenith, azimuth, ecf, delta, eot


def _sunpos_three_dimensional(unixtime, longitude, latitude, with_refraction):

    gmst, ra, delta, ecf, eot = time_dependent_calculations(unixtime)

    Nt = len(unixtime)
    Ny = len(latitude)
    Nx = len(longitude)

    zenith = np.full((Nt, Ny, Nx), np.nan)  # solar zenith angle, radians
    azimuth = np.full((Nt, Ny, Nx), np.nan)  # solar azimuth angle, radians

    # local coordinates
    lmst = gmst[:, None, None]*15 + longitude[None, None, :]  # local mean sidereal time, degrees
    hour_angle = np.radians(lmst) - ra[:, None, None]  # hour angle, radians

    latrad = np.radians(latitude[None, :, None])
    sinlat = np.sin(latrad)
    coslat = np.cos(latrad)

    sinhour = np.sin(hour_angle)
    coshour = np.cos(hour_angle)

    sindec = np.sin(delta[:, None, None])
    cosdec = np.cos(delta[:, None, None])
    tandec = np.tan(delta[:, None, None])

    acos_zenith = coslat*coshour*cosdec + sindec*sinlat
    zenith = np.arccos(np.round(acos_zenith, 6))  # radians
    zenith = zenith + EMRAU * np.sin(zenith)  # parallax correction

    if with_refraction is True:
        # adapted from the NREL's SPA..
        sun_elev = HALFPI - zenith
        A = 1  # assumed A=1: A = (pressure / 1010.) * (283 / (273+temperature))
        refr_corr = A * 1.7802e-2/(60*np.tan(sun_elev + 3.1376e-3 / (sun_elev + 8.9186e-2)))
        zenith = HALFPI - np.where(sun_elev >= -1.4545e-2, sun_elev + refr_corr, sun_elev)

    azimuth = np.arctan2(-sinhour, tandec*coslat - sinlat*coshour)
    azimuth = np.where(azimuth < 0., azimuth + PI, azimuth - PI)

    # zenith, in radians; azimuth, in radians (0 south);
    # ecf, in AU; delta, in radians; eot, in minutes
    return zenith, azimuth, ecf, delta, eot
