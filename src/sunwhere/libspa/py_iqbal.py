
import numpy as np


DTYPE_TIME = np.int32
DTYPE_DATA = np.float64

PI = np.pi
TWOPI = 2.*PI
HALFPI = 0.5*PI
RAD = PI / 180.


def sunpos(unixtime, longitude, latitude, ndim=1, with_refraction=False):

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


def day_angle(year, month, day, std_time):
    FIRST_DAY_OF_MONTH = np.array([0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365])
    # day of year (1 for 1st January) and day angle
    doy = day + FIRST_DAY_OF_MONTH[month - 1]
    bisextile = (((year % 4) == 0) & ((year % 100) != 0)) | ((year % 400) == 0)
    doy = np.where(bisextile & (month > 2.), doy + 1, doy)
    return TWOPI * (doy + (std_time/24.)) / 365.


def eccentricity_correction_factor(dayangl):
    # # eccentricity correction factor, Eq (1.2.1), Iqbal
    # return 1.000110
    #      + 0.034221*cos(dayangl) + 0.001280*sin(dayangl)
    #      + 0.000719*cos(2*dayangl) + 0.000077*sin(2*dayangl);

    # refit to match NREL's ecf (Sep 2021)
    return (1.00012825e+00
            + 3.33572253e-02*np.cos(dayangl) + 1.99417932e-03*np.sin(dayangl)
            + 6.97434604e-04*np.cos(2*dayangl) + 1.14828177e-04*np.sin(2*dayangl)
            + 1.74028675e-05*np.cos(3*dayangl) + 1.13778023e-05*np.sin(3*dayangl))


def equation_of_time(dayangl):
    # # equation of time, hours, Eq (1.4.1), Iqbal
    # return (0.000075
    #       + 0.001868*cos(dayangl) - 0.032077*sin(dayangl)
    #       - 0.014615*cos(2*dayangl) - 0.04089*sin(2*dayangl)
    #       )*229.18 / 60.;

    # equation of time, hours, refit to match NREL's eot
    return (0.00986571
            + 0.58688718*np.cos(dayangl) - 7.34538133*np.sin(dayangl)
            - 3.31493999*np.cos(2*dayangl) - 9.35366541*np.sin(2*dayangl)
            - 0.08151750*np.cos(3*dayangl) - 0.30892409*np.sin(3*dayangl)
            - 0.13532889*np.cos(4*dayangl) - 0.17336220*np.sin(4*dayangl))/60.


def solar_declination(dayangl):
    # sun declination, radians, Eq (1.3.1), Iqbal
    return (0.006918
            - 0.399912*np.cos(dayangl) + 0.070257*np.sin(dayangl)
            - 0.006758*np.cos(2*dayangl) + 0.000907*np.sin(2*dayangl)
            - 0.002697*np.cos(3*dayangl) + 0.001480*np.sin(3*dayangl))


def _sunpos_one_dimensional(unixtime, longitude, latitude, with_refraction):

    times = unixtime.astype('datetime64[s]')
    Y, M, D, h, m, _ = [times.astype(f'M8[{t}]') for t in 'YMDhms']
    year = (1970+Y).astype(np.int32)
    month = ((M-Y)+1).astype(np.int32)
    day = ((D-M)+1).astype(np.int32)
    hour = ((times-D).astype('m8[h]')).astype(np.int32)
    minute = ((times-h).astype('m8[m]')).astype(np.int32)
    second = ((times-m).astype('m8[s]')).astype(np.int32)

    std_time = hour + (minute + second / 60.) / 60.  # time, hours
    dayangl = day_angle(year, month, day, std_time)  # day angle, radians
    ecf = eccentricity_correction_factor(dayangl)  # eccentricity correction factor
    eot = equation_of_time(dayangl)  # equation of time, hours
    declination = solar_declination(dayangl)  # radians

    Nr = len(year)

    sza_rad = np.full((Nr,), np.nan)
    azi_rad = np.full((Nr,), np.nan)

    sinDecl = np.sin(declination)
    cosDecl = np.cos(declination)
    tanDecl = np.tan(declination)

    lat_rad = latitude*RAD
    cosLat = np.cos(lat_rad)
    sinLat = np.sin(lat_rad)

    LAT = (std_time + longitude/15. + eot)  # hours, Eq (1.4.2)
    LAT = np.where(LAT < 0., LAT + 24., LAT)
    LAT = np.where(LAT >= 24., LAT - 24., LAT)
    hour_angle = (LAT - 12.)*15.*RAD
    cosHour = np.cos(hour_angle)

    sza_rad = np.arccos(sinDecl*sinLat + cosDecl*cosLat*cosHour)

    if with_refraction is True:
        # Simple refraction correction from Rigollier et al. (2000),
        # doi:10.1016/S0038-092X(99)00055-9
        alpha = HALFPI - sza_rad
        alpha2 = alpha*alpha
        refcor = (0.061359*(0.1594+(1.1230*alpha)
                  + (0.065656*alpha2))/(1.+(28.9344*alpha)
                  + (277.3971*alpha2)))
        sza_rad = HALFPI - (alpha + refcor)

    dy = -np.sin(hour_angle)
    dx = tanDecl*cosLat - sinLat*cosHour
    azi_rad = np.arctan2(dy, dx)
    azi_rad = np.where(azi_rad < 0., azi_rad + PI, azi_rad - PI)

    eot = eot * 60  # minutes

    return sza_rad, azi_rad, ecf, declination, eot


def _sunpos_two_dimensional(unixtime, longitude, latitude, with_refraction):

    times = unixtime.astype('datetime64[s]')
    Y, M, D, h, m, _ = [times.astype(f'M8[{t}]') for t in 'YMDhms']
    year = (1970+Y).astype(np.int32)
    month = ((M-Y)+1).astype(np.int32)
    day = ((D-M)+1).astype(np.int32)
    hour = ((times-D).astype('m8[h]')).astype(np.int32)
    minute = ((times-h).astype('m8[m]')).astype(np.int32)
    second = ((times-m).astype('m8[s]')).astype(np.int32)

    std_time = hour + (minute + second / 60.) / 60.  # time, hours
    dayangl = day_angle(year, month, day, std_time)  # day angle, radians
    ecf = eccentricity_correction_factor(dayangl)  # eccentricity correction factor
    eot = equation_of_time(dayangl)  # equation of time, hours
    declination = solar_declination(dayangl)  # radians

    Nt = len(year)
    Nr = len(latitude)

    sza_rad = np.full((Nt, Nr), np.nan)
    azi_rad = np.full((Nt, Nr), np.nan)

    sinDecl = np.sin(declination[:, None])
    cosDecl = np.cos(declination[:, None])
    tanDecl = np.tan(declination[:, None])

    lat_rad = latitude[None, :]*RAD
    cosLat = np.cos(lat_rad)
    sinLat = np.sin(lat_rad)

    LAT = (std_time[:, None] + longitude[None, :]/15. + eot[:, None])  # hours, Eq (1.4.2)
    LAT = np.where(LAT < 0., LAT + 24., LAT)
    LAT = np.where(LAT >= 24., LAT - 24., LAT)
    hour_angle = (LAT - 12.)*15.*RAD
    cosHour = np.cos(hour_angle)

    sza_rad = np.arccos(sinDecl*sinLat + cosDecl*cosLat*cosHour)

    if with_refraction is True:
        # Simple refraction correction from Rigollier et al. (2000),
        # doi:10.1016/S0038-092X(99)00055-9
        alpha = HALFPI - sza_rad
        alpha2 = alpha*alpha
        refcor = (0.061359*(0.1594+(1.1230*alpha)
                  + (0.065656*alpha2))/(1.+(28.9344*alpha)
                  + (277.3971*alpha2)))
        sza_rad = HALFPI - (alpha + refcor)

    dy = -np.sin(hour_angle)
    dx = tanDecl*cosLat - sinLat*cosHour
    azi_rad = np.arctan2(dy, dx)
    azi_rad = np.where(azi_rad < 0., azi_rad + PI, azi_rad - PI)

    eot = eot * 60  # minutes

    return sza_rad, azi_rad, ecf, declination, eot


def _sunpos_three_dimensional(unixtime, longitude, latitude, with_refraction):

    times = unixtime.astype('datetime64[s]')
    Y, M, D, h, m, _ = [times.astype(f'M8[{t}]') for t in 'YMDhms']
    year = (1970+Y).astype(np.int32)
    month = ((M-Y)+1).astype(np.int32)
    day = ((D-M)+1).astype(np.int32)
    hour = ((times-D).astype('m8[h]')).astype(np.int32)
    minute = ((times-h).astype('m8[m]')).astype(np.int32)
    second = ((times-m).astype('m8[s]')).astype(np.int32)

    std_time = hour + (minute + second / 60.) / 60.  # time, hours
    dayangl = day_angle(year, month, day, std_time)  # day angle, radians
    ecf = eccentricity_correction_factor(dayangl)  # eccentricity correction factor
    eot = equation_of_time(dayangl)  # equation of time, hours
    declination = solar_declination(dayangl)  # radians

    Nt = len(year)
    Ny = len(latitude)
    Nx = len(longitude)

    sza_rad = np.full((Nt, Ny, Nx), np.nan)
    azi_rad = np.full((Nt, Ny, Nx), np.nan)

    sinDecl = np.sin(declination[:, None, None])
    cosDecl = np.cos(declination[:, None, None])
    tanDecl = np.tan(declination[:, None, None])

    lat_rad = latitude[None, :, None]*RAD
    cosLat = np.cos(lat_rad)
    sinLat = np.sin(lat_rad)

    LAT = (std_time[:, None, None] + longitude[None, None, :]/15.
           + eot[:, None, None])  # hours, Eq (1.4.2)
    LAT = np.where(LAT < 0., LAT + 24., LAT)
    LAT = np.where(LAT >= 24., LAT - 24., LAT)
    hour_angle = (LAT - 12.)*15.*RAD
    cosHour = np.cos(hour_angle)

    sza_rad = np.arccos(sinDecl*sinLat + cosDecl*cosLat*cosHour)

    if with_refraction is True:
        # Simple refraction correction from Rigollier et al. (2000),
        # doi:10.1016/S0038-092X(99)00055-9
        alpha = HALFPI - sza_rad
        alpha2 = alpha*alpha
        refcor = (0.061359*(0.1594+(1.1230*alpha)
                  + (0.065656*alpha2))/(1.+(28.9344*alpha)
                  + (277.3971*alpha2)))
        sza_rad = HALFPI - (alpha + refcor)

    dy = -np.sin(hour_angle)
    dx = tanDecl*cosLat - sinLat*cosHour
    azi_rad = np.arctan2(dy, dx)
    azi_rad = np.where(azi_rad < 0., azi_rad + PI, azi_rad - PI)

    eot = eot * 60  # minutes

    return sza_rad, azi_rad, ecf, declination, eot
