
# flake8: noqa: F841
# pylint: disable=unused-variable,unused-argument,too-many-locals,too-many-statements

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#  NOTE: The `numexpr` reduction functions `sum` and `prod` are very slow.  #
#  It is much faster to perform such reductions directly in numpy !!!!      #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import numpy as np
import numexpr as ne


HALFPI = 0.5 * np.pi
DEG2RAD = np.radians(1)
RAD2DEG = np.degrees(1)
DELTA_T = 67.0

ELEVATION = 0.  # observer's elevation, meters
PRESSURE = 1013.25  # atmospheric pressure, hPa
TEMPERATURE = 12.  # air temperature, degC
ATM_REFR = (PRESSURE/1010.0) * (283./(273+TEMPERATURE))


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
    longitude = np.array(longitude, ndmin=1, dtype=np.float64).ravel()
    latitude = np.array(latitude, ndmin=1, dtype=np.float64).ravel()

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

    # julian day
    jd = nex('unixtime / 86400 + 2440587.5')

    # julian ephemeris day
    jde = nex('jd + DELTA_T / 86400')

    # julian century
    jc = nex('(jd - 2451545) / 36525')

    # julian ephemeris century
    jce = nex('(jde - 2451545) / 36525')

    # julian ephemeris millenium
    jme = nex('jce / 10')
    jme_ = jme[None, :]

    # heliocentric radius vector
    r0, r1, r2, r3, r4 = np.sum(nex('C_R_0*cos(C_R_1+C_R_2*jme_)'), axis=1)
    R = nex('(r0 + r1 * jme + r2 * jme**2 + r3 * jme**3 + r4 * jme**4)/10**8')

    # heliocentric longitude, degrees
    l0, l1, l2, l3, l4, l5 = np.sum(nex('C_LONG_0*cos(C_LONG_1+C_LONG_2*jme_)'), axis=1)
    L = nex('(RAD2DEG*(l0 + l1*jme + l2*jme**2 + l3*jme**3 + l4*jme**4 + l5*jme**5)/10**8) % 360')

    # heliocentric latitude, degrees
    b0, b1 = np.sum(nex('C_LAT_0*cos(C_LAT_1+C_LAT_2*jme_)'), axis=1)
    B = nex('RAD2DEG*(b0 + b1 * jme)/10**8')

    # geocentric longitude
    theta = nex('(L + 180.) % 360')

    # geocentric latitude, radians
    beta = nex('-DEG2RAD * B')

    # mean elongation
    x0 = nex('297.85036 + 445267.111480*jce - 0.0019142*jce**2 + jce**3 / 189474')
    x0_ = x0[None, :]

    # mean anomaly sun
    x1 = nex('357.52772 + 35999.050340*jce - 0.0001603*jce**2 - jce**3 / 300000')
    x1_ = x1[None, :]

    # mean anomaly moon
    x2 = nex('134.96298 + 477198.867398*jce + 0.0086972*jce**2 + jce**3 / 56250')
    x2_ = x2[None, :]

    # moon argument latitude
    x3 = nex('93.27191 + 483202.017538*jce - 0.0036825*jce**2 + jce**3 / 327270')
    x3_ = x3[None, :]

    # moon ascending longitude
    x4 = nex('125.04452 - 1934.136261*jce + 0.0020708*jce**2 + jce**3 / 450000')
    x4_ = x4[None, :]

    # longitude nutation
    delta_psi = np.sum(
        nex(
            '(C_NUT_A + C_NUT_B*jme_)*'
            'sin('
            '  (C_NUT_Y_0*x0_ + C_NUT_Y_1*x1_ + '
            '   C_NUT_Y_2*x2_ + C_NUT_Y_3*x3_ + '
            '   C_NUT_Y_4*x4_)*DEG2RAD'
            ')/36000000'
        ), axis=0
    )

    # obliquity nutation
    delta_eps = np.sum(
        nex(
            '(C_NUT_C + C_NUT_D*jme_)*'
            'cos('
            '  (C_NUT_Y_0*x0_ + C_NUT_Y_1*x1_ + '
            '   C_NUT_Y_2*x2_ + C_NUT_Y_3*x3_ + '
            '   C_NUT_Y_4*x4_)*DEG2RAD'
            ')/36000000'
        ), axis=0
    )

    # mean ecliptic obliquity
    epsilon0 = nex(
        '84381.448 - 4680.93*(jme/10.) - 1.55*(jme/10.)**2 + 1999.25*(jme/10.)**3 - '
        '51.38*(jme/10.)**4 - 249.67*(jme/10.)**5 - 39.05*(jme/10.)**6 + 7.12*(jme/10.)**7 + '
        '27.87*(jme/10.)**8 + 5.79*(jme/10.)**9 + 2.45*(jme/10.)**10')

    # true ecliptic obliquity, radians
    epsilon = nex('(epsilon0/3600. + delta_eps)*DEG2RAD')

    # aberration correction
    delta_tau = nex('-20.4898 / (3600. * R)')

    # apparent sun longitude, radians
    lamd = nex('(theta + delta_psi + delta_tau)*DEG2RAD')

    # mean sidereal time
    v0 = nex('(280.46061837 + 360.98564736629 * (jd - 2451545) + '
             '0.000387933*jc**2 - jc**3 / 38710000) % 360.')

    # apparent sidereal time
    v = nex('v0 + delta_psi * cos(epsilon)')

    # geocentric sun right ascension
    alpha = nex(
        'arctan2('
        '  sin(lamd)*cos(epsilon) - tan(beta)*sin(epsilon), cos(lamd)'
        ')*RAD2DEG % 360')

    # geocentric sun declination, radians
    delta = nex('arcsin(sin(beta)*cos(epsilon) + cos(beta)*sin(epsilon)*sin(lamd))')

    # sun mean longitude
    m = nex('280.4664567 + 360007.6982779*jme + 0.03032028*jme**2 + '
            'jme**3 / 49931 - jme**4 / 15300 - jme**5 / 2000000')

    # equation of time, minutes
    eot = nex('((m - 0.0057183 - alpha + delta_psi*cos(epsilon)) % 360)*4')
    eot = nex('where(eot < -20., eot + 1440, eot)')
    eot = nex('where(eot > 20., eot - 1440, eot)')

    return alpha, delta, eot, v, R


def _sunpos_one_dimensional(unixtime, longitude, latitude, with_refraction):

    # unixtime: seconds since 1970-01-01T00 UTC

    nex = ne.evaluate

    alpha, delta, eot, v, R = time_dependent_calculations(unixtime)

    # start of spatial calculations...

    latrad = latitude*DEG2RAD

    # local hour angle, radians
    H = nex('((v + longitude - alpha) % 360)*DEG2RAD')

    # equatorial horizontal parallax
    xi = nex('(8.794 / (3600 * R))*DEG2RAD')  # radians
    u = nex('arctan(0.99664719*tan(latrad))')
    x = nex('cos(u) + ELEVATION / 6378140. * cos(latrad)')
    y = nex('0.99664719*sin(u) + ELEVATION / 6378140. * sin(latrad)')

    # parallax sun right ascension, radians
    deltaa = nex('arctan2(-x*sin(xi)*sin(H), cos(delta) - x*sin(xi)*cos(H))')

    # topocentric sun declination, radians
    deltap = nex('arctan2((sin(delta)-y*sin(xi))*cos(deltaa), cos(delta)-x*sin(xi)*cos(H))')

    # topocentric local hour angle, radians
    Hp = nex('H - deltaa')

    # topocentric elevation angle without atmosphere, radians
    e = nex('arcsin(sin(latrad)*sin(deltap) + cos(latrad)*cos(deltap)*cos(Hp))')

    # topocentric elevation angle with refraction correction, radians
    if with_refraction is True:
        refr_corr = nex('ATM_REFR*1.7802e-2/(60*tan(e + 3.1376e-3 / (e + 8.9186e-2)))')
        e = nex('where(e >= -1.4545e-2, e + refr_corr, e)')

    zenith = HALFPI - e

    # topocentric azimuth angle, radians
    azimuth = nex('arctan2(sin(Hp), cos(Hp)*sin(latrad) - tan(deltap)*cos(latrad))')

    # sun-earth orbit eccentricity correction factor
    ecf = nex('1 / R**2')

    return zenith, azimuth, ecf, delta, eot


def _sunpos_two_dimensional(unixtime, longitude, latitude, with_refraction):

    # unixtime: seconds since 1970-01-01T00 UTC

    nex = ne.evaluate

    alpha, delta, eot, v, R  = time_dependent_calculations(unixtime)

    # start of spatial calculations...

    lon_ = longitude[None, :]
    latrad_ = latitude[None, :]*DEG2RAD
    alpha_ = alpha[:, None]
    delta_ = delta[:, None]
    v_ = v[:, None]

    # local hour angle, radians
    H = nex('((v_ + lon_ - alpha_) % 360)*DEG2RAD')

    # equatorial horizontal parallax
    xi = nex('(8.794 / (3600 * R))*DEG2RAD')[:, None]  # radians
    u = nex('arctan(0.99664719*tan(latrad_))')
    x = nex('cos(u) + ELEVATION / 6378140. * cos(latrad_)')
    y = nex('0.99664719*sin(u) + ELEVATION / 6378140. * sin(latrad_)')

    # parallax sun right ascension, radians
    deltaa = nex('arctan2(-x*sin(xi)*sin(H), cos(delta_) - x*sin(xi)*cos(H))')

    # topocentric sun declination, radians
    deltap = nex('arctan2((sin(delta_)-y*sin(xi))*cos(deltaa), cos(delta_)-x*sin(xi)*cos(H))')

    # topocentric local hour angle, radians
    Hp = nex('H - deltaa')

    # topocentric elevation angle without atmosphere, radians
    e = nex('arcsin(sin(latrad_)*sin(deltap) + cos(latrad_)*cos(deltap)*cos(Hp))')

    # topocentric elevation angle with refraction correction, radians
    if with_refraction is True:
        refr_corr = nex('ATM_REFR*1.7802e-2/(60*tan(e + 3.1376e-3 / (e + 8.9186e-2)))')
        e = nex('where(e >= -1.4545e-2, e + refr_corr, e)')

    zenith = HALFPI - e

    # topocentric azimuth angle, radians
    azimuth = nex('arctan2(sin(Hp), cos(Hp)*sin(latrad_) - tan(deltap)*cos(latrad_))')

    # sun-earth orbit eccentricity correction factor
    ecf = nex('1 / R**2')

    return zenith, azimuth, ecf, delta, eot


def _sunpos_three_dimensional(unixtime, longitude, latitude, with_refraction):

    # unixtime: seconds since 1970-01-01T00 UTC

    nex = ne.evaluate

    alpha, delta, eot, v, R = time_dependent_calculations(unixtime)

    # start of spatial calculations...

    lon_ = longitude[None, None, :]
    latrad_ = latitude[None, :, None]*DEG2RAD
    alpha_ = alpha[:, None, None]
    delta_ = delta[:, None, None]
    v_ = v[:, None, None]

    # local hour angle, radians
    H = nex('((v_ + lon_ - alpha_) % 360)*DEG2RAD')

    # equatorial horizontal parallax
    xi = nex('(8.794 / (3600 * R))*DEG2RAD')[:, None, None]  # radians
    u = nex('arctan(0.99664719*tan(latrad_))')
    x = nex('cos(u) + ELEVATION / 6378140. * cos(latrad_)')
    y = nex('0.99664719*sin(u) + ELEVATION / 6378140. * sin(latrad_)')

    # parallax sun right ascension, radians
    deltaa = nex('arctan2(-x*sin(xi)*sin(H), cos(delta_) - x*sin(xi)*cos(H))')

    # topocentric sun declination, radians
    deltap = nex('arctan2((sin(delta_)-y*sin(xi))*cos(deltaa), cos(delta_)-x*sin(xi)*cos(H))')

    # topocentric local hour angle, radians
    Hp = nex('H - deltaa')

    # topocentric elevation angle without atmosphere, radians
    e = nex('arcsin(sin(latrad_)*sin(deltap) + cos(latrad_)*cos(deltap)*cos(Hp))')

    # topocentric elevation angle with refraction correction, radians
    if with_refraction is True:
        refr_corr = nex('ATM_REFR*1.7802e-2/(60*tan(e + 3.1376e-3 / (e + 8.9186e-2)))')
        e = nex('where(e >= -1.4545e-2, e + refr_corr, e)')

    zenith = HALFPI - e

    # topocentric azimuth angle, radians
    azimuth = nex('arctan2(sin(Hp), cos(Hp)*sin(latrad_) - tan(deltap)*cos(latrad_))')

    # sun-earth orbit eccentricity correction factor
    ecf = nex('1 / R**2')

    return zenith, azimuth, ecf, delta, eot


############################################################
#  PARAMETERS TO CALCULATE THE HELIOCENTRIC RADIUS VECTOR  #
############################################################

PARAMS_R_0 = np.array(
    [[100013989.0, 0.0, 0.0],
     [1670700.0, 3.0984635, 6283.07585],
     [13956.0, 3.05525, 12566.1517],
     [3084.0, 5.1985, 77713.7715],
     [1628.0, 1.1739, 5753.3849],
     [1576.0, 2.8469, 7860.4194],
     [925.0, 5.453, 11506.77],
     [542.0, 4.564, 3930.21],
     [472.0, 3.661, 5884.927],
     [346.0, 0.964, 5507.553],
     [329.0, 5.9, 5223.694],
     [307.0, 0.299, 5573.143],
     [243.0, 4.273, 11790.629],
     [212.0, 5.847, 1577.344],
     [186.0, 5.022, 10977.079],
     [175.0, 3.012, 18849.228],
     [110.0, 5.055, 5486.778],
     [98.0, 0.89, 6069.78],
     [86.0, 5.69, 15720.84],
     [86.0, 1.27, 161000.69],
     [65.0, 0.27, 17260.15],
     [63.0, 0.92, 529.69],
     [57.0, 2.01, 83996.85],
     [56.0, 5.24, 71430.7],
     [49.0, 3.25, 2544.31],
     [47.0, 2.58, 775.52],
     [45.0, 5.54, 9437.76],
     [43.0, 6.01, 6275.96],
     [39.0, 5.36, 4694.0],
     [38.0, 2.39, 8827.39],
     [37.0, 0.83, 19651.05],
     [37.0, 4.9, 12139.55],
     [36.0, 1.67, 12036.46],
     [35.0, 1.84, 2942.46],
     [33.0, 0.24, 7084.9],
     [32.0, 0.18, 5088.63],
     [32.0, 1.78, 398.15],
     [28.0, 1.21, 6286.6],
     [28.0, 1.9, 6279.55],
     [26.0, 4.59, 10447.39]]
)
PARAMS_R_0.resize((40, 3), refcheck=False)

PARAMS_R_1 = np.array(
 [[103019.0, 1.10749, 6283.07585],
  [1721.0, 1.0644, 12566.1517],
  [702.0, 3.142, 0.0],
  [32.0, 1.02, 18849.23],
  [31.0, 2.84, 5507.55],
  [25.0, 1.32, 5223.69],
  [18.0, 1.42, 1577.34],
  [10.0, 5.91, 10977.08],
  [9.0, 1.42, 6275.96],
  [9.0, 0.27, 5486.78]]
)
PARAMS_R_1.resize((40, 3), refcheck=False)

PARAMS_R_2 = np.array(
   [[4359.0, 5.7846, 6283.0758],
    [124.0, 5.579, 12566.152],
    [12.0, 3.14, 0.0],
    [9.0, 3.63, 77713.77],
    [6.0, 1.87, 5573.14],
    [3.0, 5.47, 18849.23]]
)
PARAMS_R_2.resize((40, 3), refcheck=False)

PARAMS_R_3 = np.array(
   [[145.0, 4.273, 6283.076],
    [7.0, 3.92, 12566.15]]
)
PARAMS_R_3.resize((40, 3), refcheck=False)

PARAMS_R_4 = np.array(
   [[4.0, 2.56, 6283.08]]
)
PARAMS_R_4.resize((40, 3), refcheck=False)

PARAMS_R_i = np.stack(
    [PARAMS_R_0, PARAMS_R_1, PARAMS_R_2,
     PARAMS_R_3, PARAMS_R_4])
C_R_0 = PARAMS_R_i[:, :, 0:1]
C_R_1 = PARAMS_R_i[:, :, 1:2]
C_R_2 = PARAMS_R_i[:, :, 2:]

########################################################
#  PARAMETERS TO CALCULATE THE HELIOCENTRIC LONGITUDE  #
########################################################

PARAMS_LONG_0 = np.array(
   [[175347046.0, 0.0, 0.0],
    [3341656.0, 4.6692568, 6283.07585],
    [34894.0, 4.6261, 12566.1517],
    [3497.0, 2.7441, 5753.3849],
    [3418.0, 2.8289, 3.5231],
    [3136.0, 3.6277, 77713.7715],
    [2676.0, 4.4181, 7860.4194],
    [2343.0, 6.1352, 3930.2097],
    [1324.0, 0.7425, 11506.7698],
    [1273.0, 2.0371, 529.691],
    [1199.0, 1.1096, 1577.3435],
    [990.0, 5.233, 5884.927],
    [902.0, 2.045, 26.298],
    [857.0, 3.508, 398.149],
    [780.0, 1.179, 5223.694],
    [753.0, 2.533, 5507.553],
    [505.0, 4.583, 18849.228],
    [492.0, 4.205, 775.523],
    [357.0, 2.92, 0.067],
    [317.0, 5.849, 11790.629],
    [284.0, 1.899, 796.298],
    [271.0, 0.315, 10977.079],
    [243.0, 0.345, 5486.778],
    [206.0, 4.806, 2544.314],
    [205.0, 1.869, 5573.143],
    [202.0, 2.458, 6069.777],
    [156.0, 0.833, 213.299],
    [132.0, 3.411, 2942.463],
    [126.0, 1.083, 20.775],
    [115.0, 0.645, 0.98],
    [103.0, 0.636, 4694.003],
    [102.0, 0.976, 15720.839],
    [102.0, 4.267, 7.114],
    [99.0, 6.21, 2146.17],
    [98.0, 0.68, 155.42],
    [86.0, 5.98, 161000.69],
    [85.0, 1.3, 6275.96],
    [85.0, 3.67, 71430.7],
    [80.0, 1.81, 17260.15],
    [79.0, 3.04, 12036.46],
    [75.0, 1.76, 5088.63],
    [74.0, 3.5, 3154.69],
    [74.0, 4.68, 801.82],
    [70.0, 0.83, 9437.76],
    [62.0, 3.98, 8827.39],
    [61.0, 1.82, 7084.9],
    [57.0, 2.78, 6286.6],
    [56.0, 4.39, 14143.5],
    [56.0, 3.47, 6279.55],
    [52.0, 0.19, 12139.55],
    [52.0, 1.33, 1748.02],
    [51.0, 0.28, 5856.48],
    [49.0, 0.49, 1194.45],
    [41.0, 5.37, 8429.24],
    [41.0, 2.4, 19651.05],
    [39.0, 6.17, 10447.39],
    [37.0, 6.04, 10213.29],
    [37.0, 2.57, 1059.38],
    [36.0, 1.71, 2352.87],
    [36.0, 1.78, 6812.77],
    [33.0, 0.59, 17789.85],
    [30.0, 0.44, 83996.85],
    [30.0, 2.74, 1349.87],
    [25.0, 3.16, 4690.48]])
PARAMS_LONG_0.resize((64, 3), refcheck=False)

PARAMS_LONG_1 = np.array(
   [[628331966747.0, 0.0, 0.0],
    [206059.0, 2.678235, 6283.07585],
    [4303.0, 2.6351, 12566.1517],
    [425.0, 1.59, 3.523],
    [119.0, 5.796, 26.298],
    [109.0, 2.966, 1577.344],
    [93.0, 2.59, 18849.23],
    [72.0, 1.14, 529.69],
    [68.0, 1.87, 398.15],
    [67.0, 4.41, 5507.55],
    [59.0, 2.89, 5223.69],
    [56.0, 2.17, 155.42],
    [45.0, 0.4, 796.3],
    [36.0, 0.47, 775.52],
    [29.0, 2.65, 7.11],
    [21.0, 5.34, 0.98],
    [19.0, 1.85, 5486.78],
    [19.0, 4.97, 213.3],
    [17.0, 2.99, 6275.96],
    [16.0, 0.03, 2544.31],
    [16.0, 1.43, 2146.17],
    [15.0, 1.21, 10977.08],
    [12.0, 2.83, 1748.02],
    [12.0, 3.26, 5088.63],
    [12.0, 5.27, 1194.45],
    [12.0, 2.08, 4694.0],
    [11.0, 0.77, 553.57],
    [10.0, 1.3, 6286.6],
    [10.0, 4.24, 1349.87],
    [9.0, 2.7, 242.73],
    [9.0, 5.64, 951.72],
    [8.0, 5.3, 2352.87],
    [6.0, 2.65, 9437.76],
    [6.0, 4.67, 4690.48]])
PARAMS_LONG_1.resize((64, 3), refcheck=False)

PARAMS_LONG_2 = np.array(
   [[52919.0, 0.0, 0.0],
    [8720.0, 1.0721, 6283.0758],
    [309.0, 0.867, 12566.152],
    [27.0, 0.05, 3.52],
    [16.0, 5.19, 26.3],
    [16.0, 3.68, 155.42],
    [10.0, 0.76, 18849.23],
    [9.0, 2.06, 77713.77],
    [7.0, 0.83, 775.52],
    [5.0, 4.66, 1577.34],
    [4.0, 1.03, 7.11],
    [4.0, 3.44, 5573.14],
    [3.0, 5.14, 796.3],
    [3.0, 6.05, 5507.55],
    [3.0, 1.19, 242.73],
    [3.0, 6.12, 529.69],
    [3.0, 0.31, 398.15],
    [3.0, 2.28, 553.57],
    [2.0, 4.38, 5223.69],
    [2.0, 3.75, 0.98]])
PARAMS_LONG_2.resize((64, 3), refcheck=False)

PARAMS_LONG_3 = np.array(
   [[289.0, 5.844, 6283.076],
    [35.0, 0.0, 0.0],
    [17.0, 5.49, 12566.15],
    [3.0, 5.2, 155.42],
    [1.0, 4.72, 3.52],
    [1.0, 5.3, 18849.23],
    [1.0, 5.97, 242.73]])
PARAMS_LONG_3.resize((64, 3), refcheck=False)

PARAMS_LONG_4 = np.array(
   [[114.0, 3.142, 0.0],
    [8.0, 4.13, 6283.08],
    [1.0, 3.84, 12566.15]])
PARAMS_LONG_4.resize((64, 3), refcheck=False)

PARAMS_LONG_5 = np.array(
   [[1.0, 3.14, 0.0]])
PARAMS_LONG_5.resize((64, 3), refcheck=False)

PARAMS_LONG_i = np.stack(
    [PARAMS_LONG_0, PARAMS_LONG_1, PARAMS_LONG_2,
     PARAMS_LONG_3, PARAMS_LONG_4, PARAMS_LONG_5])
C_LONG_0 = PARAMS_LONG_i[:, :, 0:1]
C_LONG_1 = PARAMS_LONG_i[:, :, 1:2]
C_LONG_2 = PARAMS_LONG_i[:, :, 2:]

#######################################################
#  PARAMETERS TO CALCULATE THE HELIOCENTRIC LATITUDE  #
#######################################################

PARAMS_LAT_0 = np.array(
   [[280.0, 3.199, 84334.662],
    [102.0, 5.422, 5507.553],
    [80.0, 3.88, 5223.69],
    [44.0, 3.7, 2352.87],
    [32.0, 4.0, 1577.34]])
PARAMS_LAT_0.resize((5, 3), refcheck=False)

PARAMS_LAT_1 = np.array(
   [[9.0, 3.9, 5507.55],
    [6.0, 1.73, 5223.69]])
PARAMS_LAT_1.resize((5, 3), refcheck=False)

PARAMS_LAT_i = np.stack([PARAMS_LAT_0, PARAMS_LAT_1])
C_LAT_0 = PARAMS_LAT_i[:, :, 0:1]
C_LAT_1 = PARAMS_LAT_i[:, :, 1:2]
C_LAT_2 = PARAMS_LAT_i[:, :, 2:]

#####################################################
#  COEFFICIENTS A, B, C AND D TO CALCULATE NUTATION #
#####################################################

PARAMS_NUTATION_ABCD = np.array([
    [-171996, -174.2, 92025, 8.9],
    [-13187, -1.6, 5736, -3.1],
    [-2274, -0.2, 977, -0.5],
    [2062, 0.2, -895, 0.5],
    [1426, -3.4, 54, -0.1],
    [712, 0.1, -7, 0],
    [-517, 1.2, 224, -0.6],
    [-386, -0.4, 200, 0],
    [-301, 0, 129, -0.1],
    [217, -0.5, -95, 0.3],
    [-158, 0, 0, 0],
    [129, 0.1, -70, 0],
    [123, 0, -53, 0],
    [63, 0, 0, 0],
    [63, 0.1, -33, 0],
    [-59, 0, 26, 0],
    [-58, -0.1, 32, 0],
    [-51, 0, 27, 0],
    [48, 0, 0, 0],
    [46, 0, -24, 0],
    [-38, 0, 16, 0],
    [-31, 0, 13, 0],
    [29, 0, 0, 0],
    [29, 0, -12, 0],
    [26, 0, 0, 0],
    [-22, 0, 0, 0],
    [21, 0, -10, 0],
    [17, -0.1, 0, 0],
    [16, 0, -8, 0],
    [-16, 0.1, 7, 0],
    [-15, 0, 9, 0],
    [-13, 0, 7, 0],
    [-12, 0, 6, 0],
    [11, 0, 0, 0],
    [-10, 0, 5, 0],
    [-8, 0, 3, 0],
    [7, 0, -3, 0],
    [-7, 0, 0, 0],
    [-7, 0, 3, 0],
    [-7, 0, 3, 0],
    [6, 0, 0, 0],
    [6, 0, -3, 0],
    [6, 0, -3, 0],
    [-6, 0, 3, 0],
    [-6, 0, 3, 0],
    [5, 0, 0, 0],
    [-5, 0, 3, 0],
    [-5, 0, 3, 0],
    [-5, 0, 3, 0],
    [4, 0, 0, 0],
    [4, 0, 0, 0],
    [4, 0, 0, 0],
    [-4, 0, 0, 0],
    [-4, 0, 0, 0],
    [-4, 0, 0, 0],
    [3, 0, 0, 0],
    [-3, 0, 0, 0],
    [-3, 0, 0, 0],
    [-3, 0, 0, 0],
    [-3, 0, 0, 0],
    [-3, 0, 0, 0],
    [-3, 0, 0, 0],
    [-3, 0, 0, 0]])

C_NUT_A = PARAMS_NUTATION_ABCD[:, 0:1]
C_NUT_B = PARAMS_NUTATION_ABCD[:, 1:2]
C_NUT_C = PARAMS_NUTATION_ABCD[:, 2:3]
C_NUT_D = PARAMS_NUTATION_ABCD[:, 3:]

PARAMS_NUTATION_YVALUES = np.array([
    [0, 0, 0, 0, 1],
    [-2, 0, 0, 2, 2],
    [0, 0, 0, 2, 2],
    [0, 0, 0, 0, 2],
    [0, 1, 0, 0, 0],
    [0, 0, 1, 0, 0],
    [-2, 1, 0, 2, 2],
    [0, 0, 0, 2, 1],
    [0, 0, 1, 2, 2],
    [-2, -1, 0, 2, 2],
    [-2, 0, 1, 0, 0],
    [-2, 0, 0, 2, 1],
    [0, 0, -1, 2, 2],
    [2, 0, 0, 0, 0],
    [0, 0, 1, 0, 1],
    [2, 0, -1, 2, 2],
    [0, 0, -1, 0, 1],
    [0, 0, 1, 2, 1],
    [-2, 0, 2, 0, 0],
    [0, 0, -2, 2, 1],
    [2, 0, 0, 2, 2],
    [0, 0, 2, 2, 2],
    [0, 0, 2, 0, 0],
    [-2, 0, 1, 2, 2],
    [0, 0, 0, 2, 0],
    [-2, 0, 0, 2, 0],
    [0, 0, -1, 2, 1],
    [0, 2, 0, 0, 0],
    [2, 0, -1, 0, 1],
    [-2, 2, 0, 2, 2],
    [0, 1, 0, 0, 1],
    [-2, 0, 1, 0, 1],
    [0, -1, 0, 0, 1],
    [0, 0, 2, -2, 0],
    [2, 0, -1, 2, 1],
    [2, 0, 1, 2, 2],
    [0, 1, 0, 2, 2],
    [-2, 1, 1, 0, 0],
    [0, -1, 0, 2, 2],
    [2, 0, 0, 2, 1],
    [2, 0, 1, 0, 0],
    [-2, 0, 2, 2, 2],
    [-2, 0, 1, 2, 1],
    [2, 0, -2, 0, 1],
    [2, 0, 0, 0, 1],
    [0, -1, 1, 0, 0],
    [-2, -1, 0, 2, 1],
    [-2, 0, 0, 0, 1],
    [0, 0, 2, 2, 1],
    [-2, 0, 2, 0, 1],
    [-2, 1, 0, 2, 1],
    [0, 0, 1, -2, 0],
    [-1, 0, 1, 0, 0],
    [-2, 1, 0, 0, 0],
    [1, 0, 0, 0, 0],
    [0, 0, 1, 2, 0],
    [0, 0, -2, 2, 2],
    [-1, -1, 1, 0, 0],
    [0, 1, 1, 0, 0],
    [0, -1, 1, 2, 2],
    [2, -1, -1, 2, 2],
    [0, 0, 3, 2, 2],
    [2, -1, 0, 2, 2]])

C_NUT_Y_0 = PARAMS_NUTATION_YVALUES[:, 0:1]
C_NUT_Y_1 = PARAMS_NUTATION_YVALUES[:, 1:2]
C_NUT_Y_2 = PARAMS_NUTATION_YVALUES[:, 2:3]
C_NUT_Y_3 = PARAMS_NUTATION_YVALUES[:, 3:4]
C_NUT_Y_4 = PARAMS_NUTATION_YVALUES[:, 4:]
