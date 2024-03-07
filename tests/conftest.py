
import numpy as np

from pytest_allclose import report_rmses


def cosr(a):
    return np.cos(np.radians(a))


def sinr(a):
    return np.sin(np.radians(a))


def get_da_for_sites(da):
    return da.to_numpy()


def get_da_for_regular_grid(da):
    if da.coords.dims == ('time',):
        return da.to_numpy()
    if da.coords.dims == ('time', 'location'):
        return da.isel(location=0).to_numpy()
    return da.isel(latitude=0, longitude=0).to_numpy()


def get_da_for_transect(da):
    return get_da_for_regular_grid(da)


def angular_distance(a1, a2):
    """Used to test the azimuth angle distance between SPAs"""
    dotp = cosr(a1)*cosr(a2) + sinr(a1)*sinr(a2)
    dotp = np.round(dotp, 6)  # to prevent nans when dotp is only slightly higher than 1
    return np.degrees(np.arccos(dotp))


def azimdist(a1, a2, thresh):
    ang = angular_distance(a1, a2)
    condition = ang < 5.
    errmsg = (f'test failed at: {a1[~condition]} and {a2[~condition]} '
              f'with values {ang[~condition]}')
    return np.all(condition), errmsg


def pytest_terminal_summary(terminalreporter):
    report_rmses(terminalreporter)
