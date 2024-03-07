
import pytest
import numpy as np
import pandas as pd

from conftest import (
    get_da_for_sites,
    get_da_for_regular_grid,
    get_da_for_transect,
    azimdist
)

import sunwhere  # pylint: disable=import-error


@pytest.fixture(scope='module', autouse=True)
def space_time_sites():
    N_PERIODS = 2*24*365
    N_LOCATIONS = 200
    times = pd.date_range('2024-01-01', periods=N_PERIODS, freq='H')
    rg = np.random.default_rng()
    lats = rg.uniform(-90., 90., N_LOCATIONS)
    lons = rg.uniform(-180., 180., N_LOCATIONS)
    return (times, lats, lons)


@pytest.fixture(scope='module', autouse=True)
def space_time_regular_grid():
    N_PERIODS = 2*24*365
    N_LATS, N_LONS = 1, 1
    times = pd.date_range('2024-01-01', periods=N_PERIODS, freq='H')
    rg = np.random.default_rng()
    lats = rg.uniform(-90., 90., N_LATS)
    lons = rg.uniform(-180., 180., N_LONS)
    return (times, lats, lons)


@pytest.fixture(scope='module', autouse=True)
def space_time_transect():
    N_STEPS = 50
    times = pd.date_range('2024-01-01', periods=N_STEPS, freq='H')
    rg = np.random.default_rng()
    lats = rg.uniform(-90., 90., N_STEPS)
    lons = rg.uniform(-180., 180., N_STEPS)
    return (times, lats, lons)


@pytest.fixture(scope='module', autouse=True)
def sites_nrel(space_time_sites):
    kwargs = {'algorithm': 'nrel', 'engine': 'numexpr', 'refraction': False}
    return sunwhere.sites(*space_time_sites, **kwargs)


@pytest.mark.filterwarnings("ignore")
def test_psa_numpy_consistency(space_time_sites, sites_nrel, allclose):
    kwargs = {'algorithm': 'psa', 'engine': 'numpy', 'refraction': False}
    nrel = sites_nrel
    get = get_da_for_sites
    psa = sunwhere.sites(*space_time_sites, **kwargs)
    assert allclose(get(nrel.sza), get(psa.sza), atol=2., rtol=0.)
    assert allclose(get(nrel.dec), get(psa.dec), atol=2., rtol=0.)
    assert allclose(get(nrel.eot), get(psa.eot), atol=10, rtol=0.)
    assert allclose(get(nrel.ecf), get(psa.ecf), atol=0., rtol=0.01)
    assert azimdist(get(nrel.saa), get(psa.saa), thresh=1.)


@pytest.mark.filterwarnings("ignore")
def test_psa_numexpr_consistency(space_time_sites, sites_nrel, allclose):
    kwargs = {'algorithm': 'psa', 'engine': 'numexpr', 'refraction': False}
    nrel = sites_nrel
    get = get_da_for_sites
    psa = sunwhere.sites(*space_time_sites, **kwargs)
    assert allclose(get(nrel.sza), get(psa.sza), atol=2., rtol=0.)
    assert allclose(get(nrel.dec), get(psa.dec), atol=2., rtol=0.)
    assert allclose(get(nrel.eot), get(psa.eot), atol=10, rtol=0.)
    assert allclose(get(nrel.ecf), get(psa.ecf), atol=0., rtol=0.01)
    assert azimdist(get(nrel.saa), get(psa.saa), thresh=1.)


@pytest.mark.filterwarnings("ignore")
def test_soltrack_numpy_consistency(space_time_sites, sites_nrel, allclose):
    kwargs = {'algorithm': 'soltrack', 'engine': 'numpy', 'refraction': False}
    nrel = sites_nrel
    get = get_da_for_sites
    soltrack = sunwhere.sites(*space_time_sites, **kwargs)
    assert allclose(get(nrel.sza), get(soltrack.sza), atol=2., rtol=0.)
    assert allclose(get(nrel.dec), get(soltrack.dec), atol=2., rtol=0.)
    assert allclose(get(nrel.eot), get(soltrack.eot), atol=10, rtol=0.)
    assert allclose(get(nrel.ecf), get(soltrack.ecf), atol=0., rtol=0.01)
    assert azimdist(get(nrel.saa), get(soltrack.saa), thresh=1.)


@pytest.mark.filterwarnings("ignore")
def test_soltrack_numexpr_consistency(space_time_sites, sites_nrel, allclose):
    kwargs = {'algorithm': 'soltrack', 'engine': 'numexpr', 'refraction': False}
    nrel = sites_nrel
    get = get_da_for_sites
    soltrack = sunwhere.sites(*space_time_sites, **kwargs)
    assert allclose(get(nrel.sza), get(soltrack.sza), atol=2., rtol=0.)
    assert allclose(get(nrel.dec), get(soltrack.dec), atol=2., rtol=0.)
    assert allclose(get(nrel.eot), get(soltrack.eot), atol=10, rtol=0.)
    assert allclose(get(nrel.ecf), get(soltrack.ecf), atol=0., rtol=0.01)
    assert azimdist(get(nrel.saa), get(soltrack.saa), thresh=1.)


@pytest.mark.filterwarnings("ignore")
def test_iqbal_numpy_consistency(space_time_sites, sites_nrel, allclose):
    kwargs = {'algorithm': 'iqbal', 'engine': 'numpy', 'refraction': False}
    nrel = sites_nrel
    get = get_da_for_sites
    iqbal = sunwhere.sites(*space_time_sites, **kwargs)
    assert allclose(get(nrel.sza), get(iqbal.sza), atol=2., rtol=0.)
    assert allclose(get(nrel.dec), get(iqbal.dec), atol=2., rtol=0.)
    assert allclose(get(nrel.eot), get(iqbal.eot), atol=10, rtol=0.)
    assert allclose(get(nrel.ecf), get(iqbal.ecf), atol=0., rtol=0.01)
    assert azimdist(get(nrel.saa), get(iqbal.saa), thresh=1.)


@pytest.mark.filterwarnings("ignore")
def test_iqbal_numexpr_consistency(space_time_sites, sites_nrel, allclose):
    kwargs = {'algorithm': 'iqbal', 'engine': 'numexpr', 'refraction': False}
    nrel = sites_nrel
    get = get_da_for_sites
    iqbal = sunwhere.sites(*space_time_sites, **kwargs)
    assert allclose(get(nrel.sza), get(iqbal.sza), atol=2., rtol=0.)
    assert allclose(get(nrel.dec), get(iqbal.dec), atol=2., rtol=0.)
    assert allclose(get(nrel.eot), get(iqbal.eot), atol=10, rtol=0.)
    assert allclose(get(nrel.ecf), get(iqbal.ecf), atol=0., rtol=0.01)
    assert azimdist(get(nrel.saa), get(iqbal.saa), thresh=1.)


@pytest.mark.filterwarnings("ignore")
def test_nrel_numexpr_regular_grid_consistency(space_time_regular_grid, allclose):
    kwargs = {'algorithm': 'nrel', 'engine': 'numexpr', 'refraction': False}
    get = get_da_for_regular_grid
    sites = sunwhere.sites(*space_time_regular_grid, **kwargs)
    rgrid = sunwhere.regular_grid(*space_time_regular_grid, **kwargs)
    assert allclose(get(sites.sza), get(rgrid.sza), atol=2., rtol=0.)
    assert allclose(get(sites.dec), get(rgrid.dec), atol=2., rtol=0.)
    assert allclose(get(sites.eot), get(rgrid.eot), atol=10., rtol=0.)
    assert allclose(get(sites.ecf), get(rgrid.ecf), atol=0., rtol=0.01)
    assert azimdist(get(sites.saa), get(rgrid.saa), thresh=1.)


@pytest.mark.filterwarnings("ignore")
def test_psa_numpy_regular_grid_consistency(space_time_regular_grid, allclose):
    kwargs = {'algorithm': 'psa', 'engine': 'numpy', 'refraction': False}
    get = get_da_for_regular_grid
    sites = sunwhere.sites(*space_time_regular_grid, **kwargs)
    rgrid = sunwhere.regular_grid(*space_time_regular_grid, **kwargs)
    assert allclose(get(sites.sza), get(rgrid.sza), atol=2., rtol=0.)
    assert allclose(get(sites.dec), get(rgrid.dec), atol=2., rtol=0.)
    assert allclose(get(sites.eot), get(rgrid.eot), atol=10., rtol=0.)
    assert allclose(get(sites.ecf), get(rgrid.ecf), atol=0., rtol=0.01)
    assert azimdist(get(sites.saa), get(rgrid.saa), thresh=1.)


@pytest.mark.filterwarnings("ignore")
def test_psa_numexpr_regular_grid_consistency(space_time_regular_grid, allclose):
    kwargs = {'algorithm': 'psa', 'engine': 'numexpr', 'refraction': False}
    get = get_da_for_regular_grid
    sites = sunwhere.sites(*space_time_regular_grid, **kwargs)
    rgrid = sunwhere.regular_grid(*space_time_regular_grid, **kwargs)
    assert allclose(get(sites.sza), get(rgrid.sza), atol=2., rtol=0.)
    assert allclose(get(sites.dec), get(rgrid.dec), atol=2., rtol=0.)
    assert allclose(get(sites.eot), get(rgrid.eot), atol=10., rtol=0.)
    assert allclose(get(sites.ecf), get(rgrid.ecf), atol=0., rtol=0.01)
    assert azimdist(get(sites.saa), get(rgrid.saa), thresh=1.)


@pytest.mark.filterwarnings("ignore")
def test_soltrack_numpy_regular_grid_consistency(space_time_regular_grid, allclose):
    kwargs = {'algorithm': 'soltrack', 'engine': 'numpy', 'refraction': False}
    get = get_da_for_regular_grid
    sites = sunwhere.sites(*space_time_regular_grid, **kwargs)
    rgrid = sunwhere.regular_grid(*space_time_regular_grid, **kwargs)
    assert allclose(get(sites.sza), get(rgrid.sza), atol=2., rtol=0.)
    assert allclose(get(sites.dec), get(rgrid.dec), atol=2., rtol=0.)
    assert allclose(get(sites.eot), get(rgrid.eot), atol=10., rtol=0.)
    assert allclose(get(sites.ecf), get(rgrid.ecf), atol=0., rtol=0.01)
    assert azimdist(get(sites.saa), get(rgrid.saa), thresh=1.)


@pytest.mark.filterwarnings("ignore")
def test_soltrack_numexpr_regular_grid_consistency(space_time_regular_grid, allclose):
    kwargs = {'algorithm': 'soltrack', 'engine': 'numexpr', 'refraction': False}
    get = get_da_for_regular_grid
    sites = sunwhere.sites(*space_time_regular_grid, **kwargs)
    rgrid = sunwhere.regular_grid(*space_time_regular_grid, **kwargs)
    assert allclose(get(sites.sza), get(rgrid.sza), atol=2., rtol=0.)
    assert allclose(get(sites.dec), get(rgrid.dec), atol=2., rtol=0.)
    assert allclose(get(sites.eot), get(rgrid.eot), atol=10., rtol=0.)
    assert allclose(get(sites.ecf), get(rgrid.ecf), atol=0., rtol=0.01)
    assert azimdist(get(sites.saa), get(rgrid.saa), thresh=1.)


@pytest.mark.filterwarnings("ignore")
def test_iqbal_numpy_regular_grid_consistency(space_time_regular_grid, allclose):
    kwargs = {'algorithm': 'iqbal', 'engine': 'numpy', 'refraction': False}
    get = get_da_for_regular_grid
    sites = sunwhere.sites(*space_time_regular_grid, **kwargs)
    rgrid = sunwhere.regular_grid(*space_time_regular_grid, **kwargs)
    assert allclose(get(sites.sza), get(rgrid.sza), atol=2., rtol=0.)
    assert allclose(get(sites.dec), get(rgrid.dec), atol=2., rtol=0.)
    assert allclose(get(sites.eot), get(rgrid.eot), atol=10., rtol=0.)
    assert allclose(get(sites.ecf), get(rgrid.ecf), atol=0., rtol=0.01)
    assert azimdist(get(sites.saa), get(rgrid.saa), thresh=1.)


@pytest.mark.filterwarnings("ignore")
def test_iqbal_numexpr_regular_grid_consistency(space_time_regular_grid, allclose):
    kwargs = {'algorithm': 'iqbal', 'engine': 'numexpr', 'refraction': False}
    get = get_da_for_regular_grid
    sites = sunwhere.sites(*space_time_regular_grid, **kwargs)
    rgrid = sunwhere.regular_grid(*space_time_regular_grid, **kwargs)
    assert allclose(get(sites.sza), get(rgrid.sza), atol=2., rtol=0.)
    assert allclose(get(sites.dec), get(rgrid.dec), atol=2., rtol=0.)
    assert allclose(get(sites.eot), get(rgrid.eot), atol=10., rtol=0.)
    assert allclose(get(sites.ecf), get(rgrid.ecf), atol=0., rtol=0.01)
    assert azimdist(get(sites.saa), get(rgrid.saa), thresh=1.)


@pytest.mark.filterwarnings("ignore")
def test_nrel_numexpr_transect_consistency(space_time_transect, allclose):
    kwargs = {'algorithm': 'nrel', 'engine': 'numexpr', 'refraction': False}
    get = get_da_for_transect
    for time, lat, lon in zip(*space_time_transect):
        sites = sunwhere.sites(time, lat, lon, **kwargs)
        trans = sunwhere.transect(time, lat, lon, **kwargs)
        assert allclose(get(sites.sza), get(trans.sza), atol=2., rtol=0.)
        assert allclose(get(sites.dec), get(trans.dec), atol=2., rtol=0.)
        assert allclose(get(sites.eot), get(trans.eot), atol=10., rtol=0.)
        assert allclose(get(sites.ecf), get(trans.ecf), atol=0., rtol=0.01)
        assert azimdist(get(sites.saa), get(trans.saa), thresh=1.)


@pytest.mark.filterwarnings("ignore")
def test_psa_numpy_transect_consistency(space_time_transect, allclose):
    kwargs = {'algorithm': 'psa', 'engine': 'numpy', 'refraction': False}
    get = get_da_for_transect
    for time, lat, lon in zip(*space_time_transect):
        sites = sunwhere.sites(time, lat, lon, **kwargs)
        trans = sunwhere.transect(time, lat, lon, **kwargs)
        assert allclose(get(sites.sza), get(trans.sza), atol=2., rtol=0.)
        assert allclose(get(sites.dec), get(trans.dec), atol=2., rtol=0.)
        assert allclose(get(sites.eot), get(trans.eot), atol=10., rtol=0.)
        assert allclose(get(sites.ecf), get(trans.ecf), atol=0., rtol=0.01)
        assert azimdist(get(sites.saa), get(trans.saa), thresh=1.)


@pytest.mark.filterwarnings("ignore")
def test_psa_numexpr_transect_consistency(space_time_transect, allclose):
    kwargs = {'algorithm': 'psa', 'engine': 'numexpr', 'refraction': False}
    get = get_da_for_transect
    for time, lat, lon in zip(*space_time_transect):
        sites = sunwhere.sites(time, lat, lon, **kwargs)
        trans = sunwhere.transect(time, lat, lon, **kwargs)
        assert allclose(get(sites.sza), get(trans.sza), atol=2., rtol=0.)
        assert allclose(get(sites.dec), get(trans.dec), atol=2., rtol=0.)
        assert allclose(get(sites.eot), get(trans.eot), atol=10., rtol=0.)
        assert allclose(get(sites.ecf), get(trans.ecf), atol=0., rtol=0.01)
        assert azimdist(get(sites.saa), get(trans.saa), thresh=1.)


@pytest.mark.filterwarnings("ignore")
def test_soltrack_numpy_transect_consistency(space_time_transect, allclose):
    kwargs = {'algorithm': 'soltrack', 'engine': 'numpy', 'refraction': False}
    get = get_da_for_transect
    for time, lat, lon in zip(*space_time_transect):
        sites = sunwhere.sites(time, lat, lon, **kwargs)
        trans = sunwhere.transect(time, lat, lon, **kwargs)
        assert allclose(get(sites.sza), get(trans.sza), atol=2., rtol=0.)
        assert allclose(get(sites.dec), get(trans.dec), atol=2., rtol=0.)
        assert allclose(get(sites.eot), get(trans.eot), atol=10., rtol=0.)
        assert allclose(get(sites.ecf), get(trans.ecf), atol=0., rtol=0.01)
        assert azimdist(get(sites.saa), get(trans.saa), thresh=1.)


@pytest.mark.filterwarnings("ignore")
def test_soltrack_numexpr_transect_consistency(space_time_transect, allclose):
    kwargs = {'algorithm': 'soltrack', 'engine': 'numexpr', 'refraction': False}
    get = get_da_for_transect
    for time, lat, lon in zip(*space_time_transect):
        sites = sunwhere.sites(time, lat, lon, **kwargs)
        trans = sunwhere.transect(time, lat, lon, **kwargs)
        assert allclose(get(sites.sza), get(trans.sza), atol=2., rtol=0.)
        assert allclose(get(sites.dec), get(trans.dec), atol=2., rtol=0.)
        assert allclose(get(sites.eot), get(trans.eot), atol=10., rtol=0.)
        assert allclose(get(sites.ecf), get(trans.ecf), atol=0., rtol=0.01)
        assert azimdist(get(sites.saa), get(trans.saa), thresh=1.)


@pytest.mark.filterwarnings("ignore")
def test_iqbal_numpy_transect_consistency(space_time_transect, allclose):
    kwargs = {'algorithm': 'iqbal', 'engine': 'numpy', 'refraction': False}
    get = get_da_for_transect
    for time, lat, lon in zip(*space_time_transect):
        sites = sunwhere.sites(time, lat, lon, **kwargs)
        trans = sunwhere.transect(time, lat, lon, **kwargs)
        assert allclose(get(sites.sza), get(trans.sza), atol=2., rtol=0.)
        assert allclose(get(sites.dec), get(trans.dec), atol=2., rtol=0.)
        assert allclose(get(sites.eot), get(trans.eot), atol=10., rtol=0.)
        assert allclose(get(sites.ecf), get(trans.ecf), atol=0., rtol=0.01)
        assert azimdist(get(sites.saa), get(trans.saa), thresh=1.)


@pytest.mark.filterwarnings("ignore")
def test_iqbal_numexpr_transect_consistency(space_time_transect, allclose):
    kwargs = {'algorithm': 'iqbal', 'engine': 'numexpr', 'refraction': False}
    get = get_da_for_transect
    for time, lat, lon in zip(*space_time_transect):
        sites = sunwhere.sites(time, lat, lon, **kwargs)
        trans = sunwhere.transect(time, lat, lon, **kwargs)
        assert allclose(get(sites.sza), get(trans.sza), atol=2., rtol=0.)
        assert allclose(get(sites.dec), get(trans.dec), atol=2., rtol=0.)
        assert allclose(get(sites.eot), get(trans.eot), atol=10., rtol=0.)
        assert allclose(get(sites.ecf), get(trans.ecf), atol=0., rtol=0.01)
        assert azimdist(get(sites.saa), get(trans.saa), thresh=1.)
