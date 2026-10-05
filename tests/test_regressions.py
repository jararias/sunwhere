"""Regression tests for bugs fixed in v1.5.0."""

import numpy as np
import pandas as pd
import pytest
import xarray as xr

import sunwhere


@pytest.fixture(params=["UTC", "Europe/Madrid", None])
def times(request):
    return pd.date_range("2025-01-01", "2025-12-31 23:00", freq="h", tz=request.param)


class TestTimezoneAwareTimes:
    """Outputs of tz-aware inputs must work with the whole xarray API."""

    def test_time_coordinate_is_naive_utc(self, times):
        result = sunwhere.sites(times, [53.14, 40.0], [8.21, -3.0])
        expected = pd.to_datetime(times, utc=True).tz_localize(None)
        assert result.sza.indexes["time"].tz is None
        assert result.sza.time.dtype == np.dtype("datetime64[ns]")
        np.testing.assert_array_equal(result.sza.indexes["time"], expected)

    def test_times_property_keeps_timezone(self, times):
        result = sunwhere.sites(times, 40.0, -3.0)
        assert result.times.tz == times.tz

    @pytest.mark.parametrize("usecase", ["sites", "regular_grid", "transect"])
    def test_xarray_operations(self, times, usecase):
        n = len(times)
        if usecase == "transect":
            result = sunwhere.transect(times, np.zeros(n), np.zeros(n))
        else:
            result = getattr(sunwhere, usecase)(times, [10.0, 20.0], [0.0, 5.0])
        sza = result.sza
        sza.groupby("time.month").mean()
        sza.resample(time="D").mean()
        sza.interp(time=sza.time[:3] + np.timedelta64(30, "m"))
        sza.sel(time=sza.time[5].values)
        xr.Dataset({"sza": sza, "saa": result.saa}).to_netcdf()
        naive = xr.DataArray(np.ones(n), coords={"time": sza.time.values}, dims="time")
        assert (sza + naive).sizes["time"] == n

    def test_same_results_regardless_of_timezone(self):
        utc = pd.date_range("2025-06-21", periods=24, freq="h", tz="UTC")
        ref = sunwhere.sites(utc, 40.0, -3.0).sza
        for tz in ("Europe/Madrid", "America/New_York"):
            other = sunwhere.sites(utc.tz_convert(tz), 40.0, -3.0).sza
            xr.testing.assert_allclose(ref, other)


class TestSubsecondTimes:
    @pytest.mark.parametrize("usecase", ["sites", "regular_grid", "transect"])
    def test_subsecond_resolution_is_kept(self, usecase):
        times = pd.to_datetime(["2025-06-21 12:00:00.0", "2025-06-21 12:00:00.5"])
        if usecase == "transect":
            result = sunwhere.transect(times, [0.0, 0.0], [0.0, 0.0])
        else:
            result = getattr(sunwhere, usecase)(times, 0.0, 0.0)
        assert np.diff(result.sza.values.ravel())[0] != 0.0


class TestSunriseSunset:
    @pytest.mark.parametrize("units", ["deg", "rad"])
    def test_sunset_angle_is_minus_sunrise_angle(self, units):
        times = pd.date_range("2025-06-21", periods=24, freq="h")
        result = sunwhere.sites(times, [40.0, -33.9], [-3.0, 151.2])
        xr.testing.assert_allclose(result.sunset(units), -result.sunrise(units))

    def test_sunset_after_sunrise(self):
        times = pd.date_range("2025-06-21", periods=24, freq="h")
        result = sunwhere.sites(times, 40.0, -3.0)
        assert (result.sunset("utc") > result.sunrise("utc")).all()


class TestAzimuth:
    def test_azimuth_range_zero_south(self):
        times = pd.date_range("2025-01-01", "2025-12-31", freq="h")
        saa = sunwhere.sites(times, 53.14, 8.21).saa
        assert saa.min() >= -180.0 and saa.max() <= 180.0
        assert saa.min() < -90.0 and saa.max() > 90.0

    def test_azimuth_north(self):
        times = pd.date_range("2025-06-21", periods=24, freq="h")
        result = sunwhere.sites(times, 40.0, -3.0)
        az_n = result.azimuth_north
        assert az_n.min() >= 0.0 and az_n.max() < 360.0
        xr.testing.assert_allclose((az_n - 180.0).rename("azimuth"), result.azimuth,
                                   check_dim_order=True)

    def test_azimuth_north_matches_pvlib(self):
        pvlib = pytest.importorskip("pvlib")
        times = pd.date_range("2025-06-21", periods=24, freq="h", tz="UTC")
        sw = sunwhere.sites(times, 40.0, -3.0, algorithm="nrel")
        pv = pvlib.solarposition.get_solarposition(times, 40.0, -3.0, method="nrel_numpy")
        diff = (sw.azimuth_north.values[:, 0] - pv["azimuth"].values + 180) % 360 - 180
        assert np.abs(diff).max() < 0.01


class TestMisc:
    def test_incidence_dims_match_zenith(self):
        times = pd.date_range("2025-06-21", periods=24, freq="h")
        result = sunwhere.sites(times, [40.0, 50.0], [-3.0, 8.0])
        assert result.incidence(30.0, 0.0).dims == result.zenith.dims

    def test_incidence_horizontal_equals_cosz(self):
        times = pd.date_range("2025-06-21", periods=24, freq="h")
        # without refraction, as incidence() uses the geometric hour angle
        result = sunwhere.sites(times, [40.0, 50.0], [-3.0, 8.0], refraction=False)
        np.testing.assert_allclose(result.incidence(0.0, 0.0), result.cosz, atol=1e-3)

    def test_universal_time_coordinated(self):
        times_utc = pd.date_range("2025-06-21", periods=24, freq="h")
        result = sunwhere.sites(times_utc, 40.0, -3.0)
        tst = result.true_solar_time.values[:, 0]
        back = sunwhere.universal_time_coordinated(tst, -3.0)
        err = np.abs((back - times_utc.values).astype("timedelta64[ns]").astype(float))
        assert err.max() < 1e9  # less than one second


class TestTransectScalars:
    def test_scalar_longitude_is_broadcast(self):
        times = pd.date_range("2024-06-21 12:00", periods=181, freq="h")
        lats = np.linspace(-90, 90, 181)
        result = sunwhere.transect(times, lats, 0.0)
        assert result.sza.shape == (181,)
        np.testing.assert_array_equal(result.sza.lon.values, 0.0)

    def test_scalar_latitude_and_longitude(self):
        times = pd.date_range("2024-06-21", periods=24, freq="h")
        transect = sunwhere.transect(times, 40.0, -3.0).sza.values
        sites = sunwhere.sites(times, 40.0, -3.0).sza.values[:, 0]
        np.testing.assert_allclose(transect, sites)


def test_public_algorithms():
    assert set(sunwhere.SPA_ALGORITHMS) == {"psa", "iqbal", "nrel"}
