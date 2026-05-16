"""Tests for edge cases and boundary conditions."""

import numpy as np
import pandas as pd
import sunwhere


class TestGeographicBoundaries:
    """Test calculations at geographic boundaries."""

    def test_north_pole(self):
        """Test solar position at North Pole."""
        time = pd.date_range("2024-06-21 12:00", periods=1, freq="h", tz="UTC")
        result = sunwhere.sites(time, 90.0, 0.0)
        assert result.sza is not None
        assert result.saa is not None
        assert not np.isnan(result.sza.values).any()

    def test_south_pole(self):
        """Test solar position at South Pole."""
        time = pd.date_range("2024-12-21 12:00", periods=1, freq="h", tz="UTC")
        result = sunwhere.sites(time, -90.0, 0.0)
        assert result.sza is not None
        assert result.saa is not None
        assert not np.isnan(result.sza.values).any()

    def test_equator(self):
        """Test solar position at equator."""
        time = pd.date_range("2024-03-20 12:00", periods=1, freq="h", tz="UTC")
        result = sunwhere.sites(time, 0.0, 0.0)
        assert result.sza is not None
        assert result.saa is not None
        # At equinox, solar zenith angle should be close to 0 at noon
        assert result.sza.values[0] < 10.0  # Within 10 degrees

    def test_antimeridian(self):
        """Test solar position at antimeridian (±180°)."""
        time = pd.date_range("2024-01-01 12:00", periods=1, freq="h", tz="UTC")
        # Longitude must be < 180, so use 179.99 instead
        result_pos = sunwhere.sites(time, 0.0, 179.99)
        result_neg = sunwhere.sites(time, 0.0, -180.0)
        # Results should be very close near antimeridian
        assert np.allclose(result_pos.sza.values, result_neg.sza.values, atol=0.1)

    def test_prime_meridian(self):
        """Test solar position at prime meridian (0°)."""
        time = pd.date_range("2024-01-01 12:00", periods=1, freq="h", tz="UTC")
        result = sunwhere.sites(time, 0.0, 0.0)
        assert result.sza is not None
        assert not np.isnan(result.sza.values).any()


class TestTemporalBoundaries:
    """Test calculations at temporal boundaries."""

    def test_summer_solstice(self):
        """Test solar position at summer solstice."""
        time = pd.date_range("2024-06-21", periods=24, freq="h", tz="UTC")
        result = sunwhere.sites(time, 40.0, 0.0)
        # Declination should be close to 23.44°
        assert np.max(result.dec.values) > 23.0
        assert np.max(result.dec.values) < 24.0

    def test_winter_solstice(self):
        """Test solar position at winter solstice."""
        time = pd.date_range("2024-12-21", periods=24, freq="h", tz="UTC")
        result = sunwhere.sites(time, 40.0, 0.0)
        # Declination should be close to -23.44°
        assert np.min(result.dec.values) < -23.0
        assert np.min(result.dec.values) > -24.0

    def test_vernal_equinox(self):
        """Test solar position at vernal equinox."""
        time = pd.date_range("2024-03-20", periods=24, freq="h", tz="UTC")
        result = sunwhere.sites(time, 0.0, 0.0)
        # Declination should be close to 0°
        assert np.abs(np.mean(result.dec.values)) < 2.0

    def test_autumnal_equinox(self):
        """Test solar position at autumnal equinox."""
        time = pd.date_range("2024-09-22", periods=24, freq="h", tz="UTC")
        result = sunwhere.sites(time, 0.0, 0.0)
        # Declination should be close to 0°
        assert np.abs(np.mean(result.dec.values)) < 2.0

    def test_midnight_sun(self):
        """Test midnight sun phenomenon in Arctic summer."""
        # North of Arctic Circle on summer solstice
        time = pd.date_range("2024-06-21", periods=24, freq="h", tz="UTC")
        result = sunwhere.sites(time, 70.0, 0.0)
        # Sun should stay above horizon all day
        assert np.all(result.sza.values < 90.0)

    def test_polar_night(self):
        """Test polar night phenomenon in Arctic winter."""
        # North of Arctic Circle on winter solstice
        time = pd.date_range("2024-12-21", periods=24, freq="h", tz="UTC")
        result = sunwhere.sites(time, 70.0, 0.0)
        # Sun should stay below horizon all day
        assert np.all(result.sza.values > 90.0)

    def test_year_2000(self):
        """Test leap year 2000."""
        time = pd.date_range("2000-02-29", periods=1, freq="h", tz="UTC")
        result = sunwhere.sites(time, 40.0, 0.0)
        assert result.sza is not None
        assert not np.isnan(result.sza.values).any()

    def test_year_1900(self):
        """Test non-leap year 1900."""
        time = pd.date_range("1900-03-01", periods=1, freq="h", tz="UTC")
        result = sunwhere.sites(time, 40.0, 0.0)
        assert result.sza is not None
        assert not np.isnan(result.sza.values).any()

    def test_far_future(self):
        """Test calculation for year 2100."""
        time = pd.date_range("2100-01-01", periods=1, freq="h", tz="UTC")
        result = sunwhere.sites(time, 40.0, 0.0)
        assert result.sza is not None
        assert not np.isnan(result.sza.values).any()


class TestTimezones:
    """Test calculations with different timezones."""

    def test_utc_timezone(self):
        """Test with UTC timezone."""
        time = pd.date_range("2024-01-01 12:00", periods=1, freq="h", tz="UTC")
        result = sunwhere.sites(time, 0.0, 0.0)
        assert result.sza is not None

    def test_local_timezone(self):
        """Test with local timezone (Europe/Madrid)."""
        time = pd.date_range("2024-01-01 12:00", periods=1, freq="h", tz="Europe/Madrid")
        result = sunwhere.sites(time, 40.0, -3.0)
        assert result.sza is not None

    def test_us_timezone(self):
        """Test with US timezone (America/New_York)."""
        time = pd.date_range("2024-01-01 12:00", periods=1, freq="h", tz="America/New_York")
        result = sunwhere.sites(time, 40.0, -74.0)
        assert result.sza is not None

    def test_asia_timezone(self):
        """Test with Asia timezone (Asia/Tokyo)."""
        time = pd.date_range("2024-01-01 12:00", periods=1, freq="h", tz="Asia/Tokyo")
        result = sunwhere.sites(time, 35.0, 139.0)
        assert result.sza is not None


class TestNumericalStability:
    """Test numerical stability with extreme inputs."""

    def test_very_high_temporal_resolution(self):
        """Test with very high temporal resolution (1 minute)."""
        time = pd.date_range("2024-01-01", periods=60, freq="min", tz="UTC")
        result = sunwhere.sites(time, 40.0, 0.0)
        assert not np.isnan(result.sza.values).any()
        # Check temporal continuity
        sza_diff = np.diff(result.sza.values)
        assert np.all(np.abs(sza_diff) < 1.0)  # Change less than 1° per minute

    def test_very_low_temporal_resolution(self):
        """Test with very low temporal resolution (1 month)."""
        time = pd.date_range("2024-01-01", periods=12, freq="MS", tz="UTC")
        result = sunwhere.sites(time, 40.0, 0.0)
        assert not np.isnan(result.sza.values).any()

    def test_large_number_of_sites(self):
        """Test with many sites (1000 sites)."""
        time = pd.date_range("2024-01-01", periods=1, freq="h", tz="UTC")
        lat = np.random.uniform(-90, 90, 1000)
        lon = np.random.uniform(-180, 180, 1000)
        result = sunwhere.sites(time, lat, lon)
        assert result.sza.shape == (1, 1000)
        assert not np.isnan(result.sza.values).any()

    def test_large_regular_grid(self):
        """Test with large regular grid (100x100)."""
        time = pd.date_range("2024-01-01", periods=1, freq="h", tz="UTC")
        lat = np.linspace(-90, 90, 100)
        # Longitude must be < 180, so use endpoint=False
        lon = np.linspace(-180, 180, 100, endpoint=False)
        result = sunwhere.regular_grid(time, lat, lon)
        assert result.sza.shape == (1, 100, 100)
        assert not np.isnan(result.sza.values).any()


class TestRefractionEffects:
    """Test refraction parameter."""

    def test_refraction_enabled(self):
        """Test with refraction enabled."""
        time = pd.date_range("2024-01-01 06:00", periods=1, freq="h", tz="UTC")
        result = sunwhere.sites(time, 40.0, 0.0, refraction=True)
        assert result.sza is not None

    def test_refraction_disabled(self):
        """Test with refraction disabled."""
        time = pd.date_range("2024-01-01 06:00", periods=1, freq="h", tz="UTC")
        result = sunwhere.sites(time, 40.0, 0.0, refraction=False)
        assert result.sza is not None

    def test_refraction_difference(self):
        """Test that refraction makes a difference near horizon."""
        time = pd.date_range("2024-01-01 06:00", periods=1, freq="h", tz="UTC")
        result_with = sunwhere.sites(time, 40.0, 0.0, refraction=True)
        result_without = sunwhere.sites(time, 40.0, 0.0, refraction=False)
        # Refraction should reduce zenith angle
        assert result_with.sza.values[0] <= result_without.sza.values[0]
