"""Tests for input validation."""

import numpy as np
import pandas as pd
import pytest
import sunwhere


class TestSitesValidation:
    """Test input validation for sites use case."""

    def test_invalid_latitude_range(self):
        """Test that latitude outside [-90, 90] raises ValueError."""
        time = pd.date_range("2024-01-01", periods=1, freq="h", tz="UTC")
        with pytest.raises(ValueError, match="latitude out of bounds"):
            sunwhere.sites(time, 95.0, 0.0)
        with pytest.raises(ValueError, match="latitude out of bounds"):
            sunwhere.sites(time, -95.0, 0.0)

    def test_invalid_longitude_range(self):
        """Test that longitude outside [-180, 180] raises ValueError."""
        time = pd.date_range("2024-01-01", periods=1, freq="h", tz="UTC")
        with pytest.raises(ValueError, match="longitude out of bounds"):
            sunwhere.sites(time, 0.0, 185.0)
        with pytest.raises(ValueError, match="longitude out of bounds"):
            sunwhere.sites(time, 0.0, -185.0)

    def test_invalid_latitude_type(self):
        """Test that latitude with wrong dimensionality raises ValueError."""
        time = pd.date_range("2024-01-01", periods=1, freq="h", tz="UTC")
        with pytest.raises(ValueError, match="wrong dimensions"):
            sunwhere.sites(time, [[40.0]], 0.0)

    def test_invalid_longitude_type(self):
        """Test that longitude with wrong dimensionality raises ValueError."""
        time = pd.date_range("2024-01-01", periods=1, freq="h", tz="UTC")
        with pytest.raises(ValueError, match="wrong dimensions"):
            sunwhere.sites(time, 0.0, [[0.0]])

    def test_latitude_longitude_shape_mismatch(self):
        """Test that mismatched latitude/longitude shapes raise ValueError."""
        time = pd.date_range("2024-01-01", periods=1, freq="h", tz="UTC")
        lat = np.array([40.0, 41.0])
        lon = np.array([0.0, 1.0, 2.0])
        with pytest.raises(ValueError, match="shape mismatch"):
            sunwhere.sites(time, lat, lon)

    def test_invalid_algorithm(self):
        """Test that invalid algorithm name raises ValueError."""
        time = pd.date_range("2024-01-01", periods=1, freq="h", tz="UTC")
        with pytest.raises(ValueError, match="missing algorithm"):
            sunwhere.sites(time, 0.0, 0.0, algorithm="invalid")

    def test_invalid_engine(self):
        """Test that invalid engine name raises ValueError."""
        time = pd.date_range("2024-01-01", periods=1, freq="h", tz="UTC")
        with pytest.raises(ValueError, match="missing engine"):
            sunwhere.sites(time, 0.0, 0.0, engine="invalid")


class TestRegularGridValidation:
    """Test input validation for regular_grid use case."""

    def test_latitude_not_1d(self):
        """Test that 2D latitude raises ValueError."""
        time = pd.date_range("2024-01-01", periods=1, freq="h", tz="UTC")
        lat = np.array([[40.0, 41.0], [42.0, 43.0]])
        lon = np.array([0.0, 1.0])
        with pytest.raises(ValueError, match="wrong dimensions"):
            sunwhere.regular_grid(time, lat, lon)

    def test_longitude_not_1d(self):
        """Test that 2D longitude raises ValueError."""
        time = pd.date_range("2024-01-01", periods=1, freq="h", tz="UTC")
        lat = np.array([40.0, 41.0])
        lon = np.array([[0.0, 1.0], [2.0, 3.0]])
        with pytest.raises(ValueError, match="wrong dimensions"):
            sunwhere.regular_grid(time, lat, lon)

    def test_invalid_latitude_range(self):
        """Test that latitude outside [-90, 90] raises ValueError."""
        time = pd.date_range("2024-01-01", periods=1, freq="h", tz="UTC")
        with pytest.raises(ValueError, match="latitude out of bounds"):
            sunwhere.regular_grid(time, np.array([95.0]), np.array([0.0]))

    def test_invalid_longitude_range(self):
        """Test that longitude outside [-180, 180] raises ValueError."""
        time = pd.date_range("2024-01-01", periods=1, freq="h", tz="UTC")
        with pytest.raises(ValueError, match="longitude out of bounds"):
            sunwhere.regular_grid(time, np.array([0.0]), np.array([185.0]))


class TestTransectValidation:
    """Test input validation for transect use case."""

    def test_time_latitude_shape_mismatch(self):
        """Test that mismatched time/latitude shapes raise ValueError."""
        time = pd.date_range("2024-01-01", periods=2, freq="h", tz="UTC")
        lat = np.array([40.0, 41.0, 42.0])
        lon = np.array([0.0, 1.0])
        with pytest.raises(ValueError, match="shape mismatch"):
            sunwhere.transect(time, lat, lon)

    def test_time_longitude_shape_mismatch(self):
        """Test that mismatched time/longitude shapes raise ValueError."""
        time = pd.date_range("2024-01-01", periods=2, freq="h", tz="UTC")
        lat = np.array([40.0, 41.0])
        lon = np.array([0.0, 1.0, 2.0])
        with pytest.raises(ValueError, match="shape mismatch"):
            sunwhere.transect(time, lat, lon)

    def test_invalid_latitude_range(self):
        """Test that latitude outside [-90, 90] raises ValueError."""
        time = pd.date_range("2024-01-01", periods=1, freq="h", tz="UTC")
        with pytest.raises(ValueError, match="latitude out of bounds"):
            sunwhere.transect(time, np.array([95.0]), np.array([0.0]))

    def test_invalid_longitude_range(self):
        """Test that longitude outside [-180, 180] raises ValueError."""
        time = pd.date_range("2024-01-01", periods=1, freq="h", tz="UTC")
        with pytest.raises(ValueError, match="longitude out of bounds"):
            sunwhere.transect(time, np.array([0.0]), np.array([185.0]))
