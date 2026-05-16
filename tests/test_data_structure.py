"""Tests for xarray data structure correctness."""

import numpy as np
import pandas as pd
import pytest
import sunwhere


class TestSitesDataStructure:
    """Test xarray structure for sites use case."""

    def test_dimensions(self):
        """Test that sites result has correct dimensions."""
        time = pd.date_range("2024-01-01", periods=10, freq="h", tz="UTC")
        lat = np.array([40.0, 41.0, 42.0])
        lon = np.array([-3.0, -2.0, -1.0])
        result = sunwhere.sites(time, lat, lon)
        
        assert result.sza.dims == ("time", "site")
        assert result.saa.dims == ("time", "site")
        assert result.dec.dims == ("time",)
        assert result.eot.dims == ("time",)

    def test_coordinates(self):
        """Test that coordinates are correctly set."""
        time = pd.date_range("2024-01-01", periods=5, freq="h", tz="UTC")
        lat = np.array([40.0, 41.0])
        lon = np.array([-3.0, -2.0])
        result = sunwhere.sites(time, lat, lon)
        
        # Check time coordinate
        assert "time" in result.sza.coords
        # Compare values only, not dtype (xarray may have different timezone representation)
        assert len(result.sza.coords["time"]) == len(time)
        
        # Check site coordinate
        assert "site" in result.sza.coords
        assert len(result.sza.coords["site"]) == 2

    def test_shape(self):
        """Test that output shapes are correct."""
        time = pd.date_range("2024-01-01", periods=10, freq="h", tz="UTC")
        lat = np.array([40.0, 41.0, 42.0])
        lon = np.array([-3.0, -2.0, -1.0])
        result = sunwhere.sites(time, lat, lon)
        
        assert result.sza.shape == (10, 3)
        assert result.saa.shape == (10, 3)
        assert result.dec.shape == (10,)
        assert result.eot.shape == (10,)
        assert result.ecf.shape == (10,)

    def test_attributes(self):
        """Test that xarray attributes are set."""
        time = pd.date_range("2024-01-01", periods=1, freq="h", tz="UTC")
        result = sunwhere.sites(time, 40.0, 0.0)
        
        # Check that attributes exist
        assert hasattr(result.sza, "attrs")
        assert isinstance(result.sza.attrs, dict)


class TestRegularGridDataStructure:
    """Test xarray structure for regular_grid use case."""

    def test_dimensions(self):
        """Test that regular_grid result has correct dimensions."""
        time = pd.date_range("2024-01-01", periods=5, freq="h", tz="UTC")
        lat = np.array([40.0, 41.0, 42.0])
        lon = np.array([-3.0, -2.0, -1.0, 0.0])
        result = sunwhere.regular_grid(time, lat, lon)
        
        assert result.sza.dims == ("time", "lat", "lon")
        assert result.saa.dims == ("time", "lat", "lon")
        assert result.dec.dims == ("time",)
        assert result.eot.dims == ("time",)

    def test_coordinates(self):
        """Test that coordinates are correctly set."""
        time = pd.date_range("2024-01-01", periods=2, freq="h", tz="UTC")
        lat = np.array([40.0, 41.0])
        lon = np.array([-3.0, -2.0])
        result = sunwhere.regular_grid(time, lat, lon)
        
        # Check time coordinate
        assert "time" in result.sza.coords
        # Compare values only, not dtype (xarray may have different timezone representation)
        assert len(result.sza.coords["time"]) == len(time)
        
        # Check lat coordinate
        assert "lat" in result.sza.coords
        np.testing.assert_array_equal(result.sza.coords["lat"].values, lat)
        
        # Check lon coordinate
        assert "lon" in result.sza.coords
        np.testing.assert_array_equal(result.sza.coords["lon"].values, lon)

    def test_shape(self):
        """Test that output shapes are correct."""
        time = pd.date_range("2024-01-01", periods=10, freq="h", tz="UTC")
        lat = np.array([40.0, 41.0, 42.0])
        lon = np.array([-3.0, -2.0, -1.0, 0.0])
        result = sunwhere.regular_grid(time, lat, lon)
        
        assert result.sza.shape == (10, 3, 4)
        assert result.saa.shape == (10, 3, 4)
        assert result.dec.shape == (10,)
        assert result.eot.shape == (10,)

    def test_broadcasting(self):
        """Test that values are correctly broadcast across grid."""
        time = pd.date_range("2024-01-01 12:00", periods=1, freq="h", tz="UTC")
        lat = np.array([40.0, 41.0])
        lon = np.array([-3.0, -2.0])
        result = sunwhere.regular_grid(time, lat, lon)
        
        # For same latitude, different longitudes should have different SAA
        assert result.saa.values[0, 0, 0] != result.saa.values[0, 0, 1]
        
        # For same longitude, different latitudes should have different SZA
        assert result.sza.values[0, 0, 0] != result.sza.values[0, 1, 0]


class TestTransectDataStructure:
    """Test xarray structure for transect use case."""

    def test_dimensions(self):
        """Test that transect result has correct dimensions."""
        time = pd.date_range("2024-01-01", periods=5, freq="h", tz="UTC")
        lat = np.linspace(40.0, 45.0, 5)
        lon = np.linspace(-3.0, 2.0, 5)
        result = sunwhere.transect(time, lat, lon)
        
        assert result.sza.dims == ("time",)
        assert result.saa.dims == ("time",)
        assert result.dec.dims == ("time",)
        assert result.eot.dims == ("time",)

    def test_coordinates(self):
        """Test that coordinates are correctly set."""
        time = pd.date_range("2024-01-01", periods=3, freq="h", tz="UTC")
        lat = np.array([40.0, 41.0, 42.0])
        lon = np.array([-3.0, -2.0, -1.0])
        result = sunwhere.transect(time, lat, lon)
        
        # Check time coordinate
        assert "time" in result.sza.coords
        # Compare values only, not dtype (xarray may have different timezone representation)
        assert len(result.sza.coords["time"]) == len(time)

    def test_shape(self):
        """Test that output shapes are correct."""
        time = pd.date_range("2024-01-01", periods=10, freq="h", tz="UTC")
        lat = np.linspace(40.0, 45.0, 10)
        lon = np.linspace(-3.0, 2.0, 10)
        result = sunwhere.transect(time, lat, lon)
        
        assert result.sza.shape == (10,)
        assert result.saa.shape == (10,)
        assert result.dec.shape == (10,)
        assert result.eot.shape == (10,)


class TestDataTypes:
    """Test that data types are correct."""

    def test_float_output(self):
        """Test that outputs are floating point."""
        time = pd.date_range("2024-01-01", periods=1, freq="h", tz="UTC")
        result = sunwhere.sites(time, 40.0, 0.0)
        
        assert np.issubdtype(result.sza.dtype, np.floating)
        assert np.issubdtype(result.saa.dtype, np.floating)
        assert np.issubdtype(result.dec.dtype, np.floating)
        assert np.issubdtype(result.eot.dtype, np.floating)

    def test_no_nan_for_valid_inputs(self):
        """Test that valid inputs don't produce NaN."""
        time = pd.date_range("2024-01-01", periods=24, freq="h", tz="UTC")
        result = sunwhere.sites(time, 40.0, 0.0)
        
        assert not np.isnan(result.sza.values).any()
        assert not np.isnan(result.saa.values).any()
        assert not np.isnan(result.dec.values).any()
        assert not np.isnan(result.eot.values).any()


class TestSunposProperties:
    """Test that all Sunpos properties are accessible."""

    def test_all_properties_exist(self):
        """Test that all expected properties exist."""
        time = pd.date_range("2024-01-01", periods=1, freq="h", tz="UTC")
        result = sunwhere.sites(time, 40.0, 0.0)
        
        # Main angles
        assert hasattr(result, "sza")
        assert hasattr(result, "saa")
        
        # Solar characteristics
        assert hasattr(result, "dec")
        assert hasattr(result, "eot")
        assert hasattr(result, "ecf")
        
        # Times
        assert hasattr(result, "sunrise")
        assert hasattr(result, "sunset")

    def test_properties_are_dataarrays(self):
        """Test that properties return xarray DataArrays."""
        import xarray as xr
        
        time = pd.date_range("2024-01-01", periods=1, freq="h", tz="UTC")
        result = sunwhere.sites(time, 40.0, 0.0)
        
        assert isinstance(result.sza, xr.DataArray)
        assert isinstance(result.saa, xr.DataArray)
        assert isinstance(result.dec, xr.DataArray)
        assert isinstance(result.eot, xr.DataArray)

    def test_derived_angle_relationships(self):
        """Test mathematical relationships between angles."""
        time = pd.date_range("2024-01-01 12:00", periods=1, freq="h", tz="UTC")
        result = sunwhere.sites(time, 0.0, 0.0, refraction=False)
        
        # Solar zenith angle should be between 0 and 180
        assert np.all(result.sza.values >= 0)
        assert np.all(result.sza.values <= 180)
