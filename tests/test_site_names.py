"""Tests for site_names parameter in sites() function"""

import numpy as np
import pandas as pd
import pytest

import sunwhere


class TestSiteNames:
    """Test suite for site_names feature in sites() function"""

    def test_sites_with_names(self):
        """Test sites() with custom site names"""
        times = pd.date_range('2024-06-21', periods=24, freq='h', tz='UTC')
        lats = [40.4168, 40.7128, 35.6762]
        lons = [-3.7038, -74.0060, 139.6503]
        names = ['Madrid', 'New York', 'Tokyo']

        result = sunwhere.sites(times, lats, lons, site_names=names)

        # Check that site coordinate has the correct names
        assert list(result.sza.coords['site'].values) == names
        assert list(result.azimuth.coords['site'].values) == names

        # Check that latitude and longitude arrays have the correct coords
        assert list(result.latitude.coords['site'].values) == names
        assert list(result.longitude.coords['site'].values) == names

    def test_sites_without_names(self):
        """Test sites() without site_names (backward compatibility)"""
        times = pd.date_range('2024-06-21', periods=24, freq='h', tz='UTC')
        lats = [40.4168, 40.7128, 35.6762]
        lons = [-3.7038, -74.0060, 139.6503]

        result = sunwhere.sites(times, lats, lons)

        # Check that site coordinate uses numeric indices
        assert list(result.sza.coords['site'].values) == [0, 1, 2]

    def test_site_selection_by_name(self):
        """Test selection of data by site name"""
        times = pd.date_range('2024-06-21', periods=24, freq='h', tz='UTC')
        lats = [40.4168, 40.7128]
        lons = [-3.7038, -74.0060]
        names = ['Madrid', 'NYC']

        result = sunwhere.sites(times, lats, lons, site_names=names)

        # Select by site name
        madrid_data = result.sza.sel(site='Madrid')
        nyc_data = result.sza.sel(site='NYC')

        assert madrid_data.shape == (24,)
        assert nyc_data.shape == (24,)

        # Verify they are different
        assert not np.allclose(madrid_data.values, nyc_data.values)

    def test_single_site_with_name(self):
        """Test single site with custom name"""
        times = pd.date_range('2024-06-21', periods=24, freq='h', tz='UTC')
        sunpos = sunwhere.sites(times, 40.0, -3.0, site_names=['MyLocation'])

        assert list(sunpos.sza.coords['site'].values) == ['MyLocation']
        assert sunpos.sza.shape == (24, 1)

    def test_site_names_length_mismatch(self):
        """Test that mismatched site_names length raises ValueError"""
        times = pd.date_range('2024-06-21', periods=24, freq='h', tz='UTC')
        lats = [40.0, 41.0, 42.0]
        lons = [-3.0, -4.0, -5.0]
        names = ['Site1', 'Site2']  # Only 2 names for 3 sites

        with pytest.raises(ValueError, match='site_names must have the same length'):
            sunwhere.sites(times, lats, lons, site_names=names)

    def test_site_names_as_list(self):
        """Test site_names provided as list"""
        times = pd.date_range('2024-06-21', periods=24, freq='h', tz='UTC')
        lats = [40.0, 41.0]
        lons = [-3.0, -4.0]
        names = ['Site1', 'Site2']

        result = sunwhere.sites(times, lats, lons, site_names=names)
        assert list(result.sza.coords['site'].values) == names

    def test_site_names_as_numpy_array(self):
        """Test site_names provided as numpy array"""
        times = pd.date_range('2024-06-21', periods=24, freq='h', tz='UTC')
        lats = [40.0, 41.0]
        lons = [-3.0, -4.0]
        names = np.array(['Site1', 'Site2'])

        result = sunwhere.sites(times, lats, lons, site_names=names)
        assert list(result.sza.coords['site'].values) == ['Site1', 'Site2']

    def test_site_names_with_unicode(self):
        """Test site_names with Unicode characters"""
        times = pd.date_range('2024-06-21', periods=24, freq='h', tz='UTC')
        lats = [40.4168, 35.6762]
        lons = [-3.7038, 139.6503]
        names = ['Madrid', '東京']  # Tokyo in Japanese

        result = sunwhere.sites(times, lats, lons, site_names=names)
        assert list(result.sza.coords['site'].values) == names

    def test_site_names_with_spaces(self):
        """Test site_names with spaces"""
        times = pd.date_range('2024-06-21', periods=24, freq='h', tz='UTC')
        lats = [40.7128, 51.5074]
        lons = [-74.0060, -0.1278]
        names = ['New York City', 'Greater London']

        result = sunwhere.sites(times, lats, lons, site_names=names)
        assert list(result.sza.coords['site'].values) == names

    def test_site_names_None_explicit(self):
        """Test site_names=None explicitly"""
        times = pd.date_range('2024-06-21', periods=24, freq='h', tz='UTC')
        lats = [40.0, 41.0]
        lons = [-3.0, -4.0]

        result = sunwhere.sites(times, lats, lons, site_names=None)
        assert list(result.sza.coords['site'].values) == [0, 1]

    def test_all_algorithms_with_site_names(self):
        """Test that site_names works with all algorithms"""
        times = pd.date_range('2024-06-21', periods=5, freq='h', tz='UTC')
        lats = [40.0]
        lons = [-3.0]
        names = ['TestSite']

        for algorithm in ['psa', 'nrel', 'iqbal']:
            result = sunwhere.sites(times, lats, lons, 
                                   algorithm=algorithm, site_names=names)
            assert list(result.sza.coords['site'].values) == names
