# User Guide

This guide covers the three main use cases of Sunwhere and shows you how to calculate solar positions for your specific application.

## Core Concepts

### Time Handling

Sunwhere accepts various time formats:

- Python `datetime.datetime` objects
- Pandas `DatetimeIndex`
- NumPy `datetime64` arrays
- Any format accepted by `pandas.to_datetime()`

**Important**: Timezone-naive times are assumed to be UTC.

```python
import pandas as pd
from datetime import datetime

# All of these work:
times1 = pd.date_range('2024-01-01', periods=24, freq='h', tz='UTC')
times2 = pd.to_datetime(['2024-01-01 12:00', '2024-01-01 13:00'], utc=True)
times3 = [datetime(2024, 1, 1, 12, 0), datetime(2024, 1, 1, 13, 0)]
```

### Geographic Coordinates

- **Latitude**: -90° (South Pole) to +90° (North Pole)
- **Longitude**: -180° to +180° (exclusive), positive east

```python
# Valid coordinates
madrid = (40.4, -3.7)       # 40.4°N, 3.7°W
tokyo = (35.7, 139.7)       # 35.7°N, 139.7°E
sydney = (-33.9, 151.2)     # 33.9°S, 151.2°E
```

### Algorithm Selection

Choose an algorithm based on your accuracy needs:

```python
import sunwhere

# Default: PSA (fast, accurate)
result = sunwhere.sites(times, lat, lon)

# Maximum accuracy: NREL
result = sunwhere.sites(times, lat, lon, algorithm='nrel')

# Educational/testing: Iqbal
result = sunwhere.sites(times, lat, lon, algorithm='iqbal')
```

### Engine Selection

```python
# NumExpr engine (default, faster for large datasets)
result = sunwhere.sites(times, lat, lon, engine='numexpr')

# NumPy engine (may be faster for small calculations)
result = sunwhere.sites(times, lat, lon, engine='numpy')
```

## Use Case 1: Multiple Sites

Use `sites()` when you have multiple fixed locations sharing a common time grid.

### Basic Example: Three Cities

```python
import pandas as pd
import sunwhere

# Define time range
times = pd.date_range('2024-06-21', periods=24, freq='h', tz='UTC')

# Define locations
cities = ['Madrid', 'Tokyo', 'New York']
lats = [40.4, 35.7, 40.7]
lons = [-3.7, 139.7, -74.0]

# Calculate solar position
result = sunwhere.sites(times, lats, lons)

# Access results
print(result.sza)        # Solar zenith angle (time × site)
print(result.saa)        # Solar azimuth angle (time × site)
print(result.elevation)  # Solar elevation angle (time × site)
```

### Output Structure

Results have dimensions `(time, site)`:

```python
print(result.sza.shape)  # (24, 3) - 24 hours × 3 cities
print(result.sza.dims)   # ('time', 'site')

# Access specific values
print(result.sza.sel(time='2024-06-21 12:00'))  # Zenith at noon for all cities
print(result.sza.isel(site=0))  # Time series for first site (Madrid)
```

### Using Site Names

For easier data exploration, you can provide custom names for each site using the `site_names` parameter:

```python
# Define time range
times = pd.date_range('2024-06-21', periods=24, freq='h', tz='UTC')

# Define locations with names
cities = ['Madrid', 'Tokyo', 'New York']
lats = [40.4, 35.7, 40.7]
lons = [-3.7, 139.7, -74.0]

# Calculate with site names
result = sunwhere.sites(times, lats, lons, site_names=cities)

# Now you can select by name instead of index
madrid_sza = result.sza.sel(site='Madrid')
tokyo_sza = result.sza.sel(site='Tokyo')

# The site coordinate will show the names
print(result.sza.coords['site'].values)  # ['Madrid' 'Tokyo' 'New York']

# Find solar noon for each city
solar_noon_times = result.sza.idxmin(dim='time')
print(f"Solar noon in Madrid: {solar_noon_times.sel(site='Madrid')}")
```

**Note**: Without `site_names`, sites are indexed numerically (0, 1, 2, ...) and you must use `.isel(site=0)` instead of `.sel(site='Madrid')`.

### Single Location

For a single location, you can use scalars:

```python
# Solar position for Madrid
times = pd.date_range('2024-01-01', periods=8760, freq='h', tz='UTC')
result = sunwhere.sites(times, 40.4, -3.7)

print(result.sza.shape)  # (8760, 1) - one year hourly
```

### Available Variables

All returned variables are xarray DataArrays:

```python
result.sza         # Solar zenith angle (degrees)
result.saa         # Solar azimuth angle (degrees, 0=North, 90=East)
result.zenith      # Alias for sza
result.azimuth     # Alias for saa
result.elevation   # Solar elevation angle = 90 - sza (degrees)
result.cosz        # Cosine of zenith angle
result.dec         # Solar declination (degrees)
result.eot         # Equation of time (minutes)
result.ecf         # Earth-sun distance correction factor
```

### Practical Example: Finding Sunrise/Sunset

```python
import numpy as np
import matplotlib.pyplot as plt

# One day in Madrid
times = pd.date_range('2024-06-21', periods=24*60, freq='min', tz='UTC')
result = sunwhere.sites(times, 40.4, -3.7)

# Find when sun is above horizon (elevation > 0)
daylight = result.elevation > 0

# Find sunrise and sunset
sunrise_idx = np.where(daylight[:-1] != daylight[1:])[0][0]
sunset_idx = np.where(daylight[:-1] != daylight[1:])[0][1]

print(f"Sunrise: {times[sunrise_idx]}")
print(f"Sunset: {times[sunset_idx]}")

# Plot solar elevation throughout the day
result.elevation.plot()
plt.axhline(0, color='k', linestyle='--', label='Horizon')
plt.ylabel('Solar Elevation (°)')
plt.title('Solar Elevation - Madrid - Summer Solstice')
plt.legend()
plt.show()
```

## Use Case 2: Regular Grids

Use `regular_grid()` for latitude-longitude grids - ideal for climate data, maps, and gridded analyses.

### Basic Example: European Domain

```python
import numpy as np
import sunwhere

# Define grid
lats = np.arange(35, 71, 1.0)   # 35°N to 70°N, 1° resolution
lons = np.arange(-10, 41, 1.0)  # 10°W to 40°E, 1° resolution
times = pd.date_range('2024-06-21 12:00', periods=1, freq='h', tz='UTC')

# Calculate solar position across grid
result = sunwhere.regular_grid(times, lats, lons)

print(result.sza.shape)  # (1, 36, 51) - time × lat × lon
```

### Creating Solar Maps

```python
import matplotlib.pyplot as plt

# Summer solstice noon across Europe
result = sunwhere.regular_grid(times, lats, lons)

# Plot solar zenith angle
fig, ax = plt.subplots(figsize=(10, 8))
result.sza[0].plot(x='lon', y='lat', cmap='RdYlBu_r', 
                    cbar_kwargs={'label': 'Solar Zenith Angle (°)'})
plt.title('Solar Zenith Angle - Summer Solstice Noon UTC')
plt.xlabel('Longitude (°)')
plt.ylabel('Latitude (°)')
plt.show()
```

### Time Series of Gridded Data

```python
# Full year, daily resolution
times = pd.date_range('2024-01-01', '2024-12-31', freq='D', tz='UTC')
lats = np.linspace(35, 70, 36)
lons = np.linspace(-10, 40, 51)

result = sunwhere.regular_grid(times, lats, lons)

print(result.sza.shape)  # (365, 36, 51)

# Compute annual mean solar zenith angle
annual_mean = result.sza.mean(dim='time')
annual_mean.plot(x='lon', y='lat')
plt.title('Annual Mean Solar Zenith Angle')
plt.show()
```

### High-Resolution Analysis

```python
# High-resolution grid for Spain
lats = np.linspace(36, 44, 321)   # 0.025° resolution
lons = np.linspace(-10, 4, 561)
times = pd.date_range('2024-06-21 06:00', '2024-06-21 18:00', freq='h', tz='UTC')

result = sunwhere.regular_grid(times, lats, lons)

# Save to NetCDF for later analysis
result.sza.to_netcdf('solar_zenith_spain_summer.nc')
```

### Comparison with Sites

For equivalent calculations, `sites()` and `regular_grid()` produce identical results:

```python
# Using sites (flattens the grid)
lat_grid, lon_grid = np.meshgrid(lats, lons, indexing='ij')
result_sites = sunwhere.sites(times, lat_grid.ravel(), lon_grid.ravel())

# Using regular_grid (preserves 2D structure)
result_grid = sunwhere.regular_grid(times, lats, lons)

# Results are the same, just different shapes
print(result_sites.sza.shape)  # (time, lat*lon)
print(result_grid.sza.shape)   # (time, lat, lon)
```

## Use Case 3: Transects (Moving Observers)

Use `transect()` for applications where position changes with time: satellites, aircraft, ships, etc.

### Basic Example: Satellite Ground Track

```python
import numpy as np
import sunwhere

# Simplified ISS-like orbit (~90 minute period)
times = pd.date_range('2024-01-15 00:00', periods=100, freq='54s', tz='UTC')

# Simplified sinusoidal ground track
lats = 51.6 * np.sin(2 * np.pi * np.arange(100) / 100)
lons = np.linspace(-180, 180, 100, endpoint=False)

# Calculate solar position along track
result = sunwhere.transect(times, lats, lons)

print(result.sza.shape)  # (100,) - one value per time/position

# Check when satellite is in sunlight
in_sunlight = result.elevation > 0
print(f"Sunlit for {in_sunlight.sum().item()} of {len(times)} points")
```

### Flight Path Example

```python
# NYC to Tokyo flight (~10 hours)
times = pd.date_range('2024-03-20 10:00', periods=121, freq='5min', tz='UTC')

# Linear interpolation (simplified - not actual great circle)
lats = np.linspace(40.7, 35.7, 121)  # NYC to Tokyo
lons = np.linspace(-74.0, 139.7, 121)

result = sunwhere.transect(times, lats, lons, algorithm='nrel')

# Plot solar elevation throughout flight
import matplotlib.pyplot as plt
result.elevation.plot()
plt.axhline(0, color='k', linestyle='--', label='Horizon')
plt.ylabel('Solar Elevation (°)')
plt.title('Sun Position During NYC-Tokyo Flight')
plt.legend()
plt.show()
```

### Fixed Latitude or Longitude

You can keep one coordinate constant:

```python
# Moving north along Prime Meridian
times = pd.date_range('2024-06-21 12:00', periods=181, freq='h', tz='UTC')
lats = np.linspace(-90, 90, 181)  # South Pole to North Pole
lon = 0.0  # Prime Meridian (scalar)

result = sunwhere.transect(times, lats, lon)

# Zenith angle variation from pole to pole
result.sza.plot()
plt.xlabel('Time / Latitude')
plt.ylabel('Solar Zenith Angle (°)')
plt.show()
```

### Ship Route

```python
# Vessel from Mediterranean to Atlantic
times = pd.date_range('2024-07-01', periods=200, freq='h', tz='UTC')

# Route coordinates (example path)
lats = np.array([36.0, 36.5, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0, 43.0])
lons = np.array([-5.0, -6.0, -7.0, -8.0, -9.0, -9.5, -9.8, -10.0, -10.0])

# Interpolate to hourly positions
from scipy.interpolate import interp1d
voyage_lats = interp1d(np.linspace(0, len(times)-1, len(lats)), lats)(np.arange(len(times)))
voyage_lons = interp1d(np.linspace(0, len(times)-1, len(lons)), lons)(np.arange(len(times)))

result = sunwhere.transect(times, voyage_lats, voyage_lons)
```

## Advanced Topics

### Atmospheric Refraction

By default, atmospheric refraction correction is enabled. You can disable it:

```python
# Without refraction correction
result = sunwhere.sites(times, lat, lon, refraction=False)

# Compare with and without refraction
result_ref = sunwhere.sites(times, lat, lon, refraction=True)
result_noref = sunwhere.sites(times, lat, lon, refraction=False)

diff = result_ref.sza - result_noref.sza
print(f"Max refraction effect: {diff.max().item():.3f}°")
```

Refraction correction is most significant near the horizon (~0.5° at horizon).

### Working with xarray DataArrays

Sunwhere returns xarray DataArrays, which provide powerful functionality:

```python
# Select specific times
noon_values = result.sza.sel(time='2024-06-21 12:00')

# Select time ranges
summer = result.sza.sel(time=slice('2024-06-01', '2024-08-31'))

# Statistical operations
daily_max = result.sza.resample(time='D').max()
monthly_mean = result.sza.resample(time='ME').mean()

# Coordinate-based selection
madrid_data = result.sza.sel(site='Madrid')  # By name (if site_names provided)
madrid_data = result.sza.isel(site=0)  # By numeric index

# Plotting
result.sza.plot()  # Simple 1D/2D plots
result.sza.plot.contourf()  # Contour plots for 2D data

# Export to various formats
result.sza.to_netcdf('solar_zenith.nc')  # NetCDF
result.sza.to_pandas()  # Pandas DataFrame
result.sza.values  # Raw NumPy array
```

### Performance Tips

1. **Use numexpr engine** for large datasets (default)
2. **Batch calculations** when possible instead of looping
3. **Choose appropriate algorithm**: PSA for speed, NREL for accuracy
4. **Pre-allocate time arrays** with consistent frequency

```python
# Good: Single call for all locations
result = sunwhere.sites(times, all_lats, all_lons)

# Avoid: Looping over locations
for lat, lon in zip(lats, lons):
    result = sunwhere.sites(times, lat, lon)  # Slow!
```

### Handling Missing Data (NaT)

Sunwhere handles NaT (Not-a-Time) values:

```python
import numpy as np

times = pd.DatetimeIndex(['2024-01-01', pd.NaT, '2024-01-03'])
result = sunwhere.sites(times, 40.0, 0.0)

# NaT times produce NaN in output
print(result.sza)  # [value, nan, value]
```

## Next Steps

- Explore the [API Reference](api.md) for complete documentation
- Check the GitHub repository for more examples
- Read the scientific references for algorithm details

## Troubleshooting

### ValueError: Latitude out of range

Ensure latitudes are within \[-90, 90\]:

```python
# Invalid
sunwhere.sites(times, 91.0, 0.0)  # Error!

# Valid
sunwhere.sites(times, 90.0, 0.0)  # OK
```

### ValueError: Longitude out of range

Ensure longitudes are within [-180, 180):

```python
# Invalid
sunwhere.sites(times, 0.0, 180.0)  # Error! (exactly 180 is invalid)

# Valid
sunwhere.sites(times, 0.0, 179.9)  # OK
sunwhere.sites(times, 0.0, -180.0)  # OK
```

### Shape mismatch errors

For `sites()`, lat and lon must have the same length:

```python
# Invalid
sunwhere.sites(times, [40, 35], [0, 5, 10])  # Error! Different lengths

# Valid
sunwhere.sites(times, [40, 35, 30], [0, 5, 10])  # OK
```

For `transect()`, time, lat, and lon must have the same length (unless lat/lon are scalars).
