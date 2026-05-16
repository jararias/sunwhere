"""
Example demonstrating the site_names feature in sunwhere.sites()

This shows how to use custom site names for more intuitive data access.
"""

import pandas as pd
import sunwhere

# Example 1: Multiple cities with names
print("=" * 60)
print("Example 1: Multiple cities with named sites")
print("=" * 60)

times = pd.date_range('2024-06-21', periods=24, freq='h', tz='UTC')
cities = ['Madrid', 'Tokyo', 'New York', 'Sydney']
lats = [40.4168, 35.6762, 40.7128, -33.8688]
lons = [-3.7038, 139.6503, -74.0060, 151.2093]

result = sunwhere.sites(times, lats, lons, site_names=cities)

print(f"\nSite names: {list(result.sza.coords['site'].values)}")
print(f"Data shape: {result.sza.shape} (time × sites)")

# Access data by site name
madrid_sza = result.sza.sel(site='Madrid')
tokyo_sza = result.sza.sel(site='Tokyo')

print(f"\nMadrid noon SZA: {madrid_sza.values[12]:.2f}°")
print(f"Tokyo noon SZA: {tokyo_sza.values[12]:.2f}°")

# Find solar noon for each city
solar_noon_times = result.sza.idxmin(dim='time')
print("\nSolar noon times:")
for city in cities:
    print(f"  {city}: {solar_noon_times.sel(site=city).values}")


# Example 2: Energy production sites
print("\n" + "=" * 60)
print("Example 2: Solar farm monitoring sites")
print("=" * 60)

sites = ['Plant-A', 'Plant-B', 'Plant-C']
lats = [37.5, 38.2, 39.1]
lons = [-6.5, -5.8, -5.1]

times = pd.date_range('2024-03-15 06:00', periods=13, freq='h', tz='UTC')
result = sunwhere.sites(times, lats, lons, site_names=sites)

# Calculate daylight hours (elevation > 0)
daylight_hours = (result.elevation > 0).sum(dim='time')

print("\nDaylight hours on 2024-03-15:")
for site in sites:
    hours = daylight_hours.sel(site=site).values
    print(f"  {site}: {hours} hours")


# Example 3: Comparison with and without names
print("\n" + "=" * 60)
print("Example 3: With vs without site_names")
print("=" * 60)

times = pd.date_range('2024-01-01 12:00', periods=1, tz='UTC')
lats = [40.0, 41.0]
lons = [-3.0, -4.0]

# Without names (numeric indices)
result_numeric = sunwhere.sites(times, lats, lons)
print(f"\nWithout names - site coordinate: {list(result_numeric.sza.coords['site'].values)}")
print("Access by index: result.sza.isel(site=0)")

# With names
result_named = sunwhere.sites(times, lats, lons, site_names=['North', 'South'])
print(f"\nWith names - site coordinate: {list(result_named.sza.coords['site'].values)}")
print("Access by name: result.sza.sel(site='North')")
print(f"Value: {result_named.sza.sel(site='North').values[0]:.2f}°")

print("\n" + "=" * 60)
print("✓ All examples completed successfully!")
print("=" * 60)
