## Algorithm Details

### PSA (Plataforma Solar de Almería)

**Accuracy**: Average error ~0.002° (zenith and azimuth)  
**Speed**: Fast  
**Valid range**: Optimized for 2020-2050  
**References**: 
- Blanco-Muriel, M. et al., 2001. Computing the solar vector. Solar Energy, 70(5), 431-441.
- Blanco, M. et al., 2020. Updating the PSA sun position algorithm. Solar Energy, 212, 339-341.

**When to use**: Default choice for most applications. Excellent balance of speed and accuracy.

### NREL SPA

**Accuracy**: ±0.0003° (highest accuracy)  
**Speed**: Slower than PSA  
**Valid range**: Years -2000 to 6000  
**Reference**: Reda, I. and Andreas, A., 2003. Solar Position Algorithm for Solar Radiation Applications. NREL Report No. TP-560-34302.

**When to use**: When maximum accuracy is required, especially for long-term studies or high-precision applications.

### Iqbal

**Accuracy**: Lower than PSA and NREL (errors can exceed 1°)  
**Speed**: Similar to PSA  
**Valid range**: General use  
**Reference**: Iqbal, M., 1983. An introduction to solar radiation. Academic Press.

**When to use**: Educational purposes or applications where approximate values are sufficient.

## Engine Selection

### NumExpr Engine (default)

- **Best for**: Large datasets (>1000 points)
- **Performance**: Uses multi-threading and optimized evaluation
- **Memory**: Efficient memory usage
- **Availability**: Works with all algorithms

### NumPy Engine

- **Best for**: Small datasets or when numexpr is unavailable
- **Performance**: Standard NumPy operations
- **Memory**: Standard NumPy memory usage
- **Availability**: Works with all algorithms except NREL

## Coordinate Conventions

### Latitude

- Range: -90° to +90°
- Positive: North
- Negative: South
- Examples: 
  - 40.0° = 40°N (Madrid)
  - -33.9° = 33.9°S (Sydney)

### Longitude

- Range: -180° to +180° (exclusive upper bound)
- Positive: East
- Negative: West
- Zero: Prime Meridian (Greenwich)
- Examples:
  - -3.7° = 3.7°W (Madrid)
  - 139.7° = 139.7°E (Tokyo)

### Solar Azimuth Angle

sunwhere returns the solar azimuth angle (`azimuth`, alias `saa`) in the range **\[-180°, 180°\], with zero south**:

- 0° = South
- -90° = East (negative values before solar noon)
- +90° = West (positive values after solar noon)
- ±180° = North

pvlib and most PV software measure azimuth **clockwise from north in [0°, 360°)** (0° = North, 90° = East, 180° = South, 270° = West). To convert, add 180°, or use the `azimuth_north` property:

```python
result.azimuth_north  # == (result.azimuth + 180) % 360, zero north, clockwise
```

### Time Zones

- Timezone-naive input times are assumed to be **UTC**.
- Timezone-aware input times are converted to UTC.
- The `time` coordinate of all output DataArrays is **naive UTC** (`datetime64[ns]`), whatever the input timezone. xarray does not fully support timezone-aware coordinates.
- The original input times, with their timezone, are kept in `result.times` and `result.timezone`.

```python
times = pd.date_range('2024-06-21', periods=24, freq='h', tz='Europe/Madrid')
result = sunwhere.sites(times, 40.4, -3.7)
result.sza.sel(time='2024-06-21 10:00')                  # 10:00 UTC (12:00 in Madrid)
result.sza.sel(time=times[12].tz_convert(None))           # select with a tz-aware timestamp
result.sza.to_pandas().tz_localize('UTC').tz_convert('Europe/Madrid')  # back to local time
```

### Sunrise and Sunset Hour Angles

`sunrise('deg')` returns a positive hour angle and `sunset('deg')` returns a negative one, because the hour angle is positive in the morning. Both have the same magnitude.

### Solar Zenith Angle

- 0° = Sun directly overhead (zenith)
- 90° = Sun at horizon
- \>90° = Sun below horizon

Related to elevation: zenith = 90° - elevation

## Error Handling

sunwhere validates inputs and raises appropriate errors:

### ValueError: Invalid latitude

```python
sunwhere.sites(times, 95.0, 0.0)  # Latitude must be in [-90, 90]
```

### ValueError: Invalid longitude

```python
sunwhere.sites(times, 0.0, 185.0)  # Longitude must be in [-180, 180)
```

### ValueError: Invalid algorithm

```python
sunwhere.sites(times, 0.0, 0.0, algorithm='invalid')
# Valid algorithms are: psa, iqbal, nrel
```

### ValueError: Shape mismatch

```python
sunwhere.sites(times, [40, 35], [0, 5, 10])
# Latitude and longitude must have same length
```

### ValueError: Invalid engine

```python
sunwhere.sites(times, 0.0, 0.0, algorithm='nrel', engine='numpy')
# NREL algorithm only supports numexpr engine
```
