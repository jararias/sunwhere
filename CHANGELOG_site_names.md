# Changelog - site_names Feature

## New Feature: Named Sites for sites() Function

**Date**: 2024  
**Version**: Development

### Summary

Added optional `site_names` parameter to `sunwhere.sites()` function, enabling users to assign custom names to observation sites for more intuitive data access.

### Changes

#### API Changes

**New Parameter**: `site_names` (optional)
- **Type**: array-like of str
- **Default**: None
- **Description**: Custom names for each site. Must match the length of latitude/longitude arrays.

**Example**:
```python
import sunwhere
import pandas as pd

times = pd.date_range('2024-06-21', periods=24, freq='h', tz='UTC')
cities = ['Madrid', 'Tokyo', 'New York']
lats = [40.4168, 35.6762, 40.7128]
lons = [-3.7038, 139.6503, -74.0060]

# With site names
result = sunwhere.sites(times, lats, lons, site_names=cities)

# Access by name
madrid_sza = result.sza.sel(site='Madrid')
```

#### Benefits

1. **Improved readability**: Use descriptive names instead of numeric indices
2. **Intuitive selection**: `result.sza.sel(site='Madrid')` vs `result.sza.isel(site=0)`
3. **Better documentation**: Self-documenting code with meaningful site labels
4. **Backward compatible**: Existing code without `site_names` continues to work unchanged

#### Files Modified

- `src/sunwhere/usecases.py`: Added parameter and validation
- `src/sunwhere/_base.py`: Updated Sunpos constructor and xarray coordinate creation
- `docs/user-guide.md`: Added documentation section with examples
- `tests/test_site_names.py`: Comprehensive test suite (11 tests)
- `examples/site_names_demo.py`: Demonstration examples

#### Validation

The implementation includes robust validation:
- Checks that `site_names` length matches latitude/longitude
- Verifies 1-dimensional input
- Converts various input types (list, numpy array) to consistent format
- Supports Unicode characters and spaces in names

#### Backward Compatibility

✅ **Fully backward compatible**
- Without `site_names` parameter, sites are indexed numerically (0, 1, 2, ...)
- All existing code continues to work without modification
- `site_names=None` is the default behavior

#### Testing

- 11 new unit tests added
- All tests pass (100%)
- Coverage includes:
  - Basic functionality with names
  - Backward compatibility without names
  - Data selection by name
  - Input validation
  - Edge cases (Unicode, spaces, single site)
  - All algorithms (psa, nrel, iqbal)

### Migration Guide

No migration needed. This is an optional enhancement that doesn't break existing code.

**Before** (still works):
```python
result = sunwhere.sites(times, lats, lons)
first_site = result.sza.isel(site=0)  # Numeric index
```

**After** (new capability):
```python
result = sunwhere.sites(times, lats, lons, site_names=['Madrid', 'Tokyo'])
madrid = result.sza.sel(site='Madrid')  # Named selection
```

### Notes

- Site names are stored as xarray coordinate for the 'site' dimension
- Names must be unique for proper xarray selection
- Names are converted to strings internally
