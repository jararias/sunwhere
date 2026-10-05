# Changelog

## v1.5.0 (2026-10-05)

### Changed

- The `time` coordinate of the output DataArrays is now always **naive UTC** (`datetime64[ns]`), also for timezone-aware input times. Before, it kept the input timezone, which broke xarray operations such as `groupby('time.month')`, `interp`, `to_netcdf` or alignment with naive data (`TypeError: Cannot interpret 'datetime64[us, UTC]' as a data type`). The input times, with their timezone, are still available in `Sunpos.times`. To select with a tz-aware timestamp, convert it first: `ts.tz_convert(None)`.
- `sunwhere.SPA_ALGORITHMS` lists only the public algorithms (`psa`, `iqbal`, `nrel`).

### Added

- `Sunpos.azimuth_north`: solar azimuth clockwise from north in [0°, 360°), as in pvlib.
- `transect()` accepts a scalar latitude or longitude, used for all times (as documented).
- Documentation: Benchmarks page, time zone handling, azimuth convention, and a clearer license statement.

### Fixed

- Sub-second input times were truncated to whole seconds.
- `Sunpos.sunset('deg')` and `Sunpos.sunset('rad')` raised `UnboundLocalError`.
- `sunwhere.universal_time_coordinated()` raised `AttributeError`.
- `Sunpos.incidence()` returned dims transposed with respect to `zenith`.
- Iqbal and SolTrack algorithms applied the refraction correction below the horizon (the SolTrack formula has a pole near -5° elevation).
- Latitudes and longitudes were cast to float32 inside the algorithms; they are now float64.
- CLI: any explicit time raised `TypeError`, and "Time (UTC)" always showed the current time. A clear message is shown if the `cli` extra is not installed.
- Documentation: wrong azimuth convention in the Quick Reference and User Guide, `NameError` in the README `regular_grid` example and in the "Fixed Latitude" example, speed claim on the home page (about 20x at 10 sites, not 100x), and the DOI badge now uses the concept DOI (always the latest version).
