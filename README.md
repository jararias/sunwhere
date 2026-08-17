![python versions](https://img.shields.io/badge/python-3.11%2C3.12%2C3.13-blue.svg) ![tests-badge](https://raw.githubusercontent.com/jararias/sunwhere/main/docs/images/tests-badge.svg) ![coverage-badge](https://raw.githubusercontent.com/jararias/sunwhere/main/docs/images/coverage-badge.svg) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21975085.svg)](https://doi.org/10.5281/zenodo.21975085)


# sunwhere. Solar position calculation for solar resource assessment

sunwhere is a Python library designed for fast and accurate calculations of solar position for solar resource applications 🌞.

sunwhere provides solar zenith and azimuth angles, sun-earth distance correction factor, and secondary parameters such as solar declination, equation of time, and many more. It's optimized for typical workflows in solar energy research and engineering.

It is optimized for three common use cases that cover most practical scenarios:

- sites() - Multiple arbitrary locations with common time grid

- regular_grid() - Lat-lon regular grids (ideal for spatio-temporal data)

- transect() - Moving observers (satellites, aircraft, ships)

## Installation

```sh
pip install sunwhere
```

or using [uv](https://docs.astral.sh/uv/):

```sh
uv add sunwhere
```

## Quick examples

```python
import sunwhere
lat, lon = 36.9369, -3.7943
times = pd.date_range("2026-06-02", periods=24*60, freq="min", tz="Europe/Madrid")
solpos = sunwhere.sites(times, lat, lon, site_names=("Jayena",))
solar_zenith = solpos.zenith  # xarray's dataarray with dims (time, site)
solar_zenith = solpos.zenith.sel(site="Jayena").to_pandas()  # pandas series
```

```python
import sunwhere
lats = np.arange(30, 60.1, 0.5)
lons = np.arange(-10, 50.1, 0.5)
times = pd.date_range("2026-06-02", periods=24, freq="h", tz="Europe/Madrid")
solpos = sunwhere.regular_grid(times, lat, lon)
solar_zenith = solpos.zenith  # xarray's dataarray with dims (time, lat, lon)
```


## Documentation

Full documentation — installation guide, user guide, quick reference, and API reference — is available at:

**<https://jararias.github.io/sunwhere/>**

## Citation

If you use sunwhere in your research, please cite:

```bibtex
@software{sunwhere2024,
  author = {Ruiz-Arias, Jose A.},
  title = {sunwhere: Solar position for solar resource assessment},
  year = {2024},
  url = {https://github.com/jararias/sunwhere}
}
```

## License

[CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/) — free for non-commercial use with attribution.
