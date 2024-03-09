![python versions](https://img.shields.io/badge/python-3.8%2C3.9%2C3.10-blue.svg)

TODO:
- [x] review docstrings in usecases and core
- [x] add header image (sun track)
- [x] move benchmark to main.py
- [ ] add documentation to README.md (images with folding code in usage showcases)
- [ ] add solar chart with optional analemas in main.py

# Solar position for solar resource assessment

![header](assets/headerfig.png)

*sunwhere* is tailored for typical applications in solar resource assessment. It provides the solar zenith and azimuth angles, the sun-earth's distance correction factor, and secondary parameters such as solar declination, equation of time, and sunrise and sunset hours, among others.

## Main features
*sunwhere* optionally uses the NREL[^1], Plataforma Solar de Almería (PSA)[^2], SolTrack[^3] or Iqbal[^4] solar position algorithms, which provide three different levels of speed and accuracy to adapt the operation to the requirements of different applications.

Unlike other packages, *sunwhere* focuses on usage cases. Three of them are specifically considered:

- [_sites_](./api-docs/usecases.md#function-sites), intended for calculations across a number of arbitrary locations throughout a common time grid.

- [_regular_grid_](./api-docs/usecases.md#function-regular_grid), intended for lon-lat regular grids of locations throughout a common time grid

- [_transect_](./api-docs/usecases.md#function-transect), intended for changing locations over time (e.g., scanning path of a satellite sensor).

Contrarily, other packages consider only the single-location case, having to iterate the calculations over each location in multi-site evaluations.

*sunwhere* wraps its outputs in usage-case-tailored [xarray](https://docs.xarray.dev/) [DataArrays](https://docs.xarray.dev/en/latest/user-guide/data-structures.html) to provide a clear distinction between usage cases and easy access to the rich xarray's [API](https://docs.xarray.dev/en/latest/api.html#).

[^1]: Reda I and Andreas A, 2003. Solar Position Algorithm for Solar Radiation Applications. 55 pp.; NREL Report No. TP-560-34302, Revised January 2008 [pdf](http://www.nrel.gov/docs/fy08osti/34302.pdf) [url](https://midcdmz.nrel.gov/spa/).
[^2]: Blanco, M. et al. 2020. Updating the PSA sun position algorithm. Solar Energy, Vol. 212, pp. 339-341, and Blanco-Muriel, M. et al. 2001. Computing the solar vector. Solar Energy, Vol. 70, pp. 431-441 [url](https://doi.org/10.1016/j.solener.2020.10.084).
[^3]: van der Sluys M and van Kan P, 2022. SolTrack: a free, fast and accurate routine to compute the position of the Sun [url](https://doi.org/10.48550/arXiv.2209.01557) [code](https://github.com/MarcvdSluys/SolTrack-Python)
[^4]: Iqbal, M., An introduction to solar radiation. Academic Press. 1983 [url](https://www.sciencedirect.com/book/9780123737502/an-introduction-to-solar-radiation)

## Installation notes

```sh
python3 -m pip install git+https://github.com/jararias/sunwhere@main
```

If you want to run the *sunwhere*'s benchmark script, first, be sure that all dependences are installed:

```sh
python3 -m pip install "sunwhere[benchmark] @ git+https://github.com/jararias/sunwhere@main"
```

> [!NOTE]
> See examples in https://pip.pypa.io/en/stable/cli/pip_install/#examples. See also [PEP 440](https://peps.python.org/pep-0440/#direct-references) for usage of @.

Or, from a cloned local copy:

```sh
python3 -m pip install <path_to_local_copy>/.[benchmark]
```



## Usage showcases

### sunwhere.sites

```python
import numpy as np
import pylab as pl

import sunwhere

# select n_sites random locations...
n_sites = 6
rg = np.random.RandomState(2024)
lats = -90. + 180*rg.random(n_sites)
lons = -180 + 360*rg.random(n_sites)

# ...and the time grid
times = np.arange('2023-01-15', '2023-01-18', np.timedelta64(1, 'm'), dtype='datetime64[ns]')
# or: times = pd.date_range('2023-01-15', '2023-01-18', freq='1min', inclusive='left')

# run sunwhere.sites... len(lats) must be equal to len(lons)!!
sw = sunwhere.sites(times, lats, lons)  # algorithm='psa' refraction=True

print(sw.sza.head())

# plot the results...
pl.figure(layout='constrained')
for k in range(len(sw.sza.location)):
    pl.subplot(2, 3, k+1)
    sw.sza.isel(location=k).plot()
    sw.elevation.isel(location=k).plot()
pl.show()
```

### sunwhere.regular_grid

```python
import numpy as np
import pylab as pl
import cartopy.crs as ccrs

import sunwhere

# select the spatial grid...
lats = np.arange(-89.5, 90, 1)
lons = np.arange(-179.5, 180, 1)

# ...and the time grid
times = np.arange('2023-01-15', '2023-01-16', np.timedelta64(4, 'h'), dtype='datetime64[ns]')
# or: times = pd.date_range('2023-01-15', '2023-01-18', freq='1min', inclusive='left')

# run sunwhere.sites...
sw = sunwhere.regular_grid(times, lats, lons)  # algorithm='psa' refraction=True
print(sw.sza.head())

# plot the results...
w, h = pl.figaspect(0.70)
pl.figure(figsize=(w, h), layout='constrained')
for k in range(len(sw.sza.time)):
    ax = pl.subplot(3, 2, k+1, projection=ccrs.PlateCarree())
    sw.sza.isel(time=k).plot(ax=ax, cmap='jet')
    ax.coastlines()
pl.show()
```

### sunwhere.transect

```python
```

### Speed benchmark
:construction_worker:

### Accuracy benchmark

![accuracy benchmark](assets/accuracy_benchmark.png)

### References

