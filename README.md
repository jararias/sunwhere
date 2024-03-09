![python versions](https://img.shields.io/badge/python-3.8%2C3.9%2C3.10-blue.svg)

TODO:
- [x] review docstrings in usecases and core
- [x] add header image (sun track)
- [x] move benchmark to main.py
- [ ] add documentation to README.md (images with folding code in usage showcases)
- [ ] add solar chart with optional analemas in main.py

# Solar position for solar resource assessment

![header](assets/headerfig.png)

*sunwhere* is tailored for typical applications in solar resource assessment. It provides the solar zenith and azimuth angles, the sun-earth's distance correction factor, and secondary parameters such as solar declination, equation of time, and sunrise and sunset times, among others.

## Main features
*sunwhere* optionally uses the NREL[^1], Plataforma Solar de Almería (PSA)[^2], SolTrack[^3] or Iqbal[^4] solar position algorithms, which provide alternative levels of speed and accuracy for each application's requirements.

*sunwhere* focuses on usage cases to optimize the computing performance. Three cases are specifically considered that hopefully encompass most practical situations:

- [_sites_](./api-docs/usecases.md#function-sites), intended to evaluate the solar position across any number of arbitrary locations throughout a common time grid at once.

- [_regular_grid_](./api-docs/usecases.md#function-regular_grid), intended for lon-lat regular grids throughout a common time grid.

- [_transect_](./api-docs/usecases.md#function-transect), intended for changing locations over time (e.g., scanning path of a satellite sensor).

Conversely, other packages only consider single-location calculations, having to iterate over each location in multi-site evaluations.

*sunwhere* returns a [Sunpos instance](api-docs/base.md) that provides access to the calculation results via xarray [DataArrays](https://docs.xarray.dev/en/stable/generated/xarray.DataArray.html) that adapt to each specific usage case and provide access to the rich xarray's [API](https://docs.xarray.dev/en/latest/api.html#).

[^1]: Reda I and Andreas A, 2003. Solar Position Algorithm for Solar Radiation Applications. 55 pp.; NREL Report No. TP-560-34302, Revised January 2008 [pdf](http://www.nrel.gov/docs/fy08osti/34302.pdf) [url](https://midcdmz.nrel.gov/spa/).
[^2]: Blanco, M. et al. 2020. Updating the PSA sun position algorithm. Solar Energy, Vol. 212, pp. 339-341, and Blanco-Muriel, M. et al. 2001. Computing the solar vector. Solar Energy, Vol. 70, pp. 431-441 [url](https://doi.org/10.1016/j.solener.2020.10.084).
[^3]: van der Sluys M and van Kan P, 2022. SolTrack: a free, fast and accurate routine to compute the position of the Sun [url](https://doi.org/10.48550/arXiv.2209.01557) [code](https://github.com/MarcvdSluys/SolTrack-Python)
[^4]: Iqbal, M., An introduction to solar radiation. Academic Press. 1983 [url](https://www.sciencedirect.com/book/9780123737502/an-introduction-to-solar-radiation)

## Installation notes

```sh
python3 -m pip install git+https://github.com/jararias/sunwhere@main
```

<!-- If you want to run the *sunwhere*'s benchmark script, first, be sure that all dependences are installed:

```sh
python3 -m pip install "sunwhere[benchmark] @ git+https://github.com/jararias/sunwhere@main"
```

> [!NOTE]
> See examples in https://pip.pypa.io/en/stable/cli/pip_install/#examples. See also [PEP 440](https://peps.python.org/pep-0440/#direct-references) for usage of @.

Or, from a cloned local copy:

```sh
python3 -m pip install <path_to_local_copy>/.[benchmark]
``` -->

## Usage showcases

### Case 1. sunwhere.sites

It requires a 1-dim sequence of times (it they are not timezone aware, UTC is assumed) and 1-dim arrays of latitude and longitude. They can also be scalars for single-location calculations. The latitude and longitude arrays must have exactly the same length. They are the geographic coordinates of the locations where solar position will be evaluated. The following image shows the solar zenith and elevation angles produced at 5 sites randomly selected.

![case1: sites](assets/case1_sites.png)

<details>

<summary>Python's code</summary>

```python
import numpy as np
import pylab as pl
import pandas as pd
import cartopy.crs as ccrs

import sunwhere

# select locations randomly...
n_sites = 5
rg = np.random.RandomState(20240307)
lats = -90. + 180*rg.random(n_sites)
lons = -180 + 360*rg.random(n_sites)

# ...and the time grid
times = pd.date_range('2023-01-15', '2023-01-18', freq='1min', inclusive='left', tz='UTC')
# or, using numpy datetime64:
# times = np.arange('2023-01-15', '2023-01-18', np.timedelta64(1, 'm'), dtype='datetime64[ns]')
# However, numpy does not provide good support for time zones if they are needed.

# run sunwhere.sites... len(lats) must be equal to len(lons)!!
sw = sunwhere.sites(times, lats, lons)  # algorithm='psa' refraction=True

print(sw.sza.head())

# draw some results...
pl.rcParams['axes.titlesize'] = 'small'
pl.rcParams['axes.labelsize'] = 'small'
pl.rcParams['xtick.labelsize'] = 'small'
pl.rcParams['ytick.labelsize'] = 'small'
pl.figure(figsize=(12, 6), layout='constrained')
for k in range(len(sw.sza.location)):
    pl.subplot(2, 3, k+2)
    sw.sza.isel(location=k).plot(label='zenith')
    sw.elevation.isel(location=k).plot(label='elevation')
pl.legend()

# draw the locations in a map...
ax = pl.subplot(231, projection=ccrs.PlateCarree())
location = sw.sza.coords['location'].to_numpy()
latitude = sw.sza.coords['latitude'].to_numpy()
longitude = sw.sza.coords['longitude'].to_numpy()
for loc, lon, lat in zip(location, longitude, latitude):
    ax.plot(lon, lat, 'r.', ms=8)
    ax.text(lon, lat, loc, ha='left', va='bottom')
ax.coastlines(lw=0.5, color='0.5')
ax.set_global()

pl.show()
```
</details>

### Case 2. sunwhere.regular_grid

As `sunwhere.sites`, `sunwhere.regular_grid` requires a 1-dim sequence of times (it they are not timezone aware, UTC is assumed) and 1-dim arrays of latitude and longitude. However, the length of the latitude and longitude arrays does not have to be necessarily the same. Now, they represent the rectangular coordinates of the lon-lat regular grid. The following image shows the solar zenith angle (left column), the cosine of solar zenith angle, which is important for the evaluation of solar irradiance (middle column) and the cosine of the incidence angle for a plane of array with an inclination of 30&#x00b0; and an azimuth of 60&#x00b0;. The calculations are performed over a regular grid with 1-deg cellsize.

![case2: sites](assets/case2_regular_grid.png)

<details>

<summary>Python's code</summary>

```python
import numpy as np
import pylab as pl
import cartopy.crs as ccrs

import sunwhere

# select the spatial grid...
lats = np.arange(-89.5, 90, 1)
lons = np.arange(-179.5, 180, 1)

# ...and the time grid
times = np.arange('2023-01-15', '2023-01-16', np.timedelta64(6, 'h'), dtype='datetime64[ns]')
# or: times = pd.date_range('2023-01-15', '2023-01-18', freq='1min', inclusive='left')

# run sunwhere.sites...
sw = sunwhere.regular_grid(times, lats, lons)  # algorithm='psa' refraction=True

print(sw.sza.head())

# draw some results...
pl.rcParams['axes.titlesize'] = 'small'
pl.rcParams['axes.labelsize'] = 'small'
pl.rcParams['xtick.labelsize'] = 'small'
pl.rcParams['ytick.labelsize'] = 'small'
pl.figure(figsize=(14, 8), layout='constrained')
for k in range(len(sw.sza.time)):
    ax = pl.subplot(4, 3, 3*k+1, projection=ccrs.PlateCarree())
    sw.sza.isel(time=k).plot(ax=ax, cmap='inferno_r')
    ax.coastlines(lw=0.3, color='w')
    ax = pl.subplot(4, 3, 3*k+2, projection=ccrs.PlateCarree())
    sw.cosz.isel(time=k).plot(ax=ax, cmap='RdBu_r')
    ax.coastlines(lw=0.5, color='0.3')
    ax = pl.subplot(4, 3, 3*k+3, projection=ccrs.PlateCarree())
    sw.incidence(30, 60).isel(time=k).plot(ax=ax, cmap='RdBu_r')
    ax.coastlines(lw=0.5, color='0.3')

pl.show()
```
</details>

### Case 3. sunwhere.transect

```python
```

### Case 4. Command line interface

```sh
sunwhere at --help  # to calculate solar position
sunwhere chart --help  # to produce solar charts
```

## Benchmark: which SPA to choose?

que criterio? evaluamos las efemerides? o quiza es mejor la radiacion? como hacerlo? con que comparamos las efemerides? que pegas tiene usar un modelo de cielo despejado?
importa la velocidad del algoritmo? que paquetes python hay disponibles contra los que compararar?

### Speed benchmark
:construction_worker:

### Accuracy benchmark

![accuracy benchmark](assets/accuracy_benchmark.png)

### References

