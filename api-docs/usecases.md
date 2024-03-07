<!-- markdownlint-disable -->

<a href="https://github.com/jararias/sunwhere/blob/master/sunwhere/usecases.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square" /></a>

# <kbd>module</kbd> `sunwhere.usecases`





---

<a href="https://github.com/jararias/sunwhere/blob/master/sunwhere/usecases.py#L22"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square" /></a>

## <kbd>function</kbd> `sites`

```python
sites(
    times,
    latitude,
    longitude,
    algorithm='psa',
    refraction=True,
    engine='numexpr'
)
```

Solar position across an arbitrary number of locations throughout a common time grid. 

Latitude and longitude *must* have the same size. 



**Args:**
 
 - <b>`times`</b> (1-dim array of datetime64, or similar, of size K):  Time  instants at which solar position is evaluated. Similar types  are, for instance, datetime and pandas DatetimeIndex. 
 - <b>`latitude`</b> (scalar or 1-dim array like of floats of size J):  Latitudes  (degrees) where solar position is evaluated. It must be in the  range [-90, 90]. 
 - <b>`longitude`</b> (scalar or 1-dim array like of floats of size J):  Longitudes  (degrees) where solar position is evaluated. It must be in the  range [-180, 180]. 
 - <b>`algorithm`</b> ({_psa_, _nrel_, _iqbal_}):  Solar position algorithm.  _nrel_ is for the NREL's Solar Position Algorithm [1], valid  for the years -2000 to 6000. Its expected uncertainty is  +/- 0.0003 degrees. _psa_ is for the Plataforma Solar de  Almería's (PSA) algorithm [2] with updated coefficients for  the period 2020-2050 [3]. Its expected average error is 0.002  degrees, poorer than that of _nrel_ but it is much faster.  _iqbal_ is for the algorithm described in Iqbal, M. [4], which  has lower accuracy than the former algorithms, but it is  slightly faster than the _psa_. 
 - <b>`refraction`</b> (bool):  Whether atmospheric refraction must be considered. 
 - <b>`engine`</b> ({_numpy_, _numexpr_}):  Baseline code implementation to perform  the solar position calculations. _numexpr_ is expected to be faster  than numpy, especially for big spatial and/or temporal grids. 



**Returns:**
 A Sunpos instance in which the shape of the variables that are not location-dependent (e.g., declination) is (K,) and that of the location-dependent variables (e.g., zenith) is (K, J). 



**References:**
 

[1] Reda I and Andreas A, 2003. Solar Position Algorithm for Solar Radiation Applications. 55 pp.; NREL Report No. TP-560-34302, Revised January 2008 [pdf](http://www.nrel.gov/docs/fy08osti/34302.pdf) [url](https://midcdmz.nrel.gov/spa/). 

[2] Blanco-Muriel, M. et al. 2001. Computing the solar vector. Solar Energy, Vol. 70, pp. 431-441 doi: [10.1016/S0038-092X(00)00156-0](https://doi.org/10.1016/S0038-092X(00)00156-0). 

[3] Blanco, M. et al. 2020. Updating the PSA sun position algorithm. Solar Energy, Vol. 212, pp. 339-341 doi: [10.1016/j.solener.2020.10.084](https://doi.org/10.1016/j.solener.2020.10.084). 

[4] Iqbal, M., An introduction to solar radiation. Academic Press. 1983 [url](https://www.sciencedirect.com/book/9780123737502/an-introduction-to-solar-radiation) 


---

<a href="https://github.com/jararias/sunwhere/blob/master/sunwhere/usecases.py#L112"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square" /></a>

## <kbd>function</kbd> `regular_grid`

```python
regular_grid(
    times,
    latitude,
    longitude,
    algorithm='psa',
    refraction=True,
    engine='numexpr'
)
```

Solar position across a lon-lat regular grid throughout a common time grid. 



**Args:**
 
 - <b>`times`</b> (1-dim array of datetime64, or similar, of size K):  Time  instants at which solar position is evaluated. Similar types  are, for instance, datetime and pandas DatetimeIndex. 
 - <b>`latitude`</b> (scalar or 1-dim array like of floats of size J):  Latitudes  (degrees) where solar position is evaluated. It must be in the  range [-90, 90]. 
 - <b>`longitude`</b> (scalar or 1-dim array like of floats of size I):  Longitudes  (degrees) where solar position is evaluated. It must be in the  range [-180, 180]. 
 - <b>`algorithm`</b> ({_psa_, _nrel_, _iqbal_}):  Solar position algorithm.  _nrel_ is for the NREL's Solar Position Algorithm [1], valid  for the years -2000 to 6000. Its expected uncertainty is  +/- 0.0003 degrees. _psa_ is for the Plataforma Solar de  Almería's (PSA) algorithm [2] with updated coefficients for  the period 2020-2050 [3]. Its expected average error is 0.002  degrees, poorer than that of _nrel_ but it is much faster.  _iqbal_ is for the algorithm described in Iqbal, M. [4], which  has lower accuracy than the former algorithms, but it is  slightly faster than the _psa_. 
 - <b>`refraction`</b> (bool):  Whether atmospheric refraction must be considered. 
 - <b>`engine`</b> ({_numpy_, _numexpr_}):  Baseline code implementation to perform  the solar position calculations. _numexpr_ is expected to be faster  than numpy, especially for big spatial and/or temporal grids. 



**Returns:**
 A Sunpos instance in which the shape of the variables that are not location-dependent (e.g., declination) is (K,) and that of the location-dependent variables (e.g., zenith) is (K, J, I). 



**References:**
 

[1] Reda I and Andreas A, 2003. Solar Position Algorithm for Solar Radiation Applications. 55 pp.; NREL Report No. TP-560-34302, Revised January 2008 [pdf](http://www.nrel.gov/docs/fy08osti/34302.pdf) [url](https://midcdmz.nrel.gov/spa/). 

[2] Blanco-Muriel, M. et al. 2001. Computing the solar vector. Solar Energy, Vol. 70, pp. 431-441 doi: [10.1016/S0038-092X(00)00156-0](https://doi.org/10.1016/S0038-092X(00)00156-0). 

[3] Blanco, M. et al. 2020. Updating the PSA sun position algorithm. Solar Energy, Vol. 212, pp. 339-341 doi: [10.1016/j.solener.2020.10.084](https://doi.org/10.1016/j.solener.2020.10.084). 

[4] Iqbal, M., An introduction to solar radiation. Academic Press. 1983 [url](https://www.sciencedirect.com/book/9780123737502/an-introduction-to-solar-radiation) 


---

<a href="https://github.com/jararias/sunwhere/blob/master/sunwhere/usecases.py#L194"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square" /></a>

## <kbd>function</kbd> `transect`

```python
transect(
    times,
    latitude,
    longitude,
    algorithm='psa',
    refraction=True,
    engine='numexpr'
)
```

Solar position across a transect. 

In a transect, the observer's position changes throughout time. 



**Args:**
 
 - <b>`times`</b> (1-dim array of datetime64, or similar, of size K):  Time  instants at which solar position is evaluated. Similar types  are, for instance, datetime and pandas DatetimeIndex. 
 - <b>`latitude`</b> (scalar or 1-dim array like of floats of size K):  Latitudes  (degrees) where solar position is evaluated. It must be in the  range [-90, 90]. 
 - <b>`longitude`</b> (scalar or 1-dim array like of floats of size K):  Longitudes  (degrees) where solar position is evaluated. It must be in the  range [-180, 180]. 
 - <b>`algorithm`</b> ({_psa_, _nrel_, _iqbal_}):  Solar position algorithm.  _nrel_ is for the NREL's Solar Position Algorithm [1], valid  for the years -2000 to 6000. Its expected uncertainty is  +/- 0.0003 degrees. _psa_ is for the Plataforma Solar de  Almería's (PSA) algorithm [2] with updated coefficients for  the period 2020-2050 [3]. Its expected average error is 0.002  degrees, poorer than that of _nrel_ but it is much faster.  _iqbal_ is for the algorithm described in Iqbal, M. [4], which  has lower accuracy than the former algorithms, but it is  slightly faster than the _psa_. 
 - <b>`refraction`</b> (bool):  Whether atmospheric refraction must be considered. 
 - <b>`engine`</b> ({_numpy_, _numexpr_}):  Baseline code implementation to perform  the solar position calculations. _numexpr_ is expected to be faster  than numpy, especially for big spatial and/or temporal grids. 



**Returns:**
 A Sunpos instance in which the shape of all variables is (K,). 



**References:**
 

[1] Reda I and Andreas A, 2003. Solar Position Algorithm for Solar Radiation Applications. 55 pp.; NREL Report No. TP-560-34302, Revised January 2008 [pdf](http://www.nrel.gov/docs/fy08osti/34302.pdf) [url](https://midcdmz.nrel.gov/spa/). 

[2] Blanco-Muriel, M. et al. 2001. Computing the solar vector. Solar Energy, Vol. 70, pp. 431-441 doi: [10.1016/S0038-092X(00)00156-0](https://doi.org/10.1016/S0038-092X(00)00156-0). 

[3] Blanco, M. et al. 2020. Updating the PSA sun position algorithm. Solar Energy, Vol. 212, pp. 339-341 doi: [10.1016/j.solener.2020.10.084](https://doi.org/10.1016/j.solener.2020.10.084). 

[4] Iqbal, M., An introduction to solar radiation. Academic Press. 1983 [url](https://www.sciencedirect.com/book/9780123737502/an-introduction-to-solar-radiation) 




---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
