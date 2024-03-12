<!-- markdownlint-disable -->

<a href="https://github.com/jararias/sunwhere/blob/master/sunwhere/_base.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square" /></a>

# <kbd>module</kbd> `_base.py`






---

## <kbd>class</kbd> `Sunpos`
Container class for [sunwhere](https://github.com/jararias/sunwhere) outputs. 

Do not instantiate!! 

<a href="https://github.com/jararias/sunwhere/blob/master/sunwhere/_base.py#L17"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square" /></a>

### <kbd>function</kbd> `__init__`

```python
__init__(
    times,
    latitude,
    longitude,
    algorithm,
    engine,
    refraction,
    usecase,
    ecf,
    eot,
    declination,
    zenith,
    azimuth
)
```

Creates a Sunpos' instance. 

Provides access to the solar geometry parameters, generally as [xarray's DataArrays](https://docs.xarray.dev/en/stable/generated/xarray.DataArray.html). 



**Args:**
 
 - <b>`times`</b> (sequence of numpy datetime64, or convertible to it):  times  where solar geometry is evaluated. 
 - <b>`latitude`</b> (sequence of floats):  latitude of the locations where  solar geometry is evaluated. Must be in the range [-90, 90]. 
 - <b>`longitude`</b> (sequence of floats):  longitude of the locations where  solar geometry is evaluated. Must be in the range [-180, 180). 
 - <b>`algorithm`</b> (str):  solar position algorithm: nrel, psa or  iqbal. 
 - <b>`engine`</b> (str):  code implementation: numpy or numexpr 
 - <b>`refraction`</b> (bool):  whether atmospheric refraction is considered 
 - <b>`usecase`</b> (str):  sites, regular_grid or transect 
 - <b>`ecf`</b> (1-dim array-like of floats):  sun-earth distance (eccentricity)  correction factor. 
 - <b>`eot`</b> (1-dim array-like of floats):  equation of time, minutes 
 - <b>`declination`</b> (1-dim array-like of floats):  solar declination,  degrees 
 - <b>`zenith`</b> (array of floats):  solar zenith angle, degrees (1-dim for  transect, 2-dim for sites and 3-dim for regular_grid) 
 - <b>`azimuth`</b> (array of floats):  solar azimuth angle, degrees (1-dim for  transect, 2-dim for sites and 3-dim for regular_grid) 



**Raises:**
 
 - <b>`ValueError`</b>:  if the inputs are not of the proper type or shape 


---

#### <kbd>property</kbd> algorithm

str: Solar position algorithm. 

---

#### <kbd>property</kbd> azimuth

DataArray: solar azimuth angle, in degrees [-180°, 180°], zero south. 

---

#### <kbd>property</kbd> cosz

DataArray: cosine of solar zenith angle. 

---

#### <kbd>property</kbd> dec

DataArray: solar declination, in degrees. Alias for `declination`. 

---

#### <kbd>property</kbd> declination

DataArray: solar declination, in degrees. 

---

#### <kbd>property</kbd> ecf

DataArray: sun-earth orbit's eccentricity correction factor. 

---

#### <kbd>property</kbd> elevation

DataArray: solar elevation angle, in degrees [-90, 90]. 

---

#### <kbd>property</kbd> engine

str: calculation engine. 

---

#### <kbd>property</kbd> eot

DataArray: equation of time, in minutes. 

---

#### <kbd>property</kbd> has_refraction

bool: whether solar zenith angle consideres atmospheric refraction. 

---

#### <kbd>property</kbd> latitude

DataArray: Latitudes where solar geometry is evaluated, degrees. 

---

#### <kbd>property</kbd> local_standard_time

Array of datetime64: Times where solar geometry is evaluated. 

---

#### <kbd>property</kbd> longitude

DataArray: Longitudes where solar geometry is evaluated, degrees. 

---

#### <kbd>property</kbd> saa

DataArray: solar azimuth angle, in degrees [-180°, 180°], zero south. Alias for `azimuth`. 

---

#### <kbd>property</kbd> sza

DataArray: solar zenith angle, in degrees [0, 180]. Alias for `zenith`. 

---

#### <kbd>property</kbd> times

Array of datetime64: Times where solar geometry is evaluated. 

---

#### <kbd>property</kbd> times_utc

Array of datetime64: Universal coordinated times where solar geometry is evaluated. 

---

#### <kbd>property</kbd> timezone

str: times time zone. 

---

#### <kbd>property</kbd> true_solar_time

DataArray: True solar time, known also to as local apparent time. 

---

#### <kbd>property</kbd> usecase

str: usage case. 

---

#### <kbd>property</kbd> zenith

DataArray: solar zenith angle, in degrees [0, 180]. 



---

<a href="https://github.com/jararias/sunwhere/blob/master/sunwhere/_base.py#L362"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square" /></a>

### <kbd>function</kbd> `airmass`

```python
airmass(parameterization='gueymard_2003')
```

Atmosphere relative optical air mass. 

Calculates the relative optical air mass of the atmosphere according to various parameterizations 



**Args:**
 
 - <b>`parameterization`</b> (str):  `kasten_1965` [1], `kasten_and_young_1989` [2], `gueymard_2001` [3] and `gueymard_2003` [4] 



**Returns:**
 A xarray's DataArray 



**Raises:**
 
 - <b>`ValueError`</b>:  if the parameterization is unknown 



**References:**
 

[1] Table V in Kasten, F. (1965) A new table and approximation formula for the relative optical air mass. CRREL Tech. Report 136 

[2] Kasten, F. and Young, A.T. (1989) Revised optical air mass tables and approximation formula. Applied Optics. Vol. 28(22), pp. 4735-4738. doi 10.1364/AO.28.004735 

[3] Table A.1 in Gueymard, C.A. (2001) Parameterized transmittance model for direct bean and circumsolar spectral irradiance. Sol Energy, Vol. 71(5), pp. 325-346. doi 10.1016/S0038-092X(01)00054-8 

[4] Eq. B.8 in Gueymard, C.A. (2003) Direct solar transmittance and irradiance predictions with broadband models. Part I - Detailed theoretical performance assessment. Sol Energy, Vol. 74, pp. 355-379 doi 10.1016/S0038-092X(03)00195-6 

---

<a href="https://github.com/jararias/sunwhere/blob/master/sunwhere/_base.py#L532"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square" /></a>

### <kbd>function</kbd> `daylight_length`

```python
daylight_length()
```

Daylight length, in hours. 



**References:**
 

 [1] Eq. (1.5.5) in Iqbal, M., An Introduction to Solar Radiation,  Academic Press, 1983. 

---

<a href="https://github.com/jararias/sunwhere/blob/master/sunwhere/_base.py#L306"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square" /></a>

### <kbd>function</kbd> `eth`

```python
eth(ISC=1361.1, am_correction='none')
```

Extraterrestrial horizontal solar irradiance. 

Calculates the extraterrestrial horizontal solar irradiance, in W m^-2, with optional corrections to account for solar eclipse obscuration and for the effect of high optical airmass if eth is used to evaluate clearness index 



**Args:**
 
 - <b>`ISC`</b> (float):  Solar constant, in W m<supersript>-2</superscript> 
 - <b>`am_correction`</b> (str):  whether use airmass correction. Useful if eth is  intended to evaluate the clearness index. Allowed values are `none`,  `perez_1990` [1], `skarveith_1998` [2] and  `gonzalez_and_calbo_1999` [3] 



**Returns:**
 A xarray's DataArray 



**Raises:**
 
 - <b>`ValueError`</b>:  if am_correction is unknown 



**References:**
 

[1] Perez et al. (1990) Making full use of the clearness index for parameterizing hourly insolation conditions. Sol Energy, Vol. 45(2), pp 111-114. doi 10.1016/0038-092X(90)90036-C 

[2] Skartveit et al. (1998) An hourly diffuse fraction model with correction for variability and surface albedo. Sol Energy, Vol. 63(3), pp. 173-183. doi 10.1016/S0038-092X(98)00067-X 

[3] González and Calbó (1999) Influence of the global radiation variability on the hourly diffuse fraction correlations. Sol Energy, Vol. 65(2), pp. 119-131. doi 10.1016/S0038-092X(98)00121-2 

---

<a href="https://github.com/jararias/sunwhere/blob/master/sunwhere/_base.py#L545"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square" /></a>

### <kbd>function</kbd> `incidence`

```python
incidence(sfc_slope, sfc_azimuth)
```

Cosine of the solar incidence angle. 

Calculates the cosine of the solar incidence angle. 



**Args:**
 
 - <b>`sfc_slope`</b> (float):  surface slope, in degrees 
 - <b>`sfc_azimuth`</b> (float):  surface orientation, in degrees east positive 



**References:**
 

[1] Eq. (1.6.5a) in Iqbal, M., An Introduction to Solar Radiation, Academic Press, 1983. 

---

<a href="https://github.com/jararias/sunwhere/blob/master/sunwhere/_base.py#L423"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square" /></a>

### <kbd>function</kbd> `sunrise`

```python
sunrise(units='deg')
```

Sunrise. 

Calculates the sunrise angle or sunrise time 



**Args:**
 
 - <b>`units`</b> (str):  `deg` for sunrise angle in degrees, `rad` for radians, `tst` for true solar time sunrise, `utc` for UTC and `local` for local time sunrise 



**Returns:**
 A xarray's DataArray 



**Raises:**
 
 - <b>`ValueError`</b>:  if units are unknown 



**References:**
 

[1] Eq. (1.5.4) in Iqbal, M., An Introduction to Solar Radiation, Academic Press, 1983. 

---

<a href="https://github.com/jararias/sunwhere/blob/master/sunwhere/_base.py#L492"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square" /></a>

### <kbd>function</kbd> `sunset`

```python
sunset(units='deg')
```

Sunset. 

Calculates the sunset angle or sunset time 



**Args:**
 
 - <b>`units`</b> (str):  `deg` for sunset angle in degrees, `rad` for radians, `tst` for true solar time sunset, `utc` for UTC and `local` for local time sunset 



**Returns:**
 A xarray's DataArray 



**Raises:**
 
 - <b>`ValueError`</b>:  if units are unknown 



**References:**
 

[1] Eq. (1.5.4) in Iqbal, M., An Introduction to Solar Radiation, Academic Press, 1983. 




---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
