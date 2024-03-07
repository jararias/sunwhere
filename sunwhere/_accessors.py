
import numpy as np
import pandas as pd
import xarray as xr


class SeriesAccessor:
    def __init__(self, sunpos_obj):
        self._sunpos_obj = sunpos_obj

    def __getattr__(self, item):
        times = self._sunpos_obj.times_utc
        value = getattr(self._sunpos_obj, item)
        if callable(value):
            value = value()
        return pd.Series(data=value, index=times, name=item)


class DataArrayAccessor:
    def __init__(self, sunpos_obj):
        self._sunpos_obj = sunpos_obj

    def __getattr__(self, item):
        sp = self._sunpos_obj

        data = getattr(sp, item)
        if callable(data):
            data = data()

        times = sp.times.to_numpy()

        if sp.simulation_grid == 'sites':
            if data.ndim == 1:
                data = data[:, None]
            dims = ['time', 'location']
            coords = {
                'time': times,
                'location': np.arange(len(sp.latitude)),
                'latitude': ('location', sp.latitude),
                'longitude': ('location', sp.longitude)
            }

        if sp.simulation_grid == 'regular_grid':
            dims = ['time', 'latitude', 'longitude']
            coords = {
                'time': times,
                'latitude': sp.latitude,
                'longitude': sp.longitude
            }

        # TODO: async_regular_grid and async_nonregular_grid !!

        attrs = {
            'algorithm': sp.algorithm,
            'simulation_grid': sp.simulation_grid,
            'engine': sp.engine}

        return xr.DataArray(data, dims=dims, coords=coords,
                            name=item, attrs=attrs)
