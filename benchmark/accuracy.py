
import gc
import warnings
import itertools as itt
from timeit import Timer
from pathlib import Path

import urllib3
import dateutil.parser as dtparser
from loguru import logger
from tqdm import tqdm

import numpy as np
import pylab as pl
import pandas as pd
import seaborn as sns

import pvlib.solarposition as pvsol
import sunwhere
from soltrack import SolTrack
# !python3 -m pip install git+https://github.com/jararias/SolTrack
from pysparta import SPARTA


# NOTE: keep an eye on https://github.com/MarcvdSluys/SolTrack-Python


def get_jpl_ephemeris(out_filename, start_time, stop_time, site_lon, site_lat, site_elev=0., step_size=1):
    # site_elev in km
    # step_size in minutes

    server = 'https://ssd.jpl.nasa.gov/api/horizons.api'
    format = 'format=text'
    _options = {
        'COMMAND': '10',                # the Sun
        'OBJ_DATA': 'YES',
        'MAKE_EPHEM': 'YES',
        'EPHEM_TYPE': 'OBSERVER',
        'CENTER': 'coord@399',
        # 'START_TIME': '2024-01-01',   # format: yyyy-mm-dd hh:mm
        # 'STOP_TIME': '2024-01-02',    # format: yyyy-mm-dd hh:mm
        # 'STEP_SIZE': '1 MINUTES',
        'QUANTITIES': '4,20',           # solar azimuth, elevation, and sun-earth distance (i.e., range)
        'COORD_TYPE': 'GEODETIC',       # to provide the location coordinates as lon and lat
        # 'SITE_COORD': '-3.5,37.5,0',  # location's lon, lat and altitude
        'REF_SYSTEM': 'ICRF',
        'CAL_FORMAT': 'BOTH',
        'CAL_TYPE': 'M',
        'TIME_DIGITS': 'SECONDS',
        'ANG_FORMAT': 'DEG',            # units for angles
        'APPARENT': 'AIRLESS',          # airless for NO atmospheric refraction
        'RANGE_UNITS': 'AU',            # units for the sun-earth distance
        'SUPPRESS_RANGE_RATE': 'YES',
        'SKIP_DAYLT': 'NO',
        'SOLAR_ELONG': '0,180',
        'EXTRA_PREC': 'NO',
        'R_T_S_ONLY': 'NO',
        'CSV_FORMAT': 'NO'
    }

    options = _options | {
        'START_TIME': pd.to_datetime(start_time).strftime('%Y-%m-%d %H:%M'),
        'STOP_TIME': pd.to_datetime(stop_time).strftime('%Y-%m-%d %H:%M'),
        'SITE_COORD': f'{site_lon:.6f},{site_lat:.6f},{site_elev:.3f}',
        'STEP_SIZE': f'{step_size:.0f} MINUTES'
    }

    max_lines = 90024  # max. allowed request lines !!!
    n_times = len(pd.date_range(options['START_TIME'], options['STOP_TIME'],
                                freq=f"{options['STEP_SIZE'].split()[0]}min"))
    if n_times > max_lines:
        raise ValueError(f'your request involves {n_times} data lines, which '
                         f'exceeds the max. allowed request lines of {max_lines}')

    url = server + '?' + format + '&' + '&'.join([f"{k}='{v}'" for k, v in options.items()])

    http = urllib3.PoolManager()
    response = http.request('GET', url)
    data = response.data.decode('utf-8')
    http.clear()

    path = Path(out_filename)
    if not path.parent.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
    path.open(mode='w').write(data)

    return path


def read_jpl_ephemeris_file(filename):

    def parse_line(line):
        dt64 = 'datetime64[ns]'
        dt = np.datetime64(dtparser.parse(line[1:21])).astype(dt64).astype('float64')
        jd = float(line[22:39])
        # sun_moon_presence = line[40:42]
        az = float(line[43:54])
        el = float(line[55:66])
        r = float(line[67:84])  # in AU
        return dt, jd, az, el, r

    with open(filename, 'r') as f:
        lines = f.readlines()
    first_data_line = lines.index('$$SOE\n')+1
    last_data_line = lines.index('$$EOE\n', first_data_line)
    a = np.array([parse_line(line) for line in lines[first_data_line: last_data_line]])
    return (pd.DataFrame(a)
            .set_index(pd.to_datetime(a[:, 0])).drop(columns=0)
            .set_axis(['julian_day', 'azimuth0-360', 'elevation', 'R'], axis=1)
            .rename_axis('time_utc', axis=0)
            .assign(
                azimuth=lambda df: df['azimuth0-360']-180,
                zenith=lambda df: 90.-df['elevation'],
                ecf=lambda df: 1/df['R']**2)
            .drop(columns=['azimuth0-360', 'elevation', 'R']))


def load_jpl_ephemeris(year=2024, site_lon=-3.5, site_lat=37.5):
    out_dir = Path('data')
    if not out_dir.exists():
        out_dir.mkdir(parents=True, exist_ok=True)

    whole_data_filename = out_dir / f'jpl_ephemeris_highres_{year}.parquet'
    if not whole_data_filename.exists():
        data = []
        filenames = []
        for month in range(1, 13):
            out_filename = out_dir / f'jpl_ephemeris_highres_{year}-{month:02d}.txt'
            if not out_filename.exists():
                logger.info(f'downloading high-res ephemeris to file {out_filename}')
                start_time = pd.to_datetime(f'{year}-{month:02d}-01 00:00')
                stop_time = (start_time + pd.Timedelta('35days')).replace(day=1) - pd.Timedelta('1min')
                get_jpl_ephemeris(out_filename, start_time, stop_time, site_lon, site_lat, step_size=12)
                filenames.append(out_filename)
                logger.success('data downloaded')
            logger.info(f'loading high-res ephemeris from file {out_filename}')
            data.append(read_jpl_ephemeris_file(out_filename).drop(columns='julian_day'))
        df = pd.concat(data, axis=0)
        df.insert(0, 'site_lat', site_lat)
        df.insert(0, 'site_lon', site_lon)

        df.to_parquet(whole_data_filename)

        for filename in filenames:
            filename.unlink()

    ephemeris = pd.read_parquet(whole_data_filename)
    return ephemeris.set_index(ephemeris.index.tz_localize('UTC'))


def add_level(df, label, level_names=None):
    new_columns = [[label], df.columns]
    names = level_names or ['method', 'variable']
    columns = pd.MultiIndex.from_product(new_columns, names=names)
    return df.set_axis(columns, axis=1)


def pvlib_ephemeris(times, lat, lon, methods=None):
    default_methods = ['nrel_numba', 'nrel_numpy', 'ephemeris', 'pyephem']
    methods = methods or default_methods

    def get_ephemeris(method):
        with warnings.catch_warnings():
            warnings.simplefilter(action='ignore', category=UserWarning)
            pv = pvsol.get_solarposition(times, lat, lon, method=method)
        return (pv.get(['zenith', 'azimuth'])
                .assign(azimuth=lambda df: df.azimuth - 180, ecf=float('nan')))

    return pd.concat([get_ephemeris(method) for method in methods],
                     keys=[f'pvlib.{method}' for method in methods],
                     names=['method', 'variable'], axis=1)


def sunwhere_ephemeris(times, lat, lon, methods=None, **kwargs):
    default_methods = ['nrel_numexpr', 'soltrack_numexpr', 'soltrack_numpy',
                       'psa_numexpr', 'psa_numpy', 'iqbal_numexpr', 'iqbal_numpy']
    methods = methods or default_methods

    def get_ephemeris(method, **kwargs):
        algorithm, engine = method.split('_')
        kwargs.setdefault('algorithm', algorithm)
        kwargs.setdefault('engine', engine)
        kwargs.setdefault('refraction', False)

        def get_variable(obj, name):
            df = getattr(obj, name).to_dataframe()
            return df.droplevel('location') if 'location' in df.index.names else df

        names = ['zenith', 'azimuth', 'ecf']
        sw = sunwhere.sites(times, lat, lon, **kwargs)
        return (pd.concat([get_variable(sw, name) for name in names], axis=1)
                .drop(columns=['latitude', 'longitude'])
                .rename_axis('time_utc', axis=0)
                .tz_localize('UTC'))

    return pd.concat([get_ephemeris(method, **kwargs) for method in methods],
                     keys=[f'sunwhere.{method}' for method in methods],
                     names=['method', 'variable'], axis=1)


def soltrack_ephemeris(times, lat, lon, **kwargs):
    # soltrack's original implementation ephemeris data...
    kwargs.setdefault('use_degrees', True)
    kwargs.setdefault('with_refraction', False)
    st = SolTrack(lon, lat, **kwargs)
    st.set_date_time(times)
    st.compute_position()
    st_ephemeris = {'zenith': 90. - st.altitude, 'azimuth': st.azimuth, 'ecf': 1/st.distance**2}
    return add_level(pd.DataFrame(index=times, data=st_ephemeris), 'soltrack.SolTrack')


def add_clearsky_irradiance(ephemeris, **kwargs):
    df = (ephemeris.stack(level='method')
          .swaplevel(axis=0).sort_index().reset_index())
    cosz = np.cos(np.radians(df.zenith))
    return (df.assign(**SPARTA(cosz=cosz, ecf=df.ecf, as_dict=True, **kwargs))
            .pivot(index='time_utc', columns='method')
            .swaplevel(axis=1).sort_index(axis=1))


def calc_residues(data, reference, relative=False):

    def calc_error(m, o):
        return (m - o) / o if relative is True else m - o

    residue = data.copy()
    methods = residue.columns.get_level_values('method')
    for method in methods.unique().drop(reference):
        residue[method] = calc_error(residue[method], residue[reference])
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        residue = residue.drop(columns=reference)
    return residue


def get_bulk_clearsky_errors(ephemeris, **kwargs):
    methods = {
        'sunwhere.nrel_numexpr': 'nrel',
        'sunwhere.psa_numexpr': 'psa',
        'sunwhere.iqbal_numexpr': 'iqbal',
        'sunwhere.soltrack_numexpr': 'soltrack',
        'jpl': 'jpl'
    }

    residue = calc_residues(
        # keep only the relevant sunwhere's methods...
        ephemeris.loc[:, methods.keys()].rename(columns=methods)
        .drop('azimuth', level='variable', axis=1)  # drop azimuth
        .loc[ephemeris[('jpl', 'zenith')] < 90.]  # drop night time
        .pipe(add_clearsky_irradiance, **kwargs),  # add SPARTA
        reference='jpl', relative=False
    )

    # calculate error scores...
    def named_function(f, name):
        f.__name__ = name
        return f

    scores = residue.agg(
        [
            named_function(lambda r: r.mean(), 'mbe'),
            named_function(lambda r: r.abs().mean(), 'mae'),
            named_function(lambda r: r.abs().mean(), 'std'),
            named_function(lambda r: r.pow(2).mean()**0.5, 'rmse'),
            named_function(lambda r: (r-r.mean()).abs().quantile(0.66), 'p66'),
            named_function(lambda r: (r-r.mean()).abs().quantile(0.90), 'p90')
        ]
    )

    # re-format error scores table...
    scores = (
        scores.stack(level='method')
        .rename_axis(['score', 'method'], axis=0).reset_index()
        .pivot(index='method', columns='score')
        .swaplevel(axis=1).sort_index(axis=1)
    )

    # append and re-format error scores table...
    variables = ['ghi', 'dni', 'dif']
    scores = pd.concat(
        [scores.xs(variable, level='variable', axis=1)
         .loc[['nrel', 'psa', 'soltrack', 'iqbal']] for variable in variables],
        keys=variables, names=['variable', 'score'], axis=1)

    return (scores.stack(level='variable').swaplevel(axis=0)
            .rename_axis(['variable', 'algorithm'], axis=0)
            .sort_index(level='variable', sort_remaining=False, axis=0))


def get_bulk_ephemeris_errors(ephemeris):

    residue = (calc_residues(ephemeris, reference='jpl', relative=False)
               .loc[ephemeris[('jpl', 'zenith')] < 90.])

    # calculate error scores...
    def named_function(f, name):
        f.__name__ = name
        return f

    scores = residue.agg(
        [
            named_function(lambda r: r.mean(), 'mbe'),
            named_function(lambda r: r.abs().mean(), 'mae'),
            named_function(lambda r: r.abs().mean(), 'std'),
            named_function(lambda r: r.pow(2).mean()**0.5, 'rmse'),
            named_function(lambda r: r.abs().quantile(0.66), 'p66'),
            named_function(lambda r: r.abs().quantile(0.90), 'p90')
        ]
    )

    # re-format error scores table...
    scores = (
        scores.stack(level='method')
        .rename_axis(['score', 'method'], axis=0).reset_index()
        .pivot(index='method', columns='score')
        .swaplevel(axis=1).sort_index(axis=1)
    )

    # append and re-format error scores table...
    variables = ['zenith', 'ecf', 'azimuth']
    scores = pd.concat(
        [scores.xs(variable, level='variable', axis=1) for variable in variables],
        keys=variables, names=['variable', 'score'], axis=1)

    return (scores.stack(level='variable').swaplevel(axis=0)
            .rename_axis(['variable', 'algorithm'], axis=0)
            .sort_index(level='variable', sort_remaining=False, axis=0))


def plot_accuracy(ephemeris, sort_by=None):

    df = (
        calc_residues(ephemeris, reference='jpl', relative=False)
        .stack(level='method').reset_index()
        .drop(columns=['time_utc', 'ecf'])
        .assign(
            zenith=lambda df: df.zenith.abs(),
            azimuth=lambda df: df.azimuth.abs()
        )
    )

    if sort_by is not None:
        df = (
            df.assign(mae=df.method.map(lambda m: sort_by[m]))
            .sort_values(by='mae')
            .drop(columns='mae')
        )

    methods_to_remove = ['sunwhere.nrel_numpy', 'sunwhere.psa_numpy',
                         'sunwhere.iqbal_numpy', 'sunwhere.soltrack_numpy',
                         'pvlib.nrel_numpy', 'soltrack.SolTrack']
    df = df.loc[~df.method.isin(methods_to_remove)]

    df['method'] = df.method.map(
        lambda name: {
            'pvlib.nrel_numba': "pvlib's NREL",
            'sunwhere.nrel_numexpr': "sunwhere's NREL",
            'pvlib.pyephem': "pvlib's PYEPHEM",
            'pvlib.ephemeris': "pvlib's EPHEMERIS",
            'sunwhere.iqbal_numexpr': "sunwhere's IQBAL",
            'sunwhere.psa_numexpr': "sunwhere's PSA",
            'sunwhere.soltrack_numexpr': "sunwhere's SOLTRACK",
        }.get(name, name))

    pl.rcParams['axes.labelsize'] = 'xx-large'
    sns.set_theme(style='white', rc={'axes.facecolor': (0, 0, 0, 0)})

    methods = df.method.unique()
    pal1 = sns.cubehelix_palette(len(methods), rot=-.25, light=.7)
    pal2 = sns.cubehelix_palette(len(methods))

    grid = sns.FacetGrid(df, row='method', hue='method',
                         palette=pal1, height=1.2, aspect=15)

    grid.map(sns.kdeplot, 'azimuth', log_scale=True, bw_adjust=.8, clip_on=True,
             fill=True, alpha=0.4, linewidth=0.5, warn_singular=False,
             hue=df['method'], palette=pal2)

    grid.map(sns.kdeplot, 'zenith', log_scale=True, bw_adjust=.8, clip_on=True,
             fill=True, alpha=0.4, linewidth=0.5, warn_singular=False,
             hue=df['method'], palette=pal1)

    grid.refline(y=0, linewidth=0.8, linestyle='-', color=None, clip_on=False)

    def update_label(x, color, label):
        ax = pl.gca()
        ax.text(0, .25, label, fontweight='bold', color=color,
                ha='left', va='center', transform=ax.transAxes,
                fontsize='large')

    grid.map(update_label, 'zenith')
    grid.set(yticks=[], xlim=(5e-7, 1e0), ylabel='', title='')
    grid.despine(bottom=True, left=True)

    grid.tick_params(axis='x', which='minor')
    grid.figure.subplots_adjust(hspace=-.25)

    pl.tick_params(bottom=True, labelsize='large')
    pl.xlabel('Absolute error (º)', fontsize='x-large')
    pl.suptitle('Absolute solar zenith and azimuth angle differences (bluish and reddish, '
                f'respectively) against JPL Horizons System ephemeris ({len(df)} samples)',
                fontsize='x-large')
    pl.savefig('../assets/accuracy_benchmark.png')


def speed_meter(library, algorithm, n_times, n_locs=1, n_repeats=3, number=1, quiet=True):
    """
    library: {pvlib, sunwhere}
    algorithm: the method for pvlib, or algorithm for sunwhere
    n_times: number of time steps to simulate
    n_locs: can be int or str. If int, it refers to the number of locations (as
      in sites in sunwhere). If str, it must be `intxint`, where the first integer
      refers to the number of lats and the 2nd to the number of lons to make a
      regular grid (as in regular_grid in sunwhere)
    n_repeats: number of timeit repetitions
    number: number of runs per repeat
    """
    if not quiet:
        print(f'{library}.{algorithm}({n_times}, {n_locs}): ', end='', flush=True)
    rg = np.random.default_rng()
    lag_days = rg.integers(0, 1000, size=1).item()
    t0 = np.datetime64('2024-01-01 00:00:00') + np.timedelta64(lag_days, 'D')
    times = pd.date_range(t0, freq='1min', periods=n_times)

    if isinstance(n_locs, int):
        site_lons = rg.uniform(-180., 180., n_locs)
        site_lats = rg.uniform(-90., 90., n_locs)
        if n_locs == 1:
            site_lons = site_lons.item()
            site_lats = site_lats.item()
        n_total_locs = n_locs

    if isinstance(n_locs, str):
        n_lats = int(n_locs.split('x')[0])
        n_lons = int(n_locs.split('x')[1])
        lon0 = rg.uniform(-175, 170, 1).item()
        lat0 = rg.uniform(-85, 80, 1).item()
        site_lons = np.linspace(lon0, lon0+5, n_lons)
        site_lats = np.linspace(lat0, lat0+5, n_lats)
        n_total_locs = n_lats*n_lons

    if library == 'sunwhere':
        algo, engine = algorithm.split('_')
        setup = f'kwargs=dict(engine="{engine}")'
        if isinstance(n_locs, int):
            predicate = (f'sunwhere.sites(times, lats, lons, '
                         f'refraction=False, algorithm="{algo}", **kwargs)')
        if isinstance(n_locs, str):
            predicate = (f'sunwhere.regular_grid(times, lats, lons, '
                         f'refraction=False, algorithm="{algo}", **kwargs)')

    if library == 'pvlib':
        setup = ('kwargs=dict()' if algorithm != 'nrel_numba'
                 else f'kwargs=dict(numthreads={min(8, n_times)})')
        if isinstance(n_locs, int):
            if n_locs > 1:
                predicate = (f'[pvsol.get_solarposition(times, lat, lon, '
                             f'method="{algorithm}", **kwargs) for lat, lon in zip(lats, lons)]')
            else:
                predicate = (f'pvsol.get_solarposition(times, lats, '
                             f'lons, method="{algorithm}", **kwargs)')
        if isinstance(n_locs, str):
            predicate = (f'[[pvsol.get_solarposition(times, lat, lon, '
                         f'method="{algorithm}", **kwargs) for lat in lats] for lon in lons]')

    variables = {'times': times, 'lats': site_lats, 'lons': site_lons}
    timer = Timer(predicate, globals=variables | globals(), setup=setup)

    # n_loops, time_taken = timer.autorange()
    # mean_time = time_taken / (n_loops*n_times*n_total_locs)*1e6
    # std_time = 0.

    lapses = np.array(timer.repeat(n_repeats, number))/(number*n_times*n_total_locs)
    mean_time = np.mean(lapses)*1e6
    std_time = np.std(lapses)*1e6

    del timer
    gc.collect()

    if not quiet:
        print(f'{mean_time:.3f} \u00B1 {std_time:.3f} us (coef. var:{std_time/mean_time:.1%})')
    return mean_time, std_time  # us / sim-step


def load_speed_tests():

    out_dir = Path('data')
    if not out_dir.exists():
        out_dir.mkdir(parents=True, exist_ok=True)

    library_and_algorithms = {
        'sunwhere': ['nrel_numexpr',
                     'soltrack_numexpr', 'soltrack_numpy',
                     'psa_numexpr', 'psa_numpy',
                     'iqbal_numexpr', 'iqbal_numpy'],
        'pvlib': ['nrel_numba', 'nrel_numpy']  # , 'pyephem', 'ephemeris']
    }  # pyephem and ephemeris are too much slow as to be considered !!

    time_steps = [1, 10, 100, 1000, 10000, 50000, 100000, 500000]
    locations = [1, 10, 100, '10x10', '50x50', '100x100']

    # time_steps = [1, 10]
    # locations = [1, 10, '4x4']

    def iter_dataframe(df):
        for idx in df.index:
            for col in df.columns:
                if (idx >= 10000) and (col in ('50x50', '100x100')):
                    continue
                yield int(idx), col

    keys = []
    exec_time = []
    with warnings.catch_warnings():
        warnings.simplefilter(action='ignore', category=UserWarning)
        for library, algorithms in library_and_algorithms.items():
            for algorithm in algorithms:
                this_filename = out_dir / f'exec_time_us_{library}_{algorithm}.parquet'
                if not this_filename.exists():
                    df = pd.DataFrame(index=time_steps, columns=locations)
                    kwargs = {'desc': f'{library}.{algorithm}',
                              'total': len(list(iter_dataframe(df)))}
                    pbar = tqdm(iter_dataframe(df), **kwargs)
                    for n_times, n_locs in pbar:
                        pbar.set_postfix({'n_times': n_times, 'n_locs': n_locs})
                        mean_time, std_time = speed_meter(
                            library, algorithm, n_times, n_locs)
                        df.loc[n_times][n_locs] = mean_time
                    df.to_parquet(this_filename)
                keys.append(f'{library}.{algorithm}')
                exec_time.append(pd.read_parquet(this_filename))

        return (pd.concat(exec_time, keys=keys, axis=1, names=['algorithm', 'n_locs'])
                .rename_axis('n_times', axis=0))


def plot_exec_time_and_accuracy(ephemeris, cs_errors):

    error_metrics = get_bulk_ephemeris_errors(ephemeris)

    exec_time = load_speed_tests().mul(1e-6)  # seconds
    exec_time = (
        exec_time.stack(level=['n_locs', 'algorithm'])
        .rename('t_rel (s)').reset_index())
    exec_time = exec_time.loc[exec_time.n_locs.isin(['1', '100', '10x10', '50x50'])]
    exec_time = exec_time.loc[exec_time.n_times.isin([1, 1000, 50000, 500000])]
    n_locs = exec_time['n_locs'].str.split('x').map(lambda seq: np.array(seq, dtype=int).prod())
    exec_time['t_tot (s)'] = exec_time['t_rel (s)']*exec_time['n_times']*n_locs

    mae = error_metrics.xs('zenith', level='variable', axis=0)['mae']
    ghi_errors = cs_errors.xs('ghi', level='variable', axis=0)

    def get_algorithm_ci(name):
        if 'nrel' in name.lower():
            return ghi_errors.loc['nrel'][['mbe', 'p90']].tolist()
        if 'psa' in name.lower():
            return ghi_errors.loc['psa'][['mbe', 'p90']].tolist()
        if 'soltrack' in name.lower():
            return ghi_errors.loc['soltrack'][['mbe', 'p90']].tolist()
        if 'iqbal' in name.lower():
            return ghi_errors.loc['iqbal'][['mbe', 'p90']].tolist()
        return 0, 0

    benchmark = exec_time.assign(
        zenith_mae=exec_time.algorithm.map(lambda x: mae.loc[x]),
        # cs_ci=exec_time.algorithm.map(get_algorithm_class)
    )
    print(benchmark)

    def iterate_dataframe(df):
        list_of_times = df.n_times.unique()
        list_of_locs = df.n_locs.unique()
        for k, (n_times, n_locs) in enumerate(itt.product(list_of_times, list_of_locs)):
            subset = df.loc[(df.n_times == n_times) & (df.n_locs == n_locs)]
            yield k, n_times, n_locs, subset

    pl.rcParams['axes.titlesize'] = 11
    pl.rcParams['axes.titlepad'] = 2
    pl.rcParams['axes.labelsize'] = 12

    sns.set_style('darkgrid')
    w, h = pl.figaspect(1)
    fig, axes = pl.subplots(4, 4, figsize=(w*3, h*3))
    pl.subplots_adjust(0.05, 0.04, 0.99, 0.98, wspace=0.25, hspace=0.25)

    algorithms = benchmark.algorithm.unique().tolist()
    colors = itt.cycle(sns.color_palette('Set1', len(algorithms)))
    for k, n_times, n_locs, subset in iterate_dataframe(benchmark):
        ax = axes[k // 4, k % 4]
        if subset.size == 0:
            ax.set_visible(False)
            continue

        for algorithm in algorithms:
            ci = get_algorithm_ci(algorithm)
            ax.scatter(data=subset.loc[subset.algorithm == algorithm],
                       x='t_tot (s)', y='zenith_mae', marker='o',
                       s=100, color=next(colors), edgecolors='white',
                       label=f'{algorithm} ({ci[0]:.2f} \u00b1 {ci[1]:.2f})')

        ax.set(xscale='log', yscale='log',
               title=f'n_times={n_times}, n_locs={n_locs}')

    ymin, ymax = np.inf, -np.inf
    for row_of_axes in axes:
        # xmin = min([ax.get_xlim()[0] for ax in row_of_axes])
        # xmax = max([ax.get_xlim()[1] for ax in row_of_axes])
        for ax in row_of_axes:
            ax.set_xlim(2e-3, 6e2)  # xmin, xmax)
    ymin = min([ax.get_ylim()[0] for ax in axes.flatten()])
    ymax = max([ax.get_ylim()[1] for ax in axes.flatten()])
    for ax in axes.flatten():
        ax.set_ylim(max(1e-5, ymin), ymax)

    for ax in axes[-1, :]:
        ax.set_xlabel('Execution time (s)')
    for ax in axes[:, 0]:
        ax.set_ylabel('Zenith Mean Abs. Error (º)')

    axes[2, 2].legend(bbox_to_anchor=(1.1, 0.5), loc='center left',
                      title='Algorithm (GHI MBE \u00b1 GHI CI 90%) W/m$^2$:')

    pl.savefig('../assets/accu_speed_benchmark.png')


if __name__ == '__main__':

    # load JPL Horizons ephemeris...
    ephemeris = load_jpl_ephemeris()
    times = ephemeris.index
    site_lon = ephemeris.pop('site_lon').unique().item()
    site_lat = ephemeris.pop('site_lat').unique().item()

    # append ephemeris calculated with pv-lib, sunwhere and SolTrack...
    ephemeris = (
        ephemeris.pipe(add_level, 'jpl')
        .join(pvlib_ephemeris(times, site_lat, site_lon))
        .join(sunwhere_ephemeris(times, site_lat, site_lon))
        .join(soltrack_ephemeris(times, site_lat, site_lon))
    )

    # ephemeris.to_parquet('kk.parquet')
    # ephemeris = pd.read_parquet('kk.parquet')

    # calculate ephemeris error scores...
    error_metrics = get_bulk_ephemeris_errors(ephemeris)
    print(error_metrics, '\n')

    # plot the zenith and azimuth absolute error distributions...
    plot_accuracy(ephemeris, sort_by=error_metrics.loc[('zenith', 'mae')])

    # calculate the clear-sky irradiance difference scores...
    cs_bulk_errors = get_bulk_clearsky_errors(ephemeris)  # in W/m2
    print(cs_bulk_errors, '\n')

    plot_exec_time_and_accuracy(ephemeris, cs_bulk_errors)
