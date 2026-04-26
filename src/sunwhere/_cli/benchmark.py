
import gc
import sys
import tempfile
import warnings
from timeit import Timer
from pathlib import Path

from tqdm import tqdm
from loguru import logger
from tabulate import tabulate

import numpy as np
import pylab as pl
import pandas as pd
import seaborn as sns

import pvlib.solarposition as pvsol
# !python3 -m pip install git+https://github.com/jararias/SolTrack
from soltrack import SolTrack
from sg2 import sun_position as sg2_solpos
from pysparta import SPARTA

import sunwhere
from sunwhere.utils.jpl_horizons import (
    get_jpl_ephemeris,
    read_jpl_ephemeris_file
)


logger.remove()
logger.add(sys.stdout, format="<lvl>{message}</lvl>", level="INFO")


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


def sg2_ephemeris(times, lats, lons):
    _lons = np.array(lons, ndmin=1, dtype='float64')
    _lats = np.array(lats, ndmin=1, dtype='float64')
    geopoints = np.c_[_lons, _lats, np.zeros(len(_lats))]
    times = np.array(times.tz_localize(None), ndmin=1, dtype='datetime64[ns]')
    outputs = ["geoc.EOT", "geoc.R", "topoc.delta", "topoc.gamma_S0", "topoc.alpha_S"]
    sp = sg2_solpos(geopoints, times, outputs)
    # eot = (sp.geoc.EOT/np.pi - np.floor(0.5 + sp.geoc.EOT/np.pi))*12*60
    return pd.DataFrame(
        index=times,
        data={
            'zenith': 90. - np.degrees(sp.topoc.gamma_S0.T).ravel(),
            'azimuth': np.degrees(sp.topoc.alpha_S.T).ravel() - 180,
            # 'dec': np.degrees(sp.topoc.delta.T).ravel(),
            'ecf': 1 / (sp.geoc.R**2).ravel(),
            # 'eot': eot.ravel(),  # minutes
        }
    ).rename_axis('time_utc', axis=0).tz_localize('UTC')


def add_clearsky_irradiance(ephemeris, **kwargs):
    df = (ephemeris.stack(level='method')
          .swaplevel(axis=0).sort_index().reset_index())
    cosz = np.cos(np.radians(df.zenith))
    return (df.assign(**SPARTA(cosz=cosz, ecf=df.ecf, as_dict=True, **kwargs))
            .pivot(index='time_utc', columns='method')
            .swaplevel(axis=1).sort_index(axis=1))


def calc_residues(data, reference, relative=False):

    def calc_error(m, o):
        diff = m - o
        if 'azimuth' in diff.columns:
            diff['azimuth'] = diff['azimuth'].where(diff['azimuth'] > -330, -360 - diff['azimuth'])
            diff['azimuth'] = diff['azimuth'].where(diff['azimuth'] < 330, 360 - diff['azimuth'])
        return diff / o if relative is True else diff

    residue = data.copy()
    methods = residue.columns.get_level_values('method')
    for method in methods.unique().drop(reference):
        residue[method] = calc_error(data[method], data[reference])
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        residue = residue.drop(columns=reference)
    return residue


def get_bulk_ephemeris_errors(ephemeris):

    residue = (calc_residues(ephemeris, reference='jpl', relative=False)
               .loc[ephemeris[('jpl', 'zenith')] < 180.])

    # calculate error scores...
    def named_function(f, name):
        f.__name__ = name
        return f

    scores = residue.agg(
        [
            named_function(lambda r: r.mean(), 'mbe'),
            named_function(lambda r: r.std(), 'std'),
            named_function(lambda r: r.abs().mean(), 'mae'),
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
    variables = ['zenith', 'ecf', 'azimuth']
    scores = pd.concat(
        [scores.xs(variable, level='variable', axis=1) for variable in variables],
        keys=variables, names=['variable', 'score'], axis=1)

    return (scores.stack(level='variable').swaplevel(axis=0)
            .rename_axis(['variable', 'algorithm'], axis=0)
            .sort_index(level='variable', sort_remaining=False, axis=0))


def get_bulk_clearsky_errors(ephemeris, **kwargs):
    required_methods_mapping = {
        'sunwhere.nrel_numexpr': 'nrel',
        'sunwhere.psa_numexpr': 'psa',
        'sunwhere.iqbal_numexpr': 'iqbal',
        'sunwhere.soltrack_numexpr': 'soltrack',
        'sg2.sg2_c++': 'sg2',
        'jpl': 'jpl'
    }

    available_methods = ephemeris.columns.get_level_values('method').unique()
    methods = available_methods.intersection(list(required_methods_mapping.keys()))
    residue = calc_residues(
        # keep only the relevant sunwhere's methods...
        ephemeris.loc[:, methods]
        .rename(columns=required_methods_mapping)
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
            named_function(lambda r: r.std(), 'std'),
            named_function(lambda r: r.abs().mean(), 'mae'),
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
    available_methods = scores.index.intersection(
        ['nrel', 'psa', 'soltrack', 'iqbal', 'sg2'])
    scores = pd.concat(
        [scores.xs(variable, level='variable', axis=1)
         .loc[available_methods] for variable in variables],
        keys=variables, names=['variable', 'score'], axis=1)

    return (scores.stack(level='variable').swaplevel(axis=0)
            .rename_axis(['variable', 'algorithm'], axis=0)
            .sort_index(level='variable', sort_remaining=False, axis=0))


def plot_accuracy(ephemeris, sort_by=None, filename=None):

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

    methods_to_remove = [
        'sunwhere.nrel_numpy', 'sunwhere.psa_numpy',
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
            'sg2.sg2_c++': "SG2"
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
    if filename is not None:
        pl.savefig(filename)  # '../assets/accuracy_benchmark.png')


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

    if library == 'sg2':
        setup = ('lons1d = np.array(lons, ndmin=1, dtype="float64");'
                 'lats1d = np.array(lats, ndmin=1, dtype="float64");')
        if isinstance(n_locs, int):
            predicate = ('sg2_solpos(np.c_[lons1d, lats1d, np.zeros(len(lats1d))],'
                         'np.array(times.tz_localize(None), dtype="datetime64[ns]"),'
                         '["geoc.R", "topoc.gamma_S0", "topoc.alpha_S"])')
        else:
            setup = ('xlons, xlats = np.meshgrid(lons, lats);'
                     'lons1d = xlons.reshape(-1);'
                     'lats1d = xlats.reshape(-1);')
            predicate = ('sg2_solpos(np.c_[lons1d, lats1d, np.zeros(len(lats1d))],'
                         'np.array(times.tz_localize(None), dtype="datetime64[ns]"),'
                         '["geoc.R", "topoc.gamma_S0", "topoc.alpha_S"])')

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


def load_speed_tests(out_dir=None):

    out_dir = Path(out_dir or '.')
    if not out_dir.exists():
        out_dir.mkdir(parents=True, exist_ok=True)

    library_and_algorithms = {
        'sunwhere': ['nrel_numexpr',
                     'soltrack_numexpr', 'soltrack_numpy',
                     'psa_numexpr', 'psa_numpy',
                     'iqbal_numexpr', 'iqbal_numpy'],
        'pvlib': ['nrel_numba', 'nrel_numpy'],  # pyephem and ephemeris are too slow !!
        'sg2': ['sg2_c++']
    }

    time_steps = [1, 10, 100, 1000, 10000, 50000, 100000, 500000]
    locations = [1, 10, 100, '10x10', '50x50', '100x100']

    def iter_dataframe(df):
        for idx in df.index:
            for col in df.columns:
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
                        if (n_times > 1000) and (n_locs in ('50x50', '100x100')):
                            continue
                        mean_time, std_time = speed_meter(
                            library, algorithm, n_times, n_locs)
                        df.loc[n_times][n_locs] = mean_time
                    df.to_parquet(this_filename)
                keys.append(f'{library}.{algorithm}')
                exec_time.append(pd.read_parquet(this_filename))

        return (pd.concat(exec_time, keys=keys, axis=1, names=['algorithm', 'n_locs'])
                .rename_axis('n_times', axis=0))


def plot_exec_time(filename=None):

    # out_dir = '/home/jararias/code/devel/sunwhere/benchmark/data'
    # exec_time = (load_speed_tests()  # out_dir)
    #              .stack(level=['n_locs', 'algorithm'])
    #              .rename('t_rel (us)').reset_index())
    # exec_time.to_csv('/home/jararias/code/devel/sunwhere/assets/exec_time_us.csv')

    inp_filename = '/home/jararias/code/devel/sunwhere/assets/exec_time_us.csv'
    exec_time = pd.read_csv(inp_filename, index_col=0)  # in microseconds

    fig, axes = pl.subplots(1, 3, figsize=(16, 6), layout='constrained')
    for k, n_locs in enumerate(('1', '100', '10x10')):
        df = (exec_time.loc[exec_time.n_locs == n_locs].drop(columns='n_locs')
              .pivot(index='n_times', columns='algorithm')
              .droplevel(0, axis=1))  # microseconds / ephemeris

        df.columns = pd.MultiIndex.from_tuples(
            [(e[0], *e[1].split('_')) for e in df.columns.str.split('.')],
            names=['library', 'algorithm', 'engine'])

        n_tot_locs = np.prod(n_locs.split('x'), dtype=float) if 'x' in n_locs else float(n_locs)
        total_sec = df.apply(lambda x: x * x.index * n_tot_locs * 1e-6)  # total seconds

        ax = axes[k]
        sns.heatmap(total_sec, annot=True, fmt=".2f", linewidths=0.5,
                    norm=pl.cm.colors.LogNorm(), cmap='magma_r',
                    annot_kws={'fontsize': 8}, ax=ax)
        ax.set_xlabel(None)
        ax.set_title(f'n_locs={n_locs}', fontsize=10)
        ax.tick_params(axis='y', labelrotation=0)
        fig.suptitle('Total execution time, seconds', fontsize=13)

        # print(f'\nTOTAL EXECUTION TIME FOR {n_locs} LOCATIONS (in seconds)')
        # print(
        #     tabulate(
        #         (total_sec.set_axis(['\n'.join(col) for col in total_sec.columns], axis=1)
        #         .rename_axis('\n\nn_times', axis=0)), headers='keys', floatfmt=".3f"
        #     ), end='\n'
        # )

        # rate = 1 / df  # Mill. ephemerides / second
        # print(
        #     tabulate(
        #         (rate.set_axis(['\n'.join(col) for col in rate.columns], axis=1)
        #          .rename_axis('\n\nn_times', axis=0)), headers='keys', floatfmt=".3f"
        #     ), end='\n'
        # )

    if filename is not None:
        pl.savefig(filename)  # '../assets/exec_time_benchmark.png')


def print_bulk_errors(df, variables=None, reformat_name=False):
    df = (df.get(['mbe', 'std', 'p66', 'p90', 'mae', 'rmse'])
          .rename(columns={'p66': '\u00b1CI66', 'p90': '\u00b1CI90'}))

    def reformat(name):
        if not reformat_name:
            return name
        package, method = name.split('.')
        spa, engine = method.split('_') if '_' in method else (method, None)
        new_name = f"{package}'s {spa.upper()}"
        return new_name if engine is None else new_name + f" ({engine})"

    for variable in variables or df.index.get_level_values('variable').unique():
        this_var = (df.loc[variable].sort_values(by='rmse').rename_axis(variable, axis=0)
                    .pipe(lambda x: x.set_index(x.index.map(reformat))))
        print(tabulate(this_var, headers='keys', floatfmt='.3f', tablefmt='mixed_outline'), end='\n')


def run_benchmark(year, site_lat, site_lon, show_accuracy=False, show_exec_time=False):

    logger.info(f'Performing benchmark tests for {year} at lat={site_lat} lon={site_lon}')

    sign = lambda x: int(0.5*(1 + x / abs(x)))  # noqa: E731
    lat_str = f'{abs(site_lat):07.4f}{("S", "N")[sign(site_lat)]}'.replace('.', 'p')
    lon_str = f'{abs(site_lon):08.4f}{("W", "E")[sign(site_lon)]}'.replace('.', 'p')
    out_filename = (Path(tempfile.gettempdir()) /
                    f'jpl_ephemerides_{year}_{lat_str}_{lon_str}.txt')
    logger.debug(f'JPL filename: {out_filename}')

    if not out_filename.exists():
        get_jpl_ephemeris(out_filename, f'{year}-01-01T00:00:00',
                          f'{year+1}-01-01T00:00:00', site_lon,
                          site_lat, site_elev=0., step_size=12)
        logger.success(f'JPL data downloaded to {out_filename}')

    # load JPL Horizons ephemeris file
    ephemeris = read_jpl_ephemeris_file(out_filename)
    times = ephemeris.index

    # append pvlib ephemeris
    logger.info('Calculating and appending sunwhere and pvlib ephemerides')
    ephemeris = (
        ephemeris.pipe(add_level, 'jpl')
        .join(pvlib_ephemeris(times, site_lat, site_lon))
        .join(sg2_ephemeris(times, site_lat, site_lon).pipe(add_level, 'sg2.sg2_c++'))
        .join(sunwhere_ephemeris(times, site_lat, site_lon))
    )

    # ephemeris.to_parquet('/home/jararias/kk.parquet')
    # ephemeris = pd.read_parquet('/home/jararias/kk.parquet')

    logger.info('EPHEMERIDES ERRORS (zenith and azimuth in arc-sec; ecf is unitless):')
    logger.info('  The errors are evaluated against ephemerides from the the JPL Horizons')
    logger.info('  service at https://ssd.jpl.nasa.gov/horizons/')
    ephemeris_errors = get_bulk_ephemeris_errors(ephemeris)
    # zenith and azimuth errors to arc-sec and ecf to milli...
    ephemeris_errors_arcsec = (ephemeris_errors.unstack(level='variable')
                               .swaplevel(axis=1).sort_index(axis=1))
    ephemeris_errors_arcsec['zenith'] *= 3600
    ephemeris_errors_arcsec['azimuth'] *= 3600
    ephemeris_errors_arcsec['ecf'] *= 1e3
    ephemeris_errors_arcsec = ephemeris_errors_arcsec.rename(
        columns={'ecf': 'ecf (x10\u207b\u00b3)'})
    ephemeris_errors_arcsec = (
        ephemeris_errors_arcsec.stack(level=0)
        .swaplevel(axis=0).loc[['zenith', 'azimuth', 'ecf (x10\u207b\u00b3)']])
    print_bulk_errors(ephemeris_errors_arcsec, reformat_name=True,
                      variables=('zenith', 'azimuth', 'ecf (x10\u207b\u00b3)'))

    logger.info('CLEAR-SKY IRRADIANCE ERRORS (in W m\u207b\u00b2):')
    logger.info('  The errors are evaluated using the SPARTA clear-sky solar')
    logger.info('  irradiance model with a default (standard) atmosphere and')
    logger.info('  solar geometry from JPL Horizons (true reference) and ')
    logger.info('  alternatively calculated with the solar position algorithms')
    clearsky_errors = get_bulk_clearsky_errors(ephemeris)  # in W/m2
    print_bulk_errors(clearsky_errors, variables=('ghi', 'dni', 'dif'))

    # plot the zenith and azimuth absolute error distributions...
    if show_accuracy is True:
        plot_accuracy(
            ephemeris, sort_by=ephemeris_errors.loc[('zenith', 'mae')],
            filename='/home/jararias/code/devel/sunwhere/assets/accuracy_benchmark.png'
        )

    # plot_exec_time('/home/jararias/code/devel/sunwhere/assets/exec_time_benchmark.png')
    if show_exec_time is True:
        plot_exec_time()
