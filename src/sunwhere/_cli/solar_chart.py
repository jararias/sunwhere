
import itertools as itt
from datetime import datetime

import numpy as np
import pylab as pl
import pandas as pd

import sunwhere


pl.rcParams['font.family'] = 'arial'
pl.rcParams['axes.labelsize'] = 17
pl.rcParams['axes.labelpad'] = 12
pl.rcParams['xtick.labelsize'] = 14
pl.rcParams['ytick.labelsize'] = 14


class SolarPosition:
    def __init__(self, latitude, longitude=0):
        self.latitude = latitude
        assert -90 <= latitude <= 90
        self.longitude = longitude
        assert -180 <= longitude < 180

    def __call__(self, times):
        sw = sunwhere.sites(times, self.latitude, self.longitude)
        return {'azimuth': sw.saa.isel(location=0).to_numpy(),
                'altitude': sw.elevation.isel(location=0).to_numpy()}


class SolarChart(object):

    color_solstice = 'purple'
    color_equinox = 'purple'
    color_solar_grid = 'mediumpurple'
    color_solar_patch = 'lightyellow'
    color_user_solar_path = 'red'

    def __init__(self, latitude, longitude=0, user_datetime=None,
                 timezone=None, polar=True):
        self.latitude = float(latitude)
        assert -90 < self.latitude < 90
        self.longitude = float(longitude)
        assert -180 <= self.longitude <= 180
        self._user_datetime = user_datetime
        self._timezone = timezone or 'UTC'
        self.sunpath = SolarPosition(latitude, longitude)
        self.projection = ('rectilinear', 'polar')[polar]
        figsize = ((12, 9), (10, 10))[polar]
        self.fig, self.ax = pl.subplots(
            1, 1, figsize=figsize, dpi=150, layout='constrained',
            subplot_kw={'projection': self.projection})

        if self.is_polar:
            # choose radial tick label positions
            self.ax.set_rlabel_position(180)
            # rotate South down
            self.ax.set_theta_offset(-np.pi / 2.)
            self.ax.set_ylim([90., 0.1])
        else:
            self.ax.set_xlim([-180., 180.])
            self.ax.set_ylim([0., 90.])

    @property
    def is_polar(self):
        return self.projection == 'polar'

    @property
    def timezone(self):
        return self._timezone

    def _date_range(self, start, n_periods, time_step):
        # print(start, n_periods, time_step)
        # tstamp = pd.Timestamp(start).tz_localize(self.timezone, ambiguous='NaT')
        # print(tstamp)
        times = pd.date_range(start, periods=n_periods, freq=time_step)
        times = (pd.to_datetime(times).tz_localize(None)
                 .tz_localize(self.timezone, ambiguous='NaT', nonexistent='NaT'))
        return times[~times.isna()]

    def plot_sunpath(self, sunpath, **kwargs):
        x = sunpath['azimuth']
        y = sunpath['altitude']
        self.ax.plot(np.radians(x) if self.is_polar else x, y, **kwargs)

    def text(self, azimuth, altitude, txt, **kwargs):
        az = np.radians(azimuth) if self.is_polar else azimuth
        self.ax.text(az, altitude, txt, **kwargs)

    def draw_solar_grid(self, days=None, months=None, hours=None, minutes=None,
                        include_solstice=True, include_equinox=True):

        ts = pd.Timestamp(datetime.now())
        if self._user_datetime is not None:
            ts = pd.Timestamp(self._user_datetime)
        ts = ts.tz_localize(self.timezone)
        time_ref = ts.tzname()

        zone_name = self.timezone
        if zone_name == 'UTC':
            zone_name = 'Universal Coordinated Time'

        time_step_minutes = 2
        n_steps = 1440 // time_step_minutes
        time_step = pd.Timedelta(time_step_minutes, 'min')
        kwargs = dict(ls='-', marker='')

        # solsticio de verano
        if include_solstice is True:
            times = self._date_range(f'21 Jun {ts.year}', n_steps, time_step)
            solstice = self.sunpath(times.tz_convert('UTC'))
            kwargs.update({'color': SolarChart.color_solstice})
            self.plot_sunpath(solstice, **kwargs)

        # equinoccio de invierno
        if include_equinox is True:
            times = self._date_range(f'21 Dec {ts.year}', n_steps, time_step)
            equinox = self.sunpath(times.tz_convert('UTC'))
            kwargs.update({'color': SolarChart.color_equinox})
            self.plot_sunpath(equinox, **kwargs)

        # fill space between the equinox and the solstice
        if include_solstice is True and include_equinox is True:
            x = np.linspace(-180, 180, 361)
            y1 = np.interp(x, solstice['azimuth'], solstice['altitude'], period=360)
            y2 = np.interp(x, equinox['azimuth'], equinox['altitude'], period=360)
            self.ax.fill_between(np.radians(x) if self.is_polar else x, y1,
                                 y2, color=SolarChart.color_solar_patch, zorder=-100)

        # DRAW INTRA-DAILY GRID AT hours:minutes

        color = SolarChart.color_solar_grid

        # draw grid
        hours_and_minutes = list(
            itt.product(hours or range(0, 24), minutes or (0,))
        )
        for ho, mn in hours_and_minutes:
            # draw analemmas
            date_s = pd.Timestamp(f'1 Jan {ts.year} {ho:02d}:{mn:02d}')
            times = self._date_range(date_s, 366 if ts.is_leap_year else 365, '1D')
            sunpath = self.sunpath(pd.to_datetime(times).tz_localize(None))
            self.plot_sunpath(sunpath, ls='-', marker='', lw=0.7, color=color)

            # draw time labels
            solpos = self.sunpath(
                [pd.Timestamp(f'21 Jun {ts.year} {ho:d}:{mn:d}', tz=self.timezone)])
            az, al = solpos['azimuth'], solpos['altitude']
            if al >= 0:
                daz, dal = 0., 1.
                ha = 'center' if abs(az) < 60 else ('left' if az < 0 else 'right')
                if not self.is_polar:
                    daz, dal = np.sign(az) * 1.5, 0.
                    ha = 'center' if az == 0 else ('right' if az < 0 else 'left')

                txt = f'{ho}' if mn == 0 else f'{ho}:{mn:02d}'
                bbox = dict(facecolor='white', edgecolor='none', pad=.1)
                self.text(az + daz, al + dal, txt + f' {time_ref}', ha=ha,
                          va='bottom', color=color, fontsize=8, bbox=bbox)

        # DRAW DAILY GRID AT days/months between the summer solstice and winter equinox

        days = list(itt.product(days or (15,), months or np.r_[1:6, 7:12]))
        for day in ((21, 6), (21, 12)):
            if day not in days:
                days.append(day)

        for day, month in days:
            # draw sun path
            sunpath = self.sunpath(
                self._date_range(f'{day}/{month}/{ts.year}', 288, '5min'))
            self.plot_sunpath(
                sunpath, ls='-', lw=0.5, marker='',
                color=SolarChart.color_solar_grid)

            az = 0.
            # if (self._timezone is None) and len(hours_and_minutes) > 1:
            #     sp = self.sunpath([datetime(2013, month, day, ho, mn)
            #                        for ho, mn in hours_and_minutes])
            #     ix = np.digitize(0, sp['azimuth'])
            #     az = 0.5 * (sp['azimuth'][ix-1] + sp['azimuth'][ix])

            if (day, month) == (21, 6):
                color = SolarChart.color_solstice
                weight = 'bold'
            elif (day, month) == (21, 12):
                color = SolarChart.color_equinox
                weight = 'bold'
            else:
                color = SolarChart.color_solar_grid
                weight = 'normal'

            altitude = np.interp(az, sunpath['azimuth'],
                                 sunpath['altitude'], period=360)
            bbox = dict(fc=SolarChart.color_solar_patch, ec='none', pad=1)
            self.text(az, altitude, f'{day:02d}/{month:02d}', fontsize=6,
                      ha='center', va='center', bbox=bbox, color=color,
                      weight=weight)

        # ADD OTHER LABELS AND INFORMATION

        def add_text(x, y, txt, **kwargs):
            kwargs.update({'transform': self.ax.transAxes})
            self.ax.text(x, y, txt, **kwargs)

        add_text(-.05, -.08, f'{time_ref}: {zone_name}', fontsize=10)

        if self.is_polar:
            add_text(-.08, 1.02, r'Latitude=$%.2f^{\circ}$' % self.latitude,
                     bbox=dict(fc='white', ec='none', pad=4), fontsize=14)
            add_text(-.08, .99, r'Longitude=$%.2f^{\circ}$' % self.longitude,
                     bbox=dict(fc='white', ec='none', pad=3), fontsize=14)
            add_text(0.5, 1.08, 'N', ha='center', fontsize=16)
            add_text(0.5, -.08, 'S', ha='center', fontsize=16)
            add_text(-.09, .5, 'W', va='center', fontsize=16)
            add_text(1.09, .5, 'E', va='center', fontsize=16)
        else:
            add_text(.01, .97, r'Latitude=$%.2f^{\circ}$' % self.latitude,
                     bbox=dict(fc='white', ec='none', pad=4), fontsize=14)
            add_text(.01, .94, r'Longitude=$%.2f^{\circ}$' % self.longitude,
                     bbox=dict(fc='white', ec='none', pad=4), fontsize=14)
            self.ax.set_xlabel('Solar Azimuth')
            self.ax.set_ylabel('Solar Altitude')
            add_text(0., 1.06, 'North', fontsize=16)
            add_text(1., 1.06, 'NortH', ha='right', fontsize=14)
            add_text(0.5, 1.06, 'South', ha='center', fontsize=14)
            add_text(0.25, 1.06, 'West', ha='center', fontsize=14)
            add_text(0.75, 1.06, 'East', ha='center', fontsize=14)

        self.draw_background_grid()

    def draw_background_grid(
            self, major_altitude_ticks=None,
            minor_altitude_ticks=None, major_azimuth_ticks=None,
            minor_azimuth_ticks=None):

        def get_polar_ticks(key, value):
            ticks = {
                'altitude_major': np.arange(10, 90, 10),
                'altitude_minor': np.arange(1, 90, 1),
                'azimuth_major': np.arange(-160., 181., 20),
                'azimuth_minor': np.arange(-175., 180., 5.)
            }
            return value if value is not None else ticks[key]

        def get_rect_ticks(key, value):
            ticks = {
                'altitude_major': np.arange(0, 90, 10),
                'altitude_minor': np.arange(2, 90, 2),
                'azimuth_major': np.arange(-150, 180, 30),
                'azimuth_minor': np.arange(-175, 180, 5)
            }
            return value if value is not None else ticks[key]

        def transform(azimuth):
            return np.where(azimuth >= 0., azimuth, azimuth + 360.)

        if self.is_polar:
            # major solar altitude labels
            ticks = get_polar_ticks('altitude_major', major_altitude_ticks)
            labels = [f'{y:.0f}' + r'$^{\circ}$' for y in ticks]
            self.ax.yaxis.set_ticks(ticks)
            self.ax.yaxis.set_ticklabels(labels)
            # minor solar altitude labels
            ticks = get_polar_ticks('altitude_minor', minor_altitude_ticks)
            self.ax.yaxis.set_ticks(ticks, minor=True)
            pl.setp(pl.getp(self.ax, 'yticklabels'), ha='center')
            self.ax.tick_params(axis='y', which='minor', labelright=False)
            # major solar azimuth labels
            ticks = get_polar_ticks('azimuth_major', major_azimuth_ticks)
            labels = [f'{x:.0f}' + r'$^{\circ}$' for x in ticks]
            self.ax.xaxis.set_ticks(np.radians(transform(ticks)))
            self.ax.xaxis.set_ticklabels(labels)
            self.ax.tick_params(axis='x', which='major', pad=7)
            # minor solar azimuth labels
            ticks = get_polar_ticks('azimuth_minor', minor_azimuth_ticks)
            self.ax.xaxis.set_ticks(np.radians(transform(ticks)), minor=True)
            self.ax.tick_params(axis='x', which='minor', labelbottom=False)
        else:
            # solar altitude ticks and labels
            ticks = get_rect_ticks('altitude_major', major_altitude_ticks)
            labels = [f'{x:.0f}' + r'$^{\circ}$' for x in ticks]
            self.ax.yaxis.set_ticks(ticks)
            self.ax.yaxis.set_ticklabels(labels)
            ticks = get_rect_ticks('altitude_minor', minor_altitude_ticks)
            self.ax.yaxis.set_ticks(ticks, minor=True)
            # solar azimuth ticks and labels
            ticks = get_rect_ticks('azimuth_major', major_azimuth_ticks)
            labels = [f'{x:.0f}' + r'$^{\circ}$' for x in ticks]
            self.ax.xaxis.set_ticks(ticks)
            self.ax.xaxis.set_ticklabels(labels)
            ticks = get_rect_ticks('azimuth_minor', minor_azimuth_ticks)
            self.ax.xaxis.set_ticks(ticks, minor=True)
            self.ax.tick_params(axis='x', which='minor', labelbottom=False)
            self.ax.tick_params(which='both', top=True, labeltop=True)

        self.ax.xaxis.grid(dashes=(20, 0), lw=.5, color='0.7', zorder=-100)
        self.ax.xaxis.grid(ls='--', lw=.5, color='0.8', which='minor', zorder=-100)
        self.ax.yaxis.grid(dashes=(20, 0), lw=.5, color='0.7', zorder=-100)
        self.ax.yaxis.grid(ls='--', lw=.5, color='0.8', which='minor', zorder=-100)

    def draw_solar_path(self, user_datetime):

        color = SolarChart.color_user_solar_path

        ts = pd.Timestamp(user_datetime).tz_localize(self.timezone)

        # diurnal sunpath
        times = self._date_range(f'{ts.year}/{ts.month}/{ts.day}', 288, '5min')  # ts.date(), 1440, '1min')
        sunpath = self.sunpath(times)  # pd.to_datetime(times).tz_localize(None))
        self.plot_sunpath(sunpath, ls='-', marker='', lw=1, color=color)

        # analemma
        date_s = ts.tz_convert('UTC').replace(month=1, day=1)
        times = self._date_range(date_s, 366 if ts.is_leap_year else 365, '1D')
        sunpath = self.sunpath(pd.to_datetime(times).tz_localize(None))
        self.plot_sunpath(sunpath, ls='-', marker='', lw=0.7, color=color)

        # solpos = self.sunpath([ts.tz_localize(None).tz_localize(self.timezone)])
        solpos = self.sunpath([ts.tz_convert('UTC')])
        self.plot_sunpath(solpos, ls='', marker='.', ms=8, color=color)

        txt = f'User path: {ts}'
        x = -.08 if self.is_polar else .01
        y = .977 if self.is_polar else .925
        self.ax.text(x, y, txt, transform=self.ax.transAxes,
                     ha='left', va='top', fontsize=14, color=color,
                     bbox=dict(facecolor='white', edgecolor='none', pad=3))

    def savefig(self, filename, dpi=150, **kwargs):
        self.fig.set_dpi(dpi)
        self.fig.savefig(filename, dpi=dpi, **kwargs)


def make_solar_chart(site_lat, site_lon, user_datetime=None, timezone='UTC',
                     polar=True, transparent=False, filename=None):

    solar_chart = SolarChart(site_lat, site_lon, user_datetime, timezone, polar)
    solar_chart.draw_solar_grid()
    solar_chart.draw_solar_path(user_datetime)
    if filename is not None:
        solar_chart.savefig(filename, transparent=transparent)
