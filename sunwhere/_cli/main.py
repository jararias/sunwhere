
import typer

import pylab as pl

from .position import calculate_position
# from .benchmark import run_benchmark
from .solar_chart import make_solar_chart


MY_SITE_LAT = 36.949
MY_SITE_LON = -3.822


app = typer.Typer(
    help="sunwhere's CLI",
    add_completion=False
)


def parse_latitude(lat_str):
    hemisphere = lat_str[-1].casefold()
    assert hemisphere in ('n', 's')
    latitude = (-1, 1)[hemisphere == 'n'] * float(lat_str[:-1])
    assert -90 <= latitude <= 90
    return latitude


def parse_longitude(lon_str):
    hemisphere = lon_str[-1].casefold()
    assert hemisphere in ('e', 'w')
    longitude = (-1, 1)[hemisphere == 'e'] * float(lon_str[:-1])
    assert -180 <= longitude < 180
    return longitude


@app.command(help="computes solar position in a time instant at a given location")
def at(
    time: str = typer.Argument('now', metavar='[TIME]', help="timestamp"),
    site_lat: str = typer.Argument("36.949N", parser=parse_latitude, help="latitude [90S, 90N]"),
    site_lon: str = typer.Argument("3.822W", parser=parse_longitude, help="longitude [180W, 180E)"),
    algorithm: str = typer.Option('psa', "-a", "--algorithm", help="solar position algorithm"),
    refraction: bool = typer.Option(True, help="consider atmospheric refraction"),
    timezone: str = typer.Option('UTC', "-z", "--timezone", help="time zone"),
):

    try:
        calculate_position(time, site_lat, site_lon, algorithm, refraction, timezone)
    except Exception as exc:
        typer.echo(str(exc))
        raise typer.Exit(code=1)


@app.command(help="performs a benchmark against other solar position packages [NOT AVAILABLE YET]")
def benchmark(
    year: int = typer.Argument(2024, help="benchmark year"),
    site_lat: float = typer.Argument(MY_SITE_LAT, min=-90, max=90, help="latitude [-90, 90]"),
    site_lon: float = typer.Argument(MY_SITE_LON, min=-180, max=180, help="longitude [-180, 180)"),
    plot_accuracy: bool = typer.Option(False, help="show plot of solar position algorithms accuracy"),
    plot_exec_time: bool = typer.Option(False, help="show plot of total execution times")
):
    # try:
    #     run_benchmark(year, site_lat, site_lon, plot_accuracy, plot_exec_time)
    #     if plot_accuracy or plot_exec_time:
    #         pl.show()
    # except Exception as exc:
    #     typer.echo(str(exc))
    #     raise typer.Exit(code=1)
    pass


@app.command(help="plots a solar chart [NOT AVAILABLE YET]")
def chart(
    site_lat: str = typer.Argument("36.949N", parser=parse_latitude, help="latitude [90S, 90N]"),
    site_lon: str = typer.Argument("3.822W", parser=parse_longitude, help="longitude [180W, 180E)"),
    datetime: str = typer.Option(None, '-t', '--time', metavar='[DATETIME]', help="timestamp"),
    timezone: str = typer.Option('UTC', "-z", "--timezone", help="time zone"),
    polar: bool = typer.Option(True, help="polar chart"),
    transparent: bool = typer.Option(False, help="transparent background"),
    filename: str = typer.Option(None, "-f", "--file", help="output filename"),
):

    make_solar_chart(site_lat, site_lon, datetime,
                     timezone, polar, transparent, filename)
    if filename is None:
        pl.show()


if __name__ == '__main__':

    app()
