
import typer

import pylab as pl

from .position import calculate_position
from .benchmark import run_benchmark


MY_SITE_LAT = 36.949
MY_SITE_LON = -3.822


app = typer.Typer()


@app.command()
def at(
    time: str = typer.Argument('now', metavar='[TIME]', help="timestamp"),
    site_lat: float = typer.Argument(MY_SITE_LAT, min=-90, max=90, help="latitude [-90, 90]"),
    site_lon: float = typer.Argument(MY_SITE_LON, min=-180, max=180, help="longitude [-180, 180)"),
    algorithm: str = typer.Option('psa', help="solar position algorithm"),
    refraction: bool = typer.Option(True, help="consider atmospheric refraction"),
    timezone: str = typer.Option('UTC', help="time zone"),
):

    try:
        calculate_position(time, site_lat, site_lon, algorithm, refraction, timezone)
    except Exception as exc:
        typer.echo(str(exc))
        raise typer.Exit(code=1)


@app.command()
def benchmark(
    year: int = typer.Argument(2024, help="benchmark year"),
    site_lat: float = typer.Argument(MY_SITE_LAT, min=-90, max=90, help="latitude [-90, 90]"),
    site_lon: float = typer.Argument(MY_SITE_LON, min=-180, max=180, help="longitude [-180, 180)"),
    plot_accuracy: bool = typer.Option(False, help="show plot of solar position algorithms accuracy"),
    plot_exec_time: bool = typer.Option(False, help="show plot of total execution times")
):
    try:
        run_benchmark(year, site_lat, site_lon, plot_accuracy, plot_exec_time)
        if plot_accuracy or plot_exec_time:
            pl.show()
    except Exception as exc:
        typer.echo(str(exc))
        raise typer.Exit(code=1)


@app.command()
def chart():
    print('solar chart')


if __name__ == '__main__':

    app()
