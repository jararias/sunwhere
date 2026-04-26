
import importlib.metadata

from .usecases import sites, regular_grid, transect

from ._core import __ALGORITHMS__
from .utils.datetime import universal_time_coordinated

try:
    __version__ = importlib.metadata.version("sunwhere")
except importlib.metadata.PackageNotFoundError:
    __version__ = "0.0.0"

SPA_ALGORITHMS = tuple(__ALGORITHMS__.keys())

__all__ = ["sites", "regular_grid", "transect", "universal_time_coordinated", "__version__", "SPA_ALGORITHMS"]
