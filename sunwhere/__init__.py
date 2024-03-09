
# flake8: noqa: F401
# pylint: disable=import-error

from .version import version as __version__

from .usecases import (
    sites,
    regular_grid,
    transect
)

from ._core import __ALGORITHMS__
from .utils.datetime import universal_time_coordinated

SPA_ALGORITHMS = tuple(__ALGORITHMS__.keys())
