"""Common bioinformatics utilities with lazy top-level exports."""

from importlib import import_module
from importlib.metadata import PackageNotFoundError, version
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .calculators import *
    from .operators import *
    from .general import *
    from .convertors import *
    from .fetch import *
    from .metrics import *

try:
    __version__ = version("biocrate")
except PackageNotFoundError:
    __version__ = "0.0.0"

_EXPORT_MODULES = (
    "calculators",
    "operators",
    "general",
    "convertors",
    "fetch",
    "metrics",
)

__all__ = ["__version__"]


def __getattr__(name: str):
    for module_name in _EXPORT_MODULES:
        try:
            module = import_module(f"{__name__}.{module_name}")
        except ModuleNotFoundError:
            continue
        try:
            value = getattr(module, name)
        except AttributeError:
            continue
        globals()[name] = value
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
