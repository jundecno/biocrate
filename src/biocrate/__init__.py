"""Common bioinformatics utilities.

The top-level package is intentionally lightweight. Import concrete helpers
from their submodules, for example ``biocrate.general`` or ``biocrate.metrics``.
This keeps optional chemistry/structure dependencies from being imported unless
the caller actually needs them.
"""

from importlib import import_module
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("biocrate")
except PackageNotFoundError:
    __version__ = "0.0.0"

_EXPORT_MODULES = (
    "general",
    "metrics",
    "constants",
    "convertors",
    "operators",
    "fetch",
    "calculators",
)

__all__ = ["__version__"]


def __getattr__(name: str):
    for module_name in _EXPORT_MODULES:
        module = import_module(f"{__name__}.{module_name}")
        try:
            value = getattr(module, name)
        except AttributeError:
            continue
        globals()[name] = value
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
