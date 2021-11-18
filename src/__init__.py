from importlib.metadata import PackageNotFoundError, version

from . import registries

try:
    __version__ = version("src")
except PackageNotFoundError:
    __version__ = "unknown"

__all__ = ("registries",)
