from importlib.metadata import PackageNotFoundError, version

from . import registries

try:
    __version__ = version("nsv")
except PackageNotFoundError:
    __version__ = "unknown"

__all__ = ("registries",)
