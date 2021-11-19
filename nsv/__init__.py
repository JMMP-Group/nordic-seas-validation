from importlib.metadata import PackageNotFoundError, version

from .standardizer import Standardizer

try:
    __version__ = version("nsv")
except PackageNotFoundError:
    __version__ = "unknown"

__all__ = ("Standardizer",)
