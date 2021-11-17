from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("src")
except PackageNotFoundError:
    # package is not installed
    pass
