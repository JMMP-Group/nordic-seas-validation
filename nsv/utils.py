from functools import wraps

import cf_xarray  # noqa: F401
from gsw import SA_from_SP, pt0_from_t
from xarray import Dataset


def final_cleanup_before_returning(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        ds = func(*args, **kwargs)

        # Set coordinates
        coords = ["longitude", "latitude", "distance", "time"]
        ds = ds.set_coords(
            [var for var in ds.data_vars for coord in coords if coord in var]
        )
        ds = ds.cf.guess_coord_axis()

        # Add cf attributes
        ds = ds.cf.add_canonical_attributes()

        # Sort stations by longitute
        ds = ds.cf.sortby("longitude")
        if "time" in ds:
            ds = ds.sortby("time")

        # Don't use strings for stations
        if "station" in ds.variables:
            ds["station_id"] = ds["station"].astype(str)
            ds = ds.drop("station")

        # Attributes
        ds.attrs[
            "description"
        ] = f"Standardized {func.__name__.upper()}{args[1:] or ''}{kwargs or ''}"
        ds.attrs["featureType"] = "timeSeries"
        return ds

    return wrapper


def add_attributes_and_rename_variables(ds: Dataset, attrs_dict: dict) -> Dataset:

    for var, da in ds.data_vars.items():
        da.encoding.pop("coordinates", None)
        ds[var] = da

    for var, attrs in attrs_dict.items():
        ds[var].attrs = {**ds[var].attrs, **attrs}

    for var, da in ds.variables.items():
        if "standard_name" in da.attrs:
            ds = ds.rename(
                {
                    var: da.attrs["standard_name"]
                    .replace("sea_water_", "")
                    .replace("in_sea_water", "")
                }
            )

    return ds


def compute_pt0(ds: Dataset) -> Dataset:

    sa = SA_from_SP(
        ds.cf["sea_water_practical_salinity"],
        ds.cf["depth"],
        ds.cf["longitude"],
        ds.cf["latitude"],
    )
    pt0 = pt0_from_t(sa, ds.cf["sea_water_temperature"], ds.cf["depth"])
    pt0.attrs = {"standard_name": "sea_water_potential_temperature"}

    return pt0.cf.add_canonical_attributes()


def dms2d(degrees, minutes, seconds=0):
    return degrees + (minutes / 60.0) + (seconds / 3600.0)
