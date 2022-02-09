from functools import wraps

import cf_xarray  # noqa: F401
from geopy.distance import great_circle
from gsw import SA_from_SP, pt0_from_t
from xarray import DataArray, Dataset


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

        # Assign coordinates to variables (better integration with cf-xarray)
        for varname, variable in ds.data_vars.items():
            coordinates = []
            for coord in sum(ds.cf.coordinates.values(), []):
                if set(ds[coord].dims) <= set(variable.dims):
                    coordinates.append(coord)
            if coordinates:
                variable.attrs["coordinates"] = " ".join(coordinates)
            else:
                variable.attrs.pop("coordinates", None)

        # Sort
        if "distance" in ds:
            sortvar = "distance"
        elif "time" in ds and "time" not in ds.dims:
            sortvar = "time"
        else:
            sortvar = (
                "longitude"
                if ds.cf["longitude"].std() > ds.cf["latitude"].std()
                else "latitude"
            )
        ds = ds.sortby(sortvar)

        # Compute distance along transect
        if "distance" not in ds:
            lons = ds.cf["longitude"]
            lats = ds.cf["latitude"]
            distance = [0]
            for i in range(1, len(lons)):
                coords = ((lats[ind], lons[ind]) for ind in (i - 1, i))
                distance.append(distance[-1] + great_circle(*coords).km)
            distance = DataArray(
                distance,
                dims=lons.dims,
                attrs={"long_name": "distance along transect", "units": "km"},
            )
            ds = ds.assign_coords(distance=distance)

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
