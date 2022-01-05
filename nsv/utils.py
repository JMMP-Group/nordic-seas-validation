import cf_xarray  # noqa: F401
from gsw import SA_from_SP, pt0_from_t
from xarray import Dataset


def add_cf_attributes(ds: Dataset) -> Dataset:
    coords = ["longitude", "latitude", "distance"]
    ds = ds.set_coords(
        [var for var in ds.data_vars for coord in coords if coord in var]
    )
    return ds.cf.guess_coord_axis().cf.add_canonical_attributes()


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

def dms2d(Deg,Min,Sec):
    return Deg + (Min/60.) + (Sec/3600.)
